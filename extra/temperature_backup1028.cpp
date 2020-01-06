#include "stdafx.h"
#include "temperature.h"
#include "surface.h"
#include "laplacian.h"
#include "advection.h"
#include "fluxes.h"
#include "boundaries.h"
#include "interpolate.h"
#include "damping.h"

#define VAR(i,j,k) var[INDEX(i,j,k)]
#define INDEXB(i,j,k) ((i)*NY*NZ+(j)*NZ+(k))

#define INDEXT(i,j,k) ((i)*fNYfNZ+(j)*fNZ+(k))
#define SUBINDEX(i,j,k) ((i)*subNYNZ+(j)*NZ+(k))

#define DIVG(i,j,k) divg[(i)*fNYfNZ+(j)*fNZ+(k)]
#define ZETA(i,j,k) vort[(i)*fNYfNZ+(j)*fNZ+(k)]
#define ZETABAR(i,j,k) VORT[(i)*fNYfNZ+(j)*fNZ+(k)]

#define UFRIC_TEND(i,j,k) vort_ufric[(i)*fNYfNZ+(j)*fNZ+(k)]
#define VFRIC_TEND(i,j,k) vort_vfric[(i)*fNYfNZ+(j)*fNZ+(k)]

#define ALLOC(size) ((double*) calloc(size,sizeof(double)))

#define GET_VORT(u,v,i,j,k) ( (v(i,j,k)-v(i-1,j,k)) * one_d_dx - (u(i,j,k) - u(i,j-1,k)) * one_d_dy )
#define VORT_AT_SCALAR(z,i,j,k) 0.25*(z(i,j+1,k)+z(i+1,j+1,k)+z(i,j,k)+z(i+1,j,k)) 

#define TILT_AT_SCALAR(z,i,j,k) (0.125 * (z[INDEXT(i,j+1,k  )]+z[INDEXT(i+1,j+1,k  )]+ \
										  z[INDEXT(i,j  ,k  )]+z[INDEXT(i+1,j  ,k  )]+ \
										  z[INDEXT(i,j+1,k+1)]+z[INDEXT(i+1,j+1,k+1)]+ \
										  z[INDEXT(i,j  ,k+1)]+z[INDEXT(i+1,j  ,k+1)]) \
								)

#define TILT(u,v,w) -0.25* ( \
		(\
		 (v(i-1,j,k)+v(i,j,k)) - (v(i-1,j,k-1)+v(i,j,k-1))\
		) * DTZW(k)\
		*\
        (\
		 (w(i,j-1,k)+w(i,j,k)) - (w(i-1,j-1,k)+w(i-1,j,k))\
		) * one_d_dx\
	-\
		(\
		 (u(i,j-1,k)+u(i,j,k)) - (u(i,j-1,k-1)+u(i,j,k-1))\
		) * DTZW(k)\
		*\
        (\
		 (w(i-1,j,k)+w(i,j,k)) - (w(i-1,j-1,k)+w(i,j-1,k))\
		) * one_d_dy\
	   )

#define DOT_PRODUCT(var1,var2) ( var1.i*var2.i + var1.j*var2.j + var1.k*var2.k )

//---------------------------------------------------------------------------------------------------------

#define V_AT_W(v,i,j,k) 0.25*( v(i,j,k)+v(i,j+1,k)+v(i,j,k-1)+v(i,j+1,k-1) )
#define W_AT_V(v,i,j,k) 0.25*( v(i,j,k)+v(i,j,k+1)+v(i,j-1,k)+v(i,j-1,k+1) )

#define U_AT_W(v,i,j,k) 0.25*( v(i,j,k)+v(i+1,j,k)+v(i,j,k-1)+v(i+1,j,k-1) )
#define W_AT_U(v,i,j,k) 0.25*( v(i,j,k)+v(i,j,k+1)+v(i-1,j,k)+v(i-1,j,k+1) )

#define V_AT_U(v,i,j,k) 0.25 * (v(i,j,k)+v(i,j+1,k)+v(i-1,j,k)+v(i-1,j+1,k))
#define U_AT_V(v,i,j,k) 0.25 * (v(i,j,k)+v(i+1,j,k)+v(i,j-1,k)+v(i+1,j-1,k))

//---------------------------------------------------------------------------------------------------------

#define V_AT_W2(v,i,j,k) 0.25 * ( v[INDEX(i,j,k)]+v[INDEX(i,j+1,k)]+v[INDEX(i,j,k-1)]+v[INDEX(i,j+1,k-1)] )
#define W_AT_V2(v,i,j,k) 0.25 * ( v[INDEX(i,j,k)]+v[INDEX(i,j,k+1)]+v[INDEX(i,j-1,k)]+v[INDEX(i,j-1,k+1)] )

#define U_AT_W2(v,i,j,k) 0.25 * ( v[INDEX(i,j,k)]+v[INDEX(i+1,j,k)]+v[INDEX(i,j,k-1)]+v[INDEX(i+1,j,k-1)] )
#define W_AT_U2(v,i,j,k) 0.25 * ( v[INDEX(i,j,k)]+v[INDEX(i,j,k+1)]+v[INDEX(i-1,j,k)]+v[INDEX(i-1,j,k+1)] )

#define V_AT_U2(v,i,j,k) 0.25 * ( v[INDEX(i,j,k)]+v[INDEX(i,j+1,k)]+v[INDEX(i-1,j,k)]+v[INDEX(i-1,j+1,k)] )
#define U_AT_V2(v,i,j,k) 0.25 * ( v[INDEX(i,j,k)]+v[INDEX(i+1,j,k)]+v[INDEX(i,j-1,k)]+v[INDEX(i+1,j-1,k)] )

//---------------------------------------------------------------------------------------------------------

#define V_AT_I(v,i,j,k) 0.125 * ( v[INDEX(i,j,k  )]+v[INDEX(i,j+1,k  )]+v[INDEX(i-1,j,k  )]+v[INDEX(i-1,j+1,k  )] +	\
								  v[INDEX(i,j,k-1)]+v[INDEX(i,j+1,k-1)]+v[INDEX(i-1,j,k-1)]+v[INDEX(i-1,j+1,k-1)] )
#define W_AT_I(v,i,j,k) 0.125 * ( v[INDEX(i,j,k  )]+v[INDEX(i,j-1,k  )]+v[INDEX(i-1,j,k  )]+v[INDEX(i-1,j-1,k  )] +	\
								  v[INDEX(i,j,k+1)]+v[INDEX(i,j-1,k+1)]+v[INDEX(i-1,j,k+1)]+v[INDEX(i-1,j-1,k+1)] )

#define U_AT_J(v,i,j,k) 0.25 * ( v[INDEX(i,j,k)]+v[INDEX(i+1,j,k)]+v[INDEX(i,j,k-1)]+v[INDEX(i+1,j,k-1)] )
#define W_AT_J(v,i,j,k) 0.25 * ( v[INDEX(i,j,k)]+v[INDEX(i,j,k+1)]+v[INDEX(i-1,j,k)]+v[INDEX(i-1,j,k+1)] )

#define V_AT_K(v,i,j,k) 0.25 * ( v[INDEX(i,j,k)]+v[INDEX(i,j+1,k)]+v[INDEX(i-1,j,k)]+v[INDEX(i-1,j+1,k)] )
#define U_AT_K(v,i,j,k) 0.25 * ( v[INDEX(i,j,k)]+v[INDEX(i+1,j,k)]+v[INDEX(i,j-1,k)]+v[INDEX(i+1,j-1,k)] )

//---------------------------------------------------------------------------------------------------------

//#define DIABATIC(i,j,k) ( cond[INDEX(i,j,k)] + evap[INDEX(i,j,k)])
#define DIABATIC(i,j,k) ( evap[INDEX(i,j,k)] )

#define HCONV(u,v,i,j,k) ((u(i+1,j,k)-u(i,j,k))*one_d_dx + (v(i,j+1,k)-v(i,j,k))*one_d_dy) 

#define VCONV(w,i,j,k) ( rhow[k+1] * w(i,j,k+1) - rhow[k] * w(i,j,k) ) * ONE_D_DZ(k) * one_d_rhou[k]

#define UFLUX(u,var,i,j,k) 0.5 * ( u[INDEX(i+1,j,k)] * (var[INDEX(i+1,j,k)]+var[INDEX(i  ,j,k)]) \
						         - u[INDEX(i  ,j,k)] * (var[INDEX(i  ,j,k)]+var[INDEX(i-1,j,k)]) ) * dtx

#define VFLUX(v,var,i,j,k) 0.5 * ( v[INDEX(i,j+1,k)] * (var[INDEX(i,j+1,k)]+var[INDEX(i,j,k)]) -  \
							       v[INDEX(i,j,k)] * (var[INDEX(i,j,k)]+var[INDEX(i,j-1,k)]) ) * dty
																				   
#define WFLUX(w,var,i,j,k) 0.5 * ( rhow[k+1] * w[INDEX(i,j,k+1)] * (var[INDEX(i,j,k+1)]+var[INDEX(i,j,k  )]) 	\
						   		 - rhow[k  ] * w[INDEX(i,j,k  )] * (var[INDEX(i,j,k  )]+var[INDEX(i,j,k-1)]) ) * one_d_rhou[k] * DTZ(k)

#define UDIFF(i,j,k) u_friction[INDEX(i,j,k)]
#define VDIFF(i,j,k) v_friction[INDEX(i,j,k)]
#define WDIFF(i,j,k) w_friction[INDEX(i,j,k)]

void stretched_grid(double*,double*,double*,double*,double,int);

int grid_ratio = 1;

struct vector {

	double i;
	double j;
	double k;

};

//---------------------------------------------------
// Heat budget
//---------------------------------------------------
char heat_budget_filename[len];

const int tvar_count = 8;
const char tvar_names[tvar_count][20] = {"pp_hor_adv","pp_vert_adv","pb_hor_adv","pb_vert_adv",
										 "bp_hor_adv","bp_vert_adv","tdiff","m_diabatic"};

//const int tvar_count = 1;
//const char tvar_names[1][20] = {"cond"};

double *pp_hor_adv,*pp_vert_adv,*pb_hor_adv,*pb_vert_adv,*bp_hor_adv,*bp_vert_adv,*tdiff,*m_diabatic;

double *cond,*evap;

//---------------------------------------------------
// Potential energy budget
//---------------------------------------------------
char pe_budget_filename[len];

const int pe_var_count = 8;
const char pe_var_names[pe_var_count][20] = {"pe_pp_hor_adv","pe_pp_vert_adv","pe_pb_hor_adv","pe_pb_vert_adv",
											 "pe_bp_hor_adv","pe_bp_vert_adv","pe_tdiff","pe_m_diabatic"};

double *pe_pp_hor_adv,*pe_pp_vert_adv,*pe_pb_hor_adv,*pe_pb_vert_adv,*pe_bp_hor_adv,*pe_bp_vert_adv,*pe_tdiff,*pe_m_diabatic;
double *pe_cond;

//---------------------------------------------------
// Moisture budget
//---------------------------------------------------
char moisture_budget_filename[len];

const int qvar_count = 8;
const char qvar_names[qvar_count][20] = {"qpp_hor_adv","qpp_vert_adv","qpb_hor_adv","qpb_vert_adv",
										 "qbp_hor_adv","qbp_vert_adv","qdiff","qm_diabatic"};

double *qpp_hor_adv,*qpp_vert_adv,*qpb_hor_adv,*qpb_vert_adv,*qbp_hor_adv,*qbp_vert_adv,*qdiff,*qm_diabatic;

//---------------------------------------------------
// Omega equation
//---------------------------------------------------
double *omega,*Q1,*Q2;

const int ovar_count = 3;

const char ovar_names[ovar_count][20] = {"omega","Q1","Q2"};

//---------------------------------------------------
// Vorticity budget
//---------------------------------------------------
const int vvar_count = 17;

const char vvar_names[vvar_count][20] = {"VORT","vort",
					"vort_hor_adv","vort_HOR_ADV","VORT_hor_adv",
					"vort_ver_adv","VORT_ver_adv","vort_VER_ADV",
					"vort_hor_con","cor_hor_con" ,"VORT_hor_con",
					"vort_tilt","VORT_tilt" ,"vort_TILT",
					"cor_hor_adv","vort_fric","divg"};
					


double *vort,*VORT,
		*vort_hor_adv,*vort_HOR_ADV,*VORT_hor_adv,
		*vort_ver_adv,*VORT_ver_adv,*vort_VER_ADV,
		*vort_hor_con,*cor_hor_con,*VORT_hor_con,
		*vort_tilt,*VORT_tilt,*vort_TILT,
		*cor_hor_adv,*vort_ufric,*vort_vfric,
		*divg,*fdivg;

//---------------------------------------------------
// PV budget
//---------------------------------------------------
const int pvvar_count = 7;

double *pv,*pv_diabatic,*pv_advect,*PV_advect,*pv_diff,*pv_ADVECT;
double *pv_base;
struct vector *avortVector;

const char pvvar_names[pvvar_count][20] = {"pv","pv_diabatic","pv_advect","PV_advect","pv_diff","pv_ADVECT","PV_base"};
//double *pvvars[pvvar_count] = 			  { pv , pv_diabatic , pv_advect , PV_advect };



//---------------------------------------------------
// PV tracer
//---------------------------------------------------
const int pv_tracer_var_count = 3;

double *pv_trace_diabatic1,*pv_trace_diabatic2,*pv_trace_diabatic3;
double *pv_trace_advect1,*pv_trace_advect2,*pv_trace_advect3;
double *pv_trace_initial1,*pv_trace_initial2,*pv_trace_initial3;
struct cell *pvcell,*pbcell;

double *u_bound,*v_bound,*w_bound,*t_bound;
double *base_advect;

const char pv_tracer_names[pv_tracer_var_count][20] = {"pv","PV_advect","pv_diabatic"};

//---------------------------------------------------
// Shared
//---------------------------------------------------
double *smallArray;
bool smallArrayIsInitialized = false;

int subNX = NX/grid_ratio;
int subNY = NY/grid_ratio;
int subNYNZ = subNY*NZ;

const double onetwelfth = 1.0/12.0;
const double one_d_60 = 1.0/60.0;
const double one_d_12 = 1.0/12.0;

double average_abs_diff(double *var1,double *var2,int il,int ih,int jl,int jh);
void write_tvar_to_file(const char *,const char *,double *,bool);
void initialize_basic_state_PV(int il,int ih,int jl,int jh);
void initialize_vorticity();
void initialize_moisture(int);
void initialize_temperature(int);
void initialize_pe(int);
void initialize_pv();
void initialize_pv_tracer();
void basic_state_vorticity(const char*,int,int,int,int);

void calculate_heat_budget(int,int,int,int);
void calculate_moisture_budget(int,int,int,int);
void calculate_vorticity_budget(int,int,int,int);
void calculate_pv_budget(int,int,int,int);

void write_pevars_to_file(const char *);
void write_tvars_to_file(const char *);
void write_qvars_to_file(const char *);
void write_vvars_to_file(const char *);
void write_PVvars_to_file(const char *);
void write_pv_tracer_to_file(const char *);

void pv_tracer_advect(double step,int il,int ih,int jl,int jh,int stage,int last_stage);

void create_budget_file(const char *,int);

/*********************************************************************
* Creates an individual file for each requested budget and allocates
* all required memory.
**********************************************************************/
void initialize_budgets(){

	/*********************************************
	* If heat budget
	**********************************************/
	if(HEAT_BUDGET){

		if(argcount==4){ add_commandline_arguments(heat_budget_filename,out_heat_budget_filename);} 
		else {			 strcpy(heat_budget_filename,out_heat_budget_filename);}
	
		if(rank==0){
		
			printf("Heat budget file name = %s\n",heat_budget_filename);
			
			create_budget_file(heat_budget_filename,0);
		}
	
		initialize_temperature(fNX*fNY*fNZ);
	}

	/*********************************************
	* If potential energy budget
	**********************************************/
	if(PE_BUDGET){

		if(argcount==4){ add_commandline_arguments(pe_budget_filename,out_pe_budget_filename);} 
		else {			 strcpy(pe_budget_filename,out_pe_budget_filename);}
	
		if(rank==0){
		
			printf("Potential energy budget file name = %s\n",pe_budget_filename);
			
			create_budget_file(pe_budget_filename,6);
		}
	
		initialize_pe(fNX*fNY*fNZ);
	}

	/*********************************************
	* If moisture budget
	**********************************************/
	if(MOISTURE_BUDGET){

		if(argcount==4){ add_commandline_arguments(moisture_budget_filename,out_moisture_budget_filename);} 
		else {			 strcpy(moisture_budget_filename,out_moisture_budget_filename);}
	
		if(rank==0){
						
			printf("Moisture budget file name = %s\n",moisture_budget_filename);
			
			create_budget_file(moisture_budget_filename,5);
		}
	
		initialize_moisture(fNX*fNY*fNZ);
	}

	/*********************************************
	* If vorticity budget
	**********************************************/
	if(VORTICITY_BUDGET){
	
		if(rank==0){ create_budget_file(vorticity_budget_filename,1);}
	
		initialize_vorticity();
		basic_state_vorticity(vorticity_budget_filename,1,fNX-1,1,fNY-1);
	}

	/*********************************************
	* If heat budget
	**********************************************/
	if(PV_BUDGET){ 
	
		if(rank==0){ create_budget_file(pv_budget_filename,2);}
	
		initialize_pv();
	}

	/*********************************************
	* If PV tracer
	**********************************************/
	if(PV_TRACER){ 
	
		if(rank==0){ create_budget_file(pv_tracer_filename,4);}
	
		initialize_pv_tracer();
	}
}

/*********************************************************************
*
*
**********************************************************************/
void calculate_budgets(int substep,double *steps){

	if(substep==2 && (HEAT_BUDGET || PE_BUDGET) ){ calculate_heat_budget(1,fNX-1,1,fNY-1);}
	if(substep==2 && MOISTURE_BUDGET){ calculate_moisture_budget(1,fNX-1,1,fNY-1);}
	if(substep==2 && VORTICITY_BUDGET){ calculate_vorticity_budget(1,fNX-1,1,fNY-1);}
	if(substep==0 && PV_BUDGET){ calculate_pv_budget(1,fNX-1,1,fNY-1);}
	if(PV_TRACER){ pv_tracer_advect(steps[substep],3,fNX-3,3,fNY-3,substep,2);}
}

/*********************************************************************
*
*
**********************************************************************/
void write_budgets_to_file(){

	if(PE_BUDGET && bigcounter % (outfilefreq) == 0){ write_pevars_to_file(pe_budget_filename);}	
	if(HEAT_BUDGET && bigcounter % (outfilefreq) == 0){ write_tvars_to_file(heat_budget_filename);}
	if(MOISTURE_BUDGET && bigcounter % (outfilefreq) == 0){ write_qvars_to_file(moisture_budget_filename);}
	if(VORTICITY_BUDGET && bigcounter % (outfilefreq) == 0){ write_vvars_to_file(vorticity_budget_filename);}
	if(PV_BUDGET && bigcounter % (outfilefreq) == 0){ write_PVvars_to_file(pv_budget_filename);}
	if(PV_TRACER && bigcounter % (outfilefreq) == 0){ write_pv_tracer_to_file(pv_tracer_filename);}
	
}

/*********************************************************************
*
*
**********************************************************************/
void initialize_pv(){
	
	int size = fNX*fNY*fNZ;
	
	//for(int i=0;i<pvvar_count;i++){ pvvars[i] = ALLOC(size);}
	
	pv = ALLOC(size);
	pv_diff = ALLOC(size);
	pv_diabatic = ALLOC(size);
	pv_advect = ALLOC(size);
	PV_advect = ALLOC(size);	
	pv_ADVECT = ALLOC(size);
	pv_base = ALLOC(size);
	//cond = ALLOC(size);
	evap = ALLOC(size);
	
	initialize_basic_state_PV(1,fNX-1,1,fNY-1);
	
	exchange(pv_base);
	
	write_tvar_to_file(pv_budget_filename,"PV_base",pv_base,false);
	
}

/*********************************************************************
* Initialize variables used in temperature budget
*
**********************************************************************/
void initialize_temperature(int size){
	
	pp_hor_adv  = ALLOC(size);
	pp_vert_adv = ALLOC(size);
	pb_hor_adv  = ALLOC(size);
	pb_vert_adv = ALLOC(size);
	bp_hor_adv  = ALLOC(size);
	bp_vert_adv = ALLOC(size);
	tdiff 		= ALLOC(size);	
	m_diabatic 	= ALLOC(size);

#if 0

	if(!smallArrayIsInitialized){
	
		smallArray = (double*) calloc((NX/grid_ratio)*(NY/grid_ratio)*NZ,sizeof(double));
		
		smallArrayIsInitialized = true;
	}
#endif
	
}

/*********************************************************************
* Initialize variables used in temperature budget
*
**********************************************************************/
void initialize_pe(int size){
	
	pe_pp_hor_adv  = ALLOC(size);
	pe_pp_vert_adv = ALLOC(size);
	pe_pb_hor_adv  = ALLOC(size);
	pe_pb_vert_adv = ALLOC(size);
	pe_bp_hor_adv  = ALLOC(size);
	pe_bp_vert_adv = ALLOC(size);
	pe_tdiff 	   = ALLOC(size);	
	pe_m_diabatic  = ALLOC(size);
	
}

/*********************************************************************
* Initialize variables used in temperature budget
*
**********************************************************************/
void initialize_moisture(int size){
	
	qpp_hor_adv  = ALLOC(size);
	qpp_vert_adv = ALLOC(size);
	qpb_hor_adv  = ALLOC(size);
	qpb_vert_adv = ALLOC(size);
	qbp_hor_adv  = ALLOC(size);
	qbp_vert_adv = ALLOC(size);
	qdiff 		= ALLOC(size);	
	qm_diabatic = ALLOC(size);

}

/*********************************************************************
*
*
**********************************************************************/
void initialize_vorticity(){
	
	int size = fNX*fNY*fNZ;
	
	vort = ALLOC(size);
	VORT = ALLOC(size);
	
	vort_hor_adv = ALLOC(size);
	VORT_hor_adv = ALLOC(size);
	vort_HOR_ADV = ALLOC(size);
	
	vort_ver_adv = ALLOC(size);
	VORT_ver_adv = ALLOC(size);
	vort_VER_ADV = ALLOC(size);
	cor_hor_adv = ALLOC(size);
	
	vort_hor_con = ALLOC(size);
	cor_hor_con = ALLOC(size);
	VORT_hor_con = ALLOC(size);
		
	vort_tilt = ALLOC(size);
	VORT_tilt = ALLOC(size);
	vort_TILT = ALLOC(size);

	vort_ufric = ALLOC(size);
	vort_vfric = ALLOC(size);

	divg = ALLOC(size);
	//cond = ALLOC(size);
	
	if(rank==0){
		
		if(PERIODIC_BOUNDARIES){ 
			init_laplacian_solver_periodic_EW_zerograd_NS(NX,NY);
		} else {
			init_laplacian_solver_real(NX,NY);
		}
		
	} 
	
	if(false && !smallArrayIsInitialized){
		
		smallArray = (double*) calloc((NX/grid_ratio)*(NY/grid_ratio)*NZ,sizeof(double));
		
		smallArrayIsInitialized = true;
	}
}

/*********************************************************************
*
*
**********************************************************************/
void initialize_pv_tracer(){
	
	int size = fNX*fNY*fNZ;
	
	pv_trace_diabatic1 = ALLOC(size);
	pv_trace_diabatic2 = ALLOC(size);
	pv_trace_diabatic3 = ALLOC(size);
	pv_trace_advect1 = ALLOC(size);
	pv_trace_advect2 = ALLOC(size);
	pv_trace_advect3 = ALLOC(size);	
	pv_trace_initial1 = ALLOC(size);
	pv_trace_initial2 = ALLOC(size);
	pv_trace_initial3 = ALLOC(size);
	
	u_bound = ALLOC(size);
	v_bound = ALLOC(size);
	w_bound = ALLOC(size);
	t_bound = ALLOC(size);
	
	base_advect = ALLOC(size);
	
	avortVector = (vector*)calloc(size,sizeof(vector));
	
	if(!PV_BUDGET){
		
		pv = ALLOC(size);
		evap = ALLOC(size);
		pv_base = ALLOC(size);
	}
	
	
	pvcell  = (cell*)calloc(fNY*fNZ,sizeof(cell));
	pbcell  = (cell*)calloc(fNY*fNZ,sizeof(cell));
	
	initialize_basic_state_PV(1,fNX-1,1,fNY-1);
	
	exchange(pv_base);
	
	p_mirror_boundaries(pv_base,east,west,north,south,halo_buffer,fNX,fNY);

}

/*********************************************************************
*
*
**********************************************************************/
void initialize_basic_state_PV(int il,int ih,int jl,int jh){
	
	vector Zeta,dTheta;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		//----------------------------------------------------------
		// Potential temperature gradient vector
		//----------------------------------------------------------		
		dTheta.i = 0.5 * (THBAR(i+1,j,k) - THBAR(i-1,j,k)) * one_d_dx;
		dTheta.j = 0.5 * (THBAR(i,j+1,k) - THBAR(i,j-1,k)) * one_d_dy;				   
		dTheta.k = 0.5 * (THBAR(i,j,k+1) - THBAR(i,j,k-1) + tb[k+1] - tb[k-1] ) * ONE_D_DZW(k);

		//----------------------------------------------------------
		// Absolute vorticity vector
		//----------------------------------------------------------
		Zeta.i = (W_AT_V(WBAR,i,j+1,k) - W_AT_V(WBAR,i,j,k)) * one_d_dy 	- (V_AT_W(VBAR,i,j,k+1) - V_AT_W(VBAR,i,j,k)) * ONE_D_DZ(k);
		Zeta.j = (U_AT_W(UBAR,i,j,k+1) - U_AT_W(UBAR,i,j,k)) * ONE_D_DZ(k)  - (W_AT_U(WBAR,i+1,j,k) - W_AT_U(WBAR,i,j,k)) * one_d_dx;
		Zeta.k = (V_AT_U(VBAR,i+1,j,k) - V_AT_U(VBAR,i,j,k)) * one_d_dx 	- (U_AT_V(UBAR,i,j+1,k) - U_AT_V(UBAR,i,j,k)) * one_d_dy + FC(j);	

		//----------------------------------------------------------
		// Basic state potential vorticity
		//----------------------------------------------------------
		pv_base[INDEX(i,j,k)] = DOT_PRODUCT(Zeta,dTheta) * one_d_rhou[k];
		
	}}}
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
		
		pv_base[INDEX(i,j,0)] = pv_base[INDEX(i,j,1)];
		pv_base[INDEX(i,j,NZ-1)] = pv_base[INDEX(i,j,NZ-2)];		
	}}
	
}

/*********************************************************************
*
*
**********************************************************************/
void basic_state_vorticity(const char *myfilename,int il,int ih,int jl,int jh){
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
	
		VORT[INDEXT(i,j,k)] =  GET_VORT(UBAR,VBAR,i,j,k);
		
	}}}
	
	exchange(VORT);
	
	for(int i=1;i<fNX-1;i++){
	for(int j=1;j<fNY-1;j++){
	for(int k=1;k<NZ-1;k++){
		
		vort_hor_adv[INDEXT(i,j,k)] = VORT_AT_SCALAR(ZETABAR,i,j,k);
	}}}
	
	write_tvar_to_file(myfilename,"VORT",vort_hor_adv,true);
}

/*********************************************************************
*
*
**********************************************************************/
void set_diabatic_heating(double * condensation,double * evaporation){
	
	int size = fNX*fNY*fNZ;
	
	for(int i=0;i<size;i++){
	
		cond[i] += condensation[i];
		evap[i] += evaporation[i];
		
	}
	
}

/*********************************************************************
*
*
**********************************************************************/
double budget_scalar_diff_6th(double * var,int i,int j,int k){
	
	double u_mid = 0.5*( UBAR(i,j,k) + UBAR(i+1,j,k) + U(i,j,k) + U(i+1,j,k) );
	double v_mid = 0.5*( VBAR(i,j,k) + VBAR(i,j+1,k) + V(i,j,k) + V(i,j+1,k) );
	double w_mid = 0.5*( WBAR(i,j,k) + WBAR(i,j,k+1) + W(i,j,k) + W(i,j,k+1) );
	
	double diff_z;
	
	if(k>1 && k<NZ-2){
		diff_z = -diffuse_k_4th(var,i,j,k);
	} else {
		diff_z = 0;
	}

	return
		one_d_60*fabs(u_mid) * one_d_dx    * diffuse_i_6th(var,i,j,k) +
		one_d_60*fabs(v_mid) * one_d_dy	   * diffuse_j_6th(var,i,j,k) +
		one_d_12*fabs(w_mid) * ONE_D_DZ(k) * diff_z
	;
}

/*********************************************************************
* Calculate terms in temperature budget
*
* il,ih,jl,jh - array indices
**********************************************************************/
void calculate_heat_budget(int il,int ih,int jl,int jh){
	
	int ind;

	double pe_mult;
	double pp_hor_adv0,bp_hor_adv0,bp_vert_adv0,pb_hor_adv0,pb_vert_adv0,pp_vert_adv0,tdiff0;

	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-2;k++){
		
		ind = INDEX(i,j,k);
		//--------------------------------------------------------------
		// Hoizontal advection of perturbation temperature by perturbation wind,
		// i.e. the non-linear term
		//--------------------------------------------------------------
		pp_hor_adv0 = - ( UFLUX(us,ths,i,j,k) + VFLUX(vs,ths,i,j,k) ) + TH(i,j,k) * HCONV(U,V,i,j,k) * dt;

		//--------------------------------------------------------------
		// Horizontal advection of perturbation temperature by basic state wind
		//--------------------------------------------------------------		
		bp_hor_adv0 = - ( UFLUX(m_ubar,ths,i,j,k) + VFLUX(m_vbar,ths,i,j,k) ) + TH(i,j,k) * HCONV(UBAR,VBAR,i,j,k) * dt;

		//--------------------------------------------------------------
		// Horizontal advection of basic state temperature by perturbation wind
		//--------------------------------------------------------------		
		pb_hor_adv0 = - ( UFLUX(us,m_thbar,i,j,k) + VFLUX(vs,m_thbar,i,j,k) ) + THBAR(i,j,k) * HCONV(U,V,i,j,k) * dt;

		//--------------------------------------------------------------
		// Vertical advection of perturbation temperature by basic state wind
		//--------------------------------------------------------------		
		bp_vert_adv0 = - WFLUX(m_wbar,ths,i,j,k) + TH(i,j,k) * VCONV(WBAR,i,j,k) * dt;

		//--------------------------------------------------------------
		// Vertical advection of basic state temperature by perturbation wind
		//--------------------------------------------------------------	
		pb_vert_adv0 = -0.5*( rhow[k+1]*W(i,j,k+1)*(tb[k+1]-tb[k  ]) 
								 + rhow[k  ]*W(i,j,k  )*(tb[k  ]-tb[k-1])) 
							* one_d_rhou[k] * DTZ(k)
												 
							- WFLUX(ws,m_thbar,i,j,k) + THBAR(i,j,k) * VCONV(W,i,j,k) * dt;

		//--------------------------------------------------------------
		// Vertical advection of perturbation temperature by perturbation wind
		//--------------------------------------------------------------	
		pp_vert_adv0 = -WFLUX(ws,ths,i,j,k) + TH(i,j,k) * VCONV(W,i,j,k) * dt;
		
		//--------------------------------------------------------------
		// Diffusion from turbulence
		//--------------------------------------------------------------	
		tdiff0 = ( t_diffusion[ind] + budget_scalar_diff_6th(ths,i,j,k) ) * dt;

		//--------------------------------------------------------------
		// Add each term to output arrays for heat budget
		//--------------------------------------------------------------
		#if HEAT_BUDGET
			pp_hor_adv[ind] += 	pp_hor_adv0;
			bp_hor_adv[ind] += 	bp_hor_adv0;
			pb_hor_adv[ind] += pb_hor_adv0;
			bp_vert_adv[ind] += bp_vert_adv0;
			pb_vert_adv[ind] += pb_vert_adv0;
			pp_vert_adv[ind] += pp_vert_adv0;
			tdiff[ind] += tdiff0;
		#endif
		//--------------------------------------------------------------
		// Add each term to output arrays for potential energy budget
		//--------------------------------------------------------------		
		#if PE_BUDGET		
			pe_mult = grav * (ths[ind] / tb[k]) * (ZW(k+1)-ZW(k)) / (tbw[k+1]-tbw[k]);
			
			pe_pp_hor_adv[ind] += pp_hor_adv0 * pe_mult;
			pe_bp_hor_adv[ind] += bp_hor_adv0 * pe_mult;
			pe_pb_hor_adv[ind] += pb_hor_adv0 * pe_mult;
			pe_bp_vert_adv[ind] += bp_vert_adv0 * pe_mult;
			pe_pb_vert_adv[ind] += pb_vert_adv0 * pe_mult;
			pe_pp_vert_adv[ind] += pp_vert_adv0 * pe_mult;
			pe_tdiff[ind] += tdiff0 * pe_mult;
		#endif
	}}}
	
}

/*********************************************************************
* Calculate terms in temperature budget
*
* il,ih,jl,jh - array indices
**********************************************************************/
void calculate_moisture_budget(int il,int ih,int jl,int jh){
	
	int ind;

	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		
		ind = INDEX(i,j,k);
		
		//--------------------------------------------------------------
		// Hoizontal advection of perturbation temperature by perturbation wind,
		// i.e. the non-linear term
		//--------------------------------------------------------------
		qpp_hor_adv[ind] += - ( UFLUX(us,qvs,i,j,k) + VFLUX(vs,qvs,i,j,k) ) + QV(i,j,k) * HCONV(U,V,i,j,k) * dt;

		//--------------------------------------------------------------
		// Horizontal advection of perturbation temperature by basic state wind
		//--------------------------------------------------------------		
		qbp_hor_adv[ind] += - ( UFLUX(m_ubar,qvs,i,j,k) + VFLUX(m_vbar,qvs,i,j,k) ) + QV(i,j,k) * HCONV(UBAR,VBAR,i,j,k) * dt;

		//--------------------------------------------------------------
		// Vertical advection of perturbation temperature by basic state wind
		//--------------------------------------------------------------		
		qbp_vert_adv[ind] += - WFLUX(m_wbar,qvs,i,j,k) + QV(i,j,k) * VCONV(WBAR,i,j,k) * dt;

		//--------------------------------------------------------------
		// Horizontal advection of basic state temperature by perturbation wind
		//--------------------------------------------------------------		
		qpb_hor_adv[ind] += - ( UFLUX(us,m_qbar,i,j,k) + VFLUX(vs,m_qbar,i,j,k) ) + QBAR(i,j,k) * HCONV(U,V,i,j,k) * dt;

		//--------------------------------------------------------------
		// Vertical advection of basic state temperature by perturbation wind
		//--------------------------------------------------------------	
		qpb_vert_adv[ind] += -0.5*(rhow[k+1]*W(i,j,k+1)*(qb[k+1]-qb[k  ]) 
								 + rhow[k  ]*W(i,j,k  )*(qb[k  ]-qb[k-1])) 
							* one_d_rhou[k] * DTZ(k)
												 
							- WFLUX(ws,m_qbar,i,j,k) + QBAR(i,j,k) * VCONV(W,i,j,k) * dt;

		//--------------------------------------------------------------
		// Vertical advection of perturbation temperature by perturbation wind
		//--------------------------------------------------------------	
		qpp_vert_adv[ind] += -WFLUX(ws,qvs,i,j,k) + QV(i,j,k) * VCONV(W,i,j,k) * dt;
		
		//--------------------------------------------------------------
		// Diffusion from turbulence
		//--------------------------------------------------------------	
		qdiff[ind] += ( qv_diffusion[ind] + budget_scalar_diff_6th(qvs,i,j,k) ) * dt;

	}}}
	
}

/*********************************************************************
*
*
**********************************************************************/
void calculate_vorticity_budget(int il,int ih,int jl,int jh){
	
	double vort_at_scalar,VORT_at_scalar,hconv;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		vort[INDEXT(i,j,k)] = GET_VORT(U,V,i,j,k);
		divg[INDEXT(i,j,k)] = HCONV(U,V,i,j,k);
	}}}
	
	exchange(vort);
	
	int ind = 0;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
	
		//-------------------------------------------------------------
		// Vorticity
		//-------------------------------------------------------------
		
		vort_at_scalar = VORT_AT_SCALAR(ZETA,i,j,k);
		
		VORT_at_scalar = VORT_AT_SCALAR(ZETABAR,i,j,k);
		
		ind = INDEXT(i,j,k);
		
		//-------------------------------------------------------------
		// Horizontal vorticity flux convergence
		//-------------------------------------------------------------		
		vort_hor_adv[ind] += - ((U(i+1,j,k)*(ZETA(i+1,j+1,k)+ZETA(i+1,j,k)) - U(i,j,k)*(ZETA(i,j+1,k)+ZETA(i,j,k)))
							 +  (V(i,j+1,k)*(ZETA(i,j+1,k)+ZETA(i+1,j+1,k)) - V(i,j,k)*(ZETA(i,j,k)+ZETA(i+1,j,k))))*0.5*dtx;
		
		vort_HOR_ADV[ind] += - ((UBAR(i+1,j,k)*(ZETA(i+1,j+1,k)+ZETA(i+1,j,k)) - UBAR(i,j,k)*(ZETA(i,j+1,k)+ZETA(i,j,k)))
							 +  (VBAR(i,j+1,k)*(ZETA(i,j+1,k)+ZETA(i+1,j+1,k)) - VBAR(i,j,k)*(ZETA(i,j,k)+ZETA(i+1,j,k))))*0.5*dtx;

		VORT_hor_adv[ind] += - ((U(i+1,j,k)*(ZETABAR(i+1,j+1,k)+ZETABAR(i+1,j  ,k)) - U(i,j,k)*(ZETABAR(i,j+1,k)+ZETABAR(i  ,j,k)))
							 +  (V(i,j+1,k)*(ZETABAR(i  ,j+1,k)+ZETABAR(i+1,j+1,k)) - V(i,j,k)*(ZETABAR(i,j  ,k)+ZETABAR(i+1,j,k))))*0.5*dtx;

		//-------------------------------------------------------------
		// Vertical vorticity flux convergence
		//-------------------------------------------------------------		
		vort_ver_adv[ind] += -0.5*DTZ(k) * (  rhow[k+1] * W(i,j,k+1) * ( vort_at_scalar + VORT_AT_SCALAR(ZETA,i,j,k+1) )
										   -  rhow[k  ] * W(i,j,k  ) * ( vort_at_scalar + VORT_AT_SCALAR(ZETA,i,j,k-1) )
											) * one_d_rhou[k];

		vort_VER_ADV[ind] += -0.5*DTZ(k) * (  rhow[k+1] * WBAR(i,j,k+1) * ( vort_at_scalar + VORT_AT_SCALAR(ZETA,i,j,k+1) )
										   -  rhow[k  ] * WBAR(i,j,k  ) * ( vort_at_scalar + VORT_AT_SCALAR(ZETA,i,j,k-1) )
											) * one_d_rhou[k];
		
		VORT_ver_adv[ind] += -0.5*DTZ(k) * (  rhow[k+1] * W(i,j,k+1) * ( VORT_at_scalar + VORT_AT_SCALAR(ZETABAR,i,j,k+1) )
										   -  rhow[k  ] * W(i,j,k  ) * ( VORT_at_scalar + VORT_AT_SCALAR(ZETABAR,i,j,k-1) )
											) * one_d_rhou[k];
		
		//-------------------------------------------------------------
		// Vorticity convergence (i.e. stretching)
		//-------------------------------------------------------------
		hconv = HCONV(U,V,i,j,k);
		
		vort_hor_con[ind] += -vort_at_scalar * hconv * dt;

		vort_HOR_ADV[ind] += -vort_at_scalar * HCONV(UBAR,VBAR,i,j,k) * dt;

		VORT_hor_con[ind] += -VORT_at_scalar * hconv * dt;
		
		cor_hor_con[ind] += -FC(j) * hconv * dt;
				
		//-------------------------------------------------------------
		// Vorticity tilting
		//-------------------------------------------------------------
		vort_tilt[ind] += TILT(U,V,W);
		
		VORT_tilt[ind] += TILT(UBAR,VBAR,W);
		
		vort_TILT[ind] += TILT(U,V,WBAR);
		
		//-------------------------------------------------------------
		// Advection of planetary vorticity
		//-------------------------------------------------------------
		cor_hor_adv[ind] += -0.5 * (V(i,j,k)+V(i,j+1,k)) * DFDY(j) * dt;
			
	}}}
	
}

/*********************************************************************
*
* 
**********************************************************************/
double u_at_v(double *var,int i,int j,int k){
	
	double p[4][4];
	//if(i-1<0 || i+2 > fNX-1 || j+1>fNY-1 || j-2<0){printf("here!\n");}	
	for(int a=0;a<4;a++){
	
		p[a][0] = VAR(i-1+a,j-2,k);
		p[a][1] = VAR(i-1+a,j-1,k);
		p[a][2] = VAR(i-1+a,j  ,k);	
		p[a][3] = VAR(i-1+a,j+1,k);
	}
	
	return BicubicInterpolate2(p,0.5,0.5);
}

/*********************************************************************
*
* 
**********************************************************************/
double v_at_u(double *var,int i,int j,int k){
	
	double p[4][4];
	//if(i-2<0 || i+1 > fNX-1 || j+2>fNY-1 || j-1<0){printf("here!\n");}
	for(int a=0;a<4;a++){
	
		p[a][0] = VAR(i-2+a,j-1,k);
		p[a][1] = VAR(i-2+a,j  ,k);
		p[a][2] = VAR(i-2+a,j+1,k);	
		p[a][3] = VAR(i-2+a,j+2,k);
	}
	
	return BicubicInterpolate2(p,0.5,0.5);
}

/*********************************************************************
* At scalar points
*
**********************************************************************/
void vorticity_gradient_vector(double *u,double *v,double *w,int i,int j,int k,struct vector *vec){
	
	vec->i = (W_AT_V2(w,i,j+1,k) - W_AT_V2(w,i,j,k)) * one_d_dy 	 - 	(V_AT_W2(v,i,j,k+1) - V_AT_W2(v,i,j,k)) * ONE_D_DZ(k);
	vec->j = (U_AT_W2(u,i,j,k+1) - U_AT_W2(u,i,j,k)) * ONE_D_DZ(k)   -  (W_AT_U2(w,i+1,j,k) - W_AT_U2(w,i,j,k)) * one_d_dx;
	vec->k = (V_AT_U2(v,i+1,j,k) - V_AT_U2(v,i,j,k)) * one_d_dx 	 - 	(U_AT_V2(u,i,j+1,k) - U_AT_V2(u,i,j,k)) * one_d_dy;	
}


/*********************************************************************
*
*
**********************************************************************/
void scalar_gradient_vector(double *var,int i,int j,int k,struct vector *vec){
	
	vec->i = 0.5 * ( var[INDEX(i+1,j,k)] - var[INDEX(i-1,j,k)]) * one_d_dx;
	vec->j = 0.5 * ( var[INDEX(i,j+1,k)] - var[INDEX(i,j-1,k)]) * one_d_dy;
	vec->k = 0.5 * ( var[INDEX(i,j,k+1)] - var[INDEX(i,j,k-1)]) * ONE_D_DZ(k);
}

/*********************************************************************
*
*
**********************************************************************/
void scalar_gradient_vector2(double *var,int i,int j,int k,struct vector *vec){
	
	vec->i = onetwelfth * ( var[INDEX(i-2,j,k)] - 8.0*var[INDEX(i-1,j,k)] + 8.0*var[INDEX(i+1,j,k)] - var[INDEX(i+2,j,k)]) * one_d_dx;
	vec->j = onetwelfth * ( var[INDEX(i,j-2,k)] - 8.0*var[INDEX(i,j-1,k)] + 8.0*var[INDEX(i,j+1,k)] - var[INDEX(i,j+2,k)]) * one_d_dy;
	vec->k = onetwelfth * ( var[INDEX(i,j,k-2)] - 8.0*var[INDEX(i,j,k-1)] + 8.0*var[INDEX(i,j,k+1)] - var[INDEX(i,j,k+2)]) * ONE_D_DZ(k);
}


/*********************************************************************
*
*
**********************************************************************/
void calculate_pv_budget(int il,int ih,int jl,int jh){
	
	vector zeta,Zeta,zetaF;
	vector dtheta,dTheta;
	vector dthetaDot;
	
	int ind;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		ind = INDEX(i,j,k);

		//----------------------------------------------------------
		// Diabatic heating gradient vector
		//----------------------------------------------------------		
		//dthetaDot.i = 0.5 * ( DIABATIC(i+1,j,k) - DIABATIC(i-1,j,k)) * one_d_dx;
		//dthetaDot.j = 0.5 * ( DIABATIC(i,j+1,k) - DIABATIC(i,j-1,k)) * one_d_dy;
		//dthetaDot.k = 0.5 * ( DIABATIC(i,j,k+1) - DIABATIC(i,j,k-1)) * ONE_D_DZ(k);

		scalar_gradient_vector(&DIABATIC(0,0,0),i,j,k,&dthetaDot);
		
		scalar_gradient_vector(&TH(0,0,0),i,j,k,&dtheta);

		scalar_gradient_vector(&THBAR(0,0,0),i,j,k,&dTheta);

		dTheta.k += 0.5 * (tb[k+1] - tb[k-1] ) * ONE_D_DZ(k);

		//----------------------------------------------------------
		// Potential temperature gradient vector
		//----------------------------------------------------------		
		//dtheta.i = 0.5 * (TH(i+1,j,k) - TH(i-1,j,k)) * one_d_dx;
		//dtheta.j = 0.5 * (TH(i,j+1,k) - TH(i,j-1,k)) * one_d_dy;
		//dtheta.k = 0.5 * (TH(i,j,k+1) - TH(i,j,k-1)) * ONE_D_DZ(k);

		//dTheta.i = 0.5 * (THBAR(i+1,j,k) - THBAR(i-1,j,k)) * one_d_dx;
		//dTheta.j = 0.5 * (THBAR(i,j+1,k) - THBAR(i,j-1,k)) * one_d_dy;				   
		//dTheta.k = 0.5 * (THBAR(i,j,k+1) - THBAR(i,j,k-1) + tb[k+1] - tb[k-1] ) * ONE_D_DZ(k);

		//----------------------------------------------------------
		// Absolute vorticity vector
		//----------------------------------------------------------
		zeta.i = (W_AT_V(W,i,j+1,k) - W_AT_V(W,i,j,k)) * one_d_dy 	 - 	(V_AT_W(V,i,j,k+1) - V_AT_W(V,i,j,k)) * ONE_D_DZ(k);
		zeta.j = (U_AT_W(U,i,j,k+1) - U_AT_W(U,i,j,k)) * ONE_D_DZ(k) -  (W_AT_U(W,i+1,j,k) - W_AT_U(W,i,j,k)) * one_d_dx;
		zeta.k = (V_AT_U(V,i+1,j,k) - V_AT_U(V,i,j,k)) * one_d_dx 	 - 	(U_AT_V(U,i,j+1,k) - U_AT_V(U,i,j,k)) * one_d_dy;

		Zeta.i = (W_AT_V(WBAR,i,j+1,k) - W_AT_V(WBAR,i,j,k)) * one_d_dy 	- (V_AT_W(VBAR,i,j,k+1) - V_AT_W(VBAR,i,j,k)) * ONE_D_DZ(k);
		Zeta.j = (U_AT_W(UBAR,i,j,k+1) - U_AT_W(UBAR,i,j,k)) * ONE_D_DZ(k)  - (W_AT_U(WBAR,i+1,j,k) - W_AT_U(WBAR,i,j,k)) * one_d_dx;
		Zeta.k = (V_AT_U(VBAR,i+1,j,k) - V_AT_U(VBAR,i,j,k)) * one_d_dx 	- (U_AT_V(UBAR,i,j+1,k) - U_AT_V(UBAR,i,j,k)) * one_d_dy + FC(j);

		//----------------------------------------------------------
		//  Diffusional vorticity vector
		//----------------------------------------------------------
		zetaF.i = (W_AT_V(WDIFF,i,j+1,k) - W_AT_V(WDIFF,i,j,k)) * one_d_dy 	 
			  -   (V_AT_W(VDIFF,i,j,k+1) - V_AT_W(VDIFF,i,j,k)) * ONE_D_DZ(k);
		
		zetaF.j = (U_AT_W(UDIFF,i,j,k+1) - U_AT_W(UDIFF,i,j,k)) * ONE_D_DZ(k) 
			  -   (W_AT_U(WDIFF,i+1,j,k) - W_AT_U(WDIFF,i,j,k)) * one_d_dx;
		
		zetaF.k = (V_AT_U(VDIFF,i+1,j,k) - V_AT_U(VDIFF,i,j,k)) * one_d_dx 	 
			  -   (U_AT_V(UDIFF,i,j+1,k) - U_AT_V(UDIFF,i,j,k)) * one_d_dy;

		//----------------------------------------------------------
		// Potential vorticity
		//----------------------------------------------------------
		pv[ind] = (DOT_PRODUCT(zeta,dtheta) + DOT_PRODUCT(Zeta,dtheta) + DOT_PRODUCT(zeta,dTheta)) * one_d_rhou[k];
		
		//----------------------------------------------------------
		// Diabatic generation of potential vorticity
		//----------------------------------------------------------
		pv_diabatic[ind] += ( DOT_PRODUCT(zeta,dthetaDot) + DOT_PRODUCT(Zeta,dthetaDot) ) * one_d_rhou[k];

		//----------------------------------------------------------
		// Potential vorticity diffusion
		//----------------------------------------------------------
		pv_diff[ind] += ( DOT_PRODUCT(zetaF,dTheta) + DOT_PRODUCT(zetaF,dtheta) ) * one_d_rhou[k] * dt;
	
		//----------------------------------------------------------
		// Advective generation of potential vorticity
		//----------------------------------------------------------
		PV_advect[ind] -= (0.5 * (U(i+1,j,k)*(pv_base[ind]+pv_base[INDEX(i+1,j,k)]) - U(i,j,k)*(pv_base[ind]+pv_base[INDEX(i-1,j,k)])) * one_d_dx +
						   0.5 * (V(i,j+1,k)*(pv_base[ind]+pv_base[INDEX(i,j+1,k)]) - V(i,j,k)*(pv_base[ind]+pv_base[INDEX(i,j-1,k)])) * one_d_dy +
						   0.5 * (rhow[k+1]*W(i,j,k+1)*(pv_base[ind]+pv_base[INDEX(i,j,k+1)]) 
							   -  rhow[k  ]*W(i,j,k  )*(pv_base[ind]+pv_base[INDEX(i,j,k-1)])) * ONE_D_DZ(k) * one_d_rhou[k]
						  ) * dt;
						
	}}}
	
	exchange(pv);
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		ind = INDEX(i,j,k);
		
		//----------------------------------------------------------
		// Advection of potential vorticity by perturbation wind
		//----------------------------------------------------------
		pv_advect[ind] -= (0.5 * (U(i+1,j,k)*(pv[ind]+pv[INDEX(i+1,j,k)]) - U(i,j,k)*(pv[ind]+pv[INDEX(i-1,j,k)])) * one_d_dx +
						   0.5 * (V(i,j+1,k)*(pv[ind]+pv[INDEX(i,j+1,k)]) - V(i,j,k)*(pv[ind]+pv[INDEX(i,j-1,k)])) * one_d_dy +
						   0.5 * (rhow[k+1]*W(i,j,k+1)*(pv[ind]+pv[INDEX(i,j,k+1)]) 
							   -  rhow[k  ]*W(i,j,k  )*(pv[ind]+pv[INDEX(i,j,k-1)])) * ONE_D_DZ(k) * one_d_rhou[k]
						  ) * dt;
		
		//----------------------------------------------------------
		// Advection of potential vorticity by basic state wind
		//----------------------------------------------------------
		pv_ADVECT[ind] -= (0.5 * (UBAR(i+1,j,k)*(pv[ind]+pv[INDEX(i+1,j,k)]) - UBAR(i,j,k)*(pv[ind]+pv[INDEX(i-1,j,k)])) * one_d_dx +
						   0.5 * (VBAR(i,j+1,k)*(pv[ind]+pv[INDEX(i,j+1,k)]) - VBAR(i,j,k)*(pv[ind]+pv[INDEX(i,j-1,k)])) * one_d_dy +
						   0.5 * (rhow[k+1]*WBAR(i,j,k+1)*(pv[ind]+pv[INDEX(i,j,k+1)]) 
							   -  rhow[k  ]*WBAR(i,j,k  )*(pv[ind]+pv[INDEX(i,j,k-1)])) * ONE_D_DZ(k) * one_d_rhou[k]
						  ) * dt;
	}}}
	

}

/*********************************************************************
*
*
**********************************************************************/
void smooth(double * var,int il,int ih,int jl,int jh,int spacing){
	
	double tempvar;
	int count;
	int smooth_range = 1;
	

	for(int i=ibs[rank]%spacing+il;i<ih;i+=spacing){
	for(int j=jbs[rank]%spacing+jl;j<jh;j+=spacing){
	for(int k=1;k<NZ-1;k++){
#if 0
		if(k==2 && rank==2){
			printf("%d %d | ",ibs[rank]+i-3,jbs[rank]+j-3);
			fflush(stdout);
		}
#endif
		tempvar = 0;
		count = 0;
		
		for(int l=-smooth_range;l<=smooth_range;l++){
		for(int m=-smooth_range;m<=smooth_range;m++){
			
			 tempvar += var[INDEXT(i+l,j+m,k)];
			 count++;
		}}
		
		var[INDEXT(i,j,k)] = tempvar / (double)count;
		
	}}}
	
//	exit(0);
}

/*********************************************************************
*
*
**********************************************************************/
void write_tvar_to_file(const char *myfilename,const char *varname,double *var,bool zero_out){
	
	/*
	if(grid_ratio>1){ smooth(var,3,fNX-3,3,fNY-3,grid_ratio);}
	
	gatherArrays3(var,&u[0][0][0]);
	
	if(rank==0){
		
		if(grid_ratio>1){
		
			for(int i=0;i<subNX;i++){
			for(int j=0;j<subNY;j++){
			for(int k=0;k<NZ;k++){

				smallArray[SUBINDEX(i,j,k)] = u[i*grid_ratio][j*grid_ratio][k];

			}}}
		
			write_pvar_to_file(myfilename,varname,smallArray,subNX,subNY,file_time_counter);
			
		} else {
			
			write_pvar_to_file(myfilename,varname,&u[0][0][0],subNX,subNY,file_time_counter);	
		}
	}
	*/
	
	parallel_write_pvar_to_file(myfilename,varname,var,file_time_counter);
	
	if(zero_out){ for(int i=0;i<fNX*fNY*NZ;i++){ var[i] = 0; }}
	
}

/*********************************************************************
*
*
**********************************************************************/
void write_tvar_to_file_laplacian(const char *myfilename,const char *varname,double *var,bool zero_out){
	
	gatherArrays2(var,u);
		
	if(rank==0){
		
		double input[NX][NY];
		
		
		for(int k=0;k<NZ;k++){
		
			for(int i=0;i<NX;i++){
			for(int j=0;j<NY;j++){
		
				input[i][j] = u[i][j][k];				
			}}
		
			if(PERIODIC_BOUNDARIES){
				run_laplacian_solver(input);
			} else {
				run_laplacian_solver_real(input);
			}
		
			for(int i=0;i<NX;i++){
			for(int j=0;j<NY;j++){
		
				 u[i][j][k] = input[i][j];				
			}}
		
		}
		
		
		write_pvar_to_file(myfilename,varname,&u[0][0][0],NX,NY,file_time_counter);
	}
	
	if(zero_out){ for(int i=0;i<fNX*fNY*NZ;i++){ var[i] = 0; }}
	
}

/*********************************************************************
*
*
**********************************************************************/
void write_PVvars_to_file(const char *myfilename){
	
	write_tvar_to_file(myfilename,"pv",pv,true);
	write_tvar_to_file(myfilename,"pv_diabatic",pv_diabatic,true);
	write_tvar_to_file(myfilename,"PV_advect",PV_advect,true);
	write_tvar_to_file(myfilename,"pv_advect",pv_advect,true);
	write_tvar_to_file(myfilename,"pv_ADVECT",pv_ADVECT,true);
	write_tvar_to_file(myfilename,"pv_diff",pv_diff,true);
}

/*********************************************************************
*
*
**********************************************************************/
void write_pv_tracer_to_file(const char *myfilename){
	
	write_tvar_to_file(myfilename,"pv",pv_trace_initial3,false);
	write_tvar_to_file(myfilename,"pv_diabatic",pv_trace_diabatic3,false);
	write_tvar_to_file(myfilename,"PV_advect",pv_trace_advect3,false);
}

/*********************************************************************
*
*
**********************************************************************/
void write_pevars_to_file(const char *myfilename){
	
	write_tvar_to_file(myfilename,"pe_pp_hor_adv",pe_pp_hor_adv,true);
	write_tvar_to_file(myfilename,"pe_bp_hor_adv",pe_bp_hor_adv,true);
	write_tvar_to_file(myfilename,"pe_pb_hor_adv",pe_pb_hor_adv,true);
	write_tvar_to_file(myfilename,"pe_pp_vert_adv",pe_pp_vert_adv,true);
	write_tvar_to_file(myfilename,"pe_bp_vert_adv",pe_bp_vert_adv,true);
	write_tvar_to_file(myfilename,"pe_pb_vert_adv",pe_pb_vert_adv,true);
	write_tvar_to_file(myfilename,"pe_m_diabatic",pe_m_diabatic,true);
	write_tvar_to_file(myfilename,"pe_tdiff",pe_tdiff,true);
}

/*********************************************************************
*
*
**********************************************************************/
void write_tvars_to_file(const char *myfilename){
	
	write_tvar_to_file(myfilename,"pp_hor_adv",pp_hor_adv,true);
	write_tvar_to_file(myfilename,"bp_hor_adv",bp_hor_adv,true);
	write_tvar_to_file(myfilename,"pb_hor_adv",pb_hor_adv,true);
	write_tvar_to_file(myfilename,"pp_vert_adv",pp_vert_adv,true);
	write_tvar_to_file(myfilename,"bp_vert_adv",bp_vert_adv,true);
	write_tvar_to_file(myfilename,"pb_vert_adv",pb_vert_adv,true);
	write_tvar_to_file(myfilename,"m_diabatic",m_diabatic,true);
	write_tvar_to_file(myfilename,"tdiff",tdiff,true);
}

/*********************************************************************
*
*
**********************************************************************/
void write_qvars_to_file(const char *myfilename){
	
	write_tvar_to_file(myfilename,"qpp_hor_adv",qpp_hor_adv,true);
	write_tvar_to_file(myfilename,"qbp_hor_adv",qbp_hor_adv,true);
	write_tvar_to_file(myfilename,"qpb_hor_adv",qpb_hor_adv,true);
	write_tvar_to_file(myfilename,"qpp_vert_adv",qpp_vert_adv,true);
	write_tvar_to_file(myfilename,"qbp_vert_adv",qbp_vert_adv,true);
	write_tvar_to_file(myfilename,"qpb_vert_adv",qpb_vert_adv,true);
	write_tvar_to_file(myfilename,"qm_diabatic",qm_diabatic,true);
	write_tvar_to_file(myfilename,"qdiff",qdiff,true);
}

/*********************************************************************
*
*
**********************************************************************/
void write_vvars_to_file(const char *myfilename){
	
	write_tvar_to_file_laplacian(myfilename,"vort_hor_adv",vort_hor_adv,true);
	
	for(int i=1;i<fNX-1;i++){
	for(int j=1;j<fNY-1;j++){
	for(int k=1;k<NZ-1;k++){
		
		vort_hor_adv[INDEXT(i,j,k)] = VORT_AT_SCALAR(ZETA,i,j,k);
	}}}
	
	write_tvar_to_file_laplacian(myfilename,"vort",vort_hor_adv,true);
	
	write_tvar_to_file_laplacian(myfilename,"vort_hor_con",vort_hor_con,true);
	write_tvar_to_file_laplacian(myfilename,"cor_hor_con" ,cor_hor_con ,true);
	write_tvar_to_file_laplacian(myfilename,"vort_HOR_ADV",vort_HOR_ADV,true);
	write_tvar_to_file_laplacian(myfilename,"VORT_hor_adv",VORT_hor_adv,true);
	write_tvar_to_file_laplacian(myfilename,"VORT_hor_con",VORT_hor_con,true);
	write_tvar_to_file_laplacian(myfilename,"vort_ver_adv",vort_ver_adv,true);
	write_tvar_to_file_laplacian(myfilename,"VORT_ver_adv",VORT_ver_adv,true);
	write_tvar_to_file_laplacian(myfilename,"vort_VER_ADV",vort_VER_ADV,true);
	write_tvar_to_file_laplacian(myfilename,"cor_hor_adv",cor_hor_adv,true);
	
	for(int i=1;i<fNX-1;i++){
	for(int j=1;j<fNY-1;j++){
	for(int k=1;k<NZ-1;k++){
		
		vort_hor_adv[INDEXT(i,j,k)] = TILT_AT_SCALAR(vort_tilt,i,j,k);
		VORT_hor_adv[INDEXT(i,j,k)] = TILT_AT_SCALAR(VORT_tilt,i,j,k);
		vort_HOR_ADV[INDEXT(i,j,k)] = TILT_AT_SCALAR(vort_TILT,i,j,k);

	}}}
	
	write_tvar_to_file_laplacian(myfilename,"vort_tilt",vort_hor_adv,true);
	write_tvar_to_file_laplacian(myfilename,"VORT_tilt",VORT_hor_adv,true);
	write_tvar_to_file_laplacian(myfilename,"vort_TILT",vort_HOR_ADV,true);
	
	for(int i=1;i<fNX-1;i++){
	for(int j=1;j<fNY-1;j++){
	for(int k=1;k<NZ-1;k++){
		
		ZETA(i,j,k) = GET_VORT(UFRIC_TEND,VFRIC_TEND,i,j,k);

	}}}
	
	for(int i=1;i<fNX-1;i++){
	for(int j=1;j<fNY-1;j++){
	for(int k=1;k<NZ-1;k++){
		
		vort_hor_adv[INDEXT(i,j,k)] = VORT_AT_SCALAR(ZETA,i,j,k);

	}}}
	
	
	write_tvar_to_file_laplacian(myfilename,"vort_fric",vort_hor_adv,true);

	write_tvar_to_file_laplacian(myfilename,"divg",divg,true);
	//write_tvar_to_file_laplacian(myfilename,"cond",cond,true);
	
	for(int i=0;i<fNX*fNY*fNZ;i++){
		
		vort_tilt[i] = 0;
		VORT_tilt[i] = 0;
		vort_TILT[i] = 0;
		vort_ufric[i] = 0;
		vort_vfric[i] = 0;
	}
			
}

/*********************************************************************
*
*
**********************************************************************/
void calculate_Q1(double *u,double *v,double *t,double *q1){
	
	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=0;k<NZ;k++){
	
		q1[INDEXB(i,j,k)] = 
		- (
			   (u[INDEXB(i+1,j,k)] - u[INDEXB(i  ,j,k)]) * one_d_dx *
		0.50 * (t[INDEXB(i+1,j,k)] - t[INDEXB(i-1,j,k)]) * one_d_dx
			+
		0.25 * (
			(v[INDEXB(i+1,j+1,k)] + v[INDEXB(i+1,j,k)]) - (v[INDEXB(i-1,j+1,k)] + v[INDEXB(i-1,j,k)])
			)  * one_d_dx *
		0.5 * (
			t[INDEXB(i,j+1,k)] - t[INDEXB(i,j-1,k)]
			 ) * one_d_dy			
		)
			;
				
	}}}
	/*
	Q[j,i] = (
	(u[j,i+1] - u[j,i])/dx * 
	0.5*(t[j,i+1]-t[j,i-1])/dx 
	+
	0.5*(0.5*(v[j+1,i+1]+v[j,i+1]) - 0.5*(v[j+1,i-1]+v[j,i-1]))/dx *
	0.5*(t[j+1,i]-t[j-1,i])/dy
	)
	*/
}

/*********************************************************************
*
*
**********************************************************************/
void calculate_Q2(double *u,double *v,double *t,double *q2){
	
	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=0;k<NZ;k++){

		q2[INDEXB(i,j,k)] =
	
		- (
		0.25 * (
			(u[INDEXB(i,j+1,k)] + u[INDEXB(i+1,j+1,k)]) - (u[INDEXB(i,j-1,k)] + u[INDEXB(i+1,j-1,k)])
				) * one_d_dy
		     * (
				t[INDEXB(i+1,j,k)] - t[INDEXB(i-1,j,k)]
			   ) * one_d_dx
			+
		0.50 * (
			v[INDEXB(i,j+1,k)] - v[INDEXB(i,j,k)]
			)  * one_d_dy
		   * (
			t[INDEXB(i,j+1,k)] - t[INDEXB(i,j-1,k)]
			 ) * one_d_dy			
		)
			;
				
	}}}
	/*
	Q[j,i] = (
	0.5*(0.5*(u[j+1,i]+u[j+1,i+1]) - 0.5*(u[j-1,i]+u[j-1,i+1]))/dy *
	0.5*(t[j,i+1]-t[j,i-1])/dx
	+
	(v[j+1,i] - v[j,i])/dy *
	0.5*(t[j+1,i]-t[j-1,i])/dy
	)
	*/
}

/*********************************************************************
* Adds two arrays
*
**********************************************************************/
void add_array(double *x,double *y,int il,int ih,int jl,int jh){
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=0;k<NZ;k++){
	
		x[INDEX(i,j,k)] += y[INDEX(i,j,k)];
				
	}}}
}

/*********************************************************************
* Adds two arrays
*
**********************************************************************/
void add_array(double *x,double *y){
	
	add_array(x,y,0,NX,0,NY);

}

/*********************************************************************
*
*
**********************************************************************/
void calculate_omega(double *u,double *v,double *t,double *ubar,double *vbar,double *tbar,double *o,double *q1,double *q2,int ni,int nj,int nk){
	
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<NZ;k++){
	
		q1[INDEXB(i,j,k)] = 0;
		q2[INDEXB(i,j,k)] = 0;
		omega[INDEXB(i,j,k)] = 0;
				
	}}}
	
	
	calculate_Q1(u,v,tbar,omega);
	add_array(q1,omega);
	calculate_Q1(ubar,vbar,t,omega);
	add_array(q1,omega);
	calculate_Q1(u,v,t,omega);
	add_array(q1,omega);

	calculate_Q2(u,v,tbar,omega);
	add_array(q2,omega);
	calculate_Q2(ubar,vbar,t,omega);
	add_array(q2,omega);
	calculate_Q2(u,v,t,omega);
	add_array(q2,omega);
	
	
	for(int i=1;i<ni-1;i++){
	for(int j=1;j<nj-1;j++){
	for(int k=0;k<NZ;k++){
	
		omega[INDEXB(i,j,k)] = 2.0*(0.5*(one_d_dx*(q1[INDEXB(i+1,j,k)] - q1[INDEXB(i-1,j,k)]) + one_d_dy*(q2[INDEXB(i,j+1,k)] - q2[INDEXB(i,j-1,k)])));
		
		omega[INDEXB(i,j,k)] += DFDY(j) * 0.5* ( t[INDEXB(i+1,j,k)] - t[INDEXB(i-1,j,k)] ) * one_d_dx;
				
	}}}
}

/*********************************************************************
*
*
**********************************************************************/
void get_vorticity_at_scalars(double *us,double *vs,double *vort,int il,int ih,int jl,int jh){
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=0;k<NZ;k++){
	
		VORT[INDEXB(i,j,k)] = GET_VORT(U,V,i,j,k);
		
	}}}
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=0;k<NZ;k++){
		
		vort[INDEXB(i,j,k)] = 0.25*(VORT[INDEXB(i,j+1,k)]+VORT[INDEXB(i+1,j+1,k)]+VORT[INDEXB(i,j,k)]+VORT[INDEXB(i+1,j,k)]);

	}}}
}

/*********************************************************************
*
*
**********************************************************************/
void get_nondivergent_winds(double *us,double *vs,double *vort){

	double input[NX][NY];

	for(int k=0;k<NZ;k++){

		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){

			input[i][j] = vort[INDEXB(i,j,k)];
			
		}}

		run_laplacian_solver_real(input);

		for(int i=1;i<NX-1;i++){
		for(int j=1;j<NY-1;j++){

			vort[INDEXB(i,j,k)] = input[i][j];
			
			U(i,j,k) = -0.5*one_d_dy * (input[i][j+1] - input[i][j-1]);
			V(i,j,k) = 0.5*one_d_dx * (input[i+1][j] - input[i-1][j]);
			
		}}
	}
	
}

/*********************************************************************
*
*
**********************************************************************/
void calculate_omega_driver(const char *infilename,const char *outfilename,int *times,int count){
	
	//----------------------------------------------
	// Get grid data from input file
	//----------------------------------------------
	get_data(infilename,"lat",NY,&outLats[0]);
	get_data(infilename,"lon",NX,&outLons[0]);
	get_data(infilename,"tb",NZ,&tb[0]);
	
	//----------------------------------------------------------------
	// Coriolis parameter
	//----------------------------------------------------------------
	for(int j=0;j<NY;j++){ 
		
		f[j] = 2*fc*sin(outLats[j]*trigpi/180.);
		dfdy[j] = 2*fc*cos(outLats[j]*trigpi/180.) * (1./meters_per_degree) * (trigpi/180.) ;
	}
	
	stretched_grid(&zsu[0],&mu[0],&zsw[0],&mw[0],height_lowest_level,index_lowest_level);
	
	//----------------------------------------------
	// Create the output file
	//----------------------------------------------
	create_budget_file(outfilename,3);

	init_laplacian_solver_real(NX,NY);

	init_laplacian_solver3d(NX,NY,NZ);

	//----------------------------------------------
	// Memory allocation
	//----------------------------------------------	
	int size = NX*NY*NZ;
	
	m_ubar  = ALLOC(size);
	m_vbar  = ALLOC(size);
	m_thbar = ALLOC(size);

	us  = ALLOC(size);
	vs  = ALLOC(size);
	ths = ALLOC(size);
	
	omega = ALLOC(size);
	Q1	  = ALLOC(size);
	Q2	  = ALLOC(size);
	
	vort = ALLOC(size);
	VORT = ALLOC(size);
	
	
	//----------------------------------------------
	// Get basic state data
	//----------------------------------------------
	get_data(infilename,"ubar",size,m_ubar);
	get_data(infilename,"vbar",size,m_vbar);
	get_data(infilename,"thbar",size,m_thbar);

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
		THBAR(i,j,k) = grav*( THBAR(i,j,k)  ) / tb[k];
	}}}

	//----------------------------------------------
	// Calculate for each time step
	//----------------------------------------------	
	for(int i=0;i<count;i++){
		
		get_data_at_time(infilename,"u-wind",times[i],us);
		get_data_at_time(infilename,"v-wind",times[i],vs);
		get_data_at_time(infilename,"theta",times[i],ths);
		
		for(int a=0;a<NX;a++){
		for(int j=0;j<NY;j++){
		for(int k=0;k<NZ;k++){
			   TH(a,j,k) = grav*( TH(a,j,k) ) / tb[k];
		}}}
		
		get_vorticity_at_scalars(us,vs,vort,1,NX-1,1,NY-1);
		
		get_nondivergent_winds(us,vs,vort);
		write_pvar_to_file(outfilename,ovar_names[1],   vort,NX,NY,i);
		calculate_omega(us,vs,ths,m_ubar,m_vbar,m_thbar,omega,Q1,Q2,NX,NY,NZ);
		
		run_laplacian_solver3d(NX,NY,NZ,omega);
		
		//printf("%d %f %s\n",times[i],TH(NX/2,NY/2,10),infilename);
		
		write_pvar_to_file(outfilename,ovar_names[0],omega,NX,NY,i);

		write_pvar_to_file(outfilename,ovar_names[2],   vs,NX,NY,i);
	}
	
	//----------------------------------------------
	// Clean up
	//----------------------------------------------
	free(iubar); free(ivbar); free(ithbar);
	free(us); free(vs); free(ths);
}

/*********************************************************************
*
*
**********************************************************************/
void initialize_omega(){
	
}


/*********************************************************************
*
*
**********************************************************************/
void pv_tracer_advect(double step,int il,int ih,int jl,int jh,int stage,int last_stage){
	
	exchange(pv_trace_diabatic2);
	exchange(pv_trace_advect2  );
	exchange(pv_trace_initial2  );
		
	advect_scalar(step,il,ih,jl,jh,pv_trace_diabatic1,pv_trace_diabatic2,pv_trace_diabatic3,pvcell);
	
	//advect_scalar(step,il,ih,jl,jh,pv_trace_advect1,pv_trace_advect2,pv_trace_advect3,pv_base,pvcell,pbcell);

	advect_scalar(step,il,ih,jl,jh,pv_trace_advect1,pv_trace_advect2,pv_trace_advect3,pvcell);

	advect_scalar(step,il,ih,jl,jh,pv_trace_initial1,pv_trace_initial2,pv_trace_initial3,pvcell);
	
	//set_terrain(pv_trace_diabatic3,istopos,fNX*fNY*fNZ);
	//set_terrain(pv_trace_advect3  ,istopos,fNX*fNY*fNZ);
	//set_terrain(pv_trace_initial3 ,istopos,fNX*fNY*fNZ);
	
	if(stage < last_stage){
		
		switch_array(&pv_trace_diabatic2,&pv_trace_diabatic3);
		switch_array(&pv_trace_advect2,  &pv_trace_advect3  );
		switch_array(&pv_trace_initial2, &pv_trace_initial3  );	
	} 

}

/*********************************************************************
*
*
**********************************************************************/
void pv_tracer_sources(int il,int ih,int jl,int jh){
	
	vector zeta,Zeta,zetaF;
	vector zetaB,thetaB;
	vector dtheta,dTheta;
	vector dthetaDot;
	
	double advection_pv,dflux1,dflux2,boundary;
	
	int ind;
	
	//advect_scalar(1.0,il,ih,jl,jh,us,vs,ws,pv_base,base_advect,pvcell);

	exchange(ups);
	exchange(vps);
	exchange(wps);
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		
		vorticity_gradient_vector(ups,vps,wps,i,j,k,&zeta);
		vorticity_gradient_vector(m_ubar,m_vbar,m_wbar,i,j,k,&Zeta);
		
		avortVector[INDEX(i,j,k)].i = zeta.i + Zeta.i;
		avortVector[INDEX(i,j,k)].j = zeta.j + Zeta.j;
		avortVector[INDEX(i,j,k)].k = zeta.k + Zeta.k + FC(j);		
				
	}}}
	
	exchange(&DIABATIC(0,0,0));
	exchange(u_friction);
	exchange(v_friction);
	exchange(w_friction);
	exchange(u_bound);
	exchange(v_bound);
	exchange(w_bound);
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		ind = INDEX(i,j,k);

		//----------------------------------------------------------
		// Diabatic heating gradient vector
		//----------------------------------------------------------
		scalar_gradient_vector(  &DIABATIC(0,0,0),i,j,k,&dthetaDot);
		
		//----------------------------------------------------------
		// Potential temperature gradient vector
		//----------------------------------------------------------	
		scalar_gradient_vector(  &TH(0,0,0)      ,i,j,k,&dtheta   );

		scalar_gradient_vector(  &THBAR(0,0,0)   ,i,j,k,&dTheta   );

		dTheta.k += 0.5 * (tb[k+1] - tb[k-1] ) * ONE_D_DZ(k);

		scalar_gradient_vector(  t_bound   ,i,j,k,&thetaB   );
		
		//----------------------------------------------------------
		// Vorticity vectors
		//----------------------------------------------------------
		vorticity_gradient_vector(us,    vs,    ws,    i,j,k,&zeta);
		
		vorticity_gradient_vector(m_ubar,m_vbar,m_wbar,i,j,k,&Zeta);
		
		Zeta.k += FC(j);
		
		//----------------------------------------------------------
		//  Diffusional vorticity vector
		//----------------------------------------------------------
		vorticity_gradient_vector(u_friction,v_friction,w_friction,i,j,k,&zetaF);

		//----------------------------------------------------------
		//  Boundary velocity tendency vector
		//----------------------------------------------------------
		vorticity_gradient_vector(u_bound,v_bound,w_bound,i,j,k,&zetaB);
				
		//----------------------------------------------------------
		// Potential vorticity
		//----------------------------------------------------------
		pv[ind] = (DOT_PRODUCT(zeta,dtheta) + DOT_PRODUCT(Zeta,dtheta) + DOT_PRODUCT(zeta,dTheta)) * one_d_rhou[k];
		
		//----------------------------------------------------------
		// Diabatic generation of potential vorticity
		//----------------------------------------------------------
		//pv_trace_diabatic3[ind] += ( DOT_PRODUCT(zeta,dthetaDot) + DOT_PRODUCT(Zeta,dthetaDot) ) * one_d_rhou[k];

		//----------------------------------------------------------
		// Potential vorticity diffusion
		//----------------------------------------------------------
		pv_trace_diabatic3[ind] += ( DOT_PRODUCT(zetaF,dTheta) + DOT_PRODUCT(zetaF,dtheta) ) * one_d_rhou[k] * dt;

		//----------------------------------------------------------
		// Boundary
		//----------------------------------------------------------
		//pv_trace_diabatic3[ind] += ( DOT_PRODUCT(zetaB,dTheta) + DOT_PRODUCT(zetaB,dtheta) ) * one_d_rhou[k];
		//pv_trace_diabatic3[ind] += ( DOT_PRODUCT(zeta,thetaB) + DOT_PRODUCT(Zeta,thetaB) ) * one_d_rhou[k];
	
		//----------------------------------------------------------
		// Advective generation of potential vorticity
		//----------------------------------------------------------
//		advection_pv = -base_advect[ind];
#if 1
		advection_pv = 
						(0.5 * (U(i+1,j,k)*(pv_base[ind]+pv_base[INDEX(i+1,j,k)]) - U(i,j,k)*(pv_base[ind]+pv_base[INDEX(i-1,j,k)])) * one_d_dx +
						 0.5 * (V(i,j+1,k)*(pv_base[ind]+pv_base[INDEX(i,j+1,k)]) - V(i,j,k)*(pv_base[ind]+pv_base[INDEX(i,j-1,k)])) * one_d_dy +
							 0.5 * (rhow[k+1]*W(i,j,k+1)*(pv_base[ind]+pv_base[INDEX(i,j,k+1)]) 
								  -  rhow[k  ]*W(i,j,k  )*(pv_base[ind]+pv_base[INDEX(i,j,k-1)])) * ONE_D_DZ(k) * one_d_rhou[k]
						  ) * dt;
#endif
		
		
		if(k>1 && k<NZ-2){
		
			dflux1 = CubicInterpolate(DIABATIC(i,j,k-1),DIABATIC(i,j,k),DIABATIC(i,j,k+1),DIABATIC(i,j,k+2),
			zsu[k-1],zsu[k],zsu[k+1],zsu[k+2],(zsw[k+1]-zsu[k])/(zsu[k+1]-zsu[k]) );
		
			dflux2 = CubicInterpolate(DIABATIC(i,j,k-2),DIABATIC(i,j,k-1),DIABATIC(i,j,k),DIABATIC(i,j,k+1),
			zsu[k-2],zsu[k-1],zsu[k],zsu[k+1],(zsw[k]-zsu[k-1])/(zsu[k]-zsu[k-1]) );
		} else {
			dflux1 = 0.5*(DIABATIC(i,j,k),DIABATIC(i,j,k+1));
			dflux2 = 0.5*(DIABATIC(i,j,k),DIABATIC(i,j,k-1));
		}
		
		pv_trace_diabatic3[ind] +=
						(
						(
							   //(avortVector[INDEX(i,j,k)].i+avortVector[INDEX(i+1,j,k)].i)	*
CubicInterpolate2(avortVector[INDEX(i-1,j,k)].i,avortVector[INDEX(i,j,k)].i,avortVector[INDEX(i+1,j,k)].i,avortVector[INDEX(i+2,j,k)].i,0.5)								   
								   
								    * CubicInterpolate2(DIABATIC(i-1,j,k),DIABATIC(i,j,k),DIABATIC(i+1,j,k),DIABATIC(i+2,j,k),0.5) 
							-  
																//(avortVector[INDEX(i,j,k)].i+avortVector[INDEX(i-1,j,k)].i)	*
CubicInterpolate2(avortVector[INDEX(i-2,j,k)].i,avortVector[INDEX(i-1,j,k)].i,avortVector[INDEX(i,j,k)].i,avortVector[INDEX(i+1,j,k)].i,0.5)*																									 CubicInterpolate2(DIABATIC(i-2,j,k),DIABATIC(i-1,j,k),DIABATIC(i,j,k),DIABATIC(i+1,j,k),0.5)
						) * one_d_dx 
								   +
								
  						(
  							   //(avortVector[INDEX(i,j,k)].j+avortVector[INDEX(i,j+1,k)].j)	* 
				CubicInterpolate2(avortVector[INDEX(i,j-1,k)].j,avortVector[INDEX(i,j,k)].j,avortVector[INDEX(i,j+1,k)].j,avortVector[INDEX(i,j+2,k)].j,0.5) * 
				CubicInterpolate2(DIABATIC(i,j-1,k),DIABATIC(i,j,k),DIABATIC(i,j+1,k),DIABATIC(i,j+2,k),0.5)
													 //(DIABATIC(i,j,k)+DIABATIC(i,j+1,k)) 
  							-  //(avortVector[INDEX(i,j,k)].j+avortVector[INDEX(i,j-1,k)].j)	*
		CubicInterpolate2(avortVector[INDEX(i,j-2,k)].j,avortVector[INDEX(i,j-1,k)].j,avortVector[INDEX(i,j,k)].j,avortVector[INDEX(i,j+1,k)].j,0.5) *									  
				CubicInterpolate2(DIABATIC(i,j-2,k),DIABATIC(i,j-1,k),DIABATIC(i,j,k),DIABATIC(i,j+1,k),0.5)
																		   //(DIABATIC(i,j,k)+DIABATIC(i,j-1,k)) 
  						) * one_d_dy 
								+
						0.5 * (
							   (avortVector[INDEX(i,j,k)].k+avortVector[INDEX(i,j,k+1)].k)	* dflux1
																								
																				
																								
							-  (avortVector[INDEX(i,j,k)].k+avortVector[INDEX(i,j,k-1)].k)	* dflux2
						) * ONE_D_DZ(k)
																	
						  ) * one_d_rhou[k];
		
		
		if( !isSaturated[INDEX(i,j,k)] ){
			pv_trace_advect3[ind]   -= advection_pv;
		} else {
			pv_trace_diabatic3[ind] -= advection_pv;
		}
						
	}}}
#if 0
	for(int i=0;i<fNX*fNY*NZ;i++){ base_advect[i] = 0; }
	add_array(base_advect,pv_trace_diabatic2,il,ih,jl,jh);
	add_array(base_advect,pv_trace_advect2,il,ih,jl,jh);
	add_array(base_advect,pv_trace_initial2,il,ih,jl,jh);
	
	if(rank==2){
	
		printf("%e\n",average_abs_diff(base_advect,pv,il,ih,jl,jh) );
	}
#endif
	if(bigcounter==0){
		
		exchange(pv);
		memcpy(pv_trace_initial1,pv,fNX*fNY*fNZ*sizeof(double));
		memcpy(pv_trace_initial2,pv,fNX*fNY*fNZ*sizeof(double));
		memcpy(pv_trace_initial3,pv,fNX*fNY*fNZ*sizeof(double));
	}

	memcpy(pv_trace_diabatic1,pv_trace_diabatic3,fNX*fNY*fNZ*sizeof(double));
	memcpy(pv_trace_advect1,  pv_trace_advect3,  fNX*fNY*fNZ*sizeof(double));
	memcpy(pv_trace_initial1, pv_trace_initial3, fNX*fNY*fNZ*sizeof(double));
	
	switch_array(&pv_trace_diabatic2,&pv_trace_diabatic3);
	switch_array(&pv_trace_advect2,  &pv_trace_advect3  );
	switch_array(&pv_trace_initial2, &pv_trace_initial3 );

}


/*********************************************************************
*
**********************************************************************/
double average_abs_diff(double *var1,double *var2,int il,int ih,int jl,int jh){
	
	int counter = 0;
	double sum = 0;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		
		sum += fabs(var1[INDEX(i,j,k)] - var2[INDEX(i,j,k)]);
		
		counter++;
		
	}}}
	
	return sum / (double)counter;
}

/*********************************************************************
*
*
* @param type
*				0 - heat budget
*				1 - vorticity budget
*				2 - pv budget
**********************************************************************/
void create_budget_file(const char *myfilename,int type){
	
	int status;
	int ncid;

	int xID, yID, zID, tID;
	double ndx,ndy;
	
	ndx = dx * grid_ratio;
	ndy = dy * grid_ratio;
	
	int var_dimids[4];
	int var_id;

	/*************************************************
	* Create file
	**************************************************/
	status = nc_create(myfilename, NC_CLOBBER, &ncid);
	if (status != NC_NOERR) handle_error(status);

	//printf("nc id = %d\n",ncid);

	/*************************************************
	* Add global attributes
	**************************************************/
     status = nc_put_att_double (ncid, NC_GLOBAL, "dx", NC_FLOAT, 1, &ndx);
     if (status != NC_NOERR) handle_error(status);

     status = nc_put_att_double (ncid, NC_GLOBAL, "dy", NC_FLOAT, 1, &ndy);
     if (status != NC_NOERR) handle_error(status);

     status = nc_put_att_double (ncid, NC_GLOBAL, "dz", NC_FLOAT, 1, &dz);
     if (status != NC_NOERR) handle_error(status);

     status = nc_put_att_text (ncid, NC_GLOBAL, "processed_data_file",strlen(filename), filename);
     if (status != NC_NOERR) handle_error(status);

	/*************************************************
	* Create dimensions
	**************************************************/
	status = nc_def_dim(ncid, "x", (size_t)(subNX), &xID);
	if (status != NC_NOERR) handle_error(status);

	status = nc_def_dim(ncid, "y", (size_t)(subNY), &yID);
	if (status != NC_NOERR) handle_error(status);

	status = nc_def_dim(ncid, "z", (size_t)NZ, &zID);
	if (status != NC_NOERR) handle_error(status);

	status = nc_def_dim(ncid, "t", NC_UNLIMITED, &tID);
	if (status != NC_NOERR) handle_error(status);

	/*************************************************
	* Get dimension IDs
	**************************************************/
	var_dimids[0] = tID; var_dimids[1] = xID; var_dimids[2] = yID; var_dimids[3] = zID;
	
	/*************************************************
	* Define coordinate arrays
	**************************************************/
	var_dimids[0] = yID;

	status = nc_def_var (ncid,coordvar_names[0], NC_FLOAT, 1, var_dimids, &var_id);
	if (status != NC_NOERR) handle_error(status);

	var_dimids[0] = xID;

	status = nc_def_var (ncid,coordvar_names[1], NC_FLOAT, 1, var_dimids, &var_id);
	if (status != NC_NOERR) handle_error(status);

	var_dimids[0] = tID;

	status = nc_def_var (ncid,"time", NC_FLOAT, 1, var_dimids, &var_id);
	if (status != NC_NOERR) handle_error(status);

	/*************************************************
	* Write perturbation variables
	**************************************************/
	switch(type){

		/*************************************************
		* Heat budget
		**************************************************/
		case 0:
			for(int i=0;i<tvar_count;i++){

				status = nc_def_var (ncid,tvar_names[i], NC_FLOAT, 4, var_dimids, &var_id);
				if (status != NC_NOERR) handle_error(status);
			}	
		break;
		
		/*************************************************
		* Vorticity budget
		**************************************************/
		case 1:
			for(int i=0;i<vvar_count;i++){

				status = nc_def_var (ncid,vvar_names[i], NC_FLOAT, 4, var_dimids, &var_id);
				if (status != NC_NOERR) handle_error(status);
			}	
		break;
		
		/*************************************************
		* PV budget
		**************************************************/
		case 2:
			for(int i=0;i<pvvar_count;i++){

				status = nc_def_var (ncid,pvvar_names[i], NC_FLOAT, 4, var_dimids, &var_id);
				if (status != NC_NOERR) handle_error(status);
			}	
		break;
		
		/*************************************************
		* Omega equation
		**************************************************/
		case 3:
			for(int i=0;i<ovar_count;i++){

				status = nc_def_var (ncid,ovar_names[i], NC_FLOAT, 4, var_dimids, &var_id);
				if (status != NC_NOERR) handle_error(status);
			}	
		break;
		
		/*************************************************
		* PV tracer advection
		**************************************************/
		case 4:
			for(int i=0;i<pv_tracer_var_count;i++){

				status = nc_def_var (ncid,pv_tracer_names[i], NC_FLOAT, 4, var_dimids, &var_id);
				if (status != NC_NOERR) handle_error(status);
			}	
		break;
		
		/*************************************************
		* Moisture budget
		**************************************************/
		case 5:
			for(int i=0;i<qvar_count;i++){

				status = nc_def_var (ncid,qvar_names[i], NC_FLOAT, 4, var_dimids, &var_id);
				if (status != NC_NOERR) handle_error(status);
			}	
		break;
		/*************************************************
		* Potential energy budget
		**************************************************/
		case 6:
			for(int i=0;i<pe_var_count;i++){

				status = nc_def_var (ncid,pe_var_names[i], NC_FLOAT, 4, var_dimids, &var_id);
				if (status != NC_NOERR) handle_error(status);
			}	
		break;
	}

	/*************************************************
	* Write basic state variables
	**************************************************/
	if(type==10){

		var_dimids[0] = xID; var_dimids[1] = yID; var_dimids[2] = zID;

		status = nc_def_var (ncid,"VORT", NC_FLOAT, 3, var_dimids, &var_id);
		if (status != NC_NOERR) handle_error(status);
	}
	

	/*************************************************
	* Close file
	**************************************************/
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);
	

	double * sublon = (double*) calloc(subNX,sizeof(double));	
	double * sublat = (double*) calloc(subNY,sizeof(double));
	
	for(int i=0;i<subNX;i++){ sublon[i] = outLons[i*grid_ratio];}

	for(int i=0;i<subNY;i++){ sublat[i] = outLats[i*grid_ratio];}
	
	/********************************************************
	* Open file, get variable ID
	*********************************************************/
	status = nc_open(myfilename, NC_WRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid (ncid,coordvar_names[0], &var_id);
	if (status != NC_NOERR) handle_error(status);

	status = nc_put_var_double(ncid, var_id, sublat);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid (ncid,coordvar_names[1], &var_id);
	if (status != NC_NOERR) handle_error(status);

	status = nc_put_var_double(ncid, var_id, sublon);
	if (status != NC_NOERR) handle_error(status);

	/*************************************************
	* Close file
	**************************************************/
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);
	
	free(sublon); free(sublat);	
}