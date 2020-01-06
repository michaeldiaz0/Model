#include "stdafx.h"
#include "temperature.h"
#include "surface.h"

#define INDEXT(i,j,k) ((i)*fNYfNZ+(j)*fNZ+(k))
#define SUBINDEX(i,j,k) ((i)*subNYNZ+(j)*NZ+(k))

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

#define V_AT_W(v,i,j,k) 0.25*( v(i,j,k)+v(i,j+1,k)+v(i,j,k-1)+v(i,j+1,k-1) )
#define W_AT_V(v,i,j,k) 0.25*( v(i,j,k)+v(i,j,k+1)+v(i,j-1,k)+v(i,j-1,k+1) )

#define U_AT_W(v,i,j,k) 0.25*( v(i,j,k)+v(i+1,j,k)+v(i,j,k-1)+v(i+1,j,k-1) )
#define W_AT_U(v,i,j,k) 0.25*( v(i,j,k)+v(i,j,k+1)+v(i-1,j,k)+v(i-1,j,k+1) )

#define V_AT_U(v,i,j,k) (0.25 * (v(i,j,k)+v(i,j+1,k)+v(i-1,j,k)+v(i-1,j+1,k)))
#define U_AT_V(v,i,j,k) (0.25 * (v(i,j,k)+v(i+1,j,k)+v(i,j-1,k)+v(i+1,j-1,k)))

#define DIABATIC(i,j,k) ( cond[INDEX(i,j,k)] + evap[INDEX(i,j,k)] + t_diffusion[INDEX(i,j,k)]*dt )

#define HCONV(u,v,i,j,k) ((u(i+1,j,k)-u(i,j,k))*one_d_dx + (v(i,j+1,k)-v(i,j,k))*one_d_dy)

#define UDIFF(i,j,k) u_friction[INDEX(i,j,k)]
#define VDIFF(i,j,k) v_friction[INDEX(i,j,k)]
#define WDIFF(i,j,k) w_friction[INDEX(i,j,k)]

int grid_ratio = 1;

struct vector {

	double i;
	double j;
	double k;

};

//---------------------------------------------------
// Heat budget
//---------------------------------------------------
//const char tvar_names[11][20] = {"pp_hor_adv","pp_vert_adv","pb_hor_adv","pb_vert_adv","bp_hor_adv","bp_vert_adv","pp_conv","cond","evap"};
//const int tvar_count = 9;

const char tvar_names[2][20] = {"cond","evap"};
const int tvar_count = 2;

double *pp_hor_adv,*pp_vert_adv,*pb_hor_adv,*pb_vert_adv,*bp_hor_adv,*bp_vert_adv,*pp_conv,*cond,*evap;

//---------------------------------------------------
// Vorticity budget
//---------------------------------------------------
const char vvar_names[18][20] = {"VORT","vort",
					"vort_hor_adv","vort_HOR_ADV","VORT_hor_adv",
					"vort_ver_adv","VORT_ver_adv","vort_VER_ADV",
					"vort_hor_con","cor_hor_con" ,"VORT_hor_con",
					"vort_tilt","VORT_tilt" ,"vort_TILT",
					"cor_hor_adv","vort_fric"};
					
const int vvar_count = 16;

double *vort,*VORT,
		*vort_hor_adv,*vort_HOR_ADV,*VORT_hor_adv,
		*vort_ver_adv,*VORT_ver_adv,*vort_VER_ADV,
		*vort_hor_con,*cor_hor_con,*VORT_hor_con,
		*vort_tilt,*VORT_tilt,*vort_TILT,
		*cor_hor_adv,*vort_ufric,*vort_vfric;

//---------------------------------------------------
// PV budget
//---------------------------------------------------
const int pvvar_count = 4;

double *pv,*pv_diabatic,*pv_advect,*PV_advect;
double *pv_base;

const char pvvar_names[pvvar_count][20] = {"pv","pv_diabatic","pv_advect","PV_advect"};
double *pvvars[pvvar_count] = 			  { pv , pv_diabatic , pv_advect , PV_advect };


//---------------------------------------------------
// Shared
//---------------------------------------------------
double *smallArray;
bool smallArrayIsInitialized = false;

int subNX = NX/grid_ratio;
int subNY = NY/grid_ratio;
int subNYNZ = subNY*NZ;

void write_tvar_to_file(const char *,const char *,double *,bool);
void initialize_basic_state_PV(int il,int ih,int jl,int jh);

/*********************************************************************
*
*
**********************************************************************/
void initialize_pv(){
	
	int size = fNX*fNY*fNZ;
	
	//for(int i=0;i<pvvar_count;i++){ pvvars[i] = ALLOC(size);}
	
	pv = ALLOC(size);
	pv_diabatic = ALLOC(size);
	pv_advect = ALLOC(size);
	PV_advect = ALLOC(size);	
	pv_base = ALLOC(size);
	cond = ALLOC(size);
	evap = ALLOC(size);
	
	initialize_basic_state_PV(1,fNX-1,1,fNY-1);
	
}

/*********************************************************************
*
*
**********************************************************************/
void initialize_temperature(){
	
	int size = fNX*fNY*fNZ;
	
	pp_hor_adv = ALLOC(size);
	pp_vert_adv = ALLOC(size);
	pb_hor_adv = ALLOC(size);
	pb_vert_adv = ALLOC(size);
	bp_hor_adv = ALLOC(size);
	bp_vert_adv = ALLOC(size);
	pp_conv = ALLOC(size);
	cond = ALLOC(size);
	evap = ALLOC(size);
	
	if(!smallArrayIsInitialized){
	
		smallArray = (double*) calloc((NX/grid_ratio)*(NY/grid_ratio)*NZ,sizeof(double));
		
		smallArrayIsInitialized = true;
	}
	
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
	
	if(!smallArrayIsInitialized){
		
		smallArray = (double*) calloc((NX/grid_ratio)*(NY/grid_ratio)*NZ,sizeof(double));
		
		smallArrayIsInitialized = true;
	}
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
void calculate_heat_budget(int il,int ih,int jl,int jh){
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
	
		pp_hor_adv[INDEXT(i,j,k)] += - ((U(i+1,j,k)*(TH(i+1,j,k)+TH(i,j,k)) - U(i,j,k)*(TH(i,j,k)+TH(i-1,j,k)))
									 +  (V(i,j+1,k)*(TH(i,j+1,k)+TH(i,j,k)) - V(i,j,k)*(TH(i,j,k)+TH(i,j-1,k))))*0.5*dtx
									 + TH(i,j,k)*HCONV(U,V,i,j,k)
										 ;
		
		bp_hor_adv[INDEXT(i,j,k)] += - ((UBAR(i+1,j,k)*(TH(i+1,j,k)+TH(i,j,k)) - UBAR(i,j,k)*(TH(i,j,k)+TH(i-1,j,k)))
									 +  (VBAR(i,j+1,k)*(TH(i,j+1,k)+TH(i,j,k)) - VBAR(i,j,k)*(TH(i,j,k)+TH(i,j-1,k))))*0.5*dtx
									 + TH(i,j,k)*HCONV(UBAR,VBAR,i,j,k);
		
		bp_hor_adv[INDEXT(i,j,k)] += -0.5*( rhow[k+1]*WBAR(i,j,k+1)*(TH(i,j,k+1)+TH(i,j,k)) - rhow[k]*WBAR(i,j,k)*(TH(i,j,k)+TH(i,j,k-1)) ) * one_d_rhou[k] * DTZ(k);
		
		pb_hor_adv[INDEXT(i,j,k)] += - ((U(i+1,j,k)*(THBAR(i+1,j,k)+THBAR(i,j,k)) - U(i,j,k)*(THBAR(i,j,k)+THBAR(i-1,j,k)))
									 +  (V(i,j+1,k)*(THBAR(i,j+1,k)+THBAR(i,j,k)) - V(i,j,k)*(THBAR(i,j,k)+THBAR(i,j-1,k))))*0.5*dtx
									 + THBAR(i,j,k)*HCONV(U,V,i,j,k);
		
		pb_hor_adv[INDEXT(i,j,k)] += -0.5*( rhow[k+1]*W(i,j,k+1)*(THBAR(i,j,k+1)+THBAR(i,j,k)) - rhow[k]*W(i,j,k)*(THBAR(i,j,k)+THBAR(i,j,k-1)) ) * one_d_rhou[k] * DTZ(k);
		
		pb_vert_adv[INDEXT(i,j,k)] += -0.5*( rhow[k+1]*W(i,j,k+1)*(tb[k+1]-tb[k  ]) + rhow[k]*W(i,j,k)*(tb[k]-tb[k-1])) * one_d_rhou[k] * DTZ(k);
		
		pp_vert_adv[INDEXT(i,j,k)] += -0.5*( rhow[k+1]*W(i,j,k+1)*(TH(i,j,k+1)+TH(i,j,k)) - rhow[k]*W(i,j,k)*(TH(i,j,k)+TH(i,j,k-1)) ) * one_d_rhou[k] * DTZ(k);
		
		pp_conv[INDEXT(i,j,k)] += -HCONV(U,V,i,j,k) * TH(i,j,k);
	
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
	}}}
	
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
		dthetaDot.i = 0.5 * ( DIABATIC(i+1,j,k) - DIABATIC(i-1,j,k)) * one_d_dx;
		dthetaDot.j = 0.5 * ( DIABATIC(i,j+1,k) - DIABATIC(i,j-1,k)) * one_d_dy;
		dthetaDot.k = 0.5 * ( DIABATIC(i,j,k+1) - DIABATIC(i,j,k-1)) * ONE_D_DZ(k);

		//----------------------------------------------------------
		// Potential temperature gradient vector
		//----------------------------------------------------------		
		dtheta.i = 0.5 * (TH(i+1,j,k) - TH(i-1,j,k)) * one_d_dx;
		dtheta.j = 0.5 * (TH(i,j+1,k) - TH(i,j-1,k)) * one_d_dy;
		dtheta.k = 0.5 * (TH(i,j,k+1) - TH(i,j,k-1)) * ONE_D_DZ(k);

		dTheta.i = 0.5 * (THBAR(i+1,j,k) - THBAR(i-1,j,k)) * one_d_dx;
		dTheta.j = 0.5 * (THBAR(i,j+1,k) - THBAR(i,j-1,k)) * one_d_dy;				   
		dTheta.k = 0.5 * (THBAR(i,j,k+1) - THBAR(i,j,k-1) + tb[k+1] - tb[k-1] ) * ONE_D_DZ(k);

		//----------------------------------------------------------
		// Absolute vorticity vector
		//----------------------------------------------------------
		zeta.i = (W_AT_V(W,i,j+1,k) - W_AT_V(W,i,j,k)) * one_d_dy 	 - 	(V_AT_W(V,i,j,k+1) - V_AT_W(V,i,j,k)) * ONE_D_DZ(k);
		zeta.j = (U_AT_W(U,i,j,k+1) - U_AT_W(U,i,j,k)) * ONE_D_DZ(k) -  (W_AT_U(W,i+1,j,k) - W_AT_U(W,i,j,k)) * one_d_dx;
		zeta.k = (V_AT_U(V,i+1,j,k) - V_AT_U(V,i,j,k)) * one_d_dx 	 - 	(U_AT_V(U,i,j+1,k) - U_AT_V(U,i,j,k)) * one_d_dy;

		Zeta.i = (W_AT_V(WBAR,i,j+1,k) - W_AT_V(WBAR,i,j,k)) * one_d_dy 	- (V_AT_W(VBAR,i,j,k+1) - V_AT_W(VBAR,i,j,k)) * ONE_D_DZ(k);
		Zeta.j = (U_AT_W(UBAR,i,j,k+1) - U_AT_W(UBAR,i,j,k)) * ONE_D_DZ(k)  - (W_AT_U(WBAR,i+1,j,k) - W_AT_U(WBAR,i,j,k)) * one_d_dx;
		Zeta.k = (V_AT_U(VBAR,i+1,j,k) - V_AT_U(VBAR,i,j,k)) * one_d_dx 	- (U_AT_V(UBAR,i,j+1,k) - U_AT_V(UBAR,i,j,k)) * one_d_dy + FC(j);

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
		pv_diabatic[ind] += ( DOT_PRODUCT(zetaF,dTheta)*dt + DOT_PRODUCT(zetaF,dtheta)*dt ) * one_d_rhou[k];
		
		//----------------------------------------------------------
		// Advective generation of potential vorticity
		//----------------------------------------------------------
		pv_advect[ind] -= (0.5 * (U(i+1,j,k)*(pv_base[ind]+pv_base[INDEX(i+1,j,k)]) - U(i,j,k)*(pv_base[ind]+pv_base[INDEX(i-1,j,k)])) * one_d_dx +
						   0.5 * (V(i,j+1,k)*(pv_base[ind]+pv_base[INDEX(i,j+1,k)]) - V(i,j,k)*(pv_base[ind]+pv_base[INDEX(i,j-1,k)])) * one_d_dy +
						   0.5 * (rhow[k+1]*W(i,j,k+1)*(pv_base[ind]+pv_base[INDEX(i,j,k+1)]) 
							   -  rhow[k  ]*W(i,j,k  )*(pv_base[ind]+pv_base[INDEX(i,j,k-1)])) * ONE_D_DZ(k) * one_d_rhou[k]
						  );
			
			
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
	
	if(grid_ratio>1){ smooth(var,3,fNX-3,3,fNY-3,grid_ratio);}
	
	gatherArrays2(var,u);
	
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
	
	//for(int i=0;i<pvvar_count;i++){
		//write_tvar_to_file(myfilename,pvvar_names[i],pvvars[i],true);
	//}
}

/*********************************************************************
*
*
**********************************************************************/
void write_tvars_to_file(const char *myfilename){
	/*
	write_tvar_to_file(myfilename,"pp_hor_adv",pp_hor_adv,true);
	write_tvar_to_file(myfilename,"bp_hor_adv",bp_hor_adv,true);
	write_tvar_to_file(myfilename,"pb_hor_adv",pb_hor_adv,true);
	write_tvar_to_file(myfilename,"pb_vert_adv",pb_vert_adv,true);
	write_tvar_to_file(myfilename,"pp_vert_adv",pp_vert_adv,true);
	write_tvar_to_file(myfilename,"pp_conv",pp_vert_adv,true);
	*/
	write_tvar_to_file(myfilename,"cond",cond,true);
	write_tvar_to_file(myfilename,"evap",evap,true);
}

/*********************************************************************
*
*
**********************************************************************/
void write_vvars_to_file(const char *myfilename){
	
	write_tvar_to_file(myfilename,"vort_hor_adv",vort_hor_adv,true);
	
	for(int i=1;i<fNX-1;i++){
	for(int j=1;j<fNY-1;j++){
	for(int k=1;k<NZ-1;k++){
		
		vort_hor_adv[INDEXT(i,j,k)] = VORT_AT_SCALAR(ZETA,i,j,k);
	}}}
	
	write_tvar_to_file(myfilename,"vort",vort_hor_adv,true);
	
	write_tvar_to_file(myfilename,"vort_hor_con",vort_hor_con,true);
	write_tvar_to_file(myfilename,"cor_hor_con" ,cor_hor_con ,true);
	write_tvar_to_file(myfilename,"vort_HOR_ADV",vort_HOR_ADV,true);
	write_tvar_to_file(myfilename,"VORT_hor_adv",VORT_hor_adv,true);
	write_tvar_to_file(myfilename,"VORT_hor_con",VORT_hor_con,true);
	write_tvar_to_file(myfilename,"vort_ver_adv",vort_ver_adv,true);
	write_tvar_to_file(myfilename,"VORT_ver_adv",VORT_ver_adv,true);
	write_tvar_to_file(myfilename,"vort_VER_ADV",vort_VER_ADV,true);
	write_tvar_to_file(myfilename,"cor_hor_adv",cor_hor_adv,true);
	
	for(int i=1;i<fNX-1;i++){
	for(int j=1;j<fNY-1;j++){
	for(int k=1;k<NZ-1;k++){
		
		vort_hor_adv[INDEXT(i,j,k)] = TILT_AT_SCALAR(vort_tilt,i,j,k);
		VORT_hor_adv[INDEXT(i,j,k)] = TILT_AT_SCALAR(VORT_tilt,i,j,k);
		vort_HOR_ADV[INDEXT(i,j,k)] = TILT_AT_SCALAR(vort_TILT,i,j,k);

	}}}
	
	write_tvar_to_file(myfilename,"vort_tilt",vort_hor_adv,true);
	write_tvar_to_file(myfilename,"VORT_tilt",VORT_hor_adv,true);
	write_tvar_to_file(myfilename,"vort_TILT",vort_HOR_ADV,true);
	
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
	
	
	write_tvar_to_file(myfilename,"vort_fric",vort_hor_adv,true);
	
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

	printf("nc id = %d\n",ncid);

	/*************************************************
	* Add global attributes
	**************************************************/
     status = nc_put_att_double (ncid, NC_GLOBAL, "dx", NC_FLOAT, 1, &ndx);
     if (status != NC_NOERR) handle_error(status);

     status = nc_put_att_double (ncid, NC_GLOBAL, "dy", NC_FLOAT, 1, &ndy);
     if (status != NC_NOERR) handle_error(status);

     status = nc_put_att_double (ncid, NC_GLOBAL, "dz", NC_FLOAT, 1, &dz);
     if (status != NC_NOERR) handle_error(status);

     status = nc_put_att_text (ncid, NC_GLOBAL, "input_data",strlen(datafile), datafile);
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
	* Define forecasted variables
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

#if 0
		vort_hor_adv[INDEXT(i,j,k)] = 0.125* (vort_tilt[INDEXT(i,j+1,k  )]+vort_tilt[INDEXT(i+1,j+1,k  )]+
											  vort_tilt[INDEXT(i,j  ,k  )]+vort_tilt[INDEXT(i+1,j  ,k  )]+
											  vort_tilt[INDEXT(i,j+1,k+1)]+vort_tilt[INDEXT(i+1,j+1,k+1)]+
											  vort_tilt[INDEXT(i,j  ,k+1)]+vort_tilt[INDEXT(i+1,j  ,k+1)]);
		
		VORT_hor_adv[INDEXT(i,j,k)] = 0.125* (VORT_tilt[INDEXT(i,j+1,k  )]+VORT_tilt[INDEXT(i+1,j+1,k  )]+
											  VORT_tilt[INDEXT(i,j  ,k  )]+VORT_tilt[INDEXT(i+1,j  ,k  )]+
											  VORT_tilt[INDEXT(i,j+1,k+1)]+VORT_tilt[INDEXT(i+1,j+1,k+1)]+
											  VORT_tilt[INDEXT(i,j  ,k+1)]+VORT_tilt[INDEXT(i+1,j  ,k+1)]);

	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		
		DIFFUSE_VARIABLE(vort_HOR_ADV);
		DIFFUSE_VARIABLE(VORT_hor_adv);
		DIFFUSE_VARIABLE(vort_ver_adv);
		DIFFUSE_VARIABLE(VORT_ver_adv);
		DIFFUSE_VARIABLE(vort_VER_ADV);
		DIFFUSE_VARIABLE(vort_hor_con);
		DIFFUSE_VARIABLE(cor_hor_con);
		DIFFUSE_VARIABLE(VORT_hor_con);
		DIFFUSE_VARIABLE(vort_tilt);
		DIFFUSE_VARIABLE(VORT_tilt);
		DIFFUSE_VARIABLE(vort_TILT);
		DIFFUSE_VARIABLE(cor_hor_adv);
		DIFFUSE_VARIABLE(vort_ufric);
		DIFFUSE_VARIABLE(vort_vfric);
		DIFFUSE_VARIABLE(vort_hor_adv);
	}}}

const double dx2 = 1.0/(dx*dx);
const double dy2 = 1.0/(dy*dy);
const double dz2 = 1.0/(dz*dz);

#define DIFFUSE_VARIABLE(v)	v[INDEXT(i,j,k)] = v[INDEXT(i,j,k)] + dt* ( 	 \
			+ 0.003 * (v[INDEXT(i+1,j,k)] - 2.0*v[INDEXT(i,j,k)]+v[INDEXT(i-1,j,k)])*dx2 \
			+ 0.003 * (v[INDEXT(i,j+1,k)] - 2.0*v[INDEXT(i,j,k)]+v[INDEXT(i,j-1,k)])*dy2 \
			+ 0.0003 * mu[k]*(mw[k+1]*v[INDEXT(i,j,k+1)] - mw[k+1]*v[INDEXT(i,j,k)] - mw[k]*v[INDEXT(i,j,k)] + mw[k]*v[INDEXT(i,j,k-1)])*dz2 \
			)	


-0.25* (
											(
											 (VBAR(i-1,j,k)+VBAR(i,j,k)) - (VBAR(i-1,j,k-1)+VBAR(i,j,k-1)) 	// dv/dz
											) * DTZ(k)
											*
						 			        (
											 (W(i,j-1,k)+W(i,j,k)) - (W(i-1,j-1,k)+W(i-1,j,k)) 	// dw/dx
											) * one_d_dx
										-
			   								(
			   								 (UBAR(i,j-1,k)+UBAR(i,j,k)) - (UBAR(i,j-1,k-1)+UBAR(i,j,k-1)) 	// du/dz
			   								) * DTZ(k)
											*
						 			        (
											 (W(i-1,j,k)+W(i,j,k)) - (W(i-1,j-1,k)+W(i,j-1,k)) 	// dw/dy
											) * one_d_dy
										   );
#endif