#include "stdafx.h"
#include "turbulence.h"
#include "interpolate.h"
#include "budgets.h"
#include "pcomm.h"

#if PARALLEL
	#define KMIX(i,j,k)   Kmix[(i)*fNYfNZ+(j)*fNZ+(k)]
	#define KSMIX(i,j,k)  one_d_Pr*Kmix[(i)*fNYfNZ+(j)*fNZ+(k)]
	#define KHMIX(i,j,k)  KHmix[(i)*fNYfNZ+(j)*fNZ+(k)]
	#define KHSMIX(i,j,k) one_d_Pr*KHmix[(i)*fNYfNZ+(j)*fNZ+(k)]
#else
	#define KMIX(i,j,k)   Kmix[(i)*NY*NZ+(j)*NZ+(k)]
	#define KSMIX(i,j,k)  one_d_Pr*Kmix[(i)*NY*NZ+(j)*NZ+(k)]
	#define KHMIX(i,j,k)  KHmix[(i)*NY*NZ+(j)*NZ+(k)]
	#define KHSMIX(i,j,k) one_d_Pr*KHmix[(i)*NY*NZ+(j)*NZ+(k)]
#endif

#define STRESS(j,k) stress[(j)*NZ+(k)]

/********************************************************
*
*********************************************************/
struct stress_tensor {

	double tau_11,tau_22,tau_33;
	double tau_12,tau_13,tau_23;
	double tau_11_west,tau_12_west,tau_13_west;
	double tau_u_surface,tau_v_surface;
};

double *landsea;

double *u_friction,*v_friction,*w_friction,*t_diffusion;
double *qv_diffusion,*qc_diffusion,*qr_diffusion;
double *Kmix,*KHmix;
double *integrated_friction_ke,*latent_heat_flux;

bool *isSaturated;

double *thetaE;
double *qsat_test;
double *vert_mixing_length_squared;
double *vert_mix_Kmax;

double es_water;

const double cmixv_max = 0.15;	// upper bound on non-dimensional horizontal mixing coefficients
const double cmixh_max = 0.05;	// lower bound on non-dimensional horizontal mixxing coefficients

double kmixh_max;		// upper bound on horizontal mixing coefficients
double kmixv_max;		// lower bound on horizontal mixxing coefficients

double lH;			// horizontal mixing length

const double one_d_Pr = 3.0;		// inverse Prantl number
const double lv_max = 200.0;		// upper bound for vertical mixing length
const double cs = 0.25;
double water_temp;

struct stress_tensor *stress;

/********************************************************
* FUNCTION PROTOTYPES
*********************************************************/
void get_Kmix(int,int,int,int);
void get_Kmix_vertical(int,int,int,int);
void get_Kmix_dry(int il,int ih,int jl,int jh);
void turbulent_diffusion_velocity(int,int,int,int);
void turbulent_stress_vertical(int,int,int,int);
void turbulent_diffusion_scalars(int,int,int,int);
void turbulent_diffusion_theta(int il,int ih,int jl,int jh);
void calculate_diff_tend_full(int il,int ih,int jl,int jh);
void calculate_diff_tend_full_dry(int il,int ih,int jl,int jh);
void calculate_diff_tend_vertical(int il,int ih,int jl,int jh);
void test();

/********************************************************
*
* INITIALIZER FOR TURBULENCE AND SURFACE PARAMETERIZATION
*
* nx,ny,nz - number of grid points in each direction
* dz_wlev - height of u-levels (scalars levels)
*********************************************************/
void init_kmix(int nx,int ny,int nz,double *zlevs){

	int size = nx*ny*nz;

	kmixh_max = cmixh_max*dx*dx/dt;
	kmixv_max = cmixv_max*dz*dz/dt;
	lH = dx/4.0;

	Kmix  = (double*)calloc(size,sizeof(double));
	KHmix  = (double*)calloc(size,sizeof(double));

	thetaE = (double*) calloc(NZ,sizeof(double));
	qsat_test = (double*) calloc(NZ,sizeof(double));
	vert_mixing_length_squared = (double*) calloc(NZ,sizeof(double));
	vert_mix_Kmax = (double*) calloc(NZ,sizeof(double));

	//------------------------------------------------------
	// Storage for diffusion tendencies
	//------------------------------------------------------
	u_friction = (double*)calloc(size,sizeof(double));
	v_friction = (double*)calloc(size,sizeof(double));
	w_friction = (double*)calloc(size,sizeof(double));
	t_diffusion = (double*)calloc(size,sizeof(double));
	
	integrated_friction_ke = (double*)calloc(nx*ny,sizeof(double));
	latent_heat_flux = (double*)calloc(nx*ny,sizeof(double));
		
	if(USE_MICROPHYSICS && USE_TURBULENT_STRESS){
		
		qv_diffusion = (double*)calloc(size,sizeof(double));
		qc_diffusion = (double*)calloc(size,sizeof(double));
		qr_diffusion = (double*)calloc(size,sizeof(double));
	}
	//------------------------------------------------------
	// Storage for output from microphysics scheme's 
	// saturation adjustment
	//------------------------------------------------------
	isSaturated = (bool*)malloc(size*sizeof(bool));

	for(int i=0;i<size;i++){ isSaturated[i] = false;}

	//------------------------------------------------------
	// Rate of deformation tensor
	//------------------------------------------------------
	stress = (stress_tensor*)calloc(ny*nz,sizeof(stress_tensor));
	
	double deltaZ;

	//------------------------------------------------------
	// Set vertical mixing length values
	//------------------------------------------------------
	for(int k=0;k<nz;k++){
		
		deltaZ = (zlevs[k]-zlevs[k-1]);
		
		if(deltaZ > lv_max){ vert_mixing_length_squared[k] = lv_max*lv_max;  } 
		else {		 		 vert_mixing_length_squared[k] = deltaZ*deltaZ;}
		
		vert_mix_Kmax[k] = cmixv_max*vert_mixing_length_squared[k]/dt;
	}
	//------------------------------------------------------
	// Saturation vapor pressure for temperature of water surface
	//------------------------------------------------------
	water_temp = WATER_TEMP_C + 273.15;
	
	es_water = 611.2 * exp(17.67 * (water_temp-273.15) / (water_temp - 29.65) );
}

/********************************************************
*
* 
*
*********************************************************/
void calculate_diff_tend(int il,int ih,int jl,int jh){
	
	if(PARALLEL){
		
		exchangeOnePoint(ums);
		exchangeOnePoint(vms);
		exchangeOnePoint(wms);
		exchangeOnePoint(thms);
	}
	
	//--------------------------------------------------
	//--------------  WITH MOISTURE - ------------------
	//--------------------------------------------------
	if(USE_MICROPHYSICS){

		if(PARALLEL){
			
			exchangeOnePoint(qvms);
			exchangeOnePoint(qcms);
			exchangeOnePoint(qrms);
		}

		calculate_diff_tend_full(il,ih,jl,jh);
		
	//--------------------------------------------------
	//-------------  WITHOUT MOISTURE ------------------
	//--------------------------------------------------
	} else { 
		calculate_diff_tend_full_dry(il,ih,jl,jh);		
		//calculate_diff_tend_vertical(il,ih,jl,jh);
	}
}

/********************************************************
*
* 
*
*********************************************************/
void calculate_diff_tend_full(int il,int ih,int jl,int jh){

	get_Kmix(il,ih,jl,jh);

	turbulent_diffusion_velocity(il,ih,jl,jh);
	turbulent_diffusion_scalars(il,ih,jl,jh);
}

/********************************************************
*
* 
*
*********************************************************/
void calculate_diff_tend_full_dry(int il,int ih,int jl,int jh){

	get_Kmix_dry(il,ih,jl,jh);

	turbulent_diffusion_velocity(il,ih,jl,jh);
	turbulent_diffusion_theta(il,ih,jl,jh);
}

/********************************************************
*
* 
*
*********************************************************/
void calculate_diff_tend_vertical(int il,int ih,int jl,int jh){
	
	get_Kmix_vertical(il,ih,jl,jh);
	turbulent_stress_vertical(il,ih,jl,jh);
}

/********************************************************
*
* 
*
*********************************************************/
void apply_velocity_diffusion(double step,int il,int ih,int jl,int jh){
	
	int ind,kl;
	double deltaT = step*dt;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
		
		kl = HTOPO(i,j)+1;	// lowest level to apply diffusion
		ind = INDEX(i,j,kl);
	
		UP(i,j,kl) += u_friction[ind] * deltaT;
		VP(i,j,kl) += v_friction[ind] * deltaT;
		
		for(int k=kl+1;k<NZ-1;k++){
		
			ind = INDEX(i,j,k);
		
			UP(i,j,k) += u_friction[ind] * deltaT;
			VP(i,j,k) += v_friction[ind] * deltaT;
			WP(i,j,k) += w_friction[ind] * deltaT;
		}
	}}
}

/********************************************************
*
* 
*
*********************************************************/
void apply_scalar_diffusion(double step,int il,int ih,int jl,int jh){
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=HTOPO(i,j)+1;k<NZ-1;k++){
		
		THP(i,j,k) += t_diffusion[INDEX(i,j,k)] * step * dt;
	}}}
}

/********************************************************
*
* 
*
*********************************************************/
void apply_moisture_diffusion(double step,int il,int ih,int jl,int jh){
	
	double deltaT = step*dt;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=HTOPO(i,j)+1;k<NZ-1;k++){
		
		QVP(i,j,k) += qv_diffusion[INDEX(i,j,k)] * deltaT;
		QCP(i,j,k) += qc_diffusion[INDEX(i,j,k)] * deltaT;
		QRP(i,j,k) += qr_diffusion[INDEX(i,j,k)] * deltaT;
	}}}
}

/********************************************************
*
* 
*
*********************************************************/
void initialize_landsea(const char *myfilename){
	
	landsea = (double*) calloc(NX*NY,sizeof(double));

	//-------------------------------------------------------
	// If using land-sea mask from file, get it from a file
	//-------------------------------------------------------
	if(USE_LANDSEA_FROM_FILE){
		//-------------------------------------------------------
		// Get dimensions of data in files to determine how 
		// much memory to allocate
		//-------------------------------------------------------
		size_t dims[3];
	
		get_dims(myfilename,"lon","lat","time",dims);
	
		int xdim = (int)dims[0];
		int ydim = (int)dims[1];

		if(!PARALLEL || rank == 0){
			printf("-------------------------------------------------------------------\n");
			printf("Getting Land-sea mask data from file:\n %s\n",myfilename);
			printf("-------------------------------------------------------------------\n");
		}

		double * myarray = get_data2(myfilename,lsmask_name,xdim*ydim);

		horz_interpolate(myarray,landsea,xdim,ydim,1,false,false,0.25,0.25,lonoffset+0.125,latoffset+(89.875-89.463));
	
		//-------------------------------------------------------
		// Set landsea mask, 1 for ocean, 0 for land
		//-------------------------------------------------------
		if(!MERIDIONAL_CROSS_SECTION){

			for(int i=0;i<NX;i++){
			for(int j=0;j<NY;j++){
		
				if(landsea[i*NY+j] < 0.5){ landsea[i*NY+j] = 0;}
				else { landsea[i*NY+j] = 1.0;}
			}}

		//-------------------------------------------------------
		// If using meridional cross section, make landsea mask
		// zonally uniform
		//-------------------------------------------------------	
		} else {
		
			int lonpoint = get_point_from_lon(meridional_lon_section);
		
			double *cross = (double*) calloc(NY,sizeof(double));
		
			printf("lon = %d",lonpoint);
		
			for(int j=0;j<NY;j++){ cross[j] = landsea[lonpoint*NY+j]; }
		
			for(int i=0;i<NX;i++){
			for(int j=0;j<NY;j++){
		
				if(cross[j] < 0.5){ landsea[i*NY+j] = 0;}
				else { landsea[i*NY+j] = 1.0;}
			}}
		
			free(cross);
		
		}

		free(myarray);
	
	//-------------------------------------------------------
	// If not using land-sea mask from a file, set it here 
	//-------------------------------------------------------	
	} else {
		
		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
		
			landsea[i*NY+j] = 1.0;
		}}
	}

}

/*********************************************************************
* Calculate the mixing coefficient K
**********************************************************************/
void get_Kmix_vertical(int il,int ih,int jl,int jh){
	
	double dUdz,dVdz,dudz,dvdz;
	double S2;
	double dTheta,bruntv;
	double l,c;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
		
		for(int k=HTOPO(i,j)+2;k<NZ;k++){
	
			/********************************************************	
			* Shear
			*********************************************************/
			dudz = 0.5 * ( UM(i,j,k) + UM(i+1,j,k) - UM(i,j,k-1) - UM(i+1,j,k-1) ) * ONE_D_DZW(k);
			dvdz = 0.5 * ( VM(i,j,k) + VM(i,j+1,k) - VM(i,j,k-1) - VM(i,j+1,k-1) ) * ONE_D_DZW(k);
		
			dUdz = 0.5 * ( UBAR(i,j,k) + UBAR(i+1,j,k) - UBAR(i,j,k-1) - UBAR(i+1,j,k-1) ) * ONE_D_DZW(k);
			dVdz = 0.5 * ( VBAR(i,j,k) + VBAR(i,j+1,k) - VBAR(i,j,k-1) - VBAR(i,j+1,k-1) ) * ONE_D_DZW(k);

			dUdz += dudz;
			dVdz += dvdz;
		
			S2 = (dUdz*dUdz+dVdz*dVdz);

			/********************************************************	
			* Brunt-Vaisala Frequency (dry air)
			*********************************************************/
			dTheta = 	( (THM(i,j,k  )+THBAR(i,j,k  )+tbv[k ])
						- (THM(i,j,k-1)+THBAR(i,j,k-1)+tbv[k-1]));

			bruntv = ( grav / ( THM(i,j,k)+THBAR(i,j,k)+tbv[k] ) ) * dTheta * ONE_D_DZW(k);
	
			/********************************************************	
			* If there is turbulence... (or just making sure we don't
			* take the square root of a negative number)
			*********************************************************/
			if(S2-one_d_Pr*bruntv > 0){
				
				l = DZW(k);
				
				c = l*l*cs*cs;
				
				KMIX(i,j,k)  = c*sqrt(S2-one_d_Pr*bruntv);
				
				if(KMIX(i,j,k) > vert_mix_Kmax[k]){ KMIX(i,j,k) = vert_mix_Kmax[k]; }
				
			} else {
				
				KMIX(i,j,k) = 0;
			}
		}
	}}
	
}

/*********************************************************************
*
*
**********************************************************************/
void turbulent_stress_vertical(int il,int ih,int jl,int jh){

	double wind_speed;
	double tau_u_l,tau_u_h,tau_v_l,tau_v_h;
	double tau_c_l,tau_c_h,tau_r_l,tau_r_h,tau_t_l,tau_t_h;
	double tau_q_h,tau_q_l;
	double ufric,vfric;

	int k;

	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		/********************************************************	
		* Calculate surface stress
		*********************************************************/
		k = HTOPO(i,j)+1;

		wind_speed = (UM(i,j,k)+UBAR(i,j,k))*(UM(i,j,k)+UBAR(i,j,k)) + 
					 (VM(i,j,k)+VBAR(i,j,k))*(VM(i,j,k)+VBAR(i,j,k));
		
		wind_speed = sqrt(wind_speed);

		tau_q_l = 0;
		tau_c_l = 0;
		tau_r_l = 0;
		tau_t_l = 0;
		
		/********************************************************	
		* Handle land/sea contrasts
		*********************************************************/
		if(LANDSEA(i,j) < 0.5){
		
			tau_u_l = 0.003 * wind_speed * UM(i,j,k);
			tau_v_l = 0.003 * wind_speed * VM(i,j,k);
			
		} else {
			
			tau_u_l = 0.001 * wind_speed * UM(i,j,k);
			tau_v_l = 0.001 * wind_speed * VM(i,j,k);
			
			#if 0 && USE_MICROPHYSICS
			//if(j+jbs[rank]<300 && i+ibs[rank]>182 && i+ibs[rank]<NX-1 ){
            if(j+jbs[rank]>177 && j+jbs[rank]<300 && i+ibs[rank]>182+111 && i+ibs[rank]<NX-1){
				drag_coef = (1.1e-3 + 4.0e-5*wind_speed)*LANDSEA(i,j);

				tau_c_l = 0;
				tau_r_l = 0;
				tau_t_l = drag_coef*wind_speed*(tmp_surface - (THM(i,j,k) +THBAR(i,j,k) + tb[k]) );
				tau_q_l = drag_coef*wind_speed*(qvs_surface - (QVM(i,j,k) + QBAR(i,j,k) + qb[k]) );
			}
			#endif
		}

		/********************************************************	
		* Apply turbulent diffusion to wind
		*********************************************************/
		for(int k=HTOPO(i,j)+1;k<NZ-1;k++){
			
			tau_u_h = 0.5*(KMIX(i,j,k+1)+KMIX(i-1,j,k+1)) * (UM(i,j,k+1)-UM(i,j,k)) * ONE_D_DZW(k+1);
			tau_v_h = 0.5*(KMIX(i,j,k+1)+KMIX(i,j-1,k+1)) * (VM(i,j,k+1)-VM(i,j,k)) * ONE_D_DZW(k+1);

			ufric = (tau_u_h-tau_u_l) * ONE_D_DZ(k);
			vfric = (tau_v_h-tau_v_l) * ONE_D_DZ(k);

			FRICTION(i,j,k) = -ufric*UP(i,j,k) - vfric*VP(i,j,k);

			u_friction[INDEX(i,j,k)] = ufric;
			v_friction[INDEX(i,j,k)] = vfric;
			w_friction[INDEX(i,j,k)] = 0;

			tau_u_l = tau_u_h;
			tau_v_l = tau_v_h;
			
			#if VORTICITY_BUDGET
			vort_ufric[INDEX(i,j,k)] += ufric * dt;
			vort_vfric[INDEX(i,j,k)] += vfric * dt;
			#endif
		}
		/********************************************************	
		* Apply turbulent diffusion to mixing ratio
		*********************************************************/
		for(int k=HTOPO(i,j)+1;k<NZ-1;k++){
			
			if(USE_MICROPHYSICS){
				
				tau_q_h = -KSMIX(i,j,k+1)*(QV(i,j,k+1) - QV(i,j,k))*ONE_D_DZW(k+1);		// vapor
			
				tau_c_h = -KSMIX(i,j,k+1)*(QC(i,j,k+1) - QC(i,j,k))*ONE_D_DZW(k+1);		// cloud

				tau_r_h = -KSMIX(i,j,k+1)*(QR(i,j,k+1) - QR(i,j,k))*ONE_D_DZW(k+1);		// rain

				QVP(i,j,k) -= dt * (tau_q_h-tau_q_l)*ONE_D_DZ(k);

				QCP(i,j,k) -= dt * (tau_c_h-tau_c_l)*ONE_D_DZ(k);
			
				QRP(i,j,k) -= dt * (tau_r_h-tau_r_l)*ONE_D_DZ(k);
			
				tau_q_l = tau_q_h;
				tau_c_l = tau_c_h;
				tau_r_l = tau_r_h;
			}
			//THP(i,j,k) -= dt * (tau_t_h-tau_t_l)*ONE_D_DZ(k);

			tau_t_h = -KSMIX(i,j,k+1)*(TH(i,j,k+1) - TH(i,j,k))*ONE_D_DZW(k+1);		// temperature

			t_diffusion[INDEX(i,j,k)] = -(tau_t_h-tau_t_l)*ONE_D_DZ(k);

			tau_t_l = tau_t_h;			
		}
	}}
}



/**********************************************************************************
* Calculate the mixing coefficient K. Vertical mixing stored on w-levels and 
* horizontal on u-levels
*
***********************************************************************************/
void get_Kmix(int il,int ih,int jl,int jh){
	
	double A,temp,thetaVH,thetaVL,qvsatH,qvsatL,pdH,pdL,tempL,tempH,qv,dQL,eslH,eslL;
	double S2;
	double dTheta,bruntv;
	double c;
	
	int k;
	double D11,D22,D33,D12,D13,D23;
		
	double cpRd = cp/Rd;
	
	/**************************************************************************
	*
	* 			CALCULATE MIXING COEFFICIENTS
	*
	***************************************************************************/
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
			
		KMIX(i,j,HTOPO(i,j)+1) = 0;	// zero at ground
		KMIX(i,j,NZ-1) = 0;
		#if 0
		if(PARALLEL &&  j+jbs[rank] == yp && i+ibs[rank] == xp ){
			qsat_test[1] = 0;
			thetaE[1] = 0;
		}
		#endif
		
		for(k=HTOPO(i,j)+2;k<NZ-1;k++){
			/********************************************************	
			* Rate of deformation tensor
			*********************************************************/	
			//-------------------------------------------------------------------------------------
			D11 = 0.50 * (UM(i+1,j,k) +   UM(i+1,j,k-1) -   UM(i-1,j,k) -   UM(i-1,j,k-1) + 
						UBAR(i+1,j,k) + UBAR(i+1,j,k-1) - UBAR(i-1,j,k) - UBAR(i-1,j,k-1)) * one_d_dx;
			//-------------------------------------------------------------------------------------
			D22 = 0.50 * (VM(i,j+1,k) +   VM(i,j+1,k-1) -   VM(i,j-1,k) -   VM(i,j-1,k-1) + 
						VBAR(i,j+1,k) + VBAR(i,j+1,k-1) - VBAR(i,j-1,k) - VBAR(i,j-1,k-1)) * one_d_dy;
			//-------------------------------------------------------------------------------------
			D33 = 1.0  * (WM(i,j,k+1) - WM(i,j,k-1) + WBAR(i,j,k+1) - WBAR(i,j,k-1)) * ONE_D_DZW(k);
			//-------------------------------------------------------------------------------------
			D12 = 0.5*( ( UM(i,j+1,k) +   UM(i,j+1,k-1) -   UM(i,j-1,k) -   UM(i,j-1,k-1) + 
						UBAR(i,j+1,k) + UBAR(i,j+1,k-1) - UBAR(i,j-1,k) - UBAR(i,j-1,k-1) ) * one_d_dy +
		                ( VM(i+1,j,k) +   VM(i+1,j,k-1) -   VM(i-1,j,k) -   VM(i-1,j,k-1) + 
						VBAR(i+1,j,k) + VBAR(i+1,j,k-1) - VBAR(i-1,j,k) - VBAR(i-1,j,k-1) ) * one_d_dx);
			//-------------------------------------------------------------------------------------
			D13 = 0.5*( ( UM(i,j,k) +   UM(i+1,j,k) -   UM(i,j,k-1) -   UM(i+1,j,k-1) + 
						UBAR(i,j,k) + UBAR(i+1,j,k) - UBAR(i,j,k-1) - UBAR(i+1,j,k-1) ) * ONE_D_DZW(k) + 
		           	    ( WM(i+1,j,k) - WM(i-1,j,k)  + WBAR(i+1,j,k) - WBAR(i-1,j,k)  ) * one_d_dx);
			//-------------------------------------------------------------------------------------
			D23 = 0.5*( ( VM(i,j,k) +   VM(i,j+1,k) -   VM(i,j,k-1) -   VM(i,j+1,k-1) + 
						VBAR(i,j,k) + VBAR(i,j+1,k) - VBAR(i,j,k-1) - VBAR(i,j+1,k-1) ) * ONE_D_DZW(k) +
	      	             ( WM(i,j+1,k) - WM(i,j-1,k) + WBAR(i,j+1,k) - WBAR(i,j-1,k) ) * one_d_dy);
			//-------------------------------------------------------------------------------------
			S2 = D11*D11 + D22*D22 + D33*D33 + D12*D12 + D13*D13 + D23*D23;

			/********************************************************	
			* Brunt-Vaisala Frequency (unsaturated)
			*********************************************************/
			if(!isSaturated[INDEX(i,j,k)] || !isSaturated[INDEX(i,j,k-1)]){
				//-----------------------------------------------------------------
				// Calculate full virtual potential temperature
				//-----------------------------------------------------------------
				thetaVH = (THM(i,j,k) + THBAR(i,j,k) + tb[k]) * ( 1 + 0.61*(QVM(i,j,k)+QBAR(i,j,k)+qb[k]));
				
				thetaVL = (THM(i,j,k-1) + THBAR(i,j,k-1) + tb[k-1]) * (1 + 0.61*(QVM(i,j,k-1)+QBAR(i,j,k-1)+qb[k-1]));

				dTheta = thetaVH - thetaVL;
				//-----------------------------------------------------------------
				// Brunt-Vaisala Frequency
				//-----------------------------------------------------------------
				bruntv = 2.0 * grav / ( tbv[k]+tbv[k-1]) * dTheta * ONE_D_DZW(k);
				#if 0
				if(PARALLEL &&  j+jbs[rank] == yp && i+ibs[rank] == xp ){
					//printf("%d %f %f %f\n",k,thetaVH,thetaVL,TH(i,j,k) + THBAR(i,j,k) + tb[k]);
					qsat_test[k] = 0;
					thetaE[k] = 0;
				}
				#endif
			/********************************************************	
			* Brunt-Vaisala Frequency (saturated)
			*********************************************************/	
			} else {
				//-----------------------------------------------------------------
				// Calculate full potential and full actual temperature
				//-----------------------------------------------------------------
				thetaVH = THM(i,j,k) + THBAR(i,j,k) + tb[k];
				thetaVL = THM(i,j,k-1) + THBAR(i,j,k-1) + tb[k-1];
				
				tempH = thetaVH * (PBAR(i,j,k  ) + PI(i,j,k  )/(cp*tbv[k  ]));
				tempL = thetaVL * (PBAR(i,j,k-1) + PI(i,j,k-1)/(cp*tbv[k-1]));
				
				temp = 0.5 * (tempH+tempL);
				//-----------------------------------------------------------------
				// Calculate vertical gradient of water substance and full mixing 
				// ratio for w-levels
				//-----------------------------------------------------------------				
				dQL = QCM(i,j,k) + QRM(i,j,k) - QCM(i,j,k-1) - QRM(i,j,k-1) + QVM(i,j,k) - QVM(i,j,k-1) + QBAR(i,j,k) - QBAR(i,j,k-1) + qb[k] - qb[k-1];
				
				qv = 0.5*(QVM(i,j,k)+QVM(i,j,k-1)+QBAR(i,j,k)+QBAR(i,j,k-1)+qb[k]+qb[k-1]);
				//-----------------------------------------------------------------
				// Calculate saturation mixing ratio
				//-----------------------------------------------------------------						
				pdH = p0*pow(PBAR(i,j,k  ),cpRd) + PI(i,j,k  )*rhou[k  ];	// full dimensional pressure at upper level
				pdL = p0*pow(PBAR(i,j,k-1),cpRd) + PI(i,j,k-1)*rhou[k-1];	// full dimensional pressure at lower level
				
				eslH = 611.2 * exp(17.67 * (tempH-273.15) / (tempH - 29.65) );	// saturation vapor pressure at upper level
				eslL = 611.2 * exp(17.67 * (tempL-273.15) / (tempL - 29.65) );	// saturation vapor pressure at lower level

				qvsatH = 0.62197 * eslH / (pdH-eslH);	// saturation mixing ratio at upper level
				qvsatL = 0.62197 * eslL / (pdL-eslL);	// saturation mixing ratio at lower level
				//-----------------------------------------------------------------
				// Calculate equivalent potential temperature and its gradient
				//-----------------------------------------------------------------				
				thetaVH = thetaVH * (1+Lv*qvsatH / (cp*tempH));
				thetaVL = thetaVL * (1+Lv*qvsatL / (cp*tempL));

				dTheta = thetaVH - thetaVL;		
				//-----------------------------------------------------------------
				// Output for testing purposes
				//-----------------------------------------------------------------
				#if 0
				if(PARALLEL &&  j+jbs[rank] == yp && i+ibs[rank] == xp ){
					
					qsat_test[k-1] = qvsatL;
					qsat_test[k] = qvsatH;
					thetaE[k-1] = thetaVL;
					thetaE[k] = thetaVH;
				}
				#endif
				//-----------------------------------------------------------------
				// Brunt-Vaisala Frequency
				//-----------------------------------------------------------------
				A = (grav / tbw[k]) * ( 1.0 + Lv*qv / (Rd*temp)) / (1.0+0.622*Lv*Lv*qv/(cp * Rd * temp*temp));
			
				bruntv = (A * dTheta - grav*dQL) * ONE_D_DZW(k);
			}
			
			/********************************************************	
			* If there is turbulence... (or just making sure we don't
			* take the square root of a negative number). Horizontal
			* mixing lengths are interpolated to half-levels
			*********************************************************/
			if(S2-one_d_Pr*bruntv > 0){
				
				c = vert_mixing_length_squared[k]*cs*cs;
				
				KMIX(i,j,k)  = c*sqrt(S2-one_d_Pr*bruntv);
			
				KHMIX( i,j,k-1) = 0.5*( KMIX(i,j,k-1) / vert_mixing_length_squared[k-1] + KMIX(i,j,k)/vert_mixing_length_squared[k] ) * lH*lH;
				
				if(KMIX(i,j,k) >vert_mix_Kmax[k]){ KMIX(i,j,k) = vert_mix_Kmax[k]; }
				if(KHMIX(i,j,k-1) > kmixh_max){ KHMIX(i,j,k-1) = kmixh_max; }			
							
			} else {
				
				KMIX(i,j,k) = 0;
				KHMIX(i,j,k-1) = 0 + 0.5 * KMIX(i,j,k-1) / vert_mixing_length_squared[k-1] * lH*lH;
			}
		}
	}}
	
	if(PARALLEL){
		exchangeOnePoint(Kmix);
		exchangeOnePoint(KHmix);
	}
}

/**********************************************************************************
* Calculate the mixing coefficient K. Vertical mixing stored on w-levels and 
* horizontal on u-levels
*
***********************************************************************************/
void get_Kmix_dry(int il,int ih,int jl,int jh){
	
	double thetaVH,thetaVL;
	double S2;
	double dTheta,bruntv;
	double c;
	
	int k;
	double D11,D22,D33,D12,D13,D23;
	
	/**************************************************************************
	*
	* 			CALCULATE MIXING COEFFICIENTS
	*
	***************************************************************************/
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
			
		KMIX(i,j,HTOPO(i,j)+1) = 0;	// zero at ground
		KMIX(i,j,NZ-1) = 0;
		
		for(k=HTOPO(i,j)+2;k<NZ-1;k++){
			/********************************************************	
			* Rate of deformation tensor
			*********************************************************/	
			//-------------------------------------------------------------------------------------
			D11 = 0.50 * (UM(i+1,j,k) +   UM(i+1,j,k-1) -   UM(i-1,j,k) -   UM(i-1,j,k-1) + 
						UBAR(i+1,j,k) + UBAR(i+1,j,k-1) - UBAR(i-1,j,k) - UBAR(i-1,j,k-1)) * one_d_dx;
			//-------------------------------------------------------------------------------------
			D22 = 0.50 * (VM(i,j+1,k) +   VM(i,j+1,k-1) -   VM(i,j-1,k) -   VM(i,j-1,k-1) + 
						VBAR(i,j+1,k) + VBAR(i,j+1,k-1) - VBAR(i,j-1,k) - VBAR(i,j-1,k-1)) * one_d_dy;
			//-------------------------------------------------------------------------------------
			D33 = 1.0  * (WM(i,j,k+1) - WM(i,j,k-1) + WBAR(i,j,k+1) - WBAR(i,j,k-1)) * ONE_D_DZW(k);
			//-------------------------------------------------------------------------------------
			D12 = 0.5*( ( UM(i,j+1,k) +   UM(i,j+1,k-1) -   UM(i,j-1,k) -   UM(i,j-1,k-1) + 
						UBAR(i,j+1,k) + UBAR(i,j+1,k-1) - UBAR(i,j-1,k) - UBAR(i,j-1,k-1) ) * one_d_dy +
		                ( VM(i+1,j,k) +   VM(i+1,j,k-1) -   VM(i-1,j,k) -   VM(i-1,j,k-1) + 
						VBAR(i+1,j,k) + VBAR(i+1,j,k-1) - VBAR(i-1,j,k) - VBAR(i-1,j,k-1) ) * one_d_dx);
			//-------------------------------------------------------------------------------------
			D13 = 0.5*( ( UM(i,j,k) +   UM(i+1,j,k) -   UM(i,j,k-1) -   UM(i+1,j,k-1) + 
						UBAR(i,j,k) + UBAR(i+1,j,k) - UBAR(i,j,k-1) - UBAR(i+1,j,k-1) ) * ONE_D_DZW(k) + 
		           	    ( WM(i+1,j,k) - WM(i-1,j,k)  + WBAR(i+1,j,k) - WBAR(i-1,j,k)  ) * one_d_dx);
			//-------------------------------------------------------------------------------------
			D23 = 0.5*( ( VM(i,j,k) +   VM(i,j+1,k) -   VM(i,j,k-1) -   VM(i,j+1,k-1) + 
						VBAR(i,j,k) + VBAR(i,j+1,k) - VBAR(i,j,k-1) - VBAR(i,j+1,k-1) ) * ONE_D_DZW(k) +
	      	             ( WM(i,j+1,k) - WM(i,j-1,k) + WBAR(i,j+1,k) - WBAR(i,j-1,k) ) * one_d_dy);
			//-------------------------------------------------------------------------------------
			S2 = D11*D11 + D22*D22 + D33*D33 + D12*D12 + D13*D13 + D23*D23;

			/********************************************************	
			* Brunt-Vaisala Frequency (unsaturated)
			*********************************************************/
			//-----------------------------------------------------------------
			// Calculate full virtual potential temperature
			//-----------------------------------------------------------------
			thetaVH = (THM(i,j,k) + THBAR(i,j,k) + tb[k]) * ( 1 + 0.61*qb[k]);
			
			thetaVL = (THM(i,j,k-1) + THBAR(i,j,k-1) + tb[k-1]) * (1 + 0.61*qb[k-1]);

			dTheta = thetaVH - thetaVL;
			//-----------------------------------------------------------------
			// Brunt-Vaisala Frequency
			//-----------------------------------------------------------------
			bruntv = 2.0 * grav / ( tbv[k]+tbv[k-1]) * dTheta * ONE_D_DZW(k);

			/********************************************************	
			* If there is turbulence... (or just making sure we don't
			* take the square root of a negative number). Horizontal
			* mixing lengths are interpolated to half-levels
			*********************************************************/
			if(S2-one_d_Pr*bruntv > 0){
				
				c = vert_mixing_length_squared[k]*cs*cs;
				
				KMIX(i,j,k)  = c*sqrt(S2-one_d_Pr*bruntv);
			
				KHMIX( i,j,k-1) = 0.5*( KMIX(i,j,k-1) / vert_mixing_length_squared[k-1] + KMIX(i,j,k)/vert_mixing_length_squared[k] ) * lH*lH;
				
				if(KMIX(i,j,k) >vert_mix_Kmax[k]){ KMIX(i,j,k) = vert_mix_Kmax[k]; }
				if(KHMIX(i,j,k-1) > kmixh_max){ KHMIX(i,j,k-1) = kmixh_max; }			
							
			} else {
				
				KMIX(i,j,k) = 0;
				KHMIX(i,j,k-1) = 0 + 0.5 * KMIX(i,j,k-1) / vert_mixing_length_squared[k-1] * lH*lH;
			}
		}
	}}
	
	if(PARALLEL){
		exchangeOnePoint(Kmix);
		exchangeOnePoint(KHmix);
	}
}


/**************************************************************************************
* Calculate diffusion tendencies for momentum and store them in var_friction arrays
* to apply throughout RK3 loop.
*
* il,ih,jl,jh - starting and ending indicies
****************************************************************************************/
void turbulent_diffusion_velocity(int il,int ih,int jl,int jh){

	double wind_speed,wind_speed_base,drag_coef;
	double ufric,vfric,wfric;

	int k,k_begin = 0;
	bool isEdge = false;
	
	double tau_u_zl,tau_v_zl;

	//-----------------------------------------------------------------
	// Mixing coefficients
	//-----------------------------------------------------------------
	double wK23H_h,wK23H_l,wK13H_h,wK13H_l,uK11H_h,vK22H_h,wK33V_h,uK11H_l,vK22H_l;
	double wK33V_l,uK13V_h,vK23V_h,uK13V_l,vK23V_l,uK12H_h,uK12H_l,vK12H_h,vK12H_l;

	/***********************************************************************************
	*
	* 				CALCULATE DIFFUSION TENDENCIES FOR ENTIRE DOMAIN	
	*
	************************************************************************************/
	for(int i=il-1;i<ih;i++){

		for(int j=jl-1;j<jh;j++){
			//---------------------------------------------------
			// Test for edges of domain, where lowest index will 
			// overstep HTOPO and LANDSEA array
			//---------------------------------------------------
			if(PARALLEL && (i+ibs[rank]-3 < 0 || j+jbs[rank]-3 < 0) ){
					 
				isEdge = true;
				k_begin = 1;
					
			} else {
					
				isEdge = false;
				k_begin = HTOPO(i,j)+1;
			}
			
			k = k_begin;

			/***********************************************************************************
			*	
			* 							SURFACE STRESS
			*
			************************************************************************************/
			//-----------------------------------------------------------------
			// Calculate wind speed
			//-----------------------------------------------------------------
			wind_speed = sqrt((UM(i,j,k)+UBAR(i,j,k))*(UM(i,j,k)+UBAR(i,j,k)) + (VM(i,j,k)+VBAR(i,j,k))*(VM(i,j,k)+VBAR(i,j,k)));
			wind_speed_base = sqrt(UBAR(i,j,k)*UBAR(i,j,k) + VBAR(i,j,k)*VBAR(i,j,k));
			//-----------------------------------------------------------------
			// Land surface stress
			//-----------------------------------------------------------------
			if(!isEdge && LANDSEA(i,j) < 0.5){	// check for LANDSEA and ISTOPO mismatch
					
				STRESS(j,0).tau_u_surface = 0.003 * (wind_speed * UM(i,j,k) + (wind_speed - wind_speed_base) * UBAR(i,j,k) );
				STRESS(j,0).tau_v_surface = 0.003 * (wind_speed * VM(i,j,k) + (wind_speed - wind_speed_base) * VBAR(i,j,k) );
				STRESS(j,k-1).tau_33 = 0;
			//-----------------------------------------------------------------
			// Ocean surface stress
			//-----------------------------------------------------------------		
			} else {
				
				drag_coef = (1.1e-3 + 4.0e-5*wind_speed);
						
				STRESS(j,0).tau_u_surface = drag_coef * (wind_speed * UM(i,j,k) + (wind_speed - wind_speed_base) * UBAR(i,j,k) );
				STRESS(j,0).tau_v_surface = drag_coef * (wind_speed * VM(i,j,k) + (wind_speed - wind_speed_base) * VBAR(i,j,k) );
				STRESS(j,k-1).tau_33 = 0;
			}		
			/***********************************************************************************
			*	
			* 						RATE OF DEFORMATION TENSOR
			*
			************************************************************************************/
			for(int k=k_begin;k<NZ-1;k++){
				//-----------------------------------------------------------------
				// For current i iteration, set fluxes on western side of new yz
				// cross section to those of previous eastern side
				//-----------------------------------------------------------------			
				STRESS(j,k).tau_11_west = STRESS(j,k).tau_11;
				STRESS(j,k).tau_12_west = STRESS(j,k).tau_12;
				STRESS(j,k).tau_13_west = STRESS(j,k).tau_13;
				//-----------------------------------------------------------------
				// Compression
				//-----------------------------------------------------------------
				STRESS(j,k).tau_11 = 2.0*(UM(i+1,j,k)-UM(i,j,k)) * one_d_dx;		
				STRESS(j,k).tau_22 = 2.0*(VM(i,j+1,k)-VM(i,j,k)) * one_d_dy;
				STRESS(j,k).tau_33 = 2.0*(WM(i,j,k+1)-WM(i,j,k)) * ONE_D_DZ(k);
				//-----------------------------------------------------------------
				// Shear
				//-----------------------------------------------------------------
				STRESS(j,k).tau_13 = (WM(i+1,j,k+1) - WM(i,j,k+1)) * one_d_dx + (UM(i+1,j,k+1) - UM(i+1,j,k)) * ONE_D_DZW(k+1);
				STRESS(j,k).tau_23 = (WM(i,j+1,k+1) - WM(i,j,k+1)) * one_d_dy + (VM(i,j+1,k+1) - VM(i,j+1,k)) * ONE_D_DZW(k+1);
				STRESS(j,k).tau_12 = (UM(i+1,j+1,k) - UM(i+1,j,k)) * one_d_dy + (VM(i+1,j+1,k) - VM(i,j+1,k)) * one_d_dx;
			}
		}
		/***********************************************************************************	
		* 
		* 	APPLY TURBULENT STRESS TO MOMENTUM FOR A YZ-CROSS SECTION
		*
		************************************************************************************/
		if(i!=il-1){
			
			for(int j=jl;j<jh;j++){
				//-----------------------------------------------------------------
				// SURFACE STRESS
				//-----------------------------------------------------------------
				tau_u_zl = STRESS(j,0).tau_u_surface;
				tau_v_zl = STRESS(j,0).tau_v_surface;
				
				integrated_friction_ke[INDEX2D(i,j)] = 0;
				
				for(int k=HTOPO(i,j)+1;k<NZ-1;k++){
					//-----------------------------------------------------------------
					// EDDY VISCOSITY FOR ZONAL WIND
					//-----------------------------------------------------------------
					uK11H_h = KHMIX(i  ,j,k);
					uK11H_l = KHMIX(i-1,j,k);
					uK12H_h = 0.25*(KHMIX(i,j  ,k) + KHMIX(i-1,j  ,k) + KHMIX(i-1,j+1,k) + KHMIX(i,j+1,k));
					uK12H_l = 0.25*(KHMIX(i,j-1,k) + KHMIX(i-1,j-1,k) + KHMIX(i-1,j  ,k) + KHMIX(i,j  ,k));
					uK13V_h = 0.50*(KMIX(i,j,k+1)+KMIX(i-1,j,k+1));
					//-----------------------------------------------------------------
					// EDDY VISCOSITY FOR MERIDIONAL WIND
					//-----------------------------------------------------------------	
					vK12H_h = 0.25*(KHMIX(i,j-1,k) + KHMIX(i,j,k) + KHMIX(i+1,j-1,k) + KHMIX(i+1,j,k));
					vK12H_l = uK12H_l;
					vK22H_h = KHMIX(i,j  ,k);
					vK22H_l = KHMIX(i,j-1,k);
					vK23V_h = 0.5*(KMIX(i,j,k+1)+KMIX(i,j-1,k+1));
					//-----------------------------------------------------------------
					// EDDY VISCOSITY FOR VERTICAL WIND
					//-----------------------------------------------------------------			
					wK13H_h = 0.25*(KHMIX(i,j,k)+KHMIX(i,j,k-1)+KHMIX(i+1,j,k)+KHMIX(i+1,j,k-1));
					wK13H_l = 0.25*(KHMIX(i,j,k)+KHMIX(i,j,k-1)+KHMIX(i-1,j,k)+KHMIX(i-1,j,k-1));
					wK23H_h = 0.25*(KHMIX(i,j,k)+KHMIX(i,j,k-1)+KHMIX(i,j+1,k)+KHMIX(i,j+1,k-1));
					wK23H_l = 0.25*(KHMIX(i,j,k)+KHMIX(i,j,k-1)+KHMIX(i,j-1,k)+KHMIX(i,j-1,k-1));
					wK33V_h = 0.50*(KMIX(i,j,k)+KMIX(i,j,k+1));
					wK33V_l = 0.50*(KMIX(i,j,k)+KMIX(i,j,k-1));
					//-----------------------------------------------------------------
					// TURBULENT DIFFUSION
					//-----------------------------------------------------------------
					ufric = (uK11H_h*STRESS(j,k).tau_11      - uK11H_l*STRESS(j  ,k).tau_11_west  ) * one_d_dx + 
							(uK12H_h*STRESS(j,k).tau_12_west - uK12H_l*STRESS(j-1,k).tau_12_west  ) * one_d_dy + 
							(uK13V_h*STRESS(j,k).tau_13_west - tau_u_zl				 			  ) * ONE_D_DZ(k);
			
					vfric = (vK12H_h*STRESS(j-1,k).tau_12 - vK12H_l*STRESS(j-1,k).tau_12_west  )	* one_d_dx + 
							(vK22H_h*STRESS(j  ,k).tau_22 - vK22H_l*STRESS(j-1,k).tau_22	   )	* one_d_dy + 
							(vK23V_h*STRESS(j-1,k).tau_23 - tau_v_zl					  	   ) 	* ONE_D_DZ(k);
			
					wfric = (wK13H_h*STRESS(j,k-1).tau_13 - wK13H_l*STRESS(j  ,k-1).tau_13_west)	* one_d_dx + 
							(wK23H_h*STRESS(j,k-1).tau_23 - wK23H_l*STRESS(j-1,k-1).tau_23)			* one_d_dy + 
							(wK33V_h*STRESS(j,k  ).tau_33 - wK33V_l*STRESS(j  ,k-1).tau_33)			* ONE_D_DZW(k);
					//-----------------------------------------------------------------
					// STORE VALUES FOR USE DURING RK3 LOOP
					//-----------------------------------------------------------------	
					u_friction[INDEX(i,j,k)] = ufric;
					v_friction[INDEX(i,j,k)] = vfric;
					w_friction[INDEX(i,j,k)] = wfric;
					
					tau_u_zl = uK13V_h*STRESS(j,k).tau_13_west;
					tau_v_zl = vK23V_h*STRESS(j-1,k).tau_23;
					
					//-----------------------------------------------------------------
					// 	OPTIONAL OUTPUT
					//-----------------------------------------------------------------	
					integrated_friction_ke[INDEX2D(i,j)] += -(ufric*UM(i,j,k) + vfric*VM(i,j,k)) * rhou[k] * DZU(k);
					
					#if OUTPUT_FRICTION_TEND
					FRICTION(i,j,k) = -ufric*UM(i,j,k) - vfric*VM(i,j,k);
					#endif
					
					#if VORTICITY_BUDGET
					vort_ufric[INDEX(i,j,k)] += ufric * dt;
					vort_vfric[INDEX(i,j,k)] += vfric * dt;
					#endif
				}
			}
		}
	}
}

/**********************************************************************************
* Calculate diffusion tendencies for scalar variables and store them
* in var_diffusion arrays to be applied throught RK3 loop
*
* il,ih,jl,jh - staring and ending indicies
***********************************************************************************/
void turbulent_diffusion_scalars(int il,int ih,int jl,int jh){
	
	double tau_c_l,tau_c_h,tau_r_l,tau_r_h,tau_t_l,tau_t_h;
	double tau_11_c_l,tau_22_c_l,tau_11_r_l,tau_22_r_l,tau_11_q_l,tau_22_q_l,tau_11_t_l,tau_22_t_l;
	double tau_11_c_h,tau_22_c_h,tau_11_r_h,tau_22_r_h,tau_11_q_h,tau_22_q_h,tau_11_t_h,tau_22_t_h;
	double drag_coef,tau_q_h,tau_q_l;
	double wind_speed,wind_speed_base;
	double pressure;
	double KH,KL;
	double pisfc;
	double gustFactor;
	double cpRd = cp/Rd;

	double tmp_surface_base,tmp_surface_pert,qvs_surface;

	int k;

	double dp_to_surface;

	//-----------------------------------------------------------------
	// pressure difference between lowest model level and surface
	//-----------------------------------------------------------------
	if(STRETCHED_GRID){
		dp_to_surface = grav * zsu[1] / (cp * tbv[1]);
	} else {
		dp_to_surface = grav * 0.5 * dz / (cp * tbv[1]);
	}
	
	/***********************************************************************************
	*	LOOP THROUGH ALL POINTS WITHIN GIVEN RANGE
	************************************************************************************/
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		k = HTOPO(i,j)+1;

		/***********************************************************************************
		*	
		* 							SURFACE FLUXES
		*
		************************************************************************************/
		//-----------------------------------------------------------------
		// Calculate wind speed
		//-----------------------------------------------------------------
		wind_speed = sqrt((UM(i,j,k)+UBAR(i,j,k))*(UM(i,j,k)+UBAR(i,j,k)) + (VM(i,j,k)+VBAR(i,j,k))*(VM(i,j,k)+VBAR(i,j,k)));
		wind_speed_base = sqrt(UBAR(i,j,k)*UBAR(i,j,k) + VBAR(i,j,k)*VBAR(i,j,k));
	
		gustFactor = 0;

		tau_c_l = 0;
		tau_r_l = 0;
		tau_q_l = 0;
		tau_t_l = 0;

		//-----------------------------------------------------------------
		// 	Calculate surface heat flux, if necessary
		//-----------------------------------------------------------------
		if(SURFACE_HEAT_FLUX){
			//-----------------------------------------------------------------
			// 							LAND
			//-----------------------------------------------------------------
			if(LANDSEA(i,j) < 0.5){

				tau_q_l = 0;
				tau_t_l = 0;
			//-----------------------------------------------------------------
			// 							OCEAN
			//-----------------------------------------------------------------		
			} else {
			
				drag_coef = (1.1e-3 + 4.0e-5*wind_speed)*LANDSEA(i,j)*ISTOPO(i,j,1);	// surface drag coefficient

				pressure = p0*pow(PBAR(i,j,k)+dp_to_surface,cpRd) + PI(i,j,k)*rhou[k];	// full surface dimensional pressure
				pisfc = PBAR(i,j,k) + PI(i,j,k)/(cp*tbv[k]) + dp_to_surface;			// full non-dimensional pressure at surface

				qvs_surface = 0.62197 * es_water / (pressure-es_water);					// surface saturation mixing ratio
			
				tmp_surface_base = (THBAR(i,j,k)+tb[k])*pisfc;							// basic state temperature at surface
				tmp_surface_pert = THM(i,j,k)*pisfc;									// perturbation temperature at surface
			
				// perturbation surface fluxes
				tau_q_l = drag_coef * ( ( wind_speed - wind_speed_base) * (qvs_surface -  QBAR(i,j,k) - qb[k]) - wind_speed * QVM(i,j,k) );
				tau_t_l = drag_coef * ( ( wind_speed - wind_speed_base) * (water_temp - tmp_surface_base) - wind_speed * tmp_surface_pert );
			}

			latent_heat_flux[INDEX2D(i,j)] = tau_q_l * Lv * rhou[k];	// store surface latent heat flux
		}
		
		/***********************************************************************************
		*
		* 				TURBULENT DIFFUSION FOR SCALARS
		*
		************************************************************************************/
		for(int k=HTOPO(i,j)+1;k<NZ-1;k++){

			//-----------------------------------------------------------------
			// 					VERTICAL
			//-----------------------------------------------------------------
			tau_q_h = -KSMIX(i,j,k+1)*(QVM(i,j,k+1) - QVM(i,j,k))*ONE_D_DZW(k+1);	// vapor flux top
			tau_c_h = -KSMIX(i,j,k+1)*(QCM(i,j,k+1) - QCM(i,j,k))*ONE_D_DZW(k+1);	// cloud flux top
			tau_r_h = -KSMIX(i,j,k+1)*(QRM(i,j,k+1) - QRM(i,j,k))*ONE_D_DZW(k+1);	// rain flux top
			tau_t_h = -KSMIX(i,j,k+1)*(THM(i,j,k+1) - THM(i,j,k))*ONE_D_DZW(k+1);	// heat flux top

			//-----------------------------------------------------------------
			// 					MERIDIONAL
			//-----------------------------------------------------------------
			KH = 0.5 * (KHSMIX(i,j+1,k)+KHSMIX(i,j,k));
			KL = 0.5 * (KHSMIX(i,j,k)+KHSMIX(i,j-1,k));

			tau_22_q_h = -KH*(QVM(i,j+1,k) - QVM(i,j,k))*one_d_dy;	// vapor flux north
			tau_22_q_l = -KL*(QVM(i,j,k) - QVM(i,j-1,k))*one_d_dy;	// vapor flux south

			tau_22_c_h = -KH*(QCM(i,j+1,k) - QCM(i,j,k))*one_d_dy;	// cloud flux north
			tau_22_c_l = -KL*(QCM(i,j,k) - QCM(i,j-1,k))*one_d_dy;	// cloud flux south

			tau_22_r_h = -KH*(QRM(i,j+1,k) - QRM(i,j,k))*one_d_dy;	// rain flux north
			tau_22_r_l = -KL*(QRM(i,j,k) - QRM(i,j-1,k))*one_d_dy;	// rain flux south

			tau_22_t_h = -KH*(THM(i,j+1,k) - THM(i,j,k))*one_d_dy;	// heat flux north
			tau_22_t_l = -KL*(THM(i,j,k) - THM(i,j-1,k))*one_d_dy;	// heat flux south

			//-----------------------------------------------------------------
			// 						ZONAL
			//-----------------------------------------------------------------
			KH = 0.5 * (KHSMIX(i+1,j,k)+KHSMIX(i,j,k));
			KL = 0.5 * (KHSMIX(i,j,k)+KHSMIX(i-1,j,k));

			tau_11_q_h = -KH*(QVM(i+1,j,k) - QVM(i,j,k))*one_d_dx;	// vapor flux east
			tau_11_q_l = -KL*(QVM(i,j,k) - QVM(i-1,j,k))*one_d_dx;	// vapor flux west

			tau_11_c_h = -KH*(QCM(i+1,j,k) - QCM(i,j,k))*one_d_dx;	// cloud flux east
			tau_11_c_l = -KL*(QCM(i,j,k) - QCM(i-1,j,k))*one_d_dx;	// cloud flux west

			tau_11_r_h = -KH*(QRM(i+1,j,k) - QRM(i,j,k))*one_d_dx;	// rain flux east
			tau_11_r_l = -KL*(QRM(i,j,k) - QRM(i-1,j,k))*one_d_dx;	// rain flux west

			tau_11_t_h = -KH*(THM(i+1,j,k) - THM(i,j,k))*one_d_dx;	// heat flux east
			tau_11_t_l = -KL*(THM(i,j,k) - THM(i-1,j,k))*one_d_dx; 	// heat flux west

			//-----------------------------------------------------------------
			// STORE DIFFUSION TENDENCIES FOR USE DURING RK3 LOOP
			//-----------------------------------------------------------------
			qv_diffusion[INDEX(i,j,k)] = - ( (tau_11_q_h-tau_11_q_l)*one_d_dx + (tau_22_q_h-tau_22_q_l)*one_d_dy + (tau_q_h-tau_q_l)*ONE_D_DZ(k) );
			qc_diffusion[INDEX(i,j,k)] = - ( (tau_11_c_h-tau_11_c_l)*one_d_dx + (tau_22_c_h-tau_22_c_l)*one_d_dy + (tau_c_h-tau_c_l)*ONE_D_DZ(k) );		
			qr_diffusion[INDEX(i,j,k)] = - ( (tau_11_r_h-tau_11_r_l)*one_d_dx + (tau_22_r_h-tau_22_r_l)*one_d_dy + (tau_r_h-tau_r_l)*ONE_D_DZ(k) );	
			 t_diffusion[INDEX(i,j,k)] = - ( (tau_11_t_h-tau_11_t_l)*one_d_dx + (tau_22_t_h-tau_22_t_l)*one_d_dy + (tau_t_h-tau_t_l)*ONE_D_DZ(k) );

			//-----------------------------------------------------------------
			// FLUXES ON TOP OF GRID VOLUME OF THIS ITERATION EQUAL FLUXES
			// ON THE BOTTOM OF THE NEXT ITERATION
			//-----------------------------------------------------------------
			tau_q_l = tau_q_h;
			tau_c_l = tau_c_h;
			tau_r_l = tau_r_h;			
			tau_t_l = tau_t_h;		
		}}
	}
}

/**********************************************************************************
* Calculate diffusion tendencies for scalar variables and store them
* in var_diffusion arrays to be applied throught RK3 loop
*
* il,ih,jl,jh - staring and ending indicies
***********************************************************************************/
void turbulent_diffusion_theta(int il,int ih,int jl,int jh){
	
	double tau_t_l,tau_t_h,tau_11_t_l,tau_22_t_l,tau_11_t_h,tau_22_t_h;
	double wind_speed,wind_speed_base;
	double KH,KL;

	int k;

	/***********************************************************************************
	*	LOOP THROUGH ALL POINTS WITHIN GIVEN RANGE
	************************************************************************************/
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		k = HTOPO(i,j)+1;
		
		/***********************************************************************************
		*	
		* 							SURFACE FLUXES
		*
		************************************************************************************/
		//-----------------------------------------------------------------
		// Calculate wind speed
		//-----------------------------------------------------------------
		wind_speed = sqrt((UM(i,j,k)+UBAR(i,j,k))*(UM(i,j,k)+UBAR(i,j,k)) + (VM(i,j,k)+VBAR(i,j,k))*(VM(i,j,k)+VBAR(i,j,k)));
		wind_speed_base = sqrt(UBAR(i,j,k)*UBAR(i,j,k) + VBAR(i,j,k)*VBAR(i,j,k));
	
		//-----------------------------------------------------------------
		// 							LAND
		//-----------------------------------------------------------------
		if(LANDSEA(i,j) < 0.5){

			tau_t_l = 0;
		//-----------------------------------------------------------------
		// 							OCEAN
		//-----------------------------------------------------------------		
		} else {

			tau_t_l = 0;
		}

		/***********************************************************************************
		*
		* 				TURBULENT DIFFUSION FOR SCALARS
		*
		************************************************************************************/
		for(int k=HTOPO(i,j)+1;k<NZ-1;k++){

			//-----------------------------------------------------------------
			// 					VERTICAL
			//-----------------------------------------------------------------
			tau_t_h = -KSMIX(i,j,k+1)*(THM(i,j,k+1) - THM(i,j,k))*ONE_D_DZW(k+1);	// heat flux top

			//-----------------------------------------------------------------
			// 					MERIDIONAL
			//-----------------------------------------------------------------
			KH = 0.5 * (KHSMIX(i,j+1,k)+KHSMIX(i,j,k));
			KL = 0.5 * (KHSMIX(i,j,k)+KHSMIX(i,j-1,k));

			tau_22_t_h = -KH*(THM(i,j+1,k) - THM(i,j,k))*one_d_dy;	// heat flux north
			tau_22_t_l = -KL*(THM(i,j,k) - THM(i,j-1,k))*one_d_dy;	// heat flux south

			//-----------------------------------------------------------------
			// 						ZONAL
			//-----------------------------------------------------------------
			KH = 0.5 * (KHSMIX(i+1,j,k)+KHSMIX(i,j,k));
			KL = 0.5 * (KHSMIX(i,j,k)+KHSMIX(i-1,j,k));

			tau_11_t_h = -KH*(THM(i+1,j,k) - THM(i,j,k))*one_d_dx;	// heat flux east
			tau_11_t_l = -KL*(THM(i,j,k) - THM(i-1,j,k))*one_d_dx; 	// heat flux west

			//-----------------------------------------------------------------
			// STORE DIFFUSION TENDENCIES FOR USE DURING RK3 LOOP
			//-----------------------------------------------------------------
			 t_diffusion[INDEX(i,j,k)] = - ( (tau_11_t_h-tau_11_t_l)*one_d_dx + (tau_22_t_h-tau_22_t_l)*one_d_dy + (tau_t_h-tau_t_l)*ONE_D_DZ(k) );

			//-----------------------------------------------------------------
			// FLUXES ON TOP OF GRID VOLUME OF THIS ITERATION EQUAL FLUXES
			// ON THE BOTTOM OF THE NEXT ITERATION
			//-----------------------------------------------------------------		
			tau_t_l = tau_t_h;		
		}}
	}
}

/*********************************************************************
* Calculate the mixing coefficient K
**********************************************************************/
void test(){
	
	const int xp = NX/2-1;// 91;//NX/2-1;
	const int yp = NY/2-1;//100;//NY/2-1;
	
	for(int i=3;i<fNX-3;i++){
	for(int j=3;j<fNY-3;j++){
	
		if(PARALLEL &&  j+jbs[rank] == yp && i+ibs[rank] == xp ){
	
			printf("coordinates = %f N %f E\n",outLats[yp],outLons[xp]);
	
			for(int k=1;k<NZ-1;k++){
		
				printf("%2d K = %f %f %f %+f %f %f %f %f %f %f\n",k,
				zsu[k]/1000.0,
				KMIX(i,j,k),
				KHMIX(i,j,k),
				rhou[k]*PI(i,j,k)*0.01,
				sqrt(UM(i,j,k)*UM(i,j,k)+VM(i,j,k)*VM(i,j,k)), 
				THM(i,j,k)+THBAR(i,j,k)+tb[k], 
				1000.0*(QVM(i,j,k)+QBAR(i,j,k)+qb[k]),
				1000.0*QCM(i,j,k),
				qsat_test[k]*1000.0,
				thetaE[k]);
			}
		}
	}}
	
}
