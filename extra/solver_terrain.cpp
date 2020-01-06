#include "stdafx.h"
#include "advection.h"
#include "terrain.h"
#include "interpolate.h"

#include "Heating.h"
#include "energy.h"


/*
*
*
*
*/

double div_avg[NX][NY];
double sor[3][NX][NY];


double mtime = 0;

double dtx = dt/dx;
double dty = dt/dy;
double dtz = dt/dz;

double rhoavg;

const bool isLinear = true;
const bool useMicrophysics = false;
const bool output_to_file = true;
const bool print_eke_budget = false;

int file_time_counter = 0;

int outfilefreq = 36;//144;//72;//36;
int stopheating = 576;//2304;//1152;//576;

double one_d_dx = 1. / dx;
double one_d_dy = 1. / dy;

Heating heat;

void advance_inner();
void advance_outer();

/*********************************************************************
* Integrate model equations
*
* step - how many time steps to integrate (can be fractional)
* itr - which iteration within Runge-Kutta loop (used to deterime 
*       first guess for SOR algorithm
*
**********************************************************************/
void pressure_solver(double step,int itr){

	/*****************************************************************************************
	* Calculate new momentum from everything but pressure
	******************************************************************************************/
	compute_velocity_tend_t(step);
	
	// boundary condition
	for(int k=1;k<NZ-1;k++){
		for(int i=0;i<NX;i++){ up[i][0][k] = up[i][1][k]; up[i][NY-1][k] = up[i][NY-2][k];}
		for(int j=0;j<NY;j++){ vp[0][j][k] = vp[1][j][k]; vp[NX-1][j][k] = vp[NX-2][j][k];}
	}

	/***************************************************************************************
	* Column density weighted average of divergence of advective tendencies
	****************************************************************************************/
	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	
		div_avg[i][j] = 0;

		// sum each vertical level
		for(int k=1;k<NZ-1;k++)
			div_avg[i][j] = div_avg[i][j] + ( DRU_DX(up)*one_d_dx + DRV_DY(vp)*one_d_dy );

		// vertical average
		div_avg[i][j] = div_avg[i][j] / ( cp * (double)(NZ-2) * step*dt );

	}}

	/**********************************************
	* Take off average value for the lid pressure 
	* at the boundaries to constrain the potential values
	***********************************************/
	double sum = 0;

	for(int i=0;i<NX;i++)  { sum = sum + sor[itr][i][1] + sor[itr][i][NY-2];}
	for(int j=1;j<NY-1;j++){ sum = sum + sor[itr][1][j] + sor[itr][NX-2][j];}

	sum = sum / (double)(2*NX+2*NY);

	// boundary condition
	for(int i=0;i<NX;i++){ sor[itr][i][0] = sor[itr][i][1]-sum; sor[itr][i][NY-1] = sor[itr][i][NY-2]-sum;}
	for(int j=0;j<NY;j++){ sor[itr][0][j] = sor[itr][1][j]-sum; sor[itr][NX-1][j] = sor[itr][1][j]-sum;}

	for(int i=0;i<NX;i++){ div_avg[i][0] = div_avg[i][1]; div_avg[i][NY-1] = div_avg[i][NY-2];}
	for(int j=0;j<NY;j++){ div_avg[0][j] = div_avg[1][j]; div_avg[NX-1][j] = div_avg[NY-2][j];}

	//**************************************************************************************


	/**********************************************************************
	* Solve Laplacian equation for pressure field required 
	* to balance divergence tendencies from advection using SOR
	***********************************************************************/

	// over-relaxation coefficient
	double omega = 2 / (1+trigpi/(double)NX);
	double dxdy = dx*dy;
	double omega_div_4 = 0.25*omega;
	double oneminusomega = 1-omega;

	// use previous values as first guess
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		pi[i][j][NZ-1] = sor[itr][i][j];
		div_avg[i][j] =  dxdy*div_avg[i][j];

	}}

	// SOR solver
	for(int q=1;q<200;q++){
//printf("sor = %d %f ",q,sor[itr][50][50]);
		for(int i=1;i<NX-1;i++){
		for(int j=1;j<NY-1;j++){
	
			sor[itr][i][j] = oneminusomega*sor[itr][i][j] + 

				omega_div_4*(   sor[itr][i-1][j] + sor[itr][i+1][j]+ sor[itr][i][j-1] + sor[itr][i][j+1]  - div_avg[i][j]  );
		}}
	}

	/*********************************************************************
	* Calculate hydrostatic pressure
	**********************************************************************/
	double tup; double tdn;

	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){

		pip[i][j][NZ-1] = 0;
		pip[i][j][NZ-2] = -0.5*(grav/cp)*((thp[i][j][NZ-2]*(1.+0.61*qb3d[i][j][NZ-2]))/(tbv3d[i][j][NZ-2]*tbv3d[i][j][NZ-2]))*dz*Gz[i][j];

		for(int k=NZ-3;k>0;k--){
			
			tup = (thp[i][j][k+1]*(1.+0.61*qb3d[i][j][k+1]))/(tbv3d[i][j][k+1]*tbv3d[i][j][k+1]);
			tdn = (thp[i][j][k]*(1.+0.61*qb3d[i][j][k]))/(tbv3d[i][j][k]*tbv3d[i][j][k]);
			pip[i][j][k] = pip[i][j][k+1]-0.5*(grav/cp)*(tup+tdn)*dz*Gz[i][j];
		}
		
	}}

	/*****************************************************
	*	Calculate column average hydrostatic pressure
	******************************************************/
	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	
		pip[i][j][0] = 0;
	
		for(int k=1;k<NZ-1;k++)
			pip[i][j][0] = pip[i][j][0] + grhos3d[i][j][k]*tbv3d[i][j][k]*pip[i][j][k];
			
		pip[i][j][0] = pip[i][j][0]/((double)(NZ-2));

	}}

	/*****************************************************
	* Calcuate top lid pressure
	******************************************************/
	for(int i=1;i<NX-1;i++)
		for(int j=1;j<NY-1;j++)
			pip[i][j][NZ-1] = (sor[itr][i][j] - pip[i][j][0]) / grhoavg2d[i][j];

	/*****************************************************
	* Calcuate total perturbation pressure
	* total pressure = lid pressure + hydrostatic pressure
	******************************************************/
	for(int i=0;i<NX;i++)
		for(int j=0;j<NY;j++)
			for(int k=1;k<NZ-1;k++)
				pi[i][j][k] = pip[i][j][NZ-1] + pip[i][j][k];

	/*****************************************************
	* Apply lateral boundary conditions
	******************************************************/
	for(int k=1;k<NZ-1;k++){

		for(int j=0;j<NY;j++){

			pi[0][j][k] = pi[1][j][k];
			pi[NX-1][j][k] = pi[NX-2][j][k];
		}
	
		for(int i=0;i<NX;i++){

			pi[i][0][k] = pi[i][1][k];
			pi[i][NY-1][k] = pi[i][NY-2][k];
		}
	}
//printf("u = %f %f\n",up[50][50][8],pi[50][50][8]);

	/*********************************************************************
	* Apply pressure gradient force to momentum equations
	**********************************************************************/
	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=1;k<NZ-1;k++){
		
		up[i][j][k] = up[i][j][k] + step*(	-dtx*cp*tbv3d[i][j][k]*DS_DX(pi)	);
		vp[i][j][k] = vp[i][j][k] + step*(	-dty*cp*tbv3d[i][j][k]*DS_DY(pi)	);

	}}}
		//printf("u = %f\n",up[50][50][8]);
	advect_theta_t(step);

}

/*****************************************************************************************
* Calculate flux divergences for momentum equations.
* Must call an interpolation routine first.
* 
* @param i,j,k - point in array
* @param step - fraction of total time step dt to advance
******************************************************************************************/
void uv_tend_t(double step,int i,int j,int k){

	v_at_u = 0.25 * (v[i][j][k]+v[i][j+1][k]+v[i-1][j][k]+v[i-1][j+1][k]);
	u_at_v = 0.25 * (u[i][j][k]+u[i+1][j][k]+u[i][j-1][k]+u[i+1][j-1][k]);

	/******************************************************
	* Zonal velocity
	*******************************************************/
	up[i][j][k] = um[i][j][k] + step*(
		(
			//---------------------------------------------------------------------
			// advection of zonal momentum by zonal wind
			- dtx*( XG_UU_FLUX_H - XG_UU_FLUX_L )
			//---------------------------------------------------------------------
			// advection of zonal momentum by meridional wind
			- dty*( YG_UV_FLUX_H - YG_UV_FLUX_L )
			//---------------------------------------------------------------------
			// advection of zonal momentum by vertical wind
			- dtz*( ZG_UW_FLUX_H - ZG_UW_FLUX_L )

		) / grhou3d[i][j][k]
		//---------------------------------------------------------------------
		// coriolis
		+ dt* f[j] * v_at_u
		//---------------------------------------------------------------------
		// friction
		- dt * friction[i][j][k] * u[i][j][k]
		)					
		;

	/******************************************************
	* Meridional velocity
	*******************************************************/
	vp[i][j][k] = vm[i][j][k] + step*(
		(
			//---------------------------------------------------------------------
			// advection of meridional momentum by zonal wind
			- dtx*( XG_VU_FLUX_H - XG_VU_FLUX_L)
			//---------------------------------------------------------------------
			// advection of meridional momentum by meridional wind
			- dty*( YG_VV_FLUX_H - YG_VV_FLUX_L)
			//---------------------------------------------------------------------
			// advection of meridional momentum by vertical wind
			- dtz*( ZG_VW_FLUX_H - ZG_VW_FLUX_L)

		) / grhov3d[i][j][k]
		//---------------------------------------------------------------------
		// coriolis
		- dt * f[j] * u_at_v
		//---------------------------------------------------------------------
		// friction
		- dt * friction[i][j][k] * v[i][j][k]
		)
		;

}


/*****************************************************************************************
* Calculate flux divergences for potential temperature equation.
* Must call an interpolation routine first.
* 
* @param i,j,k - point in array
* @param step - fraction of total time step dt to advance 
******************************************************************************************/
void theta_tend_t(double step,int i,int j,int k){

	thp[i][j][k] = thm[i][j][k] + step* (
		(
			//---------------------------------------------------------------------
			// advection of temperature by zonal wind
			- dtx * (XG_SU_FLUX_H - XG_SU_FLUX_L)
			//---------------------------------------------------------------------
			// advection of temperature by meridional wind	
			- dty * (YG_SV_FLUX_H - YG_SV_FLUX_L)
			//---------------------------------------------------------------------
			// advection of temperature by vertical wind			
			- dtz* (ZG_SW_FLUX_H - ZG_SW_FLUX_L)

		) / grhos3d[i][j][k]

		//---------------------------------------------------------------------
		// advection of model base state temperature by vertical wind
		- dtz * ZG_ADVECT_BASE(tb3d)
		)
	;
}

/*********************************************************************
* Loop through momentum arrays, interpolate from cell
* centers to cell faces, and then calculate rate of advection using
* the flux form of the advection equation.
*
* @param step - fractional of total time step dt to advance
*
**********************************************************************/
void compute_velocity_tend_t(double step){

	/*********************************
	* For each point far enough from
	* boundary to interpolate
	*********************************/
	for(int i=3;i<NX-3;i++){
	for(int j=3;j<NY-3;j++){

		//----------------------------
		// use lower order vertical 
		// interpolation for lowest point
		//----------------------------
		uv_zinterpolate2nd(i,j,1);	// vertical interpolation
		uv_interpolate5th(i,j,1);		// horizontal interpolation
		sdot_interp(i,j,1);

		uv_tend_t(step,i,j,1);	// flux divergence for single point

		/*********************************
		* For each point far enough from
		* upper and lower boundaries to 
		* interpolate
		*********************************/
		for(int k=2;k<NZ-2;k++){

			uv_interpolate5th(i,j,k);		// horizontal interpolation
			uv_zinterpolate3rd(i,j,k);	// vertical interpolation
			sdot_interp(i,j,k);

			uv_tend_t(step,i,j,k);	// flux divergence for single point
		}

		//----------------------------
		// use lower order vertical 
		// interpolation for highest point
		//----------------------------	
		uv_interpolate5th(i,j,NZ-2);		// horizontal interpolation
		uv_zinterpolate2nd(i,j,NZ-2);	// vertical interpolation
		sdot_interp(i,j,NZ-2);

		uv_tend_t(step,i,j,NZ-2);	// flux divergence for single point

	}}
	
}

/*********************************************************************
* Loop through potential temperature array, interpolate from cell
* centers to cell faces, and then calculate rate of advection using
* the flux form of the advection equation.
*
* @param step - fraction of total time step dt to advance

**********************************************************************/
void advect_theta_t(double step){

	/*********************************
	* For each point far enough from
	* boundary to interpolate
	*********************************/
	for(int i=3;i<NX-3;i++){
	for(int j=3;j<NY-3;j++){

		//----------------------------
		// use lower order vertical 
		// interpolation for lowest point
		//----------------------------
		th_interpolate5th(i,j,1);		// horizontal interpolation
		th_zinterpolate2nd(i,j,1);	// vertical interpolation

		theta_tend_t(step,i,j,1);	// flux divergence for single point

		/*********************************
		* For each point far enough from
		* upper and lower boundaries to 
		* interpolate
		*********************************/
		for(int k=2;k<NZ-2;k++){
		
			th_interpolate5th(i,j,k);		// horizontal interpolation
			th_zinterpolate3rd(i,j,k);	// vertical interpolation

			theta_tend_t(step,i,j,k);	// flux divergence for single point
		}

		//----------------------------
		// use lower order vertical 
		// interpolation for highest point
		//----------------------------
		th_interpolate5th(i,j,NZ-2);		// horizontal interpolation
		th_zinterpolate2nd(i,j,NZ-2);	// vertical interpolation

		theta_tend_t(step,i,j,NZ-2);	// flux divergence for single point
	}}

}

/*********************************************************************
* Calcuate vertical velocity using anelastic continuity equation
* from Lipps and Hemler (1982)
**********************************************************************/
void w_velocity_LH_t(){

	double one_d_x = 1. / dx;
	double one_d_y = 1. / dy;

	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){

		sdot[i][j][0] = 0;
		sdot[i][j][1] = 0;
		sdot[i][j][NZ-1] = 0;

		for(int k=2;k<NZ-1;k++){
			
			sdot[i][j][k] = (rhow3d[i][j][k-1]/rhow3d[i][j][k]) * sdot[i][j][k-1]

			+ (1./rhow3d[i][j][k])	*	(-(up[i+1][j][k-1]*grhou3d[i+1][j][k-1] - up[i][j][k-1]*grhou3d[i][j][k-1]) * one_d_x

		      		      				- (vp[i][j+1][k-1]*grhov3d[i][j+1][k-1] - vp[i][j][k-1]*grhov3d[i][j][k-1]) * one_d_y)*dz;
		}
	}}

	//for(int k=1;k<NZ-1;k++){ printf("%d %f\n",k,sdot[88][30][k]);}

}


/*********************************************************************
* Advance model foreward one full time step using third-order Runge-Kutta
* integration. One full step consists of three smaller steps.
*
**********************************************************************/
void integrate_rk3_t(){

	/****************************
	* first step
	*****************************/
	pressure_solver(0.333333333333333,0);

	w_velocity_LH_t();
	compute_sigma_dot(up,vp);

	upper_lower_boundaries();
	sponge_boundaries();

	if(useMicrophysics && bigcounter>=stopheating){

		sp_bound_microphysics();
		microphysics_advance_inner();
	}

	advance_inner();
	//printf("u = %e\n",up[90][30][8]);		
	/****************************
	* second step
	*****************************/
	pressure_solver(0.5,1);

	w_velocity_LH_t();
	compute_sigma_dot(up,vp);

	upper_lower_boundaries();
	sponge_boundaries();

	if(useMicrophysics && bigcounter>=stopheating){

		sp_bound_microphysics();
		microphysics_advance_inner();
	}

	advance_inner();

	/****************************
	* third step
	*****************************/
	pressure_solver(1.0,2);

	w_velocity_LH_t();
	compute_sigma_dot(up,vp);

	upper_lower_boundaries();
	sponge_boundaries();

	if(useMicrophysics && bigcounter>=stopheating)
		sp_bound_microphysics();
	

	//****************************
	// post integration

//	if(isLinear)
//		diffusion_pert(1.0);
//	else
//		diffusion(1.0);
	
	if(useMicrophysics && bigcounter>=stopheating){

		//microphysics_diffusion(1.0);
		kessler_microphysics();
	}
	

	/*********************************************
	* Apply heating perturbation
	**********************************************/
	if(bigcounter<stopheating){ heat.applyHeating();}

	upper_lower_boundaries();
	sponge_boundaries();


	if(useMicrophysics){
		sp_bound_microphysics();
		microphysics_advance();
	}

	rayleigh_damping();


	advance_outer();

}



/*********************************************************************
*
*
**********************************************************************/
void run_model_t(int count,FILE *infile){

	int counter = 0;	
	double elapsed;
	clock_t start_time,finis_time;

	/*******************************
	* Step through model 'count'
	* number of times
	********************************/
	while(counter<count){

		/*******************************
		* time each time step
		********************************/
		start_time = clock();

		integrate_rk3_t();	// run model forward one time step

		finis_time = clock();

		elapsed = ((double) (finis_time - start_time)) / CLOCKS_PER_SEC;

		printf("time %0.2f hr %0.3f s\n",mtime/3600,elapsed);

		/*******************************
		* print energy budget terms
		********************************/
		if(print_eke_budget && bigcounter % 12 == 0){ print_energy_budget(30,120,10,50,1,20,infile);}

		/*******************************
		* output to file
		********************************/
		if(output_to_file && bigcounter % outfilefreq == 0){

			write_all_pvars(filename,file_time_counter);
			file_time_counter++;

		}

		mtime = mtime + dt;	// elapsed physical time
		
		bigcounter++;	// total time steps (if this function is called multiple times)
		counter++;		// time steps within this loop
	}

}

/*********************************************************************
* Substeps within RK3 loop
*
**********************************************************************/
void advance_inner(){

	size_t num_bytes = NX*NY*NZ*sizeof(double);
	
	memcpy(u,up,num_bytes);
	memcpy(v,vp,num_bytes);
	memcpy(th,thp,num_bytes);

}

/*********************************************************************
* Final step for RK3
*
**********************************************************************/
void advance_outer(){

	size_t num_bytes = NX*NY*NZ*sizeof(double);

	memcpy(um,up,num_bytes);
	memcpy(u,up,num_bytes);

	memcpy(vm,vp,num_bytes);
	memcpy(v,vp,num_bytes);

	memcpy(thm,thp,num_bytes);
	memcpy(th,thp,num_bytes);

}

/*********************************************************************
*
*
**********************************************************************/
int main(){

	initialize2();

	compute_metric_terms();

	compute_density_terms();

	heat.initialize(21.18,86.3,19.37,93.0,100000.,6.0);

	heat.printInfo();

	if(output_to_file){ outfile_init(filename);}

	run_model_t(20000);

}


