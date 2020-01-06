#include "stdafx.h"
#include "advection.h"
#include "surface.h"
#include "Heating.h"
#include "energy.h"
#include "damping.h"
#include "boundaries.h"
#include "initializer.h"
#include "pressure.h"
#include "temperature.h"

/*******************************************************************************************
* 
* TOP LEVEL SUBROUTINES FOR SOLVING MODEL EQUATIONS
*
********************************************************************************************/

size_t file_time_counter = 0;

Heating heat;

void diff_terrain_vars(int size);

/*********************************************************************
* Integrate hydrostatic model equations. Advects all wind components,
* calculates anelastics pressure field, and then applies pressure
* gradient force to momentum field.
*
* step - fractional size of time step
* itr - which iteration within Runge-Kutta loop (used to deterime 
*       first guess for SOR algorithm
* il,ih,jl,jh - starting and ending grid indices
**********************************************************************/
void integrate_hydro(double step,int itr,int il,int ih,int jl, int jh){

	int buffer,size;

	if(PARALLEL){
		buffer = 3;
		size = fNX*fNY*fNZ;
	} else {
		buffer = 1;
		size = NX*NY*NZ;
	}
	/*****************************************************************
	* Calculate new momentum from everything but pressure
	******************************************************************/
	advect_uv_velocity(step,il+3,ih-3,jl+3,jh-3);

	/*****************************************************************
	* Zero out values at or below the terrain for velocity field
	******************************************************************/
	if(USE_TERRAIN){
		set_terrain(ups,&UISTOPO(0,0,0),size);
		set_terrain(vps,&VISTOPO(0,0,0),size);
	}

	/*****************************************************************
	* Solve anelastic pressure equation
	******************************************************************/
	solve_pressure(step);

	/*****************************************************************
	* Apply pressure gradient force to momentum equations
	******************************************************************/
	for(int i=il+buffer;i<ih-buffer;i++){
	for(int j=jl+buffer;j<jh-buffer;j++){
	for(int k=1;k<NZ-1;k++){
		
		UP(i,j,k) -= step*(dtx*cp*tbv[k]*(PI(i,j,k)-PI(i-1,j,k)));
		VP(i,j,k) -= step*(dty*cp*tbv[k]*(PI(i,j,k)-PI(i,j-1,k)));
		
	}}}
}

/*********************************************************************
* Integrate non-hydrostatic model equations. Advects all wind components,
* calculates anelastic pressure field, and then applies pressure
* gradient force to velocity field.
*
* step - fractional size of time step
* il,ih,jl,jh - starting and ending grid indices
**********************************************************************/
void integrate_non_hydro(double step,int il,int ih,int jl, int jh){

	int buffer,size;

	if(PARALLEL){
		buffer = 3;
		size = fNX*fNY*fNZ;
	} else {
		buffer = 1;
		size = NX*NY*NZ;
	}
	//-----------------------------------------------------------------
	// Calculate new velocity from everything but pressure
	//-----------------------------------------------------------------
	advect_uvw_velocity(step,il+3,ih-3,jl+3,jh-3);
	
	if(USE_TURBULENT_STRESS){ apply_velocity_diffusion(step,il+3,ih-3,jl+3,jh-3);}
	if(EXTRA_DIFFUSION && !STRETCHED_GRID){ diffusion_2nd_all(step,3,NX-3,3,NY-3);}//diffusion_6th_all(step,3,NX-3,3,NY-3);}
	
	//-----------------------------------------------------------------
	// Zero out values at or below the terrain for velocity field
	//-----------------------------------------------------------------
	if(USE_TERRAIN){
		set_terrain(ups,&UISTOPO(0,0,0),size);
		set_terrain(vps,&VISTOPO(0,0,0),size);
		set_terrain(wps, &ISTOPO(0,0,0),size);	// need a WISTOPO array!
	}

	//-----------------------------------------------------------------
	// Solve anelastic pressure equation
	//-----------------------------------------------------------------
	solve_pressure(step);
	
	//-----------------------------------------------------------------
	// Apply pressure gradient force to velocity field
	//-----------------------------------------------------------------
	for(int i=il+buffer;i<ih-buffer;i++){
	for(int j=jl+buffer;j<jh-buffer;j++){

		UP(i,j,1) -= step*(dtx*(PI(i,j,1)-PI(i-1,j,1)));
		VP(i,j,1) -= step*(dty*(PI(i,j,1)-PI(i,j-1,1)));

		for(int k=2;k<NZ-1;k++){
		
			UP(i,j,k) -= step*(    dtx*(PI(i,j,k)-PI(i-1,j,k)));
			VP(i,j,k) -= step*(    dty*(PI(i,j,k)-PI(i,j-1,k)));
			WP(i,j,k) -= step*(DTZW(k)*(PI(i,j,k)-PI(i,j,k-1)));
		}
	}}
	
	//if(PV_TRACER && step > 0.99){ diff_terrain_vars(size);}
}

/*********************************************************************
* Calcuate vertical velocity using anelastic continuity equation
* from Lipps and Hemler (1982)
**********************************************************************/
void w_velocity_LH(int il,int ih,int jl,int jh){

	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		W(i,j,HTOPO(i,j)) = 0;
		W(i,j,HTOPO(i,j)+1) = 0;
		W(i,j,NZ-1) = 0;

		for(int k=2;k<NZ-1;k++){
			
			W(i,j,k) = (rhow[k-1]/rhow[k])*W(i,j,k-1)
				+ (rhou[k-1]/rhow[k])*
					( 
						-(UP(i+1,j,k-1) - UP(i,j,k-1) ) * one_d_dx	
		      			-(VP(i,j+1,k-1) - VP(i,j,k-1) ) * one_d_dy
					)*dz/mu[k];
		}
	}}
}

/*********************************************************************
* Calcuate vertical velocity using anelastic continuity equation
* from Durran (1989)
**********************************************************************/
void w_velocity(){

	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){

		W(i,j,HTOPOFULL(i,j)) = 0;	// boundary condition - zero vertical 
		W(i,j,HTOPOFULL(i,j)+1) = 0;	// velocity at top and bottom
		W(i,j,NZ-1) = 0;

		for(int k=HTOPOFULL(i,j)+2;k<NZ-1;k++){
	
			W(i,j,k) = ( (tbw[k-1]*rhow[k-1]) / (tbw[k-1]*rhow[k])  ) * W(i,j,k-1) 
				+ ( (tb[k-1]*rhou[k-1]) / (tbw[k]*rhow[k] ))*
					(
						-(UP(i+1,j,k-1) - UP(i,j,k-1)) * one_d_dx  	
		      			-(VP(i,j+1,k-1) - VP(i,j,k-1)) * one_d_dy
					)*dz;
		}
	}}

}

/*****************************************************************
* Zero out values at or below the terrain
******************************************************************/
void set_terrain(double * var,double *topo,int size){

	for(int i=0;i<size;i++){ var[i] *= topo[i];}
}

/*****************************************************************
* 
******************************************************************/
void diff_terrain(double * varIn,double * varOut,double *topo,int size){

	for(int i=0;i<size;i++){ 
		
		if(topo[i]<0.1){ varOut[i] = -varIn[i];} 
		else 		   { varOut[i] = 0; 		 }
	}
}

/*****************************************************************
* 
******************************************************************/
void diff_terrain_vars(int size){

	diff_terrain(ups,u_bound,uistopos,size);
	diff_terrain(vps,v_bound,vistopos,size);
	diff_terrain(wps,w_bound,istopos,size);
	diff_terrain(thps,w_bound,istopos,size);

}

/*********************************************************************
* Substeps within RK3 loop
*
**********************************************************************/
void switch_array(double **var1,double **var2){
	
	double *temp = *var1;
	
	*var1 = *var2;
	*var2 = temp;
}

/*********************************************************************
* Substeps within RK3 loop
*
**********************************************************************/
void advance_inner(size_t num_bytes){
	
	switch_array(&us,&ups);
	switch_array(&vs,&vps);
	switch_array(&ths,&thps);

	if(!HYDROSTATIC){ switch_array(&ws,&wps);}	
}

/*********************************************************************
* Final step for RK3
*
**********************************************************************/
void advance_outer(size_t num_bytes){

	memcpy(ums,ups,num_bytes);
	memcpy(vms,vps,num_bytes);
	memcpy(thms,thps,num_bytes);

	if(!HYDROSTATIC){ memcpy(wms,wps,num_bytes);}

	advance_inner(num_bytes);
}

/*********************************************************************
*
*
**********************************************************************/
void optional_output(FILE *infile){
	
	if(ENSEMBLE==2){
		
		if(bigcounter % 30 == 0){ get_stats(infile);}
	}
	
	//------------------------------------------------------
	// print energy budget terms
	//------------------------------------------------------
	if(PRINT_EKE_BUDGET && bigcounter % 3 == 0){
		
		if(!ISLINEAR){ //print_energy_budget(30,120,10,50,1,28,infile);
		} else { 	  //print_energy_budget(30,120,10,50,1,20,infile);
			if(MERIDIONAL_CROSS_SECTION){
				print_energy_budget(3,NX-3,5,NY-6,1,20);
			} else {
				print_energy_budget(5,NX-5,5,NY-5,1,35);
			}
		 }
	}
	
}

/*********************************************************************
* Advance model foreward one full time step using third-order Runge-Kutta
* integration. One full step consists of three smaller steps.
*
**********************************************************************/
void linear_integrate_rk3(){

	double steps[] = {1./3.,0.5,1.0};	// fractional time steps for RK

	/*******************************************************
	* Runge-Kutta Loop
	********************************************************/
	for(int s=0;s<3;s++){

		// forecast wind and pressure field
		integrate_non_hydro(steps[s],0,NX,0,NY);
		
		// forecast temperature field
		advect_theta(steps[s],3,NX-3,3,NY-3);
		
		if(USE_TERRAIN){ set_terrain(&THP(0,0,0),&ISTOPO(0,0,0),NX*NY*NZ);}
		
		if(EXTRA_DIFFUSION && !STRETCHED_GRID){ diffusion_6th_var(&THP(0,0,0),&THM(0,0,0),steps[s],3,NX-3,3,NY-3);}

		apply_boundary_condition(0);
		
		if(s<2){ advance_inner(NX*NY*NZ*sizeof(double));}

		if(USE_MICROPHYSICS){
			//if(bigcounter > 7200){ 
				advect_qv(steps[s],3,NX-3,3,NY-3);
				//}
			//set_terrain(&QVP(0,0,0),&ISTOPO(0,0,0),NX*NY*NZ);
			diffusion_6th_var(&QVP(0,0,0),&QVM(0,0,0),steps[s],3,NX-3,3,NY-3);
			apply_boundary_condition_microphysics(0);
			if(s<2){ microphysics_advance_inner(NX*NY*NZ*sizeof(double));}
		}
	
		

	}
	
	//large_scale_precip(5,NX-5,get_point_from_lat(17),get_point_from_lat(27));
	//large_scale_precip(5,NX-5,5,NY-5);
	
	/*********************************************
	* Apply heating perturbation
	**********************************************/
	if(bigcounter<stopheating){ heat.applyHeating_random();}
	
	#if FOURIER_DAMPING
		#if PERIODIC_BOUNDARIES
			fft_damp(NX-6,NY,NZ,3,&UP(0,0,0));
			fft_damp(NX-6,NY,NZ,3,&VP(0,0,0));
			fft_damp(NX-6,NY,NZ,3,&WP(0,0,0));
			fft_damp(NX-6,NY,NZ,3,&THP(0,0,0));
			
			if(USE_MICROPHYSICS){
				fft_damp(NX-6,NY,NZ,3,&QVP(0,0,0));
			}
			
		#else
			fft_damp2(NX,NY,NZ,0,&UP(0,0,0));
			fft_damp2(NX,NY,NZ,0,&VP(0,0,0));
			fft_damp2(NX,NY,NZ,0,&WP(0,0,0));
			fft_damp2(NX,NY,NZ,0,&THP(0,0,0));
		#endif
	#endif
			
	apply_boundary_condition(0);
	

	rayleigh_damping(1,NX-1,1,NY-1,raydampheight);	// should this occur before pressure solver?

	advance_outer(NX*NY*NZ*sizeof(double));
	
	if(USE_MICROPHYSICS){
		apply_boundary_condition_microphysics(0);	
		microphysics_advance(NX*NY*NZ*sizeof(double));
	}


	if(bigcounter % (12*60*60 / (int)dt) == 0 && bigcounter != 0){
	//if(bigcounter % (48) == 0 && bigcounter != 0){
	//if(bigcounter % (144) == 0 && bigcounter != 0){

		double max = find_max2(&VP(0,0,0));
		rescale_pert2(max,5.0);
		printf("Rescaling Perturbation\n");
		//double max = find_max2(&PI(0,0,0));
		//rescale_pert2(max,200.0);
	}
	
}

/*********************************************************************
* Advance model foreward one full time step using third-order Runge-Kutta
* integration. One full step consists of three smaller steps.
*
**********************************************************************/
void integrate_rk3(){

	double steps[] = {1./3.,0.5,1.0};	// fractional time steps for RK
	int size = NX*NY*NZ;

	/*******************************************************
	* Calculate frictional and diffusional tendencies to
	* be applied during RK3 loop
	********************************************************/
	if(USE_TURBULENT_STRESS){ calculate_diff_tend(1,NX-1,1,NY-1);}

	/*******************************************************
	* Runge-Kutta Loop
	********************************************************/
	for(int s=0;s<3;s++){
		
		/*******************************************************
		* Solve momentum and pressure equations using either
		* the hydrostatic or non-hydrostatic equation set
		********************************************************/
		if(HYDROSTATIC){
			integrate_hydro(steps[s],s,0,NX,0,NY);
		} else {
			integrate_non_hydro(steps[s],0,NX,0,NY);
		}
		/*******************************************************
		* Solve for new potential temperature field
		********************************************************/
		advect_theta(steps[s],3,NX-3,3,NY-3);
			
		if(USE_TURBULENT_STRESS){ apply_scalar_diffusion(steps[s],1,NX-1,1,NY-1);}
		
		if(USE_TERRAIN){ set_terrain(thps,&ISTOPO(0,0,0),NX*NY*NZ);}

		/*******************************************************
		* Advect microphysics variables
		********************************************************/
		if(USE_MICROPHYSICS){

			advect_microphysics_cell(steps[s],3,NX-3,3,NY-3);
			
			if(USE_ICE){ advect_ice_cell(steps[s],3,NX-3,3,NY-3);}
			
			if(USE_TURBULENT_STRESS){ apply_moisture_diffusion(steps[s],1,NX-1,1,NY-1);}
			
			zero_moisture(1,NX-1,1,NY-1,size);
		}
		/*******************************************************
		* If using the hydrostatic equation set, calculate the
		* velocity for the new momentum field.
		********************************************************/
		if(HYDROSTATIC){ w_velocity_LH(1,NX-1,1,NY-1);}

		/*******************************************************
		* Handle all boundary conditions
		********************************************************/
		apply_boundary_condition(0);

		if(USE_MICROPHYSICS){

			apply_boundary_condition_microphysics(0);

			if(s<2){ microphysics_advance_inner(NX*NY*NZ*sizeof(double));}
		}
		
		if(s<2){ advance_inner(NX*NY*NZ*sizeof(double));}
	}
	/***************************************************
	* POST INTEGRATION. NOW RUN PHYSICAL PARAMETERIZATIONS
	****************************************************/
	if(USE_MICROPHYSICS){
		
		run_microphysics(1,NX-1,1,NY-1);
		
		if(USE_TERRAIN){
			set_terrain(qvps,&ISTOPO(0,0,0),size);
			set_terrain(qrps,&ISTOPO(0,0,0),size);
			set_terrain(qcps,&ISTOPO(0,0,0),size);
		}
		
		apply_boundary_condition_microphysics(0);
		
		microphysics_advance(NX*NY*NZ*sizeof(double));
	}
	
	/*********************************************
	* Apply heating perturbation
	**********************************************/
	if(bigcounter<stopheating){ heat.applyHeating();}

	/*********************************************
	* BOUNDARY CONDITIONS
	**********************************************/
	apply_boundary_condition(0);

	rayleigh_damping(1,NX-1,1,NY-1,raydampheight);	// should this occur before pressure solver?

	advance_outer(NX*NY*NZ*sizeof(double));
}


/*********************************************************************
* Run model forward 'count' number of steps. Can call multiple times.
*
* @param count - number of steps to advance model
* @param infile - optional parameter to write energy budget data to file
**********************************************************************/
void run_model(int count,FILE *infile){

	int counter = 0;	
	double elapsed;
	clock_t start_time = clock(),finis_time;

	//------------------------------------------------------
	// Step through model 'count' number of times
	//------------------------------------------------------
	while(counter<count){

		//------------------------------------------------------
		// Time each time step
		//------------------------------------------------------
		finis_time = clock();

		elapsed = ((double) (finis_time - start_time)) / CLOCKS_PER_SEC;
		
		if(VERBOSE){ printf("time %0.3f hr %0.3f s\n",mtime/3600,elapsed);}

		start_time = clock();

		optional_output(infile);

		//------------------------------------------------------
		// output to file
		//------------------------------------------------------
		if(OUTPUT_TO_FILE && bigcounter % outfilefreq == 0 && ENSEMBLE==0){

			write_all_pvars(filename,file_time_counter);
			
			#if OUTPUT_FRICTION_TEND && USE_TURBULENT_STRESS && !PARALLEL
				write_pvar_to_file(filename,"fric",frictions,file_time_counter);
			#endif
			
			#if OUTPUT_DIFFUSION_TEND && !PARALLEL
				//write_pvar_to_file(filename,"diff",diff_tend,NX,NY,file_time_counter);
			#endif
			
			file_time_counter++;
		}

		//------------------------------------------------------
		// Run model forward one time step
		//------------------------------------------------------
		if(ISLINEAR){ linear_integrate_rk3();}
		else { integrate_rk3();}

		mtime += dt;	// elapsed physical time
		
		bigcounter++;	// total time steps (if this function is called multiple times)
		counter++;		// time steps within this loop
	}

}
