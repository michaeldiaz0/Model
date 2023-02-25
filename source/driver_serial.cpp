#include "stdafx.h"
#include "advection.h"
#include "turbulence.h"
#include "Heating.h"
#include "energy.h"
#include "fluxes.h"
#include "damping.h"
#include "boundaries.h"
#include "initializer.h"
#include "pressure.h"
#include "budgets.h"

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
* Top-level initializer for serial model
**********************************************************************/
void initialize_serial(){
	
	set_outfilename(filename);

	initialize_basic_state();
	
	initialize_perturbation();
	
	if(USE_TURBULENT_STRESS){ init_kmix(NX,NY,NZ,&ZU(0));}
	
	if(USE_MICROPHYSICS){ init_microphysics(NX,NY);}
	
	initialize_landsea(landseaMaskFile);

	initialize_pressure_solver();
	
	init_boundaries(iebuffer,iwbuffer,jnbuffer,jsbuffer,-1);
	
	//init_stats();
	//heat.initialize(21.18,86.3,19.37,93.0,100000.,6.0);
	//heat.initialize(15.18,-5,15.37,5,150000.,6.0);
	//heat.initialize(8,270,8,276,300000.,6.0);
	//heat.shift(1.0,-2.0);

	//heat.printInfo();

	if(OUTPUT_TO_FILE){ outfile_init(filename);}
	
	init_damping(NX,NY,NZ);
	
	init_diffusion_weights(DIFFUSION_ORDER,&ZU(0));

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
		
		if(USE_EXPLICIT_DIFFUSION){ apply_explicit_diffusion(steps[s],3,NX-3,3,NY-3);}

		apply_boundary_condition(0);
		
		if(s<2){ advance_inner(NX*NY*NZ*sizeof(double));}

		if(USE_MICROPHYSICS){
			//if(bigcounter > 7200){ 
				advect_qv(steps[s],3,NX-3,3,NY-3);
				//}
			//set_terrain(&QVP(0,0,0),&ISTOPO(0,0,0),NX*NY*NZ);
			//diffusion_6th_var(&QVP(0,0,0),&QVM(0,0,0),steps[s],3,NX-3,3,NY-3);
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
		if(PERIODIC_BOUNDARIES){
			fft_damp(NX-6,NY,NZ,3,&UP(0,0,0));
			fft_damp(NX-6,NY,NZ,3,&VP(0,0,0));
			fft_damp(NX-6,NY,NZ,3,&WP(0,0,0));
			fft_damp(NX-6,NY,NZ,3,&THP(0,0,0));
			
			if(USE_MICROPHYSICS){
				fft_damp(NX-6,NY,NZ,3,&QVP(0,0,0));
			}
			
		} else {
			fft_damp2(NX,NY,NZ,0,&UP(0,0,0));
			fft_damp2(NX,NY,NZ,0,&VP(0,0,0));
			fft_damp2(NX,NY,NZ,0,&WP(0,0,0));
			fft_damp2(NX,NY,NZ,0,&THP(0,0,0));
		}
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
void s_integrate_rk3(){

	double steps[] = {1./3.,0.5,1.0};	// fractional time steps for RK
	int size = NX*NY*NZ;

	/*******************************************************
	* Apply explicit diffusion. Do this first, because the
	* final state should satisfy the anelastic continuity
	* equation and be free of super saturation.
	********************************************************/
	if(USE_EXPLICIT_DIFFUSION){ apply_explicit_diffusion(1.0,3,NX-3,3,NY-3);}

	/*******************************************************
	* Calculate frictional and diffusional tendencies to
	* be applied during RK3 loop
	********************************************************/
	if(USE_TURBULENT_STRESS){ calculate_diff_tend(1,NX-1,1,NY-1);}

	/*******************************************************
	* Runge-Kutta Loop
	********************************************************/
	for(int s=0;s<3;s++){
		
		// sign of advection for upwind biased derivatives
		compute_sign_cells(0,NX,0,NY);
		
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
		* vertical velocity for the new momentum field.
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
void s_run_model(int count,FILE *infile){

	int counter = 0;	
	double elapsed;
	double total_cputime = 0;
	bool isFirstStep = true;
	int timer_counter = 0;
	
	clock_t start_time = clock(),finis_time;

	clock_t start_time_full = clock();

	//------------------------------------------------------
	// Step through model 'count' number of times
	//------------------------------------------------------
	while(counter<count){

		optional_output(infile);

		//------------------------------------------------------
		// output to file
		//------------------------------------------------------
		if(OUTPUT_TO_FILE && bigcounter % outfilefreq == 0 && ENSEMBLE==0){

			output_meteorological_fields_to_file(write_pvar_to_file_2d,
												 write_pvar_to_file,
												 file_time_counter);
			
			write_time_to_file(filename,file_time_counter);
			
			
		}

		//------------------------------------------------------
		// Run model forward one time step
		//------------------------------------------------------
		if(ISLINEAR){ linear_integrate_rk3();}
		else { s_integrate_rk3();}

		if(inputs.print_courant_number){ print_courant_number_serial();}

		if(OUTPUT_TO_FILE && bigcounter % outfilefreq == 0 && ENSEMBLE==0){
			
			if(!isRestartRun || !isFirstStep)
				write_pvar_to_file(filename,"pi",pis,file_time_counter);
			
			file_time_counter++;
			
			if(VERBOSE && !isFirstStep){
				print_time_estimates(total_cputime,total_cputime,timer_counter);
			}
			
			total_cputime = 0;
			timer_counter = 0;
		}

		//------------------------------------------------------
		// Time each time step
		//------------------------------------------------------
		finis_time = clock();

		elapsed = ((double) (finis_time - start_time)) / CLOCKS_PER_SEC;

		total_cputime += elapsed;
		timer_counter += 1;
			
		if(VERBOSE){ printf("Time %d | Model time %0.3f hr | CPU time %0.3f s\n",bigcounter,mtime/3600,elapsed);}
		fflush(stdout);

		start_time = clock();

		mtime += dt;	// elapsed physical time
		
		bigcounter++;	// total time steps (if this function is called multiple times)
		counter++;		// time steps within this loop
		
		isFirstStep = false;
	}

	finis_time = clock();
	elapsed = ((double) (finis_time - start_time_full)) / CLOCKS_PER_SEC;

	printf("Program successfully completed. Total runtime: %f\n",elapsed);

}

/*********************************************************************
* 
**********************************************************************/
void run_serial_model(int count,FILE *infile){

	initialize_serial();
	
	s_run_model(count,infile);
}
