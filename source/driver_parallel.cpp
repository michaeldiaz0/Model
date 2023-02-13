#include "stdafx.h"
#include "Heating.h"
#include "advection.h"
#include "energy.h"
#include "fluxes.h"
#include "turbulence.h"
#include "budgets.h"
#include "damping.h"
#include "initializer.h"
#include "boundaries.h"
#include "pressure.h"
#include "mpi_setup.h"
#include "pcomm.h"

/*******************************************************************************************
* DRIVER FOR THE PARALLEL VERSION OF THE MODEL
*
*
*
*
*
*
*
*
********************************************************************************************/

void p_run_model(int,FILE *infile=NULL);

/*********************************************************************
* Top-level initializer for parallel model
* Setup parallel environment.
*
* argc - number of command-line arguments
* argv[] - command-line arguments
**********************************************************************/
void initialize_parallel_model(){
	

	//------------------------------------------------
	// Initialize grid dimensions, indices, location
	// of neighbors, etc., and setupt data types
	// for interprocess communication
	//------------------------------------------------
	setup_grid();
 
	initAllComms();

	set_outfilename(filename);
	//------------------------------------------------
	// Initialize input data on root process
	//------------------------------------------------
	if(rank==0){ initialize_basic_state();}
	//------------------------------------------------
	// Initialize subarrays for each process
	//------------------------------------------------
	initialize_subarray(fNX,fNY,fNZ);
	//------------------------------------------------
	// Broadcast input data from the root process 
	// to all processes
	//------------------------------------------------
	broadcast_shared_data();

	initialize_landsea(landseaMaskFile);

	//if(OUTPUT_TO_FILE && rank==0){ outfile_init(filename);}
	if(OUTPUT_TO_FILE){ outfile_init(filename);}
	
	initialize_perturbation();

	if(rank==0){
		
		free(iubar);free(ivbar);free(iwbar);free(ithbar);free(iqbar); free(ipbar);
		free(iistopo);free(iuistopo);free(ivistopo);free(ifriction);//free(itopo);
	}

	if(USE_TURBULENT_STRESS){ init_kmix(fNX,fNY,fNZ,&ZU(0));}

	initialize_flux_cells(fNY,fNZ);
	initialize_microphysics_cells(fNY,fNZ);
	initialize_sign_cells(fNX,fNY,fNZ);

	initialize_pressure_solver();
	
	if(USE_MICROPHYSICS){ init_microphysics(fNX,fNY);}
	
	init_boundaries(iebuffer,iwbuffer,jnbuffer,jsbuffer,3);

	initialize_budgets();
	
	init_damping(fNX,fNY,NZ);
	
	init_diffusion_weights(DIFFUSION_ORDER,&ZU(0));
	
}

/*********************************************************************
* Advance model foreward one full time step using third-order Runge-Kutta
* integration. One full step consists of three smaller steps.
*
**********************************************************************/
void p_integrate_rk3(){

	double steps[] = {1./3.,0.5,1.0};	// fractional time steps for RK3 loop
	int size = fNX*fNY*fNZ;

	/*******************************************************
	* Apply explicit diffusion. Do this first, because the
	* final state should satisfy the anelastic continuity
	* equation and be free of super saturation.
	********************************************************/
	if(USE_EXPLICIT_DIFFUSION){ apply_explicit_diffusion(1.0,3,fNX-3,3,fNY-3);}
	
	/*******************************************************
	* Calculate frictional and diffusional tendencies to
	* be applied during RK3 loop
	********************************************************/
	if(USE_TURBULENT_STRESS){ calculate_diff_tend(3,fNX-3,3,fNY-3);}
	
	/*******************************************************
	* Runge-Kutta Loop
	********************************************************/
	for(int s=0;s<3;s++){

		/*******************************************************
		* Exchange boundary values between processes. These values
		* are needed to interpolate to the cell faces of the control
		* volumes in order to calculate advective tendencies.
		********************************************************/
		exchange(us); exchange(vs); exchange(ws); exchange(ths);
		
		calculate_budgets(s,&steps[0]);
	
		// sign of advection for upwind biased derivatives
		compute_sign_cells(0,fNX,0,fNY);

		/*******************************************************
		* Solve momentum and pressure equations using either
		* the hydrostatic or non-hydrostatic equation set
		********************************************************/
		if(HYDROSTATIC){ integrate_hydro(steps[s],s,0,fNX,0,fNY);    } 
		else { 			 integrate_non_hydro(steps[s],0,fNX,0,fNY);  }
			
		/*******************************************************
		* Solve for new potential temperature field
		********************************************************/
		advect_theta(steps[s],3,fNX-3,3,fNY-3);
		
		if(USE_TURBULENT_STRESS){ apply_scalar_diffusion(steps[s],3,fNX-3,3,fNY-3);}
		
		if(USE_TERRAIN){ set_terrain(thps,istopos,size);}

		/*******************************************************
		* Advect microphysics variables
		********************************************************/
		if(USE_MICROPHYSICS){ 

			exchange(qvs); exchange(qcs); exchange(qrs);
			advect_microphysics_cell(steps[s],3,fNX-3,3,fNY-3);
			
			if(USE_ICE){
				exchange(qss); exchange(qis);
				advect_ice_cell(steps[s],3,fNX-3,3,fNY-3);
			}
			
			if(USE_TURBULENT_STRESS){ apply_moisture_diffusion(steps[s],3,fNX-3,3,fNY-3);}
			
			zero_moisture(3,fNX-3,3,fNY-3,size);
		}
		
		/*******************************************************
		* If using the hydrostatic equation set, calculate the
		* velocity for the new momentum field.
		********************************************************/
		if(HYDROSTATIC){ 
			
			exchange(ups); exchange(vps);
			w_velocity_LH(3,fNX-3,3,fNY-3);
		}
		
		/*******************************************************
		* Handle all boundary conditions for the full domain
		********************************************************/
		apply_boundary_condition(MPI_PROC_NULL);

		/*******************************************************
		* Handle all boundary conditions for the full domain
		* for microphysics variables and then advance them
		* forward for the next time step.
		********************************************************/
		if(USE_MICROPHYSICS){

			apply_boundary_condition_microphysics(MPI_PROC_NULL);

			if(s<2){ microphysics_advance_inner(fNX*fNY*fNZ*sizeof(double));}
		}

		if(s<2){ advance_inner(size*sizeof(double));}
	}

	/*****************************************************************
	* POST INTEGRATION. NOW RUN PHYSICAL PARAMETERIZATIONS
	******************************************************************/
	//damp_var(&THP(0,0,0),3,fNX-3,3,fNY-3,1.157e-6,2.3148e-5);
	//damp_var(&THP(0,0,0),3,fNX-3,3,fNY-3,2.3148e-5,2.3148e-5);

	/*********************************************
	* Apply heating perturbation
	**********************************************/
	if(bigcounter<stopheating){ heat.p_applyHeating();}
	
	/*********************************************
	* Handle microphysics
	**********************************************/
	if(USE_MICROPHYSICS){

		run_microphysics(3,fNX-3,3,fNY-3);
		
		if(USE_TERRAIN){
			set_terrain(qvps,istopos,size);
			set_terrain(qrps,istopos,size);
			set_terrain(qcps,istopos,size);
		}

		apply_boundary_condition_microphysics(MPI_PROC_NULL);

		microphysics_advance(size*sizeof(double));
	}

	/*********************************************
	* 
	**********************************************/
	if(PV_TRACER){ pv_tracer_sources(3,fNX-3,3,fNY-3);}

	/*********************************************
	* BOUNDARY CONDITIONS
	**********************************************/
	apply_boundary_condition(MPI_PROC_NULL);

	rayleigh_damping(0,fNX,0,fNY,raydampheight);

	advance_outer(size*sizeof(double));
}

/*********************************************************************
* Primary function for running parallel model.
*
**********************************************************************/
void p_run_model(int count,FILE *infile){

	//int counter = 0;	
	double elapsed = 0;
	clock_t start_time=clock(),finis_time=0;
	double start_walltime = MPI_Wtime(),elapsed_walltime;
	double total_walltime = 0,total_cputime = 0;
	int timer_counter = 0;
	bool isFirstStep = true;
	
	clock_t start_time_full=clock();
	
	//---------------------------------------------------------------------------
	// Step through model 'count' number of times
	//---------------------------------------------------------------------------
	while(bigcounter<count){
		
		//-----------------------------------------------------------------------
		// OUTPUT TO FILE
		//-----------------------------------------------------------------------
		if(OUTPUT_TO_FILE && bigcounter % outfilefreq == 0 && (!isRestartRun || !isFirstStep) ){
			
			if(rank==0){
				
				file_output_status(MODEL_FILES_NOT_WRITTEN);
				write_time_to_file(filename,file_time_counter);
			}
			
			output_meteorological_fields_to_file(parallel_write_pvar_to_file_2d,
												 parallel_write_pvar_to_file_3d,
												 file_time_counter);
			
			if(rank==0){ file_output_status(MODEL_FILES_WRITTEN);}
			
			write_budgets_to_file();
			
			if(rank==0){ file_output_status(ALL_FILES_WRITTEN);}
		}

		p_integrate_rk3();	// run model forward one time step

		//-----------------------------------------------------------------------
		// Output pressure field and calculate time to completion
		//-----------------------------------------------------------------------
		if(OUTPUT_TO_FILE && bigcounter % outfilefreq == 0){
			
			if(!isRestartRun || !isFirstStep)
				parallel_write_pvar_to_file_3d(filename,"pi",pis,file_time_counter);
			
			file_time_counter++;
			
			if(VERBOSE && !isFirstStep && rank==0){
				print_time_estimates(total_cputime,total_walltime,timer_counter);
			}
			
			total_walltime = 0;
			total_cputime = 0;
			timer_counter = 0;
		}

		//-----------------------------------------------------------------------
		// TIME EACH TIME STEP
		//-----------------------------------------------------------------------
		if(rank==0){
	
			finis_time = clock();
	
			elapsed = ((double) (finis_time - start_time)) / CLOCKS_PER_SEC;
			elapsed_walltime = MPI_Wtime() - start_walltime;

			total_walltime += elapsed_walltime;
			total_cputime += elapsed;
			timer_counter += 1;
				

			if(VERBOSE){ printf("time %0.3f hr %0.3f s %0.3f s\n",mtime/3600,elapsed,elapsed_walltime);}
			fflush(stdout);
	
			start_time = clock();
			start_walltime = MPI_Wtime();
		}

		mtime += dt;	// elapsed physical time
		
		bigcounter++;	// total time steps (if this function is called multiple times)
		//counter++;		// time steps within this loop
		isFirstStep = false;
	}
	
	if(rank==0){
		finis_time = clock();
		elapsed = ((double) (finis_time - start_time_full)) / CLOCKS_PER_SEC;

		printf("Program successfully completed. Total runtime: %f s\n",elapsed);
	}
}

/*********************************************************************
* Run parallel model
*
**********************************************************************/
void run_parallel_model(int argc, char *argv[]){

	initialize_parallel_environment(argc,argv);
	
	initialize_parallel_model();

	p_run_model(number_of_time_steps);

	MPI_Finalize();
	
}
