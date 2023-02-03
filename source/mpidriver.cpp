#include "stdafx.h"
#include "Heating.h"
#include "advection.h"
#include "energy.h"
#include "mpi.h"
#include "fluxes.h"
#include "turbulence.h"
#include "budgets.h"
#include "damping.h"
#include "initializer.h"
#include "boundaries.h"
#include "pressure.h"

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

#define BCAST_NX(var) MPI_Bcast(&var[0],NX,MPI_DOUBLE,0,MPI_COMM_WORLD)
#define BCAST_NY(var) MPI_Bcast(&var[0],NY,MPI_DOUBLE,0,MPI_COMM_WORLD)
#define BCAST_NZ(var) MPI_Bcast(&var[0],NZ,MPI_DOUBLE,0,MPI_COMM_WORLD)

MPI_Comm comm_cart,comm_row,comm_col;

int myNX,myNY,myNZ;	// dimensions of subarrays without halo boundaries
int fNX,fNY,fNZ;	// full dimensions of subarrays
int pNX,pNY,pNZ;	// dimensions of arrays for pressure solver
int fNYfNZ;

int coords[2];

// ranks of neighboring processes
int west,east,south,north;
int nw_corner=MPI_PROC_NULL;
int ne_corner=MPI_PROC_NULL;
int sw_corner=MPI_PROC_NULL;
int se_corner=MPI_PROC_NULL;

int error_status = 0;

int rank,row_rank,col_rank;
int dims[] = {0,0};
int numtasks;
int row_size;
int col_size;

// index within the larger domain
int *big_i,*big_j;

// sizes of subarrays on each process
int *s_nx,*s_ny,*s_nz,*s_ny_p,*s_nx_p;

// starting and ending indices of each process within the full domain arrays
int *ibs,*ibe,*jbs,*jbe,*kbs,*kbe,*jbs_p,*ibs_p;

int *scounts,*rcounts,*sdisps,*rdisps;
MPI_Datatype *types;

double *pres_row,*ipres_row;
double *pres_col,*ipres_col;

double *frictions;
double *istopos,*uistopos,*vistopos;

int size,size2;


void p_integrate_hydro(double,int,int,int,int,int);
void p_run_model(int,FILE *infile=NULL);

/*********************************************************************
* Based on the number of processors along a particular axis, calculate
* how many grid cells each processor gets
*
* @param nproc number of processors along axis
* @param nd number of grid points along axis
* @param nout[nproc] calculated number of grid points for each processor
**********************************************************************/
void decompose_axis(int nproc,int nd,int * nout){

	int sub_nx;

	int nmin = nd/nproc;	// minimum number of grid points per task

	int extra_rows = nd % nproc;	// extra rows

	for(int i=0;i<nproc;i++){

		sub_nx = (i+1 <= extra_rows) ? nmin+1 : nmin;
		nout[i] = sub_nx;
	}
}

/*********************************************************************
* Calculate the starting ij-coordinates within the full domain for
* each process.
*
* @param xproc,yproc - number of processes in each direction
* @param s_nx,s_ny - arrays of sizes xproc,yproc of gridpoints
* @param ibs,jbs - output
**********************************************************************/
void get_sub_dims(int xproc,int yproc,int *s_nx,int *s_ny,int *ibs,int *jbs){

	int xcount = 0, ycount = 0;

	for(int i=0;i<xproc-1;i++){
		xcount+=s_nx[i];
		for(int j=0;j<yproc;j++){ ibs[yproc*(i+1) + j] = xcount;}
	}

	for(int i=0;i<xproc;i++){
		ycount = 0;
		for(int j=0;j<yproc-1;j++){
			ycount+=s_ny[j];
			jbs[yproc*(i) + j+1] = ycount;
		}
	}
}

/*********************************************************************
* Setup grid information specific to anelastic pressure solver
*
* @param xproc,yproc - number of processes in each direction
* @param s_nz - output for grid sizes for each processor
* @param kbs - output for beginning k index for each processor
**********************************************************************/
void pressure_dims(int xproc,int yproc,int nx,int ny,int nz){

	if(nz > yproc){
		
		s_nz = (int *)calloc(yproc,sizeof(int));
		kbs = (int *)calloc(yproc,sizeof(int));
		
		s_ny_p = (int *)calloc(xproc,sizeof(int));
		jbs_p = (int *)calloc(xproc,sizeof(int));
		
		s_nx_p = (int *)calloc(yproc,sizeof(int));
		ibs_p = (int *)calloc(yproc,sizeof(int));
		
		decompose_axis(yproc,nz-2,s_nz);
		decompose_axis(xproc,ny,s_ny_p);
		decompose_axis(yproc,nx,s_nx_p);
		
	} else {
		
		printf("Error! Vertical levels must be smaller than yproc\n");
		exit(0);
	}
	
	kbs[0] = 1;
	jbs_p[0] = 0;
	ibs_p[0] = 0;
	
	for(int k=1;k<yproc;k++){ kbs[k] = kbs[k-1]+s_nz[k-1];}//	printf("kbs %d %d\n",k,kbs[k]);}
	for(int j=1;j<xproc;j++){ jbs_p[j]=jbs_p[j-1]+s_ny_p[j-1];}//	printf("jbs %d %d\n",j,jbs[j]);}
	for(int i=1;i<yproc;i++){ ibs_p[i]=ibs_p[i-1]+s_nx_p[i-1];}//	printf("ibs %d %d\n",i,ibs[i]);}
	
	pNX = s_nx_p[col_rank];
	pNY = s_ny_p[row_rank];
	pNZ = s_nz[col_rank];
	
	
	//printf("pnx %d pny %d pnz %d\n",pNX,pNY,pNZ);
}

/*********************************************************************
*
*
**********************************************************************/
void setup_grid(){

	/**************************************
	* Set up Cartesian grid
	***************************************/
	int periods[] = {PERIODIC_BOUNDARIES,PERIODIC_BOUNDARIES_NS};

	MPI_Dims_create(numtasks,2,dims);

	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &comm_cart);

	MPI_Comm_rank(comm_cart,&rank);

	MPI_Cart_coords(comm_cart,rank,2,coords);

	/**************************************
	* Memory allocation
	***************************************/
	s_nx = (int *)calloc(dims[0],sizeof(int));
	s_ny = (int *)calloc(dims[1],sizeof(int));

	ibs = (int *)calloc(numtasks,sizeof(int));
	ibe = (int *)calloc(numtasks,sizeof(int));
	jbs = (int *)calloc(numtasks,sizeof(int));
	jbe = (int *)calloc(numtasks,sizeof(int));

	scounts = (int *)calloc(numtasks,sizeof(int));
	rcounts = (int *)calloc(numtasks,sizeof(int));
	sdisps = (int *)calloc(numtasks,sizeof(int));
	rdisps = (int *)calloc(numtasks,sizeof(int));

	types = (MPI_Datatype *)malloc(numtasks*sizeof(MPI_Datatype));

	/**************************************
	* Calculate grid sizes for each process
	***************************************/
	decompose_axis(dims[0],NX,s_nx);
	decompose_axis(dims[1],NY,s_ny);

	myNX = s_nx[coords[0]];
	myNY = s_ny[coords[1]];
	myNZ = NZ;

	fNX = s_nx[coords[0]]+2*halo_buffer;
	fNY = s_ny[coords[1]]+2*halo_buffer;
	fNZ = NZ;

	fNYfNZ = fNY*fNZ;

	size = myNX*myNY*myNZ;
	size2 = fNX*fNY*myNZ;

	get_sub_dims(dims[0],dims[1],s_nx,s_ny,ibs,jbs);

	if(VERBOSE){ printf("Process %d \t has starting i index \t %d, \t starting j index \t %d, \t and grid dimensions \t %d x %d\n",rank,ibs[rank],jbs[rank],myNX,myNY);}

	/**************************************
	* Calculate index within big domain
	***************************************/	
	big_i = (int *)calloc(fNX,sizeof(int));
	big_j = (int *)calloc(fNY,sizeof(int));
	
	for(int i=0;i<fNX;i++){ big_i[i] = i + ibs[rank] - halo_buffer;}//printf("%d ",big_i[i]);}
	//printf("\n");
	for(int j=0;j<fNY;j++){ big_j[j] = j + jbs[rank] - halo_buffer;}//printf("%d ",big_j[j]);}
	//printf("\n");
	//exit(0);
	/**************************************
	* Determine who my neighboring processes are
	***************************************/
	MPI_Cart_shift ( comm_cart, 0 , 1 , & west , & east );
	MPI_Cart_shift ( comm_cart, 1 , 1 , & south , & north );

	get_corner(1,1,coords,dims,&ne_corner);
	get_corner(-1,-1,coords,dims,&sw_corner);
	get_corner(1,-1,coords,dims,&se_corner);
	get_corner(-1,1,coords,dims,&nw_corner);

	//printf("Rank %d corners are sw %d se %d nw %d ne %d\n",rank,sw_corner,se_corner,nw_corner,ne_corner);

	/**************************************
	* Create row and column communicators
	***************************************/
	int world_size;
	int colcolor = coords[0];
	int rowcolor = coords[1];
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	MPI_Comm_split(MPI_COMM_WORLD, rowcolor, rank, &comm_row);
	MPI_Comm_split(MPI_COMM_WORLD, colcolor, rank, &comm_col);

	MPI_Comm_rank(comm_row, &row_rank);
	MPI_Comm_size(comm_row, &row_size);

	MPI_Comm_rank(comm_col, &col_rank);
	MPI_Comm_size(comm_col, &col_size);

	if(!HYDROSTATIC){ pressure_dims(dims[0],dims[1],NX,NY,NZ);}
}

/*********************************************************************
* Data which is initialized on the root process is sent to all of
* of the other processes.
**********************************************************************/
void broadcast_shared_data(){

	//----------------------------------------------
	// This data is the same on each process
	//----------------------------------------------
	BCAST_NZ(rhou); BCAST_NZ(rhow);
	BCAST_NZ(zu); BCAST_NZ(zw);
	BCAST_NZ(tb); BCAST_NZ(tbw); BCAST_NZ(tbv);
	BCAST_NZ(qb);
	BCAST_NZ(pib);
	BCAST_NZ(one_d_rhou);
	BCAST_NZ(one_d_rhow);

	BCAST_NZ(zsu); BCAST_NZ(zsw);
	BCAST_NZ(mu); BCAST_NZ(mw);
	//BCAST_NZ(kmixv_moisture);

	BCAST_NY(f); BCAST_NY(dfdy);
	BCAST_NY(outLats); BCAST_NX(outLons);

	//MPI_Bcast(landsea,NX*NY,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&RHOAVG2DFULL(0,0),NX*NY,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&HTOPOFULL(0,0),NX*NY,MPI_INT,0,MPI_COMM_WORLD);

	//----------------------------------------------
	// This data is different on each process
	//----------------------------------------------
	distributeArray_3d(m_ubar,iubar);
	distributeArray_3d(m_vbar,ivbar);
	distributeArray_3d(m_wbar,iwbar);
	distributeArray_3d(m_thbar,ithbar);
	distributeArray_3d(m_qbar,iqbar);
	distributeArray_3d(m_pbar,ipbar);
	distributeArray_3d(istopos,iistopo);
	distributeArray_3d(uistopos,iuistopo);
	distributeArray_3d(vistopos,ivistopo);
	distributeArray_3d(frictions,ifriction);

	//----------------------------------------------
	// Fill the ghost cells and lateral boundaries
	//----------------------------------------------
	exchange(m_ubar); exchange(m_vbar); exchange(m_wbar); exchange(m_thbar);  exchange(m_qbar);	exchange(m_pbar);
	exchange(istopos); exchange(uistopos); exchange(vistopos);
	exchange(frictions);

	p_mirror_boundaries(m_ubar,east,west,north,south,halo_buffer,fNX,fNY);
	p_mirror_boundaries(m_vbar,east,west,north,south,halo_buffer,fNX,fNY);
	p_mirror_boundaries(m_thbar,east,west,north,south,halo_buffer,fNX,fNY);
	p_mirror_boundaries(m_qbar,east,west,north,south,halo_buffer,fNX,fNY);
	p_mirror_boundaries(m_pbar,east,west,north,south,halo_buffer,fNX,fNY);
	p_mirror_boundaries(istopos,east,west,north,south,halo_buffer,fNX,fNY);
	p_mirror_boundaries(uistopos,east,west,north,south,halo_buffer,fNX,fNY);
	p_mirror_boundaries(vistopos,east,west,north,south,halo_buffer,fNX,fNY);
	p_mirror_boundaries(frictions,east,west,north,south,halo_buffer,fNX,fNY);
	
}

/*********************************************************************
* Broadcast data
*
* @param size - number of elements
**********************************************************************/
void broadcast_data_int(int size,int *var){
	
	MPI_Bcast(var,size,MPI_INT,0,MPI_COMM_WORLD);
}

/*********************************************************************
* Broadcast data
*
* @param size - number of elements
**********************************************************************/
void broadcast_data_bool(int size,int *var){
	
	MPI_Bcast(var,size,MPI_INT,0,MPI_COMM_WORLD);
}

/*********************************************************************
* 
*
**********************************************************************/
void exit_from_error(){
	
	error_status = -1;
	
	if(!PARALLEL || rank ==0){
	
		printf("Programm terminating");
	
	}
	
	if(PARALLEL){ MPI_Finalize();}
	
	exit(0);
}

/*********************************************************************
* 
*
**********************************************************************/
void print_time_estimates(int total_cputime,int total_walltime,int timer_counter){

	int cputime_secs = (int)((total_cputime  / (double)timer_counter) * (number_of_time_steps - bigcounter));
	int walltime_secs = (int)((total_walltime / (double)timer_counter) * (number_of_time_steps - bigcounter));

	int cpu_time_hours = cputime_secs / 3600;
	int cpu_time_min = (cputime_secs % 3600) / 60;
	int cpu_time_sec = (cputime_secs % 3600) % 60;
	int wall_time_hours = walltime_secs / 3600;
	int wall_time_min = (walltime_secs % 3600) / 60;
	int wall_time_sec = (walltime_secs % 3600) % 60;

	printf("Estimated time remaining: cpu time = %02d:%02d:%02d hours, wall time = %02d:%02d:%02d hours\n",
	cpu_time_hours,cpu_time_min,cpu_time_sec,wall_time_hours,wall_time_min,wall_time_sec);

	cputime_secs = (int)((total_cputime  / (double)timer_counter) * number_of_time_steps);
	walltime_secs = (int)((total_walltime / (double)timer_counter) * number_of_time_steps);
			
	cpu_time_hours = cputime_secs / 3600;
	cpu_time_min = (cputime_secs % 3600) / 60;
	wall_time_hours = walltime_secs / 3600;
	wall_time_min = (walltime_secs % 3600) / 60;
	cpu_time_sec = (cputime_secs % 3600) % 60;
	wall_time_sec = (walltime_secs % 3600) % 60;
		
	printf("Estimated total run time: cpu time = %02d:%02d:%02d hours, wall time = %02d:%02d:%02d hours\n",
	cpu_time_hours,cpu_time_min,cpu_time_sec,wall_time_hours,wall_time_min,wall_time_sec);
	
}
#if 0
/*********************************************************************
* 
**********************************************************************/
void output_meteorological_fields_to_file(){
	
	parallel_write_pvar_to_file_3d(filename,"u-wind",us, file_time_counter);
	parallel_write_pvar_to_file_3d(filename,"v-wind",vs, file_time_counter);
	parallel_write_pvar_to_file_3d(filename,"w-wind",ws, file_time_counter);
	parallel_write_pvar_to_file_3d(filename,"pi",    pis,file_time_counter);
	parallel_write_pvar_to_file_3d(filename,"theta", ths,file_time_counter);

	if(USE_MICROPHYSICS){
		
		parallel_write_pvar_to_file_3d(filename,"qv",qvs,file_time_counter);
		parallel_write_pvar_to_file_3d(filename,"qc",qcs,file_time_counter);
		parallel_write_pvar_to_file_3d(filename,"qr",qrs,file_time_counter);
		
		parallel_write_pvar_to_file_2d(filename,"rainfall",accRain,file_time_counter);
		
		if(USE_ICE){
			
			parallel_write_pvar_to_file_2d(filename,"snowfall",accSnow,file_time_counter);
			
			parallel_write_pvar_to_file_3d(filename,"qi",qis,file_time_counter);
			parallel_write_pvar_to_file_3d(filename,"qs",qss,file_time_counter);
		}
	}
	
	if(USE_TURBULENT_STRESS){
		
		parallel_write_pvar_to_file_2d(filename,"int_fric",integrated_friction_ke,file_time_counter);
		
		if(SURFACE_HEAT_FLUX){
			parallel_write_pvar_to_file_2d(filename,"lh_flux",latent_heat_flux,file_time_counter);
		}
	}
	
	if(EXTRA_OUTPUT){
		parallel_write_pvar_to_file_3d(filename,"rate",rate,file_time_counter);
	}

	if(OUTPUT_FRICTION_TEND){
		parallel_write_pvar_to_file_3d(filename,"fric",frictions,file_time_counter);
	}
}
#endif

/*********************************************************************
* Advance model foreward one full time step using third-order Runge-Kutta
* integration. One full step consists of three smaller steps.
*
**********************************************************************/
void p_integrate_rk3(){

	double steps[] = {1./3.,0.5,1.0};	// fractional time steps for RK3 loop
	int size = fNX*fNY*fNZ;
	
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
}

/*********************************************************************
* Top-level initializer for parallel model
* Setup parallel environment.
*
* argc - number of command-line arguments
* argv[] - command-line arguments
**********************************************************************/
void initialize_parallel_model(int argc, char *argv[]){
	
	int len; 
	char hostname[MPI_MAX_PROCESSOR_NAME];
	//------------------------------------------------
	// Set up parallel environment
	//------------------------------------------------
	MPI_Init(&argc,&argv);					// initialize MPI  
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);// get number of tasks 
	MPI_Get_processor_name(hostname, &len);
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

	initialize_pressure_solver();
	
	if(USE_MICROPHYSICS){ init_microphysics(fNX,fNY);}
	
	init_boundaries(iebuffer,iwbuffer,jnbuffer,jsbuffer,3);

	initialize_budgets();
	
}

/*********************************************************************
* Run parallel model
*
**********************************************************************/
void run_parallel_model(int argc, char *argv[]){
	
	initialize_parallel_model(argc,argv);

	p_run_model(number_of_time_steps);

	MPI_Finalize();
	
}
