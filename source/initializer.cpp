#include "stdafx.h"
#include "interpolate.h"
#include "surface.h"
#include "fluxes.h"
#include "damping.h"
#include "boundaries.h"
#include "initializer.h"
#include "data_initializer.h"
#include "pressure.h"

/******************************************************************************
* Handles higher-level initialization and focuses on initialization from model
* output files. Also contains the smaller initialization routines.
*
*******************************************************************************/

#define d(i,j,k) (xdim*ydim*(k) + xdim*(j) + i)
#define ALLOC(v,size) v = (double *)calloc(size,sizeof(double))

/***************************************************************************
* -------------------------- VARIABLES -------------------------------------
****************************************************************************/
double *output_to_file;//u[NX][NY][NZ];
//--------------------------------------
// BASIC STATE ARRAYS
//--------------------------------------
double *m_ubar,*m_vbar,*m_thbar,*m_wbar,*m_qbar,*m_pbar;
double *iubar,*ivbar,*iwbar,*ithbar,*iqbar,*ipbar;
//--------------------------------------
// PERTURBATION ARRAYS
//--------------------------------------
double *us,*vs,*ths,*pis,*ws;
double *ups,*vps,*wps,*thps;
double *ums,*vms,*wms,*thms;
double *qvps,*qvs,*qvms;
double *qcps,*qcs,*qcms;
double *qrps,*qrs,*qrms;
double *qips,*qis,*qims;
double *qsps,*qss,*qsms;
double *rate;
//--------------------------------------
// FALL SPEEDS
//--------------------------------------
double *vts,*sts,*its;
//--------------------------------------
// TOPOGRAPHY AND FRICTION ARRAYS
//--------------------------------------
double *itopo,*iistopo,*iuistopo,*ivistopo,*ifriction;
//--------------------------------------
// VERTICALLY VARYING BASIC STATE
//--------------------------------------
double *zu,*zw,*rhou,*rhow,*tb,*tbw,*tbv,*qb,*pib;
double *one_d_rhou,*one_d_rhow;
//--------------------------------------
// VERTICAL COORDINATE PARAMETERS
//--------------------------------------
double *zsu,*zsw;
double *mu,*mw;
//--------------------------------------
// OTHER ARRAYS
//--------------------------------------
//double f[NY];
double *f;
double *dfdy;

double *rhoavg2d;//[NX][NY];
double *outLats;
double *outLons;
int *htopo;//[NX][NY];

int linear_lon;
double rhoavg = 0;
double mtime = 0;
int bigcounter = 0;

int NX,NY,NZ,NYNZ;// = 1035/3;//173;//345;//173;//345;//518;
// = 777/3;//111;//222;//111;//222;

double dx;// = 15000.0;
double dy;// = 15000.0;
double dz;// = 500.0;
double dt;// = 60.0;

double dtx;
double dty;
double dtz;

double one_d_dx;
double one_d_dy;
double one_d_dz;

double lonoffset;
double latoffset;
int number_of_time_steps;
int raydampheight;
int outfilefreq;
double height_lowest_level;

int PV_BUDGET;
int HEAT_BUDGET;
int MOISTURE_BUDGET;
int VORTICITY_BUDGET;
int PE_BUDGET;

int MICROPHYSICS_OPTION;
int USE_ICE;
int USE_MICROPHYSICS;
int isRestartRun;
int SURFACE_HEAT_FLUX;

int VERBOSE;

int BASIC_STATE_OPTION;			// 0:Model output, 1:Code, 2:Reanalysis 
int PERTURBATION_OPTION;
int CREATE_NEW_OUTPUT_FILE;

int perturbationFileTime;
int startOutfileAt;

/***************************************************************************
* ---------------------- FUNCTION PROTOTYPES--------------------------------
****************************************************************************/
void initialize_from_output_serial(const char *,size_t);
void initialize_from_output_parallel(const char *,size_t);
void initialize_basic_state_from_output_file(const char*);
int get_file_start_time();

/******************************************************************************
* Takes the structure that has stored the input from the input file with grid
* specifications and initializes the appropriate variables.
*
*******************************************************************************/
void initialize_globals(){
	
	NX = inputs.nx;
	NY = inputs.ny;
	NZ = inputs.nz;
	
	NYNZ = NY*NZ;
	
	dx = inputs.dx;
	dy = inputs.dy;
	dz = inputs.dz;
	dt = inputs.dt;
	
	dtx = dt/dx;
	dty = dt/dy;
	dtz = dt/dz;

	one_d_dx = 1.0/dx;
	one_d_dy = 1.0/dy;
	one_d_dz = 1.0/dz;
	
	lonoffset = inputs.corner_lon;
	latoffset = inputs.corner_lat;
	number_of_time_steps = inputs.time_steps;
	raydampheight = inputs.rayleigh_damping_z;
	outfilefreq = inputs.output_frequency;
	height_lowest_level = inputs.height_lowest_level;
	
	isRestartRun = inputs.is_restart_run;

	BASIC_STATE_OPTION = inputs.basic_state_init_option;
	PERTURBATION_OPTION = inputs.perturbation_init_option;

	CREATE_NEW_OUTPUT_FILE = inputs.create_new_output_file;

	perturbationFileTime = inputs.perturbationFileTime;
	startOutfileAt 		 = inputs.startOutfileAt;

	PV_BUDGET 			= 	inputs.run_pv_budget;
	HEAT_BUDGET 		= 	inputs.run_heat_budget;
	MOISTURE_BUDGET 	= 	inputs.run_moisture_budget;
	VORTICITY_BUDGET	= 	inputs.run_vorticity_budget;
	PE_BUDGET 			= 	inputs.run_pe_budget;
	
	VERBOSE = inputs.verbose;
	
	MICROPHYSICS_OPTION = inputs.microphysics_option;
	
	if(MICROPHYSICS_OPTION==2){	USE_ICE = 1;} else { USE_ICE = 0;}
	if(MICROPHYSICS_OPTION > 0){ USE_MICROPHYSICS = 1;} else {USE_MICROPHYSICS = 0;}
	
	if(USE_MICROPHYSICS == 0 && MOISTURE_BUDGET == 1){ 
		MOISTURE_BUDGET = 0;
		printf("Warning: can't use moisture budget since microphysics is not active.\n");
	}
	
	if(CREATE_NEW_OUTPUT_FILE==1 && PERTURBATION_OPTION==0 && 
		(strcmp(inputs.perturbation_file,inputs.output_file) == 0 || 
		 strcmp(inputs.perturbation_file,inputs.output_file) == 0)
	){ 
		printf("Error: output file must have a name different from input files!\n");
		exit(0);

	}
	
	SURFACE_HEAT_FLUX = inputs.use_surface_heat_flux;
	
	//printf("options = %d %d %d\n",MICROPHYSICS_OPTION,USE_MICROPHYSICS,USE_ICE);
}

/******************************************************************************
*
*
*******************************************************************************/
void initialize_basic_state(){
	
	int my_basic_state_option = BASIC_STATE_OPTION;
	
	if(isRestartRun){ my_basic_state_option = 0;}
	
	switch(my_basic_state_option){	
		//----------------------------------------------------------------
		// Initialize from a model output file
		//----------------------------------------------------------------		
		case 0:
			printf("---------------------------------------------\n");
			printf("Loading basic state from a model output field\n");
			printf("---------------------------------------------\n");
		
			if(isRestartRun){
				initialize_basic_state_from_output_file(filename);
			} else {
				initialize_basic_state_from_output_file(inputBasicStateFileName);
			}
		
			break;
		//----------------------------------------------------------------
		// Initialize from code
		//----------------------------------------------------------------
		case 1:
			printf("------------------------------------------------\n");
			printf("Initializing basic state from a code subroutine \n");
			printf("------------------------------------------------\n");
					
			initialize_basic_state_from_subroutine();
			
			break;
		//----------------------------------------------------------------
		// Initialize from a reanalysis dataset (or user created file)
		//----------------------------------------------------------------
		case 2:
			printf("---------------------------------------------\n");
			printf("Loading basic state from a reanalysis\n");
			printf("---------------------------------------------\n");

			initialize_basic_state_from_reanalysis();
	
			break;		
		//----------------------------------------------------------------
		// No option selected! Exit program.
		//----------------------------------------------------------------
		default:
			printf("---------------------------------------------\n");		
			printf("Option %d in basic state initializer not supported!\n",BASIC_STATE_OPTION);
			printf("---------------------------------------------\n");
			
			exit(0);
		
			break;
	}
	
}

/******************************************************************************
*
*
*******************************************************************************/
void initialize_perturbation(){
	
	int my_perturbation_option = PERTURBATION_OPTION;
	int startTime = 0;
	
	if(isRestartRun){ my_perturbation_option = 0;}
	
	switch(my_perturbation_option){
		//----------------------------------------------------------------
		// Initialize from a model output file
		//----------------------------------------------------------------		
		case 0:
		
			if(!PARALLEL || rank==0){
				printf("----------------------------------------------\n");
				printf("Loading perturbation from a model output field\n");
				printf("----------------------------------------------\n");
				
				//--------------------------------------------------------
				// Handle start time within input file
				//--------------------------------------------------------
				startTime = get_file_start_time();

				printf("Initialized from time = %d\n",startTime);
			}
			
			if(PARALLEL){
				if(isRestartRun){
					initialize_from_output_parallel(filename,startTime);
				} else {	
					initialize_from_output_parallel(inputPerturbationFileName,startTime);
				}
			} else {
				if(isRestartRun){
					initialize_from_output_serial(filename,startTime);
				} else {
					initialize_from_output_serial(inputPerturbationFileName,startTime);
				}	
			}
		
			break;
		//----------------------------------------------------------------
		// Initialize from code
		//----------------------------------------------------------------
		case 1:
		
			if(!PARALLEL || rank==0){		
				printf("------------------------------------------------\n");
				printf("Initializing perturbation from a code subroutine\n");
				printf("------------------------------------------------\n");
			}		
					
			initialize_perturbation_from_subroutine();
			
			break;
		//----------------------------------------------------------------
		// No option selected! Exit program.
		//----------------------------------------------------------------
		default:
			if(!PARALLEL || rank==0){
				printf("----------------------------------------------------\n");		
				printf("Option %d in perturbation initializer not supported!\n",BASIC_STATE_OPTION);
				printf("----------------------------------------------------\n");
			}
			
			exit(0);
		
			break;
	}
}


/*********************************************************************
* 
**********************************************************************/
void initialize_1D_arrays(){

	outLats = (double*) calloc(NY,sizeof(double));
	outLons = (double*) calloc(NX,sizeof(double));
	f 		= (double*) calloc(NY,sizeof(double));
	dfdy 	= (double*) calloc(NY,sizeof(double));
	
	htopo 		= (int*) calloc(NX*NY,sizeof(int));
	rhoavg2d 	= (double*) calloc(NX*NY,sizeof(double));
	
	zu = (double*) calloc(NZ,sizeof(double));
	zw = (double*) calloc(NZ,sizeof(double));
	rhou = (double*) calloc(NZ,sizeof(double));
	rhow = (double*) calloc(NZ,sizeof(double));
	tb = (double*) calloc(NZ,sizeof(double));
	tbw = (double*) calloc(NZ,sizeof(double));
	tbv = (double*) calloc(NZ,sizeof(double));
	qb = (double*) calloc(NZ,sizeof(double));
	pib = (double*) calloc(NZ,sizeof(double));
	one_d_rhou = (double*) calloc(NZ,sizeof(double));
	one_d_rhow = (double*) calloc(NZ,sizeof(double));
	zsu = (double*) calloc(NZ,sizeof(double));
	zsw = (double*) calloc(NZ,sizeof(double));
	mu = (double*) calloc(NZ,sizeof(double));
	mw = (double*) calloc(NZ,sizeof(double));
}

/*********************************************************************
* Initialize all subarrays for each process.
*
* @param size - number of elements to allocate
**********************************************************************/
void initialize_subarray(int size){
	
	// this was already allocated on the root process, just
	// do it on the other processes if parallel version
	if(PARALLEL && rank != 0){ initialize_1D_arrays();}

	//--------------------------------------------------
	// basic state
	//--------------------------------------------------
	ALLOC(m_ubar,size); ALLOC(m_vbar,size); ALLOC(m_wbar,size);
	ALLOC(m_qbar,size); ALLOC(m_thbar,size); ALLOC(m_pbar,size);

	//--------------------------------------------------
	// perturbations
	//--------------------------------------------------
	ALLOC(us,size);  ALLOC(vs,size);  ALLOC(ws,size);  ALLOC(ths,size);
	ALLOC(ups,size); ALLOC(vps,size); ALLOC(wps,size); ALLOC(thps,size);
	ALLOC(ums,size); ALLOC(vms,size); ALLOC(wms,size); ALLOC(thms,size);

	ALLOC(pis,size);

	//--------------------------------------------------
	// microphysics
	//--------------------------------------------------
	if(USE_MICROPHYSICS){

		ALLOC(qvs,size);  ALLOC(qcs,size);  ALLOC(qrs,size);
		ALLOC(qvps,size); ALLOC(qcps,size); ALLOC(qrps,size);
		ALLOC(qvms,size); ALLOC(qcms,size); ALLOC(qrms,size);
		
		ALLOC(vts,size);
		
		if(USE_ICE){
			
			ALLOC(qss,size);  ALLOC(qis,size); 
			ALLOC(qsps,size); ALLOC(qips,size); 
			ALLOC(qsms,size); ALLOC(qims,size);
			
			ALLOC(sts,size); ALLOC(its,size);
		}
	}

	ALLOC(frictions,size);
	ALLOC(istopos,size); ALLOC(uistopos,size); ALLOC(vistopos,size);

	if(EXTRA_OUTPUT){ ALLOC(rate,size);}

	//--------------------------------------------------
	// pressure
	//--------------------------------------------------
	if(HYDROSTATIC){
		ALLOC(pres_col,NY*myNX); ALLOC(ipres_col,NY*myNX);
		ALLOC(pres_row,NX*myNY); ALLOC(ipres_row,NX*myNY);
	} else {
		ALLOC(pres_col,NY*myNX*NZ); ALLOC(ipres_col,NY*myNX*NZ);
		ALLOC(pres_row,NX*pNY*NZ);  ALLOC(ipres_row,NX*pNY*NZ);
	}

}


/******************************************************************************
* Allocates a subset of memory on the root process for the parallel version.
* Allocates most of the used memory for the serial version.
*
*******************************************************************************/
void setup_memory_allocation(){

	int size = NX*NY*NZ;
	
	//----------------------------------------------------------------
	// initialize common data
	//----------------------------------------------------------------
	//if(PARALLEL){ memset(&u[0][0][0],0,NX*NY*NZ*sizeof(double));}

	initialize_1D_arrays();

	itopo = (double*) calloc(size,sizeof(double));

	//----------------------------------------------------------------
	// initialize data specific to parallel version
	//----------------------------------------------------------------
	if(PARALLEL || ENERGY){
	
		if(rank==0){ ALLOC(output_to_file,size);}
	
		ALLOC(iubar,size);
		ALLOC(ivbar,size);
		ALLOC(iwbar,size);
		ALLOC(ithbar,size);
		ALLOC(ipbar,size);
		ALLOC(iqbar,size);
	
		ALLOC(iuistopo,size);
		ALLOC(ivistopo,size);
		ALLOC(iistopo,size);
		ALLOC(ifriction,size);
			
		if(ENERGY){
			
			ALLOC(us,size);
			ALLOC(vs,size);
			ALLOC(ws,size);
			ALLOC(ths,size);
			ALLOC(pis,size);

			//--------------------------------------------------
			// microphysics
			//--------------------------------------------------
			if(USE_MICROPHYSICS){

				ALLOC(qvs,size);
				ALLOC(qcs,size);
				ALLOC(qrs,size);
			}
		}
		
		if(ENERGY){ if(OUTPUT_DIFFUSION_TEND){ init_damping(NX,NY,NZ); }}
		
	//----------------------------------------------------------------
	// initialize data specific to serial version
	//----------------------------------------------------------------
	} else {
		initialize_flux_cells(NY,NZ);
		initialize_microphysics_cells(NY,NZ);

		initialize_subarray(NX*NY*NZ);

		init_damping(NX,NY,NZ);
	}
		
		
	//----------------------------------------------------------------
	// Initialize stuff for Fourier damping in linearized equations
	//----------------------------------------------------------------
	if(ISLINEAR && FOURIER_DAMPING){
	
		if(!PERIODIC_BOUNDARIES){ init_fftw(NX,NY,NZ);}
		else { init_fftw(NX-6,NY,NZ);}
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
	
	if(USE_MICROPHYSICS){ init_microphysics();}
	
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
	
}

/******************************************************************************
* Initializes the basic state from various reanalysis dataset. Each dataset
* has a specific format (e.g. different variable names, units, array 
* structure, etc.) and requires a different initialization process.
* Support for additional datasets can be implemented elsewhere and called
* in this subroutine.
*
*******************************************************************************/
void initialize_basic_state_from_reanalysis(){
	
	//----------------------------------------------------------------
	// Allocates memory for variables
	//----------------------------------------------------------------
	setup_memory_allocation();
	//----------------------------------------------------------------
	// Calculate height of each level
	//----------------------------------------------------------------
	setup_vertical_height_levels();
	//----------------------------------------------------------------
	// Choose an initialization method for the basic state. Each 
	// function should read in data from the original file, interpolate
	// data to the staggered grid, and load the result into the basic 
	// state initialization arrays (i.e. IUBAR,IVBAR,ITHBAR, etc.)
	//----------------------------------------------------------------	
	switch(REANALYSIS_INITIALIZE_OPTION){
		//------------------------------------------------------------
		// Processed file from ERA-interim reanalysis, currently 
		// composite monsoon depression basic state
		//------------------------------------------------------------		
		case 1:
			initialize_from_era(&ZU(0));	
		break;
		//------------------------------------------------------------
		// NCEP reanalysis monthly means. Currently consist of
		// individual files for each field
		//------------------------------------------------------------			
		case 2:
			initialize_from_ncar();
		break;
		//------------------------------------------------------------
		// ERA5 monthly mean. Reads file directly.
		//------------------------------------------------------------	
		case 3:
			initialize_from_ERA5(&ZU(0));
		break;
		//------------------------------------------------------------
		// Processed file from ERA5 reanalysis, currently 
		// low-passed filtered for time favoring monsoon depressions
		//------------------------------------------------------------			
		case 4:
			initialize_from_ERA5_processed(&ZU(0));
		break;
		//----------------------------------------------------------------
		// No option selected! Exit program.
		//----------------------------------------------------------------
		default:
			printf("--------------------------------------------------\n");		
			printf("Option %d in reanalysis initializer not supported!\n",REANALYSIS_INITIALIZE_OPTION);
			printf("--------------------------------------------------\n");
			
			exit(0);
		
			break;
	}
		
	//----------------------------------------------------------------
	// If requested, create zonally uniform basic state
	//----------------------------------------------------------------
	if(MERIDIONAL_CROSS_SECTION){

		linear_lon = get_point_from_lon(meridional_lon_section);

		//printf("lon = %d\n",linear_lon);

		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
		for(int k=0;k<NZ;k++){
	
			IUBAR(i,j,k) = IUBAR(linear_lon,j,k);
			IVBAR(i,j,k) = 0;
			ITHBAR(i,j,k) = ITHBAR(linear_lon,j,k);
			IQBAR(i,j,k) = IQBAR(linear_lon,j,k);
			IPBAR(i,j,k) = IPBAR(linear_lon,j,k);
		}}}

		//----------------------------------------------------------------
		// For linearized model, create meridionally uniform conditions on the
		// northern and southern boundaries	
		//----------------------------------------------------------------
		if(ISLINEAR){

			int b1 = 0;
			int b2 = 10;
	
			for(int i=0;i<NX;i++){
			for(int k=0;k<NZ;k++){

				for(int j=0;j<b1;j++){
				
					IUBAR(i,j,k) = IUBAR(i,b1,k);
					ITHBAR(i,j,k) = ITHBAR(i,b1,k);
					IQBAR(i,j,k) = IQBAR(i,b1,k);
				}

				for(int j=NY-b2;j<NY;j++){

					IUBAR(i,j,k) = IUBAR(i,NY-b2,k);
					ITHBAR(i,j,k) = ITHBAR(i,NY-b2,k);
					IQBAR(i,j,k) = IQBAR(i,NY-b2,k);
				}
			}}
			
			
		}
	}
	
	//----------------------------------------------------------------
	// For resting basic state set the wind components to zero and make
	// the temperature and moisture fields horizontally uniform
	//----------------------------------------------------------------
	if(RESTING_BASIC_STATE){
		
		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
		for(int k=0;k<NZ;k++){

			IUBAR(i,j,k) = 0;
			IVBAR(i,j,k) = 0;
			ITHBAR(i,j,k) = ITHBAR(NX/2,NY/2,k);
			IQBAR(i,j,k) = IQBAR(NX/2,NY/2,k);
			IPBAR(i,j,k) = IPBAR(NX/2,NY/2,k);		
		}}}
	}
	
	//if(USE_MICROPHYSICS){ init_microphysics();}

	//----------------------------------------------------------------
	// Upper and lower boundary conditions for base state array
	//----------------------------------------------------------------
	upper_lower_boundaries(&ITHBAR(0,0,0));
	upper_lower_boundaries(&IQBAR(0,0,0));
	upper_lower_boundaries(&IUBAR(0,0,0));	
	upper_lower_boundaries(&IVBAR(0,0,0));

	//----------------------------------------------------------------
	// Initialize vertically varying, x,y independent basic state
	//----------------------------------------------------------------	
	initialize_vertical_basic_state(NX/2,NY/2);
	
	if(VERBOSE){ print_vertical_basic_state();}
	//----------------------------------------------------------------
	// Initialize topographic array
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
		//ITOPO(i,j,k) = 0;
	}}}
	
	init_topography();

	//----------------------------------------------------------------
	// Column averaged density
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	
		RHOAVG2DFULL(i,j) = 0;

		for(int k=HTOPOFULL(i,j)+1;k<NZ-1;k++){ RHOAVG2DFULL(i,j) = RHOAVG2DFULL(i,j) + rhou[k]*tbv[k];}

		RHOAVG2DFULL(i,j) = RHOAVG2DFULL(i,j)/((double)(NZ-2-HTOPOFULL(i,j)));
	}}

	//----------------------------------------------------------------
	// Initialize arrays for SOR
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		//sor[0][i][j] = 0;
		//sor[1][i][j] = 0;
		//sor[2][i][j] = 0;
		
	}}

	//----------------------------------------------------------------
	// Remove model base state from environmental base state
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				ITHBAR(i,j,k) = ITHBAR(i,j,k) - tb[k];
				IQBAR(i,j,k) = IQBAR(i,j,k) - qb[k];
			}
		}
	}
	
	//----------------------------------------------------------------
	// Initialize friction array for linear friction
	//----------------------------------------------------------------
	init_friction();

	//----------------------------------------------------------------
	// Calculate base state vertical velocity
	//----------------------------------------------------------------
	init_basic_state_vertical_velocity();

	//----------------------------------------------------------------
	// Coriolis parameter
	//----------------------------------------------------------------
	initialize_coriolis_parameter(outLats);
	
	//----------------------------------------------------------------
	// Alter humidty, if requested
	//----------------------------------------------------------------
	#if 0 && !ENERGY
	change_humidity();
	#endif
	

}

/******************************************************************************
* 
*******************************************************************************/
int get_file_start_time(){

	if(isRestartRun){
		
		int filestatus = get_file_output_status();
		
		if(filestatus == CAN_RESTART_FROM_FILE_END || filestatus == -1){
			return -1 + (int)get_file_length(filename,"t");
		} else {
			return -2 + (int)get_file_length(filename,"t");		
		}
	}
	
	int startTime = perturbationFileTime;
	int fileSize = (int)get_file_length(inputPerturbationFileName,"t");

	if(startTime < 0){
	
		startTime = fileSize + perturbationFileTime;
	
		if(startTime<0){ printf("Error: requested restart time of %d is invalid\n",startTime); exit(0);}
		
	} else if(startTime >= fileSize){ 
			printf("Error: requested restart time index of %d exceeds input file index size of %d\n",startTime,fileSize); exit(0);
		
	}
	
	return startTime;
}

/******************************************************************************
* Initialize the Coriolis paramter given an array on latitudes
*
*******************************************************************************/
void initialize_coriolis_parameter(double *lats){
	
	for(int j=0;j<NY;j++){ 
		
		f[j] = 2*fc*sin(lats[j]*trigpi/180.);
		dfdy[j] = 2*fc*cos(lats[j]*trigpi/180.) * (1./meters_per_degree) * (trigpi/180.) ;
	}
}

/******************************************************************************
* Initializes the stretched grid. Used global variables 
* 'height_lowest_level' and 'index_lowest_level' (could add these as input
* parameters)
*******************************************************************************/
void setup_vertical_height_levels(){
	
	for(int k=0;k<NZ;k++){

		zu[k] = ((double)k-0.5)*dz;
		zw[k] = ((double)k-1)*dz;
	}

	stretched_grid(&zsu[0],&mu[0],&zsw[0],&mw[0],height_lowest_level,index_lowest_level);
}

/*********************************************************************
* Construct the rest of the base state from the base state 
* temperature (tb) and moisture (qb) profiles
* 
* @param tb - base state temperature
* @param qb - base state moisture
**********************************************************************/
void initialize_vertical_basic_state2(double *tb,double *qb){
	
	//-----------------------------------------------------
	// create base state potential temperature profile
	//-----------------------------------------------------
	tbw[0] = 0.5*(tb[0]+tb[1]);

	for(int k=1;k<NZ;k++){ tbw[k] = 0.5*(tb[k]+tb[k-1]);}

	//-----------------------------------------------------
	// create base state virtual potential temperature profile
	//-----------------------------------------------------
	for(int k=0;k<NZ;k++){ tbv[k]=tb[k]*(1.0+0.61*qb[k]);}

	//-----------------------------------------------------
	// create base state pressure profile
	//-----------------------------------------------------
	double pisfc = pow((pressfc/p0),( Rd/cp));
	
	if(STRETCHED_GRID){
		
		pib[1] = pisfc - grav * (zsu[1]-0) / (cp * tbv[1]);
		pib[0] = pisfc + grav * (0-zsu[0]) / (cp * tbv[0]);

		for(int k=2;k<NZ;k++){ pib[k] = pib[k-1]-grav*(zsu[k]-zsu[k-1])/(cp*(  0.5*(tbv[k]+tbv[k-1]) ));}
		
	} else {
		
		pib[1] = pisfc - grav * 0.5 * dz / (cp * tbv[1]);
		pib[0] = pisfc + grav * 0.5 * dz / (cp * tbv[0]);

		for(int k=2;k<NZ;k++){ pib[k] = pib[k-1]-grav*dz/(cp*(  0.5*(tbv[k]+tbv[k-1]) ));}
	}

	//-----------------------------------------------------
	// create base state density profile
	//-----------------------------------------------------
	for(int k=0;k<NZ;k++){ rhou[k] = p0*pow(pib[k],(cv/Rd)) / (Rd*tbv[k]);}

	rhow[1] = pressfc / (Rd*tbw[1]);

	for(int k=1;k<NZ;k++){ rhow[k] = 0.5*(rhou[k]+rhou[k-1]);}

	for(int k=1;k<NZ;k++){ rhoavg = rhoavg + rhou[k];}

	rhoavg = rhoavg/((double)NZ-2);

	//-----------------------------------------------------
	// set values for non-physical points
	//-----------------------------------------------------
	tbv[0]=tbv[1];
	tbv[NZ-1]=tbv[NZ-2];
	tbw[0]=tbw[1];
	tbw[NZ-1]=tbw[NZ-2];
	pib[0]=pib[1];
	pib[NZ-1]=pib[NZ-2];
	rhow[0]=rhow[1];

	for(int k=0;k<NZ;k++){

		one_d_rhow[k] = 1. / rhow[k];
		one_d_rhou[k] = 1. / rhou[k];
	}
	
}

/***********************************************************************
* Initialize vertically varying, horizontally uniform base state
*
* ibase, jbase - the i and j indices within the horizontally
* varying basic state array to use as the base state
************************************************************************/
void initialize_vertical_basic_state(int ibase,int jbase){
	
	//-----------------------------------------------------
	// create base state potential temperature profile
	//-----------------------------------------------------
	tb[0] = ITHBAR(ibase,jbase,0);

	for(int k=1;k<NZ;k++){ tb[k] = ITHBAR(ibase,jbase,k);}

	tb[0]=tb[1];
	tb[NZ-1]=tb[NZ-2];

	//-----------------------------------------------------
	// create base state specific humidity profile
	//-----------------------------------------------------
	qb[0] = IQBAR(ibase,jbase,0);

	for(int k=1;k<NZ;k++){ qb[k] = IQBAR(ibase,jbase,k);}

	qb[0]=qb[1];
	qb[NZ-1]=qb[NZ-2];


	initialize_vertical_basic_state2(tb,qb);
}

/*********************************************************************
* 
*
**********************************************************************/
void print_vertical_basic_state(){

	if(STRETCHED_GRID){
		for(int k=0;k<NZ;k++){
			printf("%d\tzu = %.4f\tzw = %.4f\tpib = %.5f\ttb = %.1f\ttbv = %.1f\trhou = %.3f\trhow = %.3f\tqb = %.3f\n",
					k,zsu[k]/1000,zsw[k]/1000,pib[k],tb[k],tbv[k],rhou[k],rhow[k],qb[k]*1000);
		}
	} else {
		for(int k=0;k<NZ;k++){
			printf("%d\tzu = %.2f\tzw = %.2f\tpib = %.3f\ttb = %.1f\ttbv = %.1f\trhou = %.3f\trhow = %.3f\tqb = %.3f\n",
					k,zu[k]/1000,zw[k]/1000,pib[k],tb[k],tbv[k],rhou[k],rhow[k],qb[k]*1000);			
		}
	}
}


/*********************************************************************
*
*
**********************************************************************/
void zero_perturbation_fields(){
	
	memset(&U(0,0,0), 0,NX*NY*NZ*sizeof(double));
	memset(&UP(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&UM(0,0,0),0,NX*NY*NZ*sizeof(double));	
	memset(&V(0,0,0), 0,NX*NY*NZ*sizeof(double));
	memset(&VP(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&VM(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&W(0,0,0), 0,NX*NY*NZ*sizeof(double));
	memset(&WP(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&WM(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&TH(0,0,0), 0,NX*NY*NZ*sizeof(double));
	memset(&THP(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&THM(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&PI(0,0,0),0,NX*NY*NZ*sizeof(double));
}

/*********************************************************************
*
*
**********************************************************************/
void const_times_array3D(double *var,double a,int size){
	
	for(int i=0;i<size;i++){ var[i] *= a;}
}

/*********************************************************************
*
*
**********************************************************************/
void reinitialize_perturbation(int time,double scale){
		
	zero_perturbation_fields();

	//if(USE_MICROPHYSICS){ init_microphysics(); }

	initialize_from_output_serial(inputPerturbationFileName,time);
	
	int size = NX*NY*NZ;
	
	const_times_array3D(&U(0,0,0),scale,size);
	const_times_array3D(&V(0,0,0),scale,size);
	const_times_array3D(&W(0,0,0),scale,size);
	const_times_array3D(&TH(0,0,0),scale,size);
	
	const_times_array3D(&UM(0,0,0),scale,size);
	const_times_array3D(&VM(0,0,0),scale,size);
	const_times_array3D(&WM(0,0,0),scale,size);
	const_times_array3D(&THM(0,0,0),scale,size);	
}

/*********************************************************************
*
*
**********************************************************************/
void reinitialize(){

	mtime = 0;
	bigcounter = 0;

	//----------------------------------------------------------------
	// Initialize arrays for SOR
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		//sor[0][i][j] = 0;
		//sor[1][i][j] = 0;
		//sor[2][i][j] = 0;
		
	}}

	//----------------------------------------------------------------
	// Reset perturbation arrays
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		U(i,j,k) = 0;
		V(i,j,k) = 0;
		TH(i,j,k) = 0;
		W(i,j,k) = 0;
		PI(i,j,k) = 0;

		UP(i,j,k) = U(i,j,k);
		UM(i,j,k) = U(i,j,k);

		VP(i,j,k) = V(i,j,k);
		VM(i,j,k) = V(i,j,k);
	#if !HYDROSTATIC
		WP(i,j,k) = W(i,j,k);
		WM(i,j,k) = W(i,j,k);
	#endif
		THP(i,j,k) = TH(i,j,k);
		THM(i,j,k) = TH(i,j,k);
	}}}
}

/********************************************************
* 
*********************************************************/
void init_basic_state_meridional_velocity(int middle_j){
	
	for(int i=0;i<NX-1;i++){
	for(int k=0;k<NZ;k++){
		
		for(int j=middle_j;j<NY-1;j++){
		
			IVBAR(i,j+1,k) = IUBAR(i,j,k) - IUBAR(i+1,j,k) + IVBAR(i,j,k);
		
		}
		
		for(int j=middle_j-1;j>=0;j--){
		
			IVBAR(i,j,k) = IUBAR(i+1,j,k) - IUBAR(i,j,k) + IVBAR(i,j+1,k);
		
		}
	}}
}

/********************************************************
* Initializes the basic state vertical velocity by solving
* the anelastic continuity equation. The basic state winds
* should be defined before executing this subroutine.
*********************************************************/
void init_basic_state_vertical_velocity(){
	
	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){

		IWBAR(i,j,HTOPOFULL(i,j)) = 0;		// boundary condition = zero vertical 
		IWBAR(i,j,HTOPOFULL(i,j)+1) = 0;	// velocity at bottom

		for(int k=HTOPOFULL(i,j)+2;k<NZ-1;k++){
			
			IWBAR(i,j,k) = (rhow[k-1]/rhow[k])*IWBAR(i,j,k-1) 
				+ (rhou[k-1]/rhow[k])*
				(
					-(IUBAR(i+1,j,k-1) - IUBAR(i,j,k-1))/dx  	
		      		-(IVBAR(i,j+1,k-1) - IVBAR(i,j,k-1))/dy
				)*DZU(k);
			
		}

		IWBAR(i,j,NZ-1) = IWBAR(i,j,NZ-2);
	}}
}

/********************************************************
* Initializes linear friction. at the lowest model level.
* If terrain option is used, sets linear friction for
* levels near the terrain.
*********************************************************/
void init_friction(){

	for(int i=0;i<NX;i++)
		for(int j=0;j<NY;j++)
			for(int k=0;k<NZ;k++)
				IFRICTION(i,j,k) = 0;

	if(USE_LINEAR_FRICTION && !USE_TERRAIN){

		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){

			if(ISLINEAR){ IFRICTION(i,j,HTOPOFULL(i,j)+1) = 2.0e-5;}
			else { IFRICTION(i,j,HTOPOFULL(i,j)+1) = 2.0e-5;}
		}}
	}
	
	if(USE_LINEAR_FRICTION && !MERIDIONAL_CROSS_SECTION && USE_TERRAIN){

		for(int i=1;i<NX-1;i++){
		for(int j=1;j<NY-1;j++){
		for(int k=HTOPOFULL(i,j)+1;k<NZ;k++){
			
			if(
				(IISTOPO(i,j+1,k)==0 ||
				IISTOPO(i,j-1,k)==0 ||
				IISTOPO(i+1,j,k)==0 ||
				IISTOPO(i-1,j,k)==0) &&
				IFRICTION(i,j,k) == 0
				
			){
				//printf("%d %d %d %f %f\n",i,j,k,outLats[j],outLons[i]);
				IFRICTION(i,j,k) = 2.0e-5;
			}
		}}}
	}
}

/********************************************************
* Initialize topography
*********************************************************/
void init_topography(){

	int height = 0;

	memset(&HTOPOFULL(0,0), 0,NX*NY*sizeof(int));

	//for(int i=NX-1;i>=0;i--){
	for(int i=NX-2;i>=0;i--){
	for(int j=NY-2;j>=0;j--){

		height = 0;

		for(int k=0;k<NZ;k++){

		#if !STRETCHED_GRID
			#if !MERIDIONAL_CROSS_SECTION
				if(zu[k] < ITOPO(i,j,k) && ITOPO(i,j,k) > 0)
			#else
				if(zu[k] < ITOPO(linear_lon,j,k) && ITOPO(linear_lon,j,k) > 70000)
			#endif
		#else
			#if !MERIDIONAL_CROSS_SECTION
				if(zsu[k] < ITOPO(i,j,k) && ITOPO(i,j,k) > 0)
			#else
				//if(i < linear_lon && (zsu[k] < ITOPO(linear_lon,j,k) && ITOPO(linear_lon,j,k) > 0) )	
				if( zsu[k] < ITOPO(linear_lon,j,k) && ITOPO(linear_lon,j,k) > 0 )
			#endif
		#endif	
			{
				
				IISTOPO(i,j,k) = 0;
				IUISTOPO(i,j,k) = 0;
				IVISTOPO(i,j,k) = 0;
				IUISTOPO(i+1,j,k) = 0;
				IVISTOPO(i,j+1,k) = 0;
				height = k;
#if 0
			} else if(i >= linear_lon && (zsu[k] < ITOPO(i,j,k) && ITOPO(i,j,k) > 0) ){
				
				IISTOPO(i,j,k) = 0;
				IUISTOPO(i,j,k) = 0;
				IVISTOPO(i,j,k) = 0;
				IUISTOPO(i+1,j,k) = 0;
				IVISTOPO(i,j+1,k) = 0;
				height = k;
#endif
			} else {

				IISTOPO(i,j,k) = 1;
				IUISTOPO(i,j,k) = 1;
				IVISTOPO(i,j,k) = 1;
			}
		}

		HTOPOFULL(i,j) = height;
		//printf("%d %d %d\n",i,j,height);
	}}

}

/********************************************************
* Make topography height equal zero
*********************************************************/
void init_topography_to_zero(){
	
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
		
		HTOPOFULL(i,j) = 0;
		
	for(int k=0;k<NZ;k++){
	
		ITOPO(i,j,k) = 0;
		IISTOPO(i,j,k) = 1;
		IUISTOPO(i,j,k) = 1;
		IVISTOPO(i,j,k) = 1;
	}}}
	
}

/*********************************************************************
* Alter basic state humidity
**********************************************************************/
void change_humidity(){

	double smr,temp,pres,rh;
	
	for(int k=0;k<NZ;k++){
		
		temp = tb[k]*pib[k];
		pres = p0*pow(pib[k],(cp/Rd));
		
		smr = get_QV_Sat(temp,pres);
		
		//printf("%d %f %f %f\n",k,smr*1000,qb[k]*1000,qb[k]/smr);
		
	}

	//initialize_landsea(landseaMaskFile);

	const double height_limit_low = 2000.;
	const double height_limit_high = 10000.;
	const double relative_humidity = 0.90;

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=1;k<NZ;k++){
		//-----------------------------------------------------------------
		// Limit region
		//-----------------------------------------------------------------
		if(outLats[j] > 18 && outLats[j] < 23 && outLons[i] > 82 && outLons[i] < 87){
			
			temp = (ITHBAR(i,j,k)+tb[k]) * IPBAR(i,j,k);	// full, actual temperature
			pres = p0*pow(IPBAR(i,j,k),(cp/Rd));			// full, dimensional pressure

			smr = get_QV_Sat(temp,pres);					// calculate saturation mixing ratio

			rh = (IQBAR(i,j,k)+qb[k])/smr;					// calculate relative humidity
		
			// change humidity
			if(zsu[k] >= height_limit_low && zsu[k] <=  height_limit_high){ 
			
				//IQBAR(i,j,k) = relative_humidity * smr - qb[k];
		
			}

		}
	}}}

	//free(landsea);
}

/*********************************************************************
* Calculate heights for stretched grid.
*
**********************************************************************/
void stretched_grid(double * zu,double *mu,double * zw,double *mw,double z0,int lev){
	
	double c1;
	double c2;
	
	int nt = NZ-1;
	double dzp = dz;
	
	c2 = ( 1 - z0/(dzp*lev)) / ( (nt-1) *dzp - dzp*lev); 
	
	c1 = z0 / (dzp*lev) - c2*dzp*lev;

	//printf("%f %f\n",c1,c2);
	
	double zh1,zh2;
	double m1,m2;
	double fullLev;
	double halfLev;
	
	for(int i=0;i<=nt;i++){
		
		fullLev = ((double)i-1.0)*dzp;
		halfLev = ((double)i-0.5)*dzp;
		
		zh1 = (c1+c2*fullLev)*fullLev;
		m1 = 1 / (c1+2*c2*fullLev);
	
		zh2 = (c1+c2*halfLev)*halfLev;
		m2 = 1 / (c1+2*c2*halfLev);
		
		//printf("%d %.0f %f %f %f\n",i,zh1,m1,zh2,m2);
		
		zu[i] = zh2;
		zw[i] = zh1;
		
		mu[i] = m2;
		mw[i] = m1;
	}
	
}

/*********************************************************************
* Calculate saturation mixing ratio.
*
* temperature - full, actual temperature field (K, not potential, not perturation)
* pressure - full, dimensional pressure (Pa)
**********************************************************************/
double get_QV_Sat(double temperature,double pressure){

	double esl = 611.2 * exp(17.67 * (temperature-273.15) / (temperature - 29.65) );

	return 0.62197 * esl / (pressure-esl);
}

/*********************************************************************
* Load data from model output file for use as an initial condition for 
* the perturbation fields. Assumes that the current model grid settings 
* (i.e. NX,dx, etc.) and the those of the input file are identical
*
* filename 	- file containing input fields for perturbation
* varname 	- name of variable within file
* var,mvar	- fully interpolated arrays as output (if mvar is a null pointer
*			  nothing will be output to it)
* time 		- time at which to initialize, negative for basic state
**********************************************************************/
void load_from_output(const char *filename,const char *varname,double *var,double *mvar,int time){
	
	//-----------------------------------------------------------------------
	// Does the variable exist in the file?
	//-----------------------------------------------------------------------
	if(fileHasVar(filename,varname)){
	
		//----------------------------------------------------------------
		// Load data on root process
		//----------------------------------------------------------------	
		if(rank==0){ 
	
			if(time >= 0){ get_data_at_time(filename,varname,time,    iubar);}
			else { 		   get_data(        filename,varname,NX*NY*NZ,iubar);}
		}
	
		distributeArray(var,iubar);	// distribute to other processes
	
		exchange(var);			// processes exchange boundaries
	
		//----------------------------------------------------------------
		// Load interpolated fields into the second array
		//----------------------------------------------------------------	
		if(mvar!=NULL){
	
			for(int i=0;i<fNX*fNY*fNZ;i++){ mvar[i] = var[i];}
		}
	//----------------------------------------------------------------------
	// If variable does not exist, do nothing but print out a warning.
	//----------------------------------------------------------------------			
	} else {
		if(!PARALLEL || rank==0){
			printf("Variable %s not in file! Initializing to zero\n",varname);
		}
	}
}

/*********************************************************************
* Compare the values of two arrays. Return true if they are the same
* within a given tolerance level, false otherwise.
*
**********************************************************************/
bool compare_values(double *zin,double *zout,int nz,double diff_tolerance){
	
	bool same = true;

	for(int i=0;i<nz;i++){
		//printf("%d %f %f\n",i,zin[i],zout[i]);
		if( fabs(zin[i]-zout[i]) > diff_tolerance ){
			
			same = false;
			break;
		}
	}

	return same;
}

/*********************************************************************
* 
**********************************************************************/
bool hasSameHeights(const char *myfilename,int zdim){
	
	const double tolerance = 1;	// in meters
	
	double *zlevs = get_data2(myfilename,"zu",zdim);

	bool isSame = compare_values(&ZU(0),zlevs,zdim,tolerance);

	free(zlevs);

	return isSame;
}

/*********************************************************************
* Load data for perturbation fields from file for use as an initial condition 
* and interpolate it to the model grid. If file does not contain variable,
* do nothing but print out warning message.
*
* filename 	- file containing input fields for perturbation
* varname 	- name of variable within file
* time 		- time at which to initialize
* varIn		- allocated memory for input field at original resolution (xdim x ydim x zdim array)
* var_interpz - allocated memory for input field interpolated to given heights height (xdim x ydim x NZ array)
* varOut   	- fully interpolated output arrays
* levs_out	- output heights
* myLons	- input longitude
* myLats	- input latitude
* zlevs		- input heights
* xdim		- input x-dimension
* ydim		- input y-dimension
* zdim		- input z-dimension
* myDX		- input dx grid spacing
* myDY		- input dy grid spacing
**********************************************************************/
void interpolate_from_output(
			const char *filename,const char *varname,
			int time,
			double *varIn,double *var_interpz,double *varOut,
			double *levs_out,
			double *myLons,double *myLats,double *zlevs,
			int xdim,int ydim,int zdim,
			double myDX,double myDY
	){
		//-----------------------------------------------------------------------
		// Does the variable exist in the file?
		//-----------------------------------------------------------------------
		if(fileHasVar(filename,varname)){
		
			for(int i=0;i<xdim*ydim*zdim;i++){ varIn[i] = 0;}
			for(int i=0;i<xdim*ydim*NZ;i++){ var_interpz[i] = 0;}
	
			//----------------------------------------------------------------
			// If time is negative, assume data has no time dimension
			//----------------------------------------------------------------
			if(time >= 0){ get_data_at_time(filename,varname,(size_t)time,varIn,xdim,ydim,zdim);}
			else { 		   get_data(        filename,varname,xdim*ydim*zdim,varIn);}
		
			vert_interpolate_1d_from_model(levs_out,zlevs, varIn, xdim,ydim,zdim,NZ,var_interpz);
			
			horz_interpolate_from_model(var_interpz,varOut,xdim,ydim,NZ,false,false,
			myDX/meters_per_degree,myDY/meters_per_degree,outLons[0]-myLons[0],outLats[0]-myLats[0]);
			
			//if(PERIODIC_BOUNDARIES){ periodic_EW_boundaries(varOut);}	// might need for serial version? probably not for parallel
			
		//----------------------------------------------------------------------
		// If variable does not exist, do nothing but print out a warning.
		//----------------------------------------------------------------------	
		} else {
			
			if(!PARALLEL || rank==0){
				printf("Variable %s not in file! Initializing to zero\n",varname);
			}
		}

}

/*********************************************************************
* Load data for perturbation fields from file for use as an initial condition 
* and interpolate it to the model grid.
*
* filename 	- file containing input fields for perturbation
* varname 	- name of variable within file
* time 		- time at which to initialize
* varIn		- allocated memory for input field at original resolution (xdim x ydim x zdim array)
* var_interpz - allocated memory for input field interpolated to given heights height (xdim x ydim x NZ array)
* var,mvar	- fully interpolated output arrays
* levs_out	- output heights
* myLons		- input longitude
* myLats		- input latitude
* zlevs		- input heights
* xdim		- input x-dimension
* ydim		- input y-dimension
* zdim		- input z-dimension
* myDX		- input dx grid spacing
* myDY		- input dy grid spacing
**********************************************************************/
void load_interpolate_from_output(
			const char *filename,const char *varname,
			size_t time,
			double *varIn,double *var_interpz,
			double *var,double *mvar,
			double *levs_out,
			double *myLons,double *myLats,double *zlevs,
			int xdim,int ydim,int zdim,
			double myDX,double myDY
	){
	//-----------------------------------------------------------------
	// Read file and interpolate on root process
	//-----------------------------------------------------------------
	if(rank==0){
		memset(iubar,0,NX*NY*NZ*sizeof(double));
		interpolate_from_output(filename,varname,time,varIn,var_interpz,iubar,levs_out,myLons,myLats,zlevs,xdim,ydim,zdim,myDX,myDY);
	}
	
	//-----------------------------------------------------------------
	// Distribute results to other processes
	//-----------------------------------------------------------------
	distributeArray(var,iubar);
	
	exchange(var);
	
	for(int i=0;i<fNX*fNY*fNZ;i++){ 
		//var[i] *= 0.2;
		mvar[i] = var[i];
	}
}

/*********************************************************************
* Get the dimensions and grid spacing from a file a tell whether or not
* they match the model's grid spacing and dimensions
*
* @param myfilename
* @param xdim		- x dimension of input file
* @param ydim		- y dimension of input file
* @param zdim		- z dimension of input file
* @param interp_dx	- x direction grid spacing of input file
* @param interp_dy	- y direction grid spacing of input file
* @param interp_dz	- z direction grid spacing of input file
* @return 1 if all parameters match the model's grid spacing and dimensions,
*		  0 otherwise
**********************************************************************/
int get_file_dims(const char *myfilename,size_t *xdim,size_t *ydim,size_t *zdim,double *interp_dx, double *interp_dy,double *interp_dz){
	
	size_t dims[3];
	float grid_spacing[3];
	int same;
	
	get_dims(myfilename,"x","y","z",dims);
	get_grid_spacing(myfilename,"dx","dy","dz",grid_spacing);

	*interp_dx = grid_spacing[0];
	*interp_dy = grid_spacing[1];
	*interp_dz = grid_spacing[2];
	
	*xdim = dims[0];
	*ydim = dims[1];
	*zdim = dims[2];

	if(*xdim == NX && *ydim == NY && *zdim == NZ && *interp_dx == dx && *interp_dy == dy && *interp_dz == dz){
		same = 1;
	} else { 
		same = 0;
	}
	
	return same;
}


/*********************************************************************
* Initialize model perturbation fields from input file. Will determine
* whether or not interpolation is required.
*
* myfilename - file containing input fields for perturbations
* time		 - time at which to initialize
**********************************************************************/
void initialize_from_output_parallel(const char *myfilename,size_t time){
	
	size_t xdim,ydim,zdim;
	double interp_dx,interp_dy,interp_dz;
	int isSame = 0;
	
	
	//------------------------------------------------------------------------------
	// Test files for same dimensions and grid spacing
	//------------------------------------------------------------------------------
	if(rank==0){
		
		isSame = get_file_dims(myfilename,&xdim,&ydim,&zdim,&interp_dx,&interp_dy,&interp_dz);
		
		// Even if the dimensions are the same, the height levels could be different
		// and thus still require interpolation. Check for this condition.
		if(isSame==1){ 		if(!hasSameHeights(myfilename,zdim)){ isSame = 0;}		}
	}

	broadcast_data_int(1,&isSame);	// tell other processors whether the dimensions/grid spacing are the same

	//------------------------------------------------------------------------------
	// If the files have the same dimensions and grid spacing, just load the data
	// directly onto the model grid.
	//------------------------------------------------------------------------------
	if(isSame==1){
		
		if(rank == 0){ printf("Files have same dimensions and grid spacing. Loading initial conditions from file %s\n",myfilename);}
		
		load_from_output(myfilename,"u-wind",us,ums,time);
		load_from_output(myfilename,"v-wind",vs,vms,time);
		load_from_output(myfilename,"w-wind",ws,wms,time);
		load_from_output(myfilename,"theta",ths,thms,time);
	
		if(OUTPUT_FRICTION_TEND){
			load_from_output(myfilename,"fric",frictions,frictions,time);
		}
	
		if(USE_MICROPHYSICS){
		
			load_from_output(myfilename,"qv",qvs,qvms,time);
			load_from_output(myfilename,"qc",qcs,qcms,time);
			load_from_output(myfilename,"qr",qrs,qrms,time);
			
			if(USE_ICE){
				
				load_from_output(myfilename,"qs",qss,qsms,time);
				load_from_output(myfilename,"qi",qis,qims,time);				
			}
		}
	//------------------------------------------------------------------------------
	// If the files have different dimensions or grid spacing, need to interpolate
	// input fields to model grid.
	//------------------------------------------------------------------------------
	} else {
		
		double *zlevs=NULL,*var=NULL,*var_interpz=NULL,*myLons=NULL,*myLats=NULL;
		//------------------------------------------------------------------------------
		// Allocate temporary memory on the root process for interpolated fields. 
		// DO NOT TRY TO ACCESS THESE VARIABLES ON THE OTHER PROCESSES!
		//------------------------------------------------------------------------------
		if(rank==0){
			printf("---------------------------------------------------\n");
			printf("Files have different dimensions or grid spacing.\n"
				"Loading initial conditions from file\n %s \n"
				"with dimension %lu x %lu x %lu and\n grid spacing dx = %f at time = %lu\n",
				myfilename,xdim,ydim,zdim,interp_dx,time);
			printf("---------------------------------------------------\n");
			
			var = (double*)calloc(xdim*ydim*zdim,sizeof(double));
			var_interpz = (double*)calloc(xdim*ydim*NZ,sizeof(double));
	
			zlevs = get_data2(myfilename,"zu",zdim);
			myLons = get_data2(myfilename,"lon",xdim);
			myLats = get_data2(myfilename,"lat",ydim);
		}
		//------------------------------------------------------------------------------
		// Get data, interpolate it, distribute it to each process
		//------------------------------------------------------------------------------
		load_interpolate_from_output(myfilename,"u-wind",time,var,var_interpz,us,ums,&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		load_interpolate_from_output(myfilename,"v-wind",time,var,var_interpz,vs,vms,&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);	
		load_interpolate_from_output(myfilename,"w-wind",time,var,var_interpz,ws,wms,&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		load_interpolate_from_output(myfilename,"theta",time,var,var_interpz,ths,thms,&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
				
		if(USE_MICROPHYSICS){
			
			load_interpolate_from_output(myfilename,"qv",time,var,var_interpz,qvs,qvms,&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
			load_interpolate_from_output(myfilename,"qc",time,var,var_interpz,qcs,qcms,&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);	
			load_interpolate_from_output(myfilename,"qr",time,var,var_interpz,qrs,qrms,&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
			
			if(USE_ICE){
				
				load_interpolate_from_output(myfilename,"qs",time,var,var_interpz,qss,qsms,&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);	
				load_interpolate_from_output(myfilename,"qi",time,var,var_interpz,qis,qims,&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
			}
		}
		

		
		if(rank==0){ free(zlevs); free(var); free(var_interpz); free(myLons); free(myLats);}
	}
	
}

/*********************************************************************
* Initialize model perturbation fields from input file. Will determine
* whether or not interpolation is required.
*
* myfilename - file containing input fields for perturbations
* time		 - time at which to initialize
**********************************************************************/
void initialize_from_output_serial(const char *myfilename,size_t time){
	
	size_t xdim,ydim,zdim;
	double interp_dx,interp_dy,interp_dz;
	int isSame = 0;
	
	//------------------------------------------------------------------------------
	// Test files for same dimensions and grid spacing
	//------------------------------------------------------------------------------	
	isSame = get_file_dims(myfilename,&xdim,&ydim,&zdim,&interp_dx,&interp_dy,&interp_dz);

	// Even if the dimensions are the same, the height levels could be different
	// and thus still require interpolation. Check for this condition.
	if(isSame==1){ 		if(!hasSameHeights(myfilename,zdim)){ isSame = 0;}		}

	//------------------------------------------------------------------------------
	// If all dimensions and grid spacing is the same, simply load input
	//------------------------------------------------------------------------------	
	if(isSame==1){
		
		printf("---------------------------------------------\n");
		printf("Files have same dimensions and grid spacing.\n"
			   "Loading initial conditions from file\n %s\n",myfilename);
		printf("---------------------------------------------\n");
	
		get_data_at_time(myfilename,"u-wind",time,&U(0,0,0));
		get_data_at_time(myfilename,"v-wind",time,&V(0,0,0));
		get_data_at_time(myfilename,"w-wind",time,&W(0,0,0));
		get_data_at_time(myfilename,"theta",time,&TH(0,0,0));
		
		if(USE_MICROPHYSICS){
		
			get_data_at_time(myfilename,"qv",time,&QV(0,0,0));
			get_data_at_time(myfilename,"qc",time,&QC(0,0,0));
			get_data_at_time(myfilename,"qr",time,&QR(0,0,0));
			
			if(USE_ICE){
				get_data_at_time(myfilename,"qi",time,&QI(0,0,0));
				get_data_at_time(myfilename,"qs",time,&QS(0,0,0));
			}
		}
		#if 0
		for(int i=0;i<NX*NY*NZ;i++){
			us[i] *= 0.2;
			vs[i] *= 0.2;
			ws[i] *= 0.2;
			ths[i] *= 0.2;
			qvs[i] *= 0.2;
		}
		#endif
	//------------------------------------------------------------------------------
	// If the files have different dimensions or grid spacing, need to interpolate
	// input fields to model grid.
	//------------------------------------------------------------------------------
	} else {
		//------------------------------------------------------------------------------
		// Allocate temporary memory on the root process for interpolated fields.
		//------------------------------------------------------------------------------			
		printf("Files have different dimensions or grid spacing.\n"
			"Loading initial conditions from file %s \n"
			"with dimension %lu x %lu x %lu and grid spacing dx = %f at time = %lu\n",
			myfilename,xdim,ydim,zdim,interp_dx,time);
		
		double *var = (double*)calloc(xdim*ydim*zdim,sizeof(double));
		double *var_interpz = (double*)calloc(xdim*ydim*NZ,sizeof(double));

		double *zlevs = get_data2(myfilename,"zu",zdim);
		double *myLons = get_data2(myfilename,"lon",xdim);
		double *myLats = get_data2(myfilename,"lat",ydim);
		
		//------------------------------------------------------------------------------
		// Get data, interpolate it
		//------------------------------------------------------------------------------		
		interpolate_from_output(myfilename,"u-wind",time,var,var_interpz,&U(0,0,0), &ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		interpolate_from_output(myfilename,"v-wind",time,var,var_interpz,&V(0,0,0), &ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);		
		interpolate_from_output(myfilename,"w-wind",time,var,var_interpz,&W(0,0,0), &ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		interpolate_from_output(myfilename,"theta", time,var,var_interpz,&TH(0,0,0),&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);

		if(USE_MICROPHYSICS){
			interpolate_from_output(myfilename,"qv",time,var,var_interpz,&QV(0,0,0),&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);		
			interpolate_from_output(myfilename,"qc",time,var,var_interpz,&QC(0,0,0),&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
			interpolate_from_output(myfilename,"qr",time,var,var_interpz,&QR(0,0,0),&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
			
			if(USE_ICE){
				interpolate_from_output(myfilename,"qi",time,var,var_interpz,&QI(0,0,0),&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
				interpolate_from_output(myfilename,"qs",time,var,var_interpz,&QS(0,0,0),&ZU(0),myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);			
			}
		}
		
		free(zlevs); free(var); free(var_interpz); free(myLons); free(myLats);
	}
	
	//------------------------------------------------------------------------------
	// Initialize ...
	//------------------------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
		
		UM(i,j,k) = U(i,j,k);
		VM(i,j,k) = V(i,j,k);
		WM(i,j,k) = W(i,j,k);
		THM(i,j,k) = TH(i,j,k);
	}}}
	
	if(USE_MICROPHYSICS){
		
		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
		for(int k=0;k<NZ;k++){
			
			QVM(i,j,k) = QV(i,j,k);
			QCM(i,j,k) = QC(i,j,k);
			QRM(i,j,k) = QR(i,j,k);
			
		}}}
		
		if(USE_ICE){
			
			for(int i=0;i<NX;i++){
			for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
			
				QIM(i,j,k) = QI(i,j,k);
				QSM(i,j,k) = QS(i,j,k);
			
			}}}
		}
	}
}

/********************************************************
* Test the model's coordinates to make sure they are within
* physical bounds. Terminate program if not.
*********************************************************/
void test_coordinates(){
	if(outLats[ 0  ] <= -90 ){ printf("Error!: Latitude value of %f degrees is out of range\n\n",outLats[0]); exit(0);}
	if(outLats[NY-1] >= 90 ){ printf("Error!: Latitude value of %f degrees is out of range\n\n",outLats[NY-1]); exit(0);}
}

/********************************************************
*
* @param leftLon - longitude at leftmost edge of input data
* @param lowerLat - latitude at bottommost edge of input data
*********************************************************/
void setup_coordinates(double leftLon,double lowerLat){

	double lon_shift;
	
	if(SHIFT_PRIME_MERIDIAN){ lon_shift = 180;} else { lon_shift = 0;}

	outLons[0] = leftLon  + lonoffset + dx/meters_per_degree - lon_shift;
	outLats[0] = lowerLat + latoffset + dy/meters_per_degree;

	for(int i=1;i<NX;i++){ outLons[i] = outLons[i-1] + dx/meters_per_degree;}//printf("%d %f\n",i,outLons[i]);}
	for(int i=1;i<NY;i++){ outLats[i] = outLats[i-1] + dy/meters_per_degree;}//printf("%d %f\n",i,outLats[i]);}
	
	test_coordinates();
}

/*********************************************************************
* From a model output file, initialized the basic state. Will perform
* interpolation if grid dimensions do not match.
*
* myfilename - file from which to get fields
**********************************************************************/
void initialize_basic_state_from_output_file(const char *myfilename){
		
	setup_memory_allocation();

	//----------------------------------------------------------------
	// Set up coordinate system
	//----------------------------------------------------------------
	setup_vertical_height_levels();
	
	setup_coordinates(0,-89.462997);
	//setup_coordinates(0,-90);

	//----------------------------------------------------------------
	// Coriolis parameter
	//----------------------------------------------------------------
	initialize_coriolis_parameter(&outLats[0]);

	//----------------------------------------------------------------
	// Initialize topographic array
	//----------------------------------------------------------------
	init_topography_to_zero();
	//----------------------------------------------------------------
	// Initialize friction array
	//----------------------------------------------------------------
	init_friction();
	//------------------------------------------------------------------------------
	// Get grid parameters and test files for same dimensions and grid spacing
	//------------------------------------------------------------------------------		
	size_t xdim,ydim,zdim;
	double interp_dx,interp_dy,interp_dz;
	int isSame = 0;
	
	isSame = get_file_dims(myfilename,&xdim,&ydim,&zdim,&interp_dx,&interp_dy,&interp_dz);
	
	// Even if the dimensions are the same, the height levels could be different
	// and thus still require interpolation. Check for this condition.
	if(isSame==1){ 		if(!hasSameHeights(myfilename,zdim)){ isSame = 0;}		}
	
	//----------------------------------------------------------------------------------
	// If the grid parameters of the input file match the current model
	// setting, no need to interpolate, just load values
	//----------------------------------------------------------------------------------
	if(isSame==1){
		
		if(VERBOSE){
			printf("---------------------------------------------\n");
			printf("Basic state has same dimensions and grid spacing as input file.\n"
				    "Loading basic state from file\n %s\n",myfilename);
			printf("---------------------------------------------\n");
		}
	
		get_data(myfilename,"ubar",NX*NY*NZ,&IUBAR(0,0,0));
		get_data(myfilename,"vbar",NX*NY*NZ,&IVBAR(0,0,0));
		get_data(myfilename,"wbar",NX*NY*NZ,&IWBAR(0,0,0));
		get_data(myfilename,"thbar",NX*NY*NZ,&ITHBAR(0,0,0));
		get_data(myfilename,"qbar",NX*NY*NZ,&IQBAR(0,0,0));
		get_data(myfilename,"pbar",NX*NY*NZ,&IPBAR(0,0,0));
		
		get_data(myfilename,"tb",NZ,&tb[0]);		
		get_data(myfilename,"qb",NZ,&qb[0]);
	
		//--------------------------------------------------
		// Topography
		//--------------------------------------------------
		if(USE_TERRAIN){
			
			double *topo2 = get_data2(myfilename,"topo",NX*NY);
	
			for(int i=0;i<NX;i++){
			for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				ITOPO(i,j,k) = topo2[j+i*NY];
			}}}
		
			free(topo2);
		}
		
	} else {			
	//------------------------------------------------------------------------------
	// Grid parameters do not match, need to interpolate
	//------------------------------------------------------------------------------
		printf("---------------------------------------------\n");
		printf("Files for basic state have different dimensions or grid spacing.\n"
			   "Loading initial conditions from file %s \n"
			   "with dimension %lu x %lu x %lu and grid spacing dx = %f\n",
				myfilename,xdim,ydim,zdim,interp_dx);
		printf("---------------------------------------------\n");
				
		//------------------------------------------------------------------------------
		// Temporary storage for interpolation
		//------------------------------------------------------------------------------	
		double *var = (double*)calloc(xdim*ydim*zdim,sizeof(double));
		double *var_interpz = (double*)calloc(xdim*ydim*NZ,sizeof(double));

		//------------------------------------------------------------------------------
		// Get coordinate information from file from which to interpolate
		//------------------------------------------------------------------------------
		double *zlevs = get_data2(myfilename,"zu",zdim);
		double *myLons = get_data2(myfilename,"lon",xdim);
		double *myLats = get_data2(myfilename,"lat",ydim);
		double *mytb = get_data2(myfilename,"tb",zdim);
		double *myqb = get_data2(myfilename,"qb",zdim);			
		double *myzu = get_data2(myfilename,"zu",zdim);

		//------------------------------------------------------------------------------
		// Realign grid
		//------------------------------------------------------------------------------		
		for(int i=0;i<NX;i++){
			//outLons[i] += (interp_dx-dx)/meters_per_degree + 1.0e-6;
			//printf("%d %f %f\n",i,myLons[i/2],outLons[i]);
		}
		
		for(int i=0;i<NY;i++){
			//outLats[i] += (interp_dy-dy)/meters_per_degree + 1.0e-6;
			//printf("%d %f %f\n",i,myLats[i/2],outLats[i]);
		}
		
		//------------------------------------------------------------------------------
		// Whether or not to interpolate to stretched grid
		//------------------------------------------------------------------------------
		double *zgrid;
		
		if(STRETCHED_GRID){ zgrid = &zsu[0];}
		else { 				zgrid = &zu[0]; }
			
		//------------------------------------------------------------------------------
		// Get data, interpolate it
		//------------------------------------------------------------------------------		
		interpolate_from_output(myfilename,"ubar", -1,var,var_interpz,&IUBAR(0,0,0), zgrid,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		interpolate_from_output(myfilename,"vbar", -1,var,var_interpz,&IVBAR(0,0,0), zgrid,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);		
		interpolate_from_output(myfilename,"wbar", -1,var,var_interpz,&IWBAR(0,0,0), zgrid,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		interpolate_from_output(myfilename,"thbar",-1,var,var_interpz,&ITHBAR(0,0,0),zgrid,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		interpolate_from_output(myfilename,"qbar", -1,var,var_interpz,&IQBAR(0,0,0), zgrid,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		interpolate_from_output(myfilename,"pbar", -1,var,var_interpz,&IPBAR(0,0,0), zgrid,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		
		vert_interpolate_1d(zdim,NZ,myzu,zgrid,mytb,tb);
		vert_interpolate_1d(zdim,NZ,myzu,zgrid,myqb,qb);
		
		tb[0]=tb[1];
		tb[NZ-1]=tb[NZ-2];

		qb[0]=qb[1];
		qb[NZ-1]=qb[NZ-2];

		for(int i=0;i<NZ;i++){

			//printf("%d %f %f %f %f %f %f\n",i,myzu[i],ZU(i),mytb[i],tb[i],myqb[i],qb[i]);
		}

		//------------------------------------------------------------------------------
		// Topography
		//------------------------------------------------------------------------------
		if(USE_TERRAIN){
			double *topo2 = get_data2(myfilename,"topo",xdim*ydim);
		
			for(int i=0;i<xdim;i++){
			for(int j=0;j<ydim;j++){
			for(int k=0;k<NZ;k++){
				var_interpz[k+j*NZ+i*NZ*ydim] = topo2[j+i*ydim];
			}}}
		
			horz_interpolate_from_model(var_interpz,&ITOPO(0,0,0),xdim,ydim,NZ,false,false,
			interp_dx/meters_per_degree,interp_dy/meters_per_degree,outLons[0]-myLons[0],outLats[0]-myLats[0]);
		
			free(topo2);
		}
		
		//------------------------------------------------------------------------------
		// Deallocate memory local to this subroutine
		//------------------------------------------------------------------------------
		free(zlevs); free(var); free(var_interpz); free(myLons); free(myLats); free(mytb); free(myqb); free(myzu);
		
	}
			
	initialize_vertical_basic_state2(tb,qb);
	
	if(VERBOSE){ print_vertical_basic_state();}
	
	if(USE_TERRAIN){ init_topography();}
			
}


/********************************************************
* Reverse the y- and z-coordinate 
*
*********************************************************/
void interpolate_terrain(const char * infile){

	/********************************************************
	* Get dimensions of data in files to determine how 
	* much memory to allocate
	*********************************************************/
	size_t dims[3];

	get_dims(infile,"lon","lat","level",dims);

	size_t xdim = dims[0]; 
	size_t ydim = dims[1];
	size_t zdim = dims[2];

	int size = xdim*ydim*zdim;

	//printf("Dimensions of input data are %zu %zu %zu\n",xdim,ydim,zdim);
	
	/********************************************************
	* Get data from file
	*********************************************************/
	double * topo2 = get_data2(infile,"zsfc", xdim*ydim);


	double * topozinterp = (double *)malloc(xdim*ydim*NZ*sizeof(double));

	for(size_t i=0;i<xdim;i++){
	for(size_t j=0;j<ydim;j++){
	for(size_t k=0;k < NZ;k++){
		topozinterp[d(i,j,k)] = topo2[d2(i,ydim-1-j)];
	}}}

	if(SHIFT_PRIME_MERIDIAN){

		flip_array(topozinterp,xdim,ydim,NZ);
	}

	/********************************************************
	* Get and process coordinates from input data
	*********************************************************/
	double * lat2  = get_data2(infile,"latitude", ydim);
	double * lon   = get_data2(infile,"longitude", xdim);
	double * levs2 = get_data2(infile,"level", zdim);

	double * lat  = (double *)malloc(ydim*sizeof(double));
	double * levs = (double *)malloc(zdim*sizeof(double));
	
	for(size_t i=0;i<ydim;i++){ lat[i] = lat2[ydim-1-i];}			// reverse y-coordinate

	/********************************************************
	* Horizontal interpolation
	*********************************************************/
	double dlat = lat[3]-lat[2];
	double dlon = lon[3]-lon[2];

	horz_interpolate(topozinterp,&ITOPO(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);

	free(levs); free(lat); free(lon); free(lat2); free(levs2); 
	free(topozinterp); free(topo2);
}

/********************************************************
* Reverse the y- and z-coordinate 
*
*********************************************************/
void reverse_yz_coord(double * var, int xdim, int ydim, int zdim){

	int size = xdim*ydim*zdim;

	double * var2 = (double *)malloc(size*sizeof(double));

	/********************************************************
	* 
	*********************************************************/
	for(int i=0;i<xdim;i++){
	for(int j=0;j<ydim;j++){
	for(int k=0;k<zdim;k++){

		var2[d(i,j,k)] = var[d(i,ydim-1-j,zdim-1-k)];

	}}}

	for(int i=0;i<size;i++){ var[i] = var2[i];}

	free(var2);
}


/********************************************************
* 
*
*********************************************************/
void flip_array(double * var, int xdim, int ydim, int zdim){

	int size = xdim*ydim*zdim;

	double * var2 = (double *)malloc(size*sizeof(double));

	int half = xdim / 2;

	/********************************************************
	* 
	*********************************************************/
	for(int i=0;i<half;i++){
	for(int j=0;j<ydim;j++){
	for(int k=0;k<zdim;k++){

		var2[d(i+half,j,k)] = var[d(i,j,k)];
		var2[d(i,j,k)] = var[d(i+half,j,k)];

	}}}

	for(int i=0;i<size;i++){ var[i] = var2[i];}

	free(var2);
}



