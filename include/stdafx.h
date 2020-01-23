#pragma once

// initialization from model output will give slightly different
// results compared to running straight through because the output
// time is at the end of the time step.

// periodic boundaries for parallel and serial not compatible, don't 
// initialize serial version with output from parallel version or vice versa.

// still need to make latoffset and lonoffset consistent amoung initialization options
// can't run moisture budget without microphysics since it requires vapor
// can't run trajectory in parallel mode (array indices)
// rain total not continued for restart runs (fixed)

/*
TO DO:
Rain fall speed
python initializer
capabilities for 2d variable output so rainfall and heat fluxes are not in friction array
maybe a file naming system that is read from a file
maybe command-line argument processor
budget for serial version
link Grabowski microphysics scheme
could initialize model from a shell script, with input parameters provided
add diffusion to snow and ice
*/

// compile parallel (1) or serial (0) version
#define PARALLEL 1

#if 1

//--------------------------------------------------------
// Grid parameters
//--------------------------------------------------------
extern int NX;
extern int NY;
extern int NZ;

extern double dx;
extern double dy;
extern double dz;
extern double dt;

extern double lonoffset;
extern double latoffset;
extern int number_of_time_steps;

extern int raydampheight;
extern int outfilefreq;

//--------------------------------------------------------
// Grid stretching
//--------------------------------------------------------
extern double height_lowest_level; // height of lowest full level in meters
const int index_lowest_level = 1;		// what is the index of this level?

//--------------------------------------
// Version to execute
//--------------------------------------
#define ENSEMBLE 0
#define ENERGY 0
//--------------------------------------
// Equation set
//--------------------------------------
#define HYDROSTATIC 0			// hydrostatic option (no longer works)
#define ISLINEAR 0				// linearize equation set
#define USE_LINEAR_FRICTION 0	// linear friction in lowest model level(s)
extern int USE_TURBULENT_STRESS;// use turbulence parameterization
#define RESTING_BASIC_STATE 0	// set all basic state terms to zero (doesn't work with some options)
//--------------------------------------
// Physics options
//--------------------------------------
extern int MICROPHYSICS_OPTION;		// 0:None 1:Kessler 2:Rutludge
#define USE_TERRAIN 0				// 0:no, 1:yes
extern int SURFACE_HEAT_FLUX;		// 0:no, 1:yes
extern double WATER_TEMP_C;			// water temperature in Celsius (for surface heat fluxes)
#define USE_LANDSEA_FROM_FILE 0		// 0:no, 1:yes
//--------------------------------------
// For linearized equation set
//--------------------------------------
#define OUTPUT_DIFFUSION_TEND 0		// 
#define EXTRA_DIFFUSION 0			// only for linearized version, no grid stretching
#define FOURIER_DAMPING 0			// use Fourier transform to specify wavenumber
#define WAVE_NUMBER 3				// the specified wavenumber
//--------------------------------------
// Grid / boundaries
//--------------------------------------
extern int PERIODIC_BOUNDARIES;
#define MERIDIONAL_CROSS_SECTION 0
#define STRETCHED_GRID 1
extern int SHIFT_PRIME_MERIDIAN;
//--------------------------------------
// Initialization Options
//--------------------------------------
extern int isRestartRun;				// using the requested output file, start where it left off
extern int BASIC_STATE_OPTION;			// 0:Model output, 1:Code, 2:Reanalysis 
extern int PERTURBATION_OPTION;			// 0:Model output, 1:Code
#define REANALYSIS_INITIALIZE_OPTION 1 	// 1:processed ERA-interim 
										// 2:NCAR-reanalysis monthly mean 
										// 3:ERA5 monthly mean 
										// 4:ERA5 processed
//--------------------------------------
// Output
//--------------------------------------
extern int VERBOSE;					// lots of output to screen?
extern int PV_BUDGET;				// calculate potential vorticity budget and write to a file
extern int HEAT_BUDGET;				// calculate heat budget and write to a file
extern int MOISTURE_BUDGET;			// calculate moisture budget and write to a file
extern int VORTICITY_BUDGET;		// calculate vorticity budget and write to a file
extern int PE_BUDGET;				// calculate potential energy budget
extern int CREATE_NEW_OUTPUT_FILE;
#define PV_TRACER 0
#define OUTPUT_FRICTION_TEND true	// output frictional tendencies
#define OUTPUT_TO_FILE true			// write output to netcdf file

#define PRINT_EKE_BUDGET 0			// print EKE budget to terminal window
#define CALCULATE_OMEGA 0
#define EXTRA_OUTPUT 0
#define FILEBASE 0					// 0 - mac, 1 - savio, 2 - NERSC

#define FFTW_FLAGS FFTW_PATIENT

#endif

//--------------------------------------
// Where is the input and output?
//--------------------------------------
#if FILEBASE == 0
	#define INPUT_FILE_PATH "/Users/michaeldiaz/Desktop/model_input/"
	#define OUTPUT_FILE_PATH "/Users/michaeldiaz/Desktop/model_output/"
#elif FILEBASE == 1
	#define INPUT_FILE_PATH "../model_input/"
	#define OUTPUT_FILE_PATH "/global/scratch/michaeldiaz/"
#elif FILEBASE == 2
	#define INPUT_FILE_PATH "/global/cscratch1/sd/michaeld/"
	#define OUTPUT_FILE_PATH "/global/cscratch1/sd/michaeld/"
#endif

extern int USE_ICE;
extern int USE_MICROPHYSICS;

//---------------------------------------------------------
// Input parameters
//---------------------------------------------------------
struct input_params {

	static const int length = 300;

	int nx,ny,nz;
	double dx,dy,dz,dt;
	double corner_lat,corner_lon;
	int shift_prime_meridian;
	int time_steps,output_frequency;
	int rayleigh_damping_z;
	double height_lowest_level;
	
	int microphysics_option;
	
	int basic_state_init_option;
	int perturbation_init_option;
	
	int perturbationFileTime;
	int startOutfileAt;
	
	char basic_state_file[length];
	char perturbation_file[length];
	char output_file[length];
	char heat_budget_filename[length];
	char moisture_budget_filename[length];
	char pe_budget_filename[length];
	char pv_budget_filename[length];
	char vorticity_budget_filename[length];

	int use_surface_heat_flux;
	int turbulence_option;
	double water_temp;

	int is_restart_run;
	
	int create_new_output_file;
	
	int periodic_ew_boundaries;
	int boundary_width_north;
	int boundary_width_south;
	int boundary_width_east;	
	int boundary_width_west;	
	
	int run_heat_budget;
	int run_moisture_budget;
	int run_pe_budget;
	int run_pv_budget;
	int run_vorticity_budget;
	
	int verbose;
	
	bool has_shear;
	double hor_shear;
	double vert_shear0;
	double vert_shear1;
	
	int vortex_initialize;
	double vortex_latitude;
	double vortex_longitude;
	double vortex_radius;
	double vortex_upper_warm_anomaly;
	double vortex_lower_cold_anomaly;
	double vortex_height_wind_max;
	double vortex_vertical_extent;
	double vortex_rh_top;
	double vortex_rh_bottom;
	double vortex_rh_prime;
	double vortex_rh_radius;
	double vortex_rh_max;
};

extern struct input_params inputs;

// C RunTime Header Files
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <stdio.h>
#include <ctime>
#include <string.h>
#include <unistd.h>
#include <ctype.h>

// Model Header Files
#include "model.h"
#include "pvars.h"
#include "microphysics.h"
#include "file.h"
#include "util.h"
#include "solver.h"

#if 0
//--------------------------------------------------------
// Grid parameters
//--------------------------------------------------------
const int NX = 144;
const int NY = 53;
const int NZ = 43*2;

const double dx = 75000;
const double dy = dx;
const double dz = 250;//500;
const double dt = 450;

const double lonoffset = 60;
const double latoffset = 90;
const int number_of_time_steps = 80*24*60*60/(int)dt;

const int raydampheight = NZ - 5;
const int outfilefreq = 12*60*60/(int)dt;
const int stopheating = -576;

//--------------------------------------------------------
// Grid stretching
//--------------------------------------------------------
const double height_lowest_level = 100; // height of lowest full level in meters
const int index_lowest_level = 1;		// what is the index of this level?

//--------------------------------------
// Version to execute
//--------------------------------------
#define PARALLEL 0
#define ENSEMBLE 0
#define ENERGY 0
//--------------------------------------
// Equation set
//--------------------------------------
#define HYDROSTATIC 0			// hydrostatic option (no longer works)
#define ISLINEAR 1				// linearize equation set
#define USE_LINEAR_FRICTION 1	// linear friction in lowest model level(s)
#define USE_TURBULENT_STRESS 0	// use turbulence parameterization
#define RESTING_BASIC_STATE 0	// set all basic state terms to zero (doesn't work with some options)
//--------------------------------------
// Physics options
//--------------------------------------
#define MICROPHYSICS_OPTION 0		// 0:None 1:Kessler 2:Rutludge
#define USE_TERRAIN 0				// 0:no, 1:yes
#define SURFACE_HEAT_FLUX 0			// 0:no, 1:yes
#define USE_LANDSEA_FROM_FILE 0		// 0:no, 1:yes
//--------------------------------------
// For linearized equation set
//--------------------------------------
#define OUTPUT_DIFFUSION_TEND 0		// 
#define EXTRA_DIFFUSION 1			// only for linearized version, no grid stretching
#define FOURIER_DAMPING 1			// use Fourier transform to specify wavenumber
#define WAVE_NUMBER 4				// the specified wavenumber
//--------------------------------------
// Grid / boundaries
//--------------------------------------
#define PERIODIC_BOUNDARIES 1
#define MERIDIONAL_CROSS_SECTION 0
#define STRETCHED_GRID 0
#define SHIFT_PRIME_MERIDIAN 0
//--------------------------------------
// Initialization Options
//--------------------------------------
#define BASIC_STATE_OPTION 1			// 0:Model output, 1:Code, 2:Reanalysis 
#define PERTURBATION_OPTION 1			// 0:Model output, 1:Code
#define REANALYSIS_INITIALIZE_OPTION 1 	// 1:processed ERA-interim 
										// 2:NCAR-reanalysis monthly mean 
										// 3:ERA5 monthly mean 
										// 4:ERA5 processed
//--------------------------------------
// Output
//--------------------------------------
#define VERBOSE 1					// lots of output to screen?
#define PV_BUDGET 0					// calculate potential vorticity budget and write to a file
#define HEAT_BUDGET 0				// calculate heat budget and write to a file
#define MOISTURE_BUDGET 0			// calculate moisture budget and write to a file
#define PE_BUDGET 1					// calculate potential energy budget
#define VORTICITY_BUDGET 0			// calculate vorticity budget and write to a file
#define PV_TRACER 0
#define OUTPUT_FRICTION_TEND true	// output frictional tendencies
#define OUTPUT_TO_FILE true			// write output to netcdf file
#define CREATE_NEW_OUTPUT_FILE true
#define PRINT_EKE_BUDGET 0			// print EKE budget to terminal window
#define CALCULATE_OMEGA 0
#define EXTRA_OUTPUT 0
#define FILEBASE 0					// 0 - mac, 1 - savio, 2 - NERSC

#define FFTW_FLAGS FFTW_PATIENT

#endif
