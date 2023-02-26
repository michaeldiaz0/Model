#define PARALLEL_IO false

#include "netcdf.h"

#if PARALLEL_IO
	#include "netcdf_par.h"
#endif

#if PARALLEL
    #include "mpi.h"
#endif

#define NC64BIT true	// whether to support 64 bit offsets (allows bigger files ~ > 2GB for each output time )

#define MODEL_FILES_NOT_WRITTEN 0
#define MODEL_FILES_WRITTEN 1
#define ALL_FILES_WRITTEN 2
#define CAN_RESTART_FROM_FILE_END ALL_FILES_WRITTEN


const int len = 150;

/*****************************************************************************
* INPUT FILES FOR INITIALIZATION AND OUTPUT
******************************************************************************/
//------------------------------------------------------------------------------
// from NCEP reanalysis
//------------------------------------------------------------------------------
const char geoheight_file[len] = INPUT_FILE_PATH"jul2014_ghgt.nc";
const char uwind_file[len] = INPUT_FILE_PATH"jul2014_uwnd.nc";
const char vwind_file[len] = INPUT_FILE_PATH"jul2014_vwnd.nc";
const char temp_file[len] = INPUT_FILE_PATH"jul2014_temp.nc";
const char topo_file[len] = INPUT_FILE_PATH"hgt.nc";
//------------------------------------------------------------------------------
// from NETCDF files converted from ERA-interim grib format
//------------------------------------------------------------------------------
const char datafile[len] = INPUT_FILE_PATH"outEOF4.nc";//"era5basic279.nc";
const char datafile2[len] = INPUT_FILE_PATH"outEOF4.nc";
const char landseaMaskFile[len] = INPUT_FILE_PATH"lsmask.oisst.v2.nc";
//------------------------------------------------------------------------------
// Perform the energy budget calculation on this file
//------------------------------------------------------------------------------
const char defaultEnergyBudgetFileName[len] = OUTPUT_FILE_PATH"outfile0417_9.nc";
//------------------------------------------------------------------------------
// Input files
//------------------------------------------------------------------------------
//const char inputPerturbationFileName[len] = OUTPUT_FILE_PATH"outfile1013_0.nc";//"outfile1207_2.nc";//"outfile0418_2.nc";
//const char inputBasicStateFileName[len] = OUTPUT_FILE_PATH"outfile0922_0.nc";//"outfile0317_1.nc";//"outfile0305_2.nc";
extern int perturbationFileTime;// = 7;//344;//12;//58;//344;//58;//344;	// what file output time to read in
extern int startOutfileAt;// = -3;

//------------------------------------------------------------------------------
// Output files
//------------------------------------------------------------------------------
//const char filename[len] = OUTPUT_FILE_PATH"outfile0417_10.nc";
//const char out_filename[len] = OUTPUT_FILE_PATH"outfile1106_0.nc";
//const char out_heat_budget_filename[len] = OUTPUT_FILE_PATH"temp_budget1106_0.nc";
//const char out_moisture_budget_filename[len] = OUTPUT_FILE_PATH"moisture_budget1106_0.nc";
//const char out_pe_budget_filename[len] = OUTPUT_FILE_PATH"pe_budget1106_0.nc";
//const char vorticity_budget_filename[len] = OUTPUT_FILE_PATH"vort_budget1028_0.nc";
//const char out_pv_budget_filename[len] = OUTPUT_FILE_PATH"pv_budget1106_0.nc";
const char pv_tracer_filename[len] = OUTPUT_FILE_PATH"pv_tracer0205_5.nc";

extern char inputPerturbationFileName[len];
extern char inputBasicStateFileName[len];

//extern char out_filename[len];
extern char filename[len];
extern char energyBudgetFileName[len];

/********************************************************
* VARIABLE NAMES AND DESCRIPTIONS
*********************************************************/
#if EXTRA_OUTPUT
const int mainVarCount = 6;	// number of primary variables
#else
const int mainVarCount = 5;	// number of primary variables
#endif

const int baseVarCount = 6;
const int mbaseVarCount = 4;
const int coordVarCount = 2;
//const int twodVarCount = 2;

const char var_names[11][10] = {"u-wind","v-wind","w-wind", "theta", "pi","rate"};
const char mp_var_names[11][10] = {"qv", "qc", "qr","qi","qs","qg"};
const char opt_var_names[11][10] = {"fric","diff"};
const char basevar_names[6][10] = {"ubar","vbar","wbar","thbar","qbar","pbar"};
const char modelbasevar_names[4][10] = {"tb","pib","qb","zu"};
const char coordvar_names[2][10] = {"lat","lon"};
const char topo_name[5] = "topo";
const char lsmask_name[10] = "lsmask";
//const char var_names_2d[11][10] = {"rainfall","LHF","int_fric"};

extern int isRestartRun;

void add_commandline_arguments(char *filename0,const char *filename1);
void set_outfilename(char *myfilename);
bool does_file_exist(const char *myfilename);
bool has_attribute(const char *myfilename,const char *att);
void set_attribute(const char *myfilename,const char *att,int value);
void file_output_status(int a);
int get_file_output_status();

/********************************************************
* FILE WRITING PROCEDURES
*********************************************************/
void output_meteorological_fields_to_file(void (*write2d)(const char*,const char*,double*,size_t),void (*write3d)(const char*,const char*,double*,size_t),int);
void create_outfile(const char *filename, bool basestate,bool modelBaseState,bool coordinates);
void write_mbvar_to_file(int ncid,const char *var_name,double *var);
void write_pvar_to_file(int ncid,const char *var_name,double *var,size_t tcount);
void write_pvar_to_file(const char *filename,const char *var_name,double *var,size_t tcount);
void write_pvar_to_file(const char *myfilename,const char *var_name,double *var,int nx,int ny,size_t tcount);
void write_pvar_to_file_2d(const char *filename,const char *var_name,double *var,size_t tcount);
void write_bvar_to_file(int ncid,const char *var_name,double *var);
void parallel_write_pvar_to_file_2d(const char *myfilename,const char *var_name,double *var,size_t tcount);
void parallel_write_pvar_to_file_3d(const char *myfilename,const char *var_name,double *var,size_t tcount);
void write_all_pvars(const char *filename,size_t tcount);
void write_all_bvars(const char *filename);
void write_all_mbvars(const char *filename);
void write_all_coordinates(const char *filename);
void write_topo_to_file(const char *filename);
void outfile_init(char *myfilename);
void write_time_to_file(const char *myfilename,size_t tcount);

/********************************************************
* FILE READING PROCEDURES
*********************************************************/
double* get_data2(const char*,const char*,int);
void get_data(const char*,const char *, int, double*);
void get_data(const char*,const char *, int, float*);
void get_data_at_time(const char*,const char*,size_t,double*);
void get_data_at_time(const char*,const char*, size_t, double *,int,int,int);
void get_dims(const char*,const char*,const char*,const char*,size_t[3]);
void get_model_data(const char*,size_t);
void get_grid_spacing(const char *myfilename,const char *name_dx,const char *name_dy,const char *name_dz,float grid_spacing[3]);
size_t get_file_length(const char *myfilename,const char *varname);

/********************************************************
* FILE ERROR HANDLING PROCEDURES
*********************************************************/
bool fileHasVar(const char *myfilename,const char *varname);
void handle_error(int);
