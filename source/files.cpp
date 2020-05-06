#include "stdafx.h"

void handle_file_error(int,const char*);
void add_attribute(int ncid,const char *att_name,const char *value);
void add_attribute(int ncid,const char *att_name,int *values,int count);
void add_attribute(int ncid,const char *att_name,int value);
void add_attribute(int ncid,const char *att_name,double value);

double * temp_storage;

char filename[len];
char inputPerturbationFileName[len];
char inputBasicStateFileName[len];

/********************************************************
* 
*********************************************************/
bool does_file_exist(const char *myfilename){
	
	int ncid,status;

	status = nc_open(myfilename, 0, &ncid);
	
	//printf("nc status = %d",status);
	
	if(status == NC_ENOTFOUND || status == 2){ return false;}
	else {
	
		if (status != NC_NOERR) handle_error(status);
	
		status = nc_close(ncid);
		if (status != NC_NOERR) handle_error(status);
		
		return true;
	}
}

/********************************************************
* Create an output file for the model data
* 
* myfilename			- name of file
* basestate 		- whether or not to define environmental base state variables
* modelBaseState 	- whether or not to define model base state variables
* coordinates		- whether or not to define coordinate variables
*********************************************************/
void create_outfile(const char *myfilename, bool basestate,bool modelBaseState,bool coordinates){

	int status;
	int ncid;

	int xID, yID, zID, tID;
	
	int var_dimids[4];
	int var_id;

	/*************************************************
	* Create file
	**************************************************/
	if(NC64BIT){
		status = nc_create(myfilename,NC_CLOBBER|NC_64BIT_OFFSET, &ncid);
	} else {
		#if !PARALLEL_IO
			status = nc_create(myfilename,NC_CLOBBER, &ncid);
		#else
			status = nc_create_par(myfilename,NC_NETCDF4 | NC_CLOBBER | NC_MPIIO, MPI_COMM_WORLD,MPI_INFO_NULL,&ncid);
		#endif
	}
	
	if (status != NC_NOERR) handle_error(status);

	//printf("nc id = %d\n",ncid);

	/*************************************************
	* Add global attributes
	**************************************************/
	add_attribute(ncid,"dx",dx);
	add_attribute(ncid,"dy",dy);
	add_attribute(ncid,"dz",dz);
	add_attribute(ncid,"dt",dt);

	add_attribute(ncid,"latoffset",latoffset);
	add_attribute(ncid,"lonoffset",lonoffset);

	add_attribute(ncid,"reanalysis_input_data",datafile);
	add_attribute(ncid,"perturbation_input_data",inputPerturbationFileName);
	add_attribute(ncid,"basic_state_input_data",inputBasicStateFileName);

	add_attribute(ncid,"hydrostatic",HYDROSTATIC);
	add_attribute(ncid,"islinear",ISLINEAR);
	add_attribute(ncid,"use_linear_friction",USE_LINEAR_FRICTION);
	add_attribute(ncid,"use_turbulent_stress",USE_TURBULENT_STRESS);
	add_attribute(ncid,"microphysics_option",MICROPHYSICS_OPTION);
	add_attribute(ncid,"use_terrain",USE_TERRAIN);
	add_attribute(ncid,"surface_heat_flux",SURFACE_HEAT_FLUX);	 
	add_attribute(ncid,"water_temperature",WATER_TEMP_C);
	add_attribute(ncid,"use_landsea_from_file",USE_LANDSEA_FROM_FILE);
	add_attribute(ncid,"extra_diffusion",EXTRA_DIFFUSION);
	add_attribute(ncid,"use_linear_friction",USE_LINEAR_FRICTION);
	add_attribute(ncid,"periodic_boundaries",PERIODIC_BOUNDARIES);
	add_attribute(ncid,"use_linear_friction",USE_LINEAR_FRICTION);
	add_attribute(ncid,"stretched_grid",STRETCHED_GRID);	 
	add_attribute(ncid,"height_lowest_level",height_lowest_level);
	add_attribute(ncid,"index_lowest_level",index_lowest_level);
	add_attribute(ncid,"basic_state_option",BASIC_STATE_OPTION);
	add_attribute(ncid,"perturbation_option",PERTURBATION_OPTION);
	add_attribute(ncid,"reanalysis_initialize_option",REANALYSIS_INITIALIZE_OPTION);
	add_attribute(ncid,"last_file_successfully_written",ALL_FILES_WRITTEN);
	 
	/*************************************************
	* Create dimensions
	**************************************************/
	status = nc_def_dim(ncid, "x", (size_t)NX, &xID);
	if (status != NC_NOERR) handle_error(status);

	status = nc_def_dim(ncid, "y", (size_t)NY, &yID);
	if (status != NC_NOERR) handle_error(status);

	status = nc_def_dim(ncid, "z", (size_t)NZ, &zID);
	if (status != NC_NOERR) handle_error(status);

	status = nc_def_dim(ncid, "t", NC_UNLIMITED, &tID);
	if (status != NC_NOERR) handle_error(status);

	/*************************************************
	* Define forecasted variables
	**************************************************/
	var_dimids[0] = tID; var_dimids[1] = xID; var_dimids[2] = yID; var_dimids[3] = zID;
	
	int pvarcount = mainVarCount;

	/*************************************************
	* Wind and temperature
	**************************************************/
	for(int i=0;i<pvarcount;i++){

		status = nc_def_var (ncid,var_names[i], NC_FLOAT, 4, var_dimids, &var_id);
		if (status != NC_NOERR) handle_error(status);
	}

	//-------------------------------------------------
	// Two-dimensional fields
	//-------------------------------------------------
	//for(int i=0;i<twodVarCount;i++){

		//status = nc_def_var (ncid,var_names_2d[i], NC_FLOAT, 3, var_dimids, &var_id);
		//if (status != NC_NOERR) handle_error(status);
		//}

	//-------------------------------------------------
	// Turbulance parameterization fields
	//-------------------------------------------------	
	if(USE_TURBULENT_STRESS){
		
		status = nc_def_var (ncid,"int_fric", NC_FLOAT, 3, var_dimids, &var_id);
		if (status != NC_NOERR) handle_error(status);
		
		if(SURFACE_HEAT_FLUX){
		
			status = nc_def_var (ncid,"lh_flux", NC_FLOAT, 3, var_dimids, &var_id);
			if (status != NC_NOERR) handle_error(status);
		}
	}
	
	/*************************************************
	* Microphysics
	**************************************************/
	if(USE_MICROPHYSICS){
	
		int mpVarCount = 3;
	
		status = nc_def_var (ncid,"rainfall", NC_FLOAT, 3, var_dimids, &var_id);
		if (status != NC_NOERR) handle_error(status);
	
		if(USE_ICE){
			
			mpVarCount = 5;
			
			status = nc_def_var (ncid,"snowfall", NC_FLOAT, 3, var_dimids, &var_id);
			if (status != NC_NOERR) handle_error(status);
		}
	
		for(int i=0;i<mpVarCount;i++){

			status = nc_def_var (ncid,mp_var_names[i], NC_FLOAT, 4, var_dimids, &var_id);
			if (status != NC_NOERR) handle_error(status);
		}
	}
	
	/*************************************************
	* Optional (friction,diffusion)
	**************************************************/
	if(OUTPUT_FRICTION_TEND){ 
	
		status = nc_def_var (ncid,opt_var_names[0], NC_FLOAT, 4, var_dimids, &var_id);
		if (status != NC_NOERR) handle_error(status);
	}
	
	if(OUTPUT_DIFFUSION_TEND){
	
		status = nc_def_var (ncid,opt_var_names[1], NC_FLOAT, 4, var_dimids, &var_id);
		if (status != NC_NOERR) handle_error(status);
	}

	/*************************************************
	* Define environmental base state variables
	**************************************************/
	if(basestate){

		var_dimids[0] = xID; var_dimids[1] = yID; var_dimids[2] = zID;

		for(int i=0;i<baseVarCount;i++){

			status = nc_def_var (ncid,basevar_names[i], NC_FLOAT, 3, var_dimids, &var_id);
			if (status != NC_NOERR) handle_error(status);
		}
	}

	/*************************************************
	* Define model base state variables
	**************************************************/
	if(modelBaseState){

		var_dimids[0] = zID;

		for(int i=0;i<mbaseVarCount;i++){

			status = nc_def_var (ncid,modelbasevar_names[i], NC_FLOAT, 1, var_dimids, &var_id);
			if (status != NC_NOERR) handle_error(status);
		}
	}

	/*************************************************
	* Define coordinate arrays
	**************************************************/
	if(coordinates){

		var_dimids[0] = yID;

		status = nc_def_var (ncid,coordvar_names[0], NC_FLOAT, 1, var_dimids, &var_id);
		if (status != NC_NOERR) handle_error(status);

		var_dimids[0] = xID;

		status = nc_def_var (ncid,coordvar_names[1], NC_FLOAT, 1, var_dimids, &var_id);
		if (status != NC_NOERR) handle_error(status);

		var_dimids[0] = tID;

		status = nc_def_var (ncid,"time", NC_FLOAT, 1, var_dimids, &var_id);
		if (status != NC_NOERR) handle_error(status);

	}

	/*************************************************
	* Define topography
	**************************************************/
	var_dimids[0] = xID; var_dimids[1] = yID;

	status = nc_def_var (ncid,topo_name, NC_FLOAT, 2, var_dimids, &var_id);
	if (status != NC_NOERR) handle_error(status);

	/*************************************************
	* Close file
	**************************************************/
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);

}

/********************************************************
*
*********************************************************/
void add_attribute(int ncid,const char *att_name,double value){

	int status = nc_put_att_double (ncid, NC_GLOBAL, att_name, NC_FLOAT, 1, &value);
	if (status != NC_NOERR) handle_error(status);
}

/********************************************************
*
*********************************************************/
void add_attribute(int ncid,const char *att_name,int value){

	int status = nc_put_att_int (ncid, NC_GLOBAL, att_name, NC_INT, 1, &value);
	if (status != NC_NOERR) handle_error(status);
}

/********************************************************
*
*********************************************************/
void add_attribute(int ncid,const char *att_name,int *values,int count){

	int status = nc_put_att_int (ncid, NC_GLOBAL, att_name, NC_INT, count, values);
	if (status != NC_NOERR) handle_error(status);
}

/********************************************************
*
*********************************************************/
void add_attribute(int ncid,const char *att_name,const char *value){

	int status = nc_put_att_text (ncid, NC_GLOBAL, att_name,strlen(value), value);
	if (status != NC_NOERR) handle_error(status);
}

/****************************************************************************
*
* 				SINGLE VARIABLE WRITING PROCEDURES
*
*****************************************************************************/

/********************************************************
* Write data for a single model base state variable to
* the netcdf file
* 
* ncid 		- netcdf file ID
* var_name 	- variable name
* var[NZ]	- variable data
*********************************************************/
void write_mbvar_to_file(int ncid,const char *var_name,double *var){

	int status;
	int var_id;

	status = nc_inq_varid (ncid, var_name, &var_id);
	if (status != NC_NOERR) handle_error(status);

	printf("writing %s %d to file\n",var_name,var_id);

	status = nc_put_var_double(ncid, var_id, &var[0]);
	if (status != NC_NOERR) handle_error(status);

}

/********************************************************
* Write data for single perturbation variable at a single 
* time to the netcdf file
* 
* ncid 				- netcdf file ID
* var_name 			- variable name
* var[NX][NY][NZ]	- variable data
* tcount			- file time variable
*********************************************************/
void write_pvar_to_file(int ncid,const char *var_name,double *var,size_t tcount){

	int status;
	int var_id;

	size_t start[] = {tcount,0,0,0};
	size_t count[] = {1,NX,NY,NZ};

	status = nc_inq_varid (ncid, var_name, &var_id);
	if (status != NC_NOERR) handle_error(status);

	printf("writing %s %d to file\n",var_name,var_id);

	status = nc_put_vara_double(ncid, var_id,start,count, var);
	if (status != NC_NOERR) handle_error(status);
	
}

/********************************************************
* Write data for single perturbation variable at a single 
* time to the netcdf file
* 
* myfilename 		- name of file
* var_name 			- variable name
* var				- variable data
* tcount			- file time variable
*********************************************************/
void write_pvar_to_file(const char *myfilename,const char *var_name,double *var,size_t tcount){

	int ncid,status;

	status = nc_open(myfilename, NC_WRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	write_pvar_to_file(ncid,var_name,var,tcount);
	
	/*************************************************
	* Close file
	**************************************************/
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);
}

/********************************************************
* Write data for single perturbation variable at a single 
* time to the netcdf file
* 
* myfilename 		- name of file
* var_name 			- variable name
* var				- variable data
* tcount			- file time variable
*********************************************************/
void write_pvar_to_file_2d(const char *myfilename,const char *var_name,double *var,size_t tcount){

	int ncid,status,var_id;

	status = nc_open(myfilename, NC_WRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	size_t start[] = {tcount,0,0};
	size_t count[] = {1,NX,NY};

	status = nc_inq_varid (ncid, var_name, &var_id);
	if (status != NC_NOERR) handle_error(status);

	printf("writing %s %d to file\n",var_name,var_id);

	status = nc_put_vara_double(ncid, var_id,start,count, var);
	if (status != NC_NOERR) handle_error(status);
	
	/*************************************************
	* Close file
	**************************************************/
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);
}

#if !PARALLEL_IO
/********************************************************
* Write data for single perturbation variable at a single 
* time to the netcdf file for parallel version.
* 
* myfilename 		- name of file
* var_name 			- variable name
* var				- variable data
* tcount			- file time variable
*********************************************************/
void parallel_write_pvar_to_file_3d(const char *myfilename,const char *var_name,double *var,size_t tcount){
	
	gatherArrays_3d(var,output_to_file_3d);
	
	if(rank==0){ write_pvar_to_file(myfilename,var_name,output_to_file_3d,tcount);}
	
}

/********************************************************
* Write data for single perturbation variable at a single 
* time to the netcdf file for parallel version.
* 
* myfilename 		- name of file
* var_name 			- variable name
* var				- variable data
* tcount			- file time variable
*********************************************************/
void parallel_write_pvar_to_file_2d(const char *myfilename,const char *var_name,double *var,size_t tcount){
	
	gatherArrays_2d(var,output_to_file_2d);
	
	if(rank==0){ write_pvar_to_file_2d(myfilename,var_name,output_to_file_2d,tcount);}
	
}
#else
/********************************************************
* 
*********************************************************/
void parallel_write_pvar_to_file(const char *myfilename,const char *var_name,double *var,size_t tcount){

	for(int i=3;i<fNX-3;i++){
	for(int j=3;j<fNY-3;j++){
	for(int k=0;k<fNZ;k++){
		
		temp_storage[(i-3)*myNY*myNZ + (j-3)*myNZ + (k)] = var[INDEX(i,j,k)];
	}}}

	int status,res;
	int var_id,ncid;

	size_t start[] = {tcount,ibs[rank],jbs[rank],0};
	size_t count[] = {1,myNX,myNY,myNZ};
	ptrdiff_t stride[] = {myNY*myNZ,myNZ,1};
	
	status = nc_open_par(myfilename, NC_WRITE|NC_MPIIO, MPI_COMM_WORLD,MPI_INFO_NULL,&ncid);
	if (status != NC_NOERR) handle_error(status);
	
	status = nc_inq_varid (ncid, var_name, &var_id);
	if (status != NC_NOERR) handle_error(status);

	status = nc_var_par_access(ncid, var_id, NC_COLLECTIVE);
	if (status != NC_NOERR) handle_error(status);
	
	if(rank==0){ printf("writing %s %d to file\n",var_name,var_id);}

	status = nc_put_vara_double(ncid, var_id,start,count, temp_storage);
	if (status != NC_NOERR) handle_error(status);

	/*************************************************
	* Close file
	**************************************************/
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);

}
#endif
/********************************************************
* Write data for single perturbation variable at a single 
* time to the netcdf file
* 
* myfilename 		- name of file
* var_name 			- variable name
* var[NX][NY][NZ]	- variable data
* tcount			- file time variable
*********************************************************/
void write_pvar_to_file(const char *myfilename,const char *var_name,double *var,int nx,int ny,size_t tcount){

	int status;
	int var_id;
	int ncid;

	size_t start[] = {tcount,0,0,0};
	size_t count[] = {1,(size_t)nx,(size_t)ny,(size_t)NZ};

	status = nc_open(myfilename, NC_WRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid (ncid, var_name, &var_id);
	if (status != NC_NOERR) handle_error(status);

	printf("writing %s %d to file\n",var_name,var_id);

	status = nc_put_vara_double(ncid, var_id,start,count, var);
	if (status != NC_NOERR) handle_error(status);
	
	/*************************************************
	* Close file
	**************************************************/
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);
}


/********************************************************
* Write data for a single environmental base state variable
* 
* ncid 				- netcdf file ID
* var_name 			- variable name
* var[NX][NY][NZ]	- variable data
*********************************************************/
void write_bvar_to_file(int ncid,const char *var_name,double *var){

	int status;
	int var_id;

	status = nc_inq_varid (ncid, var_name, &var_id);
	if (status != NC_NOERR) handle_error(status);

	printf("writing %s %d to file\n",var_name,var_id);

	status = nc_put_var_double(ncid, var_id, var);
	if (status != NC_NOERR) handle_error(status);

}

/********************************************************
* Write data for a single model base state variable to
* the netcdf file
* 
* ncid 		- netcdf file ID
* var_name 	- variable name
* var[NZ]	- variable data
*********************************************************/
void write_topo_to_file(const char *myfilename){

	int status;
	int var_id;
	int ncid;

	double outtopo[NX][NY];

	for(int i=0;i<NX;i++)
		for(int j=0;j<NY;j++)
			outtopo[i][j] = ITOPO(i,j,NZ/2);

	/********************************************************
	* Open file, get file ID, write data
	*********************************************************/
	status = nc_open(myfilename, NC_WRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid (ncid, topo_name, &var_id);
	if (status != NC_NOERR) handle_error(status);

	printf("writing %s %d to file\n",topo_name,var_id);

	status = nc_put_var_double(ncid, var_id, &outtopo[0][0]);
	if (status != NC_NOERR) handle_error(status);

	/*************************************************
	* Close file
	**************************************************/
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);

}

/********************************************************
* 
*********************************************************/
void write_time_to_file(const char *myfilename,size_t tcount){

	int ncid;
	int status;
	int var_id;
	double hours;
	size_t index[] = {tcount};

	/********************************************************
	* Open file, get file ID, write data
	*********************************************************/
	status = nc_open(myfilename, NC_WRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	/********************************************************
	* Write time variable
	*********************************************************/
	hours = mtime / 3600;

	status = nc_inq_varid (ncid,"time", &var_id);
	if (status != NC_NOERR) handle_error(status);

	status = nc_put_var1_double(ncid, var_id, index, &hours);
	if (status != NC_NOERR) handle_error(status);

	/*************************************************
	* Close file
	**************************************************/
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);

}

/****************************************************************************
*
* 				MULTIPLE VARIABLE WRITING PROCEDURES
*
*****************************************************************************/

/*********************************************************************
* Write all variables to the netcdf file
*
* @param write2d - function to write 2d data
* @param write3d - function to write 3d data
**********************************************************************/
void output_meteorological_fields_to_file(
		void (*write2d)(const char*,const char*,double*,size_t),
		void (*write3d)(const char*,const char*,double*,size_t),
		int tcount
		){
	//------------------------------------------------------
	// Core meteorological field
	//------------------------------------------------------
	write3d(filename,"u-wind",us, tcount);
	write3d(filename,"v-wind",vs, tcount);
	write3d(filename,"w-wind",ws, tcount);
	write3d(filename,"pi",    pis,tcount);
	write3d(filename,"theta", ths,tcount);

	//------------------------------------------------------
	// Microphysics specfic variables
	//------------------------------------------------------
	if(USE_MICROPHYSICS){
		
		write3d(filename,"qv",qvs,tcount);
		write3d(filename,"qc",qcs,tcount);
		write3d(filename,"qr",qrs,tcount);
		
		write2d(filename,"rainfall",accRain,tcount);
		
		if(USE_ICE){
			
			write2d(filename,"snowfall",accSnow,tcount);
			
			write3d(filename,"qi",qis,tcount);
			write3d(filename,"qs",qss,tcount);
		}
	}
	//------------------------------------------------------
	// Turbulence specfic variables
	//------------------------------------------------------	
	if(USE_TURBULENT_STRESS){
		
		write2d(filename,"int_fric",integrated_friction_ke,tcount);
		
		if(SURFACE_HEAT_FLUX){
			write2d(filename,"lh_flux",latent_heat_flux,tcount);
		}
	}
	//------------------------------------------------------
	// Extra variables
	//------------------------------------------------------	
	if(EXTRA_OUTPUT){
		write3d(filename,"rate",rate,tcount);
	}

	if(OUTPUT_FRICTION_TEND){
		write3d(filename,"fric",frictions,tcount);
	}
}
#if 0
/********************************************************
* Write data for all perturbation variables to file
* 
* filename	- name of netcdf file
* tcount	- file time variable
*********************************************************/
void write_all_pvars(const char *myfilename,size_t tcount){

	int ncid;
	int status;
	int var_id;
	double hours;
	size_t index[] = {tcount};

	/********************************************************
	* Open file, get file ID, write data
	*********************************************************/
	status = nc_open(myfilename, NC_WRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	write_pvar_to_file(ncid,var_names[0],us,tcount);
	write_pvar_to_file(ncid,var_names[1],vs,tcount);
	write_pvar_to_file(ncid,var_names[2],ws,tcount);
	write_pvar_to_file(ncid,var_names[3],ths,tcount);
	write_pvar_to_file(ncid,var_names[4],pis,tcount);

	if(USE_MICROPHYSICS){

		write_pvar_to_file(ncid,mp_var_names[0],qvs,tcount);
		write_pvar_to_file(ncid,mp_var_names[1],qcs,tcount);
		write_pvar_to_file(ncid,mp_var_names[2],qrs,tcount);
		
		if(USE_ICE){
			write_pvar_to_file(ncid,mp_var_names[3],qis,tcount);
			write_pvar_to_file(ncid,mp_var_names[4],qss,tcount);		
		}
	}

	/********************************************************
	* Write time variable
	*********************************************************/
	hours = bigcounter * dt / 3600;

	status = nc_inq_varid (ncid,"time", &var_id);
	if (status != NC_NOERR) handle_error(status);

	status = nc_put_var1_double(ncid, var_id, index, &hours);
	if (status != NC_NOERR) handle_error(status);

	/*************************************************
	* Close file
	**************************************************/
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);

}
#endif
/********************************************************
* Write data for all environmental base state variables 
* to the netcdf file
* 
* filename	- name of netcdf file
*********************************************************/
void write_all_bvars(const char *myfilename){

	int ncid;
	int status;

	/********************************************************
	* Open file, get variable ID
	*********************************************************/
	status = nc_open(myfilename, NC_WRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	write_bvar_to_file(ncid,basevar_names[0],&IUBAR(0,0,0));
	write_bvar_to_file(ncid,basevar_names[1],&IVBAR(0,0,0));
	write_bvar_to_file(ncid,basevar_names[2],&IWBAR(0,0,0));
	write_bvar_to_file(ncid,basevar_names[3],&ITHBAR(0,0,0));
	write_bvar_to_file(ncid,basevar_names[4],&IQBAR(0,0,0));
	write_bvar_to_file(ncid,basevar_names[5],&IPBAR(0,0,0));

	/*************************************************
	* Close file
	**************************************************/
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);

}

/********************************************************
* Write data for all model base state variables 
* to the netcdf file
* 
* filename	- name of netcdf file
*********************************************************/
void write_all_mbvars(const char *myfilename){

	int ncid;
	int status;

	/********************************************************
	* Open file, get variable ID
	*********************************************************/
	status = nc_open(myfilename, NC_WRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	write_mbvar_to_file(ncid,modelbasevar_names[0],tb);
	write_mbvar_to_file(ncid,modelbasevar_names[1],pib);
	write_mbvar_to_file(ncid,modelbasevar_names[2],qb);
	
	if(STRETCHED_GRID){
		write_mbvar_to_file(ncid,modelbasevar_names[3],zsu);
	} else {
		write_mbvar_to_file(ncid,modelbasevar_names[3],zu);
	}
	
	/*************************************************
	* Close file
	**************************************************/
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);

}

/********************************************************
* Write data for all model base state variables 
* to the netcdf file
* 
* filename	- name of netcdf file
*********************************************************/
void write_all_coordinates(const char *myfilename){

	int ncid;
	int status;
	int var_id;

	/********************************************************
	* Open file, get variable ID
	*********************************************************/
	status = nc_open(myfilename, NC_WRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid (ncid,coordvar_names[0], &var_id);
	if (status != NC_NOERR) handle_error(status);

	status = nc_put_var_double(ncid, var_id, &outLats[0]);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid (ncid,coordvar_names[1], &var_id);
	if (status != NC_NOERR) handle_error(status);

	status = nc_put_var_double(ncid, var_id, &outLons[0]);
	if (status != NC_NOERR) handle_error(status);

	/*************************************************
	* Close file
	**************************************************/
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);

}

/********************************************************
* 
* 
* filename	- name of netcdf file
*********************************************************/
size_t get_file_length(const char *myfilename,const char *varname){
	
	int status,ncid,timeid;
	size_t length;
	
	status = nc_open(myfilename, 0, &ncid);
	handle_file_error(status,myfilename);
	
	status = nc_inq_dimid(ncid, varname, &timeid);
	if (status != NC_NOERR) handle_error(status);
	
	status = nc_inq_dimlen(ncid, timeid, &length);
	if (status != NC_NOERR) handle_error(status);
	
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);
	
	return length;
}

/********************************************************
* 
* filename	- name of netcdf file
*********************************************************/
void handle_output_to_existing_file(const char *myfilename){
		
	//-----------------------------------------------
	// Check whether file dimensions matched
	//-----------------------------------------------
	size_t x,y,z,t;

	size_t dims[3];
	float grid_spacing[3];

	get_dims(myfilename,"x","y","z",dims);
	get_grid_spacing(myfilename,"dx","dy","dz",grid_spacing);

	x = dims[0]; y = dims[1]; z = dims[2];

	if(x != NX || y != NY || z != NZ){
	
		printf("Error: dimensions of file %s (%lu %lu %lu) don't match of current configuration (%d %d %d)!\n",myfilename,x,y,z,NX,NY,NZ);
		printf("Exiting at line %d in file %s\n",__LINE__,__FILE__);		
		exit(0);
	}
	
	if(dx != grid_spacing[0] || dy != grid_spacing[1] || dz != grid_spacing[2]){
	
		printf("Error: dimensions of file %s (%lu %lu %lu) don't match of current configuration (%d %d %d)!\n",myfilename,x,y,z,NX,NY,NZ);
		printf("Exiting at line %d in file %s\n",__LINE__,__FILE__);		
		exit(0);
	}

	int myStartOutfileAt;

	//-------------------------------------------------------
	// Determine starting index for input perturbation file
	//--------------------------------------------------------
	if(isRestartRun){
		
		int filestatus = get_file_output_status();	
		//-----------------------------------------------
		// Restart from the last time in the file
		//-----------------------------------------------
		if(filestatus == CAN_RESTART_FROM_FILE_END || filestatus == -1){
			
			myStartOutfileAt = -1;
		//-----------------------------------------------
		// Last time was incomplete, start from t - 2
		//-----------------------------------------------			
		} else {
			
			myStartOutfileAt = -2;
			
			printf("----------------------------------------------------\n");
			printf("Last output file was incomplete. Starting at time index t - 2\n");
			printf("----------------------------------------------------\n");			
		}
	//-----------------------------------------------
	// Just use the value specified in header file
	//-----------------------------------------------		
	} else { myStartOutfileAt = startOutfileAt;}

	//-----------------------------------------------
	// Open file
	//-----------------------------------------------
	int status,ncid;
			
	status = nc_open(myfilename, NC_WRITE, &ncid);
	handle_file_error(status,myfilename);

	//-----------------------------------------------
	// Check time dimension
	//-----------------------------------------------
	t = get_file_length(myfilename,"t");

	if(myStartOutfileAt < 0){
		
		if((int)t + myStartOutfileAt < 0){ printf("Error: file time index requested (%lu) is invalid!\n",t + myStartOutfileAt);}
		
		file_time_counter = t + myStartOutfileAt;
	}
	else { file_time_counter = myStartOutfileAt; }

	if(file_time_counter > t+1){
	
		printf("Error: time dimension of current configuration (%d) exceeds that of output file (%lu) %s !\n",myStartOutfileAt,t+1,myfilename);			
		exit(0);			
	}
#if 0
	status = nc_redef(ncid);
	if (status != NC_NOERR) handle_error(status);
	
	//------------------------------------------------------
	// Add restart attributes to the file
	//------------------------------------------------------
	int *prior_restart_times;
	int type;
	size_t length;

	status = nc_inq_att(ncid, NC_GLOBAL, "restart_times", &type,&length);
	
	if (status == NC_NOERR){ prior_restart_times = ((int*) calloc(length+1,sizeof(int)));}

	status = nc_get_att_int(ncid, NC_GLOBAL, "restart_times", prior_restart_times);

	//-----------------------------------------------
	// Attribute doesn't exist, create it
	//-----------------------------------------------
	if(status==NC_ENOTATT){ 
		//printf("here!\n");
		add_attribute(ncid,"restart_times",(int)file_time_counter);
	
	} else if (status != NC_NOERR){
	
		handle_error(status);

	//-----------------------------------------------
	// Attribute exists, extend it
	//-----------------------------------------------
	} else{
	
		nc_inq_attlen(ncid,NC_GLOBAL,"restart_times", &length);
		if (status != NC_NOERR) handle_error(status);

		//printf("l = %lu %d %d\n",length,prior_restart_times[0],prior_restart_times[1]);

		int *new_restart_times = ((int*) calloc(length+1,sizeof(int)));
	
		for(int i=0;i<length;i++){ new_restart_times[i] = prior_restart_times[i];}

	
		new_restart_times[length] = file_time_counter;

		add_attribute(ncid,"restart_times",new_restart_times,length+1);
	
		//free(new_restart_times);
	}
	
	status = nc_enddef(ncid);
	if (status != NC_NOERR) handle_error(status);
	
#endif
	//-----------------------------------------------
	// Initialize physical time variable
	//-----------------------------------------------
	int var_id;
	
	status = nc_inq_varid (ncid,"time", &var_id);
	if (status != NC_NOERR) handle_error(status);

	status = nc_get_var1_double(ncid, var_id,&file_time_counter, &mtime);
	if (status != NC_NOERR) handle_error(status);
	
	mtime *= 3600;	// convert from hours to seconds (input file should be in hours!)

	printf("-------------------------------------\n");

	if(isRestartRun){
		int time_steps_left = (int) ((number_of_time_steps * dt - mtime) / dt);
		
		bigcounter = number_of_time_steps - time_steps_left;
		
		printf("Restart run: %d time steps remaining out of %d total\n",time_steps_left,number_of_time_steps);
	}
	
	printf("Starting physical time at %f hours\n",mtime/3600.0);

	printf("-------------------------------------\n");
	
	//-----------------------------------------------
	// Close file
	//-----------------------------------------------	
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);
}

/********************************************************
* Adds command line arguments to the file name
*
* filename0	- output netcdf file name
* filename1	- input netcdf file name
*********************************************************/
void add_commandline_arguments(char *filename0,const char *filename1){
	
	const char * fileExtension = strrchr(filename1,'.');	// find the final period, this will begin the file extension

	strncpy(filename0,filename1,strlen(filename1)-strlen(fileExtension));

	sprintf(filename0,"%s_%1.0f_%1.0f_%1.0f%s",filename0,inputs.hor_shear*1000000,inputs.vert_shear0*-100000,inputs.vert_shear1*-100000,fileExtension);
}

/********************************************************
* 
* filename	- name of netcdf file
*********************************************************/
void set_outfilename(char *myfilename){

	strcpy(inputPerturbationFileName,inputs.perturbation_file);
	strcpy(inputBasicStateFileName,inputs.basic_state_file);
	//---------------------------------------------------
	// For command line arguments, create a file name
	// that includes these arguments
	//---------------------------------------------------
	if(inputs.has_shear==true){

		add_commandline_arguments(myfilename,inputs.output_file);

	} else {
		strcpy(myfilename,inputs.output_file);
		//strcpy(myfilename,out_filename);
	}
	if(rank==0){
		printf("------------------------------------------------\n");
		printf("Output file = %s\n",myfilename);
		printf("------------------------------------------------\n");
	}	
	//---------------------------------------------------
	// If a restart run was requested but the restart file
	// does not exist...
	//---------------------------------------------------
	if(rank == 0 && isRestartRun && does_file_exist(myfilename) == false){
		
		printf("Restart file %s does not exist. Creating it and running from initial time.\n",myfilename);
		isRestartRun = false;
	}
	if(PARALLEL){
		broadcast_data_int(1,&isRestartRun);
	}
}

/********************************************************
*  att - name of attribute
* filename	- name of netcdf file
*********************************************************/
bool has_attribute(const char *myfilename,const char *att){
	
	int status,ncid;
	bool hasAtt = false;
	//-----------------------------------------------
	// Open file
	//-----------------------------------------------
	status = nc_open(myfilename, 0, &ncid);
	handle_file_error(status,myfilename);
	
	int type;
	size_t length;
	//-----------------------------------------------
	// Check for attribute
	//-----------------------------------------------	
	status = nc_inq_att(ncid, NC_GLOBAL,att, &type,&length);
	
	if (status == NC_NOERR){ hasAtt = true;}
	//-----------------------------------------------
	// Close file
	//-----------------------------------------------	
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);
	
	return hasAtt;
}

/********************************************************
*  att - name of attribute
* filename	- name of netcdf file
*********************************************************/
void set_attribute(const char *myfilename,const char *att,int value){
	
	int status,ncid;
	//-----------------------------------------------
	// Open file
	//-----------------------------------------------
	status = nc_open(myfilename, NC_WRITE, &ncid);
	handle_file_error(status,myfilename);
	//-----------------------------------------------
	// Update value of attribute
	//-----------------------------------------------	
	status = nc_put_att_int(ncid, NC_GLOBAL,att,NC_INT, 1,&value);
	if (status != NC_NOERR) handle_error(status);
	//-----------------------------------------------
	// Close file
	//-----------------------------------------------	
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);
}

/********************************************************
* 
*********************************************************/
void file_output_status(int a){

	const char *att = "last_file_successfully_written";
	//-----------------------------------------------
	// If attribute exists, set it
	//-----------------------------------------------
	if(has_attribute(filename,att)){
		
		set_attribute(filename,att,a);
	//-----------------------------------------------
	// Otherwise, create it
	//-----------------------------------------------
	} else {
		int ncid,status;
		
		status = nc_open(filename, NC_WRITE, &ncid);
		handle_file_error(status,filename);
		
		status = nc_redef(ncid);
		if (status != NC_NOERR) handle_error(status);
		
		add_attribute(ncid,att,a);
		
		status = nc_close(ncid);
		if (status != NC_NOERR) handle_error(status);
	}
}

/********************************************************
* 
*********************************************************/
int get_file_output_status(){
	
	int ncid,status;
	int file_status = -1;
	
	if(has_attribute(filename,"last_file_successfully_written")){
	
		status = nc_open(filename, 0, &ncid);
		handle_file_error(status,filename);
	
		status = nc_get_att_int(ncid, NC_GLOBAL,"last_file_successfully_written", &file_status);
		if (status != NC_NOERR) handle_error(status);
		
		status = nc_close(ncid);
		if (status != NC_NOERR) handle_error(status);	
	}
	
	return file_status;
}

/********************************************************
* Initialize an output file with all time independent
* data
* 
* filename	- name of netcdf file
*********************************************************/
void outfile_init(char *myfilename){

	//---------------------------------------------------
	// Create a new output file
	//---------------------------------------------------
	if(CREATE_NEW_OUTPUT_FILE && !isRestartRun){

		if(PARALLEL_IO){
		
			create_outfile(myfilename,true,true,true);
			temp_storage = ((double*) calloc(myNX*myNY*myNZ,sizeof(double)));
		
		} else if(!PARALLEL || rank==0){
		
			create_outfile(myfilename,true,true,true);
		}
	
		if(!PARALLEL || rank==0){	

			write_all_mbvars(myfilename);
			write_all_bvars(myfilename);
			write_all_coordinates(myfilename);
			write_topo_to_file(myfilename);
			free(itopo);
		}

	//---------------------------------------------------
	// Write output to an existing file
	//---------------------------------------------------		
	} else {
		
		if(!PARALLEL || rank==0){ 
			handle_output_to_existing_file(myfilename);
		}
		if(PARALLEL){
			if(isRestartRun){ broadcast_data_int(1,&bigcounter);}
		}
	}
}

/****************************************************************************
*
* 					FILE READING PROCEDURES
*
*****************************************************************************/

void get_grid_spacing(const char *myfilename,const char *name_dx,const char *name_dy,const char *name_dz,float grid_spacing[3]){
	
	int ncid,var_id,status;
	
	status = nc_open(myfilename, 0, &ncid);

	handle_file_error(status,myfilename);

	if (status != NC_NOERR) handle_error(status);
	
	status = nc_get_att_float(ncid, NC_GLOBAL, "dx", &grid_spacing[0]);
	
	if(status==NC_ENOTATT){ printf("Attribute %s not found!",name_dx);}	
	else if (status != NC_NOERR) handle_error(status);
	
	status = nc_get_att_float(ncid, NC_GLOBAL, "dy", &grid_spacing[1]);
	
	if(status==NC_ENOTATT){ printf("Attribute not found!");}	
	else if (status != NC_NOERR) handle_error(status);
	
	status = nc_get_att_float(ncid, NC_GLOBAL, "dz", &grid_spacing[2]);
	
	if(status==NC_ENOTATT){ printf("Attribute not found!");}	
	else if (status != NC_NOERR) handle_error(status);


	
}

/********************************************************
* Get data from a model output file
* 
* myfilename - name of file
*
*********************************************************/
void get_model_data(const char *myfilename,size_t time){

	//get_dims(myfilename,const char *x,const char *y,const char *z,size_t dims[3])

	size_t size = NX*NY*NZ;

	get_data(myfilename,"ubar",size,&IUBAR(0,0,0));
	get_data(myfilename,"vbar",size,&IVBAR(0,0,0));
	get_data(myfilename,"wbar",size,&IWBAR(0,0,0));
	get_data(myfilename,"thbar",size,&ITHBAR(0,0,0));
	get_data(myfilename,"qbar",size,&IQBAR(0,0,0));

	get_data(myfilename,"tb",NZ,&tb[0]);
	get_data(myfilename,"pib",NZ,&pib[0]);
	get_data(myfilename,"qb",NZ,&qb[0]);

	get_data_at_time(myfilename,"u-wind",time,us);
	get_data_at_time(myfilename,"v-wind",time,vs);
	get_data_at_time(myfilename,"w-wind",time,ws);
	get_data_at_time(myfilename,"theta",time,ths);
	get_data_at_time(myfilename,   "pi",time,pis);

	if(USE_MICROPHYSICS){
		
		get_data_at_time(myfilename,"qv",time,qvs);
		get_data_at_time(myfilename,"qc",time,qcs);
		get_data_at_time(myfilename,"qr",time,qrs);
	}

	if(USE_TURBULENT_STRESS){	get_data_at_time(myfilename,"fric",time,&IFRICTION(0,0,0));}

	double * mytopo = get_data2(myfilename,"topo",NX*NY);

	for(int i=0;i<NX-1;i++){
	for(int j=0;j<NY-1;j++){
	for(int k=0;k<NZ-1;k++){

		ITOPO(i,j,k) = mytopo[j+i*NY];

	}}}

	free(mytopo);

}


/********************************************************
* Get data from netcdf file at specified time
*
* filename 		- name of file
* varname		- name of variable
* time			- time at which to fetch data
* array			- pointer to return data
*********************************************************/
void get_data_at_time(const char *myfilename,const char *varname, size_t time, double *array,int nx,int ny,int nz){

	int status, ncid, ndims, nvars, ngatts, unlimdimid;
	int var_id;

	size_t start3d[] = {time, 0, 0, 0};
	size_t count3d[] = {1, (size_t)nx, (size_t)ny, (size_t)nz};

	size_t start2d[] = {time, 0, 0};
	size_t count2d[] = {1, (size_t)nx, (size_t)ny};


	/********************************************************
	* Open file, get variable ID, get data
	*********************************************************/
	status = nc_open(myfilename, 0, &ncid);
	
	handle_file_error(status,myfilename);

	if (status != NC_NOERR) handle_error(status);

	status = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
	if (status != NC_NOERR) handle_error(status);
	
	status = nc_inq_varid (ncid, varname, &var_id);
	
	if(status == NC_ENOTVAR){
		printf("Error: variable '%s' not found in file %s\n",varname,myfilename);
		exit(0);
	}
	
	if (status != NC_NOERR) handle_error(status);

	if(nz==0){
		status = nc_get_vara_double(ncid,var_id,start2d,count2d,array);
	} else {
		status = nc_get_vara_double(ncid,var_id,start3d,count3d,array);
	}
	
	if (status != NC_NOERR) handle_error(status);

	/*********************
	* Close netCDF dataset 
	**********************/
	status = nc_close(ncid);      
	if (status != NC_NOERR) handle_error(status);

}

/********************************************************
* Does a given netcdf file have a given variable name?
*
* myfilename - netcdf file
* varname - name of variable
* return true if yes, false if no
*********************************************************/
bool fileHasVar(const char *myfilename,const char *varname){

	int status, ncid, var_id;
	bool hasVar = true;

	/********************************************************
	* Open file, get variable ID, get data
	*********************************************************/
	status = nc_open(myfilename, 0, &ncid);

	handle_file_error(status,myfilename);

	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid (ncid, varname, &var_id);
	
	if(status==NC_ENOTVAR){ hasVar = false; }
	
	/*********************
	* Close netCDF dataset 
	**********************/
	status = nc_close(ncid);      
	if (status != NC_NOERR) handle_error(status);
	
	return hasVar;
}

/********************************************************
* Get data from netcdf file at specified time
*
* filename 		- name of file
* varname		- name of variable
* time			- time at which to fetch data
* array			- pointer to return data
*********************************************************/
void get_data_at_time(const char *myfilename,const char *varname, size_t time, double *array){

	get_data_at_time(myfilename,varname,time,array,NX,NY,NZ);

}

/********************************************************
* Get data from an NCEP reanalysis netcdf file
*
* filename 		- name of file
* varname		- name of variable
* size			- number of data elements
* array			- pointer to return data
*********************************************************/
void get_data(const char *myfilename,const char *varname, int size, double *array){

	int status, ncid, ndims, nvars, ngatts, unlimdimid;
	int var_id;

	/********************************************************
	* Open file, get variable ID, get data
	*********************************************************/
	status = nc_open(myfilename, 0, &ncid);

	handle_file_error(status,myfilename);

	if (status != NC_NOERR) handle_error(status);

	status = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid (ncid, varname, &var_id);
	
	if(status == NC_ENOTVAR){
		printf("Error: variable '%s' not found in file %s\n",varname,myfilename);
		exit(0);
	}
	
	if (status != NC_NOERR) handle_error(status);

	status = nc_get_var_double(ncid,var_id,array);
	if (status != NC_NOERR) handle_error(status);
	
	/********************************************************
	* Unpack data
	* unpacked value = add_offset + ( (packed value) * scale_factor )
	*********************************************************/
	double scale_factor;
	double add_offset;
	
	// see if there is a scale factor
	status = nc_get_att_double(ncid, var_id, "scale_factor", &scale_factor);

	if(status==NC_ENOTATT){ scale_factor = 1;}
	else if (status != NC_NOERR) handle_error(status);

	// see if there is an offset
	status = nc_get_att_double(ncid, var_id, "add_offset", &add_offset);

	if(status==NC_ENOTATT){ add_offset = 0;}	
	else if (status != NC_NOERR) handle_error(status);

	// if there is a scale factor and/or offset, add them
	if(scale_factor!=1 || add_offset!=0){

		for(int i=0;i<size;i++){ array[i] = add_offset + (array[i]*scale_factor);}
	}

	//printf("\nvariable id is %d\n",var_id);
	//printf("\nfile has %d dimensions and %d variables\n",ndims,nvars);
	//printf("\ndimensions of the array are %d %d %d\n",xdim,ydim,zdim);

	/*********************
	* Close netCDF dataset 
	**********************/
	status = nc_close(ncid);      
	if (status != NC_NOERR) handle_error(status);

}

/********************************************************
* Get data from an NCEP reanalysis netcdf file
*
* @param filename - name of file
* @param varname - name of variable
* @param size	- number of data elements
* @return array - a pointer to the beginning of the data array
*********************************************************/
double* get_data2(const char *myfilename,const char *varname, int size){

	double * array = (double *)malloc(size*sizeof(double));

	get_data(myfilename,varname,size,array);

	return array;
}


/********************************************************
*
* 
*
*********************************************************/
void get_dims(const char *myfilename,const char *x,const char *y,const char *z,size_t dims[3]) {

	int status, ncid, ndims, nvars, ngatts, unlimdimid,latid,lonid,levid;

	size_t xdim, ydim, zdim;  /* dimensions*/

	//printf("%s %s %s %s\n",myfilename,x,y,z);
	
	/********************************************************
	* Open file
	*********************************************************/
	status = nc_open(myfilename, 0, &ncid);

	handle_file_error(status,myfilename);
	
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
	if (status != NC_NOERR) handle_error(status);

	/********************************************************
	* Get the dimension IDs and then inquire their sizes
	*********************************************************/

	// latitude
	status = nc_inq_dimid(ncid, y, &latid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_dimlen(ncid, latid, &ydim);
	if (status != NC_NOERR) handle_error(status);

	// longitude
	status = nc_inq_dimid(ncid, x, &lonid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_dimlen(ncid, lonid, &xdim);
	if (status != NC_NOERR) handle_error(status);

	// height
	status = nc_inq_dimid(ncid, z, &levid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_dimlen(ncid, levid, &zdim);
	if (status != NC_NOERR) handle_error(status);

	dims[0] = xdim;
	dims[1] = ydim;
	dims[2] = zdim;

	status = nc_close(ncid);       /* close netCDF dataset */
	if (status != NC_NOERR) handle_error(status);

}

/********************************************************
*
* 
*
*********************************************************/
void handle_file_error(int status,const char* myfilename) {
	
	if(status==NC_ENOTFOUND || status==2){ 
		printf("Error: File %s not found\n",myfilename);
		exit(0);
	}
}

/********************************************************
*
* 
*
*********************************************************/
void handle_error(int status) {

	if (status != NC_NOERR) {
   		fprintf(stderr, "%s\n", nc_strerror(status));
   		exit(-1);
   }
}




/********************************************************
*
* 
*
*********************************************************/
//int main(){

//	interp_era();

//	initialize2();

//	char filename[20] = "outfile.nc";

//	create_outfile(filename,TRUE,TRUE,TRUE);
//	
//	write_all_mbvars(filename);
//	write_all_bvars(filename);
//	write_all_coordinates(filename);

//	write_all_pvars(filename,0);
//	write_all_pvars(filename,1);

//}

