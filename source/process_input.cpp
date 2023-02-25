#include "stdafx.h"
#include "process_input.h"

//---------------------------------------------------------
// Input file delimiters
//---------------------------------------------------------
const char *grid_section = "&grid";
const char *file_section = "&files";
const char *settings_section = "&settings";
const char *physics_section = "&physics";
const char *boundaries_section = "&boundaries";
const char *initialization_section = "&initialization_constants";
const char *delimiter = "/";

struct input_params inputs;


//---------------------------------------------------------
// Input file labels
//---------------------------------------------------------
const char gridnames[20][20] = {
	"NX","NY","NZ",
	"dx","dy","dz","dt",
	"corner_lat","corner_lon",
	"time_steps","output_frequency",
	"rayleigh_damping_z"};

const char infilenames[20][300] = {
	"basic_state_input_file",
	"perturbation_input_file",
	"output_file",
	"heat_budget_filename",
	"moisture_budget_filename",
	"pe_budget_filename",
	"pv_budget_filename",
	"vorticity_budget_filename"
};

/*********************************************************************
*
*
**********************************************************************/
void initialize_input_struct(struct input_params inputs){
	
	inputs.has_shear = false;
}

/*********************************************************************
*
*
**********************************************************************/
void print_input_params(struct input_params inputs){
	
	printf("\n");
	printf("nx = %d ny = %d nz = %d\n",inputs.nx,inputs.ny,inputs.nz);
	printf("dx = %f dy = %f dz = %f\n",inputs.dx,inputs.dy,inputs.dz);
	printf("corner_lat = %f corner_lon = %f\n",inputs.corner_lat,inputs.corner_lon);
	printf("\n");
	printf("basic state file = %s\n",inputs.basic_state_file);
	printf("perturbation file = %s\n",inputs.perturbation_file);
	printf("\n");
	
}

/*********************************************************************
*
*
**********************************************************************/
void test_line_format(const char* buff){
	
	const char *colon_match = strchr(buff,':');
	
	if(colon_match==NULL){ 
		printf("Error: format error in parameters file. Missing colon for: %s",buff);
		exit(0);
	}
	
	int loc = colon_match - buff;
	
	if( (buff[loc-1] != ' ' && buff[loc-1] != '\t') || (buff[loc+1] != ' ' && buff[loc+1] != '\t') ){
		
		printf("Error: format error in parameters file. Color needs spacing at line: %s",buff);
		exit(0);
	}
	
}

/*********************************************************************
*
*
**********************************************************************/
void attach_file_base(char* file,const char* filebase){

	char tempstring[300];

	strcpy(tempstring,file);
	strcpy(file,filebase);
	strcat(file,tempstring);
}

/*********************************************************************
*
*
**********************************************************************/
bool set_field(const char* buff,const char* match,int *value){
	
	if(strncmp(buff,match,strlen(match))==0){

		test_line_format(buff);

		sscanf(buff,"%*s %*s %d",value);
		return true;
	}
	return false;
}

/*********************************************************************
*
*
**********************************************************************/
bool set_field(const char* buff,const char* match,double *value){

	if(strncmp(buff,match,strlen(match))==0){

		test_line_format(buff);

		sscanf(buff,"%*s %*s %lf",value);
		return true;
	}
	return false;
}

/*********************************************************************
*
*
**********************************************************************/
bool set_field(const char* buff,const char* match,char* value){
	
	if(strncmp(buff,match,strlen(match))==0){
		
		test_line_format(buff);
		
		sscanf(buff,"%*s %*s %s",value);
		return true;
	}
	return false;
}

/*********************************************************************
*
*
**********************************************************************/
struct input_params read_input_file(const char* infile){
	
	const int bufflength = 255;	
	char buff[bufflength];
			
	FILE * fp = fopen(infile, "r");
	
	if(fp==NULL){ printf("Error: file %s not found\n", infile);}
	
	while( fgets(buff,bufflength,fp) != NULL){
		//---------------------------------------------------------
		// Read grid parameters
		//---------------------------------------------------------
		if(strncmp(buff,grid_section,strlen(grid_section))==0){
			
			while( fgets(buff,bufflength,fp) != NULL){
				
				if(strncmp(buff,delimiter,1)==0){ break; }
				
				set_field(buff,gridnames[0],&inputs.nx);
				set_field(buff,gridnames[1],&inputs.ny);
				set_field(buff,gridnames[2],&inputs.nz);
				set_field(buff,gridnames[3],&inputs.dx);
				set_field(buff,gridnames[3],&inputs.dy);
				set_field(buff,gridnames[5],&inputs.dz);
				set_field(buff,gridnames[6],&inputs.dt);
				set_field(buff,gridnames[7],&inputs.corner_lat);
				set_field(buff,gridnames[8],&inputs.corner_lon);
				set_field(buff,gridnames[9],&inputs.time_steps);
				set_field(buff,gridnames[10],&inputs.output_frequency);
				set_field(buff,gridnames[11],&inputs.rayleigh_damping_z);
				set_field(buff,"height_lowest_level",&inputs.height_lowest_level);
				set_field(buff,"shift_prime_meridian",&inputs.shift_prime_meridian);
											
			}	
		}	
		//---------------------------------------------------------
		// Read settings section
		//---------------------------------------------------------
		if(strncmp(buff,settings_section,strlen(settings_section))==0){
			
			while( fgets(buff,bufflength,fp) != NULL){
				
				if(strncmp(buff,delimiter,1)==0){ break; }
				
				set_field(buff,"is_restart_run",&inputs.is_restart_run);
								
				set_field(buff,"run_heat_budget",&inputs.run_heat_budget);
				set_field(buff,"run_moisture_budget",&inputs.run_moisture_budget);
				set_field(buff,"run_pe_budget",&inputs.run_pe_budget);
				set_field(buff,"run_pv_budget",&inputs.run_pv_budget);
				set_field(buff,"run_vorticity_budget",&inputs.run_vorticity_budget);
				set_field(buff,"basic_state_init_option",&inputs.basic_state_init_option);
				set_field(buff,"perturbation_init_option",&inputs.perturbation_init_option);
				set_field(buff,"perturbationFileTime",&inputs.perturbationFileTime);
				set_field(buff,"startOutfileAt",&inputs.startOutfileAt);
				set_field(buff,"create_new_output_file",&inputs.create_new_output_file);
				set_field(buff,"verbose",&inputs.verbose);
				set_field(buff,"print_courant_number",&inputs.print_courant_number);
				
			}	
		}
		//---------------------------------------------------------
		// Read physics section
		//---------------------------------------------------------
		if(strncmp(buff,physics_section,strlen(physics_section))==0){
			
			while( fgets(buff,bufflength,fp) != NULL){
				
				if(strncmp(buff,delimiter,1)==0){ break; }

				set_field(buff,"microphysics_option",&inputs.microphysics_option);
				set_field(buff,"use_surface_heat_flux",&inputs.use_surface_heat_flux);
				set_field(buff,"turbulence_option",&inputs.turbulence_option);
				set_field(buff,"water_temp",&inputs.water_temp);
				set_field(buff,"rain_fallout",&inputs.rain_fallout);
				set_field(buff,"use_explicit_diffusion",&inputs.use_explicit_diffusion);
				set_field(buff,"diffusion_order",&inputs.diffusion_order);
				set_field(buff,"kdiffh",&inputs.kdiffh);
				set_field(buff,"kdiffv",&inputs.kdiffv);
			}	
		}

		//---------------------------------------------------------
		// Read boundaries section
		//---------------------------------------------------------
		if(strncmp(buff,boundaries_section,strlen(boundaries_section))==0){
			
			while( fgets(buff,bufflength,fp) != NULL){
				
				if(strncmp(buff,delimiter,1)==0){ break; }

				set_field(buff,"periodic_ew_boundaries",&inputs.periodic_ew_boundaries);
				
				set_field(buff,"north",&inputs.boundary_width_north);
				set_field(buff,"south",&inputs.boundary_width_south);
				set_field(buff,"east",&inputs.boundary_width_east);
				set_field(buff,"west",&inputs.boundary_width_west);
			}	
		}
		
		//---------------------------------------------------------
		// Read file names
		//---------------------------------------------------------		
		if(strncmp(buff,file_section,strlen(file_section))==0){
			
			char input_filebase[300];
			char output_filebase[300];
			bool hasInputFileBase  = false;
			bool hasOutputFileBase  = false;
			
			while( fgets(buff,bufflength,fp) != NULL){
				
				if(strncmp(buff,delimiter,1)==0){ break; }
				
				set_field(buff,infilenames[0],inputs.basic_state_file);
				set_field(buff,infilenames[1],inputs.perturbation_file);
				set_field(buff,infilenames[2],inputs.output_file);

				set_field(buff,infilenames[3],inputs.heat_budget_filename);
				set_field(buff,infilenames[4],inputs.moisture_budget_filename);
				set_field(buff,infilenames[5],inputs.pe_budget_filename);
				set_field(buff,infilenames[6],inputs.pv_budget_filename);
				set_field(buff,infilenames[7],inputs.vorticity_budget_filename);
				
				if(!hasInputFileBase){
					hasInputFileBase = set_field(buff,"input_file_base",input_filebase);
				}
				
				if(!hasOutputFileBase){
					hasOutputFileBase = set_field(buff,"output_file_base",output_filebase);
				}
			}
			//---------------------------------------------------------
			// Attach full file paths to file names
			//---------------------------------------------------------				
			if(hasInputFileBase){
				
				attach_file_base(inputs.basic_state_file,input_filebase);
				attach_file_base(inputs.perturbation_file,input_filebase);
			}
			
			if(hasOutputFileBase){
				
				attach_file_base(inputs.output_file,output_filebase);
				
				attach_file_base(inputs.heat_budget_filename,output_filebase);
				attach_file_base(inputs.moisture_budget_filename,output_filebase);
				attach_file_base(inputs.pe_budget_filename,output_filebase);
				attach_file_base(inputs.pv_budget_filename,output_filebase);
				attach_file_base(inputs.vorticity_budget_filename,output_filebase);
			}
		}
	
		//---------------------------------------------------------
		// Initialization constants section
		//---------------------------------------------------------
		if(strncmp(buff,initialization_section,strlen(initialization_section))==0){
		
			while( fgets(buff,bufflength,fp) != NULL){
			
				if(strncmp(buff,delimiter,1)==0){ break; }
			
				set_field(buff,"vortex_initialize",&inputs.vortex_initialize);
				set_field(buff,"vortex_latitude",&inputs.vortex_latitude);
				set_field(buff,"vortex_longitude",&inputs.vortex_longitude);
				set_field(buff,"vortex_radius",&inputs.vortex_radius);
				set_field(buff,"vortex_upper_warm_anomaly",&inputs.vortex_upper_warm_anomaly);
				set_field(buff,"vortex_lower_cold_anomaly",&inputs.vortex_lower_cold_anomaly);
				set_field(buff,"vortex_height_wind_max",&inputs.vortex_height_wind_max);
				set_field(buff,"vortex_vertical_extent",&inputs.vortex_vertical_extent);
				set_field(buff,"vortex_rh_top",&inputs.vortex_rh_top);
				set_field(buff,"vortex_rh_bottom",&inputs.vortex_rh_bottom);
				set_field(buff,"vortex_rh_prime",&inputs.vortex_rh_prime);
				set_field(buff,"vortex_rh_radius",&inputs.vortex_rh_radius);
				set_field(buff,"vortex_rh_max",&inputs.vortex_rh_max);

			}	
		}
	}
	
	
		//hor_shear   = atof(argv[1])*1.0e-4;
		//vert_shear0 = atof(argv[2])*1.0e-3;
		//vert_shear1 = atof(argv[3])*1.0e-3;
	//printf("%s",buff);
	//printf("nx = %d ny = %d nz = %d\n",inputs.nx,inputs.ny,inputs.nz);
	
	fclose(fp);
	
	return inputs;
}

/*********************************************************************
*
* mpirun -np 4 ./solve.exe -f input_params.txt -s "1.0 -2.0 -1.0"
**********************************************************************/
int process_command_line_args(int argc, char *argv[]){
	
	//struct input_params inputs;
	bool input_file_given = false;
	
	
	initialize_input_struct(inputs);
	
	int c;
	
	char *token0;
	char *token1;
	char *token2;
	
	while ((c = getopt (argc, argv, "f:s:")) != -1){
		switch (c){
			case 'f':
			
				input_file_given = true;
				inputs = read_input_file(optarg);
				
				break;
			case 's':
				
				token0 = strtok(optarg," ");
				token1 = strtok(NULL," ");
				token2 = strtok(NULL," ");
				
				//printf("%s %s %s\n",token0,token1,token2);
				
				inputs.has_shear = true;
				inputs.hor_shear   = atof(token0)*1.0e-4;
				inputs.vert_shear0 = atof(token1)*1.0e-3;
				inputs.vert_shear1 = atof(token2)*1.0e-3;
				
				break;
			case '?':
				if (optopt == 'f')
				  fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else if (isprint (optopt))
				  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
				  fprintf (stderr,
				           "Unknown option character `\\x%x'.\n",
				           optopt);
				return 1;
			default:
				abort ();
		}
	  }
	
	  if(!input_file_given){
		  printf("Error: no input file provided!\n");
		  exit(0);
	  }
	
	  return 1;
}
#if 0
/*********************************************************************
*
*
**********************************************************************/
int main(int argc, char *argv[]){
	
	process_command_line_args(argc,argv);
	
}
#endif