#include "stdafx.h"
#include "energy.h"
#include "initializer.h"
#include "trajectory.h"
#include "budgets.h"
#include "ensemble.h"
#include "data_initializer.h"
#include "process_input.h"


/*********************************************************************
*
*
**********************************************************************/
int main(int argc, char *argv[]){

	process_command_line_args(argc,argv);

	initialize_globals();

	//argcount = argc;

	//if(argc==4){

		//printf("argc = %d %s %s %s\n",argc,argv[1],argv[2],argv[3]);

		//hor_shear   = atof(argv[1])*1.0e-4;
		//vert_shear0 = atof(argv[2])*1.0e-3;
		//vert_shear1 = atof(argv[3])*1.0e-3;
		
		//printf("argc = %d %f %f %f\n",argc,hor_shear,vert_shear0,vert_shear1);
	//}

	//const int nxc = atoi(argv[1]);
	//printf("nxc = %d\n",nxc);
	//trajectory_from_file("../model_output/outfile0417_3.nc");
	//return 0;

	//initialize_basic_state_idealized();
	//exit(0);
	//----------------------------------------------------------------
	// Calcuate QG-omega
	//----------------------------------------------------------------	
	if(CALCULATE_OMEGA){
		
		if(PARALLEL){ printf("Must use serial version!\n"); exit(0);}
		
		const char *inputfile = OUTPUT_FILE_PATH"outfile0227_1.nc";
		
		const int tcount = 100;//100;
		
		int otimes[tcount];
		
		for(int i=0;i<tcount;i++){ otimes[i] = i;}
		
		//int otimes[] = {85,86,87,88,89};
		
		calculate_omega_driver(inputfile,"../model_output/qgomega0227_1.nc",&otimes[0],tcount);
		
		exit(0);
	}
		
	//----------------------------------------------------------------
	// Parallel version
	//----------------------------------------------------------------
	if(PARALLEL && !ENERGY && ENSEMBLE==0){

		run_parallel_model(argc,argv);

	//----------------------------------------------------------------
	// Serial version
	//----------------------------------------------------------------		
	} else if(!ENERGY && ENSEMBLE==0){
		
		if(PV_BUDGET || HEAT_BUDGET || VORTICITY_BUDGET){
			
			printf("Budgets not supported in serial mode!\n");
			
			return(0);
		}
		
		initialize_serial();

		run_model(number_of_time_steps);

	//----------------------------------------------------------------
	// Energy budget from file
	//----------------------------------------------------------------		
	} else if(ENSEMBLE==0){

		process_command_line_args_energy(argc,argv);
		
		initialize_basic_state_from_output_file(energyBudgetFileName);

		energy_budget_from_file(energyBudgetFileName);

	//----------------------------------------------------------------
	// Ensemble version
	//----------------------------------------------------------------		
	} else {
		
		if(ENSEMBLE==1){
			run_ensemble(argc,argv);
		}
		if(ENSEMBLE==2){
			run_ensemble_perturbation(argc,argv);
		}
	}
	
}
