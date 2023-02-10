#include "stdafx.h"
#include "Heating.h"
#include "mpi.h"
#include "string.h"
#include "initializer.h"
#include "pressure.h"
#include "driver_serial.h"

struct task {

	double ind_var_1;
	int ind_var_2;
};

/*********************************************************************
*
* @param numtasks number of processes
* @param totaltasks number of tasks to perform
**********************************************************************/
int* assign_number_of_tasks(int numtasks,int totaltasks){

	int * arr = (int *)malloc(numtasks*sizeof(int));

	int min_number_rows;
	int extra_rows;
	int number_rows;

	min_number_rows = totaltasks / numtasks;
	extra_rows = totaltasks % numtasks;

	for(int i=0;i<numtasks;i++){

		number_rows = (i+1 <= extra_rows) ? min_number_rows+1 : min_number_rows;

		arr[i] = number_rows;
	}
	
	return arr;
}

/*********************************************************************
*
*
**********************************************************************/
void run_ensemble(int argc, char *argv[]){

	int  numtasks, rank, len, rc; 
	char hostname[MPI_MAX_PROCESSOR_NAME];

	/*********************************************
	* Set up parallel environment
	**********************************************/

	// initialize MPI  
	MPI_Init(&argc,&argv);

	// get number of tasks 
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

	// get my rank  
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	MPI_Get_processor_name(hostname, &len);
	printf ("Number of tasks= %d My rank= %d Running on %s\n", numtasks,rank,hostname);

	/*********************************************
	* Create and open the file
	**********************************************/
	char filename[20];

	sprintf(filename, "%s%d%s", "myfile",rank,".txt");

	FILE *myfile = fopen(filename, "a");

	/*********************************************
	* Initialize model data
	**********************************************/
	initialize_basic_state();

	initialize_pressure_solver();

	heat.initialize(21.18,86.3,19.37,93.0,100000.,6.0);

	int xheat = heat.i_westlon-heat.i_eastlon;
	int yheat = heat.i_westlat-heat.i_eastlat;

	int stride = 4;

	int xl = 38, xh = 110;
	int yl = 16, yh = 48;

	int NXPROB = (xh-xl) / stride + 1;
	int NYPROB = (yh-yl) / stride + 1;

	int totaltasks = NXPROB * NYPROB;

	int * arr = assign_number_of_tasks(numtasks,totaltasks);
	
	int * icoords = (int *)malloc(arr[rank]*sizeof(double));
	int * jcoords = (int *)malloc(arr[rank]*sizeof(double));

	int sum = 0;

	for(int i=0;i<rank;i++){ sum = sum + arr[i];}

	int istart, iend;
	int jstart, jend;
	int count = 0;
	int counter = 0;


	for(int j=yl;j<=yh;j+=stride){
	for(int i=xl;i<=xh;i+=stride){

		if(count >= sum && count<=sum+arr[rank]-1){
	
			icoords[counter] = i;
			jcoords[counter] = j;
		
			counter++;
		}

		count++;

	}}

	for(int i=0;i<arr[rank];i++){

		printf("coords = %d %d \n",icoords[i],jcoords[i]);

		heat.setLocationGrid(jcoords[i],icoords[i],jcoords[i]+yheat,icoords[i]+xheat);

		heat.printInfo(myfile);
		
		run_serial_model(2016,myfile);	

		reinitialize();
	}

	fclose(myfile);

	MPI_Finalize();

}

/*********************************************************************
*
*
**********************************************************************/
void run_ensemble_perturbation(int argc, char *argv[]){

	int  numtasks, rank, len, rc; 
	char hostname[MPI_MAX_PROCESSOR_NAME];

	/*********************************************
	* Set up parallel environment
	**********************************************/

	// initialize MPI  
	MPI_Init(&argc,&argv);

	// get number of tasks 
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

	// get my rank  
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	MPI_Get_processor_name(hostname, &len);
	//printf ("Number of tasks= %d My rank= %d Running on %s\n", numtasks,rank,hostname);

	//---------------------------------------------------
	// Create and open the file to record data
	//---------------------------------------------------
	char filename[20];

	sprintf(filename, "%s%d%s", "myfile",rank,".txt");

	FILE *myfile = fopen(filename, "a");

	//---------------------------------------------------
	// Define independent variables
	//---------------------------------------------------
	const int nxvars = 5;
	const int nyvars = 7;
	
	const double xvars[nxvars] = {0.25,0.5,0.75,1.0,1.25};
	const int yvars[nyvars] = {340,344,348,352,356,360,364};

	int totalTasks = nxvars*nyvars;

	//---------------------------------------------------
	// Determine number of tasks per process
	//---------------------------------------------------
	int * taskCounts = assign_number_of_tasks(numtasks,totalTasks);

	struct task *allTasks = (task*)calloc(totalTasks,sizeof(task));

	for(int i=0;i<nxvars;i++){
		for(int j=0;j<nyvars;j++){
		
			allTasks[j+nyvars*i].ind_var_1 = xvars[i];
			allTasks[j+nyvars*i].ind_var_2 = yvars[j];
		}
	}

	int myStartIndex = 0;
	
	for(int i=0;i<rank;i++){ myStartIndex += taskCounts[i];}

	printf("Process %d starts at index %d with var1 = %f and var2 = %d\n",rank,myStartIndex,allTasks[myStartIndex].ind_var_1,allTasks[myStartIndex].ind_var_2);

	//---------------------------------------------------
	// Initialize model data
	//---------------------------------------------------
	initialize_basic_state();

	initialize_pressure_solver();

	init_stats();

	//---------------------------------------------------
	// Process loops through its tasks
	//---------------------------------------------------
	for(int i=0;i<taskCounts[rank];i++){

		double scale = allTasks[myStartIndex+i].ind_var_1;		
		int time = allTasks[myStartIndex+i].ind_var_2;
		
		printf("-----------------------------------------------\n");
		printf("Process %d doing var1 = %d and var2 = %f\n",rank,time,scale);
		printf("-----------------------------------------------\n");

		fprintf(myfile,"-----------------------------------------------\n");
		fprintf(myfile,"Process %d doing var1 = %d and var2 = %f\n",rank,time,scale);
		fprintf(myfile,"-----------------------------------------------\n");
				
		reinitialize_perturbation(time,scale);
		
		run_serial_model(5,myfile);
	}



	fclose(myfile);

	MPI_Finalize();


}