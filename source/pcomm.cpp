#include "stdafx.h"
#include "mpi.h"

extern int *scounts,*rcounts,*sdisps,*rdisps;
extern MPI_Datatype *types;

extern MPI_Comm comm_cart,comm_row,comm_col;

MPI_Datatype mysubarray[4];
MPI_Datatype mysubarray_2d[4];

MPI_Datatype *mysubarrayptr,*mysubarrayptr_2d;
MPI_Datatype nsBound,ewBound,interior,interior_2d;
MPI_Datatype nsBoundOnePoint,ewBoundOnePoint;

MPI_Request request;
MPI_Request *requests;

MPI_Status status;
MPI_Status *statuses;

//-------------------------------------------------
// A collection of information to help with sending
// and receiving data using AlltoAll MPI commands
//-------------------------------------------------
struct comm_layout {

	int *scounts;	// number of items to send to each other process
	int *rcounts;	// number of items to be received from each other process
	int *sdisps;	// displacement into send buffer to begin data to be sent
	int *rdisps;	// displacement into receive buffer to receive data
	MPI_Datatype *mytype;	// data type of send process
	MPI_Datatype *neighbortype;	// data type of receiving processes

} row_alltoall, col_alltoall, vert_alltoall,scatter_all;

//-------------------------------------------------
// function prototypes
//-------------------------------------------------
void allocateStruct(struct comm_layout*,int);
inline void nsExchange(double *);
inline void ewExchange(double *);
inline void cornerExchange(double *);
void initComm(int,int,int,int,MPI_Datatype *);
void initComm2d(int,int,int,int,int,int,MPI_Datatype *);
#if 0
/*********************************************************************
* Distribute initial data to each process
**********************************************************************/
void distributeArrays(){

	/******************************************************
	* Root process has already loaded input data to its
	* arrays. Now divide them up and distribute them
	* to each process.
	*******************************************************/
	if(rank==0){

		/******************************************************
		* Give the root process its share of each array.
		*******************************************************/
		for(int i=0;i<myNX;i++){
		for(int j=0;j<myNY;j++){
		for(int k=0;k<myNZ;k++){

			m_ubar[I2(i+halo_buffer,j+halo_buffer,k)] = IUBAR(i,j,k);
			m_vbar[I2(i+halo_buffer,j+halo_buffer,k)] = IVBAR(i,j,k);
			m_wbar[I2(i+halo_buffer,j+halo_buffer,k)] = IWBAR(i,j,k);
			m_thbar[I2(i+halo_buffer,j+halo_buffer,k)] = ITHBAR(i,j,k);
			m_qbar[I2(i+halo_buffer,j+halo_buffer,k)] = IQBAR(i,j,k);
			m_pbar[I2(i+halo_buffer,j+halo_buffer,k)] = IPBAR(i,j,k);

			istopos[I2(i+halo_buffer,j+halo_buffer,k)] = IISTOPO(i,j,k);
			uistopos[I2(i+halo_buffer,j+halo_buffer,k)] = IUISTOPO(i,j,k);
			vistopos[I2(i+halo_buffer,j+halo_buffer,k)] = IVISTOPO(i,j,k);

			frictions[I2(i+halo_buffer,j+halo_buffer,k)] = IFRICTION(i,j,k);

		}}}

		/******************************************************
		* Distribute the rest of each array to other processes
		*******************************************************/
		for(int i=1;i<numtasks;i++){
		
			MPI_Send(&IUBAR (ibs[i],jbs[i],0),1,mysubarrayptr[i],i,123,MPI_COMM_WORLD);
			MPI_Send(&IVBAR (ibs[i],jbs[i],0),1,mysubarrayptr[i],i,124,MPI_COMM_WORLD);
			MPI_Send(&IWBAR (ibs[i],jbs[i],0),1,mysubarrayptr[i],i,122,MPI_COMM_WORLD);
			MPI_Send(&ITHBAR(ibs[i],jbs[i],0),1,mysubarrayptr[i],i,125,MPI_COMM_WORLD);
			MPI_Send(&IPBAR (ibs[i],jbs[i],0),1,mysubarrayptr[i],i,131,MPI_COMM_WORLD);
			MPI_Send(&IQBAR (ibs[i],jbs[i],0),1,mysubarrayptr[i],i,130,MPI_COMM_WORLD);

			MPI_Send(&IISTOPO (ibs[i],jbs[i],0),1,mysubarrayptr[i],i,126,MPI_COMM_WORLD);
			MPI_Send(&IUISTOPO(ibs[i],jbs[i],0),1,mysubarrayptr[i],i,127,MPI_COMM_WORLD);
			MPI_Send(&IVISTOPO(ibs[i],jbs[i],0),1,mysubarrayptr[i],i,128,MPI_COMM_WORLD);

			MPI_Send(&IFRICTION(ibs[i],jbs[i],0),1,mysubarrayptr[i],i,129,MPI_COMM_WORLD);
		}
	}

	/******************************************************
	* Other processes receive array.
	*******************************************************/
	else {

		MPI_Recv(m_ubar,1,interior,0,123,MPI_COMM_WORLD,&status);
		MPI_Recv(m_vbar,1,interior,0,124,MPI_COMM_WORLD,&status);
		MPI_Recv(m_wbar,1,interior,0,122,MPI_COMM_WORLD,&status);
		MPI_Recv(m_thbar,1,interior,0,125,MPI_COMM_WORLD,&status);
		MPI_Recv(m_pbar,1,interior,0,131,MPI_COMM_WORLD,&status);
		MPI_Recv(m_qbar,1,interior,0,130,MPI_COMM_WORLD,&status);

		MPI_Recv(istopos,1,interior,0,126,MPI_COMM_WORLD,&status);
		MPI_Recv(uistopos,1,interior,0,127,MPI_COMM_WORLD,&status);
		MPI_Recv(vistopos,1,interior,0,128,MPI_COMM_WORLD,&status);

		MPI_Recv(frictions,1,interior,0,129,MPI_COMM_WORLD,&status);
	}
}
#endif
#if 0
/*********************************************************************
* Distribute initial data to each process
**********************************************************************/
void distributeArray(double *svar){

	/******************************************************
	* Root process has already loaded input data to its
	* arrays. Now divide them up and distribute them
	* to each process.
	*******************************************************/
	if(rank==0){
		/******************************************************
		* Give the root process its share of each array.
		*******************************************************/
		for(int i=0;i<myNX;i++){
		for(int j=0;j<myNY;j++){
		for(int k=0;k<myNZ;k++){

			svar[I2(i+halo_buffer,j+halo_buffer,k)] = IUBAR(i,j,k);
		}}}
		/******************************************************
		* Distribute the rest of each array to other processes
		*******************************************************/
		for(int i=1;i<numtasks;i++){
			
			MPI_Send(&IUBAR(ibs[i],jbs[i],0),1,mysubarrayptr[i],i,123,MPI_COMM_WORLD);
		}
	}
	/******************************************************
	* Other processes receive array.
	*******************************************************/
	else {
		MPI_Recv(svar,1,interior,0,123,MPI_COMM_WORLD,&status);
	}
}
#endif
/*********************************************************************
* Distribute initial data to each process
*
* svar - smaller array, local to each process
* var - bigger array, on the root process
**********************************************************************/
void distributeArray_3d(double *svar,double *var){

	/******************************************************
	* Root process has already loaded input data to its
	* arrays. Now divide them up and distribute them
	* to each process.
	*******************************************************/
	if(rank==0){
		/******************************************************
		* Give the root process its share of each array.
		*******************************************************/
		for(int i=0;i<myNX;i++){
		for(int j=0;j<myNY;j++){
		for(int k=0;k<myNZ;k++){

			svar[I2(i+halo_buffer,j+halo_buffer,k)] = var[(i)*NY*NZ+(j)*NZ+(k)];
		}}}
		/******************************************************
		* Distribute the rest of each array to other processes
		*******************************************************/
		int offset;
		
		for(int i=1;i<numtasks;i++){
			
			offset = ibs[i]*NY*NZ+jbs[i]*NZ;
			
			MPI_Send(var+offset,1,mysubarrayptr[i],i,123,MPI_COMM_WORLD);
		}
	}
	/******************************************************
	* Other processes receive array.
	*******************************************************/
	else {
		MPI_Recv(svar,1,interior,0,123,MPI_COMM_WORLD,&status);
	}
}

/*********************************************************************
* Distribute initial data to each process
*
* svar - smaller array, local to each process
* var - bigger array, on the root process
**********************************************************************/
void distributeArray_2d(double *svar,double *var){

	/******************************************************
	* Root process has already loaded input data to its
	* arrays. Now divide them up and distribute them
	* to each process.
	*******************************************************/
	if(rank==0){
		/******************************************************
		* Give the root process its share of each array.
		*******************************************************/
		for(int i=0;i<myNX;i++){
		for(int j=0;j<myNY;j++){

			svar[INDEX2D(i+halo_buffer,j+halo_buffer)] = var[(i)*NY+(j)];
		}}
		/******************************************************
		* Distribute the rest of each array to other processes
		*******************************************************/
		int offset;
		
		for(int i=1;i<numtasks;i++){
			
			offset = ibs[i]*NY + jbs[i];
			
			MPI_Send(var+offset,1,mysubarrayptr_2d[i],i,123,MPI_COMM_WORLD);
		}
	}
	/******************************************************
	* Other processes receive array.
	*******************************************************/
	else {
		MPI_Recv(svar,1,interior_2d,0,123,MPI_COMM_WORLD,&status);
	}
}
#if 0
/*********************************************************************
* All processes send their data to all other processes
*
* @param m_var - the small local buffer from which to send data
* @param var - the large array which receives the data
**********************************************************************/
void alltoall(double * m_var,double var[NX][NY][NZ]){

	for(int i=0;i<numtasks;i++){

		scounts[i] = 1;
		sdisps[i] = 0;
		rcounts[i] = 1;
		rdisps[i] = (&var[ibs[i]][jbs[i]][0] - &var[0][0][0])*sizeof(double);
		types[i] = interior;
	}

	MPI_Alltoallw(m_var,scounts,sdisps,types,&var[0][0][0],rcounts,rdisps,mysubarrayptr,MPI_COMM_WORLD);
}
#endif

/*********************************************************************
* Processes in each row exchange data amoung themselves
**********************************************************************/
void pres_row_alltoall(double * row){

	MPI_Alltoallw(
		row,
		row_alltoall.scounts,
		row_alltoall.sdisps,
		row_alltoall.mytype,
		row,
		row_alltoall.rcounts,
		row_alltoall.rdisps,
		row_alltoall.neighbortype,
		comm_row
	);
}

/*********************************************************************
* Processes in each row exchange data amoung themselves
**********************************************************************/
void pres_row_alltoall2(double * row1,double * row2){

	MPI_Alltoallw(
		row1,
		row_alltoall.scounts,
		row_alltoall.sdisps,
		row_alltoall.mytype,
		row2,
		row_alltoall.rcounts,
		row_alltoall.rdisps,
		row_alltoall.neighbortype,
		comm_row
	);
}

/*********************************************************************
* Processes in each row exchange data amoung themselves
**********************************************************************/
void pres_row_alltoall2_reverse(double * row1,double * row2){

	MPI_Alltoallw(
		row2,
		row_alltoall.rcounts,
		row_alltoall.rdisps,
		row_alltoall.neighbortype,
		row1,
		row_alltoall.scounts,
		row_alltoall.sdisps,
		row_alltoall.mytype,
		comm_row
	);
}

/*********************************************************************
* Processes in each row exchange data amoung themselves
**********************************************************************/
void pres_vert_alltoall(double * row){

	MPI_Alltoallw(
		row,
		vert_alltoall.scounts,
		vert_alltoall.sdisps,
		vert_alltoall.mytype,
		row,
		vert_alltoall.rcounts,
		vert_alltoall.rdisps,
		vert_alltoall.neighbortype,
		comm_col
	);
}

/*********************************************************************
* Processes in each row exchange data amoung themselves
**********************************************************************/
void pres_vert_alltoall_reverse(double * row){

	MPI_Alltoallw(
		row,
		vert_alltoall.rcounts,
		vert_alltoall.rdisps,
		vert_alltoall.neighbortype,
		row,
		vert_alltoall.scounts,
		vert_alltoall.sdisps,
		vert_alltoall.mytype,
		comm_col
	);
}

/*********************************************************************
* Processes in each column exchange data amoung themselves
**********************************************************************/
void pres_col_alltoall(double * col){

	MPI_Alltoallw(
		col,
		col_alltoall.scounts,
		col_alltoall.sdisps,
		col_alltoall.mytype,
		col,
		col_alltoall.rcounts,
		col_alltoall.rdisps,
		col_alltoall.neighbortype,
		comm_col
	);
}

/*********************************************************************
* Processes in each column exchange data amoung themselves
**********************************************************************/
void pres_col_alltoall2(double * col1,double * col2){

	MPI_Alltoallw(
		col1,
		col_alltoall.scounts,
		col_alltoall.sdisps,
		col_alltoall.mytype,
		col2,
		col_alltoall.rcounts,
		col_alltoall.rdisps,
		col_alltoall.neighbortype,
		comm_col
	);
}

/*********************************************************************
* Processes in each column exchange data amoung themselves
**********************************************************************/
void pres_col_alltoall_reverse(double * col){

	MPI_Alltoallw(
		col,
		col_alltoall.rcounts,
		col_alltoall.rdisps,
		col_alltoall.neighbortype,
		col,
		col_alltoall.scounts,
		col_alltoall.sdisps,
		col_alltoall.mytype,
		comm_col
	);
}
#if 0
/*********************************************************************
* 
**********************************************************************/
void wait(int r){
	
	MPI_Waitall(r,&requests[0],MPI_STATUSES_IGNORE);
}

/*********************************************************************
* Distribute a full 3D array from the root process to all
* processes
**********************************************************************/
void scatterArrays(double * m_var,double var[NX][NY][NZ]){

	for(int i=0;i<numtasks;i++){
		sdisps[i] = (&var[ibs[i]][jbs[i]][0] - &var[0][0][0])*sizeof(double);
	}

	MPI_Alltoallw(
		&var[0][0][0],
		scatter_all.scounts,
		sdisps,
		mysubarrayptr,
		m_var,
		scatter_all.rcounts,
		scatter_all.rdisps,
		scatter_all.mytype,
		MPI_COMM_WORLD
	);
}
#endif
#if 0
/*********************************************************************
* 3.53 2.39 1.85 1.61 1.29 1.32 1.03
**********************************************************************/
void gatherArrays(){

	if(rank!=0){
	
		MPI_Isend(m_vbar,1,interior,0,12,comm_cart,&request);

		MPI_Wait(&request,&status);

	} else {

		for(int i=1;i<numtasks;i++){

			MPI_Irecv(&u[ibs[i]][jbs[i]][0],1,mysubarrayptr[i],i,12,comm_cart,&requests[i-1]);
		}

		MPI_Waitall(numtasks-1,requests,statuses);
	}
}
#endif
#if 0
/*********************************************************************
* Send the interior points of the subarrays on each process to the
* root process.
*
* @param *m_var - pointer to beginning of buffer
* @param var - full size array on root to receive data
**********************************************************************/
void gatherArrays2(double * m_var,double var[NX][NY][NZ]){

	for(int i=0;i<numtasks;i++){

		scounts[i] = 0;
		sdisps[i] = 0;
		rcounts[i] = 0;
		rdisps[i] = (&var[ibs[i]][jbs[i]][0] - &var[0][0][0])*sizeof(double);
		types[i] = interior;
	}

	if(rank==0){
		for(int i=0;i<numtasks;i++){ rcounts[i] = 1;}
	}

	scounts[0] = 1;

	MPI_Alltoallw(m_var,scounts,sdisps,types,&var[0][0][0],rcounts,rdisps,mysubarrayptr,MPI_COMM_WORLD);
}
#endif
/*********************************************************************
* Send the interior points of the subarrays on each process to the
* root process.
*
* @param *m_var - pointer to beginning of buffer
* @param var - full size array on root to receive data
**********************************************************************/
void gatherArrays_3d(double * m_var,double * var){

	for(int i=0;i<numtasks;i++){

		scounts[i] = 0;
		sdisps[i] = 0;
		rcounts[i] = 0;
		rdisps[i] = ((ibs[i])*NY*NZ + (jbs[i])*NZ + (0)) * sizeof(double);
		types[i] = interior;
	}

	if(rank==0){
		for(int i=0;i<numtasks;i++){ rcounts[i] = 1;}
	}

	scounts[0] = 1;

	MPI_Alltoallw(m_var,scounts,sdisps,types,var,rcounts,rdisps,mysubarrayptr,MPI_COMM_WORLD);
}

/*********************************************************************
* Send the interior points of the subarrays on each process to the
* root process.
*
* @param *m_var - pointer to beginning of buffer
* @param var - full size array on root to receive data
**********************************************************************/
void gatherArrays_2d(double * m_var,double * var){

	for(int i=0;i<numtasks;i++){

		scounts[i] = 0;
		sdisps[i] = 0;
		rcounts[i] = 0;
		rdisps[i] = (ibs[i]*NY + jbs[i]) * sizeof(double);
		types[i] = interior_2d;
	}

	if(rank==0){
		for(int i=0;i<numtasks;i++){ rcounts[i] = 1;}
	}

	scounts[0] = 1;

	MPI_Alltoallw(m_var,scounts,sdisps,types,var,rcounts,rdisps,mysubarrayptr_2d,MPI_COMM_WORLD);
}

/*********************************************************************
* Create a new MPI derived datatype for communicating a 3D array
* among processes. The z-dimension is fixed to the model's number of
* vertical grid levels. Useful for communicating the interior points
* back to the root process for output to a file or for communicating
* the boundary points between neighboring processes.
*
* @param xstart,ystart - the starting and ending points in the full array
*					     on the local process
* @param subnx,subny - the sizes of the subarray on the local process
*					   all processes
* @param *comm	- pointer to MPI derived datatype
**********************************************************************/
void initComm(int xstart,int ystart,int subnx,int subny,MPI_Datatype *comm){

	int starts[3] = {xstart,ystart,0};	// starting location within the full sized array
	int subsizes[3] = {subnx,subny,myNZ}; // the sizes of the subarray, which excludes boundary points
	int bigsizes[3]  = {fNX,fNY,fNZ};	// the sizes of the local array, which includes any boundary points

	MPI_Type_create_subarray(3, bigsizes, subsizes, starts,MPI_ORDER_C, MPI_DOUBLE, comm);
	MPI_Type_commit(comm);
}

/*********************************************************************
* Create a new MPI derived datatype for communicating a 2D array
* among processes.
*
* @param xstart,ystart - the starting and ending points in the full array
*					     on the local process
* @param subnx,subny - the sizes of the subarray on the local process
* @param bigNX,bigNY - the sizes of the full domain distributed across
*					   all processes in communicator
* @param *comm	- pointer to MPI derived datatype
**********************************************************************/
void initComm2d(int xstart,int ystart,int subnx,int subny,int bigNX,int bigNY,MPI_Datatype *comm){

	int starts[2] = {xstart,ystart};
	int subsizes[2] = {subnx,subny};
	int bigsizes[2]  = {bigNX,bigNY};

	MPI_Type_create_subarray(2, bigsizes, subsizes, starts,MPI_ORDER_C, MPI_DOUBLE, comm);
	MPI_Type_commit(comm);
}

/*********************************************************************
* Create a new MPI derived datatype for communicating a 3D array
* among processes.
*
* @param xstart,ystart,zstart - the starting and ending points in the full array
*					     on the local process
* @param subnx,subny,subnz - the sizes of the subarray on the local process
* @param bigNX,bigNY,bigNZ - the sizes of the full domain distributed across
*					   all processes in communicator
* @param *comm	- pointer to MPI derived datatype
**********************************************************************/
void initComm3d(int xstart,int ystart,int zstart,int subnx,int subny,int subnz,int bigNX,int bigNY,int bigNZ,MPI_Datatype *comm){

	int starts[3] = {xstart,ystart,zstart};
	int subsizes[3] = {subnx,subny,subnz};
	int bigsizes[3]  = {bigNX,bigNY,bigNZ};

	MPI_Type_create_subarray(3, bigsizes, subsizes, starts,MPI_ORDER_C, MPI_DOUBLE, comm);
	MPI_Type_commit(comm);
}

/*********************************************************************
* Create a new MPI derived datatype for communicating a 3D array
* among processes.
*
* @param xstart,ystart,zstart - the starting and ending points in the full array
*					     on the local process
* @param subnx,subny,subnz - the sizes of the subarray on the local process
* @param bigNX,bigNY,bigNZ - the sizes of the full domain distributed across
*					   all processes in communicator
* @param *comm	- pointer to MPI derived datatype
**********************************************************************/
void initComm3d_complex(int xstart,int ystart,int zstart,int subnx,int subny,int subnz,int bigNX,int bigNY,int bigNZ,MPI_Datatype *comm){

	int starts[3] = {xstart,ystart,zstart};
	int subsizes[3] = {subnx,subny,subnz};
	int bigsizes[3]  = {bigNX,bigNY,bigNZ};

	MPI_Datatype complex_number;

	MPI_Type_contiguous(2,MPI_DOUBLE,&complex_number);

	MPI_Type_create_subarray(3, bigsizes, subsizes, starts,MPI_ORDER_C, complex_number, comm);
	MPI_Type_commit(comm);
}

/*********************************************************************
* 
*
**********************************************************************/
void allocateStruct(struct comm_layout *in,int size){

	in->neighbortype = (MPI_Datatype *)malloc(size*sizeof(MPI_Datatype));
	in->mytype = (MPI_Datatype *)malloc(size*sizeof(MPI_Datatype));
	in->scounts = (int *)calloc(size,sizeof(int));
	in->rcounts = (int *)calloc(size,sizeof(int));
	in->sdisps = (int *)calloc(size,sizeof(int));
	in->rdisps = (int *)calloc(size,sizeof(int));
}

/*********************************************************************
* Set up pressure communicators for 2D pressure field
*
**********************************************************************/
void presComm2d(){

	/****************************************************************
	* Set up the structure containing variables used to communication
	* amongst the processes in a single row
	*****************************************************************/
	allocateStruct(&row_alltoall,row_size);

	// each process within a row gets a copy of the MPI derived datatype
	// for each of the other processes within its row
	for(int i=0;i<row_size;i++){

		initComm2d(0,0,s_nx[i],myNY,NX,myNY,&row_alltoall.neighbortype[i]);
		row_alltoall.scounts[i] = 1;
		row_alltoall.rcounts[i] = 1;
	}

	scounts[row_rank] = 0;
	rcounts[row_rank] = 0;

	for(int i=0;i<row_size;i++){ row_alltoall.mytype[i] = row_alltoall.neighbortype[row_rank];}

	row_alltoall.rdisps[0] = 0;

	for(int i=1;i<row_size;i++){ row_alltoall.rdisps[i] = row_alltoall.rdisps[i-1] + s_nx[i-1]*myNY*sizeof(double);}

	for(int i=0;i<row_size;i++){ row_alltoall.sdisps[i] = row_alltoall.rdisps[row_rank];}

	/****************************************************************
	* Set up the structure containing variables used to communication
	* amongst the processes in a single column
	*****************************************************************/
	allocateStruct(&col_alltoall,col_size);

	// each process within a column gets a copy of the MPI derived datatype
	// for each of the other processes within its column
	for(int j=0;j<col_size;j++){

		initComm2d(0,0,myNX,s_ny[j],myNX,NY,&col_alltoall.neighbortype[j]);
		col_alltoall.scounts[j] = 1;
		col_alltoall.rcounts[j] = 1;

	}

	scounts[col_rank] = 0;
	rcounts[col_rank] = 0;

	for(int i=0;i<col_size;i++){ col_alltoall.mytype[i] = col_alltoall.neighbortype[col_rank];}

	col_alltoall.rdisps[0] = 0;

	for(int i=1;i<col_size;i++){ col_alltoall.rdisps[i] = col_alltoall.rdisps[i-1] + s_ny[i-1]*sizeof(double);}

	for(int i=0;i<col_size;i++){ col_alltoall.sdisps[i] = col_alltoall.rdisps[col_rank];}

	/****************************************************************
	* 
	*
	*****************************************************************/
	allocateStruct(&scatter_all,numtasks);

	for(int i=0;i<numtasks;i++){

		scatter_all.scounts[i] = 0;
		scatter_all.rcounts[i] = 1;
		scatter_all.rdisps[i] = 0;
		scatter_all.mytype[i] = interior;
	}

	if(rank==0){
		for(int i=0;i<numtasks;i++){ scatter_all.scounts[i] = 1;}
	}

}

/*********************************************************************
* Set up pressure communicators for 3D pressure field
*
**********************************************************************/
void presComm3d(){
	
	/****************************************************************
	* Set up the structure containing variables used to communicate
	* amongst the processes in a single column
	*****************************************************************/
	allocateStruct(&col_alltoall,col_size);
	
	/****************************************************************
	* Each process within a column gets a copy of the MPI derived 
	* datatype for each of the other processes within its column
	*****************************************************************/
	for(int j=0;j<col_size;j++){
		
		// this is a description of the data that is being received by the jth process in a column
		initComm3d(
			0,jbs[j],kbs[col_rank],	// each received piece of data is at the same vertical level (kbs[col_rank])
			myNX,				// the x dimension is the same for each piece within the vertical column
			s_ny[j],			// the y dimension can vary
			s_nz[col_rank],		// the z dimension will be the same for each piece of data received
			myNX,NY,NZ,			// these are the dimensions of the full buffer array
			&col_alltoall.neighbortype[j]
		);
		
		// this is a description of the data that is being sent to the jth process in a column
		initComm3d(
			0,jbs[rank],kbs[j], 	// since the data is being partitioned in the vertical, each process gets a piece beginning at a different k
			myNX,			// the x dimension is the same for each piece within the vertical column
			s_ny[col_rank], // the y dimension is the same for each piece within the vertical column
			s_nz[j],		// the z dimensions can vary, the upper regions can be one grid point smaller than the lower ones
			myNX,NY,NZ,		// these are the dimensions of the full buffer array
			&col_alltoall.mytype[j]
		);
		
		col_alltoall.scounts[j] = 1;
		col_alltoall.rcounts[j] = 1;
	}
	
	// columns do not exchange with their neighbors their own portion of the full domain
	scounts[col_rank] = 0;
	rcounts[col_rank] = 0;

	// the displacement into the array at which to receive data
	for(int i=0;i<col_size;i++){ col_alltoall.rdisps[i] = 0;}

	// the displacement into the array from which to send data
	for(int i=0;i<col_size;i++){ col_alltoall.sdisps[i] = col_alltoall.rdisps[col_rank];}
	
	/****************************************************************
	* Set up the structure containing variables used to communicate
	* between rows and columns
	*****************************************************************/
	allocateStruct(&row_alltoall,row_size);
	
	// each process within a row gets a copy of the MPI derived datatype
	// for each of the other processes within its row
	for(int i=0;i<row_size;i++){

		// this is a description of the data that is being received by the ith process in a row
		// the receiver is an east-west row
		initComm3d(0,0,kbs[col_rank],
			s_nx[i],
			pNY,
			pNZ,
			NX,pNY,NZ,
			&row_alltoall.neighbortype[i]
		);

		// this is a description of the data that is being sent to the ith process in a row
		// the sender is a north-south column
		initComm3d(0,0,kbs[col_rank],
			myNX,
			s_ny_p[i],
			pNZ,
			myNX,NY,NZ,
			&row_alltoall.mytype[i]
		);
	
		row_alltoall.scounts[i] = 1;
		row_alltoall.rcounts[i] = 1;
	}

	row_alltoall.rdisps[0] = 0;

	// the displacement into the array at which to receive data
	for(int i=1;i<row_size;i++){ row_alltoall.rdisps[i] = row_alltoall.rdisps[i-1] + s_nx[i-1]*pNY*NZ*sizeof(double);}

	row_alltoall.sdisps[0] = 0;

	// the displacement into the array from which to send data
	for(int i=1;i<row_size;i++){ row_alltoall.sdisps[i] = row_alltoall.sdisps[i-1] + s_ny_p[i-1]*NZ*sizeof(double);}
	
	/****************************************************************
	* Set up the structure containing variables used to communicate
	* in the vertical
	*****************************************************************/
	allocateStruct(&vert_alltoall,col_size);
	
	// each process within a row gets a copy of the MPI derived datatype
	// for each of the other processes within its row
	for(int i=0;i<col_size;i++){

		// this is a description of the data that is being received by the ith process
		// the receiver an up-down column
		initComm3d(ibs_p[col_rank],0,kbs[i],
			pNX,
			pNY,
			s_nz[i],
			NX,pNY,NZ,
			&vert_alltoall.neighbortype[i]
		);

		// this is a description of the data that the process is sending to the ith process
		// the sender is an up-down column
		initComm3d(ibs_p[i],0,kbs[col_rank],
			s_nx_p[i],
			pNY,
			pNZ,
			NX,pNY,NZ,
			&vert_alltoall.mytype[i]
		);
	
		vert_alltoall.scounts[i] = 1;
		vert_alltoall.rcounts[i] = 1;
	}

	// columns do not exchange with their neighbors their own portion of the full domain
	scounts[col_rank] = 0;
	rcounts[col_rank] = 0;

	vert_alltoall.rdisps[0] = 0;

	// the displacement into the array at which to receive data
	for(int i=1;i<col_size;i++){ vert_alltoall.rdisps[i] = 0;}

	vert_alltoall.sdisps[0] = 0;

	// the displacement into the array from which to send data
	for(int i=1;i<col_size;i++){ vert_alltoall.sdisps[i] = 0;}//vert_alltoall.sdisps[i-1] + s_nx_p[i-1]*pNY*NZ*sizeof(double);}
}

/*********************************************************************
* Construct the MPI derived data types which will be used in communication
*
**********************************************************************/
void initAllComms(){

	requests = (MPI_Request *)malloc(200*sizeof(MPI_Request));
	statuses = (MPI_Status *)malloc(200*sizeof(MPI_Status));

	// Northern/southern boundary of subarrays
	initComm(0,0,myNX,halo_buffer,&nsBound);
	initComm(0,0,myNX,1,&nsBoundOnePoint);

	// Eastern/western boundary of subarrays
	initComm(0,0,halo_buffer,myNY,&ewBound);
	initComm(0,0,1,myNY,&ewBoundOnePoint);

	// Interior of subarrays
	initComm(halo_buffer,halo_buffer,myNX,myNY,&interior);
	initComm2d(halo_buffer,halo_buffer,myNX,myNY,fNX,fNY,&interior_2d);

	if(HYDROSTATIC){ presComm2d();}
	else { presComm3d();}

	//--------------------------------------------------------------------
	// Create the data types used to break up the full domain arrays into
	// pieces for each process. The subarrays can be of four different sizes,
	// with each dimesion having the possibility of being one unit less than 
	// that of the root process's subarray.
	//--------------------------------------------------------------------
	mysubarrayptr    = (MPI_Datatype *)malloc(numtasks*sizeof(MPI_Datatype));
	mysubarrayptr_2d = (MPI_Datatype *)malloc(numtasks*sizeof(MPI_Datatype));

	int bigNX = s_nx[0];
	int bigNY = s_ny[0];

	//--------------------------------------------------------------------
	// Create and initialize each data type
	//--------------------------------------------------------------------
	for(int i=0;i<2;i++){
	for(int j=0;j<2;j++){
		
		initComm3d(0,0,0,bigNX-i,bigNY-j,myNZ,NX,NY,NZ,&mysubarray[   i*2+j]);
		initComm2d(0,0,  bigNX-i,bigNY-j,     NX,NY,   &mysubarray_2d[i*2+j]);
	}}

	//--------------------------------------------------------------------
	// Link each element of *mysubarrayptr to one of the four data types
	//--------------------------------------------------------------------
	for(int i=0;i<dims[0];i++){
	for(int j=0;j<dims[1];j++){

		int rank_index = i*dims[1]+j;

		for(int x=0;x<2;x++){
		for(int y=0;y<2;y++){

			if(s_nx[i]==bigNX-x && s_ny[j]==bigNY-y){
				
				mysubarrayptr[   rank_index] = mysubarray[   x*2+y];
				mysubarrayptr_2d[rank_index] = mysubarray_2d[x*2+y];
			}
		}}
		
		if(mysubarrayptr[rank_index]==NULL) { printf("Error\n");}

	}}
}

/*********************************************************************
* Exchange the adjacent northern and southern boundaries between
* processors
*
* @param var - buffer to exchange data
**********************************************************************/
inline void nsExchange(double * var){

	MPI_Sendrecv(
		var+I2(halo_buffer,halo_buffer,0),      1, nsBound, south, 123, 
		var+I2(halo_buffer,halo_buffer+myNY,0), 1, nsBound, north, 123,
		comm_cart, &status
	);
				 
	MPI_Sendrecv(
		var+I2(halo_buffer,myNY,0), 1, nsBound, north, 123, 
		var+I2(halo_buffer,0,0), 1,    nsBound, south, 123,
		comm_cart, &status
	);

}

/*********************************************************************
* Exchange the adjacent northern and southern boundaries between
* processors
*
* @param var - buffer to exchange data
**********************************************************************/
inline void nsExchangeOnePoint(double * var){

	MPI_Sendrecv(
		var+I2(halo_buffer,halo_buffer,0),      1, nsBoundOnePoint, south, 123, 
		var+I2(halo_buffer,halo_buffer+myNY,0), 1, nsBoundOnePoint, north, 123,
		comm_cart, &status
	);
				 
	MPI_Sendrecv(
		var+I2(halo_buffer,myNY+halo_buffer-1,0), 1, nsBoundOnePoint, north, 123, 
		var+I2(halo_buffer,halo_buffer-1,0), 1,    nsBoundOnePoint, south, 123,
		comm_cart, &status
	);

}

/*********************************************************************
* Exchange the adjacent eastern and western boundaries between
* processors
*
* @param var - buffer to exchange data
**********************************************************************/
inline void ewExchange(double * var){

	MPI_Sendrecv(
		var+I2(halo_buffer,halo_buffer,0), 		1, ewBound, west, 123, 
		var+I2(halo_buffer+myNX,halo_buffer,0), 1, ewBound, east, 123,
		comm_cart, &status
	);
		
	MPI_Sendrecv(
		var+I2(myNX,halo_buffer,0), 1, ewBound, east, 123, 
		var+I2(0,halo_buffer,0), 	1, ewBound, west, 123,
		comm_cart, &status
	);

}

/*********************************************************************
* Exchange the adjacent eastern and western boundaries between
* processors
*
* @param var - buffer to exchange data
**********************************************************************/
inline void ewExchangeOnePoint(double * var){

	MPI_Sendrecv(
		var+I2(halo_buffer,halo_buffer,0), 		1, ewBoundOnePoint, west, 123, 
		var+I2(halo_buffer+myNX,halo_buffer,0), 1, ewBoundOnePoint, east, 123,
		comm_cart, &status
	);
		
	MPI_Sendrecv(
		var+I2(myNX+halo_buffer-1,halo_buffer,0), 1, ewBoundOnePoint, east, 123, 
		var+I2(halo_buffer-1,halo_buffer,0), 	1, ewBoundOnePoint, west, 123,
		comm_cart, &status
	);

}

/*********************************************************************
* Exchange the adjacent corners between processors
*
* @param var - buffer to exchange data
**********************************************************************/
inline void cornerExchange(double * var){

	MPI_Sendrecv(
		var+I2(halo_buffer,halo_buffer,0),          myNZ,MPI_DOUBLE, sw_corner, 123, 
		var+I2(halo_buffer+myNX,halo_buffer+myNY,0),myNZ,MPI_DOUBLE, ne_corner, 123, 
		comm_cart, &status
	);
		
	MPI_Sendrecv(
		var+I2(halo_buffer+myNX-1,halo_buffer+myNY-1,0), myNZ,MPI_DOUBLE, ne_corner, 123, 
		var+I2(halo_buffer-1,halo_buffer-1,0),			 myNZ,MPI_DOUBLE, sw_corner, 123,
		comm_cart, &status
	);
	
	MPI_Sendrecv(var+I2(
		halo_buffer+myNX-1,halo_buffer,0), 		 myNZ,MPI_DOUBLE, se_corner, 123, 
		var+I2(halo_buffer-1,halo_buffer+myNY,0),myNZ,MPI_DOUBLE, nw_corner, 123,
		comm_cart, &status
	);
		
	MPI_Sendrecv(var+I2(
		halo_buffer,halo_buffer+myNY-1,0), 		 myNZ,MPI_DOUBLE, nw_corner, 123, 
		var+I2(halo_buffer+myNX,halo_buffer-1,0),myNZ,MPI_DOUBLE, se_corner, 123,
		comm_cart, &status
	);
	
}

/*********************************************************************
* Exchange all boundaries between adjacent processes
*
**********************************************************************/
void exchange(double * var){

	nsExchange(var);
	ewExchange(var);
	cornerExchange(var);
}

/*********************************************************************
* Exchange all boundaries between adjacent processes
*
**********************************************************************/
void exchangeOnePoint(double * var){

	nsExchangeOnePoint(var);
	ewExchangeOnePoint(var);
}

/*********************************************************************
*
*
**********************************************************************/
void get_corner(int xshift,int yshift,int coords[2],int dims[2],int * corner){

	int shifted[2];

	shifted[0] = coords[0] + xshift;
	shifted[1] = coords[1] + yshift;

	if(PERIODIC_BOUNDARIES){
		
		if(!PERIODIC_BOUNDARIES_NS){
		
			if(shifted[1]!=-1 && shifted[1]<dims[1]){
		
				MPI_Cart_rank(comm_cart,shifted,corner);
			}
			
		} else {
			MPI_Cart_rank(comm_cart,shifted,corner);
		}

	} else {
		
		if(shifted[0]!= -1 && shifted[1]!=-1 && shifted[0]<dims[0] && shifted[1]<dims[1]){
		
			MPI_Cart_rank(comm_cart,shifted,corner);
		}
	}
}

