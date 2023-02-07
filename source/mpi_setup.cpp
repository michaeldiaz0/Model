#include "stdafx.h"
#include "mpi.h"
#include "pcomm.h"
#include "boundaries.h"

#define BCAST_NX(var) MPI_Bcast(&var[0],NX,MPI_DOUBLE,0,MPI_COMM_WORLD)
#define BCAST_NY(var) MPI_Bcast(&var[0],NY,MPI_DOUBLE,0,MPI_COMM_WORLD)
#define BCAST_NZ(var) MPI_Bcast(&var[0],NZ,MPI_DOUBLE,0,MPI_COMM_WORLD)


// index within the larger domain
int *big_i,*big_j;

int coords[2];
int size,size2;
int error_status = 0;

/*********************************************************************
* Set up parallel environment
**********************************************************************/
void initialize_parallel_environment(int argc, char *argv[]){
	
	int len; 
	char hostname[MPI_MAX_PROCESSOR_NAME];
	
	MPI_Init(&argc,&argv);					// initialize MPI  
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);// get number of tasks 
	MPI_Get_processor_name(hostname, &len);
}


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
