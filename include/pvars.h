/*************************************************************************
*
* Variables and subroutines used in parallel solver
*
**************************************************************************/

/*********************************************
* Indexing
*********************************************/
#define I(i,j,k) ((i)*myNY*myNZ + (j)*myNZ + (k))	// within subarray interior points without boundaries
#define I2(i,j,k) ((i)*fNY*fNZ + (j)*fNZ + (k))	// subarrays with boundaries

extern int error_status;

/*********************************************
* Grid dimensions
*********************************************/
extern int myNX,myNY,myNZ;	// dimensions of subarrays without halo boundaries
extern int fNX,fNY,fNZ,pNX,pNY,pNZ;	// full dimensions of subarrays
extern int fNYfNZ;

/*********************************************
* every process will know this grid information 
* about every other process
*********************************************/
extern int *ibs,*ibe,*jbs,*jbe,*kbs,*kbe,*ibs_p,*jbs_p;	// beginning and ending i,j coordinates of each process
extern int *s_nx,*s_ny,*s_nz,*s_nx_p,*s_ny_p;	// grid dimensions of each process

/*********************************************
* 
*********************************************/
// which row or column a process is in
extern int rank,row_rank,col_rank;

const int halo_buffer = 3; // number of boundary points for each process
extern int rank,row_rank,col_rank;
extern int dims[];
extern int numtasks;	// number of processes
extern int row_size;	// processes per row
extern int col_size;	// processes per column

// indices for neighboring processes
extern int west,east,south,north;
extern int nw_corner,ne_corner,sw_corner,se_corner;


// rows and columns for pressure solver
extern double *pres_row;
extern double *pres_col;
extern double *pres_vert;
extern double *ipres_row;
extern double *ipres_col;
extern double *ipres_vert;

/*********************************************
* Subroutines
*********************************************/
void distributeArrays();
void distributeArray(double *);
void distributeArray(double *,double *);
//void alltoall(double *,double [NX][NY][NZ]);
void pres_row_alltoall(double *);
void pres_row_alltoall2(double *,double *);
void pres_row_alltoall2_reverse(double *,double *);
void pres_col_alltoall(double *);
void pres_col_alltoall2(double *,double *);
void pres_col_alltoall_reverse(double *);
void pres_vert_alltoall(double *);
void pres_vert_alltoall_reverse(double *);
void wait(int);
//void scatterArrays(double *,double [NX][NY][NZ]);
void gatherArrays();
//void gatherArrays2(double *,double [NX][NY][NZ]);
void gatherArrays3(double * m_var,double * var);
void presComm2d();
void initAllComms();
void exchange(double *);
void exchangeOnePoint(double *);
void get_corner(int,int,int [2],int [2],int *);
void run_parallel_model(int argc, char *argv[]);
void broadcast_data_int(int size,int *var);
void initialize_subarray(int size);
void exit_from_error();