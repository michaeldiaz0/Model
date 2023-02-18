extern int *scounts,*rcounts,*sdisps,*rdisps;

#if PARALLEL
extern MPI_Datatype *types;
extern MPI_Comm comm_cart,comm_row,comm_col;
#endif

extern int dims[];

// indices for neighboring processes
extern int west,east,south,north;
extern int nw_corner,ne_corner,sw_corner,se_corner;

/*********************************************
* 
*********************************************/
// which row or column a process is in
extern int row_rank,col_rank;

/*********************************************
* every process will know this grid information 
* about every other process
*********************************************/
extern int *s_nx,*s_ny,*s_nz,*s_nx_p,*s_ny_p;	// grid dimensions of each process
extern int *ibs,*ibe,*jbs,*jbe,*kbs,*kbe,*ibs_p,*jbs_p;	// beginning and ending i,j coordinates of each process

extern int numtasks;	// number of processes
extern int row_size;	// processes per row
extern int col_size;	// processes per column

extern int rank;
const int halo_buffer = 3; // number of boundary points for each process

struct val_rank { 
    double val; 
    int rank; 
};

/*********************************************
* Subroutines
*********************************************/
void send_int_from_process_to_root(int *sdata,int *rdata,int process,int count);
void get_max_mpi(val_rank *send,val_rank *receive,int count);
void distributeArray_2d(double *,double *);
void distributeArray_3d(double *,double *);
void pres_row_alltoall(double *);
void pres_row_alltoall2(double *,double *);
void pres_row_alltoall2_reverse(double *,double *);
void pres_col_alltoall(double *);
void pres_col_alltoall2(double *,double *);
void pres_col_alltoall_reverse(double *);
void pres_vert_alltoall(double *);
void pres_vert_alltoall_reverse(double *);
void gatherArrays_3d(double * m_var,double * var);
void gatherArrays_2d(double * m_var,double * var);
void presComm2d();
void initAllComms();
void exchange(double *);
void exchangeOnePoint(double *);
void get_corner(int,int,int [2],int [2],int *);
void exit_from_error();
void broadcast_data_int(int size,int *var);
void broadcast_data_bool(int size,int *var);
