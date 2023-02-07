
#include "pvars.h"



int myNX,myNY,myNZ;	// dimensions of subarrays without halo boundaries
int fNX,fNY,fNZ;	// full dimensions of subarrays
int pNX,pNY,pNZ;	// dimensions of arrays for pressure solver
int fNYfNZ;

int coords[2];


int error_status = 0;

int rank,row_rank,col_rank;
int dims[] = {0,0};
int numtasks;
int row_size;
int col_size;

// index within the larger domain
int *big_i,*big_j;

// sizes of subarrays on each process
int *s_nx,*s_ny,*s_nz,*s_ny_p,*s_nx_p;

// starting and ending indices of each process within the full domain arrays
int *ibs,*ibe,*jbs,*jbe,*kbs,*kbe,*jbs_p,*ibs_p;


double *pres_row,*ipres_row;
double *pres_col,*ipres_col;

double *frictions;
double *istopos,*uistopos,*vistopos;

int size,size2;
