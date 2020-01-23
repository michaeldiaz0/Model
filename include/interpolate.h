
#define d(i,j,k) (xdim*ydim*(k) + xdim*(j) + i)
#define d2(i,j) (xdim*(j) + i)
#define d3(i,j,k) (zdim*ydim*(i) + zdim*(j) + k)
#define d4(i,j,k) (zdim*NY*(i) + zdim*(j) + k)
#define d5(i,j,k) (zdim_out*ydim*(i) + zdim_out*(j) + k)

/********************************************************
* Interpolation procedures for initialization
*********************************************************/
void horz_interpolate(double *,double *,int,int,int,bool,bool,double,double,double,double);
void horz_interpolate_from_model(double *,double *,int,int,int,bool,bool,double,double,double,double);
double * vert_interpolate_1d(double *z, double *source_z, double *source_val, int xdim, int ydim, int zdim,int zdim_out);
//void vert_interpolate_1d(double[NX][NY][NZ], double[NZ],double[NX][NY][NZ]);
void vert_interpolate_1d_from_model(double *z, double *source_z, double *source_val, int xdim, int ydim, int zdim,int zdim_out,double *out);
void vert_interpolate_1d(int nzIn,int nzOut, double * zIn,double * zOut,double * varIn,double * varOut);
//void vert_interpolate_3d(double[NX][NY][NZ], double[NX][NY][NZ],double[NX][NY][NZ]);
double* vert_interpolate(double *, double *, double *, int, int,int,int);
double CubicInterpolate(double,double,double,double,double,double,double,double,double);
double CubicInterpolate2(double,double,double,double,double);
double BicubicInterpolate2(double[4][4],double,double);
double CosineInterpolate(double,double,double);
double LinearInterpolate(double,double,double);
void vert_interpolate_1d_linear(int nzIn,int nzOut, double * zIn,double * zOut,double * varIn,double * varOut);