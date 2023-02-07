/********************************************************
* Boundary conditions
*********************************************************/
extern int iwbuffer;	// number of points in x-direction at which damping will be applied
extern int iebuffer;
extern int jnbuffer;
extern int jsbuffer;

void apply_boundary_condition(int);
void apply_boundary_condition_microphysics(int);

void init_boundaries(int eastLength,int westLength,int northLength,int southLength,int halo);
void mirror_boundaries(double * s);
void periodic_boundaries();
void mirror_boundaries_2d(double * s);
void p_mirror_boundaries(double*,int,int,int,int,int,int,int);
void upper_lower_boundaries(double *);
void periodic_ew_sponge_ns_boundaries();
void periodic_uvw_boundaries();
void periodic_pressure_boundaries();
void zero_boundaries();