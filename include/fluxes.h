/********************************************************
* Indexing for flux control volume data structures
*********************************************************/
// velocity control volumes
#define UCELL(j,k)    ucell[NZ*(j)+k]
#define VCELL(j,k)    vcell[NZ*(j)+k]
#define WCELL(j,k)    wcell[NZ*(j)+k]
#define UBCELL(j,k)  ubcell[NZ*(j)+k]
#define VBCELL(j,k)  vbcell[NZ*(j)+k]
#define WBCELL(j,k)  wbcell[NZ*(j)+k]

// potential temperature control volumes
#define THCELL(j,k)   thcell[NZ*(j)+k]
#define THBCELL(j,k) thbcell[NZ*(j)+k]

// used in subroutine arguments
#define SCELL(j,k) scell[NZ*(j)+k]
#define BCELL(j,k) bcell[NZ*(j)+k]

// moisture variable control volumes
#define QVCELL(j,k)   qvcell[NZ*(j)+k]
#define QCCELL(j,k)   qccell[NZ*(j)+k]
#define QRCELL(j,k)   qrcell[NZ*(j)+k]
#define QVBCELL(j,k) qvbcell[NZ*(j)+k]

// sign of advecting velocity
#define SIGN_P_VEL(i,j,k) sign_p_vel[INDEX(i,j,k)]
#define SIGN_B_VEL(i,j,k) sign_b_vel[INDEX(i,j,k)]

/********************************************************
* Data structure for interpolating scalar fluxes onto the
* faces of a control volume.
*********************************************************/
struct cell {

	double top;
	double north;
	double east;
	double west;
};

/********************************************************
* Data structure for interpolating velocity fluxes onto the
* faces of a control volume.
*********************************************************/
struct vel_cell {

	// advected velocity or fluxes
	double top;
	double north;
	double east;
	double west;

	// advecting velocity
	double ua;
	double va;
	double wa;

	// normal horizontal velocity component
	double other_vel;
	
};

/********************************************************
* Data structure for sign of velocity components
*********************************************************/
struct sign_vel {

	double u;
	double v;
	double w;
};

extern struct vel_cell *ucell,*vcell,*wcell;
extern struct vel_cell *ubcell,*vbcell,*wbcell;
extern struct cell *thcell,*thbcell;
extern struct cell *qvcell,*qvbcell,*qccell,*qrcell,*ice_cell;
extern struct sign_vel *sign_p_vel,*sign_b_vel;

/********************************************************
* Interpolation and flux procedures for advection equations
*********************************************************/
void initialize_flux_cells(int ny,int nz);
void initialize_microphysics_cells(int ny,int nz);
void initialize_sign_cells(int nx, int ny,int nz);

void compute_sign_cells(int il,int ih,int jl,int jh);

void compute_fluxes(int i,int jl,int jh);
void compute_west(int i,int jl,int jh);
void compute_west_uv(int i,int jl,int jh);

void compute_fluxes_scalar(int i,int jl,int jh,double *,struct cell*);
void compute_fluxes_scalar(int i,int jl,int jh,double *,double *,struct cell*,struct cell*);
void compute_fluxes_scalar(int i,int jl,int jh,double *u,double *v,double *w,double *s,struct cell *scell);
void compute_fluxes_scalar_with_fallspeed(int i,int jl,int jh,double *s,double *fall,struct cell *scell);
void compute_west_scalar(int,int jl,int jh,double *s,double *sb,struct cell*,struct cell*);
void compute_west_scalar(int,int jl,int jh,double *s,struct cell*);
void compute_west_scalar(int i,int jl,int jh,double *u,double *s,struct cell *scell);

void compute_fluxes_moisture(int i,int jl,int jh);
void compute_west_moisture(int i,int jl,int jh);
