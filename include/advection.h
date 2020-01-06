
void advect_theta(double,int,int,int,int);
void advect_qv(double,int,int,int,int);
void advect_uv_velocity(double,int,int,int,int);
void advect_uvw_velocity(double,int,int,int,int);
void advect_microphysics_cell(double,int,int,int,int);
void advect_ice_cell(double step,int il,int ih,int jl,int jh);
void advect_scalar(double,int,int,int,int,double *,double *,double *,struct cell *);
void advect_scalar(double step,int il,int ih,int jl,int jh,double *u,double *v,double *w,double *var,double *tend,struct cell *scell);
void advect_scalar(double step,int il,int ih,int jl,int jh,double *varm,double *var,double *varp,double *varbase,struct cell *scell,struct cell *bcell);