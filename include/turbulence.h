
extern double *landsea;
extern double *u_friction,*v_friction,*w_friction,*t_diffusion;
extern double *qv_diffusion,*qc_diffusion,*qr_diffusion;

void initialize_landsea(const char*);
void init_kmix(int,int,int,double *);
void apply_velocity_diffusion(double,int il,int ih,int jl,int jh);
void apply_scalar_diffusion(double,int il,int ih,int jl,int jh);
void apply_moisture_diffusion(double,int il,int ih,int jl,int jh);
void calculate_diff_tend(int il,int ih,int jl,int jh);