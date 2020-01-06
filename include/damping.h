
extern double *u_diff_tend;
extern double *v_diff_tend;

/********************************************************
* Damping and diffusion
*********************************************************/
void init_damping(int ni,int nj,int nk);
void fft_damp(int ni,int nj,int nk,int,double *var);
void fft_damp2(int ni,int nj,int nk,int,double *var);
void asselin();
void diffusion(double step);
void diffusion_8th_all(double step,int il,int ih,int jl,int jh);
void diffusion_6th_all(double step,int il,int ih,int jl,int jh);
void diffusion_4th_all(double step,int il,int ih,int jl,int jh);
void diffusion_2nd_all(double step,int il,int ih,int jl,int jh);
void rayleigh_damping(int,int,int,int,int);
void damp_var(double *var,int il,int ih,int jl, int jh,double coef,double max);

void diffusion_6th_var(double *pvar,double *mvar,double step,int il,int ih,int jl,int jh,double *output=NULL);
void diffusion_2nd_var(double *pvar,double *mvar,double step,int il,int ih,int jl,int jh);

double diffuse_i_4th(double *var,int i,int j,int k);
double diffuse_j_4th(double *var,int i,int j,int k);
double diffuse_k_4th(double *var,int i,int j,int k);
double diffuse_i_6th(double *var,int i,int j,int k);
double diffuse_j_6th(double *var,int i,int j,int k);
double diffuse_k_6th(double *var,int i,int j,int k);

void run_fftw(int ni,int nj,int nk,int b);
void init_fftw(int ni,int nj,int nk);

void init_large_scale_precip();
void large_scale_precip(int il,int ih,int jl,int jh);