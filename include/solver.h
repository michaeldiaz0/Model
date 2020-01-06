void w_velocity_LH(int,int,int,int);
void integrate_hydro(double,int,int,int,int,int);
void integrate_non_hydro(double,int,int,int,int);
void run_model(int,FILE *infile=NULL);
void advance_inner(size_t);
void advance_outer(size_t);
void set_terrain(double *,double *,int);
void switch_array(double **var1,double **var2);