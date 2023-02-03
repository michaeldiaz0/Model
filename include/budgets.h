
extern double *m_diabatic,*q_diabatic,*cond,*evap;//,*qm_diabatic,*pe_m_diabatic;
extern double *vort_ufric,*vort_vfric;
extern double *u_bound,*v_bound,*w_bound,*t_bound;

extern char heat_budget_filename[len];
extern char moisture_budget_filename[len];
extern char pe_budget_filename[len];
extern char pv_budget_filename[len];
extern char vorticity_budget_filename[len];

void initialize_budgets();
void calculate_budgets(int,double *);
void write_budgets_to_file();

void calculate_omega_driver(const char *infilename,const char *outfilename,int *times,int count);

void pv_tracer_sources(int il,int ih,int jl,int jh);

