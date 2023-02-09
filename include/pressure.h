// rows and columns for pressure solver
extern double *pres_row;
extern double *pres_col;
extern double *ipres_row;
extern double *ipres_col;

void solve_pressure(double);
void initialize_pressure_solver();