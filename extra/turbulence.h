

extern double tke[NX][NY][NZ],tkep[NX][NY][NZ],tkem[NX][NY][NZ];

void integrate_tke(double);
void turbulent_diffusion(double);
void initialize_tke();
