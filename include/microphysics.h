
#define ALPHA(x) ( 1.0 + (x-273.15) / (tmax-tmin) )
#define SAT_VAP_WAT(t) ( 611.2 * exp(17.67 * (t-273.15) / (t - 29.65) ) )
#define SAT_VAP_ICE(t) ( 611.2 * exp(21.8745584 * (t-273.15) / (t - 7.66) ) )
#define SAT_MIX_RATIO(e,p) ( 0.62197 * e / (p-e) )

const double tmin = -30+273.15;	// temperature below which only ice
const double tmax = 0+273.15;	// temperature above which only water

// hydrostatic uses Exner pressure
// non-hydrostatic uses pressure / density
#if HYDROSTATIC
    #define CONVERT_PRESSURE(i,j,k) PI(i,j,k)
#else
    #define CONVERT_PRESSURE(i,j,k) PI(i,j,k)/(cp*tbv[k])
#endif

/********************************************************
* Perturbation values forecast by model
*
*********************************************************/
extern double *qvps,*qvs,*qvms;		// vapor (kg/kg)
extern double *qcps,*qcs,*qcms;		// cloud water (kg/kg)
extern double *qrps,*qrs,*qrms;		// rain (kg/kg)
extern double *qsps,*qss,*qsms;
extern double *qips,*qis,*qims;
extern double *qgps,*qgs,*qgms;
extern double *nips,*nis,*nims;
extern double *nrps,*nrs,*nrms;
extern double *vts,*sts,*its;		// fall speed (m/s)
extern double *accRain,*accSnow;
extern double *dBZ;
extern bool *isSaturated;

#define QV(i,j,k) qvs[INDEX(i,j,k)]
#define QC(i,j,k) qcs[INDEX(i,j,k)]
#define QR(i,j,k) qrs[INDEX(i,j,k)]
#define QI(i,j,k) qis[INDEX(i,j,k)]
#define QS(i,j,k) qss[INDEX(i,j,k)]
#define QG(i,j,k) qgs[INDEX(i,j,k)]

#define QVP(i,j,k) qvps[INDEX(i,j,k)]
#define QCP(i,j,k) qcps[INDEX(i,j,k)]
#define QRP(i,j,k) qrps[INDEX(i,j,k)]
#define QIP(i,j,k) qips[INDEX(i,j,k)]
#define QSP(i,j,k) qsps[INDEX(i,j,k)]
#define QGP(i,j,k) qgps[INDEX(i,j,k)]

#define QVM(i,j,k) qvms[INDEX(i,j,k)]
#define QCM(i,j,k) qcms[INDEX(i,j,k)]
#define QRM(i,j,k) qrms[INDEX(i,j,k)]
#define QIM(i,j,k) qims[INDEX(i,j,k)]
#define QSM(i,j,k) qsms[INDEX(i,j,k)]
#define QGM(i,j,k) qgms[INDEX(i,j,k)]

#define QBAR(i,j,k) m_qbar[INDEX(i,j,k)]

#define VT(i,j,k) vts[INDEX(i,j,k)]
#define ST(i,j,k) sts[INDEX(i,j,k)]
#define IT(i,j,k) its[INDEX(i,j,k)]

#define ACCRAIN(i,j) accRain[INDEX2D(i,j)]

#if !PARALLEL && !ENERGY
	#define IQBAR(i,j,k) QBAR(i,j,k)
#else
	#define IQBAR(i,j,k) iqbar[(i)*NY*NZ+(j)*NZ+(k)]
#endif

//extern double kmixv_moisture[NZ];

/********************************************************
* Driver functions
*
*********************************************************/
void run_microphysics(int,int,int,int);
void init_microphysics(int,int);

/********************************************************
* General shared function
*
*********************************************************/
double get_dimensional_pressure(int index,int k);
double get_qvsat_mixed(double temperature,double pd);
void microphysics_diffusion(int,int,int,int,double);
void microphysics_advance(size_t);
void microphysics_advance_inner(size_t);
void zero_moisture(int il,int ih,int jl,int jh,int size);
double get_CAPE(int i,int j,int k_p,double,double);
double get_CAPE(int i,int j,int k_p);
double get_CAPE_base(int i,int j,int k_p);
void hydrometeor_fallout(double *var,double *vel,int il,int ih,int jl,int jh,double *precip);
void precip_rate(int il,int ih,int jl,int jh,double *,double *,double *);
void calculate_rain_fall_speed_kessler(double *vt, double *qr, int il,int ih,int jl,int jh);
void calculate_rain_fall_speed_rutledge(double *, double *, int il,int ih,int jl,int jh);
void calculate_snow_fall_speed_rutledge(double *, double *, int il,int ih,int jl,int jh);
void calculate_ice_fall_speed_rutledge(double *, double *, int il,int ih,int jl,int jh);
void calculate_eulerian_fall_speed_rain(double *vt, double *qr, int il,int ih,int jl,int jh);
void calculate_eulerian_fall_speed_snow_ice(double *vts, double *qs, double *vti, double *qi, int il,int ih,int jl,int jh);

/********************************************************
* Specific schemes
*
*********************************************************/
void init_thompson_microphysics(int mpirank);
void run_rutledge_microphysics(int il,int ih,int jl,int jh);
void run_kessler_microphysics(int il,int ih,int jl,int jh);
void run_thompson_microphysics(int il,int ih, int jl, int jh, int kl, int kh,
    int nx, int ny, int nz,
    double *qv,double *qc,double *qi,double *qr,double *qs,double *qg,
    double *ni,double *nr,
    double *t3d,double *p3d,double *w3d,
    double *pptrain,double *pptsnow,bool do_radar);

