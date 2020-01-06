
extern int argcount;
extern double hor_shear;
extern double vert_shear0;
extern double vert_shear1;

/*********************************************************
*
* 			VARIABLES FOR METEOROLOGICAL FIELDS
*
**********************************************************/

//--------------------------------------------------------
// Array indices
//--------------------------------------------------------
extern int NYNZ;

#if PARALLEL && !ENERGY
	#define INDEX(i,j,k) ((i)*fNYfNZ + (j)*fNZ + (k))	// index for parallel version
#else
	#define INDEX(i,j,k) ((i)*NYNZ + (j)*NZ + (k))		// index for serial version
#endif													

#define FULL_ARRAY_INDEX(i,j,k) ((i)*NYNZ+(j)*NZ+(k))	// index for full domain arrays
#define FULL_INDEX(i,j,k) ((i)*NYNZ+(j)*NZ+(k))

//--------------------------------------------------------
// Perturbation values forecast by model
//--------------------------------------------------------
extern double *us,	*vs,	*ws,	*ths ;
extern double *ups,	*vps,	*wps,	*thps;
extern double *ums,	*vms,	*wms,	*thms;
extern double *pis;

#define UP(i,j,k) ups[INDEX(i,j,k)]		// pert. zonal wind, forecast value
#define U(i,j,k)  us [INDEX(i,j,k)]		// pert. zonal wind, middle of time step
#define UM(i,j,k) ums[INDEX(i,j,k)]		// pert. zonal wind, beginning of time step
#define VP(i,j,k) vps[INDEX(i,j,k)]		// pert. meridional wind, forecast value
#define V(i,j,k)  vs [INDEX(i,j,k)]		// pert. meridional wind, middle of time step
#define VM(i,j,k) vms[INDEX(i,j,k)]		// pert. meridional wind, beginning of time step
#define WP(i,j,k) wps[INDEX(i,j,k)]		// pert. vertical wind, forecast value
#define W(i,j,k)  ws [INDEX(i,j,k)]		// pert. vertical wind, middle of time step
#define WM(i,j,k) wms[INDEX(i,j,k)]		// pert. vertical wind, beginning of time step
#define THP(i,j,k) thps[INDEX(i,j,k)]	// pert. potential temperature, forecast value
#define TH(i,j,k)  ths [INDEX(i,j,k)]	// pert. potential temperature, middle of time step
#define THM(i,j,k) thms[INDEX(i,j,k)]	// pert. potential temperature, beginning of time step
#define PI(i,j,k) pis[INDEX(i,j,k)]		// pert. pressure, diagnostic

//--------------------------------------------------------
// Basic state while model is running
//--------------------------------------------------------
extern double *m_ubar,*m_vbar,*m_wbar,*m_thbar,*m_qbar,*m_pbar;

#define UBAR(i,j,k) m_ubar[INDEX(i,j,k)]	// basic state zonal wind
#define VBAR(i,j,k) m_vbar[INDEX(i,j,k)]	// basic state meridional wind
#define WBAR(i,j,k) m_wbar[INDEX(i,j,k)]	// basic state vertical wind
#define THBAR(i,j,k) m_thbar[INDEX(i,j,k)]	// basic state potential temperature (departure from base state (tb))
#define PBAR(i,j,k) m_pbar[INDEX(i,j,k)]	// basic state pressure (non-dimensional units, full value)

//--------------------------------------------------------
// Basic state for initialization
//--------------------------------------------------------
extern double *iubar,*ivbar,*iwbar,*ithbar,*iqbar,*ipbar;
//---------------------------------------------
// For serial version, initialize basic
// state in its final location
//---------------------------------------------
#if !PARALLEL && !ENERGY
	#define IUBAR(i,j,k) UBAR(i,j,k)
	#define IVBAR(i,j,k) VBAR(i,j,k)
	#define IWBAR(i,j,k) WBAR(i,j,k)
	#define ITHBAR(i,j,k) THBAR(i,j,k)
	#define IPBAR(i,j,k) PBAR(i,j,k)
//---------------------------------------------
// For parallel version, initialize basic
// state in a temporary location to be
// distributed to other processes
//---------------------------------------------
#else
	#define IUBAR(i,j,k)   iubar[FULL_ARRAY_INDEX(i,j,k)]
	#define IVBAR(i,j,k)   ivbar[FULL_ARRAY_INDEX(i,j,k)]
	#define IWBAR(i,j,k)   iwbar[FULL_ARRAY_INDEX(i,j,k)]
	#define ITHBAR(i,j,k) ithbar[FULL_ARRAY_INDEX(i,j,k)]
	#define IPBAR(i,j,k)   ipbar[FULL_ARRAY_INDEX(i,j,k)]
#endif

//--------------------------------------------------------
// Topography and friction
//--------------------------------------------------------
extern double *istopos,*uistopos,*vistopos;
extern int *htopo;//[NX][NY];
extern double *frictions;

//--------------------------------------------------------
// Map coordinates and Coriolis terms
//--------------------------------------------------------
extern double *outLats;//[NY];		// latitudes
extern double *outLons;//[NX];		// longitudes
extern double *f;//f[NY];			// coriolis parameter
extern double *dfdy;//[NY];			// gradient of coriolis parameter

extern double *output_to_file;//u[NX][NY][NZ];
extern double *rhoavg2d;//[NX][NY];

extern double *itopo,*iistopo,*iuistopo,*ivistopo,*ifriction;
extern double *rate;

#define ISTOPO(i,j,k) istopos[INDEX(i,j,k)]
#define UISTOPO(i,j,k) uistopos[INDEX(i,j,k)]
#define VISTOPO(i,j,k) vistopos[INDEX(i,j,k)]
#define FRICTION(i,j,k) frictions[INDEX(i,j,k)]
#define HTOPOFULL(i,j) htopo[(i)*NY+(j)]
#define RHOAVG2DFULL(i,j) rhoavg2d[(i)*NY+(j)]

//--------------------------------------------------
// Arrays for serial version
//--------------------------------------------------
#if !PARALLEL && !ENERGY
	
	#define ITOPO(i,j,k) itopo[(i)*NY*NZ+(j)*NZ+(k)]
	#define IISTOPO(i,j,k) ISTOPO(i,j,k)
	#define IUISTOPO(i,j,k) UISTOPO(i,j,k)
	#define IVISTOPO(i,j,k) VISTOPO(i,j,k)
	#define IFRICTION(i,j,k) FRICTION(i,j,k)
	#define LANDSEA(i,j) landsea[(i)*NY+(j)]
	#define HTOPO(i,j) HTOPOFULL(i,j)
	#define RHOAVG2D(i,j) RHOAVG2DFULL(i,j)
	#define FC(j) f[j]
	#define DFDY(j) dfdy[j]
	#define LATS(j) outLats[j]
	#define LONS(i) outLons[i]
//---------------------------------------------------
// Arrays for parallel version
//---------------------------------------------------
#else
	#define ITOPO(i,j,k) itopo[(i)*NY*NZ+(j)*NZ+(k)]
	#define IISTOPO(i,j,k) iistopo[(i)*NY*NZ+(j)*NZ+(k)]
	#define IUISTOPO(i,j,k) iuistopo[(i)*NY*NZ+(j)*NZ+(k)]
	#define IVISTOPO(i,j,k) ivistopo[(i)*NY*NZ+(j)*NZ+(k)]
	#define IFRICTION(i,j,k) ifriction[(i)*NY*NZ+(j)*NZ+(k)]

	#define LANDSEA(i,j) landsea[(i+ibs[rank]-3)*NY+(j+jbs[rank]-3)]
	#define HTOPO(i,j) htopo[(i+ibs[rank]-3)*NY+(j+jbs[rank]-3)]//htopo[i+ibs[rank]-3][j+jbs[rank]-3]
	#define RHOAVG2D(i,j) rhoavg2d[(i+ibs[rank]-3)*NY+(j+jbs[rank]-3)]
	#define LATS(j) outLats[j+jbs[rank]-3]
	#define LONS(i) outLons[i+ibs[rank]-3]

	#if ENERGY
		#define FC(j) f[j] 
		#define DFDY(j) dfdy[j]
	#else
		#define FC(j) f[j+jbs[rank]-3] 
		#define DFDY(j) dfdy[j+jbs[rank]-3]
	#endif
#endif

#if STRETCHED_GRID
	#define ZU(k) zsu[k]
	#define ZW(k) zsw[k]
	#define DTZ(k) (dtz*mu[k])
	#define DTZW(k) (dtz*mw[k])
	#define DZU(k) (zsw[k+1]-zsw[k])
	#define DZW(k) (zsu[k]-zsu[k-1])
	#define ONE_D_DZ(k)  (one_d_dz*mu[k])
	#define ONE_D_DZW(k) (one_d_dz*mw[k])
#else
	#define ZU(k) zu[k]
	#define ZW(k) zw[k]
	#define DTZ(k) dtz
	#define DTZW(k) dtz
	#define DZU(k) dz
	#define DZW(k) dz
	#define ONE_D_DZ(k)  one_d_dz
	#define ONE_D_DZW(k) one_d_dz
#endif

extern double dtx;// = dt/dx;
extern double dty;// = dt/dy;
extern double dtz;// = dt/dz;

extern double one_d_dx;// = 1.0/dx;
extern double one_d_dy;// = 1.0/dy;
extern double one_d_dz;// = 1.0/dz;

/********************************************************
* Model basic state
*********************************************************/
extern double *zu;		// height of u-velocity
extern double *zw;	// height of w-velocity
extern double *rhou;	// base state density at u-velocity level
extern double *rhow;	// base state density at w-velocity level
extern double *tb;	// base state potential temperature
extern double *tbv;	// base state virtual potential temperature
extern double *tbw;	// base state potential temperature at w-velocity level
extern double *qb;	// base state specific humidity
extern double *pib;	// base state non-dimensional pressure
extern double *one_d_rhou;
extern double *one_d_rhow;
extern double *zsu;
extern double *zsw;
extern double *mu;
extern double *mw;

/********************************************************
* Physical constants
*********************************************************/
const double Rd = 287.0;		// gas constant for dry air (J/kgK)
const double cp = 1004;			// specific heat at constant pressure for air (J/K)
const double cv = 717;			// specific heat at constant volume for air (J/K)
const double p0 = 100000;		// non-dimensional pressure scaling parameter (Pa)
const double grav = 9.8;		// gravitational constant
const double Lv = 2.5e6;		// latent heat of vaporization (J/kg)
const double Ls = 2.85e6;		// latent heat of sublimation (J/kg)
const double Lf = 3.50e5;		// latent heat of fusion (J/kg)
const double trigpi = 3.14159265358979324;
const double pressfc = 100200; // dimensional surface pressure
const double fc = 7.292e-5;		// coriolis parameter
const double meters_per_degree = 111000;
const double R_earth = 6.3781e6;

/********************************************************
* Other
*********************************************************/
const double meridional_lon_section = 87.0;

extern double mtime;
//extern double sor[3][NX][NY];
extern double rhoavg;
extern int bigcounter;
extern size_t file_time_counter;
