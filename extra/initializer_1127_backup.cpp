#include "stdafx.h"
#include "interpolate.h"
#include "surface.h"
#include "fluxes.h"
#include "damping.h"
#include "initializer.h"

#define d(i,j,k) (xdim*ydim*(k) + xdim*(j) + i)
#define ALLOC(v,size) v = (double *)calloc(size,sizeof(double))

/***************************************************************************
* -------------------------- VARIABLES -------------------------------------
****************************************************************************/
double u[NX][NY][NZ];
//--------------------------------------
// BASIC STATE ARRAYS
//--------------------------------------
double *m_ubar,*m_vbar,*m_thbar,*m_wbar,*m_qbar,*m_pbar;
double *iubar,*ivbar,*iwbar,*ithbar,*iqbar,*ipbar;
//--------------------------------------
// PERTURBATION ARRAYS
//--------------------------------------
double *us,*vs,*ths,*pis,*ws;
double *ups,*vps,*wps,*thps;
double *ums,*vms,*wms,*thms;
//--------------------------------------
// TOPOGRAPHY AND FRICTION ARRAYS
//--------------------------------------
double *itopo,*iistopo,*iuistopo,*ivistopo,*ifriction;
//--------------------------------------
// VERTICAL VARYING BASIC STATE
//--------------------------------------
double zu[NZ],zw[NZ],rhou[NZ],rhow[NZ],tb[NZ],tbw[NZ],tbv[NZ],qb[NZ],pib[NZ];
double one_d_rhou[NZ],one_d_rhow[NZ];
//--------------------------------------
// VERTICAL COORDINATE PARAMETERS
//--------------------------------------
double zsu[NZ],zsw[NZ];
double mu[NZ],mw[NZ];
//--------------------------------------
// OTHER ARRAYS
//--------------------------------------
double f[NY];
double dfdy[NY];

double rhoavg2d[NX][NY];
double outLats[NY];
double outLons[NX];
int htopo[NX][NY];
int bigcounter;
int linear_lon;

/***************************************************************************
* ---------------------- FUNCTION PROTOTYPES--------------------------------
****************************************************************************/
void init_friction();
void init_topography();
void initialize_from_output_serial(const char *,size_t);
void initialize_from_output(const char *,size_t);
void init_basic_state_vertical_velocity();
void init_basic_state_meridional_velocity(int);
void initialize_vertical_basic_state(int,int);
double get_QV_Sat(double temperature,double pressure);
void change_humidity();
void stretched_grid(double*,double*,double*,double*,double,int);
void reverse_yz_coord(double*, int, int, int);
void flip_array(double *, int, int, int);
void initialize_vortex(double r_end, double zm, double zt, double tmax, double tmin,int,int);
void initialize_from_ERA5(double *);

/******************************************************************************
*
*
*******************************************************************************/
void initialize_basic_state(){

	int size = NX*NY*NZ;

	rhoavg = 0;
	mtime = 0;
	bigcounter = 0;
	
	//----------------------------------------------------------------
	// initialize common data
	//----------------------------------------------------------------
	if(PARALLEL)
		memset(&u[0][0][0], 0,NX*NY*NZ*sizeof(double));

	//----------------------------------------------------------------
	// Calculate height of each level
	//----------------------------------------------------------------
	for(int k=0;k<NZ;k++){

		zu[k] = ((double)k-0.5)*dz;
		zw[k] = ((double)k-1)*dz;
	}

	stretched_grid(&zsu[0],&mu[0],&zsw[0],&mw[0],50,1);

	itopo = (double*) calloc(size,sizeof(double));

	//----------------------------------------------------------------
	// initialize data specific to parallel version
	//----------------------------------------------------------------
	#if PARALLEL || ENERGY
	
		ALLOC(iubar,size);
		ALLOC(ivbar,size);
		ALLOC(iwbar,size);
		ALLOC(ithbar,size);
		ALLOC(ipbar,size);
		ALLOC(iqbar,size);
	
		ALLOC(iuistopo,size);
		ALLOC(ivistopo,size);
		ALLOC(iistopo,size);
		ALLOC(ifriction,size);
	
		if(STRETCHED_GRID){initialize_from_era(&zsu[0]);} 
		else {			   initialize_from_era(&zu[0]);}
		
		if(ENERGY){
			
			ALLOC(us,size);
			ALLOC(vs,size);
			ALLOC(ws,size);
			ALLOC(ths,size);
			ALLOC(pis,size);

			//--------------------------------------------------
			// microphysics
			//--------------------------------------------------
			if(USE_MICROPHYSICS){

				ALLOC(qvs,size);
				ALLOC(qcs,size);
				ALLOC(qrs,size);
			}
		}
		
	//----------------------------------------------------------------
	// initialize data specific to serial version
	//----------------------------------------------------------------
	#else
		initialize_flux_cells(NY,NZ);
		initialize_microphysics_cells(NY,NZ);

		initialize_subarray(NX*NY*NZ);

		if(STRETCHED_GRID){ initialize_from_era(&zsu[0]);} 
		else { 				initialize_from_ERA5(&zu[0]);}

		initialize_landsea(landseaMaskFile);

		if(USE_TURBULENT_STRESS){ init_kmix(NX,NY,NZ);}
		if(OUTPUT_DIFFUSION_TEND){ init_damping(NX,NY,NZ);}	
	#endif
		
	//----------------------------------------------------------------
	// Initialize stuff for Fourier damping in linearized equations
	//----------------------------------------------------------------
	if(ISLINEAR && FOURIER_DAMPING){
		
		if(!PERIODIC_BOUNDARIES){ init_fftw(NX,NY,NZ);}
		else { init_fftw(NX-6,NY,NZ);}
	}
	
	//----------------------------------------------------------------
	// If requested, create zonally uniform basic state
	//----------------------------------------------------------------
	if(MERIDIONAL_CROSS_SECTION){

		//int linear_lat = get_point_from_lat(18);
		//int linear_lat = get_point_from_lat(24);
		linear_lon = get_point_from_lon(meridional_lon_section);

		printf("lon = %d\n",linear_lon);

		//for(int i=0;i<linear_lon;i++)
		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
		for(int k=0;k<NZ;k++){
	
			IUBAR(i,j,k) = IUBAR(linear_lon,j,k);
			IVBAR(i,j,k) = 0;
			ITHBAR(i,j,k) = ITHBAR(linear_lon,j,k);
			IQBAR(i,j,k) = IQBAR(linear_lon,j,k);
			IPBAR(i,j,k) = IPBAR(linear_lon,j,k);
		}}}

		//----------------------------------------------------------------
		// For linearized model, create meridionally uniform conditions on the
		// northern and southern boundaries	
		//----------------------------------------------------------------
		if(ISLINEAR){

			int b1 = 0;
			int b2 = 10;
	
			for(int i=0;i<NX;i++){
			for(int k=0;k<NZ;k++){

				for(int j=0;j<b1;j++){
				
					IUBAR(i,j,k) = IUBAR(i,b1,k);
					ITHBAR(i,j,k) = ITHBAR(i,b1,k);
					IQBAR(i,j,k) = IQBAR(i,b1,k);
				}

				for(int j=NY-b2;j<NY;j++){

					IUBAR(i,j,k) = IUBAR(i,NY-b2,k);
					ITHBAR(i,j,k) = ITHBAR(i,NY-b2,k);
					IQBAR(i,j,k) = IQBAR(i,NY-b2,k);
				}
			}}
			
			
		}
#if 0
		for(int t=0;t<10;t++){
		
			for(int i=0;i<NX;i++){
			for(int j=3;j<NY-3;j++){
			for(int k=0;k<NZ;k++){
		
				IVBAR(i,j,k) = 0.01*diffuse_j_6th(&IUBAR(0,0,0),i,j,k);
			}}}
		
			for(int i=0;i<NX;i++){
			for(int j=3;j<NY-3;j++){
			for(int k=0;k<NZ;k++){
		
				IUBAR(i,j,k) += IVBAR(i,j,k);
			}}}
			
			for(int i=0;i<NX;i++){
			for(int j=3;j<NY-3;j++){
			for(int k=0;k<NZ;k++){
		
				IVBAR(i,j,k) = 0.01*diffuse_j_6th(&ITHBAR(0,0,0),i,j,k);
			}}}
		
			for(int i=0;i<NX;i++){
			for(int j=3;j<NY-3;j++){
			for(int k=0;k<NZ;k++){
		
				ITHBAR(i,j,k) += IVBAR(i,j,k);
			}}}
		}
	
		for(int i=0;i<NX;i++){
		for(int j=3;j<NY-3;j++){
		for(int k=0;k<NZ;k++){
		
			IVBAR(i,j,k) = 0;
		}}}	
		
		for(int j=0;j<NY;j++){
			
			
			
			printf("%f %f\n",IUBAR(NX/2,j,8),ITHBAR(NX/2,j,8));
		}
		#endif
	}
	//exit(0);
	//init_large_scale_precip();
	
	//----------------------------------------------------------------
	// For resting basic state set the wind components to zero and make
	// the temperature and moisture fields horizontally uniform
	//----------------------------------------------------------------
	if(RESTING_BASIC_STATE){
		
		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
		for(int k=0;k<NZ;k++){

			IUBAR(i,j,k) = 0;
			IVBAR(i,j,k) = 0;
			ITHBAR(i,j,k) = ITHBAR(NX/2,NY/2,k);
			IQBAR(i,j,k) = IQBAR(NX/2,NY/2,k);
			IPBAR(i,j,k) = IPBAR(NX/2,NY/2,k);		
		}}}
	}
	
	//if(USE_MICROPHYSICS){ init_microphysics();}

	//----------------------------------------------------------------
	// Upper and lower boundary conditions for base state array
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		ITHBAR(i,j,0) = ITHBAR(i,j,1);
		ITHBAR(i,j,NZ-1) = ITHBAR(i,j,NZ-2);
		IQBAR(i,j,0) = IQBAR(i,j,1);
		IQBAR(i,j,NZ-1) = IQBAR(i,j,NZ-2);
		IUBAR(i,j,0) = IUBAR(i,j,1);
		IUBAR(i,j,NZ-1) = IUBAR(i,j,NZ-2);
		IVBAR(i,j,0) = IVBAR(i,j,1);
		IVBAR(i,j,NZ-1) = IVBAR(i,j,NZ-2);
	}}

	//----------------------------------------------------------------
	// Initialize vertically varying, x,y independent basic state
	//----------------------------------------------------------------	
	initialize_vertical_basic_state(NX/2,NY/2);
	
	//----------------------------------------------------------------
	// Initialize topographic array
	//----------------------------------------------------------------
	init_topography();

	//----------------------------------------------------------------
	// Column averaged density
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	
		rhoavg2d[i][j] = 0;

		for(int k=htopo[i][j]+1;k<NZ-1;k++){ rhoavg2d[i][j] = rhoavg2d[i][j] + rhou[k]*tbv[k];}

		rhoavg2d[i][j] = rhoavg2d[i][j]/((double)(NZ-2-htopo[i][j]));
	}}

	//----------------------------------------------------------------
	// Initialize arrays for SOR
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		sor[0][i][j] = 0;
		sor[1][i][j] = 0;
		sor[2][i][j] = 0;
		
	}}

	//----------------------------------------------------------------
	// Remove model base state from environmental base state
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				ITHBAR(i,j,k) = ITHBAR(i,j,k) - tb[k];
				IQBAR(i,j,k) = IQBAR(i,j,k) - qb[k];
			}
		}
	}
#if 0
	//----------------------------------------------------------------
	// Homogenize moisture field
	//----------------------------------------------------------------
	int moisture_lat1 = get_point_from_lat(10);
	int moisture_lat2 = get_point_from_lat(25);
	int moisture_lon1 = get_point_from_lon(75);
	int moisture_lon2 = get_point_from_lon(105);
	
	int moisture_lat0 = get_point_from_lat(22);
	int moisture_lon0 = get_point_from_lon(90);	

	for(int i=moisture_lon1;i<moisture_lon2;i++){
	for(int j=moisture_lat1;j<moisture_lat2;j++){
	for(int k=0;k<NZ;k++){
		
		IQBAR(i,j,k) = IQBAR(moisture_lon0,moisture_lat0,k);
		
	}}}

#endif
	
	//----------------------------------------------------------------
	// Initialize friction array
	//----------------------------------------------------------------
	init_friction();

	//----------------------------------------------------------------
	// Calculate base state vertical velocity
	//----------------------------------------------------------------
	init_basic_state_vertical_velocity();

	//----------------------------------------------------------------
	// Coriolis parameter
	//----------------------------------------------------------------
	for(int j=0;j<NY;j++){ 
		
		f[j] = 2*fc*sin(outLats[j]*trigpi/180.);
		dfdy[j] = 2*fc*cos(outLats[j]*trigpi/180.) * (1./meters_per_degree) * (trigpi/180.) ;
	}
#if 0
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		IUBAR(i,j,k) *= 0.5;
		IVBAR(i,j,k) *= 0.5;
		ITHBAR(i,j,k) *= 0.5;
		IQBAR(i,j,k) *= 0.5;//IQBAR(NX/2,NY/2,k);
		IPBAR(i,j,k) = IPBAR(NX/2,NY/2,k);	

	}}}
#endif	
	//----------------------------------------------------------------
	// Alter humidty, if requested
	//----------------------------------------------------------------
	#if 0 && !ENERGY
	change_humidity();
	#endif
	

}

/*********************************************************************
* 
**********************************************************************/
double convert_z_to_k2(double zf,double z0,double lev,double dzp){
	
	int nt = NZ-1;
	
	double c2 = ( 1 - z0/(dzp*lev)) / ( (nt-1) *dzp - dzp*lev); 
	
	double c1 = z0 / (dzp*lev) - c2*dzp*lev;

	double zi = (-c1 + sqrt(c1*c1+4*c2*zf)) / (2.0*c2);
	
	return zi/dzp + 0.5;
}

/******************************************************************************
*
*
*******************************************************************************/
void initialize_basic_state_idealized(){

	int size = NX*NY*NZ;

	rhoavg = 0;
	mtime = 0;
	bigcounter = 0;
	
	//----------------------------------------------------------------
	// initialize common data
	//----------------------------------------------------------------
	if(PARALLEL)
		memset(&u[0][0][0], 0,NX*NY*NZ*sizeof(double));

	//----------------------------------------------------------------
	// Calculate height of each level
	//----------------------------------------------------------------
	for(int k=0;k<NZ;k++){

		zu[k] = ((double)k-0.5)*dz;
		zw[k] = ((double)k-1)*dz;
	}

	stretched_grid(&zsu[0],&mu[0],&zsw[0],&mw[0],50,1);

	itopo = (double*) calloc(size,sizeof(double));

	//----------------------------------------------------------------
	// initialize data specific to parallel version
	//----------------------------------------------------------------
	#if PARALLEL || ENERGY
	
		ALLOC(iubar,size);
		ALLOC(ivbar,size);
		ALLOC(iwbar,size);
		ALLOC(ithbar,size);
		ALLOC(ipbar,size);
		ALLOC(iqbar,size);
	
		ALLOC(iuistopo,size);
		ALLOC(ivistopo,size);
		ALLOC(iistopo,size);
		ALLOC(ifriction,size);
	
		if(STRETCHED_GRID){initialize_from_era(&zsu[0]);} 
		else {			   initialize_from_era(&zu[0]);}
		
		if(ENERGY){
			
			ALLOC(us,size);
			ALLOC(vs,size);
			ALLOC(ws,size);
			ALLOC(ths,size);
			ALLOC(pis,size);

			//--------------------------------------------------
			// microphysics
			//--------------------------------------------------
			if(USE_MICROPHYSICS){

				ALLOC(qvs,size);
				ALLOC(qcs,size);
				ALLOC(qrs,size);
			}
		}
		
	//----------------------------------------------------------------
	// initialize data specific to serial version
	//----------------------------------------------------------------
	#else
		initialize_flux_cells(NY,NZ);
		initialize_microphysics_cells(NY,NZ);

		initialize_subarray(NX*NY*NZ);

		if(STRETCHED_GRID || IDEAL){ initialize_from_era(&zsu[0]);} 
		else { 				initialize_from_ERA5(&zu[0]);}

		initialize_landsea(landseaMaskFile);

		if(USE_TURBULENT_STRESS){ init_kmix(NX,NY,NZ);}
		if(OUTPUT_DIFFUSION_TEND){ init_damping(NX,NY,NZ);}	
	#endif
		
		
		//----------------------------------------------------------------
		// Initialize stuff for Fourier damping in linearized equations
		//----------------------------------------------------------------
		if(ISLINEAR && FOURIER_DAMPING){
		
			if(!PERIODIC_BOUNDARIES){ init_fftw(NX,NY,NZ);}
			else { init_fftw(NX-6,NY,NZ);}
		}
		

		int base_i = get_point_from_lon(90);
		int base_j = get_point_from_lat(20);

		double full_base_pres[NZ];
		double store_qb[NZ];
		
		for(int k=0;k<NZ;k++){
		
			full_base_pres[k] = IPBAR(base_i,base_j,k);
			//printf("%d %f %f\n",k,zsu[k],full_base_pres[k]);
		}

		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
		for(int k=0;k<NZ;k++){

			IUBAR(i,j,k) = 0;
			IVBAR(i,j,k) = 0;
			ITHBAR(i,j,k) = ITHBAR(base_i,base_j,k);
			IQBAR(i,j,k) = IQBAR(base_i,base_j,k);
			IPBAR(i,j,k) = 0;//IPBAR(NX/2,NY/2,k);	
		}}}

		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){

			ITHBAR(i,j,0) = ITHBAR(i,j,1);
			ITHBAR(i,j,NZ-1) = ITHBAR(i,j,NZ-2);
			IQBAR(i,j,0) = IQBAR(i,j,1);
			IQBAR(i,j,NZ-1) = IQBAR(i,j,NZ-2);
			IUBAR(i,j,0) = IUBAR(i,j,1);
			IUBAR(i,j,NZ-1) = IUBAR(i,j,NZ-2);
			IVBAR(i,j,0) = IVBAR(i,j,1);
			IVBAR(i,j,NZ-1) = IVBAR(i,j,NZ-2);
		}}

	//----------------------------------------------------------------
	// Initialize vertically varying, x,y independent basic state
	//----------------------------------------------------------------	
	initialize_vertical_basic_state(base_i,base_j);
	
	//----------------------------------------------------------------
	// Initialize topographic array
	//----------------------------------------------------------------
	//memset( &IISTOPO(0,0,0),1,NX*NY*NZ*sizeof(double));
	//memset(&IUISTOPO(0,0,0),1,NX*NY*NZ*sizeof(double));
	//memset(&IVISTOPO(0,0,0),1,NX*NY*NZ*sizeof(double));
	//memset(&htopo[0][0], 0,NX*NY*sizeof(int));

	/*

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
		HTOPO(i,j) = 1;
	}}
	*/
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
		//IISTOPO(i,j,k) = 1;
		//IUISTOPO(i,j,k) = 1;
		//IVISTOPO(i,j,k) = 1;
		ITOPO(i,j,k) = 0;
	}}}

	init_topography();

	double *zui;
	
	//----------------------------------------------------------------
	// Which coordinate system?
	//----------------------------------------------------------------
	if(!ISLINEAR)
		zui = &zsu[0];
	else
		zui = &zu[0];

	//----------------------------------------------------------------
	// Coriolis parameter
	//----------------------------------------------------------------
	for(int j=0;j<NY;j++){ 
		
		f[j] = 2*fc*sin(outLats[j]*trigpi/180.);
		dfdy[j] = 2*fc*cos(outLats[j]*trigpi/180.) * (1./meters_per_degree) * (trigpi/180.) ;
	}


	//----------------------------------------------------------------
	//
	// Create lower-level horizontal shear for basic state zonal wind
	//
	//----------------------------------------------------------------
	double ubar_max_lat = 25;
	double ubar_min_lat = 15;
	double ubar_north_lat = 28;
	double ubar_north_zero_lat = 35;
	double ubar_0 = 15;
	double ubar_1 = 8;
	double ubar_max_z = 7000;
	
	int ubar_max_y = get_point_from_lat(ubar_max_lat);
	int ubar_min_y = get_point_from_lat(ubar_min_lat);
	int ubar_north_y = get_point_from_lat(ubar_north_lat);
	int ubar_north_zero_y = get_point_from_lat(ubar_north_zero_lat);
	
	double ubar_north_slope = ubar_1 / (double)(ubar_north_zero_y - ubar_north_y);
	
	double zfrac = 0;
	
	double upper_profile[NY];
	
#if 1
	for(int i=0;i<NX;i++){
	for(int k=0;k<NZ;k++){
		
		if( zui[k] < ubar_max_z ){
		
			for(int j=0;j<ubar_min_y;j++){ IUBAR(i,j,k) = ubar_0; upper_profile[j] = ubar_0;}
		
			for(int j=ubar_min_y;j<ubar_max_y;j++){
			
				IUBAR(i,j,k) = IUBAR(i,j-1,k) - ((ubar_0+ubar_1) / (ubar_max_y - ubar_min_y))*(trigpi/2.0)*sin( trigpi*(j-ubar_min_y) / (ubar_max_y - ubar_min_y) )   ;
				
				upper_profile[j] = upper_profile[j-1] - ((ubar_0+ubar_0) / (ubar_max_y - ubar_min_y))*(trigpi/2.0)*sin( trigpi*(j-ubar_min_y) / (ubar_max_y - ubar_min_y) )   ;
				//if(i==NX/2 && k==10){
					//printf("%d %f\n",j,sin( trigpi*(j-ubar_min_y+1) / (ubar_max_y - ubar_min_y) ));
				//}
			}
		
			for(int j=ubar_max_y;j<ubar_north_y;j++){ IUBAR(i,j,k) = -ubar_1;}
			
			for(int j=ubar_north_y;j<NY;j++){ IUBAR(i,j,k) = -ubar_1 + ubar_north_slope * (double)(j-ubar_north_y) ;}
		}
		//if(i==40 && k==10){ printf("%d %f\n",j,IUBAR(i,j,k));}
	}}
	
	for(int k=1;k<NZ;k++){
		
		zfrac = cos( ((zui[k]-2000) / (ubar_max_z-2000)) * 0.5*trigpi );
		
		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
				
			if( zui[k] < ubar_max_z ){
				
				IUBAR(i,j,k) = IUBAR(i,j,0) * zfrac;
			
				
			} else if(j<base_j && zui[k]<12000){
				
				IUBAR(i,j,k) = upper_profile[j] * zfrac;
				
			} else if(zui[k]>12000 && j<base_j){
				
				IUBAR(i,j,k) = IUBAR(i,j,30);
			}

		}}

	}
	/*
	for(int k=1;k<NZ;k++){
		
		zfrac = cos( ((zui[k]-4000) / (ubar_max_z-0)) * 0.5*trigpi );
		
		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
				
			if( zui[k] < ubar_max_z+4000 && j > base_j ){
				
				IUBAR(i,j,k) = IUBAR(i,j,0) * zfrac;
			
			} 

		}}

	}
	*/
	
	/*
	int i_min = NX/2;
	int i_max = NX/2+NX/4;
	double i_slope = 0.5;
	
	for(int i=i_min;i<i_max;i++){
	
		zfrac = i_slope + i_slope*(double)(i_max-i) / (double)(i_max - i_min);
		//printf("%f\n",zfrac);
		for(int k=0;k<NZ;k++){
		for(int j=0;j<NY;j++){	
		
			//if(zui[k] < ubar_max_z){
		
				IUBAR(i,j,k) *= zfrac;
				//}
		
		}}
	}
	
	for(int i=i_max;i<NX;i++){
	for(int k=0;k<NZ;k++){
	for(int j=0;j<NY;j++){	
		
		//if(zui[k] < ubar_max_z){
		
			IUBAR(i,j,k) = IUBAR(i_max-1,j,k);
			//}
		
	}}}
	
	*/

	
	
	//if(i==40 && j==NY/2-9){ printf("%d %f %f %f\n",k,zfrac,zui[k] / ubar_max_z,IUBAR(i,j,k));}
#endif
	//----------------------------------------------------------------
	//
	// Create upper-level baroclinic shear for basic state zonal wind
	//
	//----------------------------------------------------------------
#if 0

	double ubar_baroclinic_bot_z = 7000;//8000;
	double ubar_baroclinic_top_z = 16000;//18000;
	double ubar_baroclinic_top_ubar = -20;
	
	double slope;
	int klev;
	
	if(!ISLINEAR)
		klev = (int)convert_z_to_k2(ubar_baroclinic_bot_z,50,1,dz);
	else
		klev = (int)(ubar_baroclinic_bot_z / dz + 0.001);
		
	for(int k=0;k<NZ;k++){
		
		if( zui[k] > ubar_baroclinic_bot_z && zui[k] < ubar_baroclinic_top_z ){
			
			for(int i=0;i<NX;i++){
			for(int j=0;j<NY;j++){
				
				slope = (ubar_baroclinic_top_ubar - IUBAR(i,j,klev)) / (ubar_baroclinic_top_z - ubar_baroclinic_bot_z);
				
				IUBAR(i,j,k) = IUBAR(i,j,k-1) + slope * (zui[k] - zui[k-1]) ;
				
				//if(i==40 && j==3){ printf("%d %f %f\n",k,slope,IUBAR(i,j,k));}
			}}
			
		} else if( zui[k] >= ubar_baroclinic_top_z) {
		
			for(int i=0;i<NX;i++){
			for(int j=0;j<NY;j++){
		
				IUBAR(i,j,k) = ubar_baroclinic_top_ubar;
			}}
			
		}

	}
#endif
	
	
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
		
		if(IUBAR(i,j,k) < -5){
			
			//UBAR(i,j,k) = -5;
		}
		
		if(zui[k] > 5000 && j > base_j){
			
			//IUBAR(i,j,k) = -5;
		}
		
	}}}
		

	
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
			
		IUBAR(i,j,0) = IUBAR(i,j,1);
	}}
	
	
	init_basic_state_meridional_velocity(base_j);
	
	
	/*
	//----------------------------------------------------------------
	// Remove model base state from environmental base state
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				//ITHBAR(i,j,k) = ITHBAR(i,j,k) - tb[k];
				//IQBAR(i,j,k) = IQBAR(i,j,k) - qb[k];
			}
		}
	}
	*/
	

	for(int i=0;i<NX-1;i++){	
	for(int k=0;k<NZ;k++){
	
		for(int j=base_j-1;j>=1;j--){
		
			IPBAR(i,j,k) = IPBAR(i,j+1,k) + 0.25*(IUBAR(i,j-1,k)+IUBAR(i+1,j-1,k)+IUBAR(i,j,k)+IUBAR(i+1,j,k)) * dy *f[j];
			//if(i==NX/2 && (k==0 || k==1 || k==2)){
				//printf("%d %d %f %f\n",j,k,IPBAR(i,j,k),IUBAR(i,j,k));
			//}
		}
		//for(int j=2;j<NY;j++){		
		for(int j=base_j+1;j<NY;j++){
		
			IPBAR(i,j,k) = IPBAR(i,j-1,k) - 0.25*(IUBAR(i,j-1,k)+IUBAR(i+1,j-1,k)+IUBAR(i,j,k)+IUBAR(i+1,j,k)) * dy * f[j];
			//if(i==NX/2 && (k==0 || k==1 || k==2)){
				//printf("%d %d %f %f\n",j,k,IPBAR(i,j,k),IUBAR(i,j,k));
			//}
			
			//if(i==NX/2 && j == NY-3){
				//printf("%d %d %f %f %f\n",j,k,zui[k],IPBAR(i,j,k),IUBAR(i,j,k));
			//}
		}
	}}
	
	for(int i=0;i<NX;i++){			
	for(int j=0;j<NY;j++){
	for(int k=1;k<NZ-1;k++){

		ITHBAR(i,j,k) = tb[k]*(IPBAR(i,j,k-1) - IPBAR(i,j,k)) / (- grav * (zui[k]-zui[k-1]) );
		
	}}}
	
	
			
	//double store_temp[NZ];
	double store_pres[NZ];


	for(int k=0;k<NZ;k++){

		//store_temp[k] = ITHBAR(base_i,base_j,k);
		store_pres[k] =  IPBAR(base_i,base_j,k);

	}

	
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		IPBAR(i,j,k) =  (IPBAR(i,j,k) - store_pres[k]) / (cp*tbv[k]) + full_base_pres[k];
		//ITHBAR(i,j,k) -= store_temp[k];

	}}}
	
	
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
			
		for(int k=NZ-2;k>0;k--){
			
			if( ITHBAR(i,j,k) + tb[k] > ITHBAR(i,j,k+1) + tb[k+1] ){
				//printf("%d %d %d %f %f\n",i,j,k,ITHBAR(i,j,k) + tb[k],ITHBAR(i,j,k+1) + tb[k+1]);
				ITHBAR(i,j,k) = ITHBAR(i,j,k+1) + tb[k+1] - tb[k] - 0.1;
				//printf("%d %d %d %f %f\n",i,j,k,ITHBAR(i,j,k) + tb[k],ITHBAR(i,j,k+1) + tb[k+1]);
			}
			
		}
			
	}}
	
	int south_thbar = get_point_from_lat(8);
	
	for(int i=0;i<NX;i++){
	for(int j=0;j<south_thbar;j++){			
	for(int k=0;k<NZ;k++){
		
		ITHBAR(i,j,k) = ITHBAR(i,south_thbar,k);
		
	}}}
			
	
	
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			
			ITHBAR(i,j,2) = ITHBAR(i,j,3);
			ITHBAR(i,j,1) = ITHBAR(i,j,2);
			ITHBAR(i,j,0) = ITHBAR(i,j,1);
	}}
	
	for(int i=0;i<NX;i++){
	for(int k=0;k<NZ;k++){
			
		ITHBAR(i,0,k) = ITHBAR(i,1,k);
	}}

	
	double smr,temp,pres,rh;
	double base_temp[NZ];
	double base_pres[NZ];
	double base_relh[NZ];
	
	//-----------------------------------------------------------------------
	// Relative humidity changes
	//-----------------------------------------------------------------------
	int rh_y_min = get_point_from_lat(5);
	int rh_y_max = get_point_from_lat(14);
	int rh_y_min2 = get_point_from_lat(30);
	int rh_y_max2 = get_point_from_lat(36);
	
	double rh_decrease = 0.2;
	double rh_decrease2 = 0.3;
	double rh_drop[NY];
	double rh_slope = rh_decrease / (double)(rh_y_max-rh_y_min);
	double rh_slope2 = rh_decrease2 / (double)(rh_y_min2-rh_y_max2);
	
	for(int j=NY;j>rh_y_max;j--){			rh_drop[j] = 0;							}
	
	for(int j=rh_y_max;j>rh_y_min;j--){		rh_drop[j] = rh_drop[j+1] + rh_slope;	}
	
	for(int j=rh_y_min;j>=0;j--){			rh_drop[j] = rh_decrease;				}

	for(int j=rh_y_min2;j<rh_y_max2;j++){	rh_drop[j] = rh_drop[j-1] - rh_slope2;	}
	
	for(int j=rh_y_max2;j<NY;j++){			rh_drop[j] = rh_decrease2;				}

	//-----------------------------------------------------------------------
	// Calculate mixing ratio
	//-----------------------------------------------------------------------	
	for(int k=0;k<NZ;k++){
		
		base_temp[k] = tb[k]*pib[k];
		base_pres[k] = p0*pow(pib[k],(cp/Rd));
		
		base_relh[k] = qb[k] / get_QV_Sat(base_temp[k],base_pres[k]);
	}


	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=1;k<NZ-1;k++){

		temp = (ITHBAR(i,j,k)+tb[k]) * IPBAR(i,j,k);	// full, actual temperature
		pres = p0*pow(IPBAR(i,j,k),(cp/Rd));			// full, dimensional pressure

		smr = get_QV_Sat(temp,pres);					// calculate saturation mixing ratio

		IQBAR(i,j,k) = (base_relh[k]-rh_drop[j]) * smr - qb[k];

		if(i==NX/2 && k==NZ-1){ printf("%d %d %f %f %f %f %f %f\n",j,k,outLats[j],pres/100.0,base_relh[k],qb[k]*1000,smr*1000,IQBAR(i,j,k)*1000);}		
		//if(i==NX/2 && j==NY/4){ printf("%d %f %f %f %f %f %f\n",k,zui[k],pres/100.0,base_relh[k],qb[k]*1000,smr*1000,IQBAR(i,j,k)*1000);}
	
	}}}

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		IQBAR(i,j,0) = IQBAR(i,j,1);
		IQBAR(i,j,NZ-1) = IQBAR(i,j,NZ-2);
	}}

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
		
		IUBAR(i,j,k) -= 5.0;
		
	}}}



	//----------------------------------------------------------------
	// Initialize friction array
	//----------------------------------------------------------------
	init_friction();

	//----------------------------------------------------------------
	// Calculate base state vertical velocity
	//----------------------------------------------------------------
	init_basic_state_vertical_velocity();




}


/*********************************************************************
* Initialize the perturbation fields
*
**********************************************************************/
void initialize_perturbation(){

	if(RESTART){
		
		if(PARALLEL){	
			initialize_from_output(restartFileName,restartFileTime);
		} else {
			initialize_from_output_serial(restartFileName,restartFileTime);			
		}
	}

	if(!RESTART){
		
		//----------------------------------------------------------------
		// Actions to perform on only the root process
		//----------------------------------------------------------------
		if(PARALLEL && rank==0){
		
			//----------------------------------------------------------------
			// Set the initial perturbation values to zero
			//----------------------------------------------------------------
			memset( &IUBAR(0,0,0),0,NX*NY*NZ*sizeof(double));
			memset( &IVBAR(0,0,0),0,NX*NY*NZ*sizeof(double));
			memset(&ITHBAR(0,0,0),0,NX*NY*NZ*sizeof(double));
			memset( &IQBAR(0,0,0),0,NX*NY*NZ*sizeof(double));
			memset( &IPBAR(0,0,0),0,NX*NY*NZ*sizeof(double));
		
			//----------------------------------------------------------------
			// Do the initialization
			//----------------------------------------------------------------
			int lat = get_point_from_lat(18);
			int lon = get_point_from_lon(89);
		
			printf("latpoint = %d\n",lat);

			//initialize_vortex(1100000,3000, 11000, -2.05, 1.5,NX/2,lat);		
			initialize_vortex(1500000,3000, 11000, -4.5, 4.05,NX/3,lat);
			//initialize_vortex(1000000,6000, 16000, 0.5, -0.5,lon,lat);
			//initialize_vortex(1500000,6000, 16000, -2.0, 2.0,NX/2-500.0/15.0,get_point_from_lat(18));
		}
	
		//----------------------------------------------------------------
		// Distribute the perturbation arrays from the root process to
		// the other processes
		//----------------------------------------------------------------
		if(PARALLEL){
	
			distributeArray(ths,&ITHBAR(0,0,0));	distributeArray(thms,&ITHBAR(0,0,0));
			distributeArray(us,  &IUBAR(0,0,0));	distributeArray(ums,  &IUBAR(0,0,0));
			distributeArray(vs,  &IVBAR(0,0,0));	distributeArray(vms,  &IVBAR(0,0,0));

			exchange(ths); exchange(thms);
			exchange(us); exchange(ums);
			exchange(vs); exchange(vms);
						
			if(USE_MICROPHYSICS){
				
				distributeArray(qvs,&IQBAR(0,0,0));	distributeArray(qvms,&IQBAR(0,0,0));
				exchange(qvs); exchange(qvms);
			}
		
		}	
	}
}

/***********************************************************************
* 
* INITIALIZE VERTICALLY VARYING BASIC STATE
*
************************************************************************/
void initialize_vertical_basic_state(int ibase,int jbase){
	
	//-----------------------------------------------------
	// create base state potential temperature profile
	//-----------------------------------------------------
	tb[0] = ITHBAR(ibase,jbase,0);
	tbw[0] = 0.5*(ITHBAR(ibase,jbase,1)+ITHBAR(ibase,jbase,0));

	for(int k=1;k<NZ;k++){ tb[k] = ITHBAR(ibase,jbase,k);}

	for(int k=1;k<NZ;k++){ tbw[k] = 0.5*(tb[k]+tb[k-1]);}

	//-----------------------------------------------------
	// create base state specific humidity profile
	//-----------------------------------------------------
	qb[0] = IQBAR(ibase,jbase,0);

	for(int k=1;k<NZ;k++){ qb[k] = IQBAR(ibase,jbase,k);}

	//-----------------------------------------------------
	// create base state virtual potential temperature profile
	//-----------------------------------------------------
	for(int k=0;k<NZ;k++){ tbv[k]=tb[k]*(1.0+0.61*qb[k]);}

	//-----------------------------------------------------
	// create base state pressure profile
	//-----------------------------------------------------
	double pisfc = pow((pressfc/p0),( Rd/cp));
	
	if(STRETCHED_GRID){
		
		pib[1] = pisfc - grav * (zsu[1]-0) / (cp * tbv[1]);
		pib[0] = pisfc + grav * (0-zsu[0]) / (cp * tbv[0]);

		for(int k=2;k<NZ;k++){ pib[k] = pib[k-1]-grav*(zsu[k]-zsu[k-1])/(cp*(  0.5*(tbv[k]+tbv[k-1]) ));}
		
	} else {
		
		pib[1] = pisfc - grav * 0.5 * dz / (cp * tbv[1]);
		pib[0] = pisfc + grav * 0.5 * dz / (cp * tbv[0]);

		for(int k=2;k<NZ;k++){ pib[k] = pib[k-1]-grav*dz/(cp*(  0.5*(tbv[k]+tbv[k-1]) ));}
	}

	//-----------------------------------------------------
	// create base state density profile
	//-----------------------------------------------------
	for(int k=0;k<NZ;k++){ rhou[k] = p0*pow(pib[k],(cv/Rd)) / (Rd*tbv[k]);}

	rhow[1] = pressfc / (Rd*tbw[1]);

	for(int k=1;k<NZ;k++){ rhow[k] = 0.5*(rhou[k]+rhou[k-1]);}

	for(int k=1;k<NZ;k++){ rhoavg = rhoavg + rhou[k];}

	rhoavg = rhoavg/((double)NZ-2);

	//-----------------------------------------------------
	// set values for non-physical points
	//-----------------------------------------------------
	tb[0]=tb[1];
	tb[NZ-1]=tb[NZ-2];
	tbv[0]=tbv[1];
	tbv[NZ-1]=tbv[NZ-2];
	tbw[0]=tbw[1];
	tbw[NZ-1]=tbw[NZ-2];
	pib[0]=pib[1];
	pib[NZ-1]=pib[NZ-2];
//	rhou[0]=rhou[1];
//	rhou[NZ-1]=rhou[NZ-2];
	rhow[0]=rhow[1];
//	rhow[NZ-1]=rhow[NZ-2];

	for(int k=0;k<NZ;k++){

		one_d_rhow[k] = 1. / rhow[k];
		one_d_rhou[k] = 1. / rhou[k];
	}

	if(VERBOSE){
		if(STRETCHED_GRID){
			for(int k=0;k<NZ;k++){
				printf("%d\tzu = %.4f\tzw = %.4f\tpib = %.5f\ttb = %.1f\ttbv = %.1f\trhou = %.3f\trhow = %.3f\tqb = %.3f\n",
						k,zsu[k]/1000,zsw[k]/1000,pib[k],tb[k],tbv[k],rhou[k],rhow[k],qb[k]*1000);
			}
		} else {
			for(int k=0;k<NZ;k++){
				printf("%d\tzu = %.2f\tzw = %.2f\tpib = %.3f\ttb = %.1f\ttbv = %.1f\trhou = %.3f\trhow = %.3f\tqb = %.3f\n",
						k,zu[k]/1000,zw[k]/1000,pib[k],tb[k],tbv[k],rhou[k],rhow[k],qb[k]*1000);			
			}
		}
	} 
	
}

/*********************************************************************
*
*
**********************************************************************/
void zero_perturbation_fields(){
	
	memset(&U(0,0,0), 0,NX*NY*NZ*sizeof(double));
	memset(&UP(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&UM(0,0,0),0,NX*NY*NZ*sizeof(double));	
	memset(&V(0,0,0), 0,NX*NY*NZ*sizeof(double));
	memset(&VP(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&VM(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&W(0,0,0), 0,NX*NY*NZ*sizeof(double));
	memset(&WP(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&WM(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&TH(0,0,0), 0,NX*NY*NZ*sizeof(double));
	memset(&THP(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&THM(0,0,0),0,NX*NY*NZ*sizeof(double));
	memset(&PI(0,0,0),0,NX*NY*NZ*sizeof(double));
}

/*********************************************************************
*
*
**********************************************************************/
void const_times_array3D(double *var,double a,int size){
	
	for(int i=0;i<size;i++){ var[i] *= a;}
}

/*********************************************************************
*
*
**********************************************************************/
void reinitialize_perturbation(int time,double scale){
	
	mtime = 0;
	bigcounter = 0;
	
	zero_perturbation_fields();

	if(USE_MICROPHYSICS){ init_microphysics(); }

	initialize_from_output_serial(restartFileName,time);
	
	int size = NX*NY*NZ;
	
	const_times_array3D(&U(0,0,0),scale,size);
	const_times_array3D(&V(0,0,0),scale,size);
	const_times_array3D(&W(0,0,0),scale,size);
	const_times_array3D(&TH(0,0,0),scale,size);
	
	const_times_array3D(&UM(0,0,0),scale,size);
	const_times_array3D(&VM(0,0,0),scale,size);
	const_times_array3D(&WM(0,0,0),scale,size);
	const_times_array3D(&THM(0,0,0),scale,size);	
}

/*********************************************************************
*
*
**********************************************************************/
void reinitialize(){

	mtime = 0;
	bigcounter = 0;

	//----------------------------------------------------------------
	// Initialize arrays for SOR
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		sor[0][i][j] = 0;
		sor[1][i][j] = 0;
		sor[2][i][j] = 0;
		
	}}

	//----------------------------------------------------------------
	// Reset perturbation arrays
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		U(i,j,k) = 0;
		V(i,j,k) = 0;
		TH(i,j,k) = 0;
		W(i,j,k) = 0;
		PI(i,j,k) = 0;

		UP(i,j,k) = U(i,j,k);
		UM(i,j,k) = U(i,j,k);

		VP(i,j,k) = V(i,j,k);
		VM(i,j,k) = V(i,j,k);
	#if !HYDROSTATIC
		WP(i,j,k) = W(i,j,k);
		WM(i,j,k) = W(i,j,k);
	#endif
		THP(i,j,k) = TH(i,j,k);
		THM(i,j,k) = TH(i,j,k);
	}}}
}

/********************************************************
* 
*********************************************************/
void init_basic_state_meridional_velocity(int middle_j){
	
	for(int i=0;i<NX-1;i++){
	for(int k=0;k<NZ;k++){
		
		for(int j=middle_j;j<NY-1;j++){
		
			IVBAR(i,j+1,k) = IUBAR(i,j,k) - IUBAR(i+1,j,k) + IVBAR(i,j,k);
		
		}
		
		for(int j=middle_j-1;j>=0;j--){
		
			IVBAR(i,j,k) = IUBAR(i+1,j,k) - IUBAR(i,j,k) + IVBAR(i,j+1,k);
		
		}
	}}
}

/********************************************************
* 
*********************************************************/
void init_basic_state_vertical_velocity(){
	
	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){

		IWBAR(i,j,htopo[i][j]) = 0;		// boundary condition = zero vertical 
		IWBAR(i,j,htopo[i][j]+1) = 0;	// velocity at bottom

		for(int k=htopo[i][j]+2;k<NZ-1;k++){
			
			IWBAR(i,j,k) = (rhow[k-1]/rhow[k])*IWBAR(i,j,k-1) 
				+ (rhou[k-1]/rhow[k])*
				(
					-(IUBAR(i+1,j,k-1) - IUBAR(i,j,k-1))/dx  	
		      		-(IVBAR(i,j+1,k-1) - IVBAR(i,j,k-1))/dy
				)*DZU(k);
			
		}

		IWBAR(i,j,NZ-1) = IWBAR(i,j,NZ-2);
	}}
}

/********************************************************
* Initialize friction
*********************************************************/
void init_friction(){

	for(int i=0;i<NX;i++)
		for(int j=0;j<NY;j++)
			for(int k=0;k<NZ;k++)
				IFRICTION(i,j,k) = 0;

	if(USE_LINEAR_FRICTION){

		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){

			if(ISLINEAR){ IFRICTION(i,j,htopo[i][j]+1) = 2.0e-5;}
			else { IFRICTION(i,j,htopo[i][j]+1) = 2.0e-5;}
		}}
	}
	
	if(USE_LINEAR_FRICTION && !MERIDIONAL_CROSS_SECTION && !IDEAL){

		for(int i=1;i<NX-1;i++){
		for(int j=1;j<NY-1;j++){
		for(int k=HTOPO(i,j)+1;k<NZ;k++){
			
			if(
				(IISTOPO(i,j+1,k)==0 ||
				IISTOPO(i,j-1,k)==0 ||
				IISTOPO(i+1,j,k)==0 ||
				IISTOPO(i-1,j,k)==0) &&
				IFRICTION(i,j,k) == 0
				
			){
				//printf("%d %d %d %f %f\n",i,j,k,outLats[j],outLons[i]);
				IFRICTION(i,j,k) = 2.0e-5;
			}
		}}}
	}
}

/********************************************************
* Initialize topography
*********************************************************/
void init_topography(){

	int height = 0;

	memset(&htopo[0][0], 0,NX*NY*sizeof(int));

	//for(int i=NX-1;i>=0;i--){
	for(int i=NX-2;i>=0;i--){
	for(int j=NY-2;j>=0;j--){

		height = 0;

		for(int k=0;k<NZ;k++){

		#if !STRETCHED_GRID
			#if !MERIDIONAL_CROSS_SECTION
				if(zu[k] < ITOPO(i,j,k) && ITOPO(i,j,k) > 0)
			#else
				if(zu[k] < ITOPO(linear_lon,j,k) && ITOPO(linear_lon,j,k) > 70000)
			#endif
		#else
			#if !MERIDIONAL_CROSS_SECTION
				if(zsu[k] < ITOPO(i,j,k) && ITOPO(i,j,k) > 0)
			#else
				//if(i < linear_lon && (zsu[k] < ITOPO(linear_lon,j,k) && ITOPO(linear_lon,j,k) > 0) )	
				if( zsu[k] < ITOPO(linear_lon,j,k) && ITOPO(linear_lon,j,k) > 0 )
			#endif
		#endif	
			{
				IISTOPO(i,j,k) = 0;
				IUISTOPO(i,j,k) = 0;
				IVISTOPO(i,j,k) = 0;
				IUISTOPO(i+1,j,k) = 0;
				IVISTOPO(i,j+1,k) = 0;
				height = k;
#if 0
			} else if(i >= linear_lon && (zsu[k] < ITOPO(i,j,k) && ITOPO(i,j,k) > 0) ){
				
				IISTOPO(i,j,k) = 0;
				IUISTOPO(i,j,k) = 0;
				IVISTOPO(i,j,k) = 0;
				IUISTOPO(i+1,j,k) = 0;
				IVISTOPO(i,j+1,k) = 0;
				height = k;
#endif
			} else {

				IISTOPO(i,j,k) = 1;
				IUISTOPO(i,j,k) = 1;
				IVISTOPO(i,j,k) = 1;
			}
		}

		htopo[i][j] = height;
		//printf("%d %d %d\n",i,j,height);
	}}

}

/*********************************************************************
* Alter basic state humidity
**********************************************************************/
void change_humidity(){

	double smr,temp,pres,rh;
	
	for(int k=0;k<NZ;k++){
		
		temp = tb[k]*pib[k];
		pres = p0*pow(pib[k],(cp/Rd));
		
		smr = get_QV_Sat(temp,pres);
		
		//printf("%d %f %f %f\n",k,smr*1000,qb[k]*1000,qb[k]/smr);
		
	}

	//initialize_landsea(landseaMaskFile);

	const double height_limit_low = 2000.;
	const double height_limit_high = 10000.;
	const double relative_humidity = 0.90;

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=1;k<NZ;k++){
		//-----------------------------------------------------------------
		// Limit region
		//-----------------------------------------------------------------
		if(outLats[j] > 18 && outLats[j] < 23 && outLons[i] > 82 && outLons[i] < 87){
			
			temp = (ITHBAR(i,j,k)+tb[k]) * IPBAR(i,j,k);	// full, actual temperature
			pres = p0*pow(IPBAR(i,j,k),(cp/Rd));			// full, dimensional pressure

			smr = get_QV_Sat(temp,pres);					// calculate saturation mixing ratio

			rh = (IQBAR(i,j,k)+qb[k])/smr;					// calculate relative humidity
		
			// change humidity
			if(zsu[k] >= height_limit_low && zsu[k] <=  height_limit_high){ 
			
				//IQBAR(i,j,k) = relative_humidity * smr - qb[k];
		
			}

		}
	}}}

	//free(landsea);
}

/*********************************************************************
* Calculate heights for stretched grid.
*
**********************************************************************/
void stretched_grid(double * zu,double *mu,double * zw,double *mw,double z0,int lev){
	
	double c1;
	double c2;
	
	int nt = NZ-1;
	double dzp = dz;
	
	c2 = ( 1 - z0/(dzp*lev)) / ( (nt-1) *dzp - dzp*lev); 
	
	c1 = z0 / (dzp*lev) - c2*dzp*lev;

	//printf("%f %f\n",c1,c2);
	
	double zh1,zh2;
	double m1,m2;
	double fullLev;
	double halfLev;
	
	for(int i=0;i<=nt;i++){
		
		fullLev = ((double)i-1.0)*dzp;
		halfLev = ((double)i-0.5)*dzp;
		
		zh1 = (c1+c2*fullLev)*fullLev;
		m1 = 1 / (c1+2*c2*fullLev);
	
		zh2 = (c1+c2*halfLev)*halfLev;
		m2 = 1 / (c1+2*c2*halfLev);
		
		//printf("%d %.0f %f %f %f\n",i,zh1,m1,zh2,m2);
		
		zu[i] = zh2;
		zw[i] = zh1;
		
		mu[i] = m2;
		mw[i] = m1;
	}
	
}

/*********************************************************************
* 
**********************************************************************/
int rindex(int r,int k,int NR){ return r+k*NR;}

/*********************************************************************
* 
**********************************************************************/
void rinterp(int xc,int yc,int x,int y,int *ind, double *frac, double offset){
	
	double dist = sqrt( (double)(x-xc)*(x-xc) + (double)(y-yc)*(y-yc) ) + offset;
	
	//printf("d = %f",dist);
	
	*ind = (int)dist;
	
	*frac = dist - (double)((int)dist);
}

/*********************************************************************
*
*
**********************************************************************/
void initialize_vortex(double r_end, double zm, double zt, double tmax, double tmin, int xpos, int ypos){
	
	int NR = (int)(r_end / dx);
	
	double *t_z = (double*) calloc(NZ,sizeof(double));
	double *t_r = (double*) calloc(NR,sizeof(double));
	double *phi = (double*) calloc(NR*NZ,sizeof(double));
	double *vr  = (double*) calloc(NR*NZ,sizeof(double));
	double *rh  = (double*) calloc(NR*NZ,sizeof(double));

	double zmid = 0.5*(zm+zt);
	double R,tmid;
	
	//-------------------------------------------------------------------------------------
	// Calculate temperature field
	//-------------------------------------------------------------------------------------	
	for(int k=0;k<NZ;k++){
		
		if(zsu[k]<=zm){ t_z[k] = tmin * exp( -(zsu[k]*zsu[k])/(zm*zm/8.0) );		
		} else { t_z[k] = tmax * exp( -( (zsu[k]-zmid)*(zsu[k]-zmid) ) / ( (zt-zm)*(zt-zm) / 32.0) );}
	}

	for(int r=0;r<NR;r++){ t_r[r] = exp( -dx*r*dx*r / (r_end*r_end/8.0) );}

	//-------------------------------------------------------------------------------------
	// Calculate moisture field
	//-------------------------------------------------------------------------------------
	double zt_rh = 11000;
	double zb_rh = -2000;
	double rh_prime = 0.20;
	double r_end_rh = 1000000;//900000;
	double rh_max = 0.95;
	
	zmid = 0.5*(zt_rh+zb_rh);
	
	//-------------------------------------------------------------------------------------
	// Anomalous relative humidity
	//-------------------------------------------------------------------------------------
	for(int k=0;k<NZ;k++){
	for(int r=0;r<NR;r++){
				
		rh[rindex(r,k,NR)] =
			rh_prime * exp( -( (zsu[k]-zmid)*(zsu[k]-zmid) ) / ( (zt_rh-zb_rh)*(zt_rh-zb_rh) / 32.0) )
			* exp( -dx*r*dx*r / (r_end_rh*r_end_rh/8.0) );
	}}
	//-------------------------------------------------------------------------------------
	// Convert to mixing ratio anomaly
	//-------------------------------------------------------------------------------------
	double temp,pres,smr,rhp;

	for(int k=0;k<NZ;k++){
		
		temp = tb[k]*pib[k];
		pres = p0*pow(pib[k],(cp/Rd));
		smr = get_QV_Sat(temp,pres);
		
		for(int r=0;r<NR;r++){ 
			
			if( rh[rindex(r,k,NR)] + qb[k]/smr < rh_max ){ rh[rindex(r,k,NR)] *= smr;} 
			else { rh[rindex(r,k,NR)] = (rh_max - qb[k]/smr) * smr;}
		
			//if(r==0){ printf("%f rh = %f %f %f %f\n",zsu[k],rh[rindex(r,k,NR)]*1000.0,qb[k]*1000.0,smr*1000.0,100*(rh[rindex(r,k,NR)]+qb[k])/smr);}
		}
	}
	//-------------------------------------------------------------------------------------
	// Calculate pressure field
	//-------------------------------------------------------------------------------------
	for(int k=NZ-2;k>0;k--){
	for(int r=0;r<NR;r++){
	
		tmid = 0.5* (t_z[k]*t_r[r]+t_z[k+1]*t_r[r]) / ((tb[k]+tb[k+1])*0.5) + 0.61*rh[rindex(r,k,NR)];
	
		phi[rindex(r,k,NR)] = phi[rindex(r,k+1,NR)] - grav * tmid * (zsu[k+1]-zsu[k]);
	

		//if(r==0){
			//printf("%d %f %f\n",k,t_z[k],phi[rindex(r,k,NR)]);
		//}

	}}
	
	//-------------------------------------------------------------------------------------
	// Caluculate gradient wind balance
	//-------------------------------------------------------------------------------------
	double fcor = 2*7.292e-5*sin(outLats[ypos]*trigpi/180.0);
	double radicand;
	
	for(int k=0;k<NZ;k++){
		
		vr[rindex(0,k,NR)] = 0;
		
		for(int r=1;r<NR-1;r++){
		
			R = r*dx;
	
			if(phi[rindex(r,k,NR)]<0){
	
				vr[rindex(r,k,NR)] = -0.5*fcor*R + 
					0.5*sqrt( fcor*fcor*R*R + 4*R * 0.5*
						fabs( phi[rindex(r+1,k,NR)]-phi[rindex(r-1,k,NR)] )
							 * one_d_dx );
			} else {
				
				radicand = fcor*fcor*R*R - 4*R * 0.5* fabs( phi[rindex(r+1,k,NR)]-phi[rindex(r-1,k,NR)] ) * one_d_dx ;
				
				if(radicand<0){ radicand = 0;}
				
				vr[rindex(r,k,NR)] = 0.5*fcor*R - 0.5*sqrt(radicand) ;
				
				vr[rindex(r,k,NR)] *= -1.0;
			}
		}
	}
	
	double iproj,jproj,windmag,frac;
	int ind;
	//-------------------------------------------------------------------
	// Loop over box containing the vortex
	//-------------------------------------------------------------------
	for(int i=xpos-NR;i<xpos+NR;i++){
	for(int j=ypos-NR;j<ypos+NR;j++){
	for(int k=1;k<NZ;k++){
		
		if(i>0 && j>0 && i<NX && j<NY){
			//---------------------------------------------------------------
			// Interpolate Cartesian to radial coordinates
			//---------------------------------------------------------------	
			rinterp(xpos,ypos,i,j,&ind,&frac,0);
			//---------------------------------------------------------------
			// Don't exceed array bounds
			//---------------------------------------------------------------		
			if(ind < NR-2){
				//-----------------------------------------------------------
				// Potential temperature
				//-----------------------------------------------------------
				ITHBAR(i,j,k) += t_z[k] * ( (1.0-frac)*t_r[ind] + frac*t_r[ind+1] );
				//-----------------------------------------------------------
				// Water vapor
				//-----------------------------------------------------------
				IQBAR(i+(int)(5*meters_per_degree/dx),j-(int)(0*meters_per_degree/dx),k) 
							= ( (1.0-frac)*rh[rindex(ind,k,NR)] + frac*rh[rindex(ind+1,k,NR)] );	// POSSIBLE ARRAY OUT OF BOUNDS!!!
				
				//IQBAR(i+(int)(0.0*meters_per_degree/dx),j-(int)(0.0*meters_per_degree/dx),k) 
				//			= ( (1.0-frac)*rh[rindex(ind,k,NR)] + frac*rh[rindex(ind+1,k,NR)] );
			
				//-----------------------------------------------------------
				// Wind components (non-staggered)
				//-----------------------------------------------------------					
				if(ind != 0){

					iproj = fabs( (double)(j-ypos)) / ((double)ind + frac);
					jproj = fabs( (double)(i-xpos)) / ((double)ind + frac);
				
					windmag = (1.0-frac)*vr[rindex(ind,k,NR)] + frac*vr[rindex(ind+1,k,NR)];

					if(j<ypos){ IUBAR(i,j,k) += +windmag * iproj; }
					else { 		IUBAR(i,j,k) += -windmag * iproj; }
		
					if(i<xpos){ IVBAR(i,j,k) += -windmag * jproj; }
					else { 		IVBAR(i,j,k) += +windmag * jproj; }
				}
			}
		}
	}}}
			
	//---------------------------------------------------------------
	// Stagger wind components
	//---------------------------------------------------------------	
	for(int i=1;i<NX;i++){
	for(int j=1;j<NY;j++){
	for(int k=0;k<NZ;k++){
	
		IWBAR(i,j,k) = 0.5*(IUBAR(i,j,k) + IUBAR(i-1,j,k));
		IPBAR(i,j,k) = 0.5*(IVBAR(i,j-1,k) + IVBAR(i,j,k));
	}}}
	
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
	
		IUBAR(i,j,k) = IWBAR(i,j,k);
		IVBAR(i,j,k) = IPBAR(i,j,k);
	}}}	
		
	free(t_z); free(t_r); free(vr); free(phi);
}



/*********************************************************************
* Calculate saturation mixing ratio.
*
* temperature - full, actual temperature field (K, not potential, not perturation)
* pressure - full, dimensional pressure (Pa)
**********************************************************************/
double get_QV_Sat(double temperature,double pressure){

	double esl = 611.2 * exp(17.67 * (temperature-273.15) / (temperature - 29.65) );

	return 0.62197 * esl / (pressure-esl);
}

/*********************************************************************
* Load data for perturbation fields from file for use as an initial condition
*
* filename 	- file containing input fields for perturbation
* varname 	- name of variable within file
* var,mvar	- fully interpolated output arrays
* time 		- time at which to initialize
**********************************************************************/
void load_from_output(const char *filename,const char *varname,double *var,double *mvar,size_t time){
	
	if(rank==0){ get_data_at_time(filename,varname,time,iubar);}
	
	distributeArray(var);
	
	exchange(var);
	
	for(int i=0;i<fNX*fNY*fNZ;i++){ mvar[i] = var[i];}
}

/*********************************************************************
* Load data for perturbation fields from file for use as an initial condition 
* and interpolate it to the model grid.
*
* filename 	- file containing input fields for perturbation
* varname 	- name of variable within file
* time 		- time at which to initialize
* varIn		- allocated memory for input field at original resolution (xdim x ydim x zdim array)
* var_interpz - allocated memory for input field interpolated to given heights height (xdim x ydim x NZ array)
* var,mvar	- fully interpolated output arrays
* levs_out	- output heights
* myLons		- input longitude
* myLats		- input latitude
* zlevs		- input heights
* xdim		- input x-dimension
* ydim		- input y-dimension
* zdim		- input z-dimension
* myDX		- input dx grid spacing
* myDY		- input dy grid spacing
**********************************************************************/
void interpolate_from_output(
			const char *filename,const char *varname,
			size_t time,
			double *varIn,double *var_interpz,double *varOut,
			double *levs_out,
			double *myLons,double *myLats,double *zlevs,
			int xdim,int ydim,int zdim,
			double myDX,double myDY
	){
		
		for(int i=0;i<xdim*ydim*zdim;i++){ varIn[i] = 0;}
		for(int i=0;i<xdim*ydim*NZ;i++){ var_interpz[i] = 0;}
	
		get_data_at_time(filename,varname,time,varIn,xdim,ydim,zdim);
		
		vert_interpolate_1d_from_model(levs_out,zlevs, varIn, xdim,ydim,zdim,NZ,var_interpz);
		
		horz_interpolate_from_model(var_interpz,varOut,xdim,ydim,NZ,false,false,
		myDX/meters_per_degree,myDY/meters_per_degree,outLons[0]-myLons[0],outLats[0]-myLats[0]);
}

/*********************************************************************
* Load data for perturbation fields from file for use as an initial condition 
* and interpolate it to the model grid.
*
* filename 	- file containing input fields for perturbation
* varname 	- name of variable within file
* time 		- time at which to initialize
* varIn		- allocated memory for input field at original resolution (xdim x ydim x zdim array)
* var_interpz - allocated memory for input field interpolated to given heights height (xdim x ydim x NZ array)
* var,mvar	- fully interpolated output arrays
* levs_out	- output heights
* myLons		- input longitude
* myLats		- input latitude
* zlevs		- input heights
* xdim		- input x-dimension
* ydim		- input y-dimension
* zdim		- input z-dimension
* myDX		- input dx grid spacing
* myDY		- input dy grid spacing
**********************************************************************/
void load_interpolate_from_output(
			const char *filename,const char *varname,
			size_t time,
			double *varIn,double *var_interpz,
			double *var,double *mvar,
			double *levs_out,
			double *myLons,double *myLats,double *zlevs,
			int xdim,int ydim,int zdim,
			double myDX,double myDY
	){
	//-----------------------------------------------------------------
	// Read file and interpolate on root process
	//-----------------------------------------------------------------
	if(rank==0){
		memset(iubar,0,NX*NY*NZ*sizeof(double));
		interpolate_from_output(filename,varname,time,varIn,var_interpz,iubar,zsu,myLons,myLats,zlevs,xdim,ydim,zdim,myDX,myDY);
	}
	
	//-----------------------------------------------------------------
	// Distribute results to other processes
	//-----------------------------------------------------------------
	distributeArray(var);
	
	exchange(var);
	
	for(int i=0;i<fNX*fNY*fNZ;i++){ 
		//var[i] *= 0.5;
		mvar[i] = var[i];
	}
}

/*********************************************************************
* Initialize model perturbation fields from input file. Will determine
* whether or not interpolation is required.
*
* myfilename - file containing input fields for perturbations
* time		 - time at which to initialize
**********************************************************************/
void initialize_from_output(const char *myfilename,size_t time){
	
	size_t dims[3],xdim,ydim,zdim;
	double interp_dx,interp_dy,interp_dz;
	float grid_spacing[3];
	int same = 0;
	
	//------------------------------------------------------------------------------
	// Test files for same dimensions and grid spacing
	//------------------------------------------------------------------------------
	if(rank==0){
	
		get_dims(myfilename,"x","y","z",dims);
		get_grid_spacing(myfilename,"dx","dy","dz",grid_spacing);
	
		interp_dx = grid_spacing[0];
		interp_dy = grid_spacing[1];
		interp_dz = grid_spacing[2];
		
		xdim = dims[0];
		ydim = dims[1];
		zdim = dims[2];
	
		if(xdim == NX && ydim == NY && zdim == NZ && interp_dx == dx && interp_dy == dy && interp_dz == dz){
			same = 1;
		} else { 
			same = 0;
		}
	}

	broadcast_data_int(1,&same);	// tell other processors whether the dimensions/grid spacing are the same

	//------------------------------------------------------------------------------
	// If the files have the same dimensions and grid spacing, just load the data
	// directly onto the model grid.
	//------------------------------------------------------------------------------
	if(same==1){
		
		if(rank == 0){ printf("Files have same dimensions and grid spacing. Loading initial conditions from file %s\n",myfilename);}
		
		load_from_output(myfilename,"u-wind",us,ums,time);
		load_from_output(myfilename,"v-wind",vs,vms,time);
		load_from_output(myfilename,"w-wind",ws,wms,time);
		load_from_output(myfilename,"theta",ths,thms,time);
	
		if(USE_MICROPHYSICS){
		
			load_from_output(myfilename,"qv",qvs,qvms,time);
			load_from_output(myfilename,"qc",qcs,qcms,time);
			load_from_output(myfilename,"qr",qrs,qrms,time);
		}
	//------------------------------------------------------------------------------
	// If the files have different dimensions or grid spacing, need to interpolate
	// input fields to model grid.
	//------------------------------------------------------------------------------
	} else {

		double *zlevs,*var,*var_interpz,*myLons,*myLats;
		//------------------------------------------------------------------------------
		// Allocate temporary memory on the root process for interpolated fields. Don't
		// try to access these variables on the other processes!
		//------------------------------------------------------------------------------
		if(rank==0){
			
			printf("Files have different dimensions or grid spacing.\nLoading initial conditions from file %s \nwith dimension %lu x %lu x %lu and grid spacing dx = %f at time = %lu\n",
				myfilename,xdim,ydim,zdim,grid_spacing[0],time);
			
			var = (double*)calloc(xdim*ydim*zdim,sizeof(double));
			var_interpz = (double*)calloc(xdim*ydim*NZ,sizeof(double));
	
			zlevs = get_data2(myfilename,"zu",zdim);
			myLons = get_data2(myfilename,"lon",xdim);
			myLats = get_data2(myfilename,"lat",ydim);
		}
		//------------------------------------------------------------------------------
		// Get data, interpolate it, distribute it to each process
		//------------------------------------------------------------------------------
		load_interpolate_from_output(myfilename,"u-wind",time,var,var_interpz,us,ums,zsu,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		load_interpolate_from_output(myfilename,"v-wind",time,var,var_interpz,vs,vms,zsu,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);	
		load_interpolate_from_output(myfilename,"w-wind",time,var,var_interpz,ws,wms,zsw,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		load_interpolate_from_output(myfilename,"theta",time,var,var_interpz,ths,thms,zsu,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		
		if(rank==0){ free(zlevs); free(var); free(var_interpz); free(myLons); free(myLats);}
	}
	
	
}

/*********************************************************************
* Initialize model perturbation fields from input file. Will determine
* whether or not interpolation is required.
*
* myfilename - file containing input fields for perturbations
* time		 - time at which to initialize
**********************************************************************/
void initialize_from_output_serial(const char *myfilename,size_t time){
	
	size_t dims[3],xdim,ydim,zdim;
	double interp_dx,interp_dy,interp_dz;
	float grid_spacing[3];
	int same = 0;
	
	//------------------------------------------------------------------------------
	// Test files for same dimensions and grid spacing
	//------------------------------------------------------------------------------	
	get_dims(myfilename,"x","y","z",dims);
	get_grid_spacing(myfilename,"dx","dy","dz",grid_spacing);

	interp_dx = grid_spacing[0];
	interp_dy = grid_spacing[1];
	interp_dz = grid_spacing[2];
	
	xdim = dims[0];
	ydim = dims[1];
	zdim = dims[2];

	//------------------------------------------------------------------------------
	// If all dimensions and grid spacing is the same, simply load input
	//------------------------------------------------------------------------------	
	if(xdim == NX && ydim == NY && zdim == NZ && interp_dx == dx && interp_dy == dy && interp_dz == dz){
		
		printf("Files have same dimensions and grid spacing. Loading initial conditions from file %s\n",myfilename);
	
		get_data_at_time(myfilename,"u-wind",time,&U(0,0,0));
		get_data_at_time(myfilename,"v-wind",time,&V(0,0,0));
		get_data_at_time(myfilename,"w-wind",time,&W(0,0,0));
		get_data_at_time(myfilename,"theta",time,&TH(0,0,0));
	
		if(USE_MICROPHYSICS){
		
			get_data_at_time(myfilename,"qv",time,&QV(0,0,0));
			get_data_at_time(myfilename,"qc",time,&QC(0,0,0));
			get_data_at_time(myfilename,"qr",time,&QR(0,0,0));
		}
	//------------------------------------------------------------------------------
	// If the files have different dimensions or grid spacing, need to interpolate
	// input fields to model grid.
	//------------------------------------------------------------------------------
	} else {

		double *zlevs,*var,*var_interpz,*myLons,*myLats;
		//------------------------------------------------------------------------------
		// Allocate temporary memory on the root process for interpolated fields.
		//------------------------------------------------------------------------------			
		printf("Files have different dimensions or grid spacing.\nLoading initial conditions from file %s \nwith dimension %lu x %lu x %lu and grid spacing dx = %f at time = %lu\n",
		myfilename,xdim,ydim,zdim,grid_spacing[0],time);
		
		var = (double*)calloc(xdim*ydim*zdim,sizeof(double));
		var_interpz = (double*)calloc(xdim*ydim*NZ,sizeof(double));

		zlevs = get_data2(myfilename,"zu",zdim);
		myLons = get_data2(myfilename,"lon",xdim);
		myLats = get_data2(myfilename,"lat",ydim);

		//------------------------------------------------------------------------------
		// Get data, interpolate it
		//------------------------------------------------------------------------------		
		interpolate_from_output(myfilename,"u-wind",time,var,var_interpz,&U(0,0,0), zsu,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		interpolate_from_output(myfilename,"v-wind",time,var,var_interpz,&V(0,0,0), zsu,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);		
		interpolate_from_output(myfilename,"w-wind",time,var,var_interpz,&W(0,0,0), zsu,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		interpolate_from_output(myfilename,"theta", time,var,var_interpz,&TH(0,0,0),zsu,myLons,myLats,zlevs,xdim,ydim,zdim,interp_dx,interp_dy);
		
		free(zlevs); free(var); free(var_interpz); free(myLons); free(myLats);
	}
	
	//------------------------------------------------------------------------------
	// Initialize ...
	//------------------------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
		
		UM(i,j,k) = U(i,j,k);
		VM(i,j,k) = V(i,j,k);
		WM(i,j,k) = W(i,j,k);
		THM(i,j,k) = TH(i,j,k);
	}}}
	
	if(USE_MICROPHYSICS){
		
		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
		for(int k=0;k<NZ;k++){
			
			QVM(i,j,k) = QV(i,j,k);
			QCM(i,j,k) = QC(i,j,k);
			QRM(i,j,k) = QR(i,j,k);
			
		}}}
	}
}

/********************************************************
*
* 
*
*********************************************************/
void initialize_from_ncar(){

	/********************************************************
	* Get dimensions of data in files to determine how 
	* much memory to allocate
	*********************************************************/
	size_t dims[3];

	get_dims(geoheight_file,"lon","lat","level",dims);

	size_t xdim = dims[0], ydim = dims[1], zdim =  dims[2];

	/********************************************************
	* Create arrays large enough to hold the variable
	*********************************************************/
	int size = xdim*ydim*zdim;

	double * hgt2 = (double *)malloc(size*sizeof(double));
	double * uz2 = (double *)malloc(size*sizeof(double));
	double * vz2 = (double *)malloc(size*sizeof(double));
	double * tz2 = (double *)malloc(size*sizeof(double));

	double * hgt = (double *)malloc(size*sizeof(double));
	double * uz = (double *)malloc(size*sizeof(double));
	double * vz = (double *)malloc(size*sizeof(double));
	double * tz = (double *)malloc(size*sizeof(double));
	
	double * lat2 = (double *)malloc(ydim*sizeof(double));
	double * lon = (double *)malloc(xdim*sizeof(double));
	double * lat = (double *)malloc(ydim*sizeof(double));
	double * levs = (double *)malloc(zdim*sizeof(double));

	double * topo2 = (double *)malloc(xdim*ydim*sizeof(double));
	double * topozinterp = (double *)malloc(xdim*ydim*NZ*sizeof(double));

	/********************************************************
	* Get data from file
	*********************************************************/
	get_data(geoheight_file,"hgt", size, hgt2);
	get_data(uwind_file,"uwnd", size, uz2);
	get_data(vwind_file,"vwnd", size, vz2);
	get_data(temp_file,"air", size, tz2);
	get_data(topo_file,"hgt", xdim*ydim, topo2);

	/********************************************************
	* Reverse the y-coordinate
	*********************************************************/
	for(size_t i=0;i<xdim;i++){
	for(size_t j=0;j<ydim;j++){

		for(size_t k=0;k<zdim;k++){

			uz[d(i,j,k)] = uz2[d(i,ydim-1-j,k)];
			vz[d(i,j,k)] = vz2[d(i,ydim-1-j,k)];
			hgt[d(i,j,k)] = hgt2[d(i,ydim-1-j,k)];
			tz[d(i,j,k)] = tz2[d(i,ydim-1-j,k)];
		}

		for(int k=0;k<NZ;k++)
			topozinterp[d(i,j,k)] = topo2[d2(i,ydim-1-j)];
	}}

	/********************************************************
	* Get and process coordinates from input data
	*********************************************************/
	get_data(geoheight_file,"lat", ydim, lat2);
	get_data(geoheight_file,"lon", xdim, lon);
	get_data(geoheight_file,"level", zdim, levs);

	for(size_t i=0;i<xdim;i++)
		if(lon[i] >= 290)
			lon[i] = lon[i] - 360;

	for(size_t i=0;i<ydim;i++){ lat[i] = lat2[ydim-1-i]; }	// reverse y-coordinate

	for(size_t k=0;k<zdim;k++){ levs[k] = levs[k]*100;}

	/********************************************************
	* Convert temperature to potential temperature
	*********************************************************/
	for(size_t k=0;k<zdim;k++)
		for(size_t j=0;j<ydim;j++)
			for(size_t i=0;i<xdim;i++)
				tz[d(i,j,k)] = tz[d(i,j,k)]*pow((100000./levs[k]),(287./1004.));

	/********************************************************
	* Vertical interpolation
	*********************************************************/
	double * z 		  = (double *)malloc(NZ*sizeof(double));

	z[0] = 0;
	for(int k=1;k<NZ;k++){ z[k] = ((double)k-0.5)*dz;}

	double * uzinterp = vert_interpolate(z,hgt,uz,xdim,ydim,zdim,NZ);
	double * vzinterp = vert_interpolate(z,hgt,vz,xdim,ydim,zdim,NZ);
	double * tzinterp = vert_interpolate(z,hgt,tz,xdim,ydim,zdim,NZ);

	/********************************************************
	* Horizontal interpolation
	*********************************************************/
	double lonoffset = 90;
	double latoffset = 1;
	double dlat = 2.5;
	double dlon = 2.5;

	horz_interpolate(uzinterp,&IUBAR(0,0,0),xdim,ydim,NZ,!false,!true,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(vzinterp,&IVBAR(0,0,0),xdim,ydim,NZ,!true,!false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(tzinterp,&ITHBAR(0,0,0),xdim,ydim,NZ,!true,!true,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(topozinterp,&ITOPO(0,0,0),xdim,ydim,NZ,!true,!true,dlat,dlon,lonoffset,latoffset);

	/********************************************************
	* Construct coordinate arrays for output data
	*********************************************************/
	outLons[0] = lon[0] + lonoffset + dx/meters_per_degree;
	outLats[0] = lat[0] + dy/meters_per_degree;	

	for(int i=1;i<NX;i++){ outLons[i] = outLons[i-1] + dx/meters_per_degree;}
	for(int i=1;i<NY;i++){ outLats[i] = outLats[i-1] + dy/meters_per_degree;}

	/********************************************************
	* Free up allocated memory
	*********************************************************/
	free(uzinterp); free(vzinterp); free(tzinterp); free(z);
	free(hgt2); free(uz2); free(vz2); free(tz2);
	free(hgt); free(uz); free(vz); free(tz);
	free(levs); free(lat); free(lon); free(lat2);
	free(topozinterp); free(topo2);
}

/********************************************************
*
* 
*
*********************************************************/
void initialize_from_era(double * z){

	/********************************************************
	* Get dimensions of data in files to determine how 
	* much memory to allocate
	*********************************************************/
	size_t dims[3];

	get_dims(datafile,"lon","lat","level",dims);

	size_t xdim = dims[0]; 
	size_t ydim = dims[1];
	size_t zdim = dims[2];

	int size = xdim*ydim*zdim;

	double *piIn = (double *)malloc(zdim*sizeof(double));

	//printf("Dimensions of input data are %zu %zu %zu\n",xdim,ydim,zdim);

	/********************************************************
	* Get data from file
	*********************************************************/
	double * hgt = get_data2(datafile,"hgt", size);
	double * uz  = get_data2(datafile,"uwnd", size);
	double * vz  = get_data2(datafile,"vwnd", size);
	double * tz  = get_data2(datafile,"temp", size);
	double * qz  = get_data2(datafile2,"qv", size);
	double * topo2 = get_data2(datafile,"zsfc", xdim*ydim);
	
	/********************************************************
	* Reverse the y- and z-coordinate
	*********************************************************/
	reverse_yz_coord(uz,xdim,ydim,zdim);
	reverse_yz_coord(vz,xdim,ydim,zdim);
	reverse_yz_coord(hgt,xdim,ydim,zdim);
	reverse_yz_coord(tz,xdim,ydim,zdim);
	reverse_yz_coord(qz,xdim,ydim,zdim);

	double * topozinterp = (double *)malloc(xdim*ydim*NZ*sizeof(double));

	for(size_t i=0;i<xdim;i++){
	for(size_t j=0;j<ydim;j++){
	for(size_t k=0;k < NZ;k++){
		topozinterp[d(i,j,k)] = topo2[d2(i,ydim-1-j)];
	}}}

	if(SHIFT_PRIME_MERIDIAN){

		flip_array(hgt,xdim,ydim,zdim);
		flip_array(uz,xdim,ydim,zdim);
		flip_array(vz,xdim,ydim,zdim);
		flip_array(tz,xdim,ydim,zdim);
		flip_array(qz,xdim,ydim,zdim);
		flip_array(topozinterp,xdim,ydim,NZ);
	}

	/********************************************************
	* Get and process coordinates from input data
	*********************************************************/
	double * lat2  = get_data2(datafile,"latitude", ydim);
	double * lon   = get_data2(datafile,"longitude", xdim);
	double * levs2 = get_data2(datafile,"level", zdim);

	double * lat  = (double *)malloc(ydim*sizeof(double));
	double * levs = (double *)malloc(zdim*sizeof(double));
	
	for(size_t i=0;i<ydim;i++){ lat[i] = lat2[ydim-1-i];}			// reverse y-coordinate

	for(size_t k=0;k<zdim;k++){ levs[k] = levs2[zdim-1-k]*100.0;}// printf("%d pres %f\n",k,levs[k]);}	// reverse z-coordinate

	for(size_t k=0;k<zdim;k++){ piIn[k] = pow((levs[k]/p0),Rd/cp);}

#if 0
	for(int i=0;i<xdim;i++){
	for(int j=0;j<ydim;j++){
		
		if(lat[j] > 19 && lat[j] < 21 && lon[i] > 88 && lon[i] < 90){
		
			for(int k=0;k<zdim;k++){
		
				printf("%d %f %f %f %f\n",k,levs[k],lat[j],lon[i],qz[d(i,j,k)]*1000);
			}	
		}
	}}
#endif
	/********************************************************
	* Convert temperature to potential temperature
	*********************************************************/
	for(size_t k=0;k<zdim;k++)
		for(size_t j=0;j<ydim;j++)
			for(size_t i=0;i<xdim;i++)
				tz[d(i,j,k)] = tz[d(i,j,k)]*pow((100000.0/levs[k]),(287.0/1004.0));

	/********************************************************
	* Vertical interpolation
	*********************************************************/
	//double * z = (double *)malloc(NZ*sizeof(double));

	// construct height coordinate
	//for(int k=0;k<NZ;k++){ z[k] = ((double)k-0.5)*dz;}

	double * uzinterp = vert_interpolate(z,hgt,uz,xdim,ydim,zdim,NZ);
	double * vzinterp = vert_interpolate(z,hgt,vz,xdim,ydim,zdim,NZ);
	double * tzinterp = vert_interpolate(z,hgt,tz,xdim,ydim,zdim,NZ);
	double * qzinterp = vert_interpolate(z,hgt,qz,xdim,ydim,zdim,NZ);
	double * pzinterp = vert_interpolate_1d(z,hgt,piIn,xdim,ydim,zdim,NZ);

	/********************************************************
	* Horizontal interpolation
	*********************************************************/
	double dlat = lat[3]-lat[2];
	double dlon = lon[3]-lon[2];
	//printf("dlat = %f dlon = %f\n",dlat,dlon);
	horz_interpolate(uzinterp,&IUBAR(0,0,0),xdim,ydim,NZ,true,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(vzinterp,&IVBAR(0,0,0),xdim,ydim,NZ,false,true,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(tzinterp,&ITHBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(qzinterp,&IQBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(pzinterp,&IPBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(topozinterp,&ITOPO(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);

	/********************************************************
	* Construct coordinate arrays for output data
	*********************************************************/
	double lon_shift;
	
	if(SHIFT_PRIME_MERIDIAN){ lon_shift = 180;} else { lon_shift = 0;}

	outLons[0] = lon[0] + lonoffset + dx/meters_per_degree - lon_shift;
	outLats[0] = lat[0] + latoffset + dy/meters_per_degree;

	for(int i=1;i<NX;i++){ outLons[i] = outLons[i-1] + dx/meters_per_degree;}
	for(int i=1;i<NY;i++){ outLats[i] = outLats[i-1] + dy/meters_per_degree;}

	/********************************************************
	* Free up allocated memory
	*********************************************************/
	free(uzinterp); free(vzinterp); free(tzinterp); free(qzinterp); free(pzinterp);
	free(hgt); free(uz); free(vz); free(tz); free(qz);
	free(levs); free(lat); free(lon); free(lat2); free(levs2); 
	free(topozinterp); free(topo2);	free(piIn);
}

/********************************************************
*
* 
*
*********************************************************/
void initialize_from_ERA5(double * z){

	/********************************************************
	* Get dimensions of data in files to determine how 
	* much memory to allocate
	*********************************************************/
	size_t dims[3];

	
	//const char datafile1[100] = "/Users/michaeldiaz/Desktop/output.nc";
	const char datafile1[100] = "/Users/michaeldiaz/Desktop/output.nc";
	//const char datafile3[100] = "/Users/michaeldiaz/Desktop/era5jan2018g.nc";
	const char datafile2[100] = "/Users/michaeldiaz/Desktop/era5surface.nc";

	get_dims(datafile1,"longitude","latitude","level",dims);

	size_t xdim = dims[0]; 
	size_t ydim = dims[1];
	size_t zdim = dims[2];

	int size = xdim*ydim*zdim;

	double *piIn = (double *)malloc(zdim*sizeof(double));

	//printf("Dimensions of input data are %zu %zu %zu\n",xdim,ydim,zdim);

	/********************************************************
	* Get data from file
	*********************************************************/
	double * hgt = get_data2(datafile1,"z", size);
	double * uz  = get_data2(datafile1,"u", size);
	double * vz  = get_data2(datafile1,"v", size);
	double * tz  = get_data2(datafile1,"t", size);
	double * qz  = get_data2(datafile1,"q", size);
	double * topo2 = get_data2(datafile2,"z", xdim*ydim);
	
	for(int i=0;i<size;i++){ hgt[i] /= 9.80665;}
	for(int i=0;i<xdim*ydim;i++){ topo2[i] /= 9.80665;}
	
	/********************************************************
	* Reverse the y- and z-coordinate
	*********************************************************/
	reverse_yz_coord(uz,xdim,ydim,zdim);
	reverse_yz_coord(vz,xdim,ydim,zdim);
	reverse_yz_coord(hgt,xdim,ydim,zdim);
	reverse_yz_coord(tz,xdim,ydim,zdim);
	reverse_yz_coord(qz,xdim,ydim,zdim);

	double * topozinterp = (double *)malloc(xdim*ydim*NZ*sizeof(double));

	for(size_t i=0;i<xdim;i++){
	for(size_t j=0;j<ydim;j++){
	for(size_t k=0;k < NZ;k++){
		topozinterp[d(i,j,k)] = topo2[d2(i,ydim-1-j)];
	}}}

	if(SHIFT_PRIME_MERIDIAN){

		flip_array(hgt,xdim,ydim,zdim);
		flip_array(uz,xdim,ydim,zdim);
		flip_array(vz,xdim,ydim,zdim);
		flip_array(tz,xdim,ydim,zdim);
		flip_array(qz,xdim,ydim,zdim);
		flip_array(topozinterp,xdim,ydim,NZ);
	}

	/********************************************************
	* Get and process coordinates from input data
	*********************************************************/
	double * lat2  = get_data2(datafile1,"latitude", ydim);
	double * lon   = get_data2(datafile1,"longitude", xdim);
	double * levs2 = get_data2(datafile1,"level", zdim);

	double * lat  = (double *)malloc(ydim*sizeof(double));
	double * levs = (double *)malloc(zdim*sizeof(double));
	
	for(size_t i=0;i<ydim;i++){ lat[i] = lat2[ydim-1-i];}			// reverse y-coordinate

	for(size_t k=0;k<zdim;k++){ levs[k] = levs2[zdim-1-k]*100.0;}// printf("%d pres %f\n",k,levs[k]);}	// reverse z-coordinate

	for(size_t k=0;k<zdim;k++){ piIn[k] = pow((levs[k]/p0),Rd/cp);}

	/********************************************************
	* Convert temperature to potential temperature
	*********************************************************/
	for(size_t k=0;k<zdim;k++)
		for(size_t j=0;j<ydim;j++)
			for(size_t i=0;i<xdim;i++)
				tz[d(i,j,k)] = tz[d(i,j,k)]*pow((100000.0/levs[k]),(287.0/1004.0));

	/********************************************************
	* Vertical interpolation
	*********************************************************/
	//double * z = (double *)malloc(NZ*sizeof(double));

	// construct height coordinate
	//for(int k=0;k<NZ;k++){ z[k] = ((double)k-0.5)*dz;}

	double * uzinterp = vert_interpolate(z,hgt,uz,xdim,ydim,zdim,NZ);
	double * vzinterp = vert_interpolate(z,hgt,vz,xdim,ydim,zdim,NZ);
	double * tzinterp = vert_interpolate(z,hgt,tz,xdim,ydim,zdim,NZ);
	double * qzinterp = vert_interpolate(z,hgt,qz,xdim,ydim,zdim,NZ);
	double * pzinterp = vert_interpolate_1d(z,hgt,piIn,xdim,ydim,zdim,NZ);

	/********************************************************
	* Horizontal interpolation
	*********************************************************/
	double dlat = lat[3]-lat[2];
	double dlon = lon[3]-lon[2];
	//printf("dlat = %f dlon = %f\n",dlat,dlon);
	horz_interpolate(uzinterp,&IUBAR(0,0,0),xdim,ydim,NZ,true,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(vzinterp,&IVBAR(0,0,0),xdim,ydim,NZ,false,true,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(tzinterp,&ITHBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(qzinterp,&IQBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(pzinterp,&IPBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(topozinterp,&ITOPO(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);

	/********************************************************
	* Construct coordinate arrays for output data
	*********************************************************/
	double lon_shift;
	
	if(SHIFT_PRIME_MERIDIAN){ lon_shift = 180;} else { lon_shift = 0;}

	outLons[0] = lon[0] + lonoffset + dx/meters_per_degree - lon_shift;
	outLats[0] = lat[0] + latoffset + dy/meters_per_degree;

	for(int i=1;i<NX;i++){ outLons[i] = outLons[i-1] + dx/meters_per_degree;}
	for(int i=1;i<NY;i++){ outLats[i] = outLats[i-1] + dy/meters_per_degree;}

	/********************************************************
	* Free up allocated memory
	*********************************************************/
	free(uzinterp); free(vzinterp); free(tzinterp); free(qzinterp); free(pzinterp);
	free(hgt); free(uz); free(vz); free(tz); free(qz);
	free(levs); free(lat); free(lon); free(lat2); free(levs2); 
	free(topozinterp); free(topo2);	free(piIn);
}


/********************************************************
* Reverse the y- and z-coordinate 
*
*********************************************************/
void reverse_yz_coord(double * var, int xdim, int ydim, int zdim){

	int size = xdim*ydim*zdim;

	double * var2 = (double *)malloc(size*sizeof(double));

	/********************************************************
	* 
	*********************************************************/
	for(int i=0;i<xdim;i++){
	for(int j=0;j<ydim;j++){
	for(int k=0;k<zdim;k++){

		var2[d(i,j,k)] = var[d(i,ydim-1-j,zdim-1-k)];

	}}}

	for(int i=0;i<size;i++){ var[i] = var2[i];}

	free(var2);
}


/********************************************************
* 
*
*********************************************************/
void flip_array(double * var, int xdim, int ydim, int zdim){

	int size = xdim*ydim*zdim;

	double * var2 = (double *)malloc(size*sizeof(double));

	int half = xdim / 2;

	/********************************************************
	* 
	*********************************************************/
	for(int i=0;i<half;i++){
	for(int j=0;j<ydim;j++){
	for(int k=0;k<zdim;k++){

		var2[d(i+half,j,k)] = var[d(i,j,k)];
		var2[d(i,j,k)] = var[d(i+half,j,k)];

	}}}

	for(int i=0;i<size;i++){ var[i] = var2[i];}

	free(var2);
}

#if 0

//			printf("%lu %lu %lu %f %f %f\n",xdim,ydim,zdim,grid_spacing[0],grid_spacing[1],grid_spacing[2]);
	
//			for(int i=0;i<zdim;i++){ printf("%f\n",zlevs[i]);}
		}
			//get_data_at_time(myfilename,"u-wind",time,var,xdim,ydim,zdim);
			
			//horz_interpolate_from_model(double *uz,double *out,int xdim,int ydim,int zdim,bool xstagger,bool ystagger,double dlat,double dlon,double lonoffset,double latoffset)
			
			//vert_interpolate_1d_from_model(zsu, zlevs, var, xdim,ydim,zdim,NZ,var_interpz);
			
			//horz_interpolate_from_model(var_interpz,iubar,xdim,ydim,NZ,false,false,grid_spacing[0]/meters_per_degree,grid_spacing[1]/meters_per_degree,myLons[0]-outLons[0],myLats[0]-outLats[0]);

			for(int i=0;i<xdim;i++){
			for(int j=0;j<ydim;j++){
	
				for(int k=0;k<zdim;k++){
				
					if(j==ydim/2 && i==xdim/2){
				
						printf("%d %f %f\n",k,zlevs[k],var[k+j*zdim+i*zdim*ydim]);
					}
				}
				
				for(int k=0;k<NZ;k++){
				
					if(j==ydim/2 && i==xdim/2){
				
						printf("%d %f %f\n",k,zsu[k],var_interpz[k+j*NZ+i*NZ*ydim]);
					}
				}
			}}
#endif
	

