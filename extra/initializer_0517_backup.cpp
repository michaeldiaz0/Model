#include "stdafx.h"
#include "interpolate.h"
#include "surface.h"
#include "fluxes.h"
#include "damping.h"

#define d(i,j,k) (xdim*ydim*(k) + xdim*(j) + i)


double u[NX][NY][NZ];
double v[NX][NY][NZ];
double w[NX][NY][NZ];
double th[NX][NY][NZ];
double pi[NX][NY][NZ];

#if !PARALLEL && !ENERGY

	double up[NX][NY][NZ]; double um[NX][NY][NZ];
	double vp[NX][NY][NZ]; double vm[NX][NY][NZ];
	double wp[NX][NY][NZ]; double wm[NX][NY][NZ];
	double thp[NX][NY][NZ]; double thm[NX][NY][NZ];

	double ubar[NX][NY][NZ];
	double vbar[NX][NY][NZ];
	double wbar[NX][NY][NZ];
	double thbar[NX][NY][NZ];
	double qbar[NX][NY][NZ];
	
	double friction[NX][NY][NZ];
	double istopo[NX][NY][NZ];
	double uistopo[NX][NY][NZ];
	double vistopo[NX][NY][NZ];
	double topo[NX][NY][NZ];

#endif

double *iubar,*ivbar,*iwbar,*ithbar,*iqbar,*ipbar;
double *itopo,*iistopo,*iuistopo,*ivistopo,*ifriction;


double zu[NZ];	// height of u-velocity
double zw[NZ];	// height of w-velocity
double rhou[NZ];// base state density at u-velocity level
double rhow[NZ];// base state density at w-velocity level
double tb[NZ];	// base state potential temperature
double tbw[NZ];	// base state potential temperature at w-velocity level
double tbv[NZ];	// base state virtual potential temperature
double qb[NZ];	// base state specific humidity
double qbs[NZ];	// base state saturation mixing ratio
double ub[NY][NZ];	// base state u-velocity
double vb[NZ];	// base state v-velocity
double pib[NZ];	// base state non-dimensional pressure
double one_d_rhou[NZ];
double one_d_rhow[NZ];

double zsu[NZ];
double zsw[NZ];
double mu[NZ];
double mw[NZ];

double f[NY];
double dfdy[NY];

double rhoavg2d[NX][NY];
double outLats[NY];
double outLons[NX];
int htopo[NX][NY];
const double meters_per_degree = 111000;
int bigcounter;
int linear_lon;

double get_QV_Sat(double temperature,double pressure);
void stretched_grid(double*,double*,double*,double*,double,int);
void reverse_yz_coord(double*, int, int, int);
void flip_array(double *, int, int, int);

/*********************************************************************
*
*
**********************************************************************/
void initialize2(){

	double xk = Rd/cp;
	double pisfc;

	int size = NX*NY*NZ;

	rhoavg = 0;
	mtime = 0;
	bigcounter = 0;

	pisfc = pow((pressfc/p0),xk);
	
#if !PERIODIC_BOUNDARIES
	init_fftw(NX,NY,NZ);
#else
	init_fftw(NX-6,NY,NZ);
#endif
	
	/********************************************
	* initialize common data
	*********************************************/
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		u[i][j][k] = 0;
		v[i][j][k] = 0;
		th[i][j][k] = 0;
		w[i][j][k] = 0;
		pi[i][j][k] = 0;
		
	}}}

	/********************************************
	* Calculate height of each level
	*********************************************/
	for(int k=0;k<NZ;k++){

		zu[k] = ((double)k-0.5)*dz;
		zw[k] = ((double)k-1)*dz;
	}

	stretched_grid(&zsu[0],&mu[0],&zsw[0],&mw[0],50,1);

	/********************************************
	* initialize data specific to parallel version
	*********************************************/
	#if PARALLEL || ENERGY
	
		iubar = (double*) calloc(size,sizeof(double));
		ivbar = (double*) calloc(size,sizeof(double));
		iwbar = (double*) calloc(size,sizeof(double));
		ithbar = (double*) calloc(size,sizeof(double));
		ipbar = (double*) calloc(size,sizeof(double));
		iqbar = (double*) calloc(size,sizeof(double));
	
		itopo = (double*) calloc(size,sizeof(double));
		iuistopo = (double*) calloc(size,sizeof(double));
		ivistopo = (double*) calloc(size,sizeof(double));
		iistopo = (double*) calloc(size,sizeof(double));
		ifriction = (double*) calloc(size,sizeof(double));
	
		if(STRETCHED_GRID){initialize_from_era(&zsu[0]);} 
		else {			   initialize_from_era(&zu[0]);}
		
		linear_lon = get_point_from_lon(90);
		
		if(MERIDIONAL_CROSS_SECTION){
		
			for(int i=0;i<NX;i++){
			for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
		
				IUBAR(i,j,k) = IUBAR(linear_lon,j,k);
				IVBAR(i,j,k) = 0;
				ITHBAR(i,j,k) = ITHBAR(linear_lon,j,k);
				IQBAR(i,j,k) = IQBAR(linear_lon,j,k);
			}}}
			
		}
		
	/********************************************
	* initialize data specific to serial version
	*********************************************/
	#else
		initialize_flux_cells(NY,NZ);
		initialize_microphysics_cells(NY,NZ);

		// basic state 3D pressure
		ipbar = (double*) calloc(size,sizeof(double));

		if(STRETCHED_GRID){ initialize_from_era(&zsu[0]);} 
		else { 				initialize_from_era(&zu[0]);}

		initialize_landsea(landseaMaskFile);

		if(USE_TURBULENT_STRESS){ init_kmix(size);}

		if(OUTPUT_DIFFUSION_TEND){ init_damping(NX,NY,NZ);}
		//---------------------------------
		// 72 b = 8 jn = 4 js = 14 convergence, wrong phase velocity direction
		// 73 b = 8 jn = 4 js = 14 wrong phase velocity direction, no convergence
		// 74 b = 8 jn = 4 js = 14 more amplitude at upper levels (compare heating-forced experiment)
		// 75 b = 8 jn = 4 js = 14 more amplitude at upper levels (compare heating-forced experiment)
		// 76 b = 8 jn = 4 js = 14	
		// 77 b = 8 jn = 4 js = 14 beginning here, more amplitude at upper levels (compare heating-forced experiment)
		// 78 b = 8 jn = 5 js = 13
		// 79 b = 8 jn = 5 js = 12 similar to farther east
		// 80 b = 9 jn = 5 js = 11
		// 81 b = 9 jn = 6 js = 10
		// 82 b = 10 jn = 6 js = 10
		// 83 b = 10 jn = 6 js = 10
		// 84 b = 10 jn = 6 js = 10
		// 85 b = 10 jn = 6 js = 10
		// 86 b = 10 jn = 6 js = 10 latoffset = 89.75
		// 87 b = 10 jn = 6 js = 10
		// 88 b = 10 jn = 6 js = 10
		// 89 b = 10 jn = 6 js = 10
		// 90 b = 10 jn = 6 js = 10
		// 90.5 b = 10 jn = 6 js = 10	periodic convergence
		// 91 b = 10 jn = 5 js = 10	convergence
		// 92 b = 10 jn = 5 js = 10 no convergence b = 11 convergence
		// 93 b = 10 jn = 5 js = 10 (or b=12 _2) (or b=11 _3)
		// 94 b = 10 jn = 5 js = 10 (or b=11 _2)
		//---------------------------------
		linear_lon = get_point_from_lon(87);

		printf("lon = %d\n",get_point_from_lon(87));

		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
		for(int k=0;k<NZ;k++){

			up[i][j][k] = u[i][j][k];
			um[i][j][k] = u[i][j][k];
			vp[i][j][k] = v[i][j][k];
			vm[i][j][k] = v[i][j][k];
			thp[i][j][k] = th[i][j][k];
			thm[i][j][k] = th[i][j][k];

		#if !HYDROSTATIC
			wp[i][j][k] = w[i][j][k];
			wm[i][j][k] = w[i][j][k];
		#endif

		#if RESTING_BASIC_STATE
			IUBAR(i,j,k) = 0;
			IVBAR(i,j,k) = 0;
			ITHBAR(i,j,k) = ITHBAR(NX/2,NY/2,k);
			IQBAR(i,j,k) = IQBAR(NX/2,NY/2,k);
		#elif MERIDIONAL_CROSS_SECTION
			IUBAR(i,j,k) = IUBAR(linear_lon,j,k);
			IVBAR(i,j,k) = 0;
			ITHBAR(i,j,k) = ITHBAR(linear_lon,j,k);
			IQBAR(i,j,k) = IQBAR(linear_lon,j,k);
		#endif
			
		}}}
			
	#endif
		
	/***********************************************************************
	* For linearized model, create meridionally uniform conditions on the
	* northern and southern boundaries	
	************************************************************************/
	if(MERIDIONAL_CROSS_SECTION && ISLINEAR){

		double llat = 15.0;//10.446912;
		double hlat = 26.663129;

		int b1 = 0;//10;//3;
		int b2 = 10;//10;//12;
#if 0		
		for(int j=1;j<NY;j++){
			
			if(outLats[j-1] <= llat && outLats[j] > llat){
				
				//printf("lat1 = %f\n",outLats[j]);
				
				if( llat - outLats[j-1] <= outLats[j] - llat ){
					
					b1 = j - 1;
				} else {
					b1 = j;
				}
			}
			
			if(outLats[j-1] <= hlat && outLats[j] > hlat){
				
				//printf("lat2 = %f\n",outLats[j]);
				
				if( hlat - outLats[j-1] <= outLats[j] - hlat ){
					
					//b2 = NY - j + 1;
				} else {
					//b2 = NY - j;
				}
			}
		}
#endif
		//b1 = 0;
		int b[NZ];
		int b3[NZ];
		//b1 = 3;
		//b2 = 0;
		double ubarmin = 100;
		int jmin = 0;
		double dU2dy;
#if 0		
		for(int k=0;k<NZ;k++){
			
			ubarmin = 100;
			printf("%d ",k);
			for(int j=1;j<NY;j++){
			
				if(ubarmin > UBAR(10,j,k)){
				
					ubarmin = UBAR(10,j,k);
					jmin = j;
				}
				
				dU2dy = 2*UBAR(10,j,k) - UBAR(10,j+1,k) - UBAR(10,j-1,k);
				printf("%0.1f ",dU2dy);
			}
			printf("\n");
			if(jmin >23 && UBAR(10,jmin,k) < 5){ b[k] = jmin;}
			else { b[k] = NY-5;}
	
			//printf("%d %f %f\n",k,outLats[jmin],UBAR(10,jmin,k));
		}
#endif
#if 0
		double dUdy1,dUdy2;
		
		for(int k=0;k<NZ;k++){
			printf("%d\n",k);
			b[k] = -1;
			b3[k] = -1;
			for(int j=1;j<NY-1;j++){
				
				dUdy1 = UBAR(10,j+1,k) - UBAR(10,j  ,k);
				dUdy2 = UBAR(10,j  ,k) - UBAR(10,j-1,k);
				
				if( (dUdy1 < 0 && dUdy2 > 0) || (dUdy1 > 0 && dUdy2 < 0) ){
					
					if(UBAR(10,j,k) < -3 && outLats[j]>20){
					
						printf("%d %f %f\n",j,outLats[j],UBAR(10,j,k) );
						b[k] = j;
					}
				}
				
				if( (dUdy1 < 0 && dUdy2 > 0) || (dUdy1 > 0 && dUdy2 < 0) ){
					
					if(UBAR(10,j,k) > -6 && outLats[j]<18 && outLats[j]>10){
					
						printf("%d %f %f\n",j,outLats[j],UBAR(10,j,k) );
						b3[k] = j;
					}
				}
				
			}
			
		}

		//for(int k=0;k<NZ;k++){	
		
		for(int k=0;k<NZ;k++){ printf("%d %d %f %f\n",k,b[k],outLats[b[k]],UBAR(10,b[k],k));}
		
		//printf("lat = %f %f\n",outLats[b1],outLats[NY-b2]);
#endif		
		for(int i=0;i<NX;i++){
		for(int k=0;k<NZ;k++){

			for(int j=0;j<b1;j++){
				
				IUBAR(i,j,k) = IUBAR(i,b1,k);
				ITHBAR(i,j,k) = ITHBAR(i,b1,k);
				IQBAR(i,j,k) = IQBAR(i,b1,k);
			}
#if 1
			for(int j=NY-b2;j<NY;j++){

				IUBAR(i,j,k) = IUBAR(i,NY-b2,k);
				ITHBAR(i,j,k) = ITHBAR(i,NY-b2,k);
				IQBAR(i,j,k) = IQBAR(i,NY-b2,k);
			}
#endif
			
#if 0
			double vort = -(IUBAR(i,b2+1,k) - IUBAR(i,b2,k)) * one_d_dy + 2.0*fc*sin(3.1416/180.0 * 0.5*(outLats[b2+1]+outLats[b2]) );
			
			for(int j=NY-b2;j<NY;j++){

				IUBAR(i,j,k) = (vort - 2.0*fc*sin(3.1416/180.0 * 0.5*(outLats[j]+outLats[j-1]))) * dy + IUBAR(i,j-1,k);
				ITHBAR(i,j,k) = ITHBAR(i,NY-b2,k);
				IQBAR(i,j,k) = IQBAR(i,NY-b2,k);
			}
#endif			

#if 0
			
			if(b3[k]<-1){
						
				for(int j=0;j<b3[k];j++){

					IUBAR(i,j,k) = IUBAR(i,b3[k],k);
					//ITHBAR(i,j,k) = ITHBAR(i,b3[k],k);
					IQBAR(i,j,k) = IQBAR(i,b3[k],k);
				}
			}
			
			if(b[k]>-1){
						
				for(int j=b[k];j<NY;j++){

					IUBAR(i,j,k) = IUBAR(i,b[k],k);
					ITHBAR(i,j,k) = ITHBAR(i,b[k],k);
					IQBAR(i,j,k) = IQBAR(i,b[k],k);
				}
				
			} else {
				
				int js = 22;
				
				for(int j=js;j<NY;j++){
					
					if(IUBAR(i,j,k) > IUBAR(10,js,k)){
					
						IUBAR(i,j,k) = IUBAR(10,js,k);
						ITHBAR(i,j,k) = ITHBAR(10,js,k);
						IQBAR(i,j,k) = IQBAR(10,js,k);
					}
				}
			}
			
			
#endif
		}}
	}

	/***********************************************************************
	* For resting basic state set the wind components to zero and make
	* the temperature and moisture fields horizontally uniform
	************************************************************************/
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
	
	if(useMicrophysics){ init_microphysics();}

	/********************************************
	* Upper and lower boundary conditions for
	* base state array
	*********************************************/
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
	/***********************************************************************
	* 
	* INITIALIZE VERTICALLY VARYING BASIC STATE
	*
	************************************************************************/

	//-----------------------------------------------------
	// create base state potential temperature profile
	//-----------------------------------------------------
	tb[0] = ITHBAR(NX/2,NY/2,0);
	tbw[0] = 0.5*(ITHBAR(NX/2,NY/2,1)+ITHBAR(NX/2,NY/2,0));

	for(int k=1;k<NZ;k++){ tb[k] = ITHBAR(NX/2,NY/2,k);}

	for(int k=1;k<NZ;k++){ tbw[k] = 0.5*(tb[k]+tb[k-1]);}

	//-----------------------------------------------------
	// create base state specific humidity profile
	//-----------------------------------------------------
	qb[0] = IQBAR(NX/2,NY/2,0);

	for(int k=1;k<NZ;k++){ qb[k] = IQBAR(NX/2,NY/2,k);}

	//-----------------------------------------------------
	// create base state virtual potential temperature profile
	//-----------------------------------------------------
	for(int k=0;k<NZ;k++){ tbv[k]=tb[k]*(1.0+0.61*qb[k]);}

	//-----------------------------------------------------
	// create base state pressure profile
	//-----------------------------------------------------
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

	/*********************************************************************
	* Initialize topographic array
	**********************************************************************/
	init_topography();

	/*********************************************************************
	* Column averaged density
	**********************************************************************/
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	
		rhoavg2d[i][j] = 0;

		for(int k=htopo[i][j]+1;k<NZ-1;k++){ rhoavg2d[i][j] = rhoavg2d[i][j] + rhou[k]*tbv[k];}

		rhoavg2d[i][j] = rhoavg2d[i][j]/((double)(NZ-2-htopo[i][j]));
	}}

	/*********************************************************************
	* Initialize arrays for SOR
	**********************************************************************/
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		sor[0][i][j] = 0;
		sor[1][i][j] = 0;
		sor[2][i][j] = 0;
		
	}}

	/*********************************************************************
	* Remove model base state from environmental base state
	**********************************************************************/
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				ITHBAR(i,j,k) = ITHBAR(i,j,k) - tb[k];
				IQBAR(i,j,k) = IQBAR(i,j,k) - qb[k];
			}
		}
	}

	for(int i=NX/8;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
	
		IUBAR(i,j,k)  *= ((double)(NX-i) / (double)NX + 0.125);
		ITHBAR(i,j,k) *= ((double)(NX-i) / (double)NX + 0.125);
		IQBAR(i,j,k)  *= ((double)(NX-i) / (double)NX + 0.125);
	
	}}}

	/*********************************
	* Initialize friction array
	**********************************/
	init_friction();
#if 0
	fft_damp(NX,NY,NZ,0,&UBAR(0,0,0));
	fft_damp(NX,NY,NZ,0,&VBAR(0,0,0));
	fft_damp(NX,NY,NZ,0,&WBAR(0,0,0));
	fft_damp(NX,NY,NZ,0,&THBAR(0,0,0));
#endif
	/*********************************
	* Calculate base state vertical velocity
	**********************************/
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


	/*********************************************************************
	* Coriolis parameter
	**********************************************************************/
	for(int j=0;j<NY;j++){ 
		
		f[j] = 2*fc*sin(outLats[j]*trigpi/180.);//2*fc*sin(outLats[j]*trigpi/180.);
		dfdy[j] = 2*fc*cos(outLats[j]*trigpi/180.) * (1./meters_per_degree) * (trigpi/180.) ;
	}
	
	double smr,temp,pres,rh;
	
	for(int k=0;k<NZ;k++){
		
		temp = tb[k]*pib[k];
		pres = p0*pow(pib[k],(cp/Rd));
		
		smr = get_QV_Sat(temp,pres);
		
		//printf("%d %f %f %f\n",k,smr*1000,qb[k]*1000,qb[k]/smr);
		
	}
	
	
	
	#if !ENERGY
	//-----------------------------------------------------------------
	// Alter basic state humidity
	//-----------------------------------------------------------------
	initialize_landsea(landseaMaskFile);
	
	const double height_limit_low = 2000.;
	const double height_limit_high = 10000.;
	const double relative_humidity = 0.95;
	
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=1;k<NZ;k++){
		//-----------------------------------------------------------------
		// Limit region
		//-----------------------------------------------------------------
		//if(outLats[j] > 5 && outLats[j] < 24 && outLons[i] > 75 && outLons[i] < 99 && LANDSEA(i+3,j+3)>0.5){
		if(outLats[j] > 15 && outLats[j] < 25 && outLons[i] > 70 && outLons[i] < 110){
				
			temp = (ITHBAR(i,j,k)+tb[k]) * IPBAR(i,j,k);	// full, actual temperature
			pres = p0*pow(IPBAR(i,j,k),(cp/Rd));			// full, dimensional pressure
	
			smr = get_QV_Sat(temp,pres);					// calculate saturation mixing ratio
	
			rh = (IQBAR(i,j,k)+qb[k])/smr;					// calculate relative humidity
			
			// change humidity
			if(zsu[k] >= height_limit_low && zsu[k] <=  height_limit_high){ 
				
				//IQBAR(i,j,k) = relative_humidity * smr - qb[k];
			
			}
	
			//if(outLats[j] > 17 && outLats[j] < 17.5 && outLons[i] > 90 && outLons[i] < 91){
				//printf("%d %f %f %f %f %f %f %f %f %f\n",k,zsu[k],outLats[j],outLons[i],smr*1000,temp,pres,(IQBAR(i,j,k)+qb[k])*1000,rh,(IQBAR(i,j,k)+qb[k])/smr);
			//}
		}
	}}}
	
	free(landsea);
	#endif
}

/*********************************************************************
*
*
**********************************************************************/
void reinitialize(){

	mtime = 0;
	bigcounter = 0;

	/*********************************************************************
	* Initialize arrays for SOR
	**********************************************************************/
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		sor[0][i][j] = 0;
		sor[1][i][j] = 0;
		sor[2][i][j] = 0;
		
	}}

	/********************************************
	* Reset perturbation arrays
	*********************************************/
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
	
	if(USE_LINEAR_FRICTION && !MERIDIONAL_CROSS_SECTION){

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
				printf("%d %d %d %f %f\n",i,j,k,outLats[j],outLons[i]);
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

	for(int i=NX-1;i>=0;i--){
	for(int j=NY-1;j>=0;j--){

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
				if(zsu[k] < ITOPO(linear_lon,j,k) && ITOPO(linear_lon,j,k) > 0)	
			#endif
		#endif	
			{
				IISTOPO(i,j,k) = 0;
				IUISTOPO(i,j,k) = 0;
				IVISTOPO(i,j,k) = 0;
				IUISTOPO(i+1,j,k) = 0;
				IVISTOPO(i,j+1,k) = 0;
				height = k;

			} else {

				IISTOPO(i,j,k) = 1;
				IUISTOPO(i,j,k) = 1;
				IVISTOPO(i,j,k) = 1;
			}
		}

		htopo[i][j] = height;
	}}

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
	
	double zmid = 0.5*(zm+zt);
	
	for(int k=0;k<NZ;k++){
		
		if(zsu[k]<=zm){ t_z[k] = tmin * exp( -(zsu[k]*zsu[k])/(zm*zm/8.0) );		
		} else { t_z[k] = tmax * exp( -( (zsu[k]-zmid)*(zsu[k]-zmid) ) / ( (zt-zm)*(zt-zm) / 32.0) );}
		
		//printf("%f %f\n",zsu[k],t_z[k]);
	}
	
	double R,tmid;
	
	for(int r=0;r<NR;r++){
		
		R = dx*r;
		
		t_r[r] = exp( -R*R / (r_end*r_end/8.0) );// printf("%d %f\n",r,t_r[r]);
	}
	
	for(int k=NZ-2;k>0;k--){
	for(int r=0;r<NR;r++){
	
		tmid = 0.5* (t_z[k]*t_r[r]+t_z[k+1]*t_r[r]) / ((tb[k]+tb[k+1])*0.5);
	
		phi[rindex(r,k,NR)] = phi[rindex(r,k+1,NR)] - grav * tmid * (zsu[k+1]-zsu[k]);
		
	}
		//printf("%d %f\n",k,phi[rindex(0,k,NR)]);
	}
	
	double fcor = 2*7.292e-5*sin(23*trigpi/180.0);//5.0e-5;
	
	
	
	for(int k=0;k<NZ;k++){
		
		vr[rindex(0,k,NR)] = 0;
		
		for(int r=1;r<NR;r++){
		
			R = r*dx;
	
			vr[rindex(r,k,NR)] = -0.5*fcor*R + 0.5*sqrt( fcor*fcor*R*R + 4*R * 0.5*( phi[rindex(r+1,k,NR)]-phi[rindex(r-1,k,NR)] ) * one_d_dx );	
		}
		//printf("%f %f\n",10*dx,vr[rindex(10,k,NR)]);
	}
	
	int ind;
	double frac;
	
	//rinterp(xpos,ypos,xpos+5,ypos+5,&ind,&frac);
	
	printf("dist = %d %d %d %d\n",xpos,ypos,NR,NY);
	if(rank==0){
		
		memset(&IUBAR(0,0,0),0,NX*NY*NZ*sizeof(double));
		
		for(int i=xpos-NR;i<xpos+NR;i++){
		for(int j=ypos-NR;j<ypos+NR;j++){
		for(int k=1;k<NZ;k++){
		
			rinterp(xpos,ypos,i,j,&ind,&frac,0);
		
			if(ind < NR-2){ IUBAR(i,j,k) = t_z[k]*( (1.0-frac)*t_r[ind]+frac*t_r[ind+1]);}
		}}}
	}
	
	distributeArray(ths);distributeArray(thms);
	exchange(ths); exchange(thms);
	
	if(rank==0){
		
		memset(&IUBAR(0,0,0),0,NX*NY*NZ*sizeof(double));
		
		for(int i=xpos-NR;i<xpos+NR;i++){
		for(int j=ypos-NR;j<ypos+NR;j++){
		for(int k=1;k<NZ;k++){
		
			rinterp(xpos,ypos,i,j,&ind,&frac,0.0);
		
			if(ind < NR-2 && ind != 0){
				//cos(trigpi*(double)(i-xpos)/(double)(j-ypos))
				if(j<ypos){ IUBAR(i,j,k) = ((1.0-frac)*vr[rindex(ind,k,NR)]+frac*vr[rindex(ind+1,k,NR)]) * fabs( (double)(j-ypos)/((double)ind+frac)) ; }
				else { IUBAR(i,j,k) = -((1.0-frac)*vr[rindex(ind,k,NR)]+frac*vr[rindex(ind+1,k,NR)]) * fabs((double)(j-ypos)/((double)ind+frac));}
			
			//if(k==15 && i == xpos)
				//printf("%d %f\n",j,IUBAR(i,j,k));
			}
		}}}
		
		for(int i=1;i<NX;i++){
		for(int j=0;j<NY;j++){
		for(int k=0;k<NZ;k++){
		
			IVBAR(i,j,k) = 0.5*(IUBAR(i,j,k) + IUBAR(i-1,j,k));
			
		}}}
		
		for(int i=0;i<NX*NY*NZ;i++){ iubar[i] = ivbar[i];}
	}
	
	distributeArray(us);distributeArray(ums);
	exchange(us); exchange(ums);
	
	if(rank==0){
		
		memset(&IUBAR(0,0,0),0,NX*NY*NZ*sizeof(double));
		
		for(int i=xpos-NR;i<xpos+NR;i++){
		for(int j=ypos-NR;j<ypos+NR;j++){
		for(int k=1;k<NZ;k++){
		
			rinterp(xpos,ypos,i,j,&ind,&frac,0.0);
		
			if(ind < NR-2 && ind != 0){
				//cos(trigpi*(double)(i-xpos)/(double)(j-ypos))
				if(i<xpos){ IUBAR(i,j,k) = -((1.0-frac)*vr[rindex(ind,k,NR)]+frac*vr[rindex(ind+1,k,NR)]) * fabs((double)(i-xpos))/((double)ind+frac) ; }
				else { IUBAR(i,j,k) = ((1.0-frac)*vr[rindex(ind,k,NR)]+frac*vr[rindex(ind+1,k,NR)]) * fabs((double)(i-xpos))/((double)ind+frac);}
			}
			//if(k==15 && i == xpos)
				//printf("%d %f\n",j,IUBAR(i,j,k));
				//}
		}}}
		
		for(int i=0;i<NX;i++){
		for(int j=1;j<NY;j++){
		for(int k=0;k<NZ;k++){
		
			IVBAR(i,j,k) = 0.5*(IUBAR(i,j-1,k) + IUBAR(i,j,k));
			
		}}}
		
		for(int i=0;i<NX*NY*NZ;i++){ iubar[i] = ivbar[i];}
		
	}
	
	distributeArray(vs);distributeArray(vms);
	exchange(vs); exchange(vms);
	
	for(int i=0;i<fNX;i++){
	for(int j=0;j<fNY;j++){
	for(int k=0;k<fNZ;k++){
		
		U(i,j,k) = -U(i,j,k);
		UM(i,j,k) = U(i,j,k);
		V(i,j,k) = -V(i,j,k);
		VM(i,j,k) = V(i,j,k);
		TH(i,j,k) = -TH(i,j,k);
		THM(i,j,k) = TH(i,j,k);
		
	}}}
	
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
	
		for(int i=0;i<xdim*ydim*zdim;i++){ varIn[i] = 0;}
		for(int i=0;i<xdim*ydim*NZ;i++){ var_interpz[i] = 0;}
	
		get_data_at_time(filename,varname,time,varIn,xdim,ydim,zdim);
		
		vert_interpolate_1d_from_model(levs_out,zlevs, varIn, xdim,ydim,zdim,NZ,var_interpz);
		
		horz_interpolate_from_model(var_interpz,iubar,xdim,ydim,NZ,false,false,
		myDX/meters_per_degree,myDY/meters_per_degree,outLons[0]-myLons[0],outLats[0]-myLats[0]);	
	}
	//-----------------------------------------------------------------
	// Distribute results to other processes
	//-----------------------------------------------------------------
	distributeArray(var);
	
	exchange(var);
	
	for(int i=0;i<fNX*fNY*fNZ;i++){ 
		//var[i] *= 0.625;
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
			
			printf("Files have different dimensions or grid spacing.\nLoading initial conditions from file %s \nwith dimension %lu x %lu x %lu and grid spacing dx = %f \n",
				myfilename,xdim,ydim,zdim,grid_spacing[0]);
			
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
	

