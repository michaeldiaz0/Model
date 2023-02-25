#include "stdafx.h"
#include "fftw3.h"
#include "boundaries.h"
#include "pcomm.h"

/***************************************************************************
* ------------------------------ MACROS ------------------------------------
****************************************************************************/
#if !PARALLEL
	#define  IN2D(i,j)  in2d[(i)*NY+(j)]
	#define OUT2D(i,j) out2d[(i)*NY+(j)]

	#define  IN3D(i,j,k)  in3d[nj*((i)+ni*(k))+(j)]	// k,i,j array
	#define OUT3D(i,j,k) out3d[nj*((i)+ni*(k))+(j)]
#else
	#define IN2D(i,j)   in2d[(i)*myNY+(j)]
	#define OUT2D(i,j) out2d[(i)*myNY+(j)]

	#define IN3D(i,j,k)   in3d[NY*((i)+myNX*(k))+(j)]	// k,i,j array
	#define OUT3D(i,j,k) out3d[NY*((i)+myNX*(k))+(j)]
#endif

// serial version
#define B2(i,j,k) bcoef[nk*((j)+nj*(i))+(k)]
#define F2(i,j,k)   rhs[nk*((j)+nj*(i))+(k)]

// parallel version
#define B(i,j,k) bcoef[NZ*((j)+pNY*(i))+(k)]
#define F(i,j,k)   rhs[NZ*((j)+pNY*(i))+(k)]

// non-hydrostatic
#define  PRES_ROW3D(i,j,k)  pres_row[(i) * pNY*NZ  + NZ*(j) + k]
#define IPRES_ROW3D(i,j,k) ipres_row[(i) * pNY*NZ  + NZ*(j) + k]
#define  PRES_COL3D(i,j,k)  pres_col[(i) *  NY*NZ  + NZ*(j) + k]
#define IPRES_COL3D(i,j,k) ipres_col[(i) *  NY*NZ  + NZ*(j) + k]

// hydrostatic
#define  PRES_COL(i,j)  pres_col[(i)*NY  +(j)]
#define  PRES_ROW(i,j)  pres_row[(i)*myNY+(j)]
#define IPRES_COL(i,j) ipres_col[(i)*NY  +(j)]
#define IPRES_ROW(i,j) ipres_row[(i)*myNY+(j)]

#define BCOEF(i,j) bcoef[(i)*NY+(j)]

// zero pressure boundary condition
#define INR3D(i,j,k)   inreal3d[nj*((i)+ni*(k))+(j)]
#define OUTR3D(i,j,k) outreal3d[nj*((i)+ni*(k))+(j)]

/***************************************************************************
* -------------------------- VARIABLES -------------------------------------
****************************************************************************/
fftw_complex *in2d = NULL, *out2d = NULL;
fftw_complex *in3d = NULL, *out3d = NULL;
double *inreal3d = NULL, *outreal3d = NULL;
fftw_plan p1,p2,fyz,byz,fxz,bxz;

fftw_iodim ns1d[1];
fftw_iodim xz_loop[2];

fftw_iodim ew1d[1];
fftw_iodim yz_loop[2];

double *kwave,*lwave;
double *acoef,*bcoef,*ccoef;
double *P,*Q;
double *rhs;

double *pres_row,*ipres_row;
double *pres_col,*ipres_col;

double *div_avg;

/***************************************************************************
* ---------------------- FUNCTION PROTOTYPES--------------------------------
****************************************************************************/
void p_solve_hydrostatic_pressure(double);
void solve_hydrostatic_pressure(double);

void init_fft2d(int,int,int,int);
void init_fft3d(int,int,int);
void p_init_fft3d(int,int,int,int,int,int);
void init_fft3d_real(int,int,int,double,double,double);

void poisson_fft3d(int,int,int,double);
void poisson_fft3d_real(int,int,int);
void p_poisson_fft3d(double);
void poisson_fft2d(double *,int,int);
void hydrostatic_pressure(int,int,int,int);
void divergence3d(double);
void divergence3d_periodic(int,int,int,int,double);
//void SOR(int,int,double[3][NX][NY],double[NX][NY]);

//------------------------------------------------------------------------------------------------

/***************************************************************************
* Based on a set of preprocessor flags, decide how to initialize the
* anelastic pressure solver
****************************************************************************/
void initialize_pressure_solver(){

	//---------------------------------------------------------------
	if(PARALLEL){	// PARALLEL
		//-----------------------------------------------------------
		if(HYDROSTATIC){ // HYDROSTATIC
		//-----------------------------------------------------------
			init_fft2d(NX,myNY,NX,NY);
		//-----------------------------------------------------------
		} else {		// NON-HYDROSTATIC
		//-----------------------------------------------------------
			p_init_fft3d(NX,NY,NZ,myNX,pNY,pNZ);
		}
	//---------------------------------------------------------------
	} else {		// SERIAL
		//-----------------------------------------------------------
		if(HYDROSTATIC){ // HYDROSTATIC
		//-----------------------------------------------------------
			if(PERIODIC_BOUNDARIES)
                init_fft2d(NX-6,NY,NX-6,NY);
            else
                init_fft2d(NX,NY,NX,NY);
		//-----------------------------------------------------------
		} else {		// NON-HYDROSTATIC
		//-----------------------------------------------------------
			if(PERIODIC_BOUNDARIES){ 	// LINEAR
		//-----------------------------------------------------------
				init_fft3d(NX-6,NY,NZ);
		//-----------------------------------------------------------
			} else {		// NON-LINEAR
		//-----------------------------------------------------------
				init_fft3d(NX,NY,NZ);
			}
		}
	}
}

/***************************************************************************
* Higher level interface to solve the anelastic pressure equation. Uses 
* a set of preprocessor flags to decide how to compile it.
*
* @param step - fractional time step in Runge-Kutta loop
****************************************************************************/
void solve_pressure(double step){

	/***********************************************************************
	* ---------------------- HYDROSTATIC VERSION ---------------------------
	************************************************************************/
	#if HYDROSTATIC
		/*******************************************************************
		* ---------------------- SERIAL VERSION ---------------------------
		*******************************************************************/
		#if !PARALLEL
			if(!PERIODIC_BOUNDARIES){
                mirror_boundaries(&UP(0,0,0));
                mirror_boundaries(&VP(0,0,0));

                solve_hydrostatic_pressure(step);

                mirror_boundaries(&PI(0,0,0));
				
            } else {

				periodic_uvw_boundaries();

                solve_hydrostatic_pressure(step);

				periodic_pressure_boundaries();

            }
		/********************************************************************
		* ---------------------- PARALLEL VERSION ---------------------------
		*********************************************************************/
		#else
			exchange(ups);
			exchange(vps);

			p_solve_hydrostatic_pressure(step);

			exchange(pis);
		#endif
	/*************************************************************************
	* -------------------- NON-HYDROSTATIC VERSION --------------------------
	**************************************************************************/
	#else
		/*********************************************************************
		* ---------------------- SERIAL VERSION ------------------------------
		**********************************************************************/
		#if !PARALLEL
			/*****************************************************************
			* ------------------ NON-LINEAR VERSION --------------------------
			******************************************************************/
			if(!PERIODIC_BOUNDARIES){
				mirror_boundaries(&UP(0,0,0));	// Enforce zero velocity gradient at lateral boundaries.
				mirror_boundaries(&VP(0,0,0));	// Need to update because they are used to calculate
				mirror_boundaries(&WP(0,0,0));	// the divergence tendency required by pressure solver.
	
				poisson_fft3d(NX,NY,NZ,step);	// solves anelastic pressure equation
	
				mirror_boundaries(&PI(0,0,0));	// enforce zero pressure gradient at lateral boundary
			/*****************************************************************
			* ------------------ LINEAR VERSION ------------------------------
			******************************************************************/
			} else {
				periodic_uvw_boundaries();
				
				poisson_fft3d(NX-6,NY,NZ,step);	// solve anelastic pressure equation	
		
				periodic_pressure_boundaries();
			}
		/************************************************************************
		* ------------------ PARALLEL VERSION / NON-LINEAR ----------------------
		*************************************************************************/
		#else
			exchangeOnePoint(ups);	// Exchange process boundaries.
			exchangeOnePoint(vps);	// Need to update because they are used to calculate
			exchangeOnePoint(wps);	// the divergence tendency required by pressure solver.

			p_poisson_fft3d(step);	// solves anelastic pressure equation

			exchangeOnePoint(pis);	// exchange pressure at process boundaries
		#endif
		
	#endif

}

/**********************************************************************
* 
***********************************************************************/
void solve_hydrostatic_pressure(double step){

	int buf,buf2;
	
	if(PERIODIC_BOUNDARIES){
		buf = 3;
		buf2 = 3;
	} else {
		buf = 1;
		buf2 = 0;
	}

	/*****************************************************************
	* Calculate vertical density-weighted average of the divergence
	* of the new momentum field
	******************************************************************/
	double zsum = 0;
	
    for(int i=buf;i<NX-buf;i++){
    for(int j=1;j<NY-1;j++){
    
        div_avg[INDEX2D(i-buf2,j)] = 0;
        zsum = 0;
        
        // sum each vertical level
        for(int k=HTOPO(i,j)+1;k<NZ-1;k++){
            div_avg[INDEX2D(i-buf2,j)] += rhou[k] * DZU(k) * ( (UP(i+1,j,k)-UP(i,j,k))*one_d_dx + (VP(i,j+1,k)-VP(i,j,k) )*one_d_dy);
            zsum += DZU(k);
        }
            
        // vertical average
        div_avg[INDEX2D(i-buf2,j)] /= ( cp * zsum * step*dt );

    }}
    

    if(!PERIODIC_BOUNDARIES)
        mirror_boundaries_2d(div_avg);
	else
		mirror_boundaries_ns_2d(div_avg,0,NX-6);
  
	/*****************************************************************
	* Solve Laplacian equation for pressure field required 
	* to balance divergence tendencies from advection
	******************************************************************/
    poisson_fft2d(div_avg,NX-2*buf2,NY);

	/*********************************************************************
	* Calculate hydrostatic pressure
	**********************************************************************/
	hydrostatic_pressure(1,NX-1,1,NY-1);

	/*****************************************************
	* Calcuate top lid pressure
	******************************************************/
    for(int i=buf;i<NX-buf;i++)
    for(int j=1;j<NY-1;j++)
        PI(i,j,NZ-1) = (div_avg[INDEX2D(i-buf2,j)] - PI(i,j,0)) / RHOAVG2D(i,j);

	/*****************************************************
	* Calcuate total perturbation pressure
	* total pressure = lid pressure + hydrostatic pressure
	******************************************************/
	for(int i=0;i<NX;i++)
	for(int j=0;j<NY;j++)
			for(int k=1;k<NZ-1;k++)
				PI(i,j,k) += PI(i,j,NZ-1);


}
#if PARALLEL
/**********************************************************************
* Parallel pressure solver for hydrostatic model
*
* @param step - fractional time step in Runge-Kutta loop
***********************************************************************/
void p_solve_hydrostatic_pressure(double step){
	
	/*****************************************************************
	* Calculate the vertical density-weighted average of the divergence
	* of the new momentum field. Have each process compute its portion
	* of the array.
	******************************************************************/
	double zsum = 0;
		
	for(int i=3;i<fNX-3;i++){
	for(int j=3;j<fNY-3;j++){
	
		PRES_ROW(i+ibs[rank]-3,j-3) = 0;
		zsum = 0;
		
		// sum each vertical level
		for(int k=HTOPO(i,j)+1;k<NZ-1;k++){
			PRES_ROW(i+ibs[rank]-3,j-3) += rhou[k] * DZU(k) * ( (UP(i+1,j,k)-UP(i,j,k))*one_d_dx + (VP(i,j+1,k)-VP(i,j,k) )*one_d_dy);
			zsum += DZU(k);
		}

		// vertical average
		PRES_ROW(i+ibs[rank]-3,j-3) = PRES_ROW(i+ibs[rank]-3,j-3) / ( cp * zsum * step*dt );
	}}
	
	/*****************************************************************
	* Each process broadcasts its portion to all other processes
	* in the same row.
	******************************************************************/
	pres_row_alltoall(&PRES_ROW(0,0));

	/*****************************************************************
	* Load entire row into input to Fourier transform routine and
	* execute Fourier transform.
	******************************************************************/
	for(int i=0;i<NX;i++){
	for(int j=0;j<myNY;j++){

		IN2D(i,j)[0] = PRES_ROW(i,j);
		IN2D(i,j)[1] = 0;
	}}

	fftw_execute(p1);

	/*****************************************************************
	* Load each process' portion of the row into its portion of the
	* column and normalize the result. Index of 0 is the real part
	* and Index of 1 is the imaginary part.
	******************************************************************/
	for(int i=0;i<myNX;i++){
	for(int j=0;j<myNY;j++){

		PRES_COL(i,j+jbs[rank]) = OUT2D(i+ibs[rank],j)[0] / NX;
		IPRES_COL(i,j+jbs[rank]) = OUT2D(i+ibs[rank],j)[1] / NX;
	}}

	/*****************************************************************
	* Each process broadcasts its portion to all other processes
	* in the same column.
	******************************************************************/
	pres_col_alltoall(&PRES_COL(0,0));
	pres_col_alltoall(&IPRES_COL(0,0));


	/*****************************************************************
	* Solve tridiagonal matrix.
	******************************************************************/
	double denom;

	for(int i=0;i<myNX;i++){

		// forward elimination
		for(int j=0;j<NY;j++){

			denom = BCOEF(i+ibs[rank],j) + acoef[j]*P[j];

			P[j+1] = -ccoef[0] / denom;
			Q[j+1] = (PRES_COL(i,j)-acoef[j]*Q[j]) / denom;

		}

		// backward substitution
		for(int j=NY-2;j>=0;j--){ PRES_COL(i,j) = P[j+1]*PRES_COL(i,j+1) + Q[j+1];}
	
		// forward elimination
		for(int j=0;j<NY;j++){

			denom = BCOEF(i+ibs[rank],j) + acoef[j]*P[j];

			P[j+1] = -ccoef[0] / denom;
			Q[j+1] = (IPRES_COL(i,j)-acoef[j]*Q[j]) / denom;

		}

		// backward substitution
		for(int j=NY-2;j>=0;j--){ IPRES_COL(i,j) = P[j+1]*IPRES_COL(i,j+1) + Q[j+1];}

	}

	for(int i=0;i<myNX;i++){
	for(int j=0;j<myNY;j++){

		PRES_ROW(i+ibs[rank],j) = PRES_COL(i,j+jbs[rank]);
		IPRES_ROW(i+ibs[rank],j) = IPRES_COL(i,j+jbs[rank]);
	}}

	pres_row_alltoall(&PRES_ROW(0,0));
	pres_row_alltoall(&IPRES_ROW(0,0));


	for(int i=0;i<NX;i++){
	for(int j=0;j<myNY;j++){

		OUT2D(i,j)[0] = PRES_ROW(i,j);
		OUT2D(i,j)[1] = IPRES_ROW(i,j);
	}}


	fftw_execute(p2);

	/*********************************************************************
	* Calculate hydrostatic pressure
	**********************************************************************/
	hydrostatic_pressure(3,fNX-3,3,fNY-3);

	/*****************************************************
	* Calcuate top lid pressure
	******************************************************/
	for(int i=3;i<fNX-3;i++){
	for(int j=3;j<fNY-3;j++){
			PI(i,j,NZ-1) = (IN2D(i+ibs[rank]-3,j-3)[0] - PI(i,j,0)) / RHOAVG2D(i,j);

			//printf("%d %d %f\n",i+ibs[rank]-3,j,IN2D(i+ibs[rank]-3,j)[0]);
	}}
	/*****************************************************
	* Calcuate total perturbation pressure
	* total pressure = lid pressure + hydrostatic pressure
	******************************************************/
	for(int i=3;i<fNX-3;i++)
	for(int j=3;j<fNY-3;j++)
			for(int k=1;k<NZ-1;k++)
				PI(i,j,k) += PI(i,j,NZ-1);
	
}
#endif
/*********************************************************************
* Calculate hydrostatic pressure
**********************************************************************/
void hydrostatic_pressure(int il,int ih,int jl,int jh){

	double tup; double tdn;

	if(USE_MICROPHYSICS){

        if(USE_ICE){
            for(int i=il;i<ih;i++){
            for(int j=jl;j<jh;j++){

                PI(i,j,NZ-1) = 0;
                PI(i,j,NZ-2) = -0.5*(grav/cp)*((TH(i,j,NZ-2)*
                (1.+0.61*QV(i,j,NZ-2)-QC(i,j,NZ-2)-QR(i,j,NZ-2)-QS(i,j,NZ-2)-QI(i,j,NZ-2)))/(tbv[NZ-2]*tbv[NZ-2]))*DZW(NZ-1);

                for(int k=NZ-3;k>0;k--){
                
                    tup = (TH(i,j,k+1)*(1.+0.61*QV(i,j,k+1)-QC(i,j,k+1)-QR(i,j,k+1)-QS(i,j,k+1)-QI(i,j,k+1)))/(tbv[k+1]*tbv[k+1]);
                    tdn = (TH(i,j,k  )*(1.+0.61*QV(i,j,k  )-QC(i,j,k  )-QR(i,j,k  )-QS(i,j,k  )-QI(i,j,k  )))/(tbv[k  ]*tbv[k  ]);
                    PI(i,j,k) = PI(i,j,k+1)-0.5*(grav/cp)*(tup+tdn)*DZW(k+1);
                }
            }}

        } else {

            for(int i=il;i<ih;i++){
            for(int j=jl;j<jh;j++){

                PI(i,j,NZ-1) = 0;
                PI(i,j,NZ-2) = -0.5*(grav/cp)*((TH(i,j,NZ-2)*(1.+0.61*QV(i,j,NZ-2)-QC(i,j,NZ-2)-QR(i,j,NZ-2)))/(tbv[NZ-2]*tbv[NZ-2]))*DZW(NZ-1);

                for(int k=NZ-3;k>0;k--){
                
                    tup = (TH(i,j,k+1)*(1.+0.61*QV(i,j,k+1)-QC(i,j,k+1)-QR(i,j,k+1)))/(tbv[k+1]*tbv[k+1]);
                    tdn = (TH(i,j,k)*(1.+0.61*QV(i,j,k)-QC(i,j,k)-QR(i,j,k)))/(tbv[k]*tbv[k]);
                    PI(i,j,k) = PI(i,j,k+1)-0.5*(grav/cp)*(tup+tdn)*DZW(k+1);
                }
            
            }}
        }


	} else {

		for(int i=il;i<ih;i++){
		for(int j=jl;j<jh;j++){

			PI(i,j,NZ-1) = 0;
			//PI(i,j,NZ-2) = -0.5*(grav/cp)*((THP(i,j,NZ-2)*(1.+0.61*QVP(i,j,NZ-2)))/(tbv[NZ-2]*tbv[NZ-2]))*DZW(NZ-1);
			PI(i,j,NZ-2) = -0.5*(grav/cp)*((TH(i,j,NZ-2)*(1))/(tbv[NZ-2]*tbv[NZ-2]))*DZW(NZ-1);

			for(int k=NZ-3;k>0;k--){
			
				tup = (TH(i,j,k+1)*(1.))/(tbv[k+1]*tbv[k+1]);
				tdn = (TH(i,j,k  )*(1.))/(tbv[k]*tbv[k]);
				//tup = (TH(i,j,k+1)*(1.+0.61*QVP(i,j,k+1)))/(tbv[k+1]*tbv[k+1]);
				//tdn = (TH(i,j,k)*(1.+0.61*QVP(i,j,k)))/(tbv[k]*tbv[k]);
				PI(i,j,k) = PI(i,j,k+1)-0.5*(grav/cp)*(tup+tdn)*DZW(k+1);
			}
		
		}}
	}

	/*****************************************************
	*	Calculate column average hydrostatic pressure
	******************************************************/
	double zsum = 0;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	
		PI(i,j,0) = 0;
		zsum = 0;
	
		for(int k=HTOPO(i,j)+1;k<NZ-1;k++){
			
			PI(i,j,0) += rhou[k]*tbv[k]*PI(i,j,k)*DZU(k);
			zsum += DZU(k);
		}
			
		PI(i,j,0) /= zsum;

	}}

}
#if 0
/**********************************************************************
* Succesive over-relaxation algorithm
***********************************************************************/
void SOR(int itr,int iterations,double p[3][NX][NY],double div_avg[NX][NY]){

	/**********************************************
	* Take off average value for the lid pressure 
	* at the boundaries to constrain the potential values
	***********************************************/
	double sum = 0;

	for(int i=0;i<NX;i++){ sum = sum + p[itr][i][1] + p[itr][i][NY-2];}
	for(int j=1;j<NY-1;j++){ sum = sum + p[itr][1][j] + p[itr][NX-2][j];}

	sum = sum / (double)(2*NX+2*NY);

	// boundary condition
	for(int i=0;i<NX;i++){ p[itr][i][0] = p[itr][i][1]-sum; p[itr][i][NY-1] = p[itr][i][NY-2]-sum;}
	for(int j=0;j<NY;j++){ p[itr][0][j] = p[itr][1][j]-sum; p[itr][NX-1][j] = p[itr][1][j]-sum;}


	double omega = 2 / (1+trigpi/(double)NX);		// over-relaxation coefficient
	double dxdy = dx*dy;
	double omega_div_4 = 0.25*omega;
	double oneminusomega = 1-omega;

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		div_avg[i][j] =  dxdy*div_avg[i][j];

	}}

	for(int q=1;q<iterations;q++){

		for(int i=1;i<NX-1;i++){
		for(int j=1;j<NY-1;j++){
	
			p[itr][i][j] = oneminusomega*p[itr][i][j] + 

				omega_div_4*(   p[itr][i-1][j] + p[itr][i+1][j]+ p[itr][i][j-1] + p[itr][i][j+1]  - div_avg[i][j]  );
		}}
	}

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		div_avg[i][j] =  p[itr][i][j];

	}}

}
#endif
/**********************************************************************
* Fast Fourier transform-based pressure solver for hydrostatic version
*
* ni,nj - grid dimensions for full domain array
***********************************************************************/
/**********************************************************************
* Fast Fourier transform-based pressure solver for hydrostatic version
*
* ni,nj - grid dimensions for full domain array
***********************************************************************/
void init_fft2d(int ni,int nj,int nx,int ny){
	
	div_avg = (double*) calloc(sizeof(double),ni*nj);
	
	in2d = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nj);
	out2d = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nj);

	int n[] = {nx};

	p1 = fftw_plan_many_dft(1, n, nj,in2d,n,nj,1,out2d,n,nj,1,FFTW_FORWARD, FFTW_MEASURE);
	p2 = fftw_plan_many_dft(1, n, nj,out2d,n,nj,1,in2d,n,nj,1,FFTW_BACKWARD,FFTW_MEASURE);

	kwave = (double*) calloc(sizeof(double),nx);
	lwave = (double*) calloc(sizeof(double),ny);

	for(int i=0;i<nx/2;i++){	kwave[i] = 2*trigpi * (double)i / nx;	}
	for(int i=nx/2;i<nx;i++){	kwave[i] = 2*trigpi * (double)(nx-i) / nx;}

	for(int i=0;i<ny/2;i++){	lwave[i] = 2*trigpi * (double)i / ny;	}
	for(int i=ny/2;i<ny;i++){	lwave[i] = 2*trigpi * (double)(ny-i) / ny;}

	acoef = (double*) malloc(sizeof(double)*ny);
	bcoef = (double*) malloc(sizeof(double)*nx*ny);
	ccoef = (double*) malloc(sizeof(double)*ny);

	P = (double*) calloc(sizeof(double),(ny+1));
	Q = (double*) calloc(sizeof(double),(ny+1));

	for(int j=0;j<ny;j++){

		acoef[j] = 1./(dy*dy);
		ccoef[j] = 1./(dy*dy);

	}

	acoef[0] = 0;	
	ccoef[ny-1] = 0;

	for(int i=0;i<nx;i++){
	for(int j=0;j<ny;j++){

		BCOEF(i,j) = 2.*(cos(kwave[i])-1.)/(dx*dx) - ccoef[j] - acoef[j];
		
		if(j==0){ BCOEF(i,j) = -rhow[2]/(dx*dx);}
	}}

}

/**********************************************************************
* Fast Fourier transform-based pressure solver
***********************************************************************/
void poisson_fft2d(double *Fxy,int ni, int nj){

	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){

		IN2D(i,j)[0] = Fxy[INDEX2D(i,j)];
		IN2D(i,j)[1] = 0;

	}}

	fftw_execute(p1);

	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){

		OUT2D(i,j)[0] = OUT2D(i,j)[0] / ni;
		OUT2D(i,j)[1] = OUT2D(i,j)[1] / ni;

	}}

	/*****************************************************************
	* Solve tridiagonal matrix.
	******************************************************************/
	double denom;

	for(int i=0;i<ni;i++){

		// forward elimination
		for(int j=0;j<nj;j++){

			denom = BCOEF(i,j) + acoef[j]*P[j];

			P[j+1] = -ccoef[j] / denom;
			Q[j+1] = (OUT2D(i,j)[0]-acoef[j]*Q[j]) / denom;

		}

		// backward substitution
		for(int j=nj-2;j>=0;j--){ OUT2D(i,j)[0] = P[j+1]*OUT2D(i,j+1)[0] + Q[j+1];}

		// forward elimination
		for(int j=0;j<nj;j++){

			denom = BCOEF(i,j) + acoef[j]*P[j];

			P[j+1] = -ccoef[j] / denom;
			Q[j+1] = (OUT2D(i,j)[1]-acoef[j]*Q[j]) / denom;

		}

		// backward substitution
		for(int j=nj-2;j>=0;j--){ OUT2D(i,j)[1] = P[j+1]*OUT2D(i,j+1)[1] + Q[j+1];}

	}

	fftw_execute(p2);
	

	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){

		Fxy[INDEX2D(i,j)] = IN2D(i,j)[0];

	}}

}

/**********************************************************************
* Fast Fourier transform-based pressure solver for serial version
*
* ni,nj,nk - grid dimensions for full domain array
***********************************************************************/
void init_fft3d(int ni,int nj,int nk){

	kwave = (double*) malloc(sizeof(double)*ni);
	lwave = (double*) malloc(sizeof(double)*nj);

	acoef = (double*) malloc(sizeof(double)*nk);
	bcoef = (double*) malloc(sizeof(double)*ni*nj*nk);
	ccoef = (double*) malloc(sizeof(double)*nk);

	in3d = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ni*nj*nk);
	out3d = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ni*nj*nk);
	rhs = (double*) malloc(sizeof(double)*ni*nj*nk);

	for(int i=0;i<ni*nj*nk;i++){ 
		
		in3d[i][0] = 0;
		in3d[i][1] = 0;
	}

	int n[] = {ni,nj};
	int *inembed = n;
	int *onembed = n;

	p1 = fftw_plan_many_dft(2, n, nk, in3d,inembed, 1,ni*nj,out3d,onembed, 1,ni*nj,FFTW_FORWARD, FFTW_MEASURE);
	p2 = fftw_plan_many_dft(2, n, nk,out3d,inembed, 1,ni*nj,in3d, onembed, 1,ni*nj,FFTW_BACKWARD, FFTW_MEASURE);


	for(int i=0;i<ni;i++){  kwave[i] = (double)i * 2*trigpi/((double)ni);}
	for(int i=0;i<nj;i++){  lwave[i] = (double)i * 2*trigpi/((double)nj);}


	if(!STRETCHED_GRID){

		for(int k=1;k<nk-1;k++){ acoef[k] = rhow[k]/(dz*dz); ccoef[k] = rhow[k+1]/(dz*dz);}

		acoef[1] = 0;
		ccoef[nk-2] = 0;

		for(int i=0;i<ni;i++){
		for(int j=0;j<nj;j++){

			if(i==0 && j==0){ B(i,j,1) = -rhow[2]/(dz*dz);}
			else {
				B2(i,j,1) = rhou[1]*2.0*(cos(kwave[i])+cos(lwave[j])-2.0)/(dx*dx) - rhow[2]/(dz*dz);
			}

			B2(i,j,nk-2) = rhou[nk-2]*2.0*(cos(kwave[i])+cos(lwave[j])-2.0)/(dx*dx) - (rhow[nk-2])/(dz*dz);
		
		for(int k=2;k<nk-2;k++){

			B2(i,j,k) = rhou[k]*2.0*(cos(kwave[i])+cos(lwave[j])-2.0)/(dx*dx) - (rhow[k+1]+rhow[k])/(dz*dz);
		}}}
		
	} else {

		for(int k=1;k<nk-1;k++){ acoef[k] = mu[k]*mw[k]*rhow[k]/(dz*dz); ccoef[k] = mu[k]*mw[k+1]*rhow[k+1]/(dz*dz);}

		acoef[1] = 0;
		ccoef[nk-2] = 0;

		for(int i=0;i<ni;i++){
		for(int j=0;j<nj;j++){

			B2(i,j,1) = rhou[1]*2.0*(cos(kwave[i])+cos(lwave[j])-2.0)/(dx*dx) - mu[1]*mw[2]*rhow[2]/(dz*dz);

			for(int k=2;k<nk-2;k++){ B2(i,j,k) = rhou[k]*2.0*(cos(kwave[i])+cos(lwave[j])-2.0)/(dx*dx) - mu[k]*(mw[k+1]*rhow[k+1]+mw[k]*rhow[k])/(dz*dz);}
		
			B2(i,j,nk-2) = rhou[nk-2]*2.0*(cos(kwave[i])+cos(lwave[j])-2.0)/(dx*dx) - (mu[nk-2]*mw[nk-2]*rhow[nk-2])/(dz*dz);
		}}
		
	}


}

/**********************************************************************
* Solve the anelastic pressure equation using Fast Fourier transforms
* in the horizontal and Gaussian elimination in the vertical
*
* @param ni,nj,nk - grid dimensions
***********************************************************************/
void poisson_fft3d(int ni,int nj,int nk,double step){

	int buffer;

	/****************************************************
	* Set input to Fourier transform subroutine
	*****************************************************/
	if(!PERIODIC_BOUNDARIES){
		buffer = 0;
		divergence3d(step);
	} else {
		buffer = 3;
		divergence3d_periodic(3,NX-3,1,NY-1,step);
	}
	/****************************************************
	* Execute Fourier transform
	*****************************************************/
	fftw_execute(p1); 

	double size2d = (double)(ni*nj);

	// normalize
	for(int i=0;i<ni*nj*nk;i++){ out3d[i][0] /= size2d; out3d[i][1] /= size2d;}

	/****************************************************
	* Solve tridiagonal matrix
	*****************************************************/
	double z;

	/**********************************************
	* Forward elimination
	***********************************************/
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){

		z = 1.0 / B2(i,j,1);

		F2(i,j,1) = ccoef[1]*z;

		OUT3D(i,j,1)[0] = OUT3D(i,j,1)[0]*z;
		OUT3D(i,j,1)[1] = OUT3D(i,j,1)[1]*z;

		for(int k=2;k<nk-1;k++){

			z = 1.0 / ( B2(i,j,k) - acoef[k]*F2(i,j,k-1));
		
			F2(i,j,k) = ccoef[k]*z;

			OUT3D(i,j,k)[0] = (OUT3D(i,j,k)[0] - acoef[k]*OUT3D(i,j,k-1)[0])*z;
			OUT3D(i,j,k)[1] = (OUT3D(i,j,k)[1] - acoef[k]*OUT3D(i,j,k-1)[1])*z;
		}
	}}
	
	OUT3D(0,0,NZ-2)[0] = 0;
	OUT3D(0,0,NZ-2)[1] = 0;

	/***********************************************
	* Back substitution
	************************************************/
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=NZ-3;k>=1;k--){

		OUT3D(i,j,k)[0] -= F2(i,j,k)*OUT3D(i,j,k+1)[0];
		OUT3D(i,j,k)[1] -= F2(i,j,k)*OUT3D(i,j,k+1)[1];

	}}}

	/****************************************************
	* Execute reverse Fourier transform
	*****************************************************/
	fftw_execute(p2);

	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){

		PI(i+buffer,j,0) = IN3D(i,j,1)[0];
		IN3D(i,j,0)[1] = 0;

		for(int k=1;k<nk-1;k++){

			PI(i+buffer,j,k) = IN3D(i,j,k)[0];
			IN3D(i,j,k)[1] = 0;
		}

		PI(i+buffer,j,NZ-1) = IN3D(i,j,NZ-2)[0];
		IN3D(i,j,NZ-1)[1] = 0;
	}}
}
#if PARALLEL
/**********************************************************************
* Set up fast Fourier transform-based pressure solver for parallel version
*
* ni,nj,nk - grid dimensions for full domain array
* snx,sny,snz - grid dimensions for pencil shaped subarrays
***********************************************************************/
void p_init_fft3d(int ni,int nj,int nk,int snx,int sny,int snz){

	rhs = (double*) malloc(sizeof(double)*NX*sny*nk);

	/****************************************************
	* Set up Fourier transform stuff
	*****************************************************/
	
	//----------------------------------------------------------
	// array description for columns
	//----------------------------------------------------------
	ns1d[0].n = nj;		ns1d[0].is = nk;	ns1d[0].os = nk;
	
	xz_loop[0].n = snx;		xz_loop[0].is = nk*nj;	xz_loop[0].os = nk*nj;
	xz_loop[1].n = snz;		xz_loop[1].is = 1;		xz_loop[1].os = 1;
	
	//----------------------------------------------------------
	// array description for rows
	//----------------------------------------------------------
	ew1d[0].n = ni;		ew1d[0].is = nk*sny;	ew1d[0].os = nk*sny;
	
	yz_loop[0].n = sny;		yz_loop[0].is = nk;		yz_loop[0].os = nk;
	yz_loop[1].n = snz;		yz_loop[1].is = 1;		yz_loop[1].os = 1;
	
	//----------------------------------------------------------
	// column (north-south) Fourier transforms	
	//----------------------------------------------------------
	fyz = fftw_plan_guru_split_dft(	
			1, ns1d,2, xz_loop, 
			&PRES_COL3D(0,0,kbs[col_rank]),&IPRES_COL3D(0,0,kbs[col_rank]),
			&PRES_COL3D(0,0,kbs[col_rank]),&IPRES_COL3D(0,0,kbs[col_rank]),
			FFTW_FLAGS
	);
	
	byz = fftw_plan_guru_split_dft(
			1, ns1d,2, xz_loop, 
			&IPRES_COL3D(0,0,kbs[col_rank]),&PRES_COL3D(0,0,kbs[col_rank]),
			&IPRES_COL3D(0,0,kbs[col_rank]),&PRES_COL3D(0,0,kbs[col_rank]),
			FFTW_FLAGS
	);
			
	//----------------------------------------------------------
	// row (east-west) Fourier transforms
	//----------------------------------------------------------
	fxz = fftw_plan_guru_split_dft(
			1, ew1d,2, yz_loop,
			&PRES_ROW3D(0,0,kbs[col_rank]),&IPRES_ROW3D(0,0,kbs[col_rank]),
			&PRES_ROW3D(0,0,kbs[col_rank]),&IPRES_ROW3D(0,0,kbs[col_rank]),
			FFTW_FLAGS
	);
			
	bxz = fftw_plan_guru_split_dft(
			1, ew1d,2, yz_loop, 
			&IPRES_ROW3D(0,0,kbs[col_rank]),&PRES_ROW3D(0,0,kbs[col_rank]),
			&IPRES_ROW3D(0,0,kbs[col_rank]),&PRES_ROW3D(0,0,kbs[col_rank]),
			FFTW_FLAGS
	);
	
	/****************************************************
	* Set up matrix eigenvalues
	*****************************************************/
	kwave = (double*) malloc(sizeof(double)*ni);
	lwave = (double*) malloc(sizeof(double)*nj);

	acoef = (double*) malloc(sizeof(double)*nk);
	bcoef = (double*) malloc(sizeof(double)*ni*sny*nk);
	ccoef = (double*) malloc(sizeof(double)*nk);
	
	for(int i=0;i<ni;i++){  kwave[i] = (double)i * 2*trigpi/((double)ni);}
	for(int i=0;i<nj;i++){  lwave[i] = (double)i * 2*trigpi/((double)nj);}


	if(!STRETCHED_GRID){

		for(int k=1;k<nk-1;k++){ acoef[k] = rhow[k]/(dz*dz); ccoef[k] = rhow[k+1]/(dz*dz);}

		acoef[1] = 0;
		ccoef[nk-2] = 0;

		for(int i=0;i<ni;i++){
		for(int j=0;j<sny;j++){

			B(i,j,1) = rhou[1]*2.0*(cos(kwave[i])+cos(lwave[j+jbs_p[row_rank]])-2.0)/(dx*dx) - rhow[2]/(dz*dz);

			for(int k=2;k<nk-2;k++){ B(i,j,k) = rhou[k]*2.0*(cos(kwave[i])+cos(lwave[j+jbs_p[row_rank]])-2.0)/(dx*dx) - (rhow[k+1]+rhow[k])/(dz*dz);}
		
			B(i,j,nk-2) = rhou[nk-2]*2.0*(cos(kwave[i])+cos(lwave[j+jbs_p[row_rank]])-2.0)/(dx*dx) - (rhow[nk-2])/(dz*dz);
		}}
	
		
	} else {

		for(int k=1;k<nk-1;k++){ acoef[k] = mu[k]*mw[k]*rhow[k]/(dz*dz); ccoef[k] = mu[k]*mw[k+1]*rhow[k+1]/(dz*dz);}

		acoef[1] = 0;
		ccoef[nk-2] = 0;

		for(int i=0;i<ni;i++){
		for(int j=0;j<sny;j++){

			B(i,j,1) = rhou[1]*2.0*(cos(kwave[i])+cos(lwave[j+jbs_p[row_rank]])-2.0)/(dx*dx) - mu[1]*mw[2]*rhow[2]/(dz*dz);

			for(int k=2;k<nk-2;k++){ B(i,j,k) = rhou[k]*2.0*(cos(kwave[i])+cos(lwave[j+jbs_p[row_rank]])-2.0)/(dx*dx) - mu[k]*(mw[k+1]*rhow[k+1]+mw[k]*rhow[k])/(dz*dz);}
		
			B(i,j,nk-2) = rhou[nk-2]*2.0*(cos(kwave[i])+cos(lwave[j+jbs_p[row_rank]])-2.0)/(dx*dx) - (mu[nk-2]*mw[nk-2]*rhow[nk-2])/(dz*dz);
		}}
		
	}
	
}

/**********************************************************************
* This is the parallel version.
* Solve the anelastic pressure equation using Fast Fourier transforms
* in the horizontal and Guassian elimination in the vertical.
*
* @param step - fraction of full time step
***********************************************************************/
void p_poisson_fft3d(double step){

	/****************************************************
	* Calculate divergence for the process' full vertical 
	* column and load it into its respective location 
	* within the north-south column 3947
	*****************************************************/
	for(int i=3;i<fNX-3;i++){
	for(int j=3;j<fNY-3;j++){
	for(int k=1;k<NZ-1;k++){
		
		PRES_COL3D(i-3,j+jbs[rank]-3,k) =
			rhou[k  ]*(UP(i+1,j,k)-		   UP(i,j,k) )*one_d_dx + 
			rhou[k  ]*(VP(i,j+1,k)-		   VP(i,j,k) )*one_d_dy + 
		  ( rhow[k+1]* WP(i,j,k+1)-rhow[k]*WP(i,j,k) )*ONE_D_DZ(k);
		
		PRES_COL3D(i-3,j+jbs[rank]-3,k) /= (step*dt);
	}}}

	/**************************************************************
	* --------------- FORWARD FOURIER TRANSFORM -------------------
	***************************************************************/

	/****************************************************
	* All processes within a column exchange their data
	*****************************************************/
	pres_col_alltoall(&PRES_COL3D(0,0,0));
	
	/****************************************************
	* Initialize imaginary part of divergence to zero
	*****************************************************/
	for(int i=0;i<myNX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<pNZ;k++){
		
		IPRES_COL3D(i,j,k+kbs[col_rank]) = 0;
	}}}

	fftw_execute(fyz);

	/****************************************************
	* All processes within a row exchange their data
	*****************************************************/
	pres_row_alltoall2(&PRES_COL3D(0,0,0) ,&PRES_ROW3D(0,0,0));	
	pres_row_alltoall2(&IPRES_COL3D(0,0,0),&IPRES_ROW3D(0,0,0));
	
	fftw_execute(fxz);
	
	double dsize = (double)(NX*NY);
	
	for(int i=0;i<NX;i++){
	for(int j=0;j<pNY;j++){
	for(int k=kbs[col_rank];k<kbs[col_rank]+pNZ;k++){
		
	 	PRES_ROW3D(i,j,k) /= dsize;
		IPRES_ROW3D(i,j,k) /= dsize;
	}}}

	/**************************************************************
	* ---------------- SOLVE TRIDIAGONAL MATRIX -------------------
	***************************************************************/
	double z;
	
	pres_vert_alltoall(&PRES_ROW3D(0,0,0));
	pres_vert_alltoall(&IPRES_ROW3D(0,0,0));

	/**********************************************
	* Forward elimination
	***********************************************/
	for(int i=ibs_p[col_rank];i<ibs_p[col_rank]+pNX;i++){
	for(int j=0;j<pNY;j++){

		z = 1.0 / B(i,j,1);

		F(i,j,1) = ccoef[1]*z;

		 PRES_ROW3D(i,j,1) *= z;
		IPRES_ROW3D(i,j,1) *= z;

		for(int k=2;k<NZ-1;k++){

			z = 1.0 / ( B(i,j,k) - acoef[k]*F(i,j,k-1));
		
			F(i,j,k) = ccoef[k]*z;

			 PRES_ROW3D(i,j,k) = ( PRES_ROW3D(i,j,k) - acoef[k]* PRES_ROW3D(i,j,k-1))*z;
			IPRES_ROW3D(i,j,k) = (IPRES_ROW3D(i,j,k) - acoef[k]*IPRES_ROW3D(i,j,k-1))*z;
		}
	}}

	if(rank==0){
	
		PRES_ROW3D(0,0,NZ-2) = 0;
		IPRES_ROW3D(0,0,NZ-2) = 0;
	}

	/***********************************************
	* Back substitution
	************************************************/
	for(int i=ibs_p[col_rank];i<ibs_p[col_rank]+pNX;i++){
	for(int j=0;j<pNY;j++){
	for(int k=NZ-3;k>=1;k--){

		 PRES_ROW3D(i,j,k) -= F(i,j,k)* PRES_ROW3D(i,j,k+1);
		IPRES_ROW3D(i,j,k) -= F(i,j,k)*IPRES_ROW3D(i,j,k+1);
	}}}

	pres_vert_alltoall_reverse(&PRES_ROW3D(0,0,0));
	pres_vert_alltoall_reverse(&IPRES_ROW3D(0,0,0));

	/**************************************************************
	* --------------- BACKWARD FOURIER TRANSFORM ------------------
	***************************************************************/
	fftw_execute(bxz);
		
	pres_row_alltoall2_reverse(&PRES_COL3D(0,0,0) ,&PRES_ROW3D(0,0,0));
	pres_row_alltoall2_reverse(&IPRES_COL3D(0,0,0),&IPRES_ROW3D(0,0,0));

	fftw_execute(byz);

	/****************************************************
	* All processes within a column exchange their data
	*****************************************************/
	pres_col_alltoall_reverse(&PRES_COL3D(0,0,0));
	
	/****************************************************
	* The result is now the anelastic pressure field.
	*****************************************************/
	for(int i=3;i<fNX-3;i++){
	for(int j=3;j<fNY-3;j++){

		PI(i,j,0) = PRES_COL3D(i-3,j-3+jbs[rank],1);	// zero gradient lower boundary

		for(int k=1;k<NZ-1;k++){ PI(i,j,k) = PRES_COL3D(i-3,j-3+jbs[rank],k);}

		PI(i,j,NZ-1) = PRES_COL3D(i-3,j-3+jbs[rank],NZ-2);	// zero gradient upper boundary
	}}

}
#endif
/**********************************************************************
* Fast Fourier transform-based pressure solver
***********************************************************************/
void init_fft3d_real(int ni,int nj,int nk,double dx,double dy,double dz){

	kwave = (double*) malloc(sizeof(double)*ni);
	lwave = (double*) malloc(sizeof(double)*nj);

	acoef = (double*) malloc(sizeof(double)*nk);
	bcoef = (double*) malloc(sizeof(double)*ni*nj*nk);
	ccoef = (double*) malloc(sizeof(double)*nk);

	P = (double*) calloc(sizeof(double),(nk+1));
	Q = (double*) calloc(sizeof(double),(nk+1));

	inreal3d = (double*) calloc(sizeof(double),ni*nj*nk);
	outreal3d = (double*) calloc(sizeof(double),ni*nj*nk);
	rhs = (double*) malloc(sizeof(double)*ni*nj*nk);

	int n[] = {ni,nj};
	int *inembed = n;
	int *onembed = n;
	fftw_r2r_kind kind1[] = {FFTW_RODFT00,FFTW_RODFT00};
	fftw_r2r_kind kind2[] = {FFTW_RODFT00,FFTW_RODFT00};

	p1 = fftw_plan_many_r2r(2, n, nk, inreal3d,inembed, 1,ni*nj,outreal3d,onembed, 1,ni*nj,kind1, FFTW_MEASURE);
	p2 = fftw_plan_many_r2r(2, n, nk,outreal3d,inembed, 1,ni*nj,inreal3d, onembed, 1,ni*nj,kind2, FFTW_MEASURE);

	double shift = 1;

	for(int i=0;i<ni;i++){  kwave[i] = ((double)(i)+shift) * trigpi/((double)(ni+1));}

	for(int i=0;i<nj;i++){  lwave[i] = ((double)(i)+shift) * trigpi/((double)(nj+1));}

	for(int k=1;k<nk-1;k++){ acoef[k] = rhow[k]/(dz*dz); ccoef[k] = rhow[k+1]/(dz*dz);}

	acoef[1] = 0;
	ccoef[nk-2] = 0;

	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){

		B(i,j,1) = rhou[1]*2.0*(cos(kwave[i])+cos(lwave[j])-2.0)/(dx*dx) - rhow[2]/(dz*dz);

		B(i,j,nk-2) = rhou[nk-2]*2.0*(cos(kwave[i])+cos(lwave[j])-2.0)/(dx*dx) - (rhow[nk-2])/(dz*dz);

	for(int k=2;k<nk-2;k++){

		B(i,j,k) = rhou[k]*2.0*(cos(kwave[i])+cos(lwave[j])-2.0)/(dx*dx) - (rhow[k+1]+rhow[k])/(dz*dz);
	}}}

}

/**********************************************************************
* 
***********************************************************************/
void poisson_fft3d_real(int ni,int nj,int nk){

	int a = 88;//ni/2;
	int b = 30;//nj/2;
	int c = 1;//nk/2;

	/****************************************************
	* Set input to Fourier transform subroutine
	*****************************************************/
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<nk;k++){

		//INR3D(i,j,k) = divg_3d[i][j][k];

	}}}

	//printf("value = %e %e\n",INR3D(a,b,c),divg_3d[a][b][c]);

	/****************************************************
	* Execute Fourier transform
	*****************************************************/
	fftw_execute(p1); 

	double size2d = (double)(4*(ni+1)*(nj+1));

	for(int i=0;i<ni*nj*nk;i++){ outreal3d[i] /= size2d;}

	INR3D(a,b,c) = 0;

	//fftw_execute(p2);

	//printf("value = %e\n",INR3D(a,b,c));
	/****************************************************
	* Solve tridiagonal matrix
	*****************************************************/
	double z;


	double testarr[NZ];
	int d = 0;
	int e = 0;

	for(int k=0;k<NZ;k++){ testarr[k] = OUTR3D(d,e,k);}

	// lowest level
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){

		z = 1.0 / B(i,j,1);

		F(i,j,1) = ccoef[1]*z;

		OUTR3D(i,j,1) = OUTR3D(i,j,1)*z;
	}}

	// interior
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){

		for(int k=2;k<nk-2;k++){

			z = 1.0 / ( B(i,j,k) - acoef[k]*F(i,j,k-1));
		
			F(i,j,k) = ccoef[k]*z;

			OUTR3D(i,j,k) = (OUTR3D(i,j,k) - acoef[k]*OUTR3D(i,j,k-1))*z;
		}
	}}
	
	// highest level
	int kmax = NZ-2;

	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){

		z =  B(i,j,kmax) - acoef[kmax]*F(i,j,kmax-1);
		
		OUTR3D(i,j,kmax) = (OUTR3D(i,j,kmax)-acoef[kmax]*OUTR3D(i,j,kmax-1)) / z;
	}}



	/****************************************************
	* Back substitution
	*****************************************************/
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=kmax-1;k>=1;k--){

		OUTR3D(i,j,k) = OUTR3D(i,j,k) - F(i,j,k)*OUTR3D(i,j,k+1);

	}}}

	/****************************************************
	* Execute reverse Fourier transform
	*****************************************************/
	fftw_execute(p2);

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		INR3D(i,j,0) = INR3D(i,j,1);
		INR3D(i,j,NZ-1) = 0;
		INR3D(i,j,NZ-1) = INR3D(i,j,NZ-2);

	}}


	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<nk;k++){

		PI(i,j,k) = INR3D(i,j,k);
	}}}

}

/**********************************************************************
* Calculate the divergence of the new velocity field and write it
* directly to the input variable to the Fourier transform. Divide it
* by the time step to convert it to a divergence tendency.
*
* @param step - fractional time step in Runge-Kutta loop
***********************************************************************/
void divergence3d(double step){

	int ni = NX;
	int nj = NY;
	int nk = NZ;

	// divergence calculation
	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=1;k<NZ-1;k++){
		
		IN3D(i,j,k)[0] = rhou[k]*( UP(i+1,j,k)-UP(i,j,k) )*one_d_dx + 
						 rhou[k]*( VP(i,j+1,k)-VP(i,j,k) )*one_d_dy + 
		  			   ( rhow[k+1]*WP(i,j,k+1)-rhow[k]*WP(i,j,k) )*ONE_D_DZ(k);

		IN3D(i,j,k)[0] /= (step*dt);
	}}}

	// mirror lateral boundary condition
	for(int k=1;k<NZ-1;k++){

		for(int j=0;j<NY;j++){

			IN3D(0,j,k)[0] = IN3D(1,j,k)[0];
			IN3D(NX-1,j,k)[0] = IN3D(NX-2,j,k)[0];
		}
	
		for(int i=0;i<NX;i++){

			IN3D(i,0,k)[0] = IN3D(i,1,k)[0];
			IN3D(i,NY-1,k)[0] = IN3D(i,NY-2,k)[0];
		}
	}

	// upper and lower boundary condition
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		IN3D(i,j,0)[0] = IN3D(i,j,1)[0];
		IN3D(i,j,NZ-1)[0] = IN3D(i,j,NZ-2)[0];
	}}
}

/**********************************************************************
* Calculate the divergence of the new velocity field and write it
* directly to the input variable to the Fourier transform. Divide it
* by the time step to convert it to a divergence tendency.
*
* @param step - fractional time step in Runge-Kutta loop
***********************************************************************/
void divergence3d_periodic(int il,int ih,int jl,int jh,double step){

	int ni = NX-6;
	int nj = NY;
	int nk = NZ;

	// divergence calculation
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		
		IN3D(i-3,j,k)[0] = rhou[k]*( UP(i+1,j,k)-UP(i,j,k) )*one_d_dx + 
						   rhou[k]*( VP(i,j+1,k)-VP(i,j,k) )*one_d_dy + 
		  			   (   rhow[k+1]*WP(i,j,k+1)-rhow[k]*WP(i,j,k) )*ONE_D_DZ(k);

		IN3D(i-3,j,k)[0] /= (step*dt);
	}}}

	for(int k=1;k<NZ-1;k++){
	
		for(int i=0;i<ni;i++){

			IN3D(i,0,k)[0] = IN3D(i,1,k)[0];
			IN3D(i,NY-1,k)[0] = IN3D(i,NY-2,k)[0];
		}
	}

	// upper and lower boundary condition
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){

		IN3D(i,j,0)[0] = IN3D(i,j,1)[0];
		IN3D(i,j,NZ-1)[0] = IN3D(i,j,NZ-2)[0];
	}}
}
