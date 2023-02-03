#include "stdafx.h"
#include "fftw3.h"
#include "boundaries.h"

/***************************************************************************
* ------------------------------ MACROS ------------------------------------
****************************************************************************/

#define  IN2D(i,j)  in2d[(i)*NY+(j)]
#define OUT2D(i,j) out2d[(i)*NY+(j)]

#define  IN3D(i,j,k)  in3d[nj*((i)+ni*(k))+(j)]	// k,i,j array
#define OUT3D(i,j,k) out3d[nj*((i)+ni*(k))+(j)]

#define BCOEF(i,j) bcoef2d[(i)*NY+(j)]

// serial version
#define B2(i,j,k) bcoef3d[nk*((j)+nj*(i))+(k)]
#define F2(i,j,k)   rhs[nk*((j)+nj*(i))+(k)]

/***************************************************************************
* -------------------------- VARIABLES -------------------------------------
****************************************************************************/
//static fftw_complex *in2d = NULL, *out2d = NULL;
static double *in2d = NULL, *out2d = NULL;
static double *in3d = NULL, *out3d = NULL;
//static fftw_complex *in3d = NULL, *out3d = NULL;
static fftw_plan p1,p2,p3,p4;
static double *rhs;

static double *P,*Q;

double *kwave2d,*lwave2d;
double *acoef2d,*bcoef2d,*ccoef2d;

double *kwave3d,*lwave3d;
double *acoef3d,*bcoef3d,*ccoef3d;
bool solver_is_initialized = false;


struct laplacian_vars {

	fftw_plan forward,reverse;
	fftw_complex *in, *out;
	
	double *P,*Q,*rhs;
	double *kwave,*lwave;
	double *acoef,*bcoef,*ccoef;

};

struct laplacian_vars *vars;

/***************************************************************************
* ---------------------- FUNCTION PROTOTYPES--------------------------------
****************************************************************************/
void init_fft2d(int,int);
//void poisson_fft2d(double [NX][NY]);

/**********************************************************************
* Fast Fourier transform-based pressure solver for hydrostatic version
*
* ni,nj - grid dimensions for full domain array
***********************************************************************/
void init_laplacian_solver_periodic_EW_zerograd_NS(int ni,int nj){

	if(!solver_is_initialized){

		solver_is_initialized = true;

		vars = (laplacian_vars*) malloc(sizeof(laplacian_vars));

		vars->in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NX*nj);
		vars->out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NX*nj);

		int n[] = {NX};

		vars->forward = fftw_plan_many_dft(1, n, nj,vars->in,n,nj,1,vars->out,n,nj,1,FFTW_FORWARD, FFTW_MEASURE);
		vars->reverse = fftw_plan_many_dft(1, n, nj,vars->out,n,nj,1,vars->in,n,nj,1,FFTW_BACKWARD,FFTW_MEASURE);

		vars->kwave = (double*) calloc(sizeof(double),NX);
		vars->lwave = (double*) calloc(sizeof(double),NY);

		for(int i=0;i<NX/2;i++){	vars->kwave[i] = 2*trigpi * (double)i / NX;	}
		for(int i=NX/2;i<NX;i++){	vars->kwave[i] = 2*trigpi * (double)(NX-i) / NX;}

		for(int i=0;i<NY/2;i++){	vars->lwave[i] = 2*trigpi * (double)i / NY;	}
		for(int i=NY/2;i<NY;i++){	vars->lwave[i] = 2*trigpi * (double)(NY-i) / NY;}

		vars->acoef = (double*) malloc(sizeof(double)*NY);
		vars->bcoef = (double*) malloc(sizeof(double)*NX*NY);
		vars->ccoef = (double*) malloc(sizeof(double)*NY);

		vars->P = (double*) calloc(sizeof(double),(NY+1));
		vars->Q = (double*) calloc(sizeof(double),(NY+1));

		for(int j=0;j<NY;j++){

			vars->acoef[j] = 1./(dy*dy);
			vars->ccoef[j] = 1./(dy*dy);

		}

		//vars->acoef[0] = 0;	
		//vars->ccoef[NY-1] = 0;

		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){

			vars->bcoef[i*NY+j] = 2.*(cos(vars->kwave[i])-1.)/(dx*dx) - vars->ccoef[j] - vars->acoef[j];
		}}
	}

}

/**********************************************************************
* Fast Fourier transform-based pressure solver
***********************************************************************/
void init_laplacian_solver_real(int ni,int nj){

	kwave2d = (double*) malloc(sizeof(double)*ni);
	lwave2d = (double*) malloc(sizeof(double)*nj);

	acoef2d = (double*) malloc(sizeof(double)*NY);
	bcoef2d = (double*) malloc(sizeof(double)*NX*NY);
	ccoef2d = (double*) malloc(sizeof(double)*NY);

	P = (double*) calloc(sizeof(double),(nj+1));
	Q = (double*) calloc(sizeof(double),(nj+1));

	in2d = (double*) calloc(sizeof(double),ni*nj);
	out2d = (double*) calloc(sizeof(double),ni*nj);

	int n[] = {NX};
	//int *inembed = n;
	//int *onembed = n;
	
	fftw_r2r_kind kind1[] = {FFTW_RODFT00};
	fftw_r2r_kind kind2[] = {FFTW_RODFT00};

	p3 = fftw_plan_many_r2r(1, n, nj, in2d, n, nj,1,out2d, n, nj,1,kind1, FFTW_MEASURE);
	p4 = fftw_plan_many_r2r(1, n, nj,out2d, n, nj,1, in2d, n, nj,1,kind2, FFTW_MEASURE);
	//------------
	//in2d = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NX*nj);
	//out2d = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NX*nj);

	//int n[] = {NX};

	//p3 = fftw_plan_many_dft(1, n, nj,in2d,  n,nj,1,out2d,n,nj,1,FFTW_FORWARD, FFTW_MEASURE);
	//p4 = fftw_plan_many_dft(1, n, nj,out2d, n,nj,1,in2d,n,nj,1,FFTW_BACKWARD,FFTW_MEASURE);
	//--------------
	double shift = 1;

	for(int i=0;i<ni;i++){  kwave2d[i] = ((double)(i)+shift) * trigpi/((double)(ni+1));}

	for(int i=0;i<nj;i++){  lwave2d[i] = ((double)(i)+shift) * trigpi/((double)(nj+1));}

	for(int j=0;j<NY;j++){

		acoef2d[j] = 1./(dy*dy);
		ccoef2d[j] = 1./(dy*dy);

	}

	//acoef2d[0] = 0;	
	//ccoef2d[NY-1] = 0;
	//ccoef2d[NY-2] = 0;

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		BCOEF(i,j) = 2.*(cos(kwave2d[i])-1.)/(dx*dx) - ccoef2d[j] - acoef2d[j];
	}}

}

/**********************************************************************
* Fast Fourier transform-based pressure solver
***********************************************************************/
void run_laplacian_solver(double *Fxy){

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		vars->in[i*NY+j][0] = Fxy[i*NY+j];
		vars->in[i*NY+j][1] = 0;

	}}

	fftw_execute(vars->forward);

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		vars->out[i*NY+j][0] = vars->out[i*NY+j][0] / NX;
		vars->out[i*NY+j][1] = vars->out[i*NY+j][1] / NX;

	}}

	/*****************************************************************
	* Solve tridiagonal matrix.
	******************************************************************/
	double denom;

	for(int i=0;i<NX;i++){

		// forward elimination
		for(int j=0;j<NY;j++){

			denom = vars->bcoef[i*NY+j] + vars->acoef[j]*vars->P[j];

			vars->P[j+1] = -vars->ccoef[j] / denom;
			vars->Q[j+1] = (vars->out[i*NY+j][0]-vars->acoef[j]*vars->Q[j]) / denom;

		}

		// backward substitution
		for(int j=NY-2;j>=0;j--){ vars->out[i*NY+j][0] = vars->P[j+1]*vars->out[i*NY+j+1][0] + vars->Q[j+1];}

		// forward elimination
		for(int j=0;j<NY;j++){

			denom = vars->bcoef[i*NY+j] + vars->acoef[j]*vars->P[j];

			vars->P[j+1] = -vars->ccoef[j] / denom;
			vars->Q[j+1] = (vars->out[i*NY+j][1] - vars->acoef[j]*vars->Q[j]) / denom;

		}

		// backward substitution
		for(int j=NY-2;j>=0;j--){ vars->out[i*NY+j][1] = vars->P[j+1]*vars->out[i*NY+j+1][1] + vars->Q[j+1];}

	}

	fftw_execute(vars->reverse);
	

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		Fxy[i*NY+j] = vars->in[i*NY+j][0];

	}}

}


/**********************************************************************
* Fast Fourier transform-based pressure solver
***********************************************************************/
void run_laplacian_solver_real(double *Fxy){

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		IN2D(i,j) = Fxy[i*NY+j];

	}}

	fftw_execute(p3);

	double size2d = (double)(2*(NX+1));

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		OUT2D(i,j) /= size2d;

	}}

	/*****************************************************************
	* Solve tridiagonal matrix.
	******************************************************************/
	double denom;

	for(int i=0;i<NX;i++){

		// forward elimination
		for(int j=0;j<NY;j++){

			denom = BCOEF(i,j) + acoef2d[j]*P[j];

			P[j+1] = -ccoef2d[j] / denom;
			Q[j+1] = (OUT2D(i,j)-acoef2d[j]*Q[j]) / denom;

		}

		// backward substitution
		for(int j=NY-2;j>=0;j--){ OUT2D(i,j) = P[j+1]*OUT2D(i,j+1) + Q[j+1];}

	}

	fftw_execute(p4);
	

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		Fxy[i*NY+j] = IN2D(i,j);

	}}

}

/**********************************************************************
* Fast Fourier transform-based pressure solver for serial version
*
* ni,nj,nk - grid dimensions for full domain array
***********************************************************************/
void init_laplacian_solver3d(int ni,int nj,int nk){

	double f0 = 2*7.292e-5*sin(20*3.1416/180.0);

	double N2[NZ];

	kwave3d = (double*) malloc(sizeof(double)*ni);
	lwave3d = (double*) malloc(sizeof(double)*nj);

	acoef3d = (double*) malloc(sizeof(double)*nk);
	bcoef3d = (double*) malloc(sizeof(double)*ni*nj*nk);
	ccoef3d = (double*) malloc(sizeof(double)*nk);

	in3d = (double*) fftw_malloc(sizeof(double)*ni*nj*nk);
	out3d = (double*) fftw_malloc(sizeof(double)*ni*nj*nk);
	rhs = (double*) malloc(sizeof(double)*ni*nj*nk);

	for(int i=0;i<ni*nj*nk;i++){ 
		
		in3d[i] = 0;
	}

	for(int k=1;k<NZ-1;k++){ N2[k] = (grav / tb[k]) * 0.5*(tb[k+1]-tb[k-1]) * ONE_D_DZ(k);}

	N2[0] = N2[1];
	N2[NZ-1] = N2[NZ-2];

	int n[] = {ni,nj};
	int *inembed = n;
	int *onembed = n;

	fftw_r2r_kind kind1[] = {FFTW_RODFT00,FFTW_RODFT00};
	fftw_r2r_kind kind2[] = {FFTW_RODFT00,FFTW_RODFT00};

	p1 = fftw_plan_many_r2r(2, n, nk, in3d,inembed, 1,ni*nj,out3d,onembed, 1,ni*nj,kind1, FFTW_MEASURE);
	p2 = fftw_plan_many_r2r(2, n, nk,out3d,inembed, 1,ni*nj,in3d, onembed, 1,ni*nj,kind2, FFTW_MEASURE);



	//for(int i=0;i<ni;i++){  kwave3d[i] = (double)i * 2*trigpi/((double)ni);}
	//for(int i=0;i<nj;i++){  lwave3d[i] = (double)i * 2*trigpi/((double)nj);}

	double shift = 1;

	for(int i=0;i<ni;i++){  kwave3d[i] = ((double)(i)+shift) * trigpi/((double)(ni+1));}

	for(int i=0;i<nj;i++){  lwave3d[i] = ((double)(i)+shift) * trigpi/((double)(nj+1));}

	if(!STRETCHED_GRID){

		for(int k=1;k<nk-1;k++){ acoef3d[k] = 1.0/(dz*dz); ccoef3d[k] = 1.0/(dz*dz);}

		for(int i=0;i<ni;i++){
		for(int j=0;j<nj;j++){

			if(i==0 && j==0){ }
			else {
				B2(i,j,1) = 1.0*2.0*(cos(kwave3d[i])+cos(lwave3d[j])-2.0)/(dx*dx) - 1.0/(dz*dz);
			}

			B2(i,j,nk-2) = 1.0*2.0*(cos(kwave3d[i])+cos(lwave3d[j])-2.0)/(dx*dx) - 1.0/(dz*dz);
		
		for(int k=2;k<nk-2;k++){

			B2(i,j,k) = 1.0*2.0*(cos(kwave3d[i])+cos(lwave3d[j])-2.0)/(dx*dx) - (1.0+1.0)/(dz*dz);
		}}}
		
	} else {

		for(int k=0;k<nk;k++){ acoef3d[k] = (f0*f0)*mu[k]*mw[k]*1.0/(dz*dz); ccoef3d[k] = (f0*f0)*mu[k]*mw[k+1]*1.0/(dz*dz);}

		for(int i=0;i<ni;i++){
		for(int j=0;j<nj;j++){

			//B2(i,j,1) = N2[1]*2.0*(cos(kwave3d[i])+cos(lwave3d[j])-2.0)/(dx*dx); //- (f0*f0)*mu[1]*mw[2]*1.0/(dz*dz);

			for(int k=0;k<nk;k++){
				B2(i,j,k) = N2[k]*2.0*(cos(kwave3d[i])+cos(lwave3d[j])-2.0)/(dx*dx) - ccoef3d[k] - acoef3d[k];
			
			}
		
			//B2(i,j,nk-2) = N2[nk-2]*2.0*(cos(kwave3d[i])+cos(lwave3d[j])-2.0)/(dx*dx); //- (f0*f0)*(mu[nk-2]*mw[nk-2]*1.0)/(dz*dz);
		}}
		/*
		B2(i,j,1) = N2[1]*2.0*(cos(kwave3d[i])+cos(lwave3d[j])-2.0)/(dx*dx) - (f0*f0)*mu[1]*mw[2]*1.0/(dz*dz);

		for(int k=2;k<nk-2;k++){ B2(i,j,k) = N2[k]*2.0*(cos(kwave3d[i])+cos(lwave3d[j])-2.0)/(dx*dx) - (f0*f0)*mu[k]*(mw[k+1]*1.0+mw[k]*1.0)/(dz*dz);}
	
		B2(i,j,nk-2) = N2[nk-2]*2.0*(cos(kwave3d[i])+cos(lwave3d[j])-2.0)/(dx*dx) - (f0*f0)*(mu[nk-2]*mw[nk-2]*1.0)/(dz*dz);
		*/
	}
//2.*(cos(kwave2d[i])-1.)/(dx*dx) - ccoef2d[j] - acoef2d[j];

}

/**********************************************************************
* Solve the anelastic pressure equation using Fast Fourier transforms
* in the horizontal and Gaussian elimination in the vertical
*
* @param ni,nj,nk - grid dimensions
***********************************************************************/
void run_laplacian_solver3d(int ni,int nj,int nk,double *input){

	/****************************************************
	* Set input to Fourier transform subroutine
	*****************************************************/
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<nk;k++){
	
		IN3D(i,j,k) = input[INDEX(i,j,k)];
	}}}

	/****************************************************
	* Execute Fourier transform
	*****************************************************/
	fftw_execute(p1); 

	double size2d = (double)(4*(ni+1)*(nj+1));

	// normalize
	for(int i=0;i<ni*nj*nk;i++){ out3d[i] /= size2d;}

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

		F2(i,j,1) = ccoef3d[1]*z;

		OUT3D(i,j,1) = OUT3D(i,j,1)*z;

		for(int k=2;k<nk-1;k++){

			z = 1.0 / ( B2(i,j,k) - acoef3d[k]*F2(i,j,k-1));
		
			F2(i,j,k) = ccoef3d[k]*z;

			OUT3D(i,j,k) = (OUT3D(i,j,k) - acoef3d[k]*OUT3D(i,j,k-1))*z;
		}
	}}
	
	OUT3D(0,0,NZ-2) = 0;

	/***********************************************
	* Back substitution
	************************************************/
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=NZ-3;k>=1;k--){

		OUT3D(i,j,k) -= F2(i,j,k)*OUT3D(i,j,k+1);

	}}}

	/****************************************************
	* Execute reverse Fourier transform
	*****************************************************/
	fftw_execute(p2);

	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<nk;k++){
	
		input[INDEX(i,j,k)] = IN3D(i,j,k);
		
	}}}
}
