#include "stdafx.h"
#include "fftw3.h"
#include "damping.h"

#define FINDEX(i,j,k) nj*((i)+ni*(k))+(j)

#if !STRETCHED_GRID

#define DIFFUSE_VARIABLE(vp,vm)	vp(i,j,k) = vp(i,j,k) + step*dt* ( 	 \
			+ kmixhx * (vm(i+1,j,k) - 2.0*vm(i,j,k)+vm(i-1,j,k))*one_d_dx2 \
			+ kmixhy * (vm(i,j+1,k) - 2.0*vm(i,j,k)+vm(i,j-1,k))*one_d_dy2 \
			+ kmixv  * (vm(i,j,k+1) - 2.0*vm(i,j,k)+vm(i,j,k-1))*one_d_dz2 \
			)
#else
				
#define DIFFUSE_VARIABLE(vp,vm)	vp(i,j,k) = vp(i,j,k) + step*dt* ( 	 \
			+ kmixhx * (vm(i+1,j,k) - 2.0*vm(i,j,k)+vm(i-1,j,k))*one_d_dx2 \
			+ kmixhy * (vm(i,j+1,k) - 2.0*vm(i,j,k)+vm(i,j-1,k))*one_d_dy2 \
			+ kmixv  * mu[k]*(mw[k+1]*vm(i,j,k+1) - mw[k+1]*vm(i,j,k) - mw[k]*vm(i,j,k) + mw[k]*vm(i,j,k-1))*one_d_dz2 \
			)		
#endif

double *weights;
double *temp_var;

const double cmixhx = 0.05;
const double cmixhy = 0.05;
const double cmixv = 0.1;

double kmixhx;
double kmixhy;
double kmixv;

const double diffh_coef_2nd = 0.005;//0.005;
const double diffv_coef_2nd = 0.0005;//0.0005;
const double diffh_coef_4th = 0.01;
const double diffv_coef_4th = 0.001;
const double diffh_coef_6th = 0.012*0.5; //0.012;
const double diffv_coef_6th = 0.0012*0.5; //0.0012;
const double diffh_coef_8th = 0.012;
const double diffv_coef_8th = 0.0012;

double kmixh2nd;
double kmixv2nd;
double kmixh4th;
double kmixv4th;
double kmixh6th;
double kmixv6th;
double kmixh8th;
double kmixv8th;

double asselincoef = 0.1;

double *kwave2,*lwave2;
double *u_diff_tend,*v_diff_tend;

double dx2; double dy2; double dz2;
double dx4; double dy4; double dz4;
double dx6; double dy6; double dz6;
double dx8; double dy8; double dz8;

double one_d_dx2; double one_d_dy2; double one_d_dz2;
double one_d_dx4; double one_d_dy4; double one_d_dz4;
double one_d_dx6; double one_d_dy6; double one_d_dz6;
double one_d_dx8; double one_d_dy8; double one_d_dz8;

const double onetwelfth = 1.0/12.0;
const double fourthirds = 4.0/3.0;
const double fivehalves = 5.0/2.0;

fftw_complex *u_fft_in,*u_fft_out;
fftw_plan u_fft_plan_forward,u_fft_plan_reverse;

/*********************************************************************
* FUNCTION PROTYPES
**********************************************************************/
double diffuse_i_2nd(double *var,int i,int j,int k);
double diffuse_j_2nd(double *var,int i,int j,int k);
double diffuse_k_2nd(double *var,int i,int j,int k);

double diffuse_i_4th(double *var,int i,int j,int k);
double diffuse_j_4th(double *var,int i,int j,int k);
double diffuse_k_4th(double *var,int i,int j,int k);

double diffuse_i_8th(double *var,int i,int j,int k);
double diffuse_j_8th(double *var,int i,int j,int k);
double diffuse_k_8th(double *var,int i,int j,int k);

void diffusion_8th_var(double *pvar,double *mvar,double step,int il,int ih,int jl,int jh,double *output=NULL);
void diffusion_6th_var(double *pvar,double *mvar,double step,int il,int ih,int jl,int jh,double *output=NULL);
void diffusion_4th_var(double *pvar,double *mvar,double step,int il,int ih,int jl,int jh,double *output=NULL);
void diffusion_2nd_var(double *pvar,double *mvar,double step,int il,int ih,int jl,int jh);

void diffusion_8th_all(double step,int il,int ih,int jl,int jh);
void diffusion_6th_all(double step,int il,int ih,int jl,int jh);
void diffusion_4th_all(double step,int il,int ih,int jl,int jh);
void diffusion_2nd_all(double step,int il,int ih,int jl,int jh);

/*********************************************************************
* INITIALIZER
**********************************************************************/
void init_damping(int ni,int nj,int nk){
	
	kmixhx = cmixhx*dx*dx/dt;
	kmixhy = cmixhy*dy*dy/dt;
	kmixv = cmixv*dz*dz/dt;

	kmixh2nd = diffh_coef_2nd*dx*dx/dt;
	kmixv2nd = diffv_coef_2nd*dz*dz/dt;
	kmixh4th = diffh_coef_4th*dx*dx*dx*dx/dt;
	kmixv4th = diffv_coef_4th*dz*dz*dz*dz/dt;
	kmixh6th = diffh_coef_6th*dx*dx*dx*dx*dx*dx/dt;
	kmixv6th = diffv_coef_6th*dz*dz*dz*dz*dz*dz/dt;
	kmixh8th = diffh_coef_8th*dx*dx*dx*dx*dx*dx*dx*dx/dt;
	kmixv8th = diffv_coef_8th*dz*dz*dz*dz*dz*dz*dz*dz/dt;

	dx2 = dx*dx;   dy2 = dy*dy;   dz2 = dz*dz;
	dx4 = dx2*dx2; dy4 = dy2*dy2; dz4 = dz2*dz2;
	dx6 = dx2*dx4; dy6 = dy2*dy4; dz6 = dz2*dz4;
	dx8 = dx4*dx4; dy8 = dy4*dy4; dz8 = dz4*dz4;

	one_d_dx2 = 1.0/dx2; one_d_dy2 = 1.0/dy2; one_d_dz2 = 1.0/dz2;
	one_d_dx4 = 1.0/dx4; one_d_dy4 = 1.0/dy4; one_d_dz4 = 1.0/dz4;
	one_d_dx6 = 1.0/dx6; one_d_dy6 = 1.0/dy6; one_d_dz6 = 1.0/dz6;
	one_d_dx8 = 1.0/dx8; one_d_dy8 = 1.0/dy8; one_d_dz8 = 1.0/dz8;
	
	if(OUTPUT_DIFFUSION_TEND){
		
		u_diff_tend = (double*)calloc(ni*nj*nk,sizeof(double));
		v_diff_tend = (double*)calloc(ni*nj*nk,sizeof(double));
	}
	
	if(EXTRA_DIFFUSION){
		//printf("%d %d %d\n",ni,nj,nk);
		 temp_var = (double*)calloc(ni*nj*nk,sizeof(double));
	}
	
}

/*********************************************************************
* INITIALIZER
**********************************************************************/
void init_diffusion_weights(int deriv_order,double *heights){
	
	int lower_order,k;
	
	double *weights = (double*)calloc((deriv_order+1)*NZ,sizeof(double));
	
	double *weights_point = (double*)calloc((deriv_order+1)*(deriv_order+1),sizeof(double));
	
	//derivative_weights(zsu[9+3],zsu+9, deriv_order, deriv_order, deriv_order,weights_point);

	for(k=deriv_order/2;k<NZ-deriv_order/2;k++){
		//printf("%d %f\n",k,zsu[k]);
		memset(weights_point,0,(deriv_order+1)*(deriv_order+1)*sizeof(double));
		
		derivative_weights(heights[k], heights-deriv_order/2+k, deriv_order, deriv_order, deriv_order, weights_point);
		
		for(int j=0;j<=deriv_order;j++){
		
			weights[index2d(deriv_order+1,k,j)] = weights_point[index2d(deriv_order+1,j,deriv_order)];
		}
	}
	
	if(deriv_order > 4){
		
		lower_order = 4;
		
		k = 2;
		
		memset(weights_point,0,(deriv_order+1)*(deriv_order+1)*sizeof(double));
		
		derivative_weights(heights[k], heights+2, lower_order, lower_order, lower_order, weights_point);
		
		for(int j=0;j<=lower_order;j++){
		
			weights[index2d(deriv_order+1,k,j)] = weights_point[index2d(lower_order+1,j,lower_order)];
			
		}
		
		k = NZ-3;
		
		memset(weights_point,0,(deriv_order+1)*(deriv_order+1)*sizeof(double));
		
		derivative_weights(heights[k], heights+2, lower_order, lower_order, lower_order, weights_point);
		
		for(int j=0;j<=lower_order;j++){
		
			weights[index2d(deriv_order+1,k,j)] = weights_point[index2d(lower_order+1,j,lower_order)];
			
		}
		
	}
	
	if(deriv_order > 2){
		
		lower_order = 2;
		
		k = 1;
		
		memset(weights_point,0,(deriv_order+1)*(deriv_order+1)*sizeof(double));
		
		derivative_weights(heights[k], heights+2, lower_order, lower_order, lower_order, weights_point);
		
		for(int j=0;j<=lower_order;j++){
		
			weights[index2d(deriv_order+1,k,j)] = weights_point[index2d(lower_order+1,j,lower_order)];
			
		}
		
		k = NZ-2;
		
		memset(weights_point,0,(deriv_order+1)*(deriv_order+1)*sizeof(double));
		
		derivative_weights(heights[k], heights+2, lower_order, lower_order, lower_order, weights_point);
		
		for(int j=0;j<=lower_order;j++){
		
			weights[index2d(deriv_order+1,k,j)] = weights_point[index2d(lower_order+1,j,lower_order)];
			
		}
		
	}
	
#if 0
	for(int j=0;j<=deriv_order;j++){
		//printf("%d\n",j);
	
		for(int i=0;i<=deriv_order;i++){
		//printf("%f  \n",weights_point[index2d(deriv_order+1,i,j)] );
	}}
	
	for(int k=1;k<NZ-1;k++){
		printf("%d\n",k);
	
		for(int i=0;i<=deriv_order;i++){
		printf("%e  \n",kmixv6th*weights[index2d(deriv_order+1,k,i)] );
	}}
#endif	
	free(weights_point);
	
}

/*********************************************************************
*
*
**********************************************************************/
void init_fftw(int ni,int nj,int nk){

	int n[] = {ni,nj};
	int *inembed = n;
	int *onembed = n;

	u_fft_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ni*nj*nk);
	u_fft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ni*nj*nk);
	
	u_fft_plan_forward = fftw_plan_many_dft(2, n, nk, u_fft_in,inembed, 1,ni*nj,u_fft_out,onembed, 1,ni*nj,FFTW_FORWARD, FFTW_MEASURE);
	u_fft_plan_reverse = fftw_plan_many_dft(2, n, nk, u_fft_out,inembed, 1,ni*nj,u_fft_in,onembed, 1,ni*nj,FFTW_BACKWARD, FFTW_MEASURE);	

	kwave2 = (double*) malloc(sizeof(double)*ni);
	lwave2 = (double*) malloc(sizeof(double)*nj);
	
	for(int i=0;i<=ni/2;i++){  kwave2[i] = (double)i * 2*trigpi/((double)ni);}
	for(int i=0;i<=nj/2;i++){  lwave2[i] = (double)i * 2*trigpi/((double)nj);}
	for(int i=ni/2+1;i<ni;i++){  kwave2[i] = ((double)(ni-i)) * 2*trigpi/((double)ni);}
	for(int i=nj/2+1;i<nj;i++){  lwave2[i] = ((double)(nj-i)) * 2*trigpi/((double)nj);}
	/*
	for(int i=0;i<=ni/2;i++){  kwave2[i] = 0.001*(ni*dx-dx*ni*i/((double)ni)); printf("%d %f\n",i,kwave2[i]);}
	for(int i=0;i<=nj/2;i++){  lwave2[i] = (double)i * 2*trigpi/((double)nj);}
	for(int i=ni/2+1;i<ni;i++){  kwave2[i] = ((double)(ni-i)) * 2*trigpi/((double)ni);}
	for(int i=nj/2+1;i<nj;i++){  lwave2[i] = ((double)(nj-i)) * 2*trigpi/((double)nj);}
	*/
	//exit(0);
}

/*********************************************************************
* Wavenumber damping using Fourier transform. This is the version
* for periodic boundaries
**********************************************************************/
void fft_damp(int ni,int nj,int nk,int b,double *var){
	
	double kh = sqrt(kwave2[WAVE_NUMBER-1]*kwave2[WAVE_NUMBER-1])+0.01;
	double kl = sqrt(kwave2[WAVE_NUMBER+1]*kwave2[WAVE_NUMBER+1])-0.01;
	
	//printf("%f %f\n",kh,kl);
	
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<nk;k++){
	
		u_fft_in[FINDEX(i,j,k)][0] = var[INDEX(i+b,j,k)];
		u_fft_in[FINDEX(i,j,k)][1] = 0;
		
	}}}
	
	fftw_execute(u_fft_plan_forward);

	// normalize	
	double size2d = (double)(ni*nj);
	
	for(int i=0;i<ni*nj*nk;i++){ u_fft_out[i][0] /= size2d; u_fft_out[i][1] /= size2d;}
	
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<nk;k++){
	
		//if(sqrt(kwave2[i]*kwave2[i]+lwave2[j]*lwave2[j]) > trigpi*0.9){
		//if(sqrt(lwave2[j]*lwave2[j]) > 1.0 && sqrt(kwave2[i]*kwave2[i]) > 1.0){	// exp8
		//if(sqrt(lwave2[j]*lwave2[j]) > 0.75 && sqrt(kwave2[i]*kwave2[i]) > 0.25){
		//if(sqrt(kwave2[i]*kwave2[i]) < 0.07){
		//if(sqrt(kwave2[i]*kwave2[i]) < 0.14 || sqrt(kwave2[i]*kwave2[i]) > 0.20){// || sqrt(lwave2[j]*lwave2[j]) > 2.5){
		if(sqrt(kwave2[i]*kwave2[i]) < kh || sqrt(kwave2[i]*kwave2[i]) > kl){
		//if(i!=2){
			u_fft_out[FINDEX(i,j,k)][0] = 0;//0.95;
			u_fft_out[FINDEX(i,j,k)][1] = 0;//0.95;
		}
	}}}
	
	fftw_execute(u_fft_plan_reverse);
	
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<nk;k++){
	
		var[INDEX(i+b,j,k)] = u_fft_in[FINDEX(i,j,k)][0];

	}}}
}

/*********************************************************************
* Wavenumber damping using Fourier transform. This is the version
* for non-periodic boundaries
**********************************************************************/
void fft_damp2(int ni,int nj,int nk,int b,double *var){

	const int kc = 30;
	const int lc = 10; 
	double lh = sqrt(lwave2[1]*lwave2[1])+0.01;
	double ll = sqrt(lwave2[NY/2-lc]*lwave2[NY/2-lc])-0.01;	
	double kh = sqrt(kwave2[1]*kwave2[1])+0.01;
	double kl = sqrt(kwave2[NX/2-kc]*kwave2[NX/2-kc])-0.01;
	
	//printf("%f %f\n",kl,ll);
	
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<nk;k++){
	
		u_fft_in[FINDEX(i,j,k)][0] = var[INDEX(i+b,j,k)];
		u_fft_in[FINDEX(i,j,k)][1] = 0;
		
	}}}
	
	fftw_execute(u_fft_plan_forward);

	// normalize	
	double size2d = (double)(ni*nj);
	
	for(int i=0;i<ni*nj*nk;i++){ u_fft_out[i][0] /= size2d; u_fft_out[i][1] /= size2d;}
	
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<nk;k++){
	
		if(sqrt(kwave2[i]*kwave2[i]) > kl || sqrt(lwave2[j]*lwave2[j]) > ll){

			u_fft_out[FINDEX(i,j,k)][0] = 0;
			u_fft_out[FINDEX(i,j,k)][1] = 0;
		}
	}}}
	
	fftw_execute(u_fft_plan_reverse);
	
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<nk;k++){
	
		var[INDEX(i+b,j,k)] = u_fft_in[FINDEX(i,j,k)][0];

	}}}
}

/*********************************************************************
*
*
**********************************************************************/
void run_fftw(int ni,int nj,int nk,int b){

	//printf("u = %f\n",UP(NX/2,NY/2,5));

	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<nk;k++){
	
		u_fft_in[FINDEX(i,j,k)][0] = PI(i+b,j,k);
		u_fft_in[FINDEX(i,j,k)][1] = 0;
		
	}}}
	
	fftw_execute(u_fft_plan_forward);

	//printf("ufft = %f\n",u_fft_out[INDEX(NX/2,NY/2,5)][0]);

	// normalize	
	double size2d = (double)(ni*nj);

	for(int i=0;i<ni*nj*nk;i++){ u_fft_out[i][0] /= size2d; u_fft_out[i][1] /= size2d;}

	//for(int i=0;i<ni*nj*nk;i++){ u_fft_in[i][0] = 0; u_fft_in[i][1] = 0;}
#if 0
	for(int i=0;i<ni;i++){
	
		//printf("k = %f %f\n",kwave2[i],u_fft_out[INDEX(i,NY/2,5)][0]);
		//printf("k = %d %f\n",i,kwave2[i]);
	}
	for(int j=0;j<nj;j++){

		//printf("l = %d %f\n",j,lwave2[j]);
	}
	
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<nk;k++){
	
		//if(sqrt(kwave2[i]*kwave2[i]+lwave2[j]*lwave2[j]) > trigpi*0.25){
	
			//u_fft_out[FINDEX(i,j,k)][0] = 0;//0.75;
			//u_fft_out[FINDEX(i,j,k)][1] = 0;//0.75;
		//}
	}}}
	
	//fftw_execute(u_fft_plan_reverse);
	
	for(int i=3;i<ni-3;i++){
	for(int j=3;j<nj-3;j++){
	for(int k=1;k<nk-1;k++){
	
		//UP(i,j,k) = u_fft_in[FINDEX(i,j,k)][0];

	}}}
	
#endif
	
	for(int i=0;i<ni;i++){
	for(int j=0;j<nj;j++){
	for(int k=0;k<nk;k++){
	
		UP(i+b,j,k) = u_fft_out[FINDEX(i,j,k)][0];
		
	}}}
	
	//printf("ufft = %f\n",u_fft_in[INDEX(NX/2,NY/2,5)][0]);
	
}

/*********************************************************************
*
*
**********************************************************************/
void asselin(){

	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=1;k<NZ-1;k++){

		U(i,j,k) += asselincoef*(UP(i,j,k)-2.0*U(i,j,k)+UM(i,j,k));
		V(i,j,k) += asselincoef*(VP(i,j,k)-2.0*V(i,j,k)+VM(i,j,k));
		W(i,j,k) += asselincoef*(WP(i,j,k)-2.0*W(i,j,k)+WM(i,j,k));
		TH(i,j,k) += asselincoef*(THP(i,j,k)-2.0*TH(i,j,k)+THM(i,j,k));
	}}}
}

/*********************************************************************
*
*
**********************************************************************/
double diffuse_i_2nd(double *var,int i,int j,int k){
	
	return var[INDEX(i+1,j,k)] - 2.0*var[INDEX(i,j,k)] + var[INDEX(i-1,j,k)];
}

/*********************************************************************
*
*
**********************************************************************/
double diffuse_j_2nd(double *var,int i,int j,int k){
	
	return var[INDEX(i,j+1,k)] - 2.0*var[INDEX(i,j,k)] + var[INDEX(i,j-1,k)];
}

/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
double diffuse_k_2nd(double *var,int i,int j,int k){
	
	return var[INDEX(i,j,k+1)] - 2.0*var[INDEX(i,j,k)] + var[INDEX(i,j,k-1)];
}

/*********************************************************************
*
*
**********************************************************************/
double diffuse_k_2nd(double *var,double *weights,int i,int j,int k){
	
	return	weights[0] * var[INDEX(i,j,k-1)] -
		 	weights[1] * var[INDEX(i,j,k  )] +
			weights[2] * var[INDEX(i,j,k+1)];
}

/*********************************************************************
*
*
**********************************************************************/
double diffuse_i_2nd_fourthorder(double *var,int i,int j,int k){
	
	return -onetwelfth*var[INDEX(i-2,j,k)] +
			fourthirds*var[INDEX(i-1,j,k)] -
			fivehalves*var[INDEX(i  ,j,k)] +
			fourthirds*var[INDEX(i+1,j,k)] -
			onetwelfth*var[INDEX(i+2,j,k)];
}

/*********************************************************************
*
*
**********************************************************************/
double diffuse_j_2nd_fourthorder(double *var,int i,int j,int k){
	
	return -onetwelfth*var[INDEX(i,j-2,k)] +
			fourthirds*var[INDEX(i,j-1,k)] -
			fivehalves*var[INDEX(i,j  ,k)] +
			fourthirds*var[INDEX(i,j+1,k)] -
			onetwelfth*var[INDEX(i,j+2,k)];
}

/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
double diffuse_k_2nd_fourthorder(double *var,int i,int j,int k){
	
	return -onetwelfth*var[INDEX(i,j,k-2)] +
			fourthirds*var[INDEX(i,j,k-1)] -
			fivehalves*var[INDEX(i,j,k  )] +
			fourthirds*var[INDEX(i,j,k+1)] -
			onetwelfth*var[INDEX(i,j,k+2)];
}

/*********************************************************************
*
*
**********************************************************************/
double diffuse_i_4th(double *var,int i,int j,int k){
	
	return 		var[INDEX(i-2,j,k)] -
			4.0*var[INDEX(i-1,j,k)] +
			6.0*var[INDEX(i  ,j,k)] -
			4.0*var[INDEX(i+1,j,k)] +
				var[INDEX(i+2,j,k)];
}

/*********************************************************************
*
*
**********************************************************************/
double diffuse_j_4th(double *var,int i,int j,int k){
	
	return 		var[INDEX(i,j-2,k)] -
			4.0*var[INDEX(i,j-1,k)] +
			6.0*var[INDEX(i,j  ,k)] -
			4.0*var[INDEX(i,j+1,k)] +
				var[INDEX(i,j+2,k)];
}

/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
double diffuse_k_4th(double *var,int i,int j,int k){
	
	return 		var[INDEX(i,j,k-2)] -
			4.0*var[INDEX(i,j,k-1)] +
			6.0*var[INDEX(i,j,k  )] -
			4.0*var[INDEX(i,j,k+1)] +
				var[INDEX(i,j,k+2)];
}

/*********************************************************************
*
*
**********************************************************************/
double diffuse_k_4th(double *var,double *weights,int i,int j,int k){
	
	return	weights[0] * var[INDEX(i,j,k-2)] -
		 	weights[1] * var[INDEX(i,j,k-1)] +
			weights[2] * var[INDEX(i,j,k  )] -
			weights[3] * var[INDEX(i,j,k+1)] +
			weights[4] * var[INDEX(i,j,k+2)];
}

/*********************************************************************
*
*
**********************************************************************/
double diffuse_i_6th(double *var,int i,int j,int k){
	
	return		 var[INDEX(i-3,j,k)] -
		 	 6.0*var[INDEX(i-2,j,k)] +
			15.0*var[INDEX(i-1,j,k)] -
			20.0*var[INDEX(i  ,j,k)] +
			15.0*var[INDEX(i+1,j,k)] -
			 6.0*var[INDEX(i+2,j,k)] +
				 var[INDEX(i+3,j,k)];
}

/*********************************************************************
*
*
**********************************************************************/
double diffuse_j_6th(double *var,int i,int j,int k){
	
	return		 var[INDEX(i,j-3,k)] -
		 	 6.0*var[INDEX(i,j-2,k)] +
			15.0*var[INDEX(i,j-1,k)] -
			20.0*var[INDEX(i,j  ,k)] +
			15.0*var[INDEX(i,j+1,k)] -
			 6.0*var[INDEX(i,j+2,k)] +
				 var[INDEX(i,j+3,k)];
}

/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
double diffuse_k_6th(double *var,int i,int j,int k){
	
	return		 var[INDEX(i,j,k-3)] -
		 	 6.0*var[INDEX(i,j,k-2)] +
			15.0*var[INDEX(i,j,k-1)] -
			20.0*var[INDEX(i,j,k  )] +
			15.0*var[INDEX(i,j,k+1)] -
			 6.0*var[INDEX(i,j,k+2)] +
				 var[INDEX(i,j,k+3)];
}

/*********************************************************************
*
*
**********************************************************************/
double diffuse_k_6th(double *var,double *weights,int i,int j,int k){
	
	return	weights[0] * var[INDEX(i,j,k-3)] -
		 	weights[1] * var[INDEX(i,j,k-2)] +
			weights[2] * var[INDEX(i,j,k-1)] -
			weights[3] * var[INDEX(i,j,k  )] +
			weights[4] * var[INDEX(i,j,k+1)] -
			weights[5] * var[INDEX(i,j,k+2)] +
			weights[6] * var[INDEX(i,j,k+3)];
}

/*********************************************************************
*
*
**********************************************************************/
double diffuse_i_8th(double *var,int i,int j,int k){
	
	return		 var[INDEX(i-4,j,k)] -
		 	 8.0*var[INDEX(i-3,j,k)] +
			28.0*var[INDEX(i-2,j,k)] -
			56.0*var[INDEX(i-1,j,k)] +
			70.0*var[INDEX(i  ,j,k)] -
			56.0*var[INDEX(i+1,j,k)] +
			28.0*var[INDEX(i+2,j,k)] -
			 8.0*var[INDEX(i+3,j,k)] +
				 var[INDEX(i+4,j,k)];
}

/*********************************************************************
*
*
**********************************************************************/
double diffuse_j_8th(double *var,int i,int j,int k){
	
	return		 var[INDEX(i,j-4,k)] -
		 	 8.0*var[INDEX(i,j-3,k)] +
			28.0*var[INDEX(i,j-2,k)] -
			56.0*var[INDEX(i,j-1,k)] +
			70.0*var[INDEX(i,j  ,k)] -
			56.0*var[INDEX(i,j+1,k)] +
			28.0*var[INDEX(i,j+2,k)] -
			 8.0*var[INDEX(i,j+3,k)] +
				 var[INDEX(i,j+4,k)];
}

/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
double diffuse_k_8th(double *var,int i,int j,int k){
	
	return		 var[INDEX(i,j,k-4)] -
		 	 8.0*var[INDEX(i,j,k-3)] +
			28.0*var[INDEX(i,j,k-2)] -
			56.0*var[INDEX(i,j,k-1)] +
			70.0*var[INDEX(i,j,k  )] -
			56.0*var[INDEX(i,j,k+1)] +
			28.0*var[INDEX(i,j,k+2)] -
			 8.0*var[INDEX(i,j,k+3)] +
				 var[INDEX(i,j,k+4)];
}


/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
void diffusion_4th_var(double *pvar,double *mvar,double step,int il,int ih,int jl,int jh,double *output){
	
	double diff_tend_var;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
		//-----------------------------------------------------
		// Use 2nd order vertical for lowest level
		//-----------------------------------------------------
		diff_tend_var = (	
			- kmixh4th * diffuse_i_4th(mvar,i,j,1) * one_d_dx4
			- kmixh4th * diffuse_j_4th(mvar,i,j,1) * one_d_dy4 	
			+ kmixv2nd * diffuse_k_2nd(mvar,i,j,1) * one_d_dz2
		);
		
		pvar[INDEX(i,j,1)] += step*dt*diff_tend_var;
		
		#if OUTPUT_DIFFUSION_TEND
		if(output != NULL)
			output[INDEX(i,j,1)] = diff_tend_var;
		#endif
		//-----------------------------------------------------
		// Middle levels
		//-----------------------------------------------------
		for(int k=2;k<NZ-2;k++){
		
			diff_tend_var = (	
				- kmixh4th * diffuse_i_4th(mvar,i,j,k) * one_d_dx4
				- kmixh4th * diffuse_j_4th(mvar,i,j,k) * one_d_dy4 	
				- kmixv4th * diffuse_k_4th(mvar,i,j,k) * one_d_dz4
			);
			
			pvar[INDEX(i,j,k)] += step*dt*diff_tend_var;
		
			#if OUTPUT_DIFFUSION_TEND
			if(output != NULL)
				output[INDEX(i,j,k)] = diff_tend_var;
			#endif
		}
		//-----------------------------------------------------
		// Use 2nd order vertical for lhighest level
		//-----------------------------------------------------
		diff_tend_var = (	
			- kmixh4th * diffuse_i_4th(mvar,i,j,NZ-2) * one_d_dx4
			- kmixh4th * diffuse_j_4th(mvar,i,j,NZ-2) * one_d_dy4 	
			+ kmixv2nd * diffuse_k_2nd(mvar,i,j,NZ-2) * one_d_dz2
		);
		
		pvar[INDEX(i,j,NZ-2)] += step*dt*diff_tend_var;
	
		#if OUTPUT_DIFFUSION_TEND
		if(output != NULL)
			output[INDEX(i,j,NZ-2)] = diff_tend_var;
		#endif
	}}
}

/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
void diffusion_6th_var(double *pvar,double *mvar,double step,int il,int ih,int jl,int jh,double *output){
	
	double diff_tend_var;
	double xdiff,ydiff,zdiff;
	int ind = DIFFUSION_ORDER + 1; // should be seven, right?
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		
		//-----------------------------------------------------
		// Use lower order vertical diffusion for uppermost
		// and lowermost heights
		//-----------------------------------------------------
#if STRETCHED_GRID
		//printf("%d ",index2d(ind,k,0));
		if(		k==1 || k==NZ-2){ zdiff = +kmixv2nd * diffuse_k_2nd(mvar,weights+index2d(ind,k,0),i,j,k);} 
		else if(k==2 || k==NZ-3){ zdiff = -kmixv4th * diffuse_k_4th(mvar,weights+index2d(ind,k,0),i,j,k);}
		else 					{ zdiff = +kmixv6th * diffuse_k_6th(mvar,weights+index2d(ind,k,0),i,j,k);}
#else
		if(		k==1 || k==NZ-2){ zdiff = +kmixv2nd * diffuse_k_2nd(mvar,i,j,k) * one_d_dz2;} 
		else if(k==2 || k==NZ-3){ zdiff = -kmixv4th * diffuse_k_4th(mvar,i,j,k) * one_d_dz4;}
		else 					{ zdiff = +kmixv6th * diffuse_k_6th(mvar,i,j,k) * one_d_dz6;}	
#endif		
		
#if MERIDIONAL_CROSS_SECTION
		xdiff = 0;
#else
		xdiff = kmixh6th * diffuse_i_6th(mvar,i,j,k) * one_d_dx6;
#endif
		ydiff = kmixh6th * diffuse_j_6th(mvar,i,j,k) * one_d_dy6;

		diff_tend_var = xdiff + ydiff + zdiff;
		
		temp_var[INDEX(i,j,k)] = step*dt*diff_tend_var;
	
		#if OUTPUT_DIFFUSION_TEND
		if(output != NULL)
			output[INDEX(i,j,k)] = diff_tend_var;
		#endif
	}}}
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		
		pvar[INDEX(i,j,k)] += temp_var[INDEX(i,j,k)];
		
	}}}
}

/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
void diffusion_8th_var(double *pvar,double *mvar,double step,int il,int ih,int jl,int jh,double *output){
	
	double diff_tend_var;
	double xdiff,ydiff,zdiff;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		
		//-----------------------------------------------------
		// Use lower order vertical diffusion for uppermost
		// and lowermost heights
		//-----------------------------------------------------
		if(		k==1 || k==NZ-2){ zdiff = +kmixv2nd * diffuse_k_2nd(mvar,i,j,k) * one_d_dz2;} 
		else if(k==2 || k==NZ-3){ zdiff = -kmixv4th * diffuse_k_4th(mvar,i,j,k) * one_d_dz4;}
		else 					{ zdiff = +kmixv6th * diffuse_k_6th(mvar,i,j,k) * one_d_dz6;}
		
		xdiff = -kmixh8th * diffuse_i_8th(mvar,i,j,k) * one_d_dx8;
		ydiff = -kmixh8th * diffuse_j_8th(mvar,i,j,k) * one_d_dy8;

		diff_tend_var = xdiff + ydiff + zdiff;
		
		pvar[INDEX(i,j,k)] += step*dt*diff_tend_var;
	
		#if OUTPUT_DIFFUSION_TEND
		if(output != NULL)
			output[INDEX(i,j,k)] = diff_tend_var;
		#endif
	}}}
}


/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
void diffusion_2nd_var(double *pvar,double *mvar,double step,int il,int ih,int jl,int jh){

	//-----------------------------------------------------
	// Middle levels
	//-----------------------------------------------------	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		
		pvar[INDEX(i,j,k)] += step*dt* (	
			+ kmixh2nd * diffuse_i_2nd(mvar,i,j,k) * one_d_dx2
			+ kmixh2nd * diffuse_j_2nd(mvar,i,j,k) * one_d_dy2 	
			+ kmixv2nd * diffuse_k_2nd(mvar,i,j,k) * one_d_dz2
			);
	}}}
}

/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
void diffusion_8th_all(double step,int il,int ih,int jl,int jh){
	
	if(OUTPUT_DIFFUSION_TEND){
		diffusion_8th_var(&UP(0,0,0),&UM(0,0,0),step,il,ih,jl,jh,u_diff_tend);
		diffusion_8th_var(&VP(0,0,0),&VM(0,0,0),step,il,ih,jl,jh,v_diff_tend);
	} else {
		diffusion_8th_var(&UP(0,0,0),&UM(0,0,0),step,il,ih,jl,jh);
		diffusion_8th_var(&VP(0,0,0),&VM(0,0,0),step,il,ih,jl,jh);
	}
	
	diffusion_8th_var(&WP(0,0,0),&WM(0,0,0),step,il,ih,jl,jh);
	diffusion_8th_var(&THP(0,0,0),&THM(0,0,0),step,il,ih,jl,jh);
}

/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
void diffusion_6th_all(double step,int il,int ih,int jl,int jh){
#if 1
	if(OUTPUT_DIFFUSION_TEND){
		diffusion_6th_var(&UM(0,0,0),&UM(0,0,0),step,il,ih,jl,jh,u_diff_tend);
		diffusion_6th_var(&VM(0,0,0),&VM(0,0,0),step,il,ih,jl,jh,v_diff_tend);
	} else {
		diffusion_6th_var(&UM(0,0,0),&UM(0,0,0),step,il,ih,jl,jh);
		diffusion_6th_var(&VM(0,0,0),&VM(0,0,0),step,il,ih,jl,jh);
	}
	
	if(!HYDROSTATIC){
		diffusion_6th_var(&WM(0,0,0),&WM(0,0,0),step,il,ih,jl,jh);
	}
	
	diffusion_6th_var(&THM(0,0,0),&THM(0,0,0),step,il,ih,jl,jh);
	
	if(USE_MICROPHYSICS)
		
		diffusion_6th_var(&QVM(0,0,0),&QVM(0,0,0),step,il,ih,jl,jh);
		diffusion_6th_var(&QCM(0,0,0),&QCM(0,0,0),step,il,ih,jl,jh);
		diffusion_6th_var(&QRM(0,0,0),&QRM(0,0,0),step,il,ih,jl,jh);
		
		if(USE_ICE){
			diffusion_6th_var(&QIM(0,0,0),&QIM(0,0,0),step,il,ih,jl,jh);
			diffusion_6th_var(&QSM(0,0,0),&QSM(0,0,0),step,il,ih,jl,jh);
		}
	#endif
}

/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
void diffusion_6th_all2(double step,int il,int ih,int jl,int jh){
	

	diffusion_6th_var(&UBAR(0,0,0),&UBAR(0,0,0),step,il,ih,jl,jh);
	diffusion_6th_var(&VBAR(0,0,0),&VBAR(0,0,0),step,il,ih,jl,jh);
	diffusion_6th_var( &WBAR(0,0,0), &WBAR(0,0,0),step,il,ih,jl,jh);
	diffusion_6th_var(&THBAR(0,0,0),&THBAR(0,0,0),step,il,ih,jl,jh);
}

/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
void diffusion_4th_all(double step,int il,int ih,int jl,int jh){
	
	if(OUTPUT_DIFFUSION_TEND){
		diffusion_4th_var(&UP(0,0,0),&UM(0,0,0),step,il,ih,jl,jh,u_diff_tend);
		diffusion_4th_var(&VP(0,0,0),&VM(0,0,0),step,il,ih,jl,jh,v_diff_tend);
	} else {
		diffusion_4th_var(&UP(0,0,0),&UM(0,0,0),step,il,ih,jl,jh);
		diffusion_4th_var(&VP(0,0,0),&VM(0,0,0),step,il,ih,jl,jh);
	}
	
	diffusion_4th_var(&WP(0,0,0),&WM(0,0,0),step,il,ih,jl,jh);
	diffusion_4th_var(&THP(0,0,0),&THM(0,0,0),step,il,ih,jl,jh);
}

/*********************************************************************
* NOT YET COMPATIBLE WITH STRETCHED GRID!!!!!!!!!
*
**********************************************************************/
void diffusion_2nd_all(double step,int il,int ih,int jl,int jh){
	
	diffusion_2nd_var(&UP(0,0,0),&UM(0,0,0),step,il,ih,jl,jh);
	diffusion_2nd_var(&VP(0,0,0),&VM(0,0,0),step,il,ih,jl,jh);
	diffusion_2nd_var(&WP(0,0,0),&WM(0,0,0),step,il,ih,jl,jh);
	diffusion_2nd_var(&THP(0,0,0),&THM(0,0,0),step,il,ih,jl,jh);
}

/*********************************************************************
* Diffuse microphysics variables
*
**********************************************************************/
void microphysics_diffusion(int il,int ih,int jl,int jh,double step){
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		DIFFUSE_VARIABLE(QVP,QVM);
		DIFFUSE_VARIABLE(QCP,QCM);	
		DIFFUSE_VARIABLE(QRP,QRM);
		
	}}}
}

/*********************************************************************
* Diffuse variables
*
**********************************************************************/
void apply_explicit_diffusion(int il,int ih,int jl,int jh,double step){
	
	diffusion_6th_all(step,3,NX-3,3,NY-3);
	
}

/*********************************************************************
*
*
**********************************************************************/
void rayleigh_damping(int il,int ih,int jl, int jh,int level){

double coefu;
double coefw;
double raydmpcoef = 0.05;
double raydmpz = zu[level];

	// begin at the lowest level at which the damping is to be applied
	for(int k=level;k<NZ;k++){
 	
		// calculate proper coefficients at the levels of the scalars and w-velocity
	 	coefu = raydmpcoef * 0.5 * (1.0-cos(trigpi*(zu[k]-raydmpz)/(zu[NZ-1]-raydmpz)));
		coefw = raydmpcoef * 0.5 * (1.0-cos(trigpi*(zw[k]-raydmpz)/(zw[NZ-1]-raydmpz)));
		
		for(int i=il;i<ih;i++){
		for(int j=jl;j<jh;j++){
	 
	 		UP(i,j,k) -= coefu*UP(i,j,k);
			VP(i,j,k) -= coefu*VP(i,j,k);
			
		#if HYDROSTATIC
			W(i,j,k) -= coefw*W(i,j,k);
		#else
			WP(i,j,k) -= coefw*WP(i,j,k);
		#endif
			
			THP(i,j,k) -= coefu*THP(i,j,k);
		}}
	}
}
#if 0
/*********************************************************************
*
*
**********************************************************************/
void init_large_scale_precip(){
	
	for(int k=0;k<NZ;k++){
		
		if(zu[k]<Z_Tropopause && zu[k]>0){
		
			R[k] = R0 + (1.0-R0)*( zu[k] / Z_Tropopause);
						
		} else {
			
			R[k] = 1;
		}
		
		printf("%d %f %f\n",k,zu[k],R[k]);
	}
	
}

/*********************************************************************
*
*
**********************************************************************/
void large_scale_precip(int il,int ih,int jl,int jh){
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		
		THP(i,j,k) += 0.004*(1.0-R[k])*WP(i,j,k) * dt;
		
	}}}
	
}
#endif
/*********************************************************************
* 2 K / day = 2.3148e-5
* 10 days = 1.157e-6
**********************************************************************/
void damp_var(double *var,int il,int ih,int jl, int jh,double coef,double max){

    double rate;
    
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ;k++){
		
		rate = coef * var[INDEX(i,j,k)];
		
		if(rate > max){ rate = max;}
        else if(rate < -max){ rate = -max;}
	
		var[INDEX(i,j,k)] = var[INDEX(i,j,k)] - rate * dt;
		
	}}}
	
}

