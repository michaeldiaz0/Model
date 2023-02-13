#include "stdafx.h"
#include "pcomm.h"

bool startOutput = false;

double *smoothVar,*smoothVarO;

/***********************************************************************************************
* -------------------------------------- MACROS ------------------------------------------------
************************************************************************************************/
#if !PARALLEL
	#define VAR(i,j,k) var[NZ*(NY*(i)+(j))+(k)]
	//int rank = 0;
#else
	#define VAR(i,j,k) var[fNZ*(fNY*(i)+(j))+(k)]
#endif

double get_min(int klev,double *var,int *ival,int *jval,int il,int ih,int jl,int jh);
static int index(int,int,int);
void smooth_var(double *var,double *svar,int klev,int times);
void smooth9(double *varIn,double *varOut,int il,int ih,int jl,int jh,int lNY);
double get_min_2D(double *var,int *ival,int *jval,int il,int ih,int jl,int jh);

/*********************************************************************
* Returns minimum value
*
**********************************************************************/
int min(int x,int y){
	
	if(x<y){ return x;}
	
	return y;
}

/****************************************************
! Bengt Fornberg (1988) algorithm for finding interpolation
! weights for derivatives
!
! Input Parameters
! z       - location where approximations are to be accurate,
! x(0:nd) - grid point locations, found in x(0:n)
! n       - one less than total number of grid points; n must not exceed the parameter nd below,
! nd      - dimension of x- and c-arrays in calling program x(0:nd) and c(0:nd,0:m), respectively, 
! m       - highest derivative for which weights are sought,
! 
! Output Parameter
! c(0:nd,0:m) weights at grid locations x(0:n) for derivatives c--- of order 0:m, found in c(0:n,0:m)
****************************************************/
void derivative_weights(double z,double *x, int n, int nd, int m,double *c){

    double c1,c2,c3,c4,c5;
    int  i, j, k, mn;

#if 0
	printf("\n");
	printf("%f\n\n",z);
	for(int i=0;i<=n;i++){ printf("%f\n",x[i]);}
	printf("\n");
#endif
	
    c1 = 1.0;
    c4 = x[0] - z;

    c[index2d(m+1,0,0)] = 1.0;

	for(int i=1;i<=n;i++){
		
        mn = min(i,m);
        c2 = 1.0;
        c5 = c4;
        c4 = x[i]-z;

		for(int j=0;j<=i-1;j++){

            c3 = x[i]-x[j];
            c2 = c2*c3;

            if (j == i-1){
				
				for(int k=mn;k>=1;k--){
					
                    c[index2d(m+1,i,k)] = c1*( k * c[index2d(m+1,i-1,k-1)] - c5 * c[index2d(m+1,i-1,k)] ) / c2;
				}

                c[index2d(m+1,i,0)] = -c1*c5*c[index2d(m+1,i-1,0)] / c2;

            }

            for(int k=mn;k>=1;k--){
				
                c[index2d(m+1,j,k)] = (c4 * c[index2d(m+1,j,k)] - k * c[index2d(m+1,j,k-1)] ) / c3;
            }

            c[index2d(m+1,j,0)]  = c4 * c[index2d(m+1,j,0)] / c3;
        }

        c1 = c2;

    }
	
	for(int j=0;j<=m;j++){
		//printf("%d\n",j);
	
		for(int i=0;i<=m;i++){
			//printf("%f  \n",c[index2d(m+1,i,j)] );
	}}
	

}


/****************************************************
* 
*****************************************************/
void init_stats(){
	
	smoothVar = (double*) calloc(NX*NY,sizeof(double));
	smoothVarO = (double*) calloc(NX*NY,sizeof(double));
}

/****************************************************
* 
*****************************************************/
double get_stats(FILE *infile){
	
	int ipos,jpos,ipos2,jpos2;
	double pmin,pmin2,thmin,pd;
	double cpRd = cp/Rd;

	int il,ih,jl,jh;

	if(infile==NULL){
		printf("time %0.3f\n",mtime/3600);
        fflush(stdout);
	} else {
		fprintf(infile,"time %d %0.3f\n",bigcounter,mtime/3600);
        fflush(stdout);
	}
	
	//---------------------------------------------------------
	// Find the lowest pressure at a single reference height
	// to be used to track the pressure minimum. Ensures that
	// we are tracking a single low pressure system.
	//---------------------------------------------------------
	const int tracking_height = 15;
	const double box_size = 500000; // in meters
	ipos = 0; jpos = 0;

	smooth_var(pis,smoothVar,tracking_height,20);
	get_min_2D(smoothVar,&ipos,&jpos,0,NX,0,NY);

	const int radius = (int)(box_size / dx);

	il = ipos - radius;
	ih = ipos + radius;
	jl = jpos - radius;
	jh = jpos + radius;
	
	if(il < 0) { il = 0;}
	if(ih > NX){ ih = NX;}
	if(jl < 0) { jl = 0;}
	if(jh > NY){ jh = NY;}

	//---------------------------------------------------------
	// For each vertical level...
	//---------------------------------------------------------
	for(int k=1;k<NZ-1;k++){
		
		ipos = 0; jpos = 0; ipos2 = 0; jpos2 = 0;
		//---------------------------------------------------------
		// Find location and magnitude of lowest pressure
		//---------------------------------------------------------		
		pmin = get_min(k,pis,&ipos,&jpos,il,ih,jl,jh);
		//---------------------------------------------------------
		// Find location and magnitude of lowest pressure for a
		// smoothed field
		//---------------------------------------------------------			
		smooth_var(pis,smoothVar,k,10);
		pmin2 = get_min_2D(smoothVar,&ipos2,&jpos2,il,ih,jl,jh);
		//---------------------------------------------------------
		// Get smoothed temperature at minimum pressure for smoothed
		// pressure fields
		//---------------------------------------------------------		
		smooth_var(ths,smoothVar,k,10);
		thmin = smoothVar[index(ipos2,jpos2,NY)];
		//---------------------------------------------------------
		// Calculate pressure in hPa
		//---------------------------------------------------------		
		pd = p0*pow(PI(ipos,jpos,k)/(cp*tbv[k]) + PBAR(ipos,jpos,k),cpRd);
		//---------------------------------------------------------
		// Output results
		//---------------------------------------------------------				
		if(infile==NULL){
			printf("%d %f %f %f %f %f %f %f %f %f %f\n",k,zsu[k]*0.001,pd*0.01,
			outLons[ipos],outLats[jpos],pmin,TH(ipos,jpos,k),outLons[ipos2],outLats[jpos2],pmin2,thmin);
		} else {
			fprintf(infile,"%d %f %f %f %f %f %f %f %f %f %f\n",k,zsu[k]*0.001,pd*0.01,
			outLons[ipos],outLats[jpos],pmin,TH(ipos,jpos,k),outLons[ipos2],outLats[jpos2],pmin2,thmin);
			fflush(stdout);
		}
	}
	
	return 0;
}

/****************************************************
* 
*****************************************************/
double get_min(int klev,double *var,int *ival,int *jval,int il,int ih,int jl,int jh){
	
	double min = 3000;
	double ipos = 0;
	double jpos = 0;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
		
		if( VAR(i,j,klev) < min ){
			
			min = VAR(i,j,klev);
			ipos = i;
			jpos = j;
		}
	}}
	
	*ival = ipos;
	*jval = jpos;
	
	return min;
}

/****************************************************
* 
*****************************************************/
double get_min_2D(double *var,int *ival,int *jval,int il,int ih,int jl,int jh){
	
	double min = 3000;
	double ipos = 0;
	double jpos = 0;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
		
		if( var[index(i,j,NY)] < min ){
			
			min = var[index(i,j,NY)];
			ipos = i;
			jpos = j;
		}
	}}
	
	*ival = ipos;
	*jval = jpos;
	
	return min;
}

/****************************************************
* 
*****************************************************/
void smooth_var(double *var,double *svar,int klev,int times){
	
	//---------------------------------------------------
	// Extract requested vertical level from input array
	//---------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	
		svar[index(i,j,NY)] = VAR(i,j,klev);
		
	}}
	
	double *temp;
	
	//---------------------------------------------------
	// Perform smoothing
	//---------------------------------------------------
	for(int i=0;i<times;i++){
		
		smooth9(svar,smoothVarO,1,NX-1,1,NY-1,NY);
		
		temp = smoothVarO;
		smoothVarO = svar;
		svar = temp;
	}
	
}

/****************************************************
* memset(smoothVarO,0,NX*NY*sizeof(double));
*****************************************************/
void smooth9(double *varIn,double *varOut,int il,int ih,int jl,int jh,int lNY){
	
	const double f = 8.0;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
		
		varOut[index(i,j,lNY)] = varIn[index(i+1,j+1,lNY)] +   varIn[index(i,j+1,lNY)] + varIn[index(i-1,j+1,lNY)] +
								 varIn[index(i+1,j  ,lNY)] + f*varIn[index(i,j  ,lNY)] + varIn[index(i-1,j  ,lNY)] +
								 varIn[index(i+1,j-1,lNY)] +   varIn[index(i,j-1,lNY)] + varIn[index(i-1,j-1,lNY)];
		
		varOut[index(i,j,lNY)] /= (f+8.0);
	}}
	
	
}

/****************************************************
* 
*****************************************************/
static int index(int i,int j,int lNY){
	return j+i*lNY;
}

/****************************************************
* 
*****************************************************/
double pivot_coriolis(double angle,int i,int j){
	
	double lat = outLats[j] + (dx/meters_per_degree) * angle/90.0 * (float)(i - NX/2);
	
	return 2*fc*sin(lat*trigpi/180.);
}

/****************************************************
* Finds the nearest grid point to a given latitude 
*
*****************************************************/
int get_point_from_lat(double lat){

	int point = 0;

	while(lat>outLats[point]){ point++;}

	if(point>=NY){
		
		printf("Error:latitude of %f not found!\n",lat);
		exit(0);
	}

	if( fabs(outLats[point-1]-lat) < fabs(outLats[point]-lat) ){

		return point-1;

	} else { return point;}
}

/****************************************************
* Finds the nearest grid point to a given longitude 
*
*****************************************************/
int get_point_from_lon(double lon){

	int point = 0;

	while(lon>outLons[point]){ point++;}

	if(point>=NX){
		
		printf("Error:longitude of %f not found!\n",lon);
		exit(0);
	}

	if( fabs(outLons[point-1]-lon) < fabs(outLons[point]-lon) ){

		return point-1;

	} else { return point;}
}

/****************************************************
* Finds the nearest grid point to a given height in meters
*
*****************************************************/
int get_point_from_height(double zlev){

	int point = 0;

	while(zlev>ZU(point)){ point++;}

	if(point>=NZ){
		
		printf("Error:height of %f not found!\n",zlev);
		exit(0);
	}

	if( fabs(ZU(point-1)-zlev) < fabs(ZU(point)-zlev) ){

		return point-1;

	} else { return point;}
}

/****************************************************
* 
*****************************************************/
void print_max_arrays(int s){
	
	double max_u=0,max_v=0,max_w=0,max_th=0,max_qv=0,max_qc=0,max_qr=0,max_k=0;
	int max_u_i = -1; int max_u_j = -1; int max_u_k = -1;
	int max_v_i = -1; int max_v_j = -1; int max_v_k = -1;
	int max_w_i = -1; int max_w_j = -1; int max_w_k = -1;
	int max_th_i = -1; int max_th_j = -1; int max_th_k = -1;
	int max_qv_i = -1; int max_qv_j = -1; int max_qv_k = -1;
	int max_qc_i = -1; int max_qc_j = -1; int max_qc_k = -1;
	int max_qr_i = -1; int max_qr_j = -1; int max_qr_k = -1;
	int max_k_i = -1; int max_k_j = -1; int max_k_k = -1;
	
	for(int i=3;i<fNX-3;i++){
	for(int j=3;j<fNY-3;j++){
	for(int k=1;k<fNZ-1;k++){

		if(fabs(U(i,j,k))  >= max_u){  max_u  = fabs(U(i,j,k));   max_u_i = i;  max_u_j = j;  max_u_k = k; }
		if(fabs(V(i,j,k))  >= max_v){  max_v  = fabs(V(i,j,k));   max_v_i = i;  max_v_j = j;  max_v_k = k; }
		if(fabs(W(i,j,k))  >= max_w){  max_w  = fabs(W(i,j,k));   max_w_i = i;  max_w_j = j;  max_w_k = k; }
		if(fabs(TH(i,j,k)) >= max_th){ max_th = fabs(TH(i,j,k));  max_th_i = i; max_th_j = j; max_th_k = k; }
		if(fabs(QV(i,j,k)) >= max_qv){ max_qv = fabs(QV(i,j,k));  max_qv_i = i; max_qv_j = j; max_qv_k = k; }
		if(fabs(QC(i,j,k)) >= max_qc){ max_qc = fabs(QC(i,j,k));  max_qc_i = i; max_qc_j = j; max_qc_k = k; }		
		if(fabs(QR(i,j,k)) >= max_qr){ max_qr = fabs(QR(i,j,k));  max_qr_i = i; max_qr_j = j; max_qr_k = k; }
		if(fabs(FRICTION(i,j,k)) >= max_k){ max_k = FRICTION(i,j,k);  max_k_i = i; max_k_j = j; max_k_k = k; }

	}}}

	if(max_qr>0.010){
		
		startOutput = true;
		
		printf("rain : %f %f %f %f %f %f %d %d %d %d %d %f %f %f %f %f %f %f\n",
			QR(max_qr_i+1,max_qr_j,max_qr_k),
			QR(max_qr_i-1,max_qr_j,max_qr_k),
			QR(max_qr_i,max_qr_j+1,max_qr_k),
			QR(max_qr_i,max_qr_j-1,max_qr_k),
			QR(max_qr_i,max_qr_j,max_qr_k+1),
			QR(max_qr_i,max_qr_j,max_qr_k-1),
			HTOPO(max_qr_i,max_qr_j),
			HTOPO(max_qr_i+1,max_qr_j),
			HTOPO(max_qr_i-1,max_qr_j),
			HTOPO(max_qr_i,max_qr_j+1),
			HTOPO(max_qr_i,max_qr_j-1),
			FRICTION(max_qr_i,max_qr_j,max_qr_k+1),
			FRICTION(max_qr_i+1,max_qr_j,max_qr_k+1),
			FRICTION(max_qr_i-1,max_qr_j,max_qr_k+1),
			FRICTION(max_qr_i,max_qr_j+1,max_qr_k+1),
			FRICTION(max_qr_i,max_qr_j-1,max_qr_k+1),
			FRICTION(max_qr_i,max_qr_j,max_qr_k+1+1),
			FRICTION(max_qr_i,max_qr_j,max_qr_k-1+1)
		);
	} else { startOutput = false;}

	printf("\nAt substep %d of rank %d\n",s,rank);
	printf("At %d %d %d of rank %d u_max = %f || ",max_u_i, max_u_j, max_u_k,rank,max_u);
	printf("At %d %d %d of rank %d v_max = %f || ",max_v_i, max_v_j, max_v_k,rank,max_v);
	printf("At %d %d %d of rank %d w_max = %f\n",max_w_i, max_w_j, max_w_k,rank,max_w);
	printf("At %d %d %d of rank %d qv_max = %f || ",max_qv_i, max_qv_j, max_qv_k,rank,max_qv*1000.0);		
	printf("At %d %d %d of rank %d qc_max = %f || ",max_qc_i, max_qc_j, max_qc_k,rank,max_qc*1000.0);		
	printf("At %d %d %d of rank %d qr_max = %f\n",max_qr_i, max_qr_j, max_qr_k,rank,max_qr*1000.0);
	printf("At %d %d %d of rank %d th_max = %f || ",max_th_i, max_th_j, max_th_k,rank,max_th);
	printf("At %d %d %d of rank %d k_max = %f\n",max_k_i, max_k_j, max_k_k,rank,max_k);		

	fflush(stdout);
	
}

#if 0
/*********************************************************************
*
*
**********************************************************************/
void debug(int a, int b, int c,double uvar[NX][NY][NZ],double vvar[NX][NY][NZ],double var[NX][NY][NZ]){

	int columns = 6;

	printf("     ");

	for(int i=0;i<columns;i++){printf("|                         ");}  printf("\n");

	printf("     ");

	for(int i=0;i<columns;i++){ printf("|        (%d,%d)          ",a+i,b+1);} printf("\n");

	printf("     ");
	
	for(int i=0;i<columns;i++){ printf(" ------ %+f  -------",vvar[a+i][b+1][c]);} printf("\n");

	printf("     ");

	for(int i=0;i<columns;i++){printf("|                         ");}  printf("\n");

	printf("  ");

	for(int i=0;i<columns;i++){printf("(%d,%d)     (%d,%d)       ",a+i,b,a+i,b);}  printf("\n");
	for(int i=0;i<columns;i++){ printf("%+f    %+f    ",uvar[a+i][b][c],var[a+i][b][c]);} printf("\n");


}

/*********************************************************************
*
*
**********************************************************************/
void zdebug(int a, int b, int c,double uvar[NX][NY][NZ],double vvar[NX][NY][NZ],double var[NX][NY][NZ]){

	int columns = 6;

	printf("     ");

	for(int i=0;i<columns;i++){printf("|                         ");}  printf("\n");

	printf("     ");

	for(int i=0;i<columns;i++){ printf("|       (%3d,%3d)         ",b+i,c+1);} printf("\n");

	printf("     ");
	
	for(int i=0;i<columns;i++){ printf(" ------ %+f  -------",vvar[a][b+i][c+1]);} printf("\n");

	printf("     ");

	for(int i=0;i<columns;i++){printf("|                         ");}  printf("\n");

	printf(" ");

	for(int i=0;i<columns;i++){printf("(%3d,%3d)   (%3d,%3d)     ",b+i,c,b+i,c);}  printf("\n");
	for(int i=0;i<columns;i++){ printf("%+f    %+f    ",uvar[a][b+i][c],var[a][b+i][c]);} printf("\n");


}
#endif
/*********************************************************************
* 
**********************************************************************/
double convert_z_to_k(double zf,double z0,double lev,double dzp){
	
	int nt = NZ-1;
	
	double c2 = ( 1 - z0/(dzp*lev)) / ( (nt-1) *dzp - dzp*lev); 
	
	double c1 = z0 / (dzp*lev) - c2*dzp*lev;

	double zi = (-c1 + sqrt(c1*c1+4*c2*zf)) / (2.0*c2);
	
	return zi/dzp + 0.5;
}

/*********************************************************************
* 
**********************************************************************/
double convert_k_to_z(double ki,double z0,double lev,double dzp){
	
	int nt = NZ-1;
	
	double c2 = ( 1 - z0/(dzp*lev)) / ( (nt-1) *dzp - dzp*lev); 
	
	double c1 = z0 / (dzp*lev) - c2*dzp*lev;

	double halfLev = (ki-0.5)*dzp;
		
	return (c1+c2*halfLev)*halfLev;
}


/*********************************************************************
* 
*
**********************************************************************/
double stable_time_step_5th(double delta,double vel){
	
	return 0.82 * delta / vel;
}

/*********************************************************************
* 
*
**********************************************************************/
double stable_time_step(double delta,double vel){
	
	return 0.9296 * delta / vel;
}

/*********************************************************************
* 
*
* @param pl - start of interval in periods
* @param ph - end of interval in periods
* @param frac - fractional distance between pl and ph
* @param shift - shift entire function up or down in units of amplitude
* @return - a number from -one to one
**********************************************************************/
double cos_profile(double pl,double ph,double frac,double shift){
	
	if(frac >= 0 && frac <= 1.0){
	
		return cos(  pl * 2*trigpi + (ph-pl)*2*trigpi*frac ) + shift;
		
	} else {
		return 0;
	}
}


/*********************************************************************
* 
**********************************************************************/
double frac_distance(int index,int low,int high){
	
	return (double)(index-low)/(double)(high-low);
}

/*********************************************************************
* 
**********************************************************************/
double frac_distance(double index,double low,double high){
	
	return (index-low)/(high-low);
}

/*********************************************************************
*
*
**********************************************************************/
void line_array(int x1,int y1,int x2,int y2,int * array){

	float slope = (float)(y2-y1)/(float)(x2-x1);

	for(int i=0;i<x1;i++)
		array[i]=y1;

	for(int i=x2+1;i<NX;i++)
		array[i]=y2;

	for(int i=x1;i<=x2;i++)
		array[i] = y1 + (int)(slope*(float)(i-x1));

}
#if 0
/*********************************************************************
*
*
**********************************************************************/
void find_max(double var[NX][NY][NZ], int * arr, double * max, double * max2){

	double mymax = 0;
	double mymax2 = 0;

	for(int i=0;i<NX;i++)
		for(int j=arr[i];j<NY;j++)
			for(int k=0;k<NZ;k++)
				if(V(i,j,k) > fabs(mymax2))
					mymax2 = fabs(var[i][j][k]);


	for(int i=0;i<NX;i++)
		for(int j=0;j<arr[i];j++)
			for(int k=0;k<NZ;k++)
				if(V(i,j,k) > fabs(mymax))			
					mymax = fabs(var[i][j][k]);

	
	*max = mymax;
	*max2 = mymax2;
}
#endif
/*********************************************************************
*
*
**********************************************************************/
double find_max2(double * var){

	double mymax = 0;

	for(int i=0;i<NX;i++)
		for(int j=0;j<NY;j++)
			for(int k=0;k<NZ;k++)
				if(var[INDEX(i,j,k)]*rhou[k] > fabs(mymax))
					mymax = fabs(var[INDEX(i,j,k)]*rhou[k]);

	return mymax;
}
#if 0
/*********************************************************************
*
*
**********************************************************************/
int find_max3(double var[NX][NY][NZ],int i_old){

	double mymax = 0;
	int max_i = 0;
	int k = 7;
	int j = 25;
	
	int start,end;
	
	start = i_old - NX/4;
	end = i_old+NX/4;
	
	if(start<0){ start = NX-NX/4; end = NX;}
	if(end>NX){ end = NX;}
	//printf("%d %d\n",start,end);
	//if(start > 3 && end < NX-3){

		for(int i=start;i<end;i++){
			//for(int j=0;j<NY;j++){
				if(var[i][j][k] < mymax){
					max_i = i;
					mymax = var[i][j][k];
				}
				//}
		}
		
		return max_i;
		
		//} else {
		//return NX/2;	
		//}

	
}
#endif
/*********************************************************************
*
*
**********************************************************************/
void rescale_pert(double max1, double max2, double scale1, double scale2, int * arr, bool upper, bool lower){

	if(lower){

		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
		for(int k=0;k<NZ;k++){

			U(i,j,k) = scale1*U(i,j,k) / max1;
			V(i,j,k) = scale1*V(i,j,k) / max1;
			W(i,j,k) = scale1*W(i,j,k) / max1;
			TH(i,j,k) = scale1*TH(i,j,k) / max1;

			UM(i,j,k) = scale1*UM(i,j,k) / max1;
			VM(i,j,k) = scale1*VM(i,j,k) / max1;
			WM(i,j,k) = scale1*WM(i,j,k) / max1;
			THM(i,j,k) = scale1*THM(i,j,k) / max1;

		}}}
	}

	if(upper){

		for(int i=0;i<40;i++){
		for(int j=arr[i];j<NY;j++){
		for(int k=0;k<NZ;k++){

			U(i,j,k) = scale2*U(i,j,k) / max2;
			V(i,j,k) = scale2*V(i,j,k) / max2;
			W(i,j,k) = scale2*W(i,j,k) / max2;
			TH(i,j,k) = scale2*TH(i,j,k) / max2;

			UM(i,j,k) = scale2*UM(i,j,k) / max2;
			VM(i,j,k) = scale2*VM(i,j,k) / max2;
			WM(i,j,k) = scale2*WM(i,j,k) / max2;
			THM(i,j,k) = scale2*THM(i,j,k) / max2;

		}}}
	}

}

/*********************************************************************
*
*
**********************************************************************/
void rescale_pert2(double max,double scale){

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		U(i,j,k) = scale*U(i,j,k) / max;
		V(i,j,k) = scale*V(i,j,k) / max;
		W(i,j,k) = scale*W(i,j,k) / max;
		TH(i,j,k) = scale*TH(i,j,k) / max;
		

		UM(i,j,k) = scale*UM(i,j,k) / max;
		VM(i,j,k) = scale*VM(i,j,k) / max;
		WM(i,j,k) = scale*WM(i,j,k) / max;
		THM(i,j,k) = scale*THM(i,j,k) / max;
		
		
#if USE_MICROPHYSICS
		QV(i,j,k) = scale*QV(i,j,k) / max;
		QVM(i,j,k) = scale*QVM(i,j,k) / max;
#endif

	}}}
}

/*********************************************************************
* 
*
**********************************************************************/
void print_time_estimates(double total_cputime,double total_walltime,int timer_counter){

	int cputime_secs = (int)((total_cputime  / (double)timer_counter) * (number_of_time_steps - bigcounter));
	int walltime_secs = (int)((total_walltime / (double)timer_counter) * (number_of_time_steps - bigcounter));

	int cpu_time_hours = cputime_secs / 3600;
	int cpu_time_min = (cputime_secs % 3600) / 60;
	int cpu_time_sec = (cputime_secs % 3600) % 60;
	int wall_time_hours = walltime_secs / 3600;
	int wall_time_min = (walltime_secs % 3600) / 60;
	int wall_time_sec = (walltime_secs % 3600) % 60;

	printf("Estimated time remaining: cpu time = %02d:%02d:%02d h:m:s, wall time = %02d:%02d:%02d hours\n",
	cpu_time_hours,cpu_time_min,cpu_time_sec,wall_time_hours,wall_time_min,wall_time_sec);

	cputime_secs = (int)((total_cputime  / (double)timer_counter) * number_of_time_steps);
	walltime_secs = (int)((total_walltime / (double)timer_counter) * number_of_time_steps);
			
	cpu_time_hours = cputime_secs / 3600;
	cpu_time_min = (cputime_secs % 3600) / 60;
	wall_time_hours = walltime_secs / 3600;
	wall_time_min = (walltime_secs % 3600) / 60;
	cpu_time_sec = (cputime_secs % 3600) % 60;
	wall_time_sec = (walltime_secs % 3600) % 60;
		
	printf("Estimated total run time: cpu time = %02d:%02d:%02d hours, wall time = %02d:%02d:%02d hours\n",
	cpu_time_hours,cpu_time_min,cpu_time_sec,wall_time_hours,wall_time_min,wall_time_sec);
	
}


/****************************************************
* Get phase speed for linear model
*****************************************************/
#if 0
	/*
	int i_max = NX-5;
	double last_time = 0;
	double vel[30];
	int vel_count = 0;
	*/

	int temp_i = i_max;
	
	i_max = find_max3(pi,i_max);
	
	if(temp_i!=i_max){
	
		double speed = dx*(double)(i_max-temp_i) / ((mtime-last_time));
		double vel_sum = 0;
	
		if(speed < 3 && speed > -9 && i_max < NX-3){

			vel[vel_count%30] = speed;			
			vel_count++;
		}
		
		for(int i=0;i<30;i++){
			vel_sum += vel[i];
		}
	
		printf("%f %d %f %f\n",mtime/3600.0,i_max,speed,vel_sum/30.0);
		last_time = mtime;
	}
#endif

			//-----------------------------------------------------------------
			// APPLY DIFFUSION
			//-----------------------------------------------------------------
			//QVP(i,j,k) -= step * dt * ( (tau_11_q_h-tau_11_q_l)*one_d_dx + (tau_22_q_h-tau_22_q_l)*one_d_dy + (tau_q_h-tau_q_l)*ONE_D_DZ(k) );
			//QCP(i,j,k) -= step * dt * ( (tau_11_c_h-tau_11_c_l)*one_d_dx + (tau_22_c_h-tau_22_c_l)*one_d_dy + (tau_c_h-tau_c_l)*ONE_D_DZ(k) );		
			//QRP(i,j,k) -= step * dt * ( (tau_11_r_h-tau_11_r_l)*one_d_dx + (tau_22_r_h-tau_22_r_l)*one_d_dy + (tau_r_h-tau_r_l)*ONE_D_DZ(k) );	
			//THP(i,j,k) -= step * dt * ( (tau_11_t_h-tau_11_t_l)*one_d_dx + (tau_22_t_h-tau_22_t_l)*one_d_dy + (tau_t_h-tau_t_l)*ONE_D_DZ(k) );


//		if(PARALLEL && j+jbs[rank]>10 && j+jbs[rank]<70 && i+ibs[rank]>30 && i+ibs[rank]<50 ){
			
			//UP(i,j,k) += 0.1 * (rand() / RAND_MAX - 0.5);
			//VP(i,j,k) += 0.1 * (rand() / RAND_MAX - 0.5);
			//WP(i,j,k) += 0.1 * (rand() / RAND_MAX - 0.5);
			
			//gustFactor = 30;
//		}


			//if(PARALLEL && j+jbs[rank]>10 && j+jbs[rank]<70 && i+ibs[rank]>30 && i+ibs[rank]<50 ){
			
				//printf("%f ",((double)rand() / (double)RAND_MAX - 0.5));
			
				//UM(i,j,k) += 1.0 * ((double)rand() / (double)RAND_MAX - 0.5);
				//THM(i,j,k) += 1.0 * ( (double)rand() / (double)RAND_MAX - 0.5 );				
				//printf("%f ",UM(i,j,k));
				
				//VM(i,j,k) += 1.0 * ((double)rand() / (double)RAND_MAX - 0.5);
				//WM(i,j,k) += 1.0 * ((double)rand() / (double)RAND_MAX - 0.5);
			
				//gustFactor = 30;
			//}


#if 0

	double water_temp_K = water_temp + 273.15;
	
	double esl = 611.2 * exp(17.67 * (water_temp_K-273.15) / (water_temp_K - 29.65) );
	

	
	qvs_surface = (380.0/pressfc) * exp(17.27*(water_temp_K-273.15)/ (water_temp_K-36.0));
	tmp_surface = water_temp_K * pow((100000.0/pressfc),0.286);
	
	double qvs = 0.62197 * esl / (pressfc-esl);

	//svp1 = 0.6112,svp2 = 17.67, spv3 = 29.65, svpt0 = 273.15

	double f = 17.67 * (273.15-29.65) * Lv / cp;

	double phi = pressfc / (pressfc-esl) * qvs * f  /  pow((tmp_surface-29.65),2);
	
	//double phi = qvs*(17.27*237*Lv)/(cp*pow(tmp_surface-36.,2.));
	
	double dqv = (0.0277-qvs) / (1.0+phi);

	double xk = Rd/cp;
	double pisfc = pow((pressfc/p0),xk);
	
	double gam = Lv / (1004.0 * pisfc);
	
	double dtheta = gam * dqv;
	
	water_temp_K += dtheta;	

	
	
	printf("qvs = %f temp = %f esl = %f qvs = %f phi = %f dqv = %f dtheta = %f\n",qvs_surface,tmp_surface,esl,qvs,phi,dqv,dtheta);
	
	esl = 611.2 * exp(17.67 * (water_temp_K-273.15) / (water_temp_K - 29.65) );
	
	qvs = 0.62197 * esl / (pressfc-esl);
	
	qvs_surface = 0.0277 - dqv;
	
	printf("esl = %f qvs = %f qvs_surface =%f\n",esl,qvs,qvs_surface);	
	fflush(stdout);
	exit(0);


					/*
					UP(i,j,k) += ufric * step * dt;
					VP(i,j,k) += vfric * step * dt;
					WP(i,j,k) += wfric * step * dt;
					*/


/*********************************************************************
*
*
**********************************************************************/
void turbulent_stress_perturbation(int il,int ih,int jl,int jh,double step){

	double wind_speed,wind_speed_base;
	double tau_u_l,tau_u_h,tau_v_l,tau_v_h;
	double tau_c_l,tau_c_h,tau_r_l,tau_r_h,tau_t_l,tau_t_h;
	double drag_coef,tau_q_h,tau_q_l;
	double ufric,vfric;

	const double qvs_surface = 0.0255;
	const double tmp_surface = 302.1;

	int k;

	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		/********************************************************	
		* Calculate surface stress
		*********************************************************/
		k = HTOPO(i,j)+1;

		wind_speed = (U(i,j,k)+UBAR(i,j,k))*(U(i,j,k)+UBAR(i,j,k)) + 
					 (V(i,j,k)+VBAR(i,j,k))*(V(i,j,k)+VBAR(i,j,k));

		wind_speed_base = UBAR(i,j,k)*UBAR(i,j,k) + VBAR(i,j,k)*VBAR(i,j,k);

		wind_speed_base = sqrt(wind_speed_base);
			
		wind_speed = sqrt(wind_speed);

		tau_q_l = 0;
		tau_c_l = 0;
		tau_r_l = 0;
		tau_t_l = 0;
		
		/********************************************************	
		* Handle land/sea contrasts
		*********************************************************/
		if(LANDSEA(i,j) < 0.5){
		
			tau_u_l = 0.003 * (wind_speed * U(i,j,k) + (wind_speed - wind_speed_base) * UBAR(i,j,k) );
			tau_v_l = 0.003 * (wind_speed * V(i,j,k) + (wind_speed - wind_speed_base) * VBAR(i,j,k) );
			
		} else {
			
			tau_u_l = 0.001 * (wind_speed * U(i,j,k) + (wind_speed - wind_speed_base) * UBAR(i,j,k) );
			tau_v_l = 0.001 * (wind_speed * V(i,j,k) + (wind_speed - wind_speed_base) * VBAR(i,j,k) );
			
			#if USE_MICROPHYSICS && USE_SURFACE_FLUX
			//if(j+jbs[rank]<300 && i+ibs[rank]>182 && i+ibs[rank]<NX-1 ){
            if(j+jbs[rank]>177 && j+jbs[rank]<300 && i+ibs[rank]>182+111 && i+ibs[rank]<NX-1){
				drag_coef = (1.1e-3 + 4.0e-5*wind_speed)*LANDSEA(i,j);

				tau_c_l = 0;
				tau_r_l = 0;
				tau_t_l = drag_coef*wind_speed*(tmp_surface - (TH(i,j,k) +THBAR(i,j,k) + tb[k]) );
				tau_q_l = drag_coef*wind_speed*(qvs_surface - (QV(i,j,k) + QBAR(i,j,k) + qb[k]) );
			}
			#endif
		}

		/********************************************************	
		* Apply turbulent diffusion to wind
		*********************************************************/
		for(int k=HTOPO(i,j)+1;k<NZ-1;k++){
			
			tau_u_h = 0.5*(KMIX(i,j,k+1)+KMIX(i-1,j,k+1)) * (U(i,j,k+1)-U(i,j,k)) * ONE_D_DZW(k+1);
			tau_v_h = 0.5*(KMIX(i,j,k+1)+KMIX(i,j-1,k+1)) * (V(i,j,k+1)-V(i,j,k)) * ONE_D_DZW(k+1);

			ufric = (tau_u_h-tau_u_l) * ONE_D_DZ(k);
			vfric = (tau_v_h-tau_v_l) * ONE_D_DZ(k);

			FRICTION(i,j,k) = -ufric*UP(i,j,k) - vfric*VP(i,j,k);

#if VORTICITY_BUDGET
			vort_ufric[INDEX(i,j,k)] += ufric * step * dt;
			vort_vfric[INDEX(i,j,k)] += vfric * step * dt;
#endif

			UP(i,j,k) += ufric * step * dt;
			VP(i,j,k) += vfric * step * dt;

			tau_u_l = tau_u_h;
			tau_v_l = tau_v_h;
		}
		
		/********************************************************	
		* Apply turbulent diffusion to mixing ratio
		*********************************************************/
		#if USE_MICROPHYSICS && USE_SURFACE_FLUX
		for(int k=HTOPO(i,j)+1;k<NZ-1;k++){

			tau_q_h = -KSMIX(i,j,k+1)*(QV(i,j,k+1) - QV(i,j,k))*ONE_D_DZW(k+1);		// vapor
			
			tau_c_h = -KSMIX(i,j,k+1)*(QC(i,j,k+1) - QC(i,j,k))*ONE_D_DZW(k+1);		// cloud

			tau_r_h = -KSMIX(i,j,k+1)*(QR(i,j,k+1) - QR(i,j,k))*ONE_D_DZW(k+1);		// rain

			tau_t_h = -KSMIX(i,j,k+1)*(TH(i,j,k+1) - TH(i,j,k))*ONE_D_DZW(k+1);		// temperature

			QVP(i,j,k) -= step * dt * (tau_q_h-tau_q_l)*ONE_D_DZ(k);

			QCP(i,j,k) -= step * dt * (tau_c_h-tau_c_l)*ONE_D_DZ(k);
			
			QRP(i,j,k) -= step * dt * (tau_r_h-tau_r_l)*ONE_D_DZ(k);
			
			THP(i,j,k) -= step * dt * (tau_t_h-tau_t_l)*ONE_D_DZ(k);
		
			tau_q_l = tau_q_h;
			tau_c_l = tau_c_h;
			tau_r_l = tau_r_h;			
			tau_t_l = tau_t_h;	
			
						
		}
		#endif
	}}
}


/*********************************************************************
*
*
**********************************************************************/
void turbulent_stress(int il,int ih,int jl,int jh,double step){

	double wind_speed,wind_speed_base;
	double tau_13_u_l,tau_13_u_h,tau_23_v_l,tau_23_v_h,tau_13_w_l,tau_13_w_h;
	double tau_11_u_l,tau_11_u_h,tau_22_v_l,tau_22_v_h,tau_33_w_l,tau_33_w_h;
	double tau_12_u_h,tau_12_u_l,tau_12_v_h,tau_12_v_l,tau_23_w_h,tau_23_w_l;
	double tau_c_l,tau_c_h,tau_r_l,tau_r_h,tau_t_l,tau_t_h;
	double tau_11_c_l,tau_22_c_l,tau_11_r_l,tau_22_r_l,tau_11_q_l,tau_22_q_l,tau_11_t_l,tau_22_t_l;
	double tau_11_c_h,tau_22_c_h,tau_11_r_h,tau_22_r_h,tau_11_q_h,tau_22_q_h,tau_11_t_h,tau_22_t_h;
	double drag_coef,tau_q_h,tau_q_l;
	double ufric,vfric,wfric;
	double KH,KL;

	double gustFactor;

	int k;

	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		/********************************************************	
		* Calculate surface stress
		*********************************************************/
		k = HTOPO(i,j)+1;

		wind_speed = sqrt((UM(i,j,k)+UBAR(i,j,k))*(UM(i,j,k)+UBAR(i,j,k)) + 
						  (VM(i,j,k)+VBAR(i,j,k))*(VM(i,j,k)+VBAR(i,j,k)));
		
		wind_speed_base = sqrt(UBAR(i,j,k)*UBAR(i,j,k) + VBAR(i,j,k)*VBAR(i,j,k));



//		if(j+jbs[rank]>47 && j+jbs[rank]<141 && i+ibs[rank]>47 && i+ibs[rank]<141 ){
		//if(j+jbs[rank]>57 && j+jbs[rank]<131 && i+ibs[rank]>57 && i+ibs[rank]<131 ){			
			//gustFactor = gustiness*(rand()/RAND_MAX);

			gustFactor = 10;//gustCorr*RED(i,j) + gustCorrSqrt * ( (   (double)rand() / (double)RAND_MAX )  );//  - 0.5*gustiness);
		
			//wind_speed += gustFactor;
			
			//wind_speed += 10;
			//} else {
			//gustFactor = 0;
			//}


		//RED(i,j) = gustFactor;

		//if(wind_speed < 0.0){ wind_speed = 0;}
			
		tau_q_l = 0;
		tau_c_l = 0;
		tau_r_l = 0;
		tau_t_l = 0;

		/********************************************************	
		* Handle land/sea contrasts
		*********************************************************/
		if(LANDSEA(i,j) < 0.5){
		
			tau_13_u_l = 0.003 * (wind_speed * UM(i,j,k) + (wind_speed - wind_speed_base) * UBAR(i,j,k) );
			tau_23_v_l = 0.003 * (wind_speed * VM(i,j,k) + (wind_speed - wind_speed_base) * VBAR(i,j,k) );
			tau_33_w_l = 0;
		
			//tau_13_u_l = 0.003 * wind_speed * U(i,j,k);
			//tau_23_v_l = 0.003 * wind_speed * V(i,j,k);
			
		} else {
			
			tau_13_u_l = 0.001 * (wind_speed * UM(i,j,k) + (wind_speed - wind_speed_base) * UBAR(i,j,k) );
			tau_23_v_l = 0.001 * (wind_speed * VM(i,j,k) + (wind_speed - wind_speed_base) * VBAR(i,j,k) );
			tau_33_w_l = 0;
			
			//tau_13_u_l = 0.001 * wind_speed * U(i,j,k);
			//tau_23_v_l = 0.001 * wind_speed * V(i,j,k);
			
			#if USE_MICROPHYSICS && USE_SURFACE_FLUX
			//if(j+jbs[rank]>0 && j+jbs[rank]<NY-1 && i+ibs[rank]>0 && i+ibs[rank]<NX-1 ){
            //if(j+jbs[rank]>177 && j+jbs[rank]<300 && i+ibs[rank]>182+111 && i+ibs[rank]<NX-1){
				drag_coef = (1.1e-3 + 4.0e-5*wind_speed)*LANDSEA(i,j)*ISTOPO(i,j,1);

				tau_q_l = drag_coef * ( (gustFactor + wind_speed - wind_speed_base) * (qvs_surface - QBAR(i,j,k) - qb[k]) - (gustFactor) * QVM(i,j,k) );
				tau_t_l = drag_coef * ( (gustFactor + wind_speed - wind_speed_base) * (tmp_surface - THBAR(i,j,k) - tb[k]) - (gustFactor) * THM(i,j,k) );

				tau_c_l = 0;
				tau_r_l = 0;
				//tau_t_l = drag_coef*wind_speed*(tmp_surface - (THBAR(i,j,k)+TH(i,j,k)+tb[k] ) );
				//tau_q_l = drag_coef*wind_speed*(qvs_surface - (QV(i,j,k) + QBAR(i,j,k) + qb[k]) );
				//}
			#endif
		}
		
		//FRICTION(i,j,NZ-1) = tau_q_l * Lv * rhou[k];
		
		if(PARALLEL && j+jbs[rank] == yp && i+ibs[rank] == xp ){
			printf("LHF = %f %f %f\n",FRICTION(i,j,NZ-1),gustFactor,gustCorrSqrt);
		}

		
#if 1
		/********************************************************	
		* Apply turbulent diffusion to wind
		*********************************************************/
		for(int k=HTOPO(i,j)+2;k<NZ-1;k++){
			
			tau_33_w_h = 0.5*(KMIX(i,j,k)+KMIX(i,j,k+1)) * (WM(i,j,k+1)-WM(i,j,k)) * ONE_D_DZ(k);
			//tau_33_w_l = 0.5*(KMIX(i,j,k)+KMIX(i,j,k-1)) * (W(i,j,k  )-W(i,j,k-1)) * ONE_D_DZ(k);
			
			tau_23_w_h = 0.25*(KHMIX(i,j,k)+KHMIX(i,j,k-1)+KHMIX(i,j+1,k)+KHMIX(i,j+1,k-1)) * (
							(WM(i,j+1,k)-WM(i,j  ,k  )) * one_d_dy +
							(VM(i,j+1,k)-VM(i,j+1,k-1)) * ONE_D_DZW(k+1)
						);
			
			tau_23_w_l = 0.25*(KHMIX(i,j,k)+KHMIX(i,j,k-1)+KHMIX(i,j-1,k)+KHMIX(i,j-1,k-1)) * (
							(WM(i,j,k)-WM(i,j-1,k  )) * one_d_dy +
							(VM(i,j,k)-VM(i,j  ,k-1)) * ONE_D_DZW(k+1)
						);
			
			tau_13_w_h = 0.25*(KHMIX(i,j,k)+KHMIX(i,j,k-1)+KHMIX(i+1,j,k)+KHMIX(i+1,j,k-1)) * (
							(WM(i+1,j,k)-WM(i  ,j,k  )) * one_d_dx +
							(UM(i+1,j,k)-UM(i+1,j,k-1)) * ONE_D_DZW(k+1)
						);
			
			tau_13_w_l = 0.25*(KHMIX(i,j,k)+KHMIX(i,j,k-1)+KHMIX(i-1,j,k)+KHMIX(i-1,j,k-1)) * (
							(WM(i,j,k)-WM(i-1,j,k  )) * one_d_dx +
							(UM(i,j,k)-UM(i  ,j,k-1)) * ONE_D_DZW(k+1)
						);
			
			wfric = (tau_13_w_h-tau_13_w_l)*one_d_dx + (tau_23_w_h-tau_23_w_l)*one_d_dy + 2.0 * (tau_33_w_h-tau_33_w_l) * ONE_D_DZ(k); 	// DZW(k)???
			
			WP(i,j,k) += wfric * step * dt;
			
			tau_33_w_l = tau_33_w_h;
			
		}
#endif

		/********************************************************	
		* Apply turbulent diffusion to wind
		*********************************************************/
		for(int k=HTOPO(i,j)+1;k<NZ-1;k++){
			
			//-----------------------------------------------------------------
			// Compression
			//-----------------------------------------------------------------
			tau_11_u_h = KHMIX(i  ,j,k) * (UM(i+1,j,k)-UM(i  ,j,k)) * one_d_dx;
			tau_11_u_l = KHMIX(i-1,j,k) * (UM(i  ,j,k)-UM(i-1,j,k)) * one_d_dx;
			
			tau_22_v_h = KHMIX(i,j  ,k) * (VM(i,j+1,k)-VM(i,j  ,k)) * one_d_dy;
			tau_22_v_l = KHMIX(i,j-1,k) * (VM(i,j  ,k)-VM(i,j-1,k)) * one_d_dy;

			//-----------------------------------------------------------------
			// Shear
			//-----------------------------------------------------------------
			tau_13_u_h = 0.5*(KMIX(i,j,k+1)+KMIX(i-1,j,k+1)) * (
							(WM(i,j,k+1)-WM(i-1,j,k+1)) * one_d_dx +
							(UM(i,j,k+1)-UM(i  ,j,k  )) * ONE_D_DZW(k+1)
						);
			
			tau_23_v_h = 0.5*(KMIX(i,j,k+1)+KMIX(i,j-1,k+1)) * (
							(WM(i,j,k+1)-WM(i,j-1,k+1)) * one_d_dy +
							(VM(i,j,k+1)-VM(i,j  ,k  )) * ONE_D_DZW(k+1)
						);

			tau_12_u_h = 0.25*(KHMIX(i,j,k) + KHMIX(i-1,j,k) + KHMIX(i-1,j+1,k) + KHMIX(i,j+1,k)) * (
							(UM(i,j+1,k) - UM(i  ,j  ,k)) * one_d_dy	+	// DU/DY at j+1 of u-cell
							(VM(i,j+1,k) - VM(i-1,j+1,k)) * one_d_dx		// DV/DX at j+1 of u-cell
						);
			
			tau_12_u_l = 0.25*(KHMIX(i,j-1,k) + KHMIX(i-1,j-1,k) + KHMIX(i-1,j,k) + KHMIX(i,j,k)) * (
							(UM(i,j,k) - UM(i,j-1,k)) * one_d_dy	+		// DU/DY at j-1 of u-cell
						 	(VM(i,j,k) - VM(i-1,j,k)) * one_d_dx			// DV/DX at j-1 of u-cell
						);

			tau_12_v_h = 0.25*(KHMIX(i,j-1,k) + KHMIX(i,j,k) + KHMIX(i+1,j-1,k) + KHMIX(i+1,j,k)) * (
							(UM(i+1,j,k) - UM(i+1,j-1,k)) * one_d_dy +	// DU/DY at i+1 of v-cell
							(VM(i+1,j,k) - VM(i  ,j  ,k)) * one_d_dx		// DV/DX at i+1 of v-cell
						);
			
			tau_12_v_l = tau_12_u_l;
/*
			tau_12_v_l = 0.25*(KHMIX(i-1,j-1,k) + KHMIX(i-1,j,k) + KHMIX(i,j-1,k) + KHMIX(i,j,k)) * (
							(U(i,j,k) - U(i,j-1,k)) * one_d_dy +	// DU/DY at i-1 of v-cell				
							(V(i,j,k) - V(i-1,j  ,k)) * one_d_dx		// DV/DX at i-1 of v-cell
						);
*/
			ufric = 2.0*(tau_11_u_h-tau_11_u_l)*one_d_dx + 	   (tau_12_u_h-tau_12_u_l)*one_d_dy + (tau_13_u_h-tau_13_u_l) * ONE_D_DZ(k);
			vfric = 	(tau_12_v_h-tau_12_v_l)*one_d_dx + 2.0*(tau_22_v_h-tau_22_v_l)*one_d_dy + (tau_23_v_h-tau_23_v_l) * ONE_D_DZ(k);

			double test = 2.0*tau_11_u_h;

			FRICTION(i,j,k) = KMIX(i,j,k);//-ufric*UP(i,j,k) - vfric*VP(i,j,k);
//			if(test != FRICTION(i,j,k)){ printf("%d %d %d %d %.20f %.20f\n",i,j,k,rank,FRICTION(i,j,k),test);exit(0);}

#if VORTICITY_BUDGET
			vort_ufric[INDEX(i,j,k)] += ufric * step * dt;
			vort_vfric[INDEX(i,j,k)] += vfric * step * dt;
#endif

			UP(i,j,k) += ufric * step * dt;
			VP(i,j,k) += vfric * step * dt;


			tau_13_u_l = tau_13_u_h;
			tau_23_v_l = tau_23_v_h;
		}
		
		/********************************************************	
		* Apply turbulent diffusion to mixing ratio
		*********************************************************/
		#if USE_MICROPHYSICS && USE_SURFACE_FLUX
		for(int k=HTOPO(i,j)+1;k<NZ-1;k++){

			//-----------------------------------------------------------------
			// DZ high
			//-----------------------------------------------------------------
			tau_q_h = -KSMIX(i,j,k+1)*(QVM(i,j,k+1) - QVM(i,j,k))*ONE_D_DZW(k+1);		// vapor
#if 1
			tau_c_h = -KSMIX(i,j,k+1)*(QCM(i,j,k+1) - QCM(i,j,k))*ONE_D_DZW(k+1);		// cloud

			tau_r_h = -KSMIX(i,j,k+1)*(QRM(i,j,k+1) - QRM(i,j,k))*ONE_D_DZW(k+1);		// rain

			tau_t_h = -KSMIX(i,j,k+1)*(THM(i,j,k+1) - THM(i,j,k))*ONE_D_DZW(k+1);		// temperature

			//-----------------------------------------------------------------
			// DY
			//-----------------------------------------------------------------
			KH = 0.5 * (KHSMIX(i,j+1,k)+KHSMIX(i,j,k));
			KL = 0.5 * (KHSMIX(i,j,k)+KHSMIX(i,j-1,k));
			
			tau_22_q_h = -KH*(QVM(i,j+1,k) - QVM(i,j,k))*one_d_dy;		// vapor
			tau_22_q_l = -KL*(QVM(i,j,k) - QVM(i,j-1,k))*one_d_dy;
			
			tau_22_c_h = -KH*(QCM(i,j+1,k) - QCM(i,j,k))*one_d_dy;		// cloud
			tau_22_c_l = -KL*(QCM(i,j,k) - QCM(i,j-1,k))*one_d_dy;

			tau_22_r_h = -KH*(QRM(i,j+1,k) - QRM(i,j,k))*one_d_dy;		// rain
			tau_22_r_l = -KL*(QRM(i,j,k) - QRM(i,j-1,k))*one_d_dy;

			tau_22_t_h = -KH*(THM(i,j+1,k) - THM(i,j,k))*one_d_dy;		// temperature
			tau_22_t_l = -KL*(THM(i,j,k) - THM(i,j-1,k))*one_d_dy;

			//-----------------------------------------------------------------
			// DX
			//-----------------------------------------------------------------
			KH = 0.5 * (KHSMIX(i+1,j,k)+KHSMIX(i,j,k));
			KL = 0.5 * (KHSMIX(i,j,k)+KHSMIX(i-1,j,k));
			
			tau_11_q_h = -KH*(QVM(i+1,j,k) - QVM(i,j,k))*one_d_dx;		// vapor
			tau_11_q_l = -KL*(QVM(i,j,k) - QVM(i-1,j,k))*one_d_dx;
			
			tau_11_c_h = -KH*(QCM(i+1,j,k) - QCM(i,j,k))*one_d_dx;		// cloud
			tau_11_c_l = -KL*(QCM(i,j,k) - QCM(i-1,j,k))*one_d_dx;

			tau_11_r_h = -KH*(QRM(i+1,j,k) - QRM(i,j,k))*one_d_dx;		// rain
			tau_11_r_l = -KL*(QRM(i,j,k) - QRM(i-1,j,k))*one_d_dx;

			tau_11_t_h = -KH*(THM(i+1,j,k) - THM(i,j,k))*one_d_dx;		// temperature
			tau_11_t_l = -KL*(THM(i,j,k) - THM(i-1,j,k))*one_d_dx;
			
#endif
			//-----------------------------------------------------------------
			// APPLY DIFFUSION
			//-----------------------------------------------------------------
			QVP(i,j,k) -= step * dt * ( (tau_11_q_h-tau_11_q_l)*one_d_dx + (tau_22_q_h-tau_22_q_l)*one_d_dy + (tau_q_h-tau_q_l)*ONE_D_DZ(k) );
#if 1
			QCP(i,j,k) -= step * dt * ( (tau_11_c_h-tau_11_c_l)*one_d_dx + (tau_22_c_h-tau_22_c_l)*one_d_dy + (tau_c_h-tau_c_l)*ONE_D_DZ(k) );
			
			QRP(i,j,k) -= step * dt * ( (tau_11_r_h-tau_11_r_l)*one_d_dx + (tau_22_r_h-tau_22_r_l)*one_d_dy + (tau_r_h-tau_r_l)*ONE_D_DZ(k) );
			
			THP(i,j,k) -= step * dt * ( (tau_11_t_h-tau_11_t_l)*one_d_dx + (tau_22_t_h-tau_22_t_l)*one_d_dy + (tau_t_h-tau_t_l)*ONE_D_DZ(k) );
#endif		
			tau_q_l = tau_q_h;
			tau_c_l = tau_c_h;
			tau_r_l = tau_r_h;			
			tau_t_l = tau_t_h;			
			
						
		}
		#endif
	}}
}
#endif


			/********************************************************	
			* Shear
			*********************************************************/
#if 0			
			dudz = 0.5 * ( U(i,j,k) + U(i+1,j,k) - U(i,j,k-1) - U(i+1,j,k-1) ) * ONE_D_DZW(k);
			dvdz = 0.5 * ( V(i,j,k) + V(i,j+1,k) - V(i,j,k-1) - V(i,j+1,k-1) ) * ONE_D_DZW(k);
		
			dUdz = 0.5 * ( UBAR(i,j,k) + UBAR(i+1,j,k) - UBAR(i,j,k-1) - UBAR(i+1,j,k-1) ) * ONE_D_DZW(k);
			dVdz = 0.5 * ( VBAR(i,j,k) + VBAR(i,j+1,k) - VBAR(i,j,k-1) - VBAR(i,j+1,k-1) ) * ONE_D_DZW(k);

			dUdz += dudz;
			dVdz += dvdz;
	
			S2 = (dUdz*dUdz+dVdz*dVdz);
#endif


	/***********************************************************************************
	*	
	* 	INITIALIZE EDGES FOR RATE OF DEFORMATION TENSOR					
	*
	************************************************************************************/
/*
	i = 3;
	
	for(int j=jl-1;j<jh;j++){
		
		if(j==jl-1){ k_begin = 1;} else { k_begin = HTOPO(i,j)+1;}

		for(int k=k_begin;k<NZ-1;k++){
		
			STRESS(j,k).tau_11_west = 2.0 * (UM(i,j,k) - UM(i-1,j,k)) * one_d_dx;
			STRESS(j,k).tau_12_west = (UM(i,j+1,k) - UM(i,j,k)) * one_d_dy + (VM(i,j+1,k) - VM(i-1,j+1,k)) * one_d_dx;
			STRESS(j,k).tau_13_west = (WM(i,j,k+1) - WM(i-1,j,k+1)) * one_d_dx + (UM(i,j,k+1)-UM(i,j,k)) * ONE_D_DZW(k+1);
		}
	}
*/

#if 0
	/*********************************************
	* SURFACE FLUX
	**********************************************/
	if(USE_TURBULENT_STRESS){
		
		/*********************************************
		* With moisture
		**********************************************/
		if(USE_MICROPHYSICS && bigcounter>stopheating){
		
			//exchangeOnePoint(ums); exchangeOnePoint(vms); exchangeOnePoint(wms);
			//exchangeOnePoint(qvms); exchangeOnePoint(qcms); exchangeOnePoint(qrms);
			//exchangeOnePoint(thms);
		
			//get_Kmix(3,fNX-3,3,fNY-3);
		
			//turbulent_stress(3,fNX-3,3,fNY-3,1.0);
			turbulent_diffusion_scalars(3,fNX-3,3,fNY-3,1.0);
		
			test();
			
		/*********************************************
		* Without moisture
		**********************************************/
		} else {
			
			get_Kmix_vertical(3,fNX-3,3,fNY-3);
			turbulent_stress_vertical(3,fNX-3,3,fNY-3,1.0);
		}
	}
#endif

#if 0
double max_u=0,max_v=0,max_w=0,max_th=0,max_qv=0,max_qc=0,max_qr=0,max_k=0;
int max_u_i = -1; int max_u_j = -1; int max_u_k = -1;
int max_v_i = -1; int max_v_j = -1; int max_v_k = -1;
int max_w_i = -1; int max_w_j = -1; int max_w_k = -1;
int max_th_i = -1; int max_th_j = -1; int max_th_k = -1;
int max_qv_i = -1; int max_qv_j = -1; int max_qv_k = -1;
int max_qc_i = -1; int max_qc_j = -1; int max_qc_k = -1;
int max_qr_i = -1; int max_qr_j = -1; int max_qr_k = -1;
int max_k_i = -1; int max_k_j = -1; int max_k_k = -1;

max_u=0,max_v=0,max_w=0,max_th=0,max_qv=0,max_qc=0,max_qr=0,max_k=0;
max_u_i = -1; max_u_j = -1; max_u_k = -1;
max_v_i = -1; max_v_j = -1; max_v_k = -1;
max_w_i = -1; max_w_j = -1; max_w_k = -1;
max_th_i = -1; max_th_j = -1; max_th_k = -1;
max_qv_i = -1; max_qv_j = -1; max_qv_k = -1;
max_qc_i = -1; max_qc_j = -1; max_qc_k = -1;
max_qr_i = -1; max_qr_j = -1; max_qr_k = -1;
max_k_i = -1; max_k_j = -1; max_k_k = -1;		

for(int i=0;i<fNX;i++){
for(int j=0;j<fNY;j++){
for(int k=0;k<fNZ;k++){
	
	if(fabs(U(i,j,k))  >= max_u){  max_u  = fabs(U(i,j,k));   max_u_i = i;  max_u_j = j;  max_u_k = k; }
	if(fabs(V(i,j,k))  >= max_v){  max_v  = fabs(V(i,j,k));   max_v_i = i;  max_v_j = j;  max_v_k = k; }
	if(fabs(W(i,j,k))  >= max_w){  max_w  = fabs(W(i,j,k));   max_w_i = i;  max_w_j = j;  max_w_k = k; }
	if(fabs(TH(i,j,k)) >= max_th){ max_th = fabs(TH(i,j,k));  max_th_i = i; max_th_j = j; max_th_k = k; }
	if(fabs(QV(i,j,k)) >= max_qv){ max_qv = fabs(QV(i,j,k));  max_qv_i = i; max_qv_j = j; max_qv_k = k; }
	if(fabs(QC(i,j,k)) >= max_qc){ max_qc = fabs(QC(i,j,k));  max_qc_i = i; max_qc_j = j; max_qc_k = k; }		
	if(fabs(QR(i,j,k)) >= max_qr){ max_qr = fabs(QR(i,j,k));  max_qr_i = i; max_qr_j = j; max_qr_k = k; }
	if(fabs(FRICTION(i,j,k)) >= max_k){ max_k = FRICTION(i,j,k);  max_k_i = i; max_k_j = j; max_k_k = k; }
	
}}}

printf("\nAt substep %d of rank %d\n",s,rank);
printf("At %d %d %d of rank %d u_max = %f || ",max_u_i, max_u_j, max_u_k,rank,max_u);
printf("At %d %d %d of rank %d v_max = %f || ",max_v_i, max_v_j, max_v_k,rank,max_v);
printf("At %d %d %d of rank %d w_max = %f\n",max_w_i, max_w_j, max_w_k,rank,max_w);
printf("At %d %d %d of rank %d qv_max = %f || ",max_qv_i, max_qv_j, max_qv_k,rank,max_qv*1000.0);		
printf("At %d %d %d of rank %d qc_max = %f || ",max_qc_i, max_qc_j, max_qc_k,rank,max_qc*1000.0);		
printf("At %d %d %d of rank %d qr_max = %f\n",max_qr_i, max_qr_j, max_qr_k,rank,max_qr*1000.0);
printf("At %d %d %d of rank %d th_max = %f || ",max_th_i, max_th_j, max_th_k,rank,max_th);
printf("At %d %d %d of rank %d k_max = %f\n",max_k_i, max_k_j, max_k_k,rank,max_k);		

fflush(stdout);
#endif

//	for(int i=0;i<NX;i++){
//	for(int j=0;j<NY;j++){
//	for(int k=0;k<NZ;k++){
//		pim[i][j][k] = pi[i][j][k]*cp*rhou[k]*tbv[k]*(1.+0.61*qb[k])*0.01;
//		//pim[i][j][k] = pip[i][j][k]*cp*rhou[k]*tbv[k]*(1.+0.61*qb[k]);
//	}}}

//	printf("%f\n",pim[30][30][1]);
	



//	ub_at_v = 0.25 * (ubar[i][j][k]+ubar[i+1][j][k]+ubar[i][j-1][k]+ubar[i+1][j-1][k]);
//	vb_at_u = 0.25 * (vbar[i][j][k]+vbar[i][j+1][k]+vbar[i-1][j][k]+vbar[i-1][j+1][k]);

//		- dt * rhou[k] * friction[i][j][k] * u[i][j][k] * sqrt( (u[i][j][k]+ubar[i][j][k])*(u[i][j][k]+ubar[i][j][k]) +
//																		(v_at_u+vb_at_u)*(v_at_u+vb_at_u)+0.1
//																	) / dz

//		- dt * rhou[k] * friction[i][j][k] * v[i][j][k] * sqrt( (v[i][j][k]+vbar[i][j][k])*(v[i][j][k]+vbar[i][j][k]) +
//																		(u_at_v+ub_at_v)*(u_at_v+ub_at_v)+0.1
//																	) / dz

//	int i = 92; int j = 29;
//	double div;
//	double avg = 0;
//	double rho2 = 0;

//	printf("topo = %d %d %d %d %d\n",htopo[i][j-1],htopo[i][j],htopo[i][j+1],htopo[i-1][j],htopo[i+1][j]);

//	for(int k=htopo[i][j]+1;k<NZ;k++){

//		div = (up[i+1][j][k] - up[i][j][k])/dx  + (vp[i][j+1][k] - vp[i][j][k])/dy;

//		avg = avg + rhou[k]*div;
//		rho2 = rho2 + rhou[k];

//		printf("div = %d %e %f %f %f %f %f\n",k,div,w[i][j][k],up[i+1][j][k],up[i][j][k],vp[i][j+1][k],vp[i][j][k]);

//	}

//	avg = avg / (double)(NZ-2-htopo[i][j]);

//	printf("avg = %e\n",avg);


//	int i = 90; int j = 52;
//	double div;
//	double avg = 0;
//	double rho2 = 0;

//	printf("topo = %d %d %d %d %d\n",htopo[i][j-1],htopo[i][j],htopo[i][j+1],htopo[i-1][j],htopo[i+1][j]);

//	for(int k=htopo[i][j]+1;k<NZ-1;k++){

//		div = (up[i+1][j][k] - up[i][j][k])/dx  + (vp[i][j+1][k] - vp[i][j][k])/dy;

//		avg = avg + rhou[k]*div;
//		rho2 = rho2 + rhou[k];

//		printf("div = %d %e %f %f %f %f %f\n",k,div,w[i][j][k],up[i+1][j][k],up[i][j][k],vp[i][j+1][k],vp[i][j][k]);

//	}

//	avg = avg / rho2;

//	printf("avg = %e\n",avg);

//	test2();

//	double ftest[NX][NY];
//	double store[NX][NY];

//	double radius;
//	double max_radius = 20;

//	for(int i=0;i<NX;i++){
//	for(int j=0;j<NY;j++){

//		radius = (double)sqrt(pow(NX/2-i,2)+pow(NY/2-j,2));

//		if(radius < max_radius){

//			ftest[i][j] = (radius-max_radius)*(radius-max_radius);
//			store[i][j] = ftest[i][j];	

//		} else {

//			ftest[i][j] = 0;
//			store[i][j] = 0;
//		}

//	}}

//	init_fft2d(NX,NY,(double)NX*dx,(double)NY*dy);
//	poisson_fft2d(ftest);

//	double test;
//	int j = NY/2;

//	for(int i=1;i<NX-1;i++){
//	//for(int j=1;j<NY-1;j++){

//		test = (ftest[i+1][j] + ftest[i][j+1] + ftest[i][j-1] + ftest[i-1][j] - 4*ftest[i][j])/(dx*dx);

//		printf("%d %d %f %f %f\n",i,j,ftest[i][j],test,store[i][j]);

//	}//}
