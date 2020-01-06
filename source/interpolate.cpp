
#include "stdafx.h"
#include "interpolate.h"

#define INTERPOLATE_VERBOSE true

/****************************************************************************
*
* 					INTERPOLATION PROCEDURES FOR INITIALIZATION
*
*****************************************************************************/

void vertical_interpolation_warnings(double *z,double *zi,int zdim,int nz);
 
double meters_per_lat = 111000.0;

/********************************************************
* Vertical interpolation
* 
* @param xdim,ydim,zdim - dimensions of input data
* @param zdim_out - output z-dimension
* @param z - array of heights of output data
* @param source_z - array of heights of input data
* @param out - newly interpolated output data
*
*********************************************************/
double * vert_interpolate(double *z, double *source_z, double *source_val, int xdim, int ydim, int zdim,int zdim_out){

	int count = 0;
	double mu;

	double *out  = (double *)malloc(xdim*ydim*zdim_out*sizeof(double));

	/********************************************************
	* Loop through input data
	*********************************************************/
	for(int i=0;i<xdim;i++){
	for(int j=0;j<ydim;j++){

		count = 0;

		/********************************************************
		* For each point in the output heights
		*********************************************************/
		for(int k=0;k<zdim_out;k++){

			/********************************************************
			* Our point is between these two points in our input data
			*********************************************************/
			if(  z[k] > source_z[d(i,j,count)] && z[k] <= source_z[d(i,j,count+1)] ){
		
				//---------------------------------------------------------
				// fractional distance between interpolating points
				//---------------------------------------------------------
				mu = (z[k] - source_z[d(i,j,count)]) / (source_z[d(i,j,count+1)] - source_z[d(i,j,count)]);

				/********************************************************
				* Interpolate points where all data needed is available
				*********************************************************/
				if(count!=0){

					out[d(i,j,k)] = CubicInterpolate(	
									source_val[d(i,j,count-1)],
									source_val[d(i,j,count)],
									source_val[d(i,j,count+1)],
									source_val[d(i,j,count+2)],
									source_z[d(i,j,count-1)],
									source_z[d(i,j,count)],
									source_z[d(i,j,count+1)],
									source_z[d(i,j,count+2)],
									mu
								);
				/********************************************************
				* Interpolate points where lowest point required for
				* Cubic interpolation is below input data
				*********************************************************/
				} else {

					out[d(i,j,k)] = CubicInterpolate(
									2*source_val[d(i,j,count)]-source_val[d(i,j,count+1)],
									source_val[d(i,j,count)],
									source_val[d(i,j,count+1)],
									source_val[d(i,j,count+2)],
									2*source_z[d(i,j,count)]-source_z[d(i,j,count+1)],
									source_z[d(i,j,count)],
									source_z[d(i,j,count+1)],
									source_z[d(i,j,count+2)],
									mu
								);
				}
		
			/********************************************************
			* Point is below input data
			*********************************************************/
			} else if(  z[k] < source_z[d(i,j,0)] ){ 

				mu = (z[k] - source_z[d(i,j,0)]) / (source_z[d(i,j,1)] - source_z[d(i,j,0)]);

				out[d(i,j,k)] = LinearInterpolate(source_val[d(i,j,0)],source_val[d(i,j,1)],mu);

			/********************************************************
			* Point is above input data, exit loop
			*********************************************************/
			} else if(  z[k] > source_z[d(i,j,zdim-1)] ){ 

				k = zdim_out;				

			/********************************************************
			* Point exactly coincides with input data height level
			*********************************************************/
			} else if(  z[k] == source_z[d(i,j,count)] ){ 

				out[d(i,j,k)] = source_val[d(i,j,count)];

			/********************************************************
			* No condition satisfied--compare next level in input data
			* Decrement index for output array since this level was 
			* not handeled.
			*********************************************************/
			} else { count++; k--;}
		}
	}}

	return out;
}

/********************************************************
* Vertical interpolation
* 
* @param xdim,ydim,zdim - dimensions of input data
* @param zdim_out - output z-dimension
* @param z - 1D array of heights of output data
* @param source_var - 1D array of input data to be interpolated
* @param source_z - 1D array of heights of input data
* @param out - 3D newly interpolated output data
*
*********************************************************/
double * vert_interpolate_1d(double *z, double *source_z, double *source_val, int xdim, int ydim, int zdim,int zdim_out){

	int count = 0;
	double mu;

	double *out  = (double *)malloc(xdim*ydim*zdim_out*sizeof(double));

	/********************************************************
	* Loop through input data
	*********************************************************/
	for(int i=0;i<xdim;i++){
	for(int j=0;j<ydim;j++){

		count = 0;

		/********************************************************
		* For each point in the output heights
		*********************************************************/
		for(int k=0;k<zdim_out;k++){

			/********************************************************
			* Our point is between these two points in our input data
			*********************************************************/
			if(  z[k] > source_z[d(i,j,count)] && z[k] <= source_z[d(i,j,count+1)] ){
		
				//---------------------------------------------------------
				// fractional distance between interpolating points
				//---------------------------------------------------------
				mu = (z[k] - source_z[d(i,j,count)]) / (source_z[d(i,j,count+1)] - source_z[d(i,j,count)]);

				/********************************************************
				* Interpolate points where all data needed is available
				*********************************************************/
				if(count!=0){

					out[d(i,j,k)] = CubicInterpolate(	
									source_val[count-1],
									source_val[count  ],
									source_val[count+1],
									source_val[count+2],
									source_z[d(i,j,count-1)],
									source_z[d(i,j,count)],
									source_z[d(i,j,count+1)],
									source_z[d(i,j,count+2)],
									mu
								);
				/********************************************************
				* Interpolate points where lowest point required for
				* Cubic interpolation is below input data
				*********************************************************/
				} else {

					out[d(i,j,k)] = CubicInterpolate(
									2*source_val[count]-source_val[count+1],
									source_val[count],
									source_val[count+1],
									source_val[count+2],
									2*source_z[d(i,j,count)]-source_z[d(i,j,count+1)],
									source_z[d(i,j,count)],
									source_z[d(i,j,count+1)],
									source_z[d(i,j,count+2)],
									mu
								);
				}
		
			/********************************************************
			* Point is below input data
			*********************************************************/
			} else if(  z[k] < source_z[d(i,j,0)] ){ 

				mu = (z[k] - source_z[d(i,j,0)]) / (source_z[d(i,j,1)] - source_z[d(i,j,0)]);

				out[d(i,j,k)] = LinearInterpolate(source_val[0],source_val[1],mu);

			/********************************************************
			* Point is above input data, exit loop
			*********************************************************/
			} else if(  z[k] > source_z[d(i,j,zdim-1)] ){ 

				k = zdim_out;				

			/********************************************************
			* Point exactly coincides with input data height level
			*********************************************************/
			} else if(  z[k] == source_z[d(i,j,count)] ){ 

				out[d(i,j,k)] = source_val[count];

			/********************************************************
			* No condition satisfied--compare next level in input data
			* Decrement index for output array since this level was 
			* not handeled.
			*********************************************************/
			} else { count++; k--;}
		}
	}}

	return out;
}

/********************************************************
* Vertical interpolation
* 
* @param xdim,ydim,zdim - dimensions of input data
* @param zdim_out - output z-dimension
* @param z - 1D array of heights of output data
* @param source_z - 1D array of heights of input data
* @param out - newly interpolated output data
*
*********************************************************/
void vert_interpolate_1d_from_model(double *z, double *source_z, double *source_val, int xdim, int ydim, int zdim,int zdim_out,double *out){

	int count = 0;
	double mu;

	vertical_interpolation_warnings(z,source_z,zdim,zdim_out);

	/********************************************************
	* Loop through input data
	*********************************************************/
	for(int i=0;i<xdim;i++){
	for(int j=0;j<ydim;j++){

		count = 0;

		/********************************************************
		* For each point in the output heights
		*********************************************************/
		for(int k=0;k<zdim_out;k++){

			/********************************************************
			* Our point is between these two points in our input data
			*********************************************************/
			if(  z[k] > source_z[count] && z[k] <= source_z[count+1] ){
		
				//---------------------------------------------------------
				// fractional distance between interpolating points
				//---------------------------------------------------------
				mu = (z[k] - source_z[count]) / (source_z[count+1] - source_z[count]);

				/********************************************************
				* Interpolate points where all data needed is available
				*********************************************************/
				if(count!=0 && count < zdim-2){

					out[d5(i,j,k)] = CubicInterpolate(	
									source_val[d3(i,j,count-1)],
									source_val[d3(i,j,count)],
									source_val[d3(i,j,count+1)],
									source_val[d3(i,j,count+2)],
									source_z[count-1],
									source_z[count],
									source_z[count+1],
									source_z[count+2],
									mu
								);
				/********************************************************
				* Interpolate points where lowest point required for
				* Cubic interpolation is below input data
				*********************************************************/
				} else if(count==0) {

					out[d5(i,j,k)] = CubicInterpolate(
									2*source_val[d3(i,j,count)]-source_val[d3(i,j,count+1)],
									source_val[d3(i,j,count)],
									source_val[d3(i,j,count+1)],
									source_val[d3(i,j,count+2)],
									2*source_z[count]-source_z[count+1],
									source_z[count],
									source_z[count+1],
									source_z[count+2],
									mu
								);
				/********************************************************
				* Interpolate points where lowest point required for
				* Cubic interpolation is above input data
				*********************************************************/					
				} else {
					//if(i==NX/2 && j==NY/2){
					//	printf("%f %f %f %f\n",source_val[count-1],source_val[count],source_val[count+1],source_val[count+2]);
					//	printf("%f %f %f %f\n",source_z[count-1],source_z[count],source_z[count+1],source_z[count+2]);
					//	printf("%f %f\n",2*source_val[d3(i,j,count+1)]-source_val[d3(i,j,count)],2*source_z[count+1]-source_z[count]);
					//}
					out[d5(i,j,k)] = CubicInterpolate(
									source_val[d3(i,j,count-1)],
									source_val[d3(i,j,count)],
									source_val[d3(i,j,count+1)],
									2*source_val[d3(i,j,count+1)]-source_val[d3(i,j,count)],
									source_z[count-1],
									source_z[count],
									source_z[count+1],
									2*source_z[count+1]-source_z[count],
									mu
								);
				}
				
			/********************************************************
			* Point is below input data
			*********************************************************/
			} else if(  z[k] < source_z[0] ){ 

				mu = (z[k] - source_z[0]) / (source_z[1] - source_z[0]);

				out[d5(i,j,k)] = LinearInterpolate(source_val[d3(i,j,0)],source_val[d3(i,j,1)],mu);

			/********************************************************
			* Point is above input data, exit loop
			*********************************************************/
			} else if(  z[k] > source_z[zdim-1] ){ 

				k = zdim_out;			

			/********************************************************
			* Point exactly coincides with input data height level
			*********************************************************/
			} else if(  z[k] == source_z[count] ){ 

				out[d5(i,j,k)] = source_val[d3(i,j,count)];

			/********************************************************
			* No condition satisfied--compare next level in input data
			* Decrement index for output array since this level was 
			* not handeled.
			*********************************************************/
			} else { count++; k--;}
		}
	}}
}

/********************************************************
* Vertical interpolation to model grid with terrain
* following coordinates
* 
* @param z   - 3D array of heights of output data
* @param var - 1D array of input uninterpolated data
*			   at model heights where terrain height
*			   equals zero
* @param out - output data interpolated to terrain
*			   following coordinates
*
*********************************************************/
void vert_interpolate_1d(double z[NX][NY][NZ], double var[NZ],double out[NX][NY][NZ]){

	double mu;
	int count = 0;

	/********************************************************
	* Loop through input data
	*********************************************************/
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		count = 0;

		/********************************************************
		* For each point in the output heights
		*********************************************************/
		for(int k=0;k<NZ;k++){

			/********************************************************
			* Our point is between these two points in our input data
			*********************************************************/
			if(  z[i][j][k] > zu[count] && z[i][j][k] <= zu[count+1] ){
		
				/*---------------------------------------------------------
				* Fractional distance between interpolating points
				*
				*	zu[count+1]		-------------+-------------
				*
				*
				*		z[i][j][k]		---------+---------
				*
				*	zu[count]		-------------+-------------
				*
				*----------------------------------------------------------*/
				mu = (z[i][j][k] - zu[count]) / (zu[count+1] - zu[count]);

				/********************************************************
				* Interpolate point
				*********************************************************/
				if(count!=0){

					out[i][j][k] = CubicInterpolate2(var[count-1],var[count],var[count+1],var[count+2],mu);

				} else {

					out[i][j][k] = CubicInterpolate2(2*var[count]-var[count+1],var[count],var[count+1],var[count+2],mu);
				}
		
			/********************************************************
			* Point is below input data
			*********************************************************/
			} else if(  z[i][j][k] < zu[0] ){ 

				mu = (z[i][j][k] - zu[0]) / (zu[1] - zu[0]);

				out[i][j][k] = LinearInterpolate(var[0],var[1],mu);

			/********************************************************
			* Point is above input data, exit loop
			*********************************************************/
			} else if(  z[i][j][k] > zu[NZ-1] ){ 

				k = NZ;				

			/********************************************************
			* Point exactly coincides with input data height level
			*********************************************************/
			} else if(  z[i][j][k] == zu[count] ){ 

				out[i][j][k] = var[count];

			/********************************************************
			* No condition satisfied--compare next level in input data
			* Decrement index for output array since this level was 
			* not handeled.
			*********************************************************/
			} else { count++; k--;}

		}
	}}

}

/********************************************************
* Vertical interpolation
* 
* @param z   - 3D array of heights of output data
* @param var - 1D array of input uninterpolated data
*			   at model heights where terrain height
*			   equals zero
* @param out - output data interpolated to terrain
*			   following coordinates
*
*********************************************************/
void vert_interpolate_1d_linear(int nzIn,int nzOut, double * zIn,double * zOut,double * varIn,double * varOut){

	double mu;
	int count = 0;

	/********************************************************
	* For each point in the output heights
	*********************************************************/
	for(int k=0;k<nzOut;k++){
		//printf("%d is %f between %f and %f?\n",k,zOut[k],zIn[count],zIn[count+1]);
		/********************************************************
		* Our point is between these two points in our input data
		*********************************************************/
		if(  zOut[k] > zIn[count] && zOut[k] <= zIn[count+1] ){
			//printf("matched %f %f %f\n",zOut[k],zIn[count],zIn[count+1]);
			/*---------------------------------------------------------
			* Fractional distance between interpolating points
			*
			*	zOut[count+1]	-------------+-------------
			*
			*
			*		zIn[k]			---------+---------
			*
			*	zOut[count]		-------------+-------------
			*
			*----------------------------------------------------------*/
			mu = (zOut[k] - zIn[count]) / (zIn[count+1] - zIn[count]);

			/********************************************************
			* Interpolate point
			*********************************************************/
			if(count!=0 && count < nzIn-2){
				//printf("interpolating %f %f %f %f %f\n",varIn[count-1],varIn[count],varIn[count+1],varIn[count+2],mu);
				varOut[k] = LinearInterpolate(varIn[count],varIn[count+1],mu);

			} else if(count==0) {
				varOut[k] = LinearInterpolate(varIn[count],varIn[count+1],mu);
				//varOut[k] = CubicInterpolate2(2*varIn[count]-varIn[count+1],varIn[count],varIn[count+1],varIn[count+2],mu);
			} else {
				varOut[k] = LinearInterpolate(varIn[count],varIn[count+1],mu);
				//varOut[k] = CubicInterpolate2(varIn[count-1],varIn[count],varIn[count+1],2*varIn[count+1]-varIn[count],mu);
			}
	
		/********************************************************
		* Point is below input data
		*********************************************************/
		} else if(  zOut[k] < zIn[0] ){ 

			mu = (zOut[k] - zIn[0]) / (zIn[1] - zIn[0]);

			varOut[k] = LinearInterpolate(varIn[0],varIn[1],mu);

		/********************************************************
		* Point is above input data, exit loop
		*********************************************************/
		} else if(  zOut[k] > zIn[nzIn-1] ){ 

			k = nzOut;				
			
			//printf("Vertical interpolation issue at line %d in file %s! Ending point outside of data range.\n",__LINE__,__FILE__);
			//exit(0);
			
		/********************************************************
		* Point exactly coincides with input data height level
		*********************************************************/
		} else if(  zOut[k] == zIn[count] ){ 

			varOut[k] = varIn[count];

		/********************************************************
		* No condition satisfied--compare next level in input data
		* Decrement index for output array since this level was 
		* not handeled.
		*********************************************************/
		} else { count++; k--;}

	}
		
}

/********************************************************
* Vertical interpolation
* 
* @param z   - 3D array of heights of output data
* @param var - 1D array of input uninterpolated data
*			   at model heights where terrain height
*			   equals zero
* @param out - output data interpolated to terrain
*			   following coordinates
*
*********************************************************/
void vert_interpolate_1d(int nzIn,int nzOut, double * zIn,double * zOut,double * varIn,double * varOut){

	double mu;
	int count = 0;

	vertical_interpolation_warnings(zOut,zIn,nzIn,nzOut);

	/********************************************************
	* For each point in the output heights
	*********************************************************/
	for(int k=0;k<nzOut;k++){
		//printf("%d is %f between %f and %f?\n",k,zOut[k],zIn[count],zIn[count+1]);
		/********************************************************
		* Our point is between these two points in our input data
		*********************************************************/
		if(  zOut[k] > zIn[count] && zOut[k] <= zIn[count+1] ){
			//printf("matched %f %f %f\n",zOut[k],zIn[count],zIn[count+1]);
			/*---------------------------------------------------------
			* Fractional distance between interpolating points
			*
			*	zOut[count+1]	-------------+-------------
			*
			*
			*		zIn[k]			---------+---------
			*
			*	zOut[count]		-------------+-------------
			*
			*----------------------------------------------------------*/
			mu = (zOut[k] - zIn[count]) / (zIn[count+1] - zIn[count]);

			/********************************************************
			* Interpolate point
			*********************************************************/
			if(count!=0 && count < nzIn-2){
				//printf("interpolating %f %f %f %f %f\n",varIn[count-1],varIn[count],varIn[count+1],varIn[count+2],mu);
				varOut[k] = CubicInterpolate2(varIn[count-1],varIn[count],varIn[count+1],varIn[count+2],mu);

			} else if(count==0) {
				
				varOut[k] = CubicInterpolate2(2*varIn[count]-varIn[count+1],varIn[count],varIn[count+1],varIn[count+2],mu);
			} else {
				
				varOut[k] = CubicInterpolate2(varIn[count-1],varIn[count],varIn[count+1],2*varIn[count+1]-varIn[count],mu);
			}
	
		/********************************************************
		* Point is below input data
		*********************************************************/
		} else if(  zOut[k] < zIn[0] ){ 

			mu = (zOut[k] - zIn[0]) / (zIn[1] - zIn[0]);

			varOut[k] = LinearInterpolate(varIn[0],varIn[1],mu);

		/********************************************************
		* Point is above input data, exit loop
		*********************************************************/
		} else if(  zOut[k] > zIn[nzIn-1] ){ 

			k = nzOut;				
			
			//printf("Vertical interpolation issue at line %d in file %s! Ending point outside of data range.\n",__LINE__,__FILE__);
			//exit(0);
			
		/********************************************************
		* Point exactly coincides with input data height level
		*********************************************************/
		} else if(  zOut[k] == zIn[count] ){ 

			varOut[k] = varIn[count];

		/********************************************************
		* No condition satisfied--compare next level in input data
		* Decrement index for output array since this level was 
		* not handeled.
		*********************************************************/
		} else { count++; k--;}

	}
		
}


/********************************************************
* Vertical interpolation to model grid with terrain
* following coordinates
* 
* @param z   - 3D array of heights of output data
* @param var - 3D array of input uninterpolated data
*			   at model heights where terrain height
*			   equals zero
* @param out - output data interpolated to terrain
*			   following coordinates
*
*********************************************************/
void vert_interpolate_3d(double z[NX][NY][NZ], double var[NX][NY][NZ],double out[NX][NY][NZ]){

	double mu;
	int count = 0;

	/********************************************************
	* Loop through input data
	*********************************************************/
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		count = 0;

		/********************************************************
		* For each point in the output heights
		*********************************************************/
		for(int k=0;k<NZ;k++){

			/********************************************************
			* Our point is between these two points in our input data
			*********************************************************/
			if(  z[i][j][k] > zu[count] && z[i][j][k] <= zu[count+1] ){
		
				/*---------------------------------------------------------
				* Fractional distance between interpolating points
				*
				*	zu[count+1]		-------------+-------------
				*
				*
				*		z[i][j][k]		---------+---------
				*
				*	zu[count]		-------------+-------------
				*
				*----------------------------------------------------------*/
				mu = (z[i][j][k] - zu[count]) / (zu[count+1] - zu[count]);

				/********************************************************
				* Interpolate point
				*********************************************************/
				if(count!=0){

					out[i][j][k] = CubicInterpolate2(var[i][j][count-1],var[i][j][count],var[i][j][count+1],var[i][j][count+2],mu);

				} else {

					out[i][j][k] = CubicInterpolate2(2*var[i][j][count]-var[i][j][count+1],var[i][j][count],var[i][j][count+1],var[i][j][count+2],mu);
				}
		
			/********************************************************
			* Point is below input data
			*********************************************************/
			} else if(  z[i][j][k] < zu[0] ){ 

				mu = (z[i][j][k] - zu[0]) / (zu[1] - zu[0]);

				out[i][j][k] = LinearInterpolate(var[i][j][0],var[i][j][1],mu);

			/********************************************************
			* Point is above input data, exit loop
			*********************************************************/
			} else if(  z[i][j][k] > zu[NZ-1] ){ 

				k = NZ;				

			/********************************************************
			* Point exactly coincides with input data height level
			*********************************************************/
			} else if(  z[i][j][k] == zu[count] ){ 

				out[i][j][k] = var[i][j][count];

			/********************************************************
			* No condition satisfied--compare next level in input data
			* Decrement index for output array since this level was 
			* not handeled.
			*********************************************************/
			} else { count++; k--;}

		}
	}}

}

/********************************************************
* Interpolate the horizontal dimensions of an array
*
* *uz - input array with dimension xdim x ydim x NZ
* out[NX][NY][NZ] - output array
* xstagger - does it have a staggered x-coordinate?
* ystagger - does it have a stagger y-coordinate?
*
*********************************************************/
void horz_interpolate(double *uz,double *out,int xdim,int ydim,int zdim,bool xstagger,bool ystagger,double dlat,double dlon,double lonoffset,double latoffset){

	int xcount = 0, ycount = 0;

	double p[4][4], mu, mux, muy;

	/********************************************************
	* Distances in meters for output interpolated array
	*********************************************************/
	double * x = (double *)malloc(NX*sizeof(double));
	double * y = (double *)malloc(NY*sizeof(double));

	x[0] = lonoffset*meters_per_lat;
	y[0] = latoffset*meters_per_lat;

	if(xstagger){ x[0] = x[0] - 0.5*dx;}

	if(ystagger){ y[0] = y[0] - 0.5*dy;}

	for(int i=1;i<NY;i++)
		y[i] = y[i-1] + dy;
	for(int i=1;i<NX;i++)
		x[i] = x[i-1] + dx;

	/********************************************************
	* Distances for input array
	*********************************************************/
	double * xi = (double *)malloc(xdim*sizeof(double));
	double * yi = (double *)malloc(ydim*sizeof(double));	

	xi[0] = 0;
	yi[0] = 0;

	for(int i=1;i<xdim;i++)
		xi[i] = xi[i-1] + dlat*meters_per_lat;
	for(int i=1;i<ydim;i++)
		yi[i] = yi[i-1] + dlon*meters_per_lat;

	/********************************************************
	* Interpolate for each vertical level
	*********************************************************/
	for(int k=0;k<zdim;k++){
	
		ycount = 0;

		/********************************************************
		* Loop through j-points, looking for 
		*********************************************************/
		for(int j=0;j<NY;j++){

			/********************************************************
			* The y-coordinate of our point lies in between these two points of the grid
			* from which we are interpolating
			*********************************************************/
			if(y[j] > yi[ycount] && y[j] <= yi[ycount+1]){

				xcount = 0;
			
				for(int i=0;i<NX;i++){

					/********************************************************
					* The x-coordinate of our point lies in between these two points of the grid
					* from which we are interpolating
					*********************************************************/
					if(x[i] > xi[xcount] && x[i] <= xi[xcount+1]){

						/******************************************************** 
						* Fill the stencil p[4][4] with the correct points 
						* for interpolation
						*********************************************************/
						//---------------------------------------------------------
						// fractional distance between interpolating points
						mux = (x[i]-xi[xcount]) / (xi[xcount+1]-xi[xcount]);
						muy = (y[j]-yi[ycount]) / (yi[ycount+1]-yi[ycount]);
						//---------------------------------------------------------
						// domain interior of input field
						if(xcount != 0 && ycount != 0){

							for(int ax=-1;ax<3;ax++)
								for(int ay=-1;ay<3;ay++)
									p[ax+1][ay+1] = uz[d(xcount+ax,ycount+ay,k)];
						//---------------------------------------------------------
						// left side of input field
						} else if(xcount == 0 && ycount != 0){

							for(int ax=0;ax<3;ax++)
								for(int ay=-1;ay<3;ay++)
									p[ax+1][ay+1] = uz[d(xcount+ax,ycount+ay,k)];

							for(int ay=-1;ay<3;ay++)
								p[0][ay+1] = 2*uz[d(0,ycount+ay,k)] - uz[d(1,ycount+ay,k)];
						//---------------------------------------------------------
						// lower edge of input field
						} else if(ycount == 0 && xcount != 0){

							for(int ax=-1;ax<3;ax++)
								for(int ay=0;ay<3;ay++)
									p[ax+1][ay+1] = uz[d(xcount+ax,ycount+ay,k)];

							for(int ax=-1;ax<3;ax++)
								p[ax+1][0] = 2*uz[d(xcount+ax,0,k)] - uz[d(xcount+ax,1,k)];
						//---------------------------------------------------------
						// lower left corner
						} else if(ycount == 0 && xcount == 0){

							for(int ax=0;ax<3;ax++)
								for(int ay=0;ay<3;ay++)
									p[ax+1][ay+1] = uz[d(xcount+ax,ycount+ay,k)];

							for(int ax=0;ax<3;ax++)
								p[ax+1][0] = 2*uz[d(xcount+ax,0,k)] - uz[d(xcount+ax,1,k)];

							for(int ay=0;ay<3;ay++)
								p[0][ay+1] = 2*uz[d(0,ycount+ay,k)] - uz[d(1,ycount+ay,k)];

							p[0][0] = 0.5*(p[0][1] + p[1][0]);
						}

						out[d4(i,j,k)] = BicubicInterpolate2(p,mux,muy);

					/********************************************************
					* Our output grid is larger than the input grid in
					* its x dimension
					*********************************************************/
					} else if(xcount > xdim-1){

						i = NX;

					/********************************************************
					* Our point coincides exactly with the x-coordinate of 
					* the grid from which we are interpolating
					*********************************************************/
					} else if(x[i]==xi[xcount]){

						// y-coordinate also exactly coincides
						if(y[j]==yi[ycount]){

							out[d4(i,j,k)] = uz[d(xcount,ycount,k)];

						} else if(ycount!=0){

							mu = (y[j]-yi[ycount]) / (yi[ycount+1]-yi[ycount]);

							out[d4(i,j,k)] = CubicInterpolate2(
								uz[d(xcount,ycount-1,k)],
								uz[d(xcount,ycount,k)],
								uz[d(xcount,ycount+1,k)],
								uz[d(xcount,ycount+2,k)],mu);

						} else if(ycount==0){

							mu = (y[j]-yi[ycount]) / (yi[ycount+1]-yi[ycount]);

							out[d4(i,j,k)] = CubicInterpolate2(2*
								uz[d(xcount,ycount,k)]-uz[d(xcount,ycount+1,k)],
								uz[d(xcount,ycount,k)],uz[d(xcount,ycount+1,k)],
								uz[d(xcount,ycount+2,k)],mu);
						}

					} else { xcount++;i--;}
				}

			} else if(ycount > ydim-1){

				j = NY;			

			/********************************************************
			* Our point coincides exactly with the y-coordinate of 
			* the grid from which we are interpolating
			*********************************************************/
			} else if(y[j]==yi[ycount]){

				xcount = 0;

				for(int i=0;i<NX;i++){

					if(x[i] > xi[xcount] && x[i] <= xi[xcount+1]){

						out[d4(i,j,k)] = CubicInterpolate2(
							uz[d(xcount-1,ycount,k)],
							uz[d(xcount,ycount,k)],
							uz[d(xcount+1,ycount,k)],
							uz[d(xcount+2,ycount,k)],mu);
						
					} else if(xcount > xdim-1){

						i = NX;

					} else if(x[i]==xi[xcount]){

						out[d4(i,j,k)] = uz[d(xcount,ycount,k)];

					} else { xcount++;i--;}
				}

			} else { ycount++;j--;}
		}
	}

	free(x); free(y); free(xi); free(yi);

}

/*******************************************************************************
* WARNING MESSAGES IF OUTPUT GRID EXTENDS BEYOND INPUT GRID FOR VERTICAL
* INTERPOLATION
********************************************************************************/
void vertical_interpolation_warnings(double *z,double *zi,int zdim,int nz){
	
	static bool already_warned = false;
		
	if(!already_warned){
		if(z[nz-1] > zi[zdim-1] && INTERPOLATE_VERBOSE){
			printf("\n---------------------------------------------------------------------------------------------\n");
			printf("Warning!: Vertical interpolation cannot be carried out near line %d in file %s.\n Ending z point outside of data range. "
				"Data will be blank for %f km in the top.\n",__LINE__,__FILE__,fabs(zi[zdim-1]-z[nz-1])*1e-3);
			already_warned = true;
			printf("---------------------------------------------------------------------------------------------\n");
		}
	}
}

/*******************************************************************************
* WARNING MESSAGES IF OUTPUT GRID EXTENDS BEYOND INPUT GRID
********************************************************************************/
void horizontal_interpolation_warnings(double *x,double *xi,double *y,double *yi,int xdim,int nx,int ydim,int ny){
	
	static bool already_warned1 = false;
	static bool already_warned2 = false;
	static bool already_warned3 = false;
	static bool already_warned4 = false;

	//-----------------------------------------------------------------------------------------------------------------------------
	if(!already_warned1){
		if(y[0] < yi[0] && INTERPOLATE_VERBOSE){
			printf("\n---------------------------------------------------------------------------------------------\n");
			printf("Warning!: Horizontal interpolation cannot be carried out near line %d in file %s.\n Starting y point outside of data range. "
				"Data will be blank for %f km in the south.\n",__LINE__,__FILE__,fabs(yi[0]-y[0])*1e-3);
			already_warned1 = true;
			printf("---------------------------------------------------------------------------------------------\n");
		}
	}
	//-----------------------------------------------------------------------------------------------------------------------------	
	if(!already_warned2){
		if(x[0] < xi[0] && INTERPOLATE_VERBOSE){
			printf("\n---------------------------------------------------------------------------------------------\n");
			printf("Warning!: Horizontal interpolation cannot be carried out near line %d in file %s.\n Starting x point outside of data range. "
				"Data will be blank for %f km in the west.\n",__LINE__,__FILE__,fabs(xi[0]-x[0])*1e-3);
			printf("---------------------------------------------------------------------------------------------\n");
			already_warned2 = true;
		}
	}
	//-----------------------------------------------------------------------------------------------------------------------------	
	if(!already_warned3){
		if(y[ny-1] > yi[ydim-1] && INTERPOLATE_VERBOSE){
			printf("\n-----------------------------------------------------\n");
			printf("Warning!:Horizontal interpolation cannot be carried out near line %d in file %s.\n Ending y point outside of data range. "
				"Data will be blank for %f km in the north.\n",__LINE__,__FILE__,fabs(yi[ydim-1]-y[ny-1])*1e-3);
			printf("-----------------------------------------------------\n");
			already_warned3 = true;
		}
	}
	//-----------------------------------------------------------------------------------------------------------------------------
	if(!already_warned4){
		if(x[nx-1] > xi[xdim-1] && INTERPOLATE_VERBOSE){
			printf("\n-----------------------------------------------------\n");
			printf("Warning!:Horizontal interpolation cannot be carried out near line %d in file %s.\n Ending x point outside of data range. " 
				"Data will be blank for %f km in the east.\n",__LINE__,__FILE__,fabs(xi[xdim-1]-x[nx-1])*1e-3);
			printf("-----------------------------------------------------\n");
			already_warned4 = true;
		}
	}
	//-----------------------------------------------------------------------------------------------------------------------------
	
}

/***********************************************************************
* 
*************************************************************************/
void construct_tiles_upper_right(int xcount,int ycount,double *uz,int xdim,int ydim,int zdim,int k,double p[4][4]){
	
	if(xcount == xdim - 2 && ycount < ydim-2){

		for(int ax=-1;ax<2;ax++)
			for(int ay=-1;ay<3;ay++)
				p[ax+1][ay+1] = uz[d3(xcount+ax,ycount+ay,k)];

		for(int ay=-1;ay<3;ay++)
			p[3][ay+1] = 2*uz[d3(xdim-2,ycount+ay,k)] - uz[d3(xdim-1,ycount+ay,k)];
	}
	
}

/***********************************************************************
* 
*************************************************************************/
void construct_tiles_lower_left(int xcount,int ycount,double *uz,int xdim,int ydim,int zdim,int k,double p[4][4]){

	//---------------------------------------------------------
	// domain interior of input field
	//---------------------------------------------------------
	if(xcount != 0 && ycount != 0){

		for(int ax=-1;ax<3;ax++)
			for(int ay=-1;ay<3;ay++)
				p[ax+1][ay+1] = uz[d3(xcount+ax,ycount+ay,k)];
	//---------------------------------------------------------
	// left side of input field
	//---------------------------------------------------------
	} else if(xcount == 0 && ycount != 0){

		for(int ax=0;ax<3;ax++)
			for(int ay=-1;ay<3;ay++)
				p[ax+1][ay+1] = uz[d3(xcount+ax,ycount+ay,k)];

		for(int ay=-1;ay<3;ay++)
			p[0][ay+1] = 2*uz[d3(0,ycount+ay,k)] - uz[d3(1,ycount+ay,k)];
	//---------------------------------------------------------
	// lower edge of input field
	//---------------------------------------------------------	
	} else if(ycount == 0 && xcount != 0){

		for(int ax=-1;ax<3;ax++)
			for(int ay=0;ay<3;ay++)
				p[ax+1][ay+1] = uz[d3(xcount+ax,ycount+ay,k)];
		
		for(int ax=-1;ax<3;ax++)
			p[ax+1][0] = 2*uz[d3(xcount+ax,0,k)] - uz[d3(xcount+ax,1,k)];				
	//---------------------------------------------------------
	// lower left corner
	//---------------------------------------------------------	
	} else if(ycount == 0 && xcount == 0){

		for(int ax=0;ax<3;ax++)
			for(int ay=0;ay<3;ay++)
				p[ax+1][ay+1] = uz[d3(xcount+ax,ycount+ay,k)];

		for(int ax=0;ax<3;ax++)
			p[ax+1][0] = 2*uz[d3(xcount+ax,0,k)] - uz[d3(xcount+ax,1,k)];

		for(int ay=0;ay<3;ay++)
			p[0][ay+1] = 2*uz[d3(0,ycount+ay,k)] - uz[d3(1,ycount+ay,k)];

		p[0][0] = 0.5*(p[0][1] + p[1][0]);
	}

}

/***********************************************************************
* Interpolate the horizontal dimensions of an array
*
* *uz - input array with dimension xdim x ydim x NZ
* out[NX][NY][NZ] - output array
* xstagger - does it have a staggered x-coordinate?
* ystagger - does it have a stagger y-coordinate?
*
*************************************************************************/
void horz_interpolate_from_model(double *uz,double *out,
								int xdim,int ydim,int zdim,
								bool xstagger,bool ystagger,
								double dlat,double dlon,
								double lonoffset,double latoffset){

	int xcount = 0, ycount = 0;

	double p[4][4], mu, mux, muy;

	bool xfail,yfail = false;

	double tolerance = dx * 1.0e-4;

	/********************************************************
	* Distances in meters for output interpolated array
	*********************************************************/
	double * x = (double *)malloc(NX*sizeof(double));
	double * y = (double *)malloc(NY*sizeof(double));

	x[0] = lonoffset*meters_per_lat;
	y[0] = latoffset*meters_per_lat;

	if(xstagger){ x[0] = x[0] - 0.5*dx;}

	if(ystagger){ y[0] = y[0] - 0.5*dy;}

	for(int i=1;i<NY;i++)
		y[i] = y[i-1] + dy;
	for(int i=1;i<NX;i++)
		x[i] = x[i-1] + dx;

	/********************************************************
	* Distances for input array
	*********************************************************/
	double * xi = (double *)malloc(xdim*sizeof(double));
	double * yi = (double *)malloc(ydim*sizeof(double));	

	xi[0] = 0;
	yi[0] = 0;

	for(int i=1;i<xdim;i++)
		xi[i] = xi[i-1] + dlat*meters_per_lat;
	for(int i=1;i<ydim;i++)
		yi[i] = yi[i-1] + dlon*meters_per_lat;


	for(int i=0;i<xdim && i<NX;i++){

		if(fabs(xi[i]-x[i]) < tolerance ){ xi[i] = x[i];}
	 }
	 

	for(int i=0;i<ydim && i<NY;i++){

		if(fabs(yi[i]-y[i]) < tolerance ){ yi[i] = y[i];}
	 }
 
	//---------------------------------------------------------------
	// WARNING MESSAGES IF OUTPUT GRID EXTENDS BEYOND INPUT GRID
	//---------------------------------------------------------------
	horizontal_interpolation_warnings(x,xi,y,yi,xdim,NX,ydim,NY);
 
	//--------------------------------------------------------------
	// If the output grid starts to the left or south of the input
	// grid, handle by resetting the index until a matching point
	// is found. The variables istart, jstart, ifound, and jfound
	// control this process.
	//--------------------------------------------------------------
	int istart = 0;
	int jstart = 0;

	bool jfound = false;
	bool ifound = false;

	/********************************************************
	* Interpolate for each vertical level
	*********************************************************/
	for(int k=0;k<zdim;k++){
		
		ycount = 0;
		jfound = false;
		jstart = 0;
		
		/********************************************************
		* Loop through j-points, looking for 
		*********************************************************/
		for(int j=jstart;j<NY;j++){
			//printf("%d y %f yi %f\n",j,y[j],yi[ycount]);
			/********************************************************
			* The y-coordinate of our point lies in between these two 
			* points of the grid from which we are interpolating
			*********************************************************/
			if(ycount+1 < ydim && (y[j] > yi[ycount] && y[j] < yi[ycount+1]) ){
				//printf("%d matched y %f yi %f %f\n",j,y[j],yi[ycount],yi[ycount+1]);
				xcount = 0;
			
				jfound = true;	// we've reached the beginning of the grid overlap in the y-direction
				ifound = false;
			
				istart = 0;
			
				for(int i=istart;i<NX;i++){
						//printf("%d %d %d %d %d %f %f %f %f %f %f %f \n ",i,j,k,xcount,ycount,out[d4(i,j,k)],y[j],yi[ycount],yi[ycount+1],x[i],xi[xcount],xi[xcount+1]);
					/********************************************************
					* The x-coordinate of our point lies in between these two 
					* points of the grid from which we are interpolating
					*********************************************************/
					if(xcount>-1 && xcount < xdim - 1 && x[i] > xi[xcount] && x[i] < xi[xcount+1]){

						ifound = true; 	// we've reached the beginning of the grid overlap in the y-direction
	//printf("%d %d %d %d %d %f %f %f %f %f %f %f \n ",i,j,k,xcount,ycount,out[d4(i,j,k)],y[j],yi[ycount],yi[ycount+1],x[i],xi[xcount],xi[xcount+1]);
						/******************************************************** 
						* Fill the stencil p[4][4] with the correct points 
						* for interpolation
						*********************************************************/
						//---------------------------------------------------------
						// fractional distance between interpolating points
						mux = (x[i]-xi[xcount]) / (xi[xcount+1]-xi[xcount]);
						muy = (y[j]-yi[ycount]) / (yi[ycount+1]-yi[ycount]);
						
						for(int ax=0;ax<4;ax++)
							for(int ay=0;ay<4;ay++)
								p[ax][ay] = 0;
						
						//------------------------------------------------------------
						// Make sure we are not on the upper or right edge of field
						//------------------------------------------------------------
						if(xcount+2 < xdim && ycount+2 < ydim){
						
							construct_tiles_lower_left(xcount,ycount,uz,xdim,ydim,zdim,k,p);
							
							
												
						} else {
								//if(xcount == xdim - 2 && ycount < ydim-2){
							//printf("%d %d %d %d %d %f %f %f %f %f %f %f \n ",i,j,k,xcount,ycount,out[d4(i,j,k)],y[j],yi[ycount],yi[ycount+1],x[i],xi[xcount],xi[xcount+1]);
									//}
							construct_tiles_upper_right(xcount,ycount,uz,xdim,ydim,zdim,k,p);
						}
						// should this depend on whether or not the previuos
						// if clauses were satisfied?
						out[d4(i,j,k)] = BicubicInterpolate2(p,mux,muy);

					} else if(xcount==xdim-1){
	//printf("%d %d %d %d %d %f %f %f %f %f %f %f \n ",i,j,k,xcount,ycount,out[d4(i,j,k)],y[j],yi[ycount],yi[ycount+1],x[i],xi[xcount],xi[xcount+1]);
						if(x[i]-xi[xdim-1]<10){ out[d4(i,j,k)] = uz[d3(xcount,ycount,k)];}
						//printf("%d %d %d %d %d %f %f %f %f %f %f %f \n ",i,j,k,xcount,ycount,out[d4(i,j,k)],y[j],yi[ycount],yi[ycount+1],x[i],xi[xcount],xi[xcount+1]);
						
						if(!ifound && i<NX){ istart++; i = istart; xcount = 0;}
						
					/********************************************************
					* Our output grid is larger than the input grid in
					* its x dimension
					*********************************************************/
					} else if(xcount > xdim-1){ // does this ever get executed?
							
						if(!ifound && i<NX){ istart++; i = istart; xcount = 0;}
						//i = NX;

					/********************************************************
					* Our point coincides exactly with the x-coordinate of 
					* the grid from which we are interpolating
					*********************************************************/
					} else if(x[i]==xi[xcount]){
						
						ifound = true;
						
						//printf("%d %d %d %f %f %f %f ||| ",i,j,k,y[j],yi[ycount],x[i],xi[xcount]);
						// y-coordinate also exactly coincides
						if(y[j]==yi[ycount]){

							out[d4(i,j,k)] = uz[d3(xcount,ycount,k)];
							
						} else if(ycount!=0 && ycount+2<ydim){

							mu = (y[j]-yi[ycount]) / (yi[ycount+1]-yi[ycount]);
							
							out[d4(i,j,k)] = CubicInterpolate2(
								uz[d3(xcount,ycount-1,k)],
								uz[d3(xcount,ycount,k)],
								uz[d3(xcount,ycount+1,k)],
								uz[d3(xcount,ycount+2,k)],mu);	// THIS COULD FAIL ON THE RIGHTMOST EDGE, CONSIDER FIXING IT!

						} else if(ycount==0 && ycount+2<ydim){

							mu = (y[j]-yi[ycount]) / (yi[ycount+1]-yi[ycount]);
							
							out[d4(i,j,k)] = CubicInterpolate2(2*
								uz[d3(xcount,ycount,k)]-uz[d3(xcount,ycount+1,k)],
								uz[d3(xcount,ycount,k)],uz[d3(xcount,ycount+1,k)],
								uz[d3(xcount,ycount+2,k)],mu);
								
						}

					} else { xcount++;i--;}
				}

			} else if(ycount > ydim-1){

				if(!jfound && j<NY){ jstart++; j = jstart; ycount = 0;}

			/********************************************************
			* Our point coincides exactly with the y-coordinate of 
			* the grid from which we are interpolating
			*********************************************************/
			} else if(y[j]==yi[ycount]){
						
				xcount = 0;
				ifound = false;
				istart = 0;

				for(int i=istart;i<NX;i++){
					
					if(xcount+2<xdim && xcount-1>-1 && x[i] > xi[xcount] && x[i] < xi[xcount+1]){
						//printf("%d %d %d %f %f %f %f %f ||| ",i,j,k,out[d4(i,j,k)],y[j],yi[ycount],x[i],xi[xcount]);
						ifound = true;
						
						out[d4(i,j,k)] = CubicInterpolate2(
							uz[d3(xcount-1,ycount,k)],
							uz[d3(xcount,ycount,k)],
							uz[d3(xcount+1,ycount,k)],
							uz[d3(xcount+2,ycount,k)],mu); 	// THIS COULD FAIL ON THE RIGHTMOST EDGE, CONSIDER FIXING IT!
						
					} else if(xcount > xdim-1){
						//printf("%d %d %d %d %d %f %f %f %f %f \n ",i,j,k,xcount,ycount,out[d4(i,j,k)],y[j],yi[ycount],x[i],xi[xcount]);
						if(!ifound && i<NX){ istart++; i = istart; xcount = 0;}
						
					} else if(x[i]==xi[xcount]){
						// printf("%d\n",k);
						out[d4(i,j,k)] = uz[d3(xcount,ycount,k)];
						
						ifound = true;
						//if(out[d4(i,j,k)]>1e-5){
						//printf("%d %d %d %d %d %f %f %f %f %f \n ",i,j,k,xcount,ycount,out[d4(i,j,k)],y[j],yi[ycount],x[i],xi[xcount]);
							//}
					} else if(xcount==0 && x[i] > xi[xcount] && x[i] < xi[xcount+1]){
						
						out[d4(i,j,k)] = CubicInterpolate2(2*
							uz[d3(xcount,ycount,k)]-uz[d3(xcount+1,ycount,k)],
							uz[d3(xcount,ycount,k)],uz[d3(xcount+1,ycount,k)],
							uz[d3(xcount+2,ycount,k)],mu);
						
						ifound = true;
						
					} else if(xcount==xdim-2 && x[i] > xi[xcount] && x[i] < xi[xcount+1]){
						//printf("%d %f %f %f\n",i,x[i],xi[xcount],xi[xcount+1]);
						out[d4(i,j,k)] = CubicInterpolate2(uz[d3(xcount-1,ycount,k)],
							uz[d3(xcount,ycount,k)],uz[d3(xcount+1,ycount,k)],
							2*uz[d3(xcount+1,ycount,k)]-uz[d3(xcount,ycount,k)],mu);
						
						ifound = true;
						
					} else { xcount++;i--;}
				}

			} else { ycount++;j--;}
		}
	}


	free(x); free(y); free(xi); free(yi);

}

/********************************************************
*
* 
*
*********************************************************/
double CubicInterpolate(double y0,double y1,double y2,double y3,double x0,double x1,double x2,double x3,double mu){
   
	double a,b,c,d,f0p,f1p;

	f0p = (y2-y0)/((x2-x0)/(x2-x1));
	f1p = (y3-y1)/((x3-x0)/(x2-x1));

	a = 2*y1 - 2*y2 +  f0p + f1p;
	b = -3*y1 + 3*y2 - 2*f0p - f1p;
	c = f0p;
	d = y1;
	
	return a*mu*mu*mu + b*mu*mu + c*mu + d;
}

/********************************************************
*
* 
*
*********************************************************/
double CubicInterpolate2(double y0,double y1,double y2,double y3,double mu){
   
	double a,b,c,d,f0p,f1p;

	f0p = (y2-y0)/2;
	f1p = (y3-y1)/2;

	a = 2*y1 - 2*y2 +  f0p + f1p;
	b = -3*y1 + 3*y2 - 2*f0p - f1p;
	c = f0p;
	d = y1;
	
	return a*mu*mu*mu + b*mu*mu + c*mu + d;
}

/********************************************************
*
* 
*
*********************************************************/
double BicubicInterpolate2(double p[4][4],double mux,double muy){
   
	double p1,p2,p3,p4;

	p1 = CubicInterpolate2(p[0][0],p[1][0],p[2][0],p[3][0],mux);
	p2 = CubicInterpolate2(p[0][1],p[1][1],p[2][1],p[3][1],mux);
	p3 = CubicInterpolate2(p[0][2],p[1][2],p[2][2],p[3][2],mux);
	p4 = CubicInterpolate2(p[0][3],p[1][3],p[2][3],p[3][3],mux);

	return CubicInterpolate2(p1,p2,p3,p4,muy);
}

/********************************************************
*
* 
*
*********************************************************/
double CosineInterpolate(double y1,double y2,double mu){
   
	double mu2;

	mu2 = (1-cos(mu*3.1416))/2;

	return(y1*(1-mu2)+y2*mu2);
}

/********************************************************
*
* 
*
*********************************************************/
double LinearInterpolate(double y1,double y2,double mu){

   return(y1*(1-mu)+y2*mu);
}

