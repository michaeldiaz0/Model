
#include "stdafx.h"
#include "Heating.h"
#include "interpolate.h"
#include "pcomm.h"

//using namespace std;

/********************************************************
*
*
*
*
*********************************************************/



// normalized heating rates in degrees per day per centimeter
double base_heating[] = {0,0.25,0.75,1.25,1.75,2.25,2.75,3.19,3.56,
						3.89,4.16,4.35,4.45,4.48,4.43,4.25,3.95,3.73,
						3.58,3.38,3.13,2.88,2.63, 2.38,2.13,1.8,1.4,
						1.125,0.975,0.8,0.6,0.425,0.275,0.175,0.125};

//double base_heating[] = {0,6.0,0.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

double z_heating[35]; 

// heating rate multiplied by precipitation rate 
// (degrees per model time step)
double *heating;//[NZ];

// end grid points of heating region
int i_eastlon,i_eastlat,i_westlon,i_westlat,i_width;

// end lat/lon points of heating region
double d_eastlon,d_eastlat,d_westlon,d_westlat,d_width;

double equiv_precip_rate;

double slope;

int * arr;

/********************************************************
* Class constructor
*********************************************************/
Heating::Heating(){

}

/********************************************************
* Class Initializer
*********************************************************/
void Heating::initialize(double lat1, double lon1, double lat2, double lon2, double width_meters, double equiv_precip_rate){

	heating = (double *)calloc(NZ,sizeof(double));

	arr = (int *)malloc(NX*sizeof(int));

	for(int k=0;k<NZ;k++){ heating[k] = 0;}

	for(int k=0;k<35;k++){
	
		z_heating[k] = ((double)k-0.5)*dz;
	}

	this->equiv_precip_rate = equiv_precip_rate;

	setHeatingRate(equiv_precip_rate);

	d_eastlon = lon1;
	d_eastlat = lat1;
	d_westlon = lon2;
	d_westlat = lat2;
	d_width = width_meters;

	// need to test outLats and outLons to make sure they're set
	
	i_eastlat = get_point_from_lat(lat1);
	i_westlat = get_point_from_lat(lat2);
	i_eastlon = get_point_from_lon(lon1);
	i_westlon = get_point_from_lon(lon2);

	slope = (d_eastlat - d_westlat) / (d_eastlon - d_westlon);

	i_width = (int)(width_meters/dy+0.1);

	line_array(i_eastlon,i_eastlat,i_westlon,i_westlat,arr);
	

}

/********************************************************
* 
*********************************************************/
void Heating::scaleHeating(double scaleFactor){

	for(int k=0;k<NZ;k++){ heating[k] *= scaleFactor;}
	
	this->equiv_precip_rate *= scaleFactor;

}

/********************************************************
* 
*********************************************************/
int Heating::changeSize(int xl,int xh,int w){


	this->i_eastlon += xh;
	this->i_westlon += xl;
	
	this->i_width += w;


	this->d_eastlat += (outLons[this->i_eastlon]-d_eastlon) * slope;
	this->d_westlat += (outLons[this->i_westlon]-d_westlon) * slope;

	//printf("slope = %f %f %f\n",slope,d_eastlat,d_westlat);

	this->d_eastlon = outLons[this->i_eastlon];
	this->d_westlon = outLons[this->i_westlon];

	this->i_eastlat = get_point_from_lat(this->d_eastlat);
	this->i_westlat = get_point_from_lat(this->d_westlat);

	line_array(i_eastlon,i_eastlat,i_westlon,i_westlat,arr);
	
	
	
	return getSizeInGridPoints();
}

/********************************************************
* 
*********************************************************/
int Heating::getSizeInGridPoints(){

	int size = 0;
	int width;

	for(int i=i_eastlon;i<i_westlon;i++){
		
		width = 0;
		
		for(int j=arr[i]-i_width;j<arr[i]+i_width;j++){
			size++;
			width++;
		}
		
		//printf("%d\n",width);
	}
	
	return size;
}

/********************************************************
* 
*********************************************************/
void Heating::setLocation(double lat1,double lon1,double lat2,double lon2){

	i_eastlat = get_point_from_lat(lat1);
	i_westlat = get_point_from_lat(lat2);
	i_eastlon = get_point_from_lon(lon1);
	i_westlon = get_point_from_lon(lon2);

	d_eastlon = lon1;
	d_eastlat = lat1;
	d_westlon = lon2;
	d_westlat = lat2;

	line_array(i_eastlon,i_eastlat,i_westlon,i_westlat,arr);
}

/********************************************************
* 
*********************************************************/
void Heating::setLocationGrid(int lat1,int lon1,int lat2,int lon2){

	i_eastlat = lat1;
	i_westlat = lat2;
	i_eastlon = lon1;
	i_westlon = lon2;

	d_eastlat = outLats[lat1];
	d_westlat = outLats[lat2];
	d_eastlon = outLons[lon1];
	d_westlon = outLons[lon2];

	line_array(i_eastlon,i_eastlat,i_westlon,i_westlat,arr);
}

/********************************************************
* 
*********************************************************/
void Heating::shift(double latshift,double lonshift){

	setLocation(d_eastlat+latshift, d_eastlon+lonshift, d_westlat+latshift, d_westlon+lonshift);

}

/********************************************************
* 
*********************************************************/
void Heating::shiftGrid(int x, int y){

	i_eastlat = i_eastlat + y;
	i_westlat = i_westlat + y;
	i_eastlon = i_eastlon + x;
	i_westlon = i_westlon + x;

	d_eastlon = outLons[i_eastlon];
	d_eastlat = outLats[i_eastlat];
	d_westlon = outLons[i_westlon];
	d_westlat = outLats[i_westlat];
}

/*********************************************************************
*
*
**********************************************************************/
void Heating::applyHeating(){

	for(int i=i_eastlon;i<i_westlon;i++){
	for(int j=arr[i]-i_width;j<arr[i]+i_width;j++){
	for(int k=0;k<NZ;k++){
		THP(i,j,k) += heating[k];
	}}}
	
	//printf("%f\n",heating[10]);

}

/*********************************************************************
*
*
**********************************************************************/
void Heating::applyHeating_random(){

	double r;

	for(int i=i_eastlon;i<i_westlon;i++){
	for(int j=arr[i]-i_width;j<arr[i]+i_width;j++){
	for(int k=0;k<NZ;k++){
		
		r = (double)rand() / (double)RAND_MAX;

		//printf("%f ",r);

		if(r>0.98){ THP(i,j,k) += heating[k];}
	}}}
}
#if PARALLEL
/*********************************************************************
*
*
**********************************************************************/
void Heating::p_applyHeating(){

	for(int i=i_eastlon;i<i_westlon;i++){
	for(int j=arr[i]-i_width;j<arr[i]+i_width;j++){
	for(int k=0;k<NZ;k++){

		if( i>=ibs[rank] && i<=ibs[rank]+myNX && j>=jbs[rank] && j<=jbs[rank]+myNY ){

			THP(i-ibs[rank]+3,j-jbs[rank]+3,k) += heating[k];	
		}
	}}}
	
}

/*********************************************************************
*
*
**********************************************************************/
void Heating::p_applyHeating_random(){

	double r;
	int sum1,sum2;
	
	sum1 = 0;
	sum2 = 0;

	//for(int i=3;i<fNX-3;i++){
	//for(int j=3;j<fNY-3;j++){
	for(int i=i_eastlon;i<i_westlon;i++){
	for(int j=arr[i]-i_width;j<arr[i]+i_width;j++){
	for(int k=0;k<NZ;k++){

		if( i>=ibs[rank] && i<=ibs[rank]+myNX && j>=jbs[rank] && j<=jbs[rank]+myNY ){

			r = (double)rand() / (double)RAND_MAX;
			
			if(r>0.9){
				//THP(i,j,k) += heating[k];
				THP(i-ibs[rank]+3,j-jbs[rank]+3,k) += heating[k];
				sum1++;
			} else {
				sum2++;
			}
		}
	}}}
	
	//printf("%d %d %f\n",sum1,sum2,(double)sum1/(double)(sum1+sum2));
	//fflush(stdout);
	
}
#endif
/*********************************************************************
*
*
**********************************************************************/
void Heating::heating_oval(int xpos,int ypos,double iradius,double jradius){

	int il,ih,jl,jh;

	double length;	

	ih = xpos + (int)(iradius/dx)+1;
	il = xpos - (int)(iradius/dx);

	jh = ypos + (int)(jradius/dy)+1;
	jl = ypos - (int)(jradius/dy);

	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		length = (double)((i-xpos)*(i-xpos)) / (iradius*iradius) + (double)((j-ypos)*(j-ypos)) / (jradius*jradius);

		if(length <= 1){

			for(int k=0;k<NZ;k++){ THP(i,j,k) = THP(i,j,k) + heating[k];}
		}
		
	}}

}

/********************************************************
* Initialize heating rate array
*********************************************************/
void Heating::setHeatingRate(double rain_rate){

	if(STRETCHED_GRID){
		
		vert_interpolate_1d(34,NZ,&z_heating[0],&zsu[0],&base_heating[0],&heating[0]);
		
		for(int k=0;k<NZ;k++){ heating[k] = rain_rate * dt * heating[k] / (24.*60.*60.);}
		
	} else {

		for(int k=0;k<35;k++){ heating[k] = rain_rate * dt * base_heating[k] / (24.*60.*60.);}

		for(int k=34;k<NZ;k++){ heating[k] = 0;}
	}

	//for(int k=0;k<NZ;k++){

//		if(k%2==0){ heating[k] = rain_rate * dt * base_heating[k/2] / (24.f*60.f*60.f);}
//		else { heating[k] = rain_rate * dt * 0.5*(base_heating[k/2]+base_heating[k/2+1]) / (24.f*60.f*60.f);}
		//printf("%d %f\n",k,heating[k]);
	//}
				
	//for(int k=0;k<NZ;k++){
	//	printf("%d %f %f\n",k,heating[k],zsu[k]);
	//}
//	
//	heating[69] = 0;

}

/********************************************************
* 
*********************************************************/
void Heating::printInfo(FILE * infile){

	if(infile==NULL){

		printf("Heating %f %f %f %f %d %f\n",outLats[i_eastlat],outLons[i_eastlon],outLats[i_westlat],outLons[i_westlon],i_width*2,equiv_precip_rate);

		printf("%d %d %d %d\n",i_eastlat,i_eastlon,i_westlat,i_westlon);

	} else {

		fprintf(infile,"Heating %f %f %f %f %d %f\n",outLats[i_eastlat],outLons[i_eastlon],outLats[i_westlat],outLons[i_westlon],i_width*2,equiv_precip_rate);

		fprintf(infile,"%d %d %d %d\n",i_eastlat,i_eastlon,i_westlat,i_westlon);
	}

}

