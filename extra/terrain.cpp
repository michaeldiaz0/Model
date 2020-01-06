#include "stdafx.h"
#include "interpolate.h"
#include "advection.h"
#include "terrain.h"

double zheight[NX][NY][NZ];	// physical height of coordinate surfaces

//double divgbar[NX][NY][NZ];	// basic state divergence
//double divg[NX][NY][NZ];	// perturbation state divergence

double sdot[NX][NY][NZ];
double sdotbar[NX][NY][NZ];

double tbv3d[NX][NY][NZ];
double tb3d[NX][NY][NZ];
double qb3d[NX][NY][NZ];

double grhoc3d[NX][NY][NZ];	// density/Gz at cell corners
double grhou3d[NX][NY][NZ];	// density/Gz at u-points
double grhov3d[NX][NY][NZ];	// density/Gz at v-points
double grhos3d[NX][NY][NZ];	// density/Gz at scalar points

double rhow3d[NX][NY][NZ];	// density at w-points

double Gx[NX][NY][NZ];
double Gy[NX][NY][NZ];

double Gz[NX][NY];
double Gzu[NX][NY];
double Gzv[NX][NY];
double one_d_Gz[NX][NY];
double one_d_Gzu[NX][NY];
double one_d_Gzv[NX][NY];

double grhoavg2d[NX][NY];

double sdot_at_uf[2];
double sdotb_at_uf[2];
double sdot_at_vf[2];
double sdotb_at_vf[2];

/*****************************************************************************************
*
******************************************************************************************/
void compute_height(double z[NZ]){

	double H = zu[NZ-1];

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		zheight[i][j][k] =  z[k]*( H - topo[i][j][1] ) / H + topo[i][j][1];

	}}}
}

/*****************************************************************************************
*
******************************************************************************************/
void compute_height_u(){

	double H = zu[NZ-1];
	double height;

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		height = TOPO_U;

		zheight[i][j][k] =  zu[k]*( H - height ) / H + height;

	}}}
}

/*****************************************************************************************
*
******************************************************************************************/
void compute_height_v(){

	double H = zu[NZ-1];
	double height;

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		height = TOPO_V;

		zheight[i][j][k] =  zu[k]*( H - height ) / H + height;

	}}}
}

/*****************************************************************************************
*
******************************************************************************************/
void compute_height_c(){

	double H = zu[NZ-1];
	double height;

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		height = TOPO_C;

		zheight[i][j][k] =  zu[k]*( H - height ) / H + height;

	}}}
}

/*****************************************************************************************
*
******************************************************************************************/
void sdot_interp(int i,int j,int k){

	sdot_at_uf[0] = 0.5* (sdot[i][j][k  ]+sdot[i-1][j][k  ]);
	sdot_at_uf[1] = 0.5* (sdot[i][j][k+1]+sdot[i-1][j][k+1]);

	sdot_at_vf[0] =  0.5*(sdot[i][j][k  ]+sdot[i][j-1][k  ]);
	sdot_at_vf[1] =  0.5*(sdot[i][j][k+1]+sdot[i][j-1][k+1]);
}

/*****************************************************************************************
*
******************************************************************************************/
void compute_sigma_dot(double u[NX][NY][NZ],double v[NX][NY][NZ]){

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

//		sdot[i][j][k] = 0.25*(u[i+1][j][k] + u[i][j][k] + u[i+1][j][k-1] + u[i][j][k-1])*Gx[i][j][k]*one_d_Gz[i][j] 
//					  + 0.25*(v[i][j+1][k] + v[i][j][k] + v[i][j+1][k-1] + v[i][j][k-1])*Gy[i][j][k]*one_d_Gz[i][j]
//					  + w[i][j][k];

		w[i][j][k] = sdot[i][j][k] - (
						0.25*(u[i+1][j][k] + u[i][j][k] + u[i+1][j][k-1] + u[i][j][k-1])*Gx[i][j][k]*one_d_Gz[i][j] 
					  + 0.25*(v[i][j+1][k] + v[i][j][k] + v[i][j+1][k-1] + v[i][j][k-1])*Gy[i][j][k]*one_d_Gz[i][j]  
					);
	}}}
}

/*****************************************************************************************
*
******************************************************************************************/
void compute_divergence(){

	

}

/*****************************************************************************************
*
******************************************************************************************/
void compute_density_terms(){

	size_t num_bytes = NX*NY*NZ*sizeof(double);

	memcpy(u,ubar,num_bytes);
	memcpy(v,vbar,num_bytes);
	memcpy(th,thbar,num_bytes);

	compute_height(zu);
	vert_interpolate_1d(zheight,rhou,grhos3d);
	vert_interpolate_1d(zheight,tbv,tbv3d);
	vert_interpolate_1d(zheight,tb,tb3d);
	vert_interpolate_1d(zheight,qb,qb3d);
	vert_interpolate_3d(zheight,th,thbar);


	compute_height(zw);
	vert_interpolate_1d(zheight,rhou,rhow3d);

	compute_height_u();
	vert_interpolate_1d(zheight,rhou,grhou3d);
	vert_interpolate_3d(zheight,u,ubar);

	compute_height_v();
	vert_interpolate_1d(zheight,rhou,grhov3d);
	vert_interpolate_3d(zheight,v,vbar);

	compute_height_c();
	vert_interpolate_1d(zheight,rhou,grhoc3d);


	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		sdotbar[i][j][k] = 0;
		sdot[i][j][k] = 0;
		u[i][j][k] = 0;
		v[i][j][k] = 0;
		th[i][j][k] = 0;
		grhos3d[i][j][k] = grhos3d[i][j][k] / Gz[i][j];
		grhou3d[i][j][k] = grhou3d[i][j][k] / Gzu[i][j];
		grhov3d[i][j][k] = grhov3d[i][j][k] / Gzv[i][j];

	}}}

	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=0;k<NZ;k++){

		grhoc3d[i][j][k] = grhoc3d[i][j][k] / (0.25*(Gz[i][j]+Gz[i-1][j]+Gz[i][j-1]+Gz[i-1][j-1]));

	}}}

	// boundary condition
	for(int k=0;k<NZ;k++){
		for(int i=0;i<NX;i++){ grhoc3d[i][0][k] = grhoc3d[i][1][k]; grhoc3d[i][NY-1][k] = grhoc3d[i][NY-2][k];}
		for(int j=0;j<NY;j++){ grhoc3d[0][j][k] = grhoc3d[1][j][k]; grhoc3d[NX-1][j][k] = grhoc3d[NX-2][j][k];}
	}

	/*********************************************************************
	* Column averaged density
	**********************************************************************/
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	
		grhoavg2d[i][j] = 0;

		for(int k=1;k<NZ-1;k++){ grhoavg2d[i][j] = grhoavg2d[i][j] + grhos3d[i][j][k]*tbv3d[i][j][k];}

		grhoavg2d[i][j] = grhoavg2d[i][j]/((double)(NZ-2));

	}}

//	memset(u,0,num_bytes);
//	memset(v,0,num_bytes);
//	memset(th,0,num_bytes);

	sdot_at_uf[0] = 0;
	sdotb_at_uf[0] = 0;
	sdot_at_vf[0] = 0;
	sdotb_at_vf[0] = 0;
	sdot_at_uf[1] = 0;
	sdotb_at_uf[1] = 0;
	sdot_at_vf[1] = 0;
	sdotb_at_vf[1] = 0;

	for(int i=0;i<NX;i++)
		for(int j=0;j<NY;j++)
			for(int k=0;k<NZ;k++)
				friction[i][j][k] = 0;

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		friction[i][j][1] = 1.0e-5;
	}}

}

/*****************************************************************************************
*
******************************************************************************************/
void compute_metric_terms(){

	double H = zu[NZ-1];

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		Gz[i][j] = H / ( H - topo[i][j][1] );
		Gzu[i][j] = H / ( H - TOPO_U);
		Gzv[i][j] = H / ( H - TOPO_V);

		one_d_Gz[i][j] = 1 / Gz[i][j];
	}}


	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){

		for(int k=0;k<NZ;k++){
	
			Gx[i][j][k] = (zw[k] - H) / ( H - topo[i][j][1] ) * ( topo[i+1][j][1] - topo[i-1][j][1] ) / (2.*dx);

			Gy[i][j][k] = (zw[k] - H) / ( H - topo[i][j][1] ) * ( topo[i][j+1][1] - topo[i][j-1][1] ) / (2.*dy);

		}
	}}

	for(int k=1;k<NZ-1;k++){
		for(int i=0;i<NX;i++){ Gx[i][0][k] = Gx[i][1][k]; Gx[i][NY-1][k] = Gx[i][NY-2][k];}
		for(int j=0;j<NY;j++){ Gx[0][j][k] = Gx[1][j][k]; Gx[NX-1][j][k] = Gx[NX-2][j][k];}
		for(int i=0;i<NX;i++){ Gy[i][0][k] = Gy[i][1][k]; Gy[i][NY-1][k] = Gy[i][NY-2][k];}
		for(int j=0;j<NY;j++){ Gy[0][j][k] = Gy[1][j][k]; Gy[NX-1][j][k] = Gy[NX-2][j][k];}
	}

	printf("Gz = %f\n",Gz[50][60]);

}



/*****************************************************************************************
* 
******************************************************************************************/
//int main(){

//	initialize2();


//	compute_metric_terms();

//	compute_density_terms();

//	compute_height_u();

//	for(int k=0;k<NZ;k++){

//		//printf("%d %f %f %f %f %f %f\n",k,zu[k],rhou[k],zheight[70][52][k],zheight[70][53][k],grhos3d[70][53][k],grhov3d[70][52][k]);

//		printf("%d %f %f %f %f\n",k,zu[k],zheight[70][52][k],u[70][52][k],ubar[70][52][k]);

//	}

//	printf("%I64u\n",NX*NY*NZ*sizeof(double));

//}

//	for(int k=0;k<NZ;k++){
//	for(int j=0;j<NY;j++){

//		printf("%f ",Gx[70][j][k]*100000);

//	}
//			printf("\n");
//	}


//	//printf("%f %f\n",topo[70][55][1],zheight[70][55][10]);


//	for(int k=0;k<NZ;k++){
//	for(int j=0;j<NY;j++){

//		printf("%f ",zheight[70][j][k]);

//	}
//			printf("\n");
//	}
