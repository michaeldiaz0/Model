#include "stdafx.h"

double shear,buoyancy,dissipation;
double bruntv;
double Kdiffhs,Kdiffvs;
double Kdiffh,Kdiffv;
double lh,lvv,l;
double C;
double Ck = .24;
double D12,D13,D23,D11,D22,D33;
double theta_diff;

double dx2;
double dy2;
double dz2;

double grid_cubicroot;

double tke[NX][NY][NZ],tkep[NX][NY][NZ],tkem[NX][NY][NZ];

/********************************************************	
* Prognostic TKE equation 
*
*********************************************************/
void integrate_tke(double step){

	grid_cubicroot = pow(dx*dy*dz,1./3.);

	printf("cubic = %f\n",grid_cubicroot);

	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=1;k<NZ-1;k++){
	
		/********************************************************	
		* Brunt-Vaisala Frequency (dry air)
		*********************************************************/
	
		theta_diff = (tb[k+1]+th[i][j][k+1]-tb[k-1]-th[i][j][k-1])/2;
	
		bruntv = grav*(1/(tb[k]+th[i][j][k]))*(theta_diff)/dz;
	
		/********************************************************	
		* Calculate length scales
		*********************************************************/
		if(theta_diff > 0){
		
			lh = dx;//fmin(   grid_cubicroot   ,   0.76 * sqrt( tke[i][j][k] ) / sqrt(bruntv)  );
			lvv = fmin(   5   ,   0.76 * sqrt( tke[i][j][k] ) / sqrt(bruntv)  );

		} else {
		//printf("here");
			lh = dx;//grid_cubicroot;
			lvv = 5;
		}
		//lvv = 4;

		if(theta_diff > 0){

			l = min( grid_cubicroot, 0.76 * sqrt( tke[i][j][k] ) / sqrt(bruntv));
		} else {

			l = grid_cubicroot;
		}

		C = 1.9*Ck + (0.93-1.9*Ck)*l / grid_cubicroot;
	
		/********************************************************	
		* Calculate eddy diffusivities
		*********************************************************/
		Kdiffh = Ck * lh * sqrt(tke[i][j][k]);
		Kdiffv = Ck * lvv * sqrt(tke[i][j][k]);//Kdiffh;
	
		Kdiffv = 0.0005;

		Kdiffhs = Kdiffh * 3;//(1 + 2*lh / grid_cubicroot);
		Kdiffvs = Kdiffv * (1 + 2*lvv / dz);
	
	
		/********************************************************	
		* Calculate rate of deformation tensors
		*********************************************************/
		D11 = (u[i+1][j][k]-u[i][j][k])/dx;
		D22 = (v[i][j+1][k]-v[i][j][k])/dy;
		D33 = (w[i][j][k+1]-w[i][j][k])/dz;
		/**************************************************************************************/ 
		D12 = .25*((u[i  ][j+1][k  ]+u[i+1][j+1][k  ] - u[i  ][j-1][k  ]-u[i+1][j-1][k  ])/dy +

	               (v[i+1][j  ][k  ]+v[i+1][j+1][k  ] - v[i-1][j  ][k  ]-v[i-1][j+1][k  ])/dx);
		/**************************************************************************************/
		D13 = .25*((u[i  ][j  ][k+1]+u[i+1][j  ][k+1] - u[i  ][j  ][k-1]-u[i+1][j  ][k-1])/dz +
	     			    
	           	   (w[i+1][j  ][k  ]+w[i+1][j  ][k+1] - w[i-1][j  ][k  ]-w[i-1][j  ][k+1])/dx);
		/**************************************************************************************/
		D23 = .25*((v[i+1][j  ][k  ]+v[i+1][j+1][k  ] - v[i-1][j  ][k  ]-v[i-1][j+1][k  ])/dz + // correct this!

      	           (w[i  ][j+1][k  ]+w[i  ][j+1][k+1] - w[i  ][j-1][k  ]-w[i  ][j-1][k+1])/dy);
		/**************************************************************************************/
	
		/********************************************************	
		* Shear production is a function of the
		* deformation of the velocity field
		*********************************************************/
		shear = Kdiffh*D11*D11 + Kdiffh*D22*D22 + Kdiffv*D33*D33 

			  + Kdiffh*D12*D12 + Kdiffv*D13*D13 + Kdiffv*D23*D23;	
	
		/********************************************************	
		* Buoyancy production is a function of the
		* Brunt-Vaisala frequency
		*********************************************************/

		buoyancy = -Kdiffv*bruntv;
	
		//if(buoyancy>0){buoyancy = 0;};

		if(lh != 0){
	
			dissipation = - C * pow(tke[i][j][k],1.5) / lh;

		} else {

			dissipation = 0;
		}
	//dissipation = 0;
		/********************************************************	
		* Prognostic TKE equation 
		*
		*********************************************************/

		tkep[i][j][k] = tkem[i][j][k] - step * dt * (
		/**************************************************************************************/
			0.5*((tke[i+1][j][k]+tke[i][j][k])*u[i+1][j][k] - (tke[i][j][k]+tke[i-1][j][k])* u[i][j][k])/dx
		/**************************************************************************************/
		+   0.5*((tke[i][j+1][k]+tke[i][j][k])*v[i][j+1][k] - (tke[i][j][k]+tke[i][j-1][k])* v[i][j][k])/dy
		/**************************************************************************************/
		+   0.5*((tke[i][j][k+1]+tke[i][j][k])*w[i][j][k+1] - (tke[i][j][k]+tke[i][j][k-1])* w[i][j][k])/dz
		/**************************************************************************************/
		-     (shear + buoyancy + dissipation)
		/**************************************************************************************/	
		);
		if(tkep[i][j][k] <= .1){ tkep[i][j][k] = .1;}

		if(i==20 && j==40 && k==5){ printf("Kdiff = %f tkem = %f tkep = %f tke = %f s = %f b = %f d = %f\n",Kdiffh,tkem[i][j][k],tkep[i][j][k],tke[i][j][k],shear,buoyancy,dissipation);}

	}}}

	printf("tke\n");
	for(int j=40;j<50;j++){

		printf("%f ",tke[20][j][5]);
		//printf("%e ",pi[i+1][49][NZ-1]-pi[i][49][NZ-1]);

	}
	printf("\n");

}

/********************************************************	
* apply diffusion to model equations
*
*********************************************************/
void turbulent_diffusion(double step){

	dx2 = dx*dx;
	dy2 = dy*dy;
	dz2 = dz*dz;

	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=1;k<NZ-1;k++){

		up[i][j][k] = up[i][j][k] + step * dt * (
						Kdiffv*((um[i][j][k+1] - 2.0*um[i][j][k]+um[i][j][k-1])/dz2) +
						Kdiffh*((um[i+1][j][k] - 2.0*um[i][j][k]+um[i-1][j][k])/dx2) +
						Kdiffh*((um[i][j+1][k] - 2.0*um[i][j][k]+um[i][j-1][k])/dy2)
						);
	
		vp[i][j][k] = vp[i][j][k] + step * dt * (
					  Kdiffv*((vm[i][j][k+1]-2.0*vm[i][j][k]+vm[i][j][k-1])/dz2) +
				      Kdiffh*((vm[i+1][j][k]-2.0*vm[i][j][k]+vm[i-1][j][k])/dx2) +
				      Kdiffh*((vm[i][j+1][k]-2.0*vm[i][j][k]+vm[i][j-1][k])/dy2)
					);
	
		      /*
		wp[i][j][k] = wp[i][j][k] + step * dt * (
					  Kdiffv*((wm[i][j][k+1]-2.0*wm[i][j][k]+wm[i][j][k-1])/dz2) +
				      Kdiffh*((wm[i+1][j][k]-2.0*wm[i][j][k]+wm[i-1][j][k])/dx2) +
				      Kdiffh*((wm[i][j+1][k]-2.0*wm[i][j][k]+wm[i][j-1][k])/dy2)
					);
	
	      
		thp[i][j][k] = thp[i][j][k] + step * dt * (
					(Kdiffvs*(thm[i][j][k+1]-2.0*thm[i][j][k]+thm[i][j][k-1])/dz2) +
					(Kdiffhs*(thm[i+1][j][k]-2.0*thm[i][j][k]+thm[i-1][j][k])/dx2)// +
					(Kdiffhs*(thm[i][j+1][k]-2.0*thm[i][j][k]+thm[i][j-1][k])/dy2) +
					//(Kdiffvs*(tb[k+1]-2.0*tb[k]+tb[k-1])/dz2)
					);
			*/		
		tkep[i][j][k] = tkep[i][j][k] + step * dt * (
					(Kdiffvs*(tkem[i][j][k+1]-2.0*tkem[i][j][k]+tkem[i][j][k-1])/dz2) +
					 (Kdiffhs*(tkem[i+1][j][k]-2.0*tkem[i][j][k]+tkem[i-1][j][k])/dx2) +
					 (Kdiffhs*(tkem[i][j+1][k]-2.0*tkem[i][j][k]+tkem[i][j-1][k])/dy2)
					);

		thp[i][j][k] = thp[i][j][k] + step*dt* (	
			+ Kdiffhs * ( (thm[i+1][j][k]-thb[i+1][j][k]) - 2.0*(thm[i][j][k]-thb[i][j][k]) + (thm[i-1][j][k]-thb[i-1][j][k] ))/(dx*dx) 
			+ Kdiffhs * ( (thm[i][j+1][k]-thb[i][j+1][k]) - 2.0*(thm[i][j][k]-thb[i][j][k]) + (thm[i][j-1][k]-thb[i][j-1][k] ))/(dy*dy) 
			+ Kdiffvs * ( (thm[i][j][k+1]-thb[i][j][k+1]) - 2.0*(thm[i][j][k]-thb[i][j][k]) + (thm[i][j][k-1]-thb[i][j][k-1] ))/(dz*dz))
			;

	}}}

}


/********************************************************	
* 
*
*********************************************************/
void initialize_tke(){

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		tke[i][j][k] = .005;
		tkem[i][j][k] = .005;
		tkep[i][j][k] = .005;
	}}}

}
