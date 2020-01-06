#include "stdafx.h"



/*********************************************************************
*
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
		
		zu[i] = zh1;
		zw[i] = zh2;
		
		mu[i] = m1;
		mw[i] = m2;
	}
	
}


/*********************************************************************
*
*
**********************************************************************/
int main(){

	double zsu[NZ];
	double zsw[NZ];
	double mu[NZ];
	double mw[NZ];
	
	stretched_grid(&zsu[0],&zsw[0],&mu[0],&mw[0],50,1);

	for(int i=0;i<NZ;i++){
		
		printf("%d %.0f %f %f %f\n",i,zsu[i],zsw[i],mu[i],mw[i]);
		
	}
	
}
