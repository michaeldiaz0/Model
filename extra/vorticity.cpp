#include "stdafx.h"
#include "vorticty.h"
#include "Heating.h"

#undef U
#undef V
#undef W
#undef TH
#undef PI

#define U(i,j,k)  u [i][j][k]
#define V(i,j,k)  v [i][j][k]
#define W(i,j,k)  w [i][j][k]
#define TH(i,j,k) th[i][j][k]
#define PI(i,j,k) pi[i][j][k]

#if STRETCHED_GRID
	#define MU(k) mu[k]
	#define MW(k) mw[k]
	#define ZU(k) zsu[k]
	#define ZW(k) zsw[k]
#else
	#define MU(k) 1.0
	#define MW(k) 1.0
	#define ZU(k) zu[k]
	#define ZW(k) zw[k]
#endif

/*********************************************************************
* Total vorticity
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
double get_vorticity(int xl, int xh, int yl, int yh, int zl, int zh){

	double eke;
	double umid, vmid;
	double eke_vavg, eke_havg;

	eke_havg = 0;

	for(int i=xl;i<xh;i++){
	for(int j=yl;j<yh;j++){

		eke_vavg = 0;

		for(int k=zl;k<zh;k++){
		
			umid = 0.5*(U(i,j,k)+U(i+1,j,k));
			vmid = 0.5*(V(i,j,k)+V(i,j+1,k));

			eke_vavg = eke_vavg + 0.5*rhou[k]*(umid*umid+vmid*vmid)*DZU(k);
		}

		eke_vavg = eke_vavg / (ZW(zh)-ZW(zl));

		eke_havg = eke_havg + eke_vavg;

	}}

	eke_havg = eke_havg / (double) ( (xh-xl)*(yh-yl) );

	return eke_havg;
}

/*********************************************************************
* Perform energy budget from a model output file.
**********************************************************************/
void vorticity_budget_from_file(const char * myfilename){

	get_model_data(myfilename,0);

	for(int i=0;i<100;i++){

		get_model_data(myfilename,i);
		
		//print_energy_budget(5,NX-5,5,NY-5,15,35);
	
		//print_energy_budget(4,NX-3,5,NY-5,1,20);
		print_energy_budget(30,120,10,50,1,28);
		mtime = mtime + 10800;
	}
}

/*********************************************************************
* Print an energy budget. 
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
void print_vorticity_budget(int xl, int xh, int yl, int yh, int zl, int zh, FILE *infile){

	double zeta;

	zeta = get_vorticity(xl,xh,yl,yh,zl,zh);

	if(infile==NULL){


	}
	
	
}

