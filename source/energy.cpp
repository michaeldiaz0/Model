#include "stdafx.h"
#include "energy.h"
#include "damping.h"
#include "Heating.h"



#if USE_MICROPHYSICS
	#define BUOYANCY(i,j,k) ( TH(i,j,k) / tbv[k] + 0.61*QV(i,j,k) - QC(i,j,k) - QR(i,j,k) )
#else
	#define BUOYANCY(i,j,k) ( TH(i,j,k) / tbv[k] )
#endif

/*********************************************************************
* Perform an eddy kinetic energy budget
**********************************************************************/

char energyBudgetFileName[len];

double myterms[6];
const double sec_per_day = 60.*60.*24.;

double get_residual(double,double,double,double,double,double,double);

double get_EKE(int,int,int,int,int,int);
double get_CKE(int,int,int,int,int,int,double* terms=NULL);
double get_PKE(int,int,int,int,int,int);
double get_FKE(int,int,int,int,int,int);
double get_pwork(int,int,int,int,int,int,int);
double get_gflux(int,int,int,int,int,int);
double get_advection(int,int,int,int,int,int);
double get_DKE(int xl, int xh, int yl, int yh, int zl, int zh);
double get_IDKE(int xl, int xh, int yl, int yh, int zl, int zh);
double get_temp_advection(int xl, int xh, int yl, int yh, int zl, int zh);

/*********************************************************************
*
*
**********************************************************************/
int process_command_line_args_energy(int argc, char *argv[]){
	
	strcpy(energyBudgetFileName,defaultEnergyBudgetFileName);
	
	//printf("%s\n",argv[2]);
	
	int c;
	
	while ((c = getopt (argc, argv, "e:s:")) != -1){
		switch (c){
			case 'e':
				if(VERBOSE){
					printf("Processing energy budget from file %s\n",optarg);
				}

				strcpy(energyBudgetFileName,optarg);

				//break;
			
			case 's':
				//printf("%s\n",optarg);
				break;
			case '?':
				if (optopt == 'e')
				  fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else if (isprint (optopt))
				  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
				  fprintf (stderr,
				           "Unknown option character `\\x%x'.\n",
				           optopt);
				return 1;
			default:
				abort ();
		}
	  }
	
	
	  return 1;
}

/********************************************************
*
* 
*
*********************************************************/
double CubicInterpolate3(double y0,double y1,double y2,double y3,double mu){
   
	double a,b,c,d,f0p,f1p;

	f0p = (y2-y0)/2;
	f1p = (y3-y1)/2;

	a = 2*y1 - 2*y2 +  f0p + f1p;
	b = -3*y1 + 3*y2 - 2*f0p - f1p;
	c = f0p;
	d = y1;
	
	return a*mu*mu*mu + b*mu*mu + c*mu + d;
}

/*********************************************************************
* Total eddy kinetic energy
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
double get_EKE(int xl, int xh, int yl, int yh, int zl, int zh){

	double umid, vmid, eke_havg;

	eke_havg = 0;

	for(int i=xl;i<xh;i++){
	for(int j=yl;j<yh;j++){
	for(int k=zl;k<zh;k++){
	
		umid = 0.5*(U(i,j,k)+U(i+1,j,k));
		vmid = 0.5*(V(i,j,k)+V(i,j+1,k));

		eke_havg += 0.5*rhou[k]*(umid*umid+vmid*vmid)*DZU(k);
		
	}}}

	eke_havg = eke_havg / (double) ( (xh-xl)*(yh-yl)*(ZW(zh)-ZW(zl)) );

	return eke_havg;
}

/*********************************************************************
* Barotropic terms (Reynold's stresses)
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
double get_CKE(int xl, int xh, int yl, int yh, int zl, int zh, double *terms){

	double umid, vmid, wmid;
	double cke_vavg, cke_havg;

	double uvdudy,uvdvdx,vvdvdy,uududx,uwdudz,vwdvdz;

	double outterms[6];

	cke_havg = 0;


	for(int i=0;i<6;i++){ outterms[i] = 0;}

	for(int i=xl;i<xh;i++){
	for(int j=yl;j<yh;j++){

		cke_vavg = 0;

		for(int k=zl;k<zh;k++){
		
			umid = 0.5*(U(i,j,k)+U(i+1,j,k));
			vmid = 0.5*(V(i,j,k)+V(i,j+1,k));
			wmid = 0.5*(W(i,j,k)+W(i,j,k+1));

			//---------------------------------------------
			// u * v * dU/dy
			//---------------------------------------------
			uvdudy =  0.25*(IUBAR(i+1,j+1,k)+IUBAR(i+1,j,k)+IUBAR(i,j+1,k)+IUBAR(i,j,k));
			uvdudy -= 0.25*(IUBAR(i+1,j-1,k)+IUBAR(i+1,j,k)+IUBAR(i,j-1,k)+IUBAR(i,j,k));

			uvdudy = -umid*vmid*uvdudy / dy;

			//---------------------------------------------
			// v * v * dV/dx
			//---------------------------------------------
			uvdvdx =  0.25*(IVBAR(i+1,j+1,k)+IVBAR(i+1,j,k)+IVBAR(i,j+1,k)+IVBAR(i,j,k));
			uvdvdx -= 0.25*(IVBAR(i-1,j+1,k)+IVBAR(i-1,j,k)+IVBAR(i,j+1,k)+IVBAR(i,j,k));

			uvdvdx = -umid*vmid*uvdvdx / dx;

			//---------------------------------------------
			// v * v * dV/dy
			//---------------------------------------------
			vvdvdy = -vmid*vmid*(IVBAR(i,j+1,k)-IVBAR(i,j,k)) / dy;

			//---------------------------------------------
			// u * u * dU/dx
			//---------------------------------------------
			uududx = -umid*umid*(IUBAR(i+1,j,k)-IUBAR(i,j,k)) / dx;

			//---------------------------------------------
			// u * w * dU/dz
			//---------------------------------------------
			uwdudz =  0.25*(IUBAR(i+1,j,k+1)+IUBAR(i,j,k+1)+IUBAR(i+1,j,k)+IUBAR(i,j,k)) * rhow[k+1];
			uwdudz -= 0.25*(IUBAR(i+1,j,k-1)+IUBAR(i,j,k-1)+IUBAR(i+1,j,k)+IUBAR(i,j,k)) * rhow[k];

			uwdudz = -one_d_rhou[k]*umid*wmid*uwdudz * ONE_D_DZ(k);

			//---------------------------------------------
			// v * w * dV/dz
			//---------------------------------------------
			vwdvdz =  0.25*(IVBAR(i,j+1,k+1)+IVBAR(i,j,k+1)+IVBAR(i,j+1,k)+IVBAR(i,j,k)) * rhow[k+1];
			vwdvdz -= 0.25*(IVBAR(i,j+1,k-1)+IVBAR(i,j,k-1)+IVBAR(i,j+1,k)+IVBAR(i,j,k)) * rhow[k];

			vwdvdz = -one_d_rhou[k]*vmid*wmid*vwdvdz * ONE_D_DZ(k);


			outterms[0] += rhou[k] * uvdudy * DZU(k);
			outterms[1] += rhou[k] * uvdvdx * DZU(k);
			outterms[2] += rhou[k] * uududx * DZU(k);
			outterms[3] += rhou[k] * vvdvdy * DZU(k);
			outterms[4] += rhou[k] * uwdudz * DZU(k);
			outterms[5] += rhou[k] * vwdvdz * DZU(k);

			cke_havg = cke_havg + rhou[k] * (uvdudy+uvdvdx+vvdvdy+uududx+uwdudz+vwdvdz) * DZU(k);

		}

		//cke_vavg = cke_vavg / (double)(zh-zl);

		//cke_havg = cke_havg + cke_vavg;

	}}

	cke_havg = cke_havg / (double) ( (xh-xl)*(yh-yl)*(ZW(zh)-ZW(zl)) );


	if(terms!=NULL){ 

		for(int i=0;i<6;i++){ terms[i] = outterms[i] / (double)( (xh-xl)*(yh-yl)*(ZW(zh)-ZW(zl)));}

	}

	return cke_havg;
}

/*********************************************************************
* Baroclinic energy production. 
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
double get_PKE(int xl, int xh, int yl, int yh, int zl, int zh){

	double wmid;
	double pke_vavg, pke_havg;

	pke_havg = 0;

	for(int i=xl;i<xh;i++){
	for(int j=yl;j<yh;j++){

		pke_vavg = 0;

		for(int k=zl;k<zh;k++){
		
			wmid = 0.5*(W(i,j,k)+W(i,j,k+1));

			pke_vavg = pke_vavg + rhou[k]*wmid*BUOYANCY(i,j,k) * DZU(k);	// could update with perturbation virtual temperature
		}

		pke_vavg = pke_vavg / (ZW(zh)-ZW(zl));// / (double)(zh-zl);

		pke_havg = pke_havg + pke_vavg;

	}}

	pke_havg = grav * pke_havg / (double) ( (xh-xl)*(yh-yl) );

	return pke_havg;
}

/*********************************************************************
* Frictional dissipation.
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
double get_FKE(int xl, int xh, int yl, int yh, int zl, int zh){

	double umid, vmid;
	double fke_vavg, fke_havg;

	fke_havg = 0;

	for(int i=xl;i<xh;i++){
	for(int j=yl;j<yh;j++){

		fke_vavg = 0;

		for(int k=zl;k<zh;k++){
		
			if(!USE_TURBULENT_STRESS){
				umid = 0.5*(U(i,j,k)+U(i+1,j,k));
				vmid = 0.5*(V(i,j,k)+V(i,j+1,k));
			
				fke_vavg = fke_vavg + rhou[k]*( IFRICTION(i,j,k)*umid*umid + IFRICTION(i,j,k)*vmid*vmid )*DZU(k);
				/*
				fke_vavg = fke_vavg + rhou[k]*( 0.5*(IFRICTION(i,j,k)*U(i,j,k)*U(i,j,k) + IFRICTION(i+1,j,k)*U(i+1,j,k)*U(i+1,j,k)) 
											  + 0.5*(IFRICTION(i,j,k)*V(i,j,k)*V(i,j,k) + IFRICTION(i,j+1,k)*V(i,j+1,k)*V(i,j+1,k) ) 
												  )*DZU(k);
				*/
			} else {
				fke_vavg = fke_vavg + rhou[k]*IFRICTION(i,j,k)*DZU(k);
			}
		}

		fke_havg = fke_havg + fke_vavg;

	}}

	fke_havg = fke_havg / (double) ( (xh-xl)*(yh-yl)*(ZW(zh)-ZW(zl)) );

	return -fke_havg;
}

/*********************************************************************
* Explicit diffusion.
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
double get_DKE(int xl, int xh, int yl, int yh, int zl, int zh){

	double dke = 0;
	double umid, vmid,udiffmid,vdiffmid;

	for(int i=xl;i<xh;i++){
	for(int j=yl;j<yh;j++){
	for(int k=zl;k<zh;k++){
		
		umid = 0.5*(U(i,j,k)+U(i+1,j,k));
		vmid = 0.5*(V(i,j,k)+V(i,j+1,k));
		
		udiffmid = 0.5*(u_diff_tend[INDEX(i,j,k)]+u_diff_tend[INDEX(i+1,j,k)]);
		vdiffmid = 0.5*(v_diff_tend[INDEX(i,j,k)]+v_diff_tend[INDEX(i,j+1,k)]);
		
		dke += rhou[k]* ( umid * udiffmid + vmid * vdiffmid ) * DZU(k);
	}}}

	return dke / (double) ( (xh-xl)*(yh-yl)*(ZW(zh)-ZW(zl)) );
}

/*********************************************************************
* Implicit diffusion
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
double get_IDKE(int xl, int xh, int yl, int yh, int zl, int zh){

	double dke = 0;
	double umid, vmid,wmid,udiff_x,vdiff_x,udiff_y,vdiff_y,udiff_z,vdiff_z,udiff,vdiff;
	double ubar_mid,vbar_mid,wbar_mid;
	double one_d_60 = 1.0/60.0;
	double one_d_12 = 1.0/12.0;
	
	double u_at_vpoint,v_at_upoint,w_at_vpoint,w_at_upoint;
	
	for(int i=xl;i<xh;i++){
	for(int j=yl;j<yh;j++){
	for(int k=zl;k<zh;k++){
		
		umid = 0.5*(U(i,j,k)+U(i+1,j,k));
		vmid = 0.5*(V(i,j,k)+V(i,j+1,k));
		wmid = 0.5*(W(i,j,k)+W(i,j,k+1));
		
		ubar_mid = 0.5*( IUBAR(i,j,k) + IUBAR(i+1,j,k) );
		vbar_mid = 0.5*( IVBAR(i,j,k) + IVBAR(i,j+1,k) );
		wbar_mid = 0.5*( IWBAR(i,j,k) + IWBAR(i,j,k+1) );
		
		udiff_x = 0.5*(diffuse_i_6th(us,i,j,k) + diffuse_i_6th(us,i+1,j,k));
		udiff_y = 0.5*(diffuse_j_6th(us,i,j,k) + diffuse_j_6th(us,i+1,j,k));
		
		vdiff_x = 0.5*(diffuse_i_6th(vs,i,j,k) + diffuse_i_6th(vs,i,j+1,k));
		vdiff_y = 0.5*(diffuse_j_6th(vs,i,j,k) + diffuse_j_6th(vs,i,j+1,k));
		
		if(k>1 && k<NZ-2){
			udiff_z = -0.5*(diffuse_k_4th(us,i,j,k) + diffuse_k_4th(us,i+1,j,k));
			vdiff_z = -0.5*(diffuse_k_4th(vs,i,j,k) + diffuse_k_4th(vs,i,j+1,k));
		} else {
			udiff_z = 0;
			vdiff_z = 0;
		}

		u_at_vpoint = 0.25 * (U(i,j,k)+U(i+1,j,k)+U(i,j-1,k)+U(i+1,j-1,k)+IUBAR(i,j,k)+IUBAR(i+1,j,k)+IUBAR(i,j-1,k)+IUBAR(i+1,j-1,k));
		v_at_upoint = 0.25 * (V(i,j,k)+V(i,j+1,k)+V(i-1,j,k)+V(i-1,j+1,k)+IVBAR(i,j,k)+IVBAR(i,j+1,k)+IVBAR(i-1,j,k)+IVBAR(i-1,j+1,k));
		
		w_at_vpoint = 0.25*(W(i,j,k+1)+W(i,j-1,k+1)+W(i,j,k)+W(i,j-1,k)+IWBAR(i,j,k+1)+IWBAR(i,j-1,k+1)+IWBAR(i,j,k)+IWBAR(i,j-1,k));
		w_at_upoint = 0.25*(W(i,j,k+1)+W(i-1,j,k+1)+W(i,j,k)+W(i-1,j,k)+IWBAR(i,j,k+1)+IWBAR(i-1,j,k+1)+IWBAR(i,j,k)+IWBAR(i-1,j,k));
		
		udiff = one_d_60*fabs(ubar_mid+umid)*one_d_dx    * udiff_x +
				one_d_60*fabs( v_at_upoint )*one_d_dy	 * udiff_y +
				one_d_12*fabs( w_at_upoint )*ONE_D_DZ(k) * udiff_z;
		
		vdiff = one_d_60*fabs( u_at_vpoint )*one_d_dx    * vdiff_x +
				one_d_60*fabs(vbar_mid+vmid)*one_d_dy    * vdiff_y +
				one_d_12*fabs( w_at_vpoint )*ONE_D_DZ(k) * vdiff_z;
		
		dke += rhou[k]* ( umid * udiff + vmid * vdiff ) * DZU(k);
	}}}

	return dke / (double) ( (xh-xl)*(yh-yl)*(ZW(zh)-ZW(zl)) );
}

/*********************************************************************
* Implicit diffusion
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
double get_IDKE_basic(int xl, int xh, int yl, int yh, int zl, int zh){

	double dke = 0;
	double umid, vmid,wmid,udiff_x,vdiff_x,udiff_y,vdiff_y,udiff_z,vdiff_z,udiff,vdiff;
	double one_d_60 = 1.0/60.0;
	double one_d_12 = 1.0/12.0;
	
	double u_at_vpoint,v_at_upoint,w_at_vpoint,w_at_upoint;
	
	for(int i=xl;i<xh;i++){
	for(int j=yl;j<yh;j++){
	for(int k=zl;k<zh;k++){
		
		umid = 0.5*(U(i,j,k)+U(i+1,j,k));
		vmid = 0.5*(V(i,j,k)+V(i,j+1,k));
		wmid = 0.5*(W(i,j,k)+W(i,j,k+1));
			
		udiff_x = 0.5*(diffuse_i_6th(iubar,i,j,k) + diffuse_i_6th(iubar,i+1,j,k));
		udiff_y = 0.5*(diffuse_j_6th(iubar,i,j,k) + diffuse_j_6th(iubar,i+1,j,k));
		
		vdiff_x = 0.5*(diffuse_i_6th(ivbar,i,j,k) + diffuse_i_6th(ivbar,i,j+1,k));
		vdiff_y = 0.5*(diffuse_j_6th(ivbar,i,j,k) + diffuse_j_6th(ivbar,i,j+1,k));
		
		if(k>1 && k<NZ-2){
			udiff_z = -0.5*(diffuse_k_4th(iubar,i,j,k) + diffuse_k_4th(iubar,i+1,j,k));
			vdiff_z = -0.5*(diffuse_k_4th(ivbar,i,j,k) + diffuse_k_4th(ivbar,i,j+1,k));
		} else {
			udiff_z = 0;
			vdiff_z = 0;
		}

		u_at_vpoint = 0.25 * (U(i,j,k)+U(i+1,j,k)+U(i,j-1,k)+U(i+1,j-1,k));
		v_at_upoint = 0.25 * (V(i,j,k)+V(i,j+1,k)+V(i-1,j,k)+V(i-1,j+1,k));
		
		w_at_vpoint = 0.25*(W(i,j,k+1)+W(i,j-1,k+1)+W(i,j,k)+W(i,j-1,k));
		w_at_upoint = 0.25*(W(i,j,k+1)+W(i-1,j,k+1)+W(i,j,k)+W(i-1,j,k));
		
		udiff = one_d_60*fabs( umid		   )*one_d_dx    * udiff_x +
				one_d_60*fabs( v_at_upoint )*one_d_dy	 * udiff_y +
				one_d_12*fabs( w_at_upoint )*ONE_D_DZ(k) * udiff_z;
		
		vdiff = one_d_60*fabs( u_at_vpoint )*one_d_dx    * vdiff_x +
				one_d_60*fabs( vmid        )*one_d_dy    * vdiff_y +
				one_d_12*fabs( w_at_vpoint )*ONE_D_DZ(k) * vdiff_z;
		
		dke += rhou[k]* ( umid * udiff + vmid * vdiff ) * DZU(k);
	}}}

	return dke / (double) ( (xh-xl)*(yh-yl)*(ZW(zh)-ZW(zl)) );
}

/*********************************************************************
* Pressure flux. 
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
double get_gflux(int xl, int xh, int yl, int yh, int zl, int zh){

	double uflux, vflux, wflux;
	double gflux_vavg, gflux_havg;

	gflux_havg = 0;

	for(int i=xl;i<xh;i++){
	for(int j=yl;j<yh;j++){

		gflux_vavg = 0;

		for(int k=zl;k<zh;k++){
		
			uflux = 0.5 *  ( U(i+1,j,k)*(PI(i,j,k)+PI(i+1,j,k)) - U(i,j,k)*(PI(i,j,k)+PI(i-1,j,k)))/dx;
			vflux = 0.5 *  ( V(i,j+1,k)*(PI(i,j,k)+PI(i,j+1,k)) - V(i,j,k)*(PI(i,j,k)+PI(i,j-1,k)))/dy;
			wflux = 0.5 *  ( W(i,j,k+1)*(PI(i,j,k)+PI(i,j,k+1)) - W(i,j,k)*(PI(i,j,k)+PI(i,j,k-1)))/dz;

			gflux_vavg = gflux_vavg + rhou[k]*cp*tb[k]*(uflux+vflux+wflux);

		}

		gflux_vavg = gflux_vavg / (double)(zh-zl);

		gflux_havg = gflux_havg + gflux_vavg;

	}}

	gflux_havg = gflux_havg / (double) ( (xh-xl)*(yh-yl) );

	return -gflux_havg;
}

/*********************************************************************
* Pressure flux. 
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
double get_temp_advection(int xl, int xh, int yl, int yh, int zl, int zh){

	double uad, vad, dtbdz;

	double total = 0;

	for(int i=xl;i<xh;i++){
	for(int j=yl;j<yh;j++){
	for(int k=zl;k<zh;k++){
		
			dtbdz = 0.5 * (tb[k+1] - tb[k-1]) * ONE_D_DZ(k);
		
			uad = -0.5 * ( U(i+1,j,k)+U(i,j,k) ) * ( ITHBAR(i+1,j,k)-ITHBAR(i,j,k) ) / dx;
			vad = -0.5 * ( V(i,j+1,k)+V(i,j,k) ) * ( ITHBAR(i,j+1,k)-ITHBAR(i,j,k) ) / dy;

			total +=  grav * TH(i,j,k) / tb[k] * (1.0/dtbdz) * (uad + vad) * rhou[k] * DZU(k);

	}}}

	total = total / (double) ( (xh-xl)*(yh-yl)*(ZW(zh)-ZW(zl)) );

	return total;
}

/*********************************************************************
* Pressure work. 
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
double get_pwork(int xl, int xh, int yl, int yh, int zl, int zh,int interp){

    double uflux, vflux;
    double pwork_havg;
    double umid, vmid;

    pwork_havg = 0;

	switch(interp){
		//--------------------------------------------------------------------------------
		// Cubic interpolate winds, 4th order pressure
		//--------------------------------------------------------------------------------
		case 0:
		    for(int i=xl;i<xh;i++){
		    for(int j=yl;j<yh;j++){
		    for(int k=zl;k<zh;k++){
				
	            umid = CubicInterpolate3(U(i-1,j,k),U(i,j,k),U(i+1,j,k),U(i+2,j,k),0.5);
	            vmid = CubicInterpolate3(V(i,j-1,k),V(i,j,k),V(i,j+1,k),V(i,j+2,k),0.5);

	            uflux = umid * (1.0*PI(i-2,j,k)-8.0*PI(i-1,j,k)+8.0*PI(i+1,j,k)-1.0*PI(i+2,j,k)) / (12.0*dx);
            
	        	vflux = vmid * (1.0*PI(i,j-2,k)-8.0*PI(i,j-1,k)+8.0*PI(i,j+1,k)-1.0*PI(i,j+2,k)) / (12.0*dy);
				
				pwork_havg += rhou[k]*(uflux+vflux)*DZU(k);
			}}}				
		break;
		//--------------------------------------------------------------------------------
		// Approximate energy in volume by actual values at south and west sides
		//--------------------------------------------------------------------------------
		case 1:		
		    for(int i=xl;i<xh;i++){
		    for(int j=yl;j<yh;j++){
		    for(int k=zl;k<zh;k++){
				
				uflux = U(i,j,k) * (PI(i,j,k)-PI(i-1,j,k)) / dx;
				vflux = V(i,j,k) * (PI(i,j,k)-PI(i,j-1,k)) / dy;
				
				pwork_havg += rhou[k]*(uflux+vflux)*DZU(k);
			}}}
		break;
		//--------------------------------------------------------------------------------
		// Linear interpolate winds to center, 2 delta centered difference for pressure
		//--------------------------------------------------------------------------------
		case 2:		
		    for(int i=xl;i<xh;i++){
		    for(int j=yl;j<yh;j++){
		    for(int k=zl;k<zh;k++){
				
				umid = 0.5*(U(i,j,k)+U(i+1,j,k));
				vmid = 0.5*(V(i,j,k)+V(i,j+1,k));
				
				uflux = 0.5 * umid * (PI(i+1,j,k)-PI(i-1,j,k)) * one_d_dx;
				vflux = 0.5 * vmid * (PI(i,j+1,k)-PI(i,j-1,k)) * one_d_dy;
				
				pwork_havg += rhou[k]*(uflux+vflux)*DZU(k);
			}}}	
		break;
		//--------------------------------------------------------------------------------
		// Approximate energy in volume by linear interpolation of actual values from all four sides
		//--------------------------------------------------------------------------------
		case 3:		
		    for(int i=xl;i<xh;i++){
		    for(int j=yl;j<yh;j++){
		    for(int k=zl;k<zh;k++){
			
				uflux = 0.5*( U(i,j,k) * (PI(i,j,k)-PI(i-1,j,k)) + U(i+1,j,k) * (PI(i+1,j,k)-PI(i,j,k))) * one_d_dx;
				vflux = 0.5*( V(i,j,k) * (PI(i,j,k)-PI(i,j-1,k)) + V(i,j+1,k) * (PI(i,j+1,k)-PI(i,j,k))) * one_d_dy;
							
				pwork_havg += rhou[k]*(uflux+vflux)*DZU(k);
			}}}			
		break;
	}


	pwork_havg = pwork_havg / (double) ( (xh-xl)*(yh-yl)*(ZW(zh)-ZW(zl)) );

	return -pwork_havg;
}

/*********************************************************************
* Advection. 
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
double get_advection(int xl, int xh, int yl, int yh, int zl, int zh){

	double umid=0,vmid=0;
	double eke1=0,eke2 = 0;
	double ns_flux=0,ew_flux=0,tb_flux = 0;
	double advection = 0;

	// northern and southern boundary fluxes
	for(int i=xl;i<xh;i++){
	for(int k=zl;k<zh;k++){

		// northern
		umid = 0.25*(U(i,yh,k)+U(i,yh-1,k)+U(i+1,yh,k)+U(i+1,yh-1,k));

		eke2 = 0.5*(umid*umid+V(i,yh,k)*V(i,yh,k));

		// southern
		umid = 0.25*(U(i,yl,k)+U(i,yl-1,k)+U(i+1,yl,k)+U(i+1,yl-1,k));

		eke1 = 0.5*(umid*umid+V(i,yl,k)*V(i,yl,k));

		// flux divergence
		#if !ISLINEAR
		ns_flux = ns_flux + rhou[k]*( (IVBAR(i,yh,k)+V(i,yh,k))*eke2 - (IVBAR(i,yl,k)+V(i,yl,k))*eke1 )*dx*DZU(k);
		#else
		ns_flux = ns_flux + rhou[k]*( (IVBAR(i,yh,k))*eke2 - (IVBAR(i,yl,k))*eke1 )*dx*DZU(k);	
		#endif

	}}

	// eastern and western boundary fluxes
	for(int j=yl;j<yh;j++){
	for(int k=zl;k<zh;k++){

		// eastern
		vmid = 0.25*(V(xh,j,k)+V(xh-1,j,k)+V(xh,j+1,k)+V(xh-1,j+1,k));

		eke2 = 0.5*(U(xh,j,k)*U(xh,j,k)+vmid*vmid);

		// western
		vmid = 0.25*(V(xl,j,k)+V(xl-1,j,k)+V(xl,j+1,k)+V(xl-1,j+1,k));

		eke1 = 0.5*(U(xl,j,k)*U(xl,j,k)+vmid*vmid);

		// flux divergence
		#if !ISLINEAR
		ew_flux = ew_flux + rhou[k]*( (IUBAR(xh,j,k)+U(xh,j,k))*eke2 - (IUBAR(xl,j,k)+U(xl,j,k))*eke1 )*dy*DZU(k);
		#else
		ew_flux = ew_flux + rhou[k]*( (IUBAR(xh,j,k))*eke2 - (IUBAR(xl,j,k))*eke1 )*dy*DZU(k);
		#endif
	}}

	// upper and lower boundary fluxes
	for(int i=xl;i<xh;i++){
	for(int j=yl;j<yh;j++){

		// upper
		umid = 0.25*(U(i,j,zh)+U(i,j,zh-1)+U(i+1,j  ,zh)+U(i+1,j  ,zh-1));
		vmid = 0.25*(V(i,j,zh)+V(i,j,zh-1)+V(i  ,j+1,zh)+V(i  ,j+1,zh-1));

		eke2 = 0.5*(umid*umid+vmid*vmid);

		// lower
		umid = 0.25*(U(i,j,zl)+U(i,j,zl-1)+U(i+1,j  ,zl)+U(i+1,j  ,zl-1));
		vmid = 0.25*(V(i,j,zl)+V(i,j,zl-1)+V(i  ,j+1,zl)+V(i  ,j+1,zl-1));

		eke1 = 0.5*(umid*umid+vmid*vmid);

		// flux divergence
		#if !ISLINEAR
		tb_flux = tb_flux + ( rhow[zh]*(IWBAR(i,j,zh)+W(i,j,zh))*eke2 - rhow[zl]*(IWBAR(i,j,zl)+W(i,j,zl))*eke1 )*dx*dy;
		#else
		tb_flux = tb_flux + ( rhow[zh]*(IWBAR(i,j,zh))*eke2 - rhow[zl]*(IWBAR(i,j,zl))*eke1 )*dx*dy;
		#endif
	}}

	advection = (-ns_flux - ew_flux - tb_flux) / (double) ( dx*dy*(xh-xl)*(yh-yl)*(ZW(zh)-ZW(zl)) );

	return advection;
}



/*********************************************************************
* Residual 
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
double get_residual(double eke1,double eke2,double dt,
					double shear,double advection,
					double pressure,double friction,
					double diffusion,double idiffusion)
{
	
	double dEKEdt = (eke2 - eke1) * 24.0 / dt;
	
	double RHS = shear + advection + pressure + friction + diffusion + idiffusion;
	
	//printf("%f %f\n",dEKEdt,RHS);
	
	return dEKEdt - RHS;
}

/*********************************************************************
* Perform energy budget from a model output file.
**********************************************************************/
void energy_budget_from_file(const char * myfilename){

	//get_model_data(myfilename,0);

	//Heating heat;
	//heat.initialize(21.18,86.3,19.37,93.0,100000.,6.0);
	//heat.printInfo();

    int start_index = 0;
	int end_index = (int)get_file_length(myfilename,"t");

	double * times = get_data2(myfilename,"time",end_index);
    
	double dhours = (times[start_index+1]-times[start_index]);
	
	mtime = (double)start_index*60.0*60.0*dhours;

	for(int i=start_index;i<end_index;i++){

		get_model_data(myfilename,i);

		print_energy_budget(5,NX-5,5,NY-5,1,35);

		mtime = mtime + dhours*60.0*60.0;
	}
}

/*********************************************************************
* Print an energy budget. 
*
* @param xh,xh - eastern and western boundary
* @param yh,yl - northern and southern boundary
* @param zh,zl - upper and lower boundary
**********************************************************************/
void print_energy_budget(int xl, int xh, int yl, int yh, int zl, int zh, FILE *infile){

	double eke=0,cke=0,pke=0,fke=0,pwork=0,advection=0,diffusion = 0,idiffusion=0,residual=0;

	static double lastTime,dHour;
	static double ekeMinusOne=-1,ekeMinusTwo=-1;
	static double shear2=0,advection2=0,pressure2=0,friction2=0,diff2=0,idiff2=0;

	eke = get_EKE(xl,xh,yl,yh,zl,zh);
	cke = get_CKE(xl,xh,yl,yh,zl,zh,&myterms[0])*sec_per_day;
	pke = get_PKE(xl,xh,yl,yh,zl,zh)*sec_per_day;
	//pke = get_temp_advection(xl,xh,yl,yh,zl,zh)*sec_per_day;
	fke = get_FKE(xl,xh,yl,yh,zl,zh)*sec_per_day;
	pwork = get_pwork(xl,xh,yl,yh,zl,zh,2)*sec_per_day;
	advection = get_advection(xl,xh,yl,yh,zl,zh)*sec_per_day;
	idiffusion = get_IDKE(xl,xh,yl,yh,zl,zh)*sec_per_day + get_IDKE_basic(xl,xh,yl,yh,zl,zh)*sec_per_day;
	
	
	if(OUTPUT_DIFFUSION_TEND){ diffusion = get_DKE(xl,xh,yl,yh,zl,zh)*sec_per_day; }
	
	if(infile==NULL){ dHour = (mtime - lastTime) / 3600; }
	else { dHour = 2.0;}
	
	//printf("%f %f %f\n",mtime,lastTime,dHour);
	if(ekeMinusOne != -1){
		
		residual = get_residual(ekeMinusOne,eke,dHour,
		0.5*(shear2+cke),
		0.5*(advection2+advection),
		0.5*(pressure2+pwork),
		0.5*(friction2+fke),
		0.5*(diff2+diffusion),
		0.5*(idiff2+idiffusion));
	} else {
		residual = -1;
	}

	lastTime = mtime;
	shear2 = cke;
	advection2 = advection;
	friction2 = fke;
	pressure2 = pwork;
	diff2 = diffusion;
	idiff2 = idiffusion;

	ekeMinusTwo = ekeMinusOne;
	ekeMinusOne = eke;

	for(int i=0;i<6;i++){ myterms[i] = myterms[i]*sec_per_day;}

	if(infile==NULL){

		printf("%.2f %f %f %f %f %f ",mtime/3600.f,eke,cke,pke,fke,pwork);

		printf("%f %f %f %f %f %f %f %f %f\n",myterms[0],myterms[1],myterms[2],myterms[3],myterms[4],myterms[5],advection,idiffusion,residual);

	} else {
	
		fprintf(infile,"%.1f %f %f %f %f %f ",mtime/3600.f,eke,cke,pke,fke,pwork);

		fprintf(infile,"%f %f %f %f %f %f %f %f\n",myterms[0],myterms[1],myterms[2],myterms[3],myterms[4],myterms[5],advection,residual);
		
		fflush(stdout);
	}
}
