
#include "stdafx.h"

double qv[NX][NY][NZ],qvp[NX][NY][NZ],qvm[NX][NY][NZ];
double qc[NX][NY][NZ],qcp[NX][NY][NZ],qcm[NX][NY][NZ];
double qr[NX][NY][NZ],qrp[NX][NY][NZ],qrm[NX][NY][NZ];
double vt[NX][NY][NZ],st[NX][NY][NZ];
double acc[NX][NY][NZ];

double qci[NX][NY][NZ],qcip[NX][NY][NZ], qcim[NX][NY][NZ];
double qs[NX][NY][NZ],qsp[NX][NY][NZ],qsm[NX][NY][NZ];

//double qvb[NX][NY][NZ];

double A,C,E,cvent,vtden,qrr,pd,theta,pressure,phi,vapor,total_convert,total_convert2;
double qvs,esw; // saturation mixing ratio / vapor pressure with respect to water
double qis,esi; // saturation mixing ratio / vapor pressure with respect to ice

double var;

/****************************************************
* Liquid accretion coefficients/variables
*****************************************************/
double N0r = 8e6;	// intercept value in raindrop size distribution
double ERC = 1;		// rain/cloud water collection efficiency
double lambdaR;		// slope of rain drop size distribution
double B;			// accretion

/****************************************************
* Ice nucleation coefficients/variables
*****************************************************/
double n0 = 1e-2;	// constant in expression for ice crystal concentration
double beta = 0.6;	// constant in ice crystal concentration
double M0 = 1e-12;	// initial mass of cloud ice crystal
double pint;		// initiation of cloud ice
double nc;			// number concentration of ice crystals 

/****************************************************
* Ice deposition coefficients/variables
*****************************************************/
double Rstar = 8.314e3;		// universal gas constant
double Mw = 18.016;			// molecular weight of water
double Ka = 2.43e-2;		// thermal conductivity of air
double chi = 2.26e-5;		// diffusivity of water vapor in air
double Ml;					// average mass of cloud ice crystal
double Dl;					// diamter of ice particles
double App;					// thermodynamic term in deposition
double Bpp;					// thermodynamic term in deposition
double pdepi;				// depositional growth

/****************************************************
* Snow autoconversion coefficients/variables
*****************************************************/
double Mmax = 9.4e-10;		// maximum mass of cloud ice crystal
double pconv;				// autoconversion rate
double qimax;				// conversion of cloud ice to snow threshhold

/****************************************************
* Snow coefficients/variables
*****************************************************/
double app = 1.139;		// fall speed coefficient
double b = 0.11;		// fall speed exponent
double lambdaS;			// slope factor
double N0S = 8e6;		// intercept value in snowflake size distribution
double rhoS = 200;		// density of snow (kg/m3)

/****************************************************
* Snow collecting cloud ice coefficients/variables
*****************************************************/
double ESI = 0.1;	// snow/cloud ice collection efficiency
double psaci;

/************************************************************
* Snow collecting cloud water coefficients/variables
*************************************************************/
double ESC = 1.f;	// snow/cloud water collection efficiency
double psacw;
double psacw_warm;	// source term for rain when above freezing
double psacw_cold;	// source term for snow when below freezing

/****************************************************
* Melting of cloud ice coefficients/variables
*****************************************************/
double psmlti;

/****************************************************
* Melting of snow coefficients/variables
*****************************************************/
double mu = 1.718e-5;
double psmlt;

/****************************************************
* Despositional growth of snow coefficients/variables
*****************************************************/
double psdep;

/****************************************************
* Evaporation of snow coefficients/variables
*****************************************************/
double pmltev;


/***********************************************************************************
* 
************************************************************************************/
void init_microphysics(){

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		qv[i][j][k] = 0; qvp[i][j][k] = 0; qvm[i][j][k] = 0;
		qc[i][j][k] = 0; qcp[i][j][k] = 0; qcm[i][j][k] = 0;
		qr[i][j][k] = 0; qrp[i][j][k] = 0; qrm[i][j][k] = 0;
		qci[i][j][k] = 0; qcip[i][j][k] = 0; qcim[i][j][k] = 0;
		qs[i][j][k] = 0; qsp[i][j][k] = 0; qsm[i][j][k] = 0;
	}}}
}

/***********************************************************************************
* Advect moisture variables
************************************************************************************/
void advect_microphysics(double step){

	double dtx = dt/dx;
	double dty = dt/dy;
	double dtz = dt/dz;

	advect_scalar(qvm,qv,qvp,qb,step);
	advect_scalar(qcm,qc,qcp,step);
	advect_scalar(qcim,qci,qcip,step);
	//advect_scalar(qsm,qs,qsp,step);

	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=1;k<NZ-1;k++){

			
		qrp[i][j][k] = qrm[i][j][k] + step * (

			- .50*dtx * (  u[i+1][j][k]*(qr[i+1][j][k]+qr[i][j][k]) - u[i][j][k]*(qr[i][j][k]+qr[i-1][j][k])  )

			- .50*dty * (  v[i][j+1][k]*(qr[i][j+1][k]+qr[i][j][k]) - v[i][j][k]*(qr[i][j][k]+qr[i][j-1][k])  )

			- .50*dtz * (  rhow[k+1]*(w[i][j][k+1]-(vt[i][j][k+1]+vt[i][j][k])/2)*(qr[i][j][k+1]+qr[i][j][k]) -

			rhow[k]*(w[i][j][k]-(vt[i][j][k]+vt[i][j][k-1])/2)*(qr[i][j][k]+qr[i][j][k-1]))/rhou[k]

		)
		;

		qsp[i][j][k] = qsm[i][j][k] + step * (

			- .50*dtx * (  u[i+1][j][k]*(qs[i+1][j][k]+qs[i][j][k]) - u[i][j][k]*(qs[i][j][k]+qs[i-1][j][k])  )

			- .50*dty * (  v[i][j+1][k]*(qs[i][j+1][k]+qs[i][j][k]) - v[i][j][k]*(qs[i][j][k]+qs[i][j-1][k])  )

			- .50*dtz * (  rhow[k+1]*(w[i][j][k+1]-(st[i][j][k+1]+st[i][j][k])/2)*(qs[i][j][k+1]+qs[i][j][k]) -

			rhow[k]*(w[i][j][k]-(st[i][j][k]+st[i][j][k-1])/2)*(qs[i][j][k]+qs[i][j][k-1]))/rhou[k]

		)
		;
			
		if(qcp[i][j][k] < 0)
			qcp[i][j][k] = 0;

		if(qvp[i][j][k]+qb[k] < 0)
			qvp[i][j][k] = -qb[k];

		if(qrp[i][j][k] < 0)
			qrp[i][j][k] = 0;

		if(qcip[i][j][k] < 0)
			qcip[i][j][k] = 0;

		if(qsp[i][j][k] < 0)
			qsp[i][j][k] = 0;

	}}}

}

/*******************************************************************************
* Perform microphysics calculations
* Handles all of the thermodynamics and phase changes
* for each grid cell
********************************************************************************/
void kessler_microphysics(){

	double cpRd = cp/Rd;
	double temperature;

	double a0=-0.267,a1=5.15e3,a2=-1.0225e6,a3=7.55e7;	// coefficients in polynomial fallspeed for rain

	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=1;k<NZ-1;k++){
	
		/********************************
		* Get actual values
		*********************************/
		theta = thp[i][j][k]+tb[k];			// actual potential temperature is base state plus perturbation
 		pressure = pi[i][j][k]+pib[k];		// actual pressure is base state plus perturbation
 		vapor = qvp[i][j][k]+qb[k];			// actual vapor is base state plus perturbation
		temperature = theta*pressure;		// temperature
		pd = p0*pow(pressure,cpRd);			// dimensional pressure

		// calculate saturation mixing ratio

		var = exp(21.87*(temperature-276.16)/ (temperature-7.66)  );
	 	qis = (380.f/pd) * var;
		esi = 6.1078 * var;
		qvs = (380.f/pd) * exp(17.27*(temperature-273.16)/ (temperature-35.86)  );
		
		if(qsp[i][j][k] != 0){ lambdaS = pow( (trigpi*rhoS*N0S) / (rhou[k]*qsp[i][j][k]),0.25); printf("%f\n",lambdaS);}
		else { lambdaS = 0;}

		/********************************
		* Autoconversion
		*********************************/
		if(qcp[i][j][k] > 1.0e-3){ A = 0.001*(qcp[i][j][k]-1.0e-3)*dt;}
		else { A = 0;}
		//A = 0;
		/****************************************************
		* Accretion - collection of cloud water by rainwater
		*****************************************************/
		if(qrp[i][j][k] > 1e-10){

			lambdaR = pow( (trigpi*1000*N0r) / (rhou[k]*qrp[i][j][k]),0.25);

			//B = 0.25*trigpi*rhou[k]*qc[i][j][k]*ERC*N0r*pow(p0/pressure,0.4) *
				//2/pow(lambdaR,3) * (a0 + a1*3/lambdaR + a2*12/pow(lambdaR,2) + a3*60/pow(lambdaR,3)) * dt;

			B = 0.25*trigpi*qc[i][j][k]*ERC*N0r*pow(p0/pd,0.4) *
				2/pow(lambdaR,3) * (a0 + a1*3/lambdaR + a2*12/pow(lambdaR,2) + a3*60/pow(lambdaR,3)) * dt;


			//B = rhou[k] * 2.2 * qcp[i][j][k] * pow(qrp[i][j][k], 0.875) * dt;

		} else { B = 0;}

		/********************************
		* Evaporation of rain water
		*********************************/
		if(qrp[i][j][k] > 1e-10){

			cvent = 1.6 + 30.39 * pow( rhou[k]*qrp[i][j][k], 0.2046);
			E = (1.f/rhou[k]) * (((1.f-vapor/qvs) * cvent * pow(rhou[k]*qrp[i][j][k],0.525)) /
				(2.03e4 + 9.58e6 / (pd*qvs)))*dt;

			E = fmin(qrp[i][j][k],E);
			E = fmax(E,0.0);

		} else { E = 0;}
	
		/************************************************************
		* Initiation of cloud ice if temperature is below freezing and 
		* air is saturated with respect to ice
		*************************************************************/
		if(temperature < 273.16 && qis < qvp[i][j][k]){

			nc = n0*exp(beta*(273.16-temperature));

			//pint = fmin(M0*nc,rhou[k]*(qvp[i][j][k]-qcip[i][j][k]));

			pint = fmin(M0*nc/rhou[k],qvp[i][j][k]-qcip[i][j][k]);
		
		} else { pint = 0; nc = 0;}

		/************************************************************
		* Depositional growth of cloud ice
		*************************************************************/
		if(temperature < 273.16 && nc!=0 && qcip[i][j][k]>1e-10){

			App = ( Lv / (Ka*temperature) ) * ( (Ls*Mw) / (Rstar*temperature) - 1);
		
			Bpp = Rstar * temperature / (chi * Mw * esi);

			Ml = qcip[i][j][k] / nc;	//Ml = qcip[i][j][k] * rhou[k] / nc;		

			Dl = 16.3 * sqrt(Ml);		

			pdepi = (4 * Dl * (qis-1)*nc) / ( App + Bpp ) * dt;

		} else { pdepi = 0;}

		/************************************************************
		* Autoconversion of cloud ice to snow
		*************************************************************/
		if(nc!=0 && qcip[i][j][k] > 1e-10){

			qimax = Mmax * nc / rhou[k];
			//printf("%d %d %d %f %f %f\n",i,j,k,qimax,qcip[i][j][k],temperature);

			if(qcip[i][j][k] > qimax)
				pconv = qcip[i][j][k]-qimax;	//pconv = rhou[k] * (qcip[i][j][k]-qimax);
			else
				pconv = 0;

		} else { pconv = 0;}
	
		/************************************************************
		* Collection of cloud ice by snow
		*************************************************************/
		if(temperature < 273.16 && qcip[i][j][k] > 1e-10 && lambdaS > 1e-15){

			psaci = (trigpi * app * qcip[i][j][k] * ESI * N0S) / 4 * pow(p0/pd,0.4) * 2.22 / pow(lambdaS,3.11)*dt;

		} else { psaci = 0;}

		/************************************************************
		* Collection of cloud water by snow
		*************************************************************/
		if(qcp[i][j][k] > 1e-10 && lambdaS > 1e-15){

			psacw = (trigpi * app * qcp[i][j][k] * ESC * N0S) / 4 * pow(p0/pd,0.4) * 2.22 / pow(lambdaS,3.11)*dt;

			if(temperature < 273.16){ psacw_cold = psacw; psacw_warm = 0;}
			else { psacw_warm = psacw; psacw_cold = 0;}

		} else { psacw_cold = 0; psacw_warm = 0;}

		/************************************************************
		* Snow
		*************************************************************/
		if(qsp[i][j][k] > 1e-10 && lambdaS > 1e-15){

			var = 0.65/(lambdaS*lambdaS) +

				0.44 * sqrt(app*rhou[k]/mu) *
 
				pow(p0/pd,0.2) * 1.383 / pow(lambdaS,2.555);

			
			/********************************************************
			* Melting of snow (negative values for melting)
			*********************************************************/
			if(temperature > 273.16){

				psmlt = -(2*trigpi * N0S / Lf) * Ka * (temperature - 273.16) * var; // kg m^-3 s^-1

				psmlt = dt * psmlt / rhou[k]; // mixing ratio for complete time step (kg/kg)

			} else { psmlt = 0;}

			/********************************************************
			* Depositional growth of snow
			*********************************************************/
			if(temperature < 273.16){

				psdep = 4 * (qis-1)*N0S / (App + Bpp) * var;	// kg m^-3 s^-1

				psdep = dt * psdep / rhou[k]; // mixing ratio for complete time step (kg/kg)


			}  else { psdep = 0;}

		} else { psmlt = 0; psdep = 0;}

		//psmlt = 0;psdep = 0;

		/************************************************************
		* Melting of cloud ice
		*************************************************************/
		if(temperature > 273.16 && qcip[i][j][k] > 1e-10){

			psmlti = qcip[i][j][k];

		} else { psmlti = 0;}

		//psaci = 0;
		//psacw_cold = 0;psacw_warm = 0;
		//pconv = 0;
		//pdepi = 0;
		//pint = 0;

		//if(pdepi<0){pdepi=0;}

		/********************************
		* Balance hydrometeors
		*********************************/
	
		qcp[i][j][k] = qcp[i][j][k] + psmlti;

		qcip[i][j][k] = qcip[i][j][k] - psmlti;

		total_convert = fmin(A + B,qcp[i][j][k]);				// ensure autoconversion+accretion does not exceed cloud water
	
		qcp[i][j][k] = qcp[i][j][k] - total_convert;			// remove rain from cloud water
	
		qrp[i][j][k] = qrp[i][j][k] + total_convert	 - E;		// add autoconversion+accretion-evaporation to rain
		
		total_convert = fmin(pint+pdepi,qvp[i][j][k]);

		total_convert = -fmin(-total_convert,qcip[i][j][k]);

		qvp[i][j][k] = qvp[i][j][k] -   total_convert + E;			// evaporation - ice initiation

		qcip[i][j][k] = qcip[i][j][k] + total_convert;

		thp[i][j][k] = thp[i][j][k] - (Lv/(cp*pib[k]))*E 
									+ (Ls/(cp*pib[k]))*(total_convert)
									+ (Lf/(cp*pib[k]))*psmlti;


		;


		total_convert = fmin(pconv + psaci,qcip[i][j][k]);

		qcip[i][j][k] = qcip[i][j][k] - total_convert;			// cloud ice balance
		
		qsp[i][j][k] =  qsp[i][j][k]  + total_convert;


		total_convert = fmin(psacw_cold,qcp[i][j][k]);

		qcp[i][j][k] = qcp[i][j][k] - total_convert;

		qsp[i][j][k] =  qsp[i][j][k] + total_convert;

		thp[i][j][k] = thp[i][j][k] + (Lf/(cp*pib[k]))*total_convert;


		total_convert = fmin(psacw_warm,qcp[i][j][k]);

		qcp[i][j][k] = qcp[i][j][k] - total_convert;

		qrp[i][j][k] =  qrp[i][j][k] + total_convert;


		total_convert = fmax(psmlt,-qsp[i][j][k]);

		qsp[i][j][k] = qsp[i][j][k] + total_convert;

		qrp[i][j][k] =  qrp[i][j][k] - total_convert;

		thp[i][j][k] = thp[i][j][k] + (Lf/(cp*pib[k]))*total_convert;


		total_convert = fmin(psdep,qvp[i][j][k]);

		total_convert = -fmin(-total_convert,qsp[i][j][k]);

		qvp[i][j][k] = qvp[i][j][k] - total_convert;

		qsp[i][j][k] =  qsp[i][j][k] + total_convert;

		thp[i][j][k] = thp[i][j][k] + (Ls/(cp*pib[k]))*total_convert;


//		if(qcip[i][j][k] < -1e-6){ printf("qcip<0 at %d %d %d %f\n",i,j,k,qcip[i][j][k]);}
//		if(qsp[i][j][k] < -1e-6){ printf("qsp<0 at %d %d %d %f\n",i,j,k,qsp[i][j][k]);} 
//		if(qvp[i][j][k] < -1e-6){ printf("qvp<0 at %d %d %d %f\n",i,j,k,qvp[i][j][k]);} 

		/*******************************************
		* Get actual values which may have changed
		********************************************/
		theta = thp[i][j][k]+tb[k];			// actual temperature is base state plus perturbation
	 	vapor = qvp[i][j][k]+qb[k];			// actual vapor is base state plus perturbation
		temperature = theta*pressure;		// temperature
	
		C = 0.0;

		// calculate saturation mixing ratio
	 	qvs = (380.f/pd) * exp(17.27*(temperature-273.0)/(temperature-36.0));
	
		phi = qvs*(17.27*237*Lv)/(cp*pow(temperature-36.,2.));	// latent heating/cooling correction term

		/********************************
		* If parcel is supersaturated
		*********************************/
		if(vapor > qvs){ 

			C = (vapor-qvs)	/(1+phi);						// condensable water
	
		/********************************
		* If parcel is subsaturated
		*********************************/	
		} else if(vapor < qvs && qcp[i][j][k] > 0){

			C = -fmin(qcp[i][j][k],(qvs-vapor)/(1+phi));	// condensable water

		}
	

		qcp[i][j][k] = qcp[i][j][k] + C;						// add condensate to cloud water
		qvp[i][j][k] = qvp[i][j][k] - C;						// remove condensed cloud from vapor

		thp[i][j][k] = thp[i][j][k] + (Lv/(cp*pib[k]))*C
									
					;		// latent heating
	
		acc[i][j][k] = A;
//		cc[i][j][k] = qvs;
//		con[i][j][k] = C;
//		evap[i][j][k] = E;
	
		/********************************
		* Calculate snow fall speed
		*********************************/
		if(qsp[i][j][k] > 0 && lambdaS > 1e-15){ st[i][j][k] = app * 1.15 * pow(lambdaS,-b) * pow(p0/pd,0.4);}
		else { 
			st[i][j][k] = 0;
		}
		


	}}}


	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=1;k<NZ-1;k++){

		vt[i][j][k] = 0.0;

		/********************************
		* Calculate fall speed
		*********************************/
		if(qrp[i][j][k] > 1.0e-12){
	
			qrr = fmax(qrp[i][j][k]*.001*rhou[k],0.0);
			vtden = sqrt(rhou[1]/rhou[k]);
			vt[i][j][k] = 36.34*(pow(qrr,0.1364)) * vtden;
		
		}

	}}}

	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){

		vt[i][j][0] = vt[i][j][1];
		st[i][j][0] = st[i][j][1];

	}}

}




/*********************************************************************
*
*
**********************************************************************/
void sp_bound_microphysics(){

	double coef;

	/*********************************
	* Eastern boundary of domain
	*********************************/
	for(int i=NX-1-iebuffer;i<NX-1;i++){
 	
		// calculate proper coefficients
 		coef = 0.5 * (1.-cos(trigpi*(double)(i+iebuffer-NX+2)/(double)iebuffer));
		//printf("%f ",coef);
		for(int j=1;j<NY-1;j++){
		for(int k=1;k<NZ-1;k++){
 		
 			qvp[i][j][k] = qvp[i][j][k] - coef*(qvp[i][j][k]-qvb[i][j][k]);
			qcp[i][j][k] = qcp[i][j][k] - coef*qcp[i][j][k];
			qrp[i][j][k] = qrp[i][j][k] - coef*qrp[i][j][k];
			qcip[i][j][k] = qcip[i][j][k] - coef*qcip[i][j][k];
			qsp[i][j][k] = qsp[i][j][k] - coef*qsp[i][j][k];
 		}}
 	}
 	//printf("\n");
	/*********************************
	* Western boundary of domain
	**********************************/
	for(int i=1;i<iwbuffer+1;i++){
 
	// calculate proper coefficients
 		coef = .5 * (1.-cos(trigpi*(double)(iwbuffer-i+1)/(double)iwbuffer));
 		//printf("%f ",coef);
		for(int j=1;j<NY-1;j++){
		for(int k=1;k<NZ-1;k++){
		
 			qvp[i][j][k] = qvp[i][j][k] - coef*(qvp[i][j][k]-qvb[i][j][k]);
			qcp[i][j][k] = qcp[i][j][k] - coef*qcp[i][j][k];
			qrp[i][j][k] = qrp[i][j][k] - coef*qrp[i][j][k];
			qcip[i][j][k] = qcip[i][j][k] - coef*qcip[i][j][k];
			qsp[i][j][k] = qsp[i][j][k] - coef*qsp[i][j][k];
 		}}
	}
	//printf("\n");
	/*********************************
	* Southern boundary of domain
	**********************************/
	for(int j=1;j<jsbuffer+1;j++){
 
	// calculate proper coefficients
		coef = .5 * (1.-cos(trigpi*(double)(jsbuffer-j+1)/(double)(jsbuffer)));
 		//printf("%f ",coef);
		for(int i=1;i<NX-1;i++){
		for(int k=1;k<NZ-1;k++){
		
 			qvp[i][j][k] = qvp[i][j][k] - coef*(qvp[i][j][k]-qvb[i][j][k]);
			qcp[i][j][k] = qcp[i][j][k] - coef*qcp[i][j][k];
			qrp[i][j][k] = qrp[i][j][k] - coef*qrp[i][j][k];
			qcip[i][j][k] = qcip[i][j][k] - coef*qcip[i][j][k];
			qsp[i][j][k] = qsp[i][j][k] - coef*qsp[i][j][k];
		}}
	}
	//printf("\n"); 
	/*********************************
	* Northern boundary of domain
	**********************************/
	for(int j=NY-1-jnbuffer;j<NY-1;j++){
 	
		// calculate proper coefficients
		coef = .5 * (1.-cos(trigpi*(double)(j+jnbuffer-NY+2)/(double)jnbuffer));
 		//printf("%f ",coef);
		for(int i=1;i<NX-1;i++){
		for(int k=1;k<NZ-1;k++){
		
 			qvp[i][j][k] = qvp[i][j][k] - coef*(qvp[i][j][k]-qvb[i][j][k]);
			qcp[i][j][k] = qcp[i][j][k] - coef*qcp[i][j][k];
			qrp[i][j][k] = qrp[i][j][k] - coef*qrp[i][j][k];
			qcip[i][j][k] = qcip[i][j][k] - coef*qcip[i][j][k];
			qsp[i][j][k] = qsp[i][j][k] - coef*qsp[i][j][k];
		}}
 	}
	//printf("\n");

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		qvp[i][j][0] = qvp[i][j][1];
		qvp[i][j][NZ-1] = qvp[i][j][NZ-2];
		qcp[i][j][0] = qcp[i][j][1];
		qcp[i][j][NZ-1] = qcp[i][j][NZ-2];
		qrp[i][j][0] = qrp[i][j][1];
		qrp[i][j][NZ-1] = qrp[i][j][NZ-2];
		qsp[i][j][0] = qsp[i][j][1];
		qsp[i][j][NZ-1] = qsp[i][j][NZ-2];
		qcip[i][j][0] = qcip[i][j][1];
		qcip[i][j][NZ-1] = qcip[i][j][NZ-2];
	}}


}


/*********************************************************************
*
*
**********************************************************************/
void microphysics_diffusion(double step){

	double dx2 = dx*dx;
	double dy2 = dy*dy;
	double dz2 = dz*dz;
	

	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=1;k<NZ-1;k++){

		qvp[i][j][k] = qvp[i][j][k] + step*dt* (	
			+ kmixhx * ((qvm[i+1][j][k]-qvb[i+1][j][k]) - 2.0*(qvm[i][j][k]-qvb[i][j][k])+(qvm[i-1][j][k]-qvb[i-1][j][k]))/dx2
			+ kmixhy * ((qvm[i][j+1][k]-qvb[i][j+1][k]) - 2.0*(qvm[i][j][k]-qvb[i][j][k])+(qvm[i][j-1][k])-qvb[i][j-1][k])/dy2	
			+  kmixv * ((qvm[i][j][k+1]-qvb[i][j][k+1]) - 2.0*(qvm[i][j][k]-qvb[i][j][k])+(qvm[i][j][k-1])-qvb[i][j][k-1])/dz2
			)
			;
		qcp[i][j][k] = qcp[i][j][k] + step*dt* (	
			+ kmixhx * (qcm[i+1][j][k] - 2.0*qcm[i][j][k]+qcm[i-1][j][k])/dx2 
			+ kmixhy * (qcm[i][j+1][k] - 2.0*qcm[i][j][k]+qcm[i][j-1][k])/dy2
			+  kmixv * (qcm[i][j][k+1] - 2.0*qcm[i][j][k]+qcm[i][j][k-1])/dz2
			)
			;	
		qrp[i][j][k] = qrp[i][j][k] + step*dt* (	
			+ kmixhx * (qrm[i+1][j][k] - 2.0*qrm[i][j][k]+qrm[i-1][j][k])/dx2
			+ kmixhy * (qrm[i][j+1][k] - 2.0*qrm[i][j][k]+qrm[i][j-1][k])/dy2 
			+  kmixv * (qrm[i][j][k+1] - 2.0*qrm[i][j][k]+qrm[i][j][k-1])/dz2
		);
		qcip[i][j][k] = qcip[i][j][k] + step*dt* (	
			+ kmixhx * (qcim[i+1][j][k] - 2.0*qcim[i][j][k]+qcim[i-1][j][k])/dx2 
			+ kmixhy * (qcim[i][j+1][k] - 2.0*qcim[i][j][k]+qcim[i][j-1][k])/dy2
			+  kmixv * (qcim[i][j][k+1] - 2.0*qcim[i][j][k]+qcim[i][j][k-1])/dz2
			)
			;	
		qsp[i][j][k] = qsp[i][j][k] + step*dt* (	
			+ kmixhx * (qsm[i+1][j][k] - 2.0*qsm[i][j][k]+qsm[i-1][j][k])/dx2
			+ kmixhy * (qsm[i][j+1][k] - 2.0*qsm[i][j][k]+qsm[i][j-1][k])/dy2 
			+  kmixv * (qsm[i][j][k+1] - 2.0*qsm[i][j][k]+qsm[i][j][k-1])/dz2
		);
	}}}
}

/*********************************************************************
*
*
**********************************************************************/
void microphysics_advance_inner(){

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		qv[i][j][k] = qvp[i][j][k];
		qc[i][j][k] = qcp[i][j][k];
		qr[i][j][k] = qrp[i][j][k];
		qci[i][j][k] = qcip[i][j][k];
		qs[i][j][k] = qsp[i][j][k];

	}}}

}

/*********************************************************************
*
*
**********************************************************************/
void microphysics_advance(){

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		qvm[i][j][k] = qvp[i][j][k];
		qcm[i][j][k] = qcp[i][j][k];
		qrm[i][j][k] = qrp[i][j][k];
		qcim[i][j][k] = qcip[i][j][k];
		qsm[i][j][k] = qsp[i][j][k];

		qv[i][j][k] = qvp[i][j][k];
		qc[i][j][k] = qcp[i][j][k];
		qr[i][j][k] = qrp[i][j][k];
		qci[i][j][k] = qcip[i][j][k];
		qs[i][j][k] = qsp[i][j][k];

	}}}

}
