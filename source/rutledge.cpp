
#include "stdafx.h"
#include "budgets.h"
#include "turbulence.h"
#include "interpolate.h"

/*******************************************************************************************
* ICE MICROPHYSICS BASED ON RUTLEDGE AND HOBBS (1983) WITH IMPROVEMENTS BY HONG ET AL. (2004).
* INCLUDES RAIN, SNOW, CLOUD WATER, AND CLOUD ICE.
*
*
********************************************************************************************/

#define ALPHA(x) ( 1.0 + (x-273.15) / (tmax-tmin) )
#define SAT_VAP_WAT(t) ( 611.2 * exp(17.67 * (t-273.15) / (t - 29.65) ) )
#define SAT_VAP_ICE(t) ( 611.2 * exp(21.8745584 * (t-273.15) / (t - 7.66) ) )
#define SAT_MIX_RATIO(e,p) ( 0.62197 * e / (p-e) )


double A,C,E,cvent,vtden,qrr,pd,theta,pressure,phi,vapor,total_convert,total_convert2;
double qvl_sat,qvi_sat; // saturation mixing ratio / vapor pressure with respect to water
double esl,esi; // saturation mixing ratio / vapor pressure with respect to ice

double var;

const double tmin = -30+273.15;	// temperature below which only ice
const double tmax = 0+273.15;		// temperature above which only water


const double minvar_snow = 1.0e-10;
const double minvar_cice = 1.0e-10;
const double minvar_rain = 1.0e-10;
const double minvar_cwat = 1.0e-10;

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
double nc0;

/****************************************************
* Ice deposition coefficients/variables
*****************************************************/
double Rstar = 8.314e3;		// universal gas constant
double Mw = 18.016;			// molecular weight of water
double Ka = 2.43e-2;		// thermal conductivity of air
double chi = 2.26e-5;		// diffusivity of water vapor in air
double MI;					// average mass of cloud ice crystal
double Dl;					// diamter of ice particles
double App;					// thermodynamic term in deposition
double Bpp;					// thermodynamic term in deposition
double pdepi;				// depositional growth
double dimax = 500.0e-6;	// maximum diameter

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
double ESC = 1.0;	// snow/cloud water collection efficiency
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
double mu_s = 1.718e-5;
double psmlt;

/****************************************************
* Despositional growth of snow coefficients/variables
*****************************************************/
double psdep;

/****************************************************
* Evaporation of snow coefficients/variables
*****************************************************/
double pmltev;
double Ap,Bp;

/****************************************************
* Evaporation of rain coefficients/variables
*****************************************************/
const double ap = 3.0e3;		// constant in linear fall speed

double pcond;
//const double Rw = 461;	// gas constant for water vapor

double pcfrz,prfrz;
double a1,a2;
//const double mn = 1.05e-18;

double get_a1(double);
double get_a2(double);
void saturation_adjustment(int,int,int,int);

/****************************************************
* Initialize moisture variables
*****************************************************/
void init_rutledge_microphysics(){


}

/*******************************************************************************
* Perform microphysics calculations
* Handles all of the thermodynamics and phase changes
* for each grid cell
********************************************************************************/
void run_rutledge_microphysics(int il,int ih,int jl,int jh){

	double cpRd = cp/Rd;
	double temperature, qvsat = 0;
	double diabatic;

	double pint_p_pdepi;

	double a0=-0.267,a1=5.15e3,a2=-1.0225e6,a3=7.55e7;	// coefficients in polynomial fallspeed for rain

	double rfall = 2115.0*pow(0.01,0.2) * 17.83786;	// rain fall speed parameter

	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
	
		/**********************************************
		* Get full values of thermodynamic variables
		***********************************************/
		theta = THP(i,j,k) + THBAR(i,j,k) + tb[k];		// full potential temperature is base state plus perturbation		
		pressure = PI(i,j,k)/(cp*tbv[k]) + PBAR(i,j,k);	// full pressure is base state plus perturbation
 		vapor = QVP(i,j,k) + QBAR(i,j,k) + qb[k];		// full vapor is base state plus perturbation
		temperature = theta*pressure;					// actual temperature
		pd = p0*pow(pressure,cpRd);						// dimensional pressure
	
		// calculate saturation mixing ratio
		esl = SAT_VAP_WAT(temperature);
		esi = SAT_VAP_ICE(temperature);

		qvl_sat = SAT_MIX_RATIO(esl,pd);
		qvi_sat = SAT_MIX_RATIO(esi,pd);

		//pcond = (vapor-qvl_sat) / ( 1 + Lv*Lv*qvl_sat / (cp*Rw*temperature*temperature) );

		//if(pcond<0){ pcond = fmax(pcond,-QCP(i,j,k));}

		pcond = 0;

		//------------------------------------------------------
		// If there is snow, calculate slope parameter
		//------------------------------------------------------		
		if(QSP(i,j,k) != 0){
			
			N0S = fmin( 2.0e6*exp( 0.12*(273.16-temperature) ) , 2.0e8)  ;
	
			lambdaS = sqrt(sqrt( (trigpi*rhoS*N0S) / (rhou[k]*QSP(i,j,k)) ));		
			//lambdaS = pow( (trigpi*rhoS*N0S) / (rhou[k]*QSP(i,j,k)),0.25);
		
		}// printf("%f\n",lambdaS);}
		else { lambdaS = 0;}

		/************************************************************
		*
		* 				RAIN PROCESSES
		*
		*************************************************************/
		A = 0; B = 0; E = 0; prfrz = 0;
		//------------------------------------------------------
		// Autoconversion
		//------------------------------------------------------
		if(QCP(i,j,k) > 1.0e-3){ A = 0.001*(QCP(i,j,k)-1.0e-3)*dt;}

		//------------------------------------------------------
		// If rain exceeds some minimum threshold...
		//------------------------------------------------------
		if(QRP(i,j,k) > minvar_rain){
			//------------------------------------------------------
			// Accretion - collection of cloud water by rainwater
			//------------------------------------------------------
			lambdaR = pow( (trigpi*1000*N0r) / (rhou[k]*QRP(i,j,k)),0.25);

			B = 0.25*trigpi*QCP(i,j,k)*ERC*N0r*pow(p0/pd,0.4) *
				2/pow(lambdaR,3) * (a0 + a1*3/lambdaR + a2*12/(lambdaR*lambdaR) + a3*60/(lambdaR*lambdaR*lambdaR)) * dt;
		
			//------------------------------------------------------
			// Evaporation of rain water
			//------------------------------------------------------
			Ap = ( Lv / (Ka*temperature) ) * ( (Lv*Mw) / (Rstar*temperature) - 1);
	
			Bp = Rstar * temperature / (chi * Mw * esl);

			var = 0.78/(lambdaR*lambdaR) + 0.31 * sqrt(ap*rhou[k]/mu_s)/(lambdaR*lambdaR*lambdaR) * 2.0 * pow((p0/pd),0.2);

			E = var * 2*trigpi*N0r* (vapor/qvl_sat-1.0) / (Ap+Bp) * dt * one_d_rhou[k];
		
			//------------------------------------------------------
			// Freezing of rain water if temperature is below freezing
			//------------------------------------------------------		
			if(temperature < 273.16){
			
				prfrz = 20*trigpi*trigpi*100.0*N0r*(1000.0/rhou[k]) * exp(0.66*(273.16-temperature) - 1 ) / (lambdaR*lambdaR*lambdaR*lambdaR*lambdaR*lambdaR*lambdaR);
				
				prfrz = fmin(prfrz * one_d_rhou[k] * dt, QRP(i,j,k) );
			}
			
			//------------------------------------------------------
			// Will end simulation if model becomes unstable
			//------------------------------------------------------
			if(E>100){

				printf("%d %d %d %f %f %f %f %f %f %f %f\n",i,j,k,E,QRP(i,j,k),qvsat,pd,vapor,theta,pressure,temperature);
				fflush(stdout);
				exit(0);
			}
		}


		/****************************************************************
		*
		* 			SOME CLOUD ICE PROCESSES
		*
		*****************************************************************/
		pint = 0; pdepi = 0; pconv = 0; psaci = 0; psmlti = 0;
		//------------------------------------------------------
		// Must be below freezing
		//------------------------------------------------------
		if(temperature < 273.16){
	
			//nc = n0*exp(beta*(273.16-temperature));	// number concentration of ice crystals

			nc0 = 1.0e3*exp(0.1*(273.16-temperature));	// number concentration of ice crystals
			
			nc = fmin(fmax( 5.38e7*pow(rhou[k]*QIP(i,j,k),0.75) ,1.0e3),1.0e6);	// formulation in Hong et al. (2004)
	
			//------------------------------------------------------
			// Initiation of cloud ice if temperature is below 
			// freezing and air is saturated with respect to ice
			//------------------------------------------------------
			if(qvi_sat < vapor){ pint = fmin(M0*nc0/rhou[k],vapor-qvi_sat);}

			//------------------------------------------------------
			// Depositional growth of cloud ice
			//------------------------------------------------------
			if(QIP(i,j,k)>minvar_cice || QSP(i,j,k)>minvar_snow){	// check for snow because it will also need App and Bpp

				App = ( Lv / (Ka*temperature) ) * ( (Ls*Mw) / (Rstar*temperature) - 1);
		
				Bpp = Rstar * temperature / (chi * Mw * esi);

				MI = rhou[k] * QIP(i,j,k) / nc;	// average mass of cloud ice particles

				Dl = 11.3 * sqrt(MI);	// average diameter of cloud ice particles

				pdepi = (4 * Dl * (vapor/qvi_sat-1)*nc) / ( App + Bpp ) * one_d_rhou[k] * dt;
			}
			//------------------------------------------------------
			// Autoconversion of cloud ice to snow
			//------------------------------------------------------
			if(QIP(i,j,k) > minvar_cice){
				//qimax = Mmax * nc / rhou[k];
				
				// critical value suggested by Hong et al. (2004) 
				// to remove strong temperature dependency, kg/kg
				if(temperature>233.16 && pd > 30000){
					qimax = 8.0e-5 * one_d_rhou[k];
				} else {
					qimax = 0.18e-3;
				}
				
				if(QIP(i,j,k) > qimax){ pconv = QIP(i,j,k)-qimax;}
			}
			//------------------------------------------------------
			// Collection of cloud ice by snow
			//------------------------------------------------------
			if(QIP(i,j,k) > minvar_cice && lambdaS > 1e-15){

				ESI = exp( 0.05*(temperature-273.16) );	// suggested by Hong et al. (2004)

				psaci = (trigpi * app * QIP(i,j,k) * ESI * N0S) / 4 * pow(p0/pd,0.4) * 2.22 / pow(lambdaS,3.11)*dt;
			}	
		//------------------------------------------------------
		// Melting of cloud ice if above freezing
		//------------------------------------------------------	
		} else { psmlti = QIP(i,j,k); }
		

		/************************************************************
		*
		* 			SOME SNOW PROCESSES
		*
		*************************************************************/
		psmlt = 0; psdep = 0; pmltev = 0;
		
		if(QSP(i,j,k) > minvar_snow && lambdaS > 1e-15){

			var = 0.65/(lambdaS*lambdaS) +

				0.44 * sqrt(app*rhou[k]/mu_s) *
 
				pow(p0/pd,0.2) * 1.383 / pow(lambdaS,2.555);

			//------------------------------------------------------
			// Melting of snow (negative values for melting)
			//------------------------------------------------------
			if(temperature > 273.16){

				psmlt = -(2*trigpi * N0S / Lf) * Ka * (temperature - 273.16) * var; // kg m^-3 s^-1

				psmlt = dt * psmlt * one_d_rhou[k]; // mixing ratio for complete time step (kg/kg)
			}
			//------------------------------------------------------
			// Depositional growth of snow
			//------------------------------------------------------
			if(temperature < 273.16){
				
				// App and Bpp already calculated with cloud ice deposition
				psdep = 4 * (vapor/qvi_sat-1)*N0S / (App + Bpp) * var;	// kg m^-3 s^-1

				psdep = dt * psdep * one_d_rhou[k]; // mixing ratio for complete time step (kg/kg)

				//printf("%d %d %d %f %f %f\n",i,j,k,psdep,vapor/qvi_sat,vapor/qvl_sat);
			}		
			//------------------------------------------------------
			// Evaporation of melting snow
			//------------------------------------------------------
			if(temperature >= 273.16){
				
				Ap = ( Lv / (Ka*temperature) ) * ( (Lv*Mw) / (Rstar*temperature) - 1);
				Bp = Rstar * temperature / (chi * Mw * esl);
				
				pmltev = 4 * (vapor/qvl_sat-1)*N0S / (Ap + Bp) * var;	// kg m^-3 s^-1

				pmltev = dt * pmltev * one_d_rhou[k]; // mixing ratio for complete time step (kg/kg)
			}
			
			//------------------------------------------------------
			// Will end simulation if model becomes unstable
			//------------------------------------------------------
			if(fabs(psmlt)>100){

				printf("%d %d %d %f %f %f %f %f %f %f %f\n",i,j,k,psmlt,QSP(i,j,k),qvsat,pd,vapor,theta,pressure,temperature);
				fflush(stdout);
				exit(0);
			}
		}

		/************************************************************
		* Collection of cloud water by snow
		*************************************************************/
		psacw_cold = 0; psacw_warm = 0;
		
		if(QCP(i,j,k) > minvar_cwat && lambdaS > 1e-15){

			psacw = (trigpi * app * QCP(i,j,k) * ESC * N0S) / 4 * pow(p0/pd,0.4) * 2.22 / pow(lambdaS,3.11)*dt;

			if(temperature < 273.16){ psacw_cold = psacw;}
			else { 					  psacw_warm = psacw;}
		}

		/************************************************************
		* Freezing of cloud water
		*************************************************************/
		pcfrz = 0;
		
		if(temperature < tmin){ pcfrz = QCP(i,j,k);}

		/************************************************************
		*
		* 			ADD UP ALL OF THE PROCESSES
		*
		*************************************************************/
		pint_p_pdepi = fmin(pint+pdepi,vapor);
		pint_p_pdepi = -fmin(-pint_p_pdepi,QIP(i,j,k));
		
		QVP(i,j,k) += -(pcond + E + psdep + pmltev + pint_p_pdepi);

		QCP(i,j,k) += pcond + psmlti - A - B - psacw_warm - psacw_cold - pcfrz;
		
		QIP(i,j,k) += pint_p_pdepi - psmlti - psaci - pconv + pcfrz;
		
		QRP(i,j,k) += E + A + B - psmlt + psacw_warm - prfrz;
		
		QSP(i,j,k) += psdep + pmltev + psaci + psmlt + psacw_cold + pconv + prfrz;

		diabatic =    Lv/(cp*pib[k]) * (pcond + E + pmltev)
				   +  Ls/(cp*pib[k]) * (pint_p_pdepi + psdep)
				   +  Lf/(cp*pib[k]) * (psmlt + psmlti + psacw_cold + pcfrz + prfrz);
			
		THP(i,j,k) += diabatic;
		
		if(HEAT_BUDGET || PE_BUDGET || PV_BUDGET){
			m_diabatic[INDEX(i,j,k)] = diabatic;
		}

		if(MOISTURE_BUDGET){
			q_diabatic[INDEX(i,j,k)] -= (pcond + E + psdep + pmltev + pint_p_pdepi);	
		}
		
		/********************************
		* Calculate snow fall speed
		*********************************/
		if(QSP(i,j,k) > minvar_snow){
							
				lambdaS = pow( (trigpi*rhoS*N0S) / (rhou[k]*QSP(i,j,k)),0.25);
				
				ST(i,j,k) = app * 1.15 * pow(lambdaS,-b) * pow(p0/pd,0.4);
		}
		else { 
			ST(i,j,k) = 0;
		}
		
		/********************************
		* Calculate rain fall speed
		*********************************/
		if(QRP(i,j,k) > 1.0e-12){
	
			lambdaR = sqrt(sqrt( (trigpi*1000*N0r) / (rhou[k]*QRP(i,j,k)) ));
		
			VT(i,j,k) = rfall / (6.0*pow(lambdaR,0.8)) * pow( p0/pd,0.4);

		} else {
			VT(i,j,k) = 0;
		}


	}}}

	
	saturation_adjustment(il,ih,jl,jh);
	
	
	/********************************
	* Calculate ice fall speed
	*********************************/
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
	
		nc = fmin(fmax( 5.38e7 * pow(rhou[k]*QIP(i,j,k),0.75) ,1.0e3),1.0e6);
	
		MI = rhou[k] * fmax(QIP(i,j,k),0) / nc;	// average mass of cloud ice particles

		Dl = fmax(fmin(11.9 * sqrt(MI),dimax),1.0e-25);
		
		IT(i,j,k) = 1.49e4 * pow(Dl,1.31);//exp(log(Dl)*1.31);

	}}}
	
			
	
	//--------------------------------------------
	// Lower boundary condition
	//--------------------------------------------
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		VT(i,j,0) = VT(i,j,1);
		ST(i,j,0) = ST(i,j,1);
		IT(i,j,0) = IT(i,j,1);

	}}
#if 0
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=0;k<NZ;k++){
		
		if(QSP(i,j,k) < 0){
			
			printf("%f %f %d %f\n",outLons[big_i[i]],outLats[big_j[j]],k,QSP(i,j,k)*1000);
		}
		
	}}}
#endif
	// maybe find out why there are negative values?
	if(PARALLEL)
		zero_moisture(il,ih,jl,jh,fNX*fNY*fNZ);
	else
		zero_moisture(il,ih,jl,jh,NX*NY*NZ);
	//--------------------------------------------
	// Semi-Lagrangian rain fallout
	//--------------------------------------------
	if(RAIN_FALLOUT==2){
		
		hydrometeor_fallout(qrps,vts,il,ih,jl,jh,accRain);
		hydrometeor_fallout(qips,its,il,ih,jl,jh,accSnow);
		hydrometeor_fallout(qsps,sts,il,ih,jl,jh,accSnow);
		
		
		// set rain fall speed to zero so that it
		// does not get advected by subsequent advection
		// calculations
		if(PARALLEL){
			memset(vts,0,fNX*fNY*fNZ*sizeof(double));
			memset(sts,0,fNX*fNY*fNZ*sizeof(double));
			memset(its,0,fNX*fNY*fNZ*sizeof(double));
		} else {
			memset(vts,0,NX*NY*NZ*sizeof(double));
			memset(sts,0,NX*NY*NZ*sizeof(double));
			memset(its,0,NX*NY*NZ*sizeof(double));
		}
	
	} else {
		
		precip_rate(il,ih,jl,jh,vts,qrps,accRain);
		precip_rate(il,ih,jl,jh,sts,qrps,accSnow);
		precip_rate(il,ih,jl,jh,sts,qips,accSnow);	
	}
	
}




/*********************************************************************
* Calculate saturation mixing ratio for water, ice, or mixed phase
* depending on the temperature.
*
* pd -> dimensional pressure (Pa)
* phi_out, phil_out, phii_out -> correction factors for saturation adjustment
* for different phases (l->liquid, i->ice)
*
**********************************************************************/
void get_qvsat(double temperature,double pd,double *qvsat_out,double *phi_out,double *phil_out,double *phii_out){

	double esl,esi,qvsat,f,fi,fl,phi,qvl_sat,qvi_sat,phil=0,phii=0;

	//------------------------------------------------
	// All rain
	//------------------------------------------------
	if(temperature > tmax){

		// calculate saturation mixing ratio
		esl = SAT_VAP_WAT(temperature);
		qvsat = SAT_MIX_RATIO(esl,pd);

		// latent heating/cooling correction term
		f = 17.67 * (273.15-29.65) * Lv / cp;
		phi = pd / (pd-esl) * qvsat * f  /  pow((temperature-29.65),2);

		*phil_out = phi;
		*phii_out = phi;
		*phi_out = phi;

	//------------------------------------------------
	// All ice
	//------------------------------------------------
	} else if(temperature < tmin) {
		// calculate saturation mixing ratio
		esl = SAT_VAP_ICE(temperature);
		qvsat = SAT_MIX_RATIO(esl,pd);

		// latent heating/cooling correction term
		f = 21.8745584 * (273.15-7.66) * Ls / cp;
		phi = pd / (pd-esl) * qvsat * f  /  pow((temperature-7.66),2);

		*phil_out = phi;
		*phii_out = phi;
		*phi_out = phi;

	//------------------------------------------------
	// Mixed phase
	//------------------------------------------------
	} else {

		// calculate saturation mixing ratio
		esl = SAT_VAP_WAT(temperature);
		esi = SAT_VAP_ICE(temperature);

		qvl_sat = SAT_MIX_RATIO(esl,pd);
		qvi_sat = SAT_MIX_RATIO(esi,pd);

		// latent heating/cooling correction term
		fl = 17.67 * (273.15-29.65) * Lv / cp;
		fi = 21.8745584 * (273.15-7.66) * Ls / cp;

		// latent heating/cooling correction term
		phil = pd / (pd-esl) * qvl_sat * fl  /  pow((temperature-29.65),2);
		phii = pd / (pd-esi) * qvi_sat * fi  /  pow((temperature-7.66),2);

		qvsat = ALPHA(temperature) * qvl_sat + (1.0-ALPHA(temperature)) * qvi_sat;
		phi = ALPHA(temperature) * phil + (1.0-ALPHA(temperature)) * phii;

		*phil_out = phil;
		*phii_out = phii;
		*phi_out = phi;
	}


	*qvsat_out = qvsat;

}

/*********************************************************************
*
*
**********************************************************************/
void saturation_adjustment(int il,int ih,int jl,int jh){

	double theta,vapor,temperature,pd,pressure,diabatic;
	double qvsat,phi,phil,phii;
	double cpRd = cp/Rd;

	double Ci,Cl;

	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		/**********************************************
		* Get full values of thermodynamic variables
		***********************************************/
		theta = THP(i,j,k) + THBAR(i,j,k) + tb[k];		// full potential temperature is base state plus perturbation		
		pressure = PI(i,j,k)/(cp*tbv[k]) + PBAR(i,j,k);	// full pressure is base state plus perturbation
 		vapor = QVP(i,j,k) + QBAR(i,j,k) + qb[k];		// full vapor is base state plus perturbation
		temperature = theta*pressure;					// actual temperature
		pd = p0*pow(pressure,cpRd);						// dimensional pressure


		get_qvsat(temperature,pd,&qvsat,&phi,&phil,&phii);

		Cl = 0;
		Ci = 0;

		//----------------------------------------------------------------
		// If parcel is supersaturated
		//----------------------------------------------------------------
		if(vapor > qvsat){		
			//-----------------------------------
			// All liquid
			//-----------------------------------
			if(temperature>tmax){

				Cl = (vapor-qvsat)/(1.0+phi);	// condensable water
				Ci = 0;			
			//-----------------------------------
			// All ice
			//-----------------------------------	
			} else if(temperature<tmin){
				
				Cl = 0;
				Ci = (vapor-qvsat)/(1.0+phi);	
			//-----------------------------------
			// Mixed phase
			//-----------------------------------	
			} else {
				
				Cl = ALPHA(temperature)			* (vapor-qvsat)/(1.0+phi);
				Ci = (1.0-ALPHA(temperature))	* (vapor-qvsat)/(1.0+phi);
			}

			if(USE_TURBULENT_STRESS){isSaturated[INDEX(i,j,k)] = true;}

		//----------------------------------------------------------------
		// If parcel is subsaturated
		//----------------------------------------------------------------
		} else if(vapor < qvsat && (QCP(i,j,k)+QIP(i,j,k)) > 0){
			//-----------------------------------
			// All liquid
			//-----------------------------------
			if(temperature>tmax){

				Cl = -fmin(QCP(i,j,k),(qvsat-vapor)/(1.0+phi));	// condensable water
				Ci = 0;
			//-----------------------------------
			// All ice
			//-----------------------------------	
			} else if(temperature<tmin){
				Cl = 0;
				Ci = -fmin(QIP(i,j,k),(qvsat-vapor)/(1.0+phi));	// condensable ice
				//printf("%d %d %d %e %f %f\n",i,j,k,Ci*1000,temperature,vapor/qvsat);	
			//-----------------------------------
			// Mixed phase
			//-----------------------------------	
			} else {
				// first remove cloud water
				Cl = -fmin(QCP(i,j,k),(qvsat-vapor)/(1.0+phil));
				// if that's not enough, remove cloud ice
				if(Cl+QCP(i,j,k)<=0){

					Ci = -fmin(QIP(i,j,k),(qvsat-vapor+Cl)/(1.0+phii));
				} else {
					Ci = 0;
				}
			}
			
			if(USE_TURBULENT_STRESS){

				if(-(Ci+Cl) < QCP(i,j,k)+QIP(i,j,k)){isSaturated[INDEX(i,j,k)] = true;}
				else {isSaturated[INDEX(i,j,k)] = false;}
			}


		} else if(USE_TURBULENT_STRESS){
			//-------------------------------------------------------
			// If parcel is exactly saturated
			//-------------------------------------------------------			
			if(vapor == qvsat){ isSaturated[INDEX(i,j,k)] = true;}
			//-------------------------------------------------------
			// If parcel is subsaturated
			//-------------------------------------------------------
			else { isSaturated[INDEX(i,j,k)] = false;}
		}

		QIP(i,j,k) += Ci;
		QCP(i,j,k) += Cl;						// add condensate to cloud water
		QVP(i,j,k) -= (Ci + Cl);				// remove condensed cloud from vapor

		diabatic = (Lv/(cp*pib[k]))*Cl + (Ls/(cp*pib[k]))*Ci;		// latent heating

		THP(i,j,k) += diabatic;

		if(HEAT_BUDGET || PE_BUDGET || PV_BUDGET){
			m_diabatic[INDEX(i,j,k)] += diabatic;
		}
				
		if(MOISTURE_BUDGET){
			q_diabatic[INDEX(i,j,k)] -= (Ci + Cl);	
		}

	}}}
	
	
	


}
