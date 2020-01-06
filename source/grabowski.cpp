
#include "stdafx.h"
#include "temperature.h"
#include "surface.h"	// for t_diffusion for PV budget

#define TAUA(x) (3.55555*(x-273.15)*(x-273.15) + 106.66666*(x-273.15) + 1000.0)
#define ALPHA(x) ( 1.0 + (x-273.15) / (tmax-tmin) )
#define SAT_VAP_WAT(t) ( 611.2 * exp(17.67 * (t-273.15) / (t - 29.65) ) )
#define SAT_VAP_ICE(t) ( 611.2 * exp(21.8745584 * (t-273.15) / (t - 7.66) ) )
#define SAT_MIX_RATIO(e,p) ( 0.62197 * e / (p-e) )
//#define SAT_VAP_WAT(temperature) ( 611.2 * exp( (2.53e6/461.5) * (1.0/273.16 - 1.0/temperature) ) )
//#define SAT_VAP_ICE(temperature) ( 611.2 * exp( (2.84e6/461.5) * (1.0/273.16 - 1.0/temperature) ) )
#define AUTO_CONV_WAT(Nd,Dd,psi) ( 1.67e-5 * psi*psi * 1.0 / (5.0 + 0.036*Nd / (Dd*psi) )  )
#define PSI_AUTO_CONV(qc,rho) (1.0e3*rho*qc)
#define AUTO_CONV2(qc) (0.001*(qc-1.0e-3))
#define LAMBDA_ICE(N0,rho,q) (pow( (2.5e-2 * N0)*2.0 / (rho * q)  , 0.333333))
#define LAMBDA_RAIN(N0,rho,q) (pow( (523.6 * N0)*6.0 / (rho * q)  , 0.25))
#define G_OF_TE(esl,t) ( 1.0e-7 / (2.2*t/esl + 2.2e2/t) )

const double N0 = 1.0e7; 	 // (m^-4)
const double Nd = 50;	 // concentration of cloud droplets (number / cm^3)
const double Dd = 0.366; // relative dispersion of cloud droplet population

const double tmin = -20+273.15;	// temperature below which only ice
const double tmax = 0+273.15;		// temperature above which only water

//double kmixv_moisture[NZ];

/****************************************************
* Initialize moisture variables
*****************************************************/
void init_microphysics(){

	//printf("%e %e\n",AUTO_CONV_WAT(Nd,Dd,PSI_AUTO_CONV( 0.002,1.0 )),AUTO_CONV2(0.002));

	for(int t=-10;t<10;t++){
		//printf("%d %f\n",t,ALPHA((double)t));
	}

	for(int t=200;t<275;t++){

		printf("%d %f\n",t,SAT_VAP_WAT( (double)t) -SAT_VAP_ICE( (double)t) ) ;
	}

	exit(0);
}

/****************************************************
* Perform microphysics calculations
*****************************************************/
void kessler_microphysics(int il,int ih,int jl,int jh){

	double A,B,C,E,cvent,vtden,qrr,qvsat,pd,theta,pressure,phi,vapor,total_convert;

	const double minvar = 1.0e-12; // threshold for microphysics calculations

	const double cpRd = cp/Rd;
	double temperature;
	double dc,de;
	double esl,f;
	double lambda,nbar,mbar,D,vt,fvent,g,dep,qs,qr;

	double esi,qvl_sat,qvi_sat,fl,fi,phil,phii,qi,ql;
	double maxDep = 0,maxSnow = 0;

	double lambda_L,lambda_I,nbar_L,nbar_I,mbar_L,mbar_I,D_L,D_I,vt_L,vt_I,fvent_L,fvent_I,g_L,g_I,dep_L,dep_I;
	
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
		
		/**************************************************************************
		* 							AUTOCONVERSION
		***************************************************************************/
		if(QCP(i,j,k) > minvar){
		
			if( temperature > tmax){
			
				A = AUTO_CONV_WAT( Nd, Dd, PSI_AUTO_CONV( QCP(i,j,k),rhou[k] ) ) * one_d_rhou[k] * dt;
				
			} else if( temperature < tmin){
			
				A = QCP(i,j,k) / TAUA(temperature) * dt;
				
			} else {
				
				A =  AUTO_CONV_WAT( Nd, Dd, PSI_AUTO_CONV( ALPHA(temperature)*QCP(i,j,k),rhou[k] ) ) * one_d_rhou[k] * dt +
					(1.0-ALPHA(temperature)) * QCP(i,j,k) / TAUA(temperature) * dt;
			}
		} else { A = 0;}

		/**************************************************************************
		* 							DEPOSITION
		***************************************************************************/
		dep = 0;
		B = 0;
		E = 0;
		
		if(QRP(i,j,k) > minvar){
			//--------------------------------------------------------
			// 							All Rain
			//--------------------------------------------------------
			if(temperature >= tmax){
			
				qr = QRP(i,j,k);
			 			
				esl = SAT_VAP_WAT(temperature);
				qvl_sat = SAT_MIX_RATIO(esl,pd);
		
				lambda_L = LAMBDA_RAIN( N0,rhou[k],qr);
		
			 	// mean concentration of precipitation particles
			 	nbar_L = N0 / lambda_L;
			
				// mean mass of a precipitation particle
				mbar_L = rhou[k] * qr / nbar_L;
		
				// diameter of mean particle
				D_L = pow( mbar_L / 523.6, 0.33333);	

				// sedimentation velocity
				vt_L = 130.0 * pow(D_L,0.5 );		

				// ventilation factor
				fvent_L = 0.78 + 0.27 * sqrt(D_L*vt_L / 2.0e-5);			

				g_L = G_OF_TE(esl,temperature);

				dep_L = nbar_L * 4.0*trigpi * D_L / 2.0 * ( vapor / qvl_sat - 1.0) * fvent_L * g_L;
			
				dep = dep_L * one_d_rhou[k] * dt;
			
				B = (nbar_L * 0.7854 * D_L*D_L * vt_L * 0.8 * 1.0 * QCP(i,j,k) ) * dt;
				
		 	}
			//--------------------------------------------------------
			// 							All ice
			//--------------------------------------------------------
			else if(temperature < tmin){ 
				
				qs = QRP(i,j,k);
						
				esi = SAT_VAP_ICE(temperature);
				qvi_sat = SAT_MIX_RATIO(esi,pd);
			
				lambda = LAMBDA_ICE(N0,rhou[k],qs);
			
				nbar = N0 / lambda; 	// mean concentration of precipitation particles
			
				mbar = rhou[k] * qs / nbar;	// mean mass of a precipitation particle
			
				D = sqrt(mbar / 2.5e-2);	// diameter of mean particle
			
				vt = 4.0 * pow(D,0.25);		// sedimentation velocity
			
				fvent = 0.65 + 0.39 * sqrt(D*vt / 2.0e-5);	// ventilation factor
			
				g = G_OF_TE(esi,temperature);
			
				dep = nbar * 4.0*trigpi * D / 3.0 * ( vapor / qvi_sat - 1.0) * fvent * g * one_d_rhou[k] * dt;
				
				B = nbar * 0.7854 * D*D * vt * 0.2 * 0.3 * QCP(i,j,k) * dt;
			}
			
			//--------------------------------------------------------
			// 							Mixed phase
			//--------------------------------------------------------
			else {
				
				qs = (1.0-ALPHA(temperature)) * QRP(i,j,k);
				qr = ALPHA(temperature) * QRP(i,j,k);
				 
 				esi = SAT_VAP_ICE(temperature);
 				qvi_sat = SAT_MIX_RATIO(esi,pd);
				
 				esl = SAT_VAP_WAT(temperature);
 				qvl_sat = SAT_MIX_RATIO(esl,pd);
			
 				lambda_L = LAMBDA_RAIN( N0,rhou[k],qr);
 				lambda_I = LAMBDA_ICE(  N0,rhou[k],qs);
			
			 	// mean concentration of precipitation particles
			 	nbar_L = N0 / lambda_L; 
 				nbar_I = N0 / lambda_I;
				
				// mean mass of a precipitation particle
 				mbar_L = rhou[k] * qr / nbar_L;			
 				mbar_I = rhou[k] * qs / nbar_I;
			
				// diameter of mean particle
				D_L = pow( mbar_L / 523.6, 0.33333);	
 				D_I = sqrt(mbar_I / 2.5e-2);	
 
				// sedimentation velocity
  				vt_L = 130.0 * pow(D_L,0.5 );				
 				vt_I =   4.0 * pow(D_I,0.25);		

				// ventilation factor
 				fvent_L = 0.78 + 0.27 * sqrt(D_L*vt_L / 2.0e-5);
 				fvent_I = 0.65 + 0.39 * sqrt(D_I*vt_I / 2.0e-5);			

 				g_L = G_OF_TE(esl,temperature);			
 				g_I = G_OF_TE(esi,temperature);

 				dep_L = nbar_L * 4.0*trigpi * D_L / 2.0 * ( vapor / qvl_sat - 1.0) * fvent_L * g_L;			
 				dep_I = nbar_I * 4.0*trigpi * D_I / 3.0 * ( vapor / qvi_sat - 1.0) * fvent_I * g_I;
				
				dep = (dep_L + dep_I) * one_d_rhou[k] * dt;
				
				B = (	
						(nbar_L * 0.7854 * D_L*D_L * vt_L * 0.8 * 1.0 * QCP(i,j,k)) + 
						(nbar_I * 0.7854 * D_I*D_I * vt_I * 0.2 * 0.3 * QCP(i,j,k))
					) * dt;
					
			 }
			 
 			if(dep>100.0 || dep<-100.0){

 				printf("%d %d %d %f %f %f %f %f %f %f %f\n",i,j,k,E,QRP(i,j,k),qvsat,pd,vapor,theta,pressure,temperature);
 				fflush(stdout);
 				exit(0);
 			}
			
			//if(dep > maxDep){ maxDep = dep; maxSnow = vapor / qvi_sat;}
		}

		E = -dep;

		/**************************************
		* Balance water elements
		***************************************/
		QVP(i,j,k) += E;						// evaporate rain into vapor
	
		total_convert = fmin(A + B,QCP(i,j,k));	// ensure autoconversion+accretion does not exceed cloud water
	
		QCP(i,j,k) -= total_convert;			// remove rain from cloud water
	
		QRP(i,j,k) += (total_convert - E);		// add autoconversion+accretion-evaporation to rain
		
		de = (Lv/(cp*pressure))*E;
		
		THP(i,j,k) -= de;		// evaporation cools air temperature
	
//		#if HEAT_BUDGET
		//#if VORTICITY_BUDGET
		//cond[INDEX(i,j,k)] -= de;	// for heat budget
//		#endif
		//evap[INDEX(i,j,k)] -= de;
		//#elif PV_BUDGET
		//cond[INDEX(i,j,k)] = de;
		//evap[INDEX(i,j,k)] = -de;	// for PV budget	
		//#endif
		
		
		/******************************************************
		*
		* 		SATURATION ADJUSTMENT
		*
		*******************************************************/
		theta = THP(i,j,k)+THBAR(i,j,k)+tb[k];			// full potential temperature is base state plus perturbation
	 	vapor = QVP(i,j,k)+QBAR(i,j,k)+qb[k];			// full vapor is base state plus perturbation
		temperature = theta*pressure;					// actual temperature

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

		//------------------------------------------------
		// All snow
		//------------------------------------------------			
		} else if(temperature < tmin) {
			// calculate saturation mixing ratio
			esl = SAT_VAP_ICE(temperature);
			qvsat = SAT_MIX_RATIO(esl,pd);
	
			// latent heating/cooling correction term
			f = 21.8745584 * (273.15-7.66) * Lv / cp;
			phi = pd / (pd-esl) * qvsat * f  /  pow((temperature-7.66),2);

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
			fi = 21.8745584 * (273.15-7.66) * Lv / cp;

			// latent heating/cooling correction term
			phil = pd / (pd-esl) * qvl_sat * fl  /  pow((temperature-29.65),2);
			phii = pd / (pd-esi) * qvi_sat * fi  /  pow((temperature-7.66),2);
			
			qvsat = ALPHA(temperature) * qvl_sat + (1.0-ALPHA(temperature)) * qvi_sat;
			phi = ALPHA(temperature) * phil + (1.0-ALPHA(temperature)) * phii;
			
		}

		C = 0.0;	// condensed vapor

		/**************************************
		* If parcel is supersaturated
		***************************************/
		if(vapor > qvsat){

			C = (vapor-qvsat)/(1.0+phi);	// condensable water
		
			if(USE_TURBULENT_STRESS){isSaturated[INDEX(i,j,k)] = true;}
		
		/*************************************
		* If parcel is subsaturated
		**************************************/
		} else if(vapor < qvsat && QCP(i,j,k) > 0){

			C = -fmin(QCP(i,j,k),(qvsat-vapor)/(1.0+phi));	// condensable water
		
			if(USE_TURBULENT_STRESS){
		
				if(-C < QCP(i,j,k)){isSaturated[INDEX(i,j,k)] = true;}
				else {isSaturated[INDEX(i,j,k)] = false;}
			}
		
		} else if(USE_TURBULENT_STRESS && vapor == qvsat){
		
			isSaturated[INDEX(i,j,k)] = true;
		
		} else if(USE_TURBULENT_STRESS){
		
			isSaturated[INDEX(i,j,k)] = false;
		}
				
		/***********************************
		* Exchanges between clouds and vapor
		************************************/	
		QCP(i,j,k) += C;					// add condensate to cloud water
		QVP(i,j,k) -= C;					// remove condensed cloud from vapor

		dc = (Lv/(cp*pressure))*C;

		THP(i,j,k) += dc;

		
//		#if HEAT_BUDGET
		//#if VORTICITY_BUDGET
		//cond[INDEX(i,j,k)] += dc;	// for heat budget
		//#elif PV_BUDGET
		//evap[INDEX(i,j,k)] += dc + t_diffusion[INDEX(i,j,k)] * dt; // for PV budget 		
		//#endif

	}}}

	//printf("%e %f\n",1000.0*maxDep/dt,100.0*maxSnow);

	double maxFall = 0;
	int kmax = 0;

	double app = 1.139;		// fall speed coefficient
	double b = 0.11;		// fall speed exponent
	double lambdaS;			// slope factor
	double N0S = 8e6;		// intercept value in snowflake size distribution
	double rhoS = 200;		// density of snow (kg/m3)
	double sv,rv,frac;

	/*******************************************
	* Calculate fall speed
	********************************************/
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		VT(i,j,k) = 0.0;

		if(QRP(i,j,k) > minvar){

			pressure = PI(i,j,k)/(cp*tbv[k]) + PBAR(i,j,k);
			pd = p0*pow(pressure,cpRd);
			theta = THP(i,j,k)+THBAR(i,j,k)+tb[k];
			temperature = theta*pressure;
			
			if(temperature>tmax){
			
				qrr = fmax(QRP(i,j,k)*0.001*rhou[k],0.0);
				vtden = sqrt(rhou[1]/rhou[k]);
				
				VT(i,j,k) = 36.34*(pow(qrr,0.1364)) * vtden;
				
			} else if(temperature > tmin && temperature <= tmax){ 
				
				qrr = fmax(ALPHA(temperature)*QRP(i,j,k)*0.001*rhou[k],0.0);
				vtden = sqrt(rhou[1]/rhou[k]);
				
				rv = 36.34*(pow(qrr,0.1364)) * vtden;
				
				lambdaS = pow( (trigpi*rhoS*N0S) / (rhou[k]*(1.0-ALPHA(temperature))*QRP(i,j,k)),0.25);
				
				sv = app * 1.15 * pow(lambdaS,-b) * pow(p0/pd,0.4);
								
				VT(i,j,k) =  (1.0-ALPHA(temperature)) * sv + ALPHA(temperature) * rv;
				
			} else {
								
				lambdaS = pow( (trigpi*rhoS*N0S) / (rhou[k]*QRP(i,j,k)),0.25);
				
				VT(i,j,k) = app * 1.15 * pow(lambdaS,-b) * pow(p0/pd,0.4);
				
				//if(VT(i,j,k) > maxFall){ maxFall = VT(i,j,k); kmax = k; }
			}

		}
		
		
	}}}

	//printf("%d %f\n",kmax,maxFall);
	/*******************************************
	* Set lower boundary
	********************************************/
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		VT(i,j,0) = VT(i,j,1);
		VT(i,j,HTOPO(i,j)) = VT(i,j,HTOPO(i,j)+1);
	}}

}

/*********************************************************************
*
*
**********************************************************************/
double get_CAPE(int i,int j,int k_p){
	
	double esl,qvsat,theta,vapor,temp,pressure,pd,f,phi,vtemp,temperature;
	double p_temp,p_vapor,p_vtemp;
	double c,buoy,cape;
	double cpRd = cp/Rd;

	//----------------------------------------------
	// Initialize parcel's temperature and
	// water vapor content
	//----------------------------------------------
	p_vapor = QBAR(i,j,k_p) + qb[k_p];	// full vapor is base state plus perturbation
	p_temp = THP(i,j,k_p) + THBAR(i,j,k_p) + tb[k_p];		
			
	//----------------------------------------------
	// Let the parcel ascend...
	//----------------------------------------------
	for(int k=k_p;k<NZ;k++){

		//----------------------------------------------
		// The ambient environment
		//----------------------------------------------
		theta = THP(i,j,k) + THBAR(i,j,k) + tb[k];		// full potential temperature is base state plus perturbation		
		pressure = PI(i,j,k)/(cp*tbv[k]) + PBAR(i,j,k);	// full pressure is base state plus perturbation			
		temperature = theta*pressure;					// actual temperature
		vapor = QBAR(i,j,k) + qb[k];
		pd = p0*pow(pressure,cpRd);

		//----------------------------------------------
		// The parcel
		//----------------------------------------------	
		esl = SAT_VAP_WAT(p_temp*pressure);		//calculate parcel's saturation vapor pressure	
		qvsat = SAT_MIX_RATIO(esl,pd);

		f = 17.67 * (273.15-29.65) * Lv / cp;				// correction to saturation vapor pressure which accounts for
		phi = pd / (pd-esl) * qvsat * f  /  ((p_temp-29.65)*(p_temp-29.65));	// latent heat release
			
		c = (p_vapor-qvsat)/(1+phi);		// remove condensate from parcel
		c = fmax(c,0.0);

		p_vapor -= c;
		p_temp  += Lv * c / (cp*pressure);

		//  calculate parcel virtual potential temperature
		p_vtemp = p_temp * (1 + 0.61*p_vapor);
		vtemp = theta * (1 + 0.61*vapor);

		// check for buoyancy
		buoy = p_vtemp - vtemp;

		// if positively buoyant then add to CAPE
		cape += grav*(zu[k]-zu[k-1])*fmax(buoy,0.0) / vtemp;

		//esl = SAT_VAP_WAT(p_temp*pressure);

		//qvsat = SAT_MIX_RATIO(esl,pd);

		//printf("%d %f %f %f %f %f\n",k,theta,p_temp,p_temp*pressure,p_vapor*1000,qvsat*1000);

	}

	return cape;
	
}

/*********************************************************************
*
*
**********************************************************************/
void zero_moisture(int il,int ih,int jl,int jh,int size){
	
	for(int i=0;i<size;i++) 
		if(qcps[i] < 0){ qcps[i] = 0;}
	
	for(int i=0;i<size;i++)
		if(qrps[i] < 0){ qrps[i] = 0;}
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
	
		// Remove negative values for moisture variables
		if(QVP(i,j,k)+QBAR(i,j,k)+qb[k] < 0){ QVP(i,j,k) = -(QBAR(i,j,k)+qb[k]);}
		//if(QCP(i,j,k) < 0){ QCP(i,j,k) = 0;}
		//if(QRP(i,j,k) < 0){ QRP(i,j,k) = 0;}
	}}}
	
}

/*********************************************************************
* Substeps within RK3 loop
**********************************************************************/
void microphysics_advance_inner(size_t num_bytes){

	switch_array(&qvs,&qvps);
	switch_array(&qcs,&qcps);
	switch_array(&qrs,&qrps);	
}

/*********************************************************************
* Final step for RK3
**********************************************************************/
void microphysics_advance(size_t num_bytes){

	memcpy(qvms,qvps,num_bytes);
	memcpy(qcms,qcps,num_bytes);
	memcpy(qrms,qrps,num_bytes);

	microphysics_advance_inner(num_bytes);
}
