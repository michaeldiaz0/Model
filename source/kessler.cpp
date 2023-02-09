
#include "stdafx.h"
#include "budgets.h"

const double minvar = 1.0e-12; // threshold for microphysics calculations


/****************************************************
* Initialize moisture variables
*****************************************************/
void init_kessler_microphysics(){


}

/****************************************************
* Perform microphysics calculations
*****************************************************/
void run_kessler_microphysics(int il,int ih,int jl,int jh){

	double A,B,C,E,cvent,qvsat,pd,theta,pressure,phi,vapor,total_convert;

	const double cpRd = cp/Rd;
	double temperature;
	double dc,de;
	double esl,f;

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

		/**********************************************
		* Autoconversion
		***********************************************/
		if(QCP(i,j,k) > 1.0e-3){ 
			
			A = 0.001*(QCP(i,j,k)-1.0e-3)*dt;}
			
		else { A = 0; }

		/**********************************************
		* Accretion
		***********************************************/
		if(QRP(i,j,k) > minvar){

			B = rhou[k] * 2.2 * QCP(i,j,k) * pow(QRP(i,j,k), 0.875) * dt;

		} else { B = 0; }	
		
		/**********************************************
		* Evaporation
		***********************************************/
		if(QRP(i,j,k) > minvar){

			// calculate saturation mixing ratio
			esl = 611.2 * exp(17.67 * (temperature-273.15) / (temperature - 29.65) );
	
			qvsat = 0.62197 * esl / (pd-esl);

			cvent = 1.6 + 30.3922 * pow( rhou[k]*QRP(i,j,k), 0.2046);

			E = (1.0/rhou[k]) * (( (1.0-vapor/qvsat) * cvent * pow(rhou[k]*QRP(i,j,k),0.525)
								) / ( 2.03e4 + 9.58e6 / (pd*qvsat) ))*dt;

			E = fmin(QRP(i,j,k),E);
			E = fmax(E,0.0);

			if(E>100){

				printf("%d %d %d %f %f %f %f %f %f %f %f\n",i,j,k,E,QRP(i,j,k),qvsat,pd,vapor,theta,pressure,temperature);
				fflush(stdout);
				exit(0);
			}

		} else { E = 0;}
	
		/**************************************
		* Balance water elements
		***************************************/
		QVP(i,j,k) += E;						// evaporate rain into vapor
	
		total_convert = fmin(A + B,QCP(i,j,k));	// ensure autoconversion+accretion does not exceed cloud water
	
		QCP(i,j,k) -= total_convert;			// remove rain from cloud water
	
		QRP(i,j,k) += (total_convert - E);		// add autoconversion+accretion-evaporation to rain
		
		de = (Lv/(cp*pressure))*E;
		
		THP(i,j,k) -= de;		// evaporation cools air temperature
	
		if(HEAT_BUDGET || PE_BUDGET || PV_BUDGET){
			m_diabatic[INDEX(i,j,k)] -= de;
		}
				
		if(MOISTURE_BUDGET){
			q_diabatic[INDEX(i,j,k)] += E;	
		}
		
		/*******************************************
		*
		* 		SATURATION ADJUSTMENT
		*
		********************************************/
		theta = THP(i,j,k)+THBAR(i,j,k)+tb[k];			// full potential temperature is base state plus perturbation
	 	vapor = QVP(i,j,k)+QBAR(i,j,k)+qb[k];			// full vapor is base state plus perturbation
		temperature = theta*pressure;					// actual temperature

		// calculate saturation mixing ratio
		esl = 611.2 * exp(17.67 * (temperature-273.15) / (temperature - 29.65) );
		qvsat = 0.62197 * esl / (pd-esl);
	
		// latent heating/cooling correction term
		f = 17.67 * (273.15-29.65) * Lv / cp;
		phi = pd / (pd-esl) * qvsat * f  /  pow((temperature-29.65),2);

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

		dc = (Lv/(cp*pressure))*C;	// latent heating

		THP(i,j,k) += dc;
		
		if(HEAT_BUDGET || PE_BUDGET || PV_BUDGET){
			m_diabatic[INDEX(i,j,k)] += dc;
		}
				
		if(MOISTURE_BUDGET){
			q_diabatic[INDEX(i,j,k)] -= C;	
		}
		
	}}}



	//--------------------------------------------
	// Semi-Lagrangian rain fallout
	//--------------------------------------------
	if(RAIN_FALLOUT==2){
		
		// need fall speed to calculate Lagrangian trajectories
		calculate_eulerian_fall_speed_rain(vts, qrps, il, ih, jl, jh);
		
		hydrometeor_fallout(qrps,vts,il,ih,jl,jh,accRain);
				
		// set rain fall speed to zero so that it
		// does not get advected by subsequent advection
		// calculations
		if(PARALLEL)
			memset(vts,0,fNX*fNY*fNZ*sizeof(double));
		else
			memset(vts,0,NX*NY*NZ*sizeof(double));
	} else {
		
		precip_rate(il,ih,jl,jh,vts,qrps,accRain);
	}

}

/****************************************************
* Calculate fall speed of rain
*
* vt -> terminal fall speed of rain (output)
* qr -> rain water mixing ratio (input)
*
*****************************************************/
void calculate_rain_fall_speed_kessler(double *vt, double *qr, int il,int ih,int jl,int jh){
	
	double vtden,qrr;
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		vt[INDEX(i,j,k)] = 0.0;

		if(qr[INDEX(i,j,k)] > minvar){

			qrr = fmax(qr[INDEX(i,j,k)]*0.001*rhou[k],0.0);
			vtden = sqrt(rhou[1]/rhou[k]);
			vt[INDEX(i,j,k)] = 36.34*(pow(qrr,0.1364)) * vtden;
		}
	
	}}}
	
	//--------------------------------------------
	// Lower boundary condition
	//--------------------------------------------
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		vt[INDEX(i,j,0)] = vt[INDEX(i,j,1)];
		vt[INDEX(i,j,HTOPO(i,j))] = vt[INDEX(i,j,HTOPO(i,j)+1)];

	}}

}

