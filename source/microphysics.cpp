#include "stdafx.h"

/*******************************************************************************
* 
********************************************************************************/
#define SAT_VAP_WAT(t) ( 611.2 * exp(17.67 * (t-273.15) / (t - 29.65) ) )
#define SAT_VAP_ICE(t) ( 611.2 * exp(21.8745584 * (t-273.15) / (t - 7.66) ) )
#define SAT_MIX_RATIO(e,p) ( 0.62197 * e / (p-e) )

//double kmixv_moisture[NZ];
double *raintotal;

/***********************************************************************************
* 
************************************************************************************/
void init_microphysics(){

	int nx,ny;

	if(PARALLEL){ nx = fNX; ny = fNY;}
	else { nx = NX; ny = NY;}

	raintotal = (double *)calloc(nx*ny,sizeof(double));
	//-------------------------------------------------------------------
	// If initializing from an output file, store the precipitation total
	//-------------------------------------------------------------------
	if(isRestartRun || (PERTURBATION_OPTION == 0 && perturbationFileTime > 0)){
		
		for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
		
			raintotal[ny*i+j] = FRICTION(i,j,NZ-2);
			
			//if(i==nx/2){ printf("%f ",raintotal[ny*i+j]);}
			
		}}
	}

	//run_microphysics = &rutledge_microphysics;
}

/******************************************************
* 
*******************************************************/
void run_microphysics(int il,int ih,int jl,int jh){
	
	switch(MICROPHYSICS_OPTION){
		
		//-------------------------------------------------
		// do nothing
		//-------------------------------------------------
		case 0:
		
			break;
		//-------------------------------------------------
		// Kessler microphysics
		//-------------------------------------------------	
		case 1:
		
			run_kessler_microphysics(il,ih,jl,jh);
			
			break;
		//-------------------------------------------------
		// Based on Rutledge and Hobbs (1983), Hong et al. (200x)
		//-------------------------------------------------			
		case 2:
		
			run_rutledge_microphysics(il,ih,jl,jh);
				
			break;
		//-------------------------------------------------
		// Unsupported option
		//-------------------------------------------------			
		default:
		
			printf("Error! Microphysics option %d not supported",MICROPHYSICS_OPTION);
			
			break;
	
	}
}

/*********************************************************************
* Accumulated rainfall (mm)
*
* il,ih,jl,jh - array index bounds
**********************************************************************/
void rain_rate(int il,int ih,int jl,int jh){
	
	int ny;
	
	if(PARALLEL){ ny = fNY;}
	else { ny = NY;}
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
			
		raintotal[ny*i+j] += rhow[1] * 0.5*(VT(i,j,0)+VT(i,j,1)) * 0.5*(QRP(i,j,0)+QRP(i,j,1)) * dt;
		
		FRICTION(i,j,NZ-2) = raintotal[ny*i+j];
	}}
}

/*********************************************************************
* Calculate fallout of hydrometeors using a forward in time 
* semi-Lagrangian method based on Juang and Hong (2010). Useful for
* longer time steps to avoid numerical instabilities.
*
* var - the hydrometeor mixing ratio field
* vel - the hydrometeor terminal velocity field
* il,ih,jl,jh - array index bounds
*
**********************************************************************/
void hydrometeor_fallout(double *var,double *vel,int il,int ih,int jl,int jh){
	
	const double DeCFL0 = 0.05;	// tunable paramter to maintain numerical stability
	const double minVar = 1e-10;
	
	double ZA_edges[NZ];
	double QA_mass[NZ];
	double QD_mass[NZ];
	double Vel_edges[NZ];
	double DeCFL;
	
	bool hasRain;
	
	//-----------------------------------------------------------
	// Initialize local variables
	//-----------------------------------------------------------
	for(int k=0;k<NZ;k++){
		
		ZA_edges[k] = ZW(k);
		Vel_edges[k] = 0;
		QD_mass[k] = 0;
		QA_mass[k] = 0;
	}
	
	//-----------------------------------------------------------
	// For each grid point in local domain...
	//-----------------------------------------------------------
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
		
		hasRain = false;
	
		for(int k=1;k<NZ;k++){ 
			if(var[INDEX(i,j,k)] > minVar){ hasRain = true; break;}
		}
		
		//---------------------------------------------------------------
		// Calculate only if some place in the column has exceeded
		// a critical value (to avoid unnecessary calculations).
		//---------------------------------------------------------------		
		if(hasRain){
		
			//-----------------------------------------------------------
			// velocity at upper and lower edges of grid cells
			//-----------------------------------------------------------
			for(int k=1;k<NZ;k++){ Vel_edges[k] = 0.5*( vel[INDEX(i,j,k)]+vel[INDEX(i,j,k-1)]);}
		
			Vel_edges[0] = Vel_edges[1];	
			//-----------------------------------------------------------
			// Modify sedimentation velocity to satisfy DeCFL criterion
			//-----------------------------------------------------------
			for(int k=NZ-2;k>=0;k--){
			
				DeCFL = ( Vel_edges[k+1] - Vel_edges[k]) * dt / DZU(k);
	
				if(DeCFL>=DeCFL0){
					Vel_edges[k] = Vel_edges[k+1] - DeCFL0 * DZU(k) / dt;
				}
			}	
			//-----------------------------------------------------------
			// Advect cell edges
			//-----------------------------------------------------------
			for(int k=0;k<NZ-1;k++){ ZA_edges[k] = ZW(k) - dt * Vel_edges[k];
				if(i==103 && j==9 && rank==1){ printf("%d %f %f %f\n",k,ZW(k),ZA_edges[k],Vel_edges[k]);}
			}
			//-----------------------------------------------------------
			// hydrometeor mass of advected parcel
			//-----------------------------------------------------------		
			for(int k=0;k<NZ-1;k++){
		
				QA_mass[k] = rhou[k]*var[INDEX(i,j,k)] * (ZW(k+1)-ZW(k)) / (ZA_edges[k+1]-ZA_edges[k]);
				
				if(i==103 && j==9 && rank==1){ printf("%d %f %f %f %f\n",k,rhou[k]*var[INDEX(i,j,k)]*1000,QA_mass[k]*1000,ZA_edges[k+1],ZA_edges[k]);}
			}	
			//-----------------------------------------------------------
			// Interpolation
			//-----------------------------------------------------------
			piecewise_interp(ZA_edges,&ZW(0),QA_mass,QD_mass,NZ);	
			//-----------------------------------------------------------
			// Divide by density for final mass
			//-----------------------------------------------------------
			for(int k=0;k<NZ;k++){ var[INDEX(i,j,k)] = QD_mass[k] / rhou[k];
				//if(i>=100 && i<=105 && j>=7 && j<=12 && rank==1){ printf("%d %d %d %f\n",i,j,k,var[INDEX(i,j,k)]*1000);
				if(i==103 && j==9 && rank==1){ printf("%f %f %f %d %d %d %f\n",LONS(i),LATS(j),ZU(k),i,j,k,var[INDEX(i,j,k)]*1000);}
			}
		
		}}
	}
}


/*********************************************************************
* Piecewise constant interpolation which preserves total mass
*
* zin - input height levels
* zout - output height levels
* qin - input field
* qout - output field
* zlevs number of levels (i.e. array length for all input arrays)
**********************************************************************/
void piecewise_interp(double *zin,double *zout,double *qin,double *qout,int zlevs){

	double factor;

	for(int i=0;i<zlevs;i++){ qout[i] = 0;}

	for(int i=1;i<zlevs-1;i++){
	for(int j=1;j<zlevs-1;j++){
		
		if(zout[i] > zin[j]){
			
			if(zout[i+1] < zin[j+1] ){ // mA3->m1
				
				qout[i] += qin[j];
				
			} else if(zout[i+1] == zin[j+1] ){
			
				factor = (zout[i] - zin[j]) / (zout[i+1]-zout[i]);
			
				qout[i] += factor * qin[j];
			}
		}
		
		if(zout[i] > zin[j] && zout[i] < zin[j+1]){
			
			 if(zout[i-1] < zin[j] ){ // mA4-> m2
			
				factor = ( zout[i] - zin[j]) / (zout[i]-zout[i-1]);
			
				qout[i-1] += factor * qin[j];	
			}
		
			 if(zout[i+1] > zin[j+1] ){ // mA4 -> m3
				
				factor = (zin[j+1] - zout[i]) / (zout[i+1]-zout[i]);
			
				qout[i] += factor * qin[j];
			} 
		}
				
		if(zout[i] == zin[j]){
		
			 if(zout[i+1] > zin[j+1] ){
			
				factor = (zin[j+1] - zout[i]) / (zout[i+1]-zout[i]);
				
				qout[i] += factor * qin[j];
				
			} else if(zout[i+1] == zin[j+1] ){
				
				qout[i] += qin[j];
			}
		}

	}}

}


/*********************************************************************
* Substeps within RK3 loop
**********************************************************************/
void microphysics_advance_inner(size_t num_bytes){

	switch_array(&qvs,&qvps);
	switch_array(&qcs,&qcps);
	switch_array(&qrs,&qrps);
	
	if(USE_ICE){
		switch_array(&qis,&qips);
		switch_array(&qss,&qsps);
	}	
}

/*********************************************************************
* Final step for RK3
**********************************************************************/
void microphysics_advance(size_t num_bytes){

	memcpy(qvms,qvps,num_bytes);
	memcpy(qcms,qcps,num_bytes);
	memcpy(qrms,qrps,num_bytes);
	
	if(USE_ICE){
		memcpy(qims,qips,num_bytes);
		memcpy(qsms,qsps,num_bytes);
	}

	microphysics_advance_inner(num_bytes);
}

/*********************************************************************
* Remove negative values for moisture variables
*
**********************************************************************/
void zero_moisture(int il,int ih,int jl,int jh,int size){
	
	for(int i=0;i<size;i++) 
		if(qcps[i] < 0){ qcps[i] = 0;}
	
	for(int i=0;i<size;i++)
		if(qrps[i] < 0){ qrps[i] = 0;}

	if(USE_ICE){

		for(int i=0;i<size;i++)
			if(qips[i] < 0){ qips[i] = 0;}
	
		for(int i=0;i<size;i++)
			if(qsps[i] < 0){ qsps[i] = 0;}
	}
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
	
		if(QVP(i,j,k)+QBAR(i,j,k)+qb[k] < 0){ QVP(i,j,k) = -(QBAR(i,j,k)+qb[k]);}

	}}}
	
}

/*********************************************************************
*
*
**********************************************************************/
double get_CAPE(int i,int j,int k_p){
	//----------------------------------------------
	// Initialize parcel's temperature and
	// water vapor content
	//----------------------------------------------
	double p_vapor = QV(i,j,k_p) + QBAR(i,j,k_p);// + qb[k_p];	// full vapor is base state plus perturbation
	double p_temp = TH(i,j,k_p) + THBAR(i,j,k_p) ;//+ tb[k_p];
	
	
	return get_CAPE(i,j,k_p,p_vapor,p_temp);
}

/*********************************************************************
*
*
**********************************************************************/
double get_CAPE(int i,int j,int k_p,double p_vapor,double p_temp){
	
	double esl,qvsat,theta,vapor,pressure,pd,f,phi,vtemp,temperature;
	//double p_temp,p_vapor
	double p_vtemp;
	double c,buoy,cape = 0;
	double cpRd = cp/Rd;

	//----------------------------------------------
	// Initialize parcel's temperature and
	// water vapor content
	//----------------------------------------------
	//p_vapor = QV(i,j,k_p) + QBAR(i,j,k_p);// + qb[k_p];	// full vapor is base state plus perturbation
	//p_temp = TH(i,j,k_p) + THBAR(i,j,k_p) ;//+ tb[k_p];		
			
	//----------------------------------------------
	// Let the parcel ascend...
	//----------------------------------------------
	for(int k=k_p;k<NZ;k++){

		//----------------------------------------------
		// The ambient environment
		//----------------------------------------------
		theta = TH(i,j,k) + THBAR(i,j,k);// + tb[k];		// full potential temperature is base state plus perturbation		
		//pressure = PI(i,j,k)/(cp*tbv[k]) + PBAR(i,j,k);	// full pressure is base state plus perturbation			
		pressure = PI(i,j,k)/(cp*tb[k]) + PBAR(i,j,k);	// full pressure is base state plus perturbation
		temperature = theta*pressure;					// actual temperature
		vapor = QV(i,j,k_p) + QBAR(i,j,k);// + qb[k];
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
double get_CAPE_base(int i,int j,int k_p){
	
	double esl,qvsat,theta,vapor,pressure,pd,f,phi,vtemp,temperature;
	double p_temp,p_vapor,p_vtemp;
	double c,buoy,cape = 0;
	double cpRd = cp/Rd;

	//----------------------------------------------
	// Initialize parcel's temperature and
	// water vapor content
	//----------------------------------------------
	p_vapor = IQBAR(i,j,k_p) + qb[k_p];	// full vapor is base state plus perturbation
	p_temp = ITHBAR(i,j,k_p) + tb[k_p];		
			
	//----------------------------------------------
	// Let the parcel ascend...
	//----------------------------------------------
	for(int k=k_p;k<NZ;k++){

		//----------------------------------------------
		// The ambient environment
		//----------------------------------------------
		theta = ITHBAR(i,j,k) + tb[k];		// full potential temperature is base state plus perturbation		
		pressure = IPBAR(i,j,k);	// full pressure is base state plus perturbation			
		temperature = theta*pressure;					// actual temperature
		vapor = IQBAR(i,j,k) + qb[k];
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