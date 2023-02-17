#include "stdafx.h"
#include "pcomm.h"

/*******************************************************************************
* 
********************************************************************************/
#define SAT_VAP_WAT(t) ( 611.2 * exp(17.67 * (t-273.15) / (t - 29.65) ) )
#define SAT_VAP_ICE(t) ( 611.2 * exp(21.8745584 * (t-273.15) / (t - 7.66) ) )
#define SAT_MIX_RATIO(e,p) ( 0.62197 * e / (p-e) )

#if HYDROSTATIC
    #define CONVERT_PRESSURE(i,j,k) PI(i,j,k)
#else
    #define CONVERT_PRESSURE(i,j,k) PI(i,j,k)/(cp*tbv[k])
#endif

double piecewise_interp(double *zin,double *zout,double *qin,double *qout,int zlevs,int,int);


/***********************************************************************************
* 
************************************************************************************/
void init_microphysics(int nx,int ny){

	//-------------------------------------------------------------------
	// If initializing from an output file, store the precipitation total
	//-------------------------------------------------------------------
	if(isRestartRun || (PERTURBATION_OPTION == 0 && perturbationFileTime > 0)){
		
		for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
		
			//raintotal[ny*i+j] = FRICTION(i,j,NZ-2);
			
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
* Accumulated precipitation (mm)
*
* il,ih,jl,jh - array index bounds
* vel - 		fall speed
* hydro_field	hydrometeor field
* output -		accumulation precipitation (mm)
*
**********************************************************************/
void precip_rate(int il,int ih,int jl,int jh,double *vel,double *hydro_field,double *output){
	
	int ny;
	
	if(PARALLEL){ ny = fNY;}
	else { ny = NY;}
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
			
		output[INDEX2D(i,j)] += 
			rhow[1] * 
			0.5 * ( vel[INDEX(i,j,0)] + vel[INDEX(i,j,1)] ) * 
			0.5 * ( hydro_field[INDEX(i,j,0)] + hydro_field[INDEX(i,j,1)] ) *
			dt;
	}}
	
}

/****************************************************
* Calculate eulerian fall speed of rain
*
* vt -> terminal fall speed of rain (output)
* qr -> rain water mixing ratio (input)
*
*****************************************************/
void calculate_eulerian_fall_speed_rain(double *vtr, double *qr, int il,int ih,int jl,int jh){
	
	if(MICROPHYSICS_OPTION==1){
	
		calculate_rain_fall_speed_kessler( vtr, qr, il, ih, jl, jh);
		
	} else if(MICROPHYSICS_OPTION==2){
		
		calculate_rain_fall_speed_rutledge(vtr, qr, il, ih, jl, jh);
	}

}

/****************************************************
* Calculate eulerian fall speed of hydrometeors
*
* vt -> terminal fall speed of rain (output)
* qr -> rain water mixing ratio (input)
*
*****************************************************/
void calculate_eulerian_fall_speed_snow_ice(double *vts, double *qs, double *vti, double *qi, int il,int ih,int jl,int jh){
	
	if(MICROPHYSICS_OPTION==2){
	
		calculate_ice_fall_speed_rutledge( vti, qi, il, ih, jl, jh);
		
		calculate_snow_fall_speed_rutledge(vts, qs, il, ih, jl, jh);
	}

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
void hydrometeor_fallout(double *var,double *vel,int il,int ih,int jl,int jh,double *precip){
	
	const double DeCFL0 = 0.05;	// tunable paramter to maintain numerical stability
	const double minVar = 1e-12;
	
	double ZA_edges[NZ];
	double QA_mass[NZ];
	double QD_mass[NZ];
	double Vel_edges[NZ];
	double DeCFL;
	
	bool hasRain;
	
	int ny;
	if(PARALLEL){ ny = fNY;}
	else { ny = NY;}
	
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
			for(int k=0;k<NZ-1;k++){ ZA_edges[k] = ZW(k) - dt * Vel_edges[k];}	
			//-----------------------------------------------------------
			// hydrometeor mass of advected parcel
			//-----------------------------------------------------------		
			for(int k=0;k<NZ-1;k++){ QA_mass[k] = rhou[k]*var[INDEX(i,j,k)] * (ZW(k+1)-ZW(k)) / (ZA_edges[k+1]-ZA_edges[k]);}
			//-----------------------------------------------------------
			// Interpolation
			//-----------------------------------------------------------
			precip[INDEX2D(i,j)] += piecewise_interp(ZA_edges,&ZW(0),QA_mass,QD_mass,NZ,i,j);
			
			//-----------------------------------------------------------
			// Divide by density for final mass
			//-----------------------------------------------------------
			for(int k=0;k<NZ;k++){ var[INDEX(i,j,k)] = QD_mass[k] * one_d_rhou[k];}
		
		}
	}}
	
}

/*********************************************************************
* Integrate linear function
**********************************************************************/
double integrate_linear(double *zin,double *qin,int p,double zl,double zh,bool debug){

	double edge_slope0 = 2.0 * (qin[p  ]-qin[p-1]) / (zin[p+1] - zin[p-1]);
	double edge_slope1 = 2.0 * (qin[p+2]-qin[p+1]) / (zin[p+2] - zin[p  ]);
	
	double mid_slope;
	
	if( (edge_slope0<0 && edge_slope1>0) || (edge_slope0>0 && edge_slope1<0) ){	

		mid_slope = 0;

	} else{

		mid_slope = 0.5*(edge_slope0 + edge_slope1);
	}

	double qp_edge = 0.5 * mid_slope * (zin[p+1]-zin[p]) + qin[p];
	double qn_edge = 2.0 * qin[p] - qp_edge;
	
	if(qp_edge < 0 || qn_edge < 0){
		
		qp_edge = qin[p];
		qn_edge = qin[p];
	}

	double t0 = (zl-zin[p]) / (zin[p+1]-zin[p]);
	double t1 = (zh-zin[p]) / (zin[p+1]-zin[p]);
	
	return (qp_edge - qn_edge) * 0.5*t1*t1 + qn_edge * t1 - ((qp_edge - qn_edge) * 0.5*t0*t0 + qn_edge * t0);
}

/*********************************************************************
* Piecewise linear interpolation which preserves total mass
*
* zin - input height levels
* zout - output height levels
* qin - input field
* qout - output field
* zlevs number of levels (i.e. array length for all input arrays)
*
* @return rainfall in kg / m^2
**********************************************************************/
double piecewise_interp(double *zin,double *zout,double *qin,double *qout,int zlevs,int ic,int jc){

	double factor;

	double rainfall = 0;

	for(int i=0;i<zlevs;i++){ qout[i] = 0;}

	int highest_negative_level = 1;

	for(int j=1;j<zlevs-2;j++){
		
		int i = 0;
		
		if(zin[j] < 0 && zin[j+1] < 0){
			
			rainfall += qin[j] * (zin[j+1]-zin[j]);
			
			highest_negative_level = j+1;
			
		}
	}
	
	if(zin[highest_negative_level]<0){
	
		rainfall += integrate_linear(zin,qin,highest_negative_level,zin[highest_negative_level],0,false) * 
			-zin[highest_negative_level] *  ((zin[highest_negative_level+1] - zin[highest_negative_level]) /
				 -zin[highest_negative_level])
				;
	}
	
	//-----------------------------------------------------------
	// Compare input and output cells to look for matches.
	// Watch out for j = 0, the integrate_linear function 
	// looks for a j - 1 value
	//-----------------------------------------------------------
	for(int i=1;i<zlevs-2;i++){
	for(int j=highest_negative_level;j<zlevs-2;j++){
		//-----------------------------------------------------------
		// Overlap (mA4 -> m2)
		//-----------------------------------------------------------
		if(zout[i] < zin[j] && zout[i+1] < zin[j+1] && zout[i+1] > zin[j]){

			qout[i] += integrate_linear(zin,qin,j,zin[j],zout[i+1],false) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}
		//-----------------------------------------------------------
		// Overlap (mA4 -> m3)
		//-----------------------------------------------------------
		if(zout[i] < zin[j+1] && zout[i+1] > zin[j+1] && zout[i] > zin[j]){

			qout[i] += integrate_linear(zin,qin,j,zout[i],zin[j+1],false) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}
		//-----------------------------------------------------------
		// Output cell falls completely within input cell (mA5 -> m4)
		//-----------------------------------------------------------
		if(zout[i] >= zin[j] && zout[i+1] <= zin[j+1]){

			qout[i] += integrate_linear(zin,qin,j,zout[i],zout[i+1],false) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}
		//-----------------------------------------------------------
		// Input cell falls completely within output cell
		//-----------------------------------------------------------		
		if(zout[i] < zin[j] && zout[i+1] > zin[j+1]){

			qout[i] += integrate_linear(zin,qin,j,zin[j],zin[j+1],false) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}	
	}}


	return rainfall;
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
double piecewise_constant_interp(double *zin,double *zout,double *qin,double *qout,int zlevs,int ic,int jc){

	double factor;

	int test_i = 13;
	int test_j = 24;


	double rainfall = 0;

	for(int i=0;i<zlevs;i++){ qout[i] = 0;}

	for(int i=1;i<zlevs-2;i++){
	for(int j=1;j<zlevs-2;j++){
		
		if(zout[i] >= zin[j] && zout[i+1] <= zin[j+1]){
			
			factor = 1.0;
			qout[i] += factor*qin[j];		
		}

		if(zout[i] < zin[j] && zout[i+1] < zin[j+1] && zout[i+1] > zin[j]){
			
			factor = (zout[i+1] - zin[j]) / (zout[i+1]-zout[i]);
			qout[i] += factor*qin[j];
		}
		
		if(zout[i] < zin[j+1] && zout[i+1] > zin[j+1] && zout[i] > zin[j]){
			
			factor = (zin[j+1] - zout[i]) / (zout[i+1]-zout[i]);
			qout[i] += factor*qin[j];
		}
	}}

	return rainfall;
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
		pressure = CONVERT_PRESSURE(i,j,k) + PBAR(i,j,k);	// full pressure is base state plus perturbation
		//pressure = PI(i,j,k)/(cp*tb[k]) + PBAR(i,j,k);	// full pressure is base state plus perturbation
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
#if 0

module piece_interpolation

    implicit none

    contains

    !*********************************************************************
    ! Calculate the values at the cell edges 
    !
    ! zin - input x values
    ! qin - input y values
    ! zl  - height of lower bound of integration
    ! zh  - height of upper bound of integration
    ! qp_edges - values on positive edge of cell
    ! qn_edges - values on negative edge of cell
    ! n - length of arrays
    !*********************************************************************/
    subroutine calculate_edge_values(zin,qin,qp_edges,qn_edges,n)

        real,dimension(:),intent(in) :: zin,qin
        real,dimension(:),intent(out) :: qp_edges,qn_edges
        integer,intent(in) :: n

        real :: edge_slope0, edge_slope1, qp_edge, qn_edge, mid_slope, t0, t1
        real :: delta_qk, delta_qk_max, delta_qk_min
        real,dimension(:),allocatable :: delta_qk_mono
        integer :: k 

        real :: q_min,q_max,q_mp,q_lc

        allocate(delta_qk_mono(n))

        !----------------------------------------------------
        ! Monotonic difference
        !----------------------------------------------------
        do k = 2,n-1

            delta_qk = 0.25 * (qin(k+1) - qin(k-1))
            delta_qk_max = max(qin(k+1),qin(k),qin(k-1)) - qin(k)
            delta_qk_min = qin(k) - min(qin(k+1),qin(k),qin(k-1))

            delta_qk_mono(k) = sign(min(abs(delta_qk),delta_qk_max,delta_qk_min),delta_qk)

        end do

        delta_qk_mono(1) = delta_qk_mono(2)
        delta_qk_mono(n) = delta_qk_mono(n)

        !----------------------------------------------------
        ! Calculate the value at the endpoints of the cell
        !----------------------------------------------------
        do k=2,n-1

            qn_edges(k) = ( qin(k-1) * (zin(k+1)-zin(k)) + qin(k) * (zin(k)-zin(k-1))) / ( zin(k+1)-zin(k-1)) - (delta_qk_mono(k) - delta_qk_mono(k-1)) / 3.0
            qp_edges(k-1) = qn_edges(k)          

        end do

        qp_edges(n) = qin(n) ! set boundary values
        qp_edges(1) = qin(1)
        qn_edges(n) = qin(n)
        qn_edges(1) = qin(1)

        !----------------------------------------------------
        ! Ensure positive definiteness 
        !----------------------------------------------------
        do k = 1,n

            qn_edges(k) = qin(k) - sign(min(abs(2.0*delta_qk_mono(k)),abs(qn_edges(k)-qin(k))),delta_qk_mono(k))
            qp_edges(k) = qin(k) + sign(min(abs(2.0*delta_qk_mono(k)),abs(qp_edges(k)-qin(k))),delta_qk_mono(k))

        end do 

        deallocate(delta_qk_mono)

    end subroutine calculate_edge_values 

    !*********************************************************************
    ! Integrate a segment assuming a parabolic approximation
    !
    ! zin - input x values
    ! qin - input y values
    ! p   - position in array
    ! zl  - height of lower bound of integration
    ! zh  - height of upper bound of integration
    ! qp - values on positive side of cell edges
    ! qn - values on negative side of cell edges
    !*********************************************************************/
    real function integrate_parabolic(qin,qp,qn,zin,zl,zh,p) result(value)

        real,intent(in) :: zl,zh
        real,dimension(:), intent(in) :: zin,qin,qp,qn
        integer,intent(in) :: p

        real :: high,low,t0,t1

        t0 = (zl-zin(p)) / (zin(p+1)-zin(p))
        t1 = (zh-zin(p)) / (zin(p+1)-zin(p))
        
        high = (qp(p)+qn(p)-2.0*qin(p))*t1*t1*t1 - (qp(p)+2.0*qn(p)-3.0*qin(p))*t1*t1 + qn(p)*t1
        low  = (qp(p)+qn(p)-2.0*qin(p))*t0*t0*t0 - (qp(p)+2.0*qn(p)-3.0*qin(p))*t0*t0 + qn(p)*t0

        value = high - low
        
    end function integrate_parabolic
    
    !*********************************************************************
    ! Integrate a segment using a linear approximation
    !
    ! zin       - input full levels
    ! qin       - input values to be interpolated at half levels
    ! zin_half  - input half levels
    ! p         - position in array
    ! zl        - height of lower bound of integration
    ! zh        - height of upper bound of integration
    !*********************************************************************/
    real function integrate_linear(zin,qin,p,zl,zh,zin_half) result(value)

        real,intent(in) :: zl,zh
        real,dimension(:), intent(in) :: zin,qin,zin_half
        integer,intent(in) :: p

        real :: edge_slope0,edge_slope1,qp_edge,qn_edge,mid_slope,t0,t1

        ! calculate the slopes at the endpoints of the cell
        edge_slope0 = (qin(p  )-qin(p-1)) / (zin_half(p ) - zin_half(p-1))
        edge_slope1 = (qin(p+1)-qin(p  )) / (zin_half(p+1) - zin_half(p  ))
       
        !print *,zl,zh
 
        
        if( (edge_slope0 .lt. 0 .and. edge_slope1 .gt. 0) .or. (edge_slope0 .gt. 0 .and. edge_slope1 .lt. 0) ) then
            mid_slope = 0
        else
            mid_slope = 0.5*(edge_slope0 + edge_slope1)
        end if

        qp_edge = 0.5 * mid_slope * (zin(p+1)-zin(p)) + qin(p)
        qn_edge = 2.0 * qin(p) - qp_edge
        
        !print *,qp_edge,qn_edge

        if(qp_edge .lt. 0 .or. qn_edge .lt. 0) then
            
            qp_edge = qin(p)
            qn_edge = qin(p)
            !printf("zero edge\n");
        end if

        t0 = (zl-zin(p)) / (zin(p+1)-zin(p))
        t1 = (zh-zin(p)) / (zin(p+1)-zin(p))
        !print *,t0,t1
        !print *,"" 
        !printf("integrate over %f to %f\n",t0,t1);

        !if(debug){
            !printf("midpoint = %f %f %f\n",1000*((qp_edge - qn_edge) * 0.5 + qn_edge),1000*qp_edge,1000*qn_edge);
        !}

        value = (qp_edge - qn_edge) * 0.5*t1*t1 + qn_edge * t1 - ((qp_edge - qn_edge) * 0.5*t0*t0 + qn_edge * t0)
        
    end function integrate_linear
    
    !*********************************************************************
    ! Piecewise linear interpolation which preserves total mass
    !
    ! zin - input height levels
    ! zout - output height levels
    ! qin - input field
    ! qout - output field
    ! zlevs number of levels (i.e. array length for all input arrays)
    !
    !*********************************************************************/
    subroutine piecewise_interp(zin,zout,zin_half,qin,qout,zlevs)

        real,dimension(:),intent(in) :: zin,zout,qin,zin_half
        real,dimension(:),intent(out) :: qout
        integer, intent(in) :: zlevs

        integer :: i,j

        qout = 0
 
        !-----------------------------------------------------------
        ! Compare input and output cells to look for matches.
        ! Watch out for j = 0, the integrate_linear function 
        ! looks for a j - 1 valueS
        !-----------------------------------------------------------
        do i=2,zlevs-2
        do j=2,zlevs-2
        
            !-----------------------------------------------------------
            ! Overlap (mA4 -> m2)
            !-----------------------------------------------------------
            if(zout(i) .lt. zin(j) .and. zout(i+1) .lt. zin(j+1) .and. zout(i+1) .gt. zin(j)) then
                
                qout(i) = qout(i) + integrate_linear(zin,qin,j,zin(j),zout(i+1),zin_half) * ((zin(j+1) - zin(j)) / (zout(i+1)-zout(i)))
                
            end if
            !-----------------------------------------------------------
            ! Overlap (mA4 -> m3)
            !-----------------------------------------------------------
            if(zout(i) .lt. zin(j+1) .and. zout(i+1) .gt. zin(j+1) .and. zout(i) .gt. zin(j)) then

                qout(i) = qout(i) + integrate_linear(zin,qin,j,zout(i),zin(j+1),zin_half) * ((zin(j+1) - zin(j)) / (zout(i+1)-zout(i)))
                
            end if
            !-----------------------------------------------------------
            ! Output cell falls completely within input cell (mA5 -> m4)
            !-----------------------------------------------------------
            if(zout(i) .ge. zin(j) .and. zout(i+1) .le. zin(j+1)) then

                qout(i) = qout(i) + integrate_linear(zin,qin,j,zout(i),zout(i+1),zin_half) * ((zin(j+1) - zin(j)) / (zout(i+1)-zout(i)))
                
            end if
            !-----------------------------------------------------------
            ! Input cell falls completely within output cell
            !-----------------------------------------------------------        
            if(zout(i) .lt. zin(j) .and. zout(i+1) .gt. zin(j+1)) then

                qout(i) = qout(i) + integrate_linear(zin,qin,j,zin(j),zin(j+1),zin_half) * ((zin(j+1) - zin(j)) / (zout(i+1)-zout(i)))
                
            end if
            
        end do
        end do

        

    end subroutine piecewise_interp
    
    !*********************************************************************
    ! Piecewise parabolic interpolation which preserves total mass
    !
    ! zin - input height levels
    ! zout - output height levels
    ! qin - input field
    ! qout - output field
    ! zlevs number of levels (i.e. array length for all input arrays)
    ! qp_edges - values on positive edge of cell
    ! qn_edges - values on negative edge of cell
    !*********************************************************************/
    subroutine piecewise_parabolic_interp(zin,zout,qp_edges,qn_edges,qin,qout,zlevs)

        real,dimension(:),intent(in) :: zin,zout,qin,qp_edges,qn_edges
        real,dimension(:),intent(out) :: qout
        integer, intent(in) :: zlevs

        integer :: i,j

        qout = 0

        !-----------------------------------------------------------
        ! Compare input and output cells to look for matches.
        ! Watch out for j = 0, the integrate_linear function 
        ! looks for a j - 1 valueS
        !-----------------------------------------------------------
        do i=1,zlevs-2
        do j=1,zlevs-2
        
            !-----------------------------------------------------------
            ! Overlap (mA4 -> m2)
            !-----------------------------------------------------------
            if(zout(i) .lt. zin(j) .and. zout(i+1) .lt. zin(j+1) .and. zout(i+1) .gt. zin(j)) then
                
                qout(i) = qout(i) + integrate_parabolic(qin,qp_edges,qn_edges,zin,zin(j),zout(i+1),j) * ((zin(j+1) - zin(j)) / (zout(i+1)-zout(i)))
                
            end if
            !-----------------------------------------------------------
            ! Overlap (mA4 -> m3)
            !-----------------------------------------------------------
            if(zout(i) .lt. zin(j+1) .and. zout(i+1) .gt. zin(j+1) .and. zout(i) .gt. zin(j)) then

                qout(i) = qout(i) + integrate_parabolic(qin,qp_edges,qn_edges,zin,zout(i),zin(j+1),j) * ((zin(j+1) - zin(j)) / (zout(i+1)-zout(i)))
                
            end if
            !-----------------------------------------------------------
            ! Output cell falls completely within input cell (mA5 -> m4)
            !-----------------------------------------------------------
            if(zout(i) .ge. zin(j) .and. zout(i+1) .le. zin(j+1)) then

                qout(i) = qout(i) + integrate_parabolic(qin,qp_edges,qn_edges,zin,zout(i),zout(i+1),j) * ((zin(j+1) - zin(j)) / (zout(i+1)-zout(i)))
                
            end if
            !-----------------------------------------------------------
            ! Input cell falls completely within output cell
            !-----------------------------------------------------------        
            if(zout(i) .lt. zin(j) .and. zout(i+1) .gt. zin(j+1)) then

                qout(i) = qout(i) + integrate_parabolic(qin,qp_edges,qn_edges,zin,zin(j),zin(j+1),j) * ((zin(j+1) - zin(j)) / (zout(i+1)-zout(i)))
                
            end if
            
        end do
        end do

    end subroutine piecewise_parabolic_interp

end module

program test

    use piece_interpolation
    use mod_adaptive_interp    

    implicit none

    integer :: a,i,n,p
    real :: b,c,d,e,f,g,h,zl,zh
    
    real,dimension(:),allocatable :: z_in,z_out,q_in,q_out,rho,z_in_f,z_out_f,qp_edges,qn_edges

    n = 88

    allocate(z_in(n))
    allocate(z_out(n))
    allocate(q_in(n))
    allocate(q_out(n))
    allocate(rho(n))
    allocate(z_in_f(n))
    allocate(z_out_f(n))
    allocate(qp_edges(n))
    allocate(qn_edges(n))

    open(unit=1,file="values.txt")
    open(unit=2,file="height_levs.txt")

    do i = 1,n
        read(1,*) a,z_out(i),z_in(i),d,e,f,q_in(i),rho(i)
    end do

    do i =1,n
        read(2,*) z_out_f(i),d
    end do

    close(1)
    close(2)    

    z_in_f(1) = 0

    do i = 2,n
        z_in_f(i) = 0.5*(z_in(i)+z_in(i-1))
    end do

    p = 28

    do i = 1,n
!        print *,z_in_f(i),z_out_f(i)
    end do

    !call piecewise_interp(z_in_f,z_out_f,z_in,q_in,q_out,n)

    do i = 1,1!5760    

        call adaptiveInterpolation1D(z_in,q_in,z_out,q_out,n,n,4,2)
!        call adaptiveInterpolation1D(z_out,q_out,z_in,q_in,n,n,4,2)

        !call piecewise_interp(z_in_f,z_out_f,z_in,q_in,q_out,n)
        !call piecewise_interp(z_out_f,z_in_f,z_out,q_out,q_in,n)
        
        !call calculate_edge_values(z_in_f,q_in,qp_edges,qn_edges,n)
        !call piecewise_parabolic_interp(z_in_f,z_out_f,qp_edges,qn_edges,q_in,q_out,n)

     !   call calculate_edge_values(z_out_f,q_out,qp_edges,qn_edges,n)
     !   call piecewise_parabolic_interp(z_out_f,z_in_f,qp_edges,qn_edges,q_out,q_in,n)
    end do


    do i = 1,n-1
        print *,i,z_in(i),z_out(i),q_in(i),q_out(i)
    end do


    zh = 0.5*(z_in(p)+z_in(p+1))
    zl = 0.5*(z_in(p-1)+z_in(p))

    b = integrate_parabolic(q_in,qp_edges,qn_edges,z_in_f,zl,zh,p) 
    !print *,zl,zh,b

    b = integrate_linear(z_in_f,q_in,p,zl,zh,z_in)
    !print *,zl,zh,b

end program


            !q_mp = qin(k) - 2.0*delta_qk_mono(k)
            !q_lc = qin(k) + 1.5*(delta_qk_mono(k+2)-delta_qk_mono(k)) - delta_qk_mono(k)

            !q_min = min(qin(k),q_mp,q_lc)
            !q_max = max(qin(k),q_mp,q_lc)

            !qn_edges(k) = min(max(qn_edges(k),q_min),q_max)

            !-----------------------------------------------------

            !q_mp = qin(k) + 2.0*delta_qk_mono(k)
            !q_lc = qin(k) + 1.5*(delta_qk_mono(k)-delta_qk_mono(k-2)) + delta_qk_mono(k)

            !q_min = min(qin(k),q_mp,q_lc)
            !q_max = max(qin(k),q_mp,q_lc)

            !qp_edges(k) = min(max(qp_edges(k),q_min),q_max)
#endif
