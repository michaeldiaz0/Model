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

const double cpRd = cp/Rd;

enum functions {
	Constant,
	Linear,
	Parabolic
} Rainfall_Interpolation_Method;

double *qp_edges;
double *qn_edges;
double *delta_qk_mono;
double *dBZ;

double (*piecewise_interp)(double *zin,double *zout,double *qin,double *qout,int zlevs);

double piecewise_constant_interp( double *zin,double *zout,double *qin,double *qout,int zlevs);
double piecewise_linear_interp(   double *zin,double *zout,double *qin,double *qout,int zlevs);
double piecewise_parabolic_interp(double *zin,double *zout,double *qin,double *qout,int zlevs);

double integrate_linear(   double *zin,double *qin,int p,double zl,double zh);
double integrate_parabolic(double *zin,double *qin,int p,double zl,double zh,double *qp,double *qn);

void calculate_edge_values(double *zin, double *qin, double *qp_edges, double *qn_edges,int n);

/***********************************************************************************
* 
************************************************************************************/
void init_microphysics(int nx,int ny){
	
	if(MICROPHYSICS_OPTION==3){
		init_thompson_microphysics(rank);
		dBZ = (double *)calloc(nx*ny*NZ,sizeof(double));
	}

	//Rainfall_Interpolation_Method = Constant;
	//Rainfall_Interpolation_Method = Linear;
	Rainfall_Interpolation_Method = Parabolic;

	if(Rainfall_Interpolation_Method == Constant){
		
		piecewise_interp = &piecewise_constant_interp;
			
	} else if(Rainfall_Interpolation_Method == Linear ){
		
		piecewise_interp = &piecewise_linear_interp;
		
	} else if(Rainfall_Interpolation_Method == Parabolic){
		
		piecewise_interp = &piecewise_parabolic_interp;
		
		qp_edges = (double *)calloc(NZ,sizeof(double));
		qn_edges = (double *)calloc(NZ,sizeof(double));
		delta_qk_mono = (double *)calloc(NZ,sizeof(double));
	}
}

/*********************************************************************
* Return the sign of x multiplied by one
*
**********************************************************************/
inline static double signof(double x){
	
	if(x>0){ return 1;}
	
	return -1;
}

/*********************************************************************
* Return the maximum value of three numbers
*
**********************************************************************/
inline double max_of_three(double x,double y,double z){
	
	return fmax(x,fmax(y,z));
}

/*********************************************************************
* Return the minimum value of three numbers
*
**********************************************************************/
inline static double min_of_three(double x,double y,double z){
	
	return fmin(x,fmin(y,z));
}

/******************************************************
* 
*******************************************************/
void run_microphysics(int il,int ih,int jl,int jh){
	
	bool do_radar = false;

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
		// Based on Rutledge and Hobbs (1983), Hong et al. (200x)
		//-------------------------------------------------			
		case 3:

			

			if(OUTPUT_TO_FILE && bigcounter % outfilefreq == 0){
				do_radar = true;
			}

			if(PARALLEL){

				zero_moisture(il,ih,jl,jh,fNX*fNY*fNZ);

				run_thompson_microphysics(
				il,ih,jl,jh,2,fNZ-1,fNX,fNY,fNZ,
				qvps,qcps,qips,qrps,qsps,qgps,
				nips,nrps,thps,pis,wps,accRain,accSnow,
				do_radar
				);
				
			} else {

				zero_moisture(il,ih,jl,jh,NX*NY*NZ);

				run_thompson_microphysics(
				il,ih,jl,jh,2,NZ-1,NX,NY,NZ,
				qvps,qcps,qips,qrps,qsps,qgps,
				nips,nrps,thps,pis,wps,accRain,accSnow,
				do_radar
				);
			}
				
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
			// hydrometeor density of advected parcel
			//-----------------------------------------------------------		
			for(int k=0;k<NZ-1;k++){ QA_mass[k] = rhou[k]*var[INDEX(i,j,k)] * (ZW(k+1)-ZW(k)) / (ZA_edges[k+1]-ZA_edges[k]);}
			//-----------------------------------------------------------
			// Interpolation
			//-----------------------------------------------------------		
	        precip[INDEX2D(i,j)] += piecewise_interp(ZA_edges,&ZW(0),QA_mass,QD_mass,NZ);				
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
double integrate_linear(double *zin,double *qin,int p,double zl,double zh){

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
* Integrate a segment assuming a parabolic approximation
*
* zin - input x values
* qin - input y values
* p   - position in array
* zl  - height of lower bound of integration
* zh  - height of upper bound of integration
* qp - values on positive side of cell edges
* qn - values on negative side of cell edges
!*********************************************************************/
double integrate_parabolic(double *zin,double *qin,int p,double zl,double zh,double *qp,double *qn){

    double high,low,t0,t1;

    t0 = (zl-zin[p]) / (zin[p+1]-zin[p]);
    t1 = (zh-zin[p]) / (zin[p+1]-zin[p]);
    
    high = (qp[p]+qn[p]-2.0*qin[p])*t1*t1*t1 - (qp[p]+2.0*qn[p]-3.0*qin[p])*t1*t1 + qn[p]*t1;
    low  = (qp[p]+qn[p]-2.0*qin[p])*t0*t0*t0 - (qp[p]+2.0*qn[p]-3.0*qin[p])*t0*t0 + qn[p]*t0;

	return high - low;
}

/*********************************************************************
* 
**********************************************************************/
double below_surface_rainfall(double *zin, double *qin, int zlevs, int *highest_negative_level){

	double rainfall;

	for(int j=1;j<zlevs-2;j++){
		
		if(zin[j] < 0 && zin[j+1] < 0){
		
			rainfall += qin[j] * (zin[j+1]-zin[j]);
		
			*highest_negative_level = j+1;
		}
	}
	
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
double piecewise_constant_interp(double *zin,double *zout,double *qin,double *qout,int zlevs){

	double factor;

	memset(qout,0,zlevs*sizeof(double));

	int highest_negative_level = 1;

	double rainfall = below_surface_rainfall(zin, qin, zlevs, &highest_negative_level);
	//-----------------------------------------------------------
	// Rainfall that has fallen between the surface and the 
	// highest edge below the surface
	//-----------------------------------------------------------	
	if(zin[highest_negative_level]<0){
	
		rainfall += integrate_linear(zin,qin,highest_negative_level,zin[highest_negative_level],0) * 
			
			-zin[highest_negative_level] *  
				
				(
				(zin[highest_negative_level+1] - zin[highest_negative_level]) /
				 -zin[highest_negative_level]
				)
				;
	}

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
double piecewise_linear_interp(double *zin,double *zout,double *qin,double *qout,int zlevs){

	memset(qout,0,zlevs*sizeof(double));

	int highest_negative_level = 1;

	double rainfall = below_surface_rainfall(zin, qin, zlevs, &highest_negative_level);
	
	//-----------------------------------------------------------
	// Rainfall that has fallen between the surface and the 
	// highest edge below the surface
	//-----------------------------------------------------------
	if(zin[highest_negative_level]<0){
	
		rainfall += integrate_linear(zin,qin,highest_negative_level,zin[highest_negative_level],0) * 
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

			qout[i] += integrate_linear(zin,qin,j,zin[j],zout[i+1]) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}
		//-----------------------------------------------------------
		// Overlap (mA4 -> m3)
		//-----------------------------------------------------------
		if(zout[i] < zin[j+1] && zout[i+1] > zin[j+1] && zout[i] > zin[j]){

			qout[i] += integrate_linear(zin,qin,j,zout[i],zin[j+1]) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}
		//-----------------------------------------------------------
		// Output cell falls completely within input cell (mA5 -> m4)
		//-----------------------------------------------------------
		if(zout[i] >= zin[j] && zout[i+1] <= zin[j+1]){

			qout[i] += integrate_linear(zin,qin,j,zout[i],zout[i+1]) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}
		//-----------------------------------------------------------
		// Input cell falls completely within output cell
		//-----------------------------------------------------------		
		if(zout[i] < zin[j] && zout[i+1] > zin[j+1]){

			qout[i] += integrate_linear(zin,qin,j,zin[j],zin[j+1]) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}	
	}}

	return rainfall;
}


/*********************************************************************
* Piecewise parabolic interpolation which preserves total mass
*
* zin - input height levels
* zout - output height levels
* qin - input field
* qout - output field
* zlevs number of levels (i.e. array length for all input arrays)
* qp_edges - values on positive edge of cell
* qn_edges - values on negative edge of cell
**********************************************************************/
double piecewise_parabolic_interp(double *zin,double *zout,double *qin,double *qout,int zlevs){

    memset(qout,0,zlevs*sizeof(double));

	int highest_negative_level = 1;

	double rainfall = below_surface_rainfall(zin, qin, zlevs, &highest_negative_level);

	calculate_edge_values(zin,qin,qp_edges,qn_edges,zlevs);
	
	if(zin[highest_negative_level]<0){
	
		rainfall += integrate_parabolic(
			zin,
			qin,
			highest_negative_level,
			zin[highest_negative_level],
			0,
			qp_edges,qn_edges
			) * 
			-zin[highest_negative_level] *  ((zin[highest_negative_level+1] - zin[highest_negative_level]) /
				 -zin[highest_negative_level])
				;
	}

    //-----------------------------------------------------------
    // Compare input and output cells to look for matches.
    // Watch out for j = 0, the integrate_linear function 
    // looks for a j - 1 valueS
    //-----------------------------------------------------------
	for(int i=0;i<zlevs-1;i++){
	for(int j=0;j<zlevs-1;j++){
        //-----------------------------------------------------------
        // Overlap (mA4 -> m2)
        //-----------------------------------------------------------
		if(zout[i] < zin[j] && zout[i+1] < zin[j+1] && zout[i+1] > zin[j]){
            
            qout[i] = qout[i] + integrate_parabolic(zin,qin,j,zin[j],zout[i+1],qp_edges,qn_edges) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));          
        }
        //-----------------------------------------------------------
        // Overlap (mA4 -> m3)
        //-----------------------------------------------------------
        if(zout[i] < zin[j+1] && zout[i+1] > zin[j+1] && zout[i] > zin[j]){

            qout[i] = qout[i] + integrate_parabolic(zin,qin,j,zout[i],zin[j+1],qp_edges,qn_edges) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));           
        }
        //-----------------------------------------------------------
        // Output cell falls completely within input cell (mA5 -> m4)
        //-----------------------------------------------------------
        if(zout[i] >= zin[j] && zout[i+1] <= zin[j+1]){

            qout[i] = qout[i] + integrate_parabolic(zin,qin,j,zout[i],zout[i+1],qp_edges,qn_edges) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));           
        }
        //-----------------------------------------------------------
        // Input cell falls completely within output cell
        //-----------------------------------------------------------        
        if(zout[i] < zin[j] && zout[i+1] > zin[j+1]){

            qout[i] = qout[i] + integrate_parabolic(zin,qin,j,zin[j],zin[j+1],qp_edges,qn_edges) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));          
        }      
    }}
	
	return rainfall;
}

/*********************************************************************
* Calculate the values at the cell edges 
*
* zin - input x values
* qin - input y values
* qp_edges - values on positive edge of cell
* qn_edges - values on negative edge of cell
* n - length of arrays
!*********************************************************************/
void calculate_edge_values(double *zin, double *qin, double *qp_edges, double *qn_edges,int n){

    double edge_slope0, edge_slope1, qp_edge, qn_edge, mid_slope, t0, t1;
    double delta_qk, delta_qk_max, delta_qk_min;
    double q_min,q_max,q_mp,q_lc;

    //----------------------------------------------------
    // Monotonic difference
    //----------------------------------------------------
	for(int k=1;k<n-1;k++){

        delta_qk = 0.25 * (qin[k+1] - qin[k-1]);
        delta_qk_max = max_of_three(qin[k+1],qin[k],qin[k-1]) - qin[k];
        delta_qk_min = qin[k] - min_of_three(qin[k+1],qin[k],qin[k-1]);

        delta_qk_mono[k] = signof(delta_qk) * min_of_three(fabs(delta_qk),delta_qk_max,delta_qk_min);

	}

    delta_qk_mono[0] = delta_qk_mono[1];
	delta_qk_mono[n-1] = delta_qk_mono[n-2];

    //----------------------------------------------------
    // Calculate the value at the endpoints of the cell
    //----------------------------------------------------
	for(int k=1;k<n-1;k++){

        qn_edges[k] = ( qin[k-1] * (zin[k+1]-zin[k]) + qin[k] * (zin[k]-zin[k-1])) / 
			( zin[k+1]-zin[k-1]) - (delta_qk_mono[k] - delta_qk_mono[k-1]) / 3.0;
		
        qp_edges[k-1] = qn_edges[k];
	}

	// set boundary values
    qp_edges[n-1] = qin[n-1];
    qp_edges[0] = qin[0];
    qn_edges[n-1] = qin[n-1];
    qn_edges[0] = qin[0];

    //----------------------------------------------------
    // Ensure positive definiteness 
    //----------------------------------------------------
	for(int k=0;k<n;k++){
		qn_edges[k] = qin[k] - signof(delta_qk_mono[k]) * fmin(fabs(2.0*delta_qk_mono[k]),abs(qn_edges[k]-qin[k]));
		qp_edges[k] = qin[k] + signof(delta_qk_mono[k]) * fmin(fabs(2.0*delta_qk_mono[k]),abs(qp_edges[k]-qin[k]));
	}

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
	
	if(MICROPHYSICS_OPTION==3){
		switch_array(&qgs,&qgps);
		switch_array(&nis,&nips);
		switch_array(&nrs,&nrps);
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

	if(MICROPHYSICS_OPTION==3){
		memcpy(qgms,qgps,num_bytes);
		memcpy(nims,nips,num_bytes);
		memcpy(nrms,nrps,num_bytes);
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

	if(MICROPHYSICS_OPTION==3){

		for(int i=0;i<size;i++)
			if(qgps[i] < 0){ qgps[i] = 0;}

		for(int i=0;i<size;i++)
			if(nrps[i] < 0){ nrps[i] = 0;}

		for(int i=0;i<size;i++)
			if(nips[i] < 0){ nips[i] = 0;}
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
double get_dimensional_pressure(int index,int k){

#if !HYDROSTATIC
	return p0*pow(m_pbar[index],cpRd) + pis[index]*rhou[k];
#else
	return p0*pow(m_pbar[index]+pis[index],cpRd);
#endif
}

/*********************************************************************
* Calculate saturation mixing ratio for water, ice, or mixed phase
* depending on the temperature.
*
* pd -> dimensional pressure (Pa)
*
* for different phases (l->liquid, i->ice)
*
**********************************************************************/
double get_qvsat_mixed(double temperature,double pd){

	double esl,esi,qvsat,qvl_sat,qvi_sat;

	//------------------------------------------------
	// All rain
	//------------------------------------------------
	if(temperature > tmax){

		// calculate saturation mixing ratio
		esl = SAT_VAP_WAT(temperature);
		qvsat = SAT_MIX_RATIO(esl,pd);

	//------------------------------------------------
	// All ice
	//------------------------------------------------
	} else if(temperature < tmin) {
		// calculate saturation mixing ratio
		esl = SAT_VAP_ICE(temperature);
		qvsat = SAT_MIX_RATIO(esl,pd);

	//------------------------------------------------
	// Mixed phase
	//------------------------------------------------
	} else {

		// calculate saturation mixing ratio
		esl = SAT_VAP_WAT(temperature);
		esi = SAT_VAP_ICE(temperature);

		qvl_sat = SAT_MIX_RATIO(esl,pd);
		qvi_sat = SAT_MIX_RATIO(esi,pd);

		qvsat = ALPHA(temperature) * qvl_sat + (1.0-ALPHA(temperature)) * qvi_sat;

	}

	return qvsat;
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