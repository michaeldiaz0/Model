#include "stdafx.h"
#include "advection.h"
#include "fluxes.h"
#include "pcomm.h"

/******************************************************************************************
*
* 							ADVECTION SUBROUTINES
*
*******************************************************************************************/

/***************************************************************************
* ------------------------------ MACROS ------------------------------------
****************************************************************************/
// Model base state advection
#define Z_ADVECT_BASE(s) ( 0.50*( rhow[k+1]*W(i,j,k+1)*(s[k+1]-s[k  ]) \
							    + rhow[k  ]*W(i,j,k  )*(s[k  ]-s[k-1])	\
							    ) * one_d_rhou[k] )

#define LOOP3D_IJK(il,ih,jl,jh,kl,kh,equation) 	for(int i=il;i<ih;i++){ \
												for(int j=jl;j<jh;j++){ \
												for(int k=kl;k<kh;k++){ \
													equation;	\
												}}}

#define LOOP2D_JK(jl,jh,kl,kh,equation) 		for(int j=jl;j<jh;j++){ \
												for(int k=kl;k<kh;k++){ \
													equation;	\
												}}

/*****************************************************************************************
* 
******************************************************************************************/
inline double buoyancy_dry(int i,int j,int k){
	return TH(i,j,k) / tbv[k];
}

/*****************************************************************************************
* 
******************************************************************************************/
inline double buoyancy_warm_micro(int i,int j,int k){
	return TH(i,j,k) / tbv[k] + 0.61*QV(i,j,k) - QC(i,j,k) - QR(i,j,k);
}

/*****************************************************************************************
* 
******************************************************************************************/
inline double buoyancy_cold_micro(int i,int j,int k){
	return TH(i,j,k) / tbv[k] + 0.61*QV(i,j,k) - QC(i,j,k) - QR(i,j,k) - QI(i,j,k) - QS(i,j,k);
}

/*****************************************************************************************
* 
******************************************************************************************/
inline double buoyancy_cold_micro_graupel(int i,int j,int k){
	return TH(i,j,k) / tbv[k] + 0.61*QV(i,j,k) - QC(i,j,k) - QR(i,j,k) - QI(i,j,k) - QS(i,j,k) - QG(i,j,k);
}

/******************************************************************************************
*
* FUNCTIONS TO CALCULATE ADVECTIVE TENDENCIES FOR A SINGLE POINT.
* MUST CALL A VERTICAL AND HORIZONTAL INTERPOLATION ROUTINE FIRST.
*
*******************************************************************************************/

/*****************************************************************************************
* Calculate zonal velocity tendencies using the flux form of the advection equation.
*
* @param i,j,k - point in array
* @return zontal momentum tendency
******************************************************************************************/
inline static double u_tend_cell(int i,int j,int k){

	return (
		//---------------------------------------------------------------------
		// advection of zonal momentum by zonal wind
		- dtx*( UCELL(j,k).east - UCELL(j,k).west )
		//---------------------------------------------------------------------
		// advection of zonal momentum by meridional wind
		- dty*( UCELL(j,k).north - UCELL(j-1,k).north )
		//---------------------------------------------------------------------
		// advection of zonal momentum by vertical wind
		- DTZ(k)*( rhow[k+1] * UCELL(j,k).top - rhow[k] * UCELL(j,k-1).top ) * one_d_rhou[k]
		//---------------------------------------------------------------------
		// coriolis
		+ dt* FC(j) * UCELL(j,k).other_vel
		//---------------------------------------------------------------------
		// friction
		#if USE_LINEAR_FRICTION
		- dt * FRICTION(i,j,k) * U(i,j,k)
		#endif
	);
}

/*****************************************************************************************
* Calculate meridional velocity tendencies using the flux form of the advection equation.
*
* @param i,j,k - point in array
******************************************************************************************/
inline static double v_tend_cell(int i,int j,int k){

	return (
		//---------------------------------------------------------------------
		// advection of meridional momentum by zonal wind
		- dtx*( VCELL(j,k).east - VCELL(j,k).west)
		//---------------------------------------------------------------------
		// advection of meridional momentum by meridional wind
		- dty*( VCELL(j,k).north - VCELL(j-1,k).north)
		//---------------------------------------------------------------------
		// advection of meridional momentum by vertical wind
		- DTZ(k)*( rhow[k+1] * VCELL(j,k).top - rhow[k] * VCELL(j,k-1).top) * one_d_rhou[k] 
		//---------------------------------------------------------------------
		// coriolis
		- dt * FC(j) * VCELL(j,k).other_vel
		//---------------------------------------------------------------------
		// friction
		#if USE_LINEAR_FRICTION
		- dt * FRICTION(i,j,k) * V(i,j,k)
		#endif
	);
}

/*****************************************************************************************
* Calculate vertical velocity tendencies using the flux form of the advection equation.
* 
* @param i,j,k - point in array
******************************************************************************************/
inline static double w_tend_cell(int i,int j,int k){

	return (
		//---------------------------------------------------------------------
		// advection of vertical momentum by zonal wind
		- dtx*( WCELL(j,k).east - WCELL(j,k).west)
		//---------------------------------------------------------------------
		// advection of vertical momentum by meridional wind
		- dty*( WCELL(j,k).north - WCELL(j-1,k).north)
		//---------------------------------------------------------------------
		// advection of vertical momentum by vertical wind
		- DTZW(k)*( rhou[k] * WCELL(j,k).top - rhou[k-1] * WCELL(j,k-1).top) * one_d_rhow[k]
		//---------------------------------------------------------------------
		// bouyancy
		//+ dt*grav*0.5*( BUOYANCY(i,j,k)+BUOYANCY(i,j,k-1) )
	);
}

/*****************************************************************************************
* Calculate potential temperature tendencies using the flux form of the advection equation.
* 
* @param i,j,k - point in array
******************************************************************************************/
inline static double theta_tend_cell(int i,int j,int k){

	return (
		//---------------------------------------------------------------------
		// advection of temperature by zonal wind
		- dtx * (THCELL(j,k).east - THCELL(j,k).west)
		//---------------------------------------------------------------------
		// advection of temperature by meridional wind	
		- dty * (THCELL(j,k).north - THCELL(j-1,k).north)
		//---------------------------------------------------------------------
		// advection of temperature by vertical wind			
		- DTZ(k)*(rhow[k+1]*THCELL(j,k).top - rhow[k]*THCELL(j,k-1).top) * one_d_rhou[k] 
		//---------------------------------------------------------------------
		// advection of model base state temperature by vertical wind
		- DTZ(k) * Z_ADVECT_BASE(tb)
	);
}

/*****************************************************************************************
* Calculate scalar tendencies using the flux form of the advection equation.
* 
* @param i,j,k - point in array
******************************************************************************************/
inline static double scalar_tend_cell(int i,int j,int k,struct cell *scell){

	return (
		//---------------------------------------------------------------------
		// advection of scalar by zonal wind
		- dtx * (SCELL(j,k).east - SCELL(j,k).west)
		//---------------------------------------------------------------------
		// advection of scalar by meridional wind	
		- dty * (SCELL(j,k).north - SCELL(j-1,k).north)
		//---------------------------------------------------------------------
		// advection of scalar by vertical wind			
		- DTZ(k)*(rhow[k+1]*SCELL(j,k).top - rhow[k]*SCELL(j,k-1).top) * one_d_rhou[k]
		//- DTZ(k)*(SCELL(j,k).top - SCELL(j,k-1).top)

	);
}

/*****************************************************************************************
* Calculate water vapor tendencies using the flux form of the advection equation.
* 
* @param i,j,k - point in array
******************************************************************************************/
inline static double qv_tend_cell(int i,int j,int k){

	return (
		//---------------------------------------------------------------------
		// advection by zonal wind
		- dtx * (QVCELL(j,k).east - QVCELL(j,k).west)
		//---------------------------------------------------------------------
		// advection by meridional wind	
		- dty * (QVCELL(j,k).north - QVCELL(j-1,k).north)
		//---------------------------------------------------------------------
		// advection by vertical wind			
		- DTZ(k) * (rhow[k+1]*QVCELL(j,k).top - rhow[k]*QVCELL(j,k-1).top) * one_d_rhou[k] 
		//---------------------------------------------------------------------
		// advection of model base state by vertical wind
		- DTZ(k) * Z_ADVECT_BASE(qb)
	);
}

/*****************************************************************************************
* Calculate cloud water tendencies using the flux form of the advection equation.
* 
* @param i,j,k - point in array
******************************************************************************************/
inline static double qc_tend_cell(int i,int j,int k){

	return (
		//---------------------------------------------------------------------
		// advection by zonal wind
		- dtx * (QCCELL(j,k).east - QCCELL(j,k).west)
		//---------------------------------------------------------------------
		// advection by meridional wind	
		- dty * (QCCELL(j,k).north - QCCELL(j-1,k).north)
		//---------------------------------------------------------------------
		// advection by vertical wind			
		- DTZ(k) * (rhow[k+1]*QCCELL(j,k).top - rhow[k]*QCCELL(j,k-1).top) * one_d_rhou[k]
	);
}

/*****************************************************************************************
* Calculate rain water tendencies using the flux form of the advection equation..
* 
* @param i,j,k - point in array
******************************************************************************************/
inline static double qr_tend_cell(int i,int j,int k){

	return (
		//---------------------------------------------------------------------
		// advection by zonal wind
		- dtx * (QRCELL(j,k).east - QRCELL(j,k).west)
		//---------------------------------------------------------------------
		// advection by meridional wind	
		- dty * (QRCELL(j,k).north - QRCELL(j-1,k).north)
		//---------------------------------------------------------------------
		// advection by vertical wind			
		- DTZ(k) * (rhow[k+1]*QRCELL(j,k).top - rhow[k]*QRCELL(j,k-1).top) * one_d_rhou[k]
	);
}

/******************************************************************************************
*
* SUBROUTINES TO LOOP THROUGH ALL POINTS IN DOMAIN AND CALCULATE ADVECTIVE TENDENCIES
*
*******************************************************************************************/

/*********************************************************************
* Loop through velocity arrays, interpolate from cell centers to cell 
* faces, and then calculate rate of advection using the flux form of 
* the advection equation.
*
* @param step - fraction of total time step dt to advance
* @param il,ih,jl,jh - beginning and ending indices
**********************************************************************/
void advect_uv_velocity(double step,int il,int ih,int jl,int jh){

	// Calculate the zonal flux on the eastern side of innermost column of boundary points.
	// This will become the flux on the western side of the leftmost interior (non-boundary) cell
	compute_west_uv(il-1,jl,jh);

	//-------------------------------------------------------
	// For each point within the east-west range
	//-------------------------------------------------------
	for(int i=il;i<ih;i++){

		// compute all fluxes for a YZ cross section
		compute_fluxes(i,jl-1,jh);

		//-------------------------------------------------------
		// For each point within the north-south range
		//-------------------------------------------------------
		for(int j=jl;j<jh;j++){
		for(int k=1;k<NZ-1;k++){

			// new values equal old values plus the rate of change
			UP(i,j,k) = UM(i,j,k) + step * u_tend_cell(i,j,k);
			VP(i,j,k) = VM(i,j,k) + step * v_tend_cell(i,j,k);
		}}
	}
}

/*********************************************************************
* Loop through velocity arrays, interpolate from cell centers to cell 
* faces, and then calculate rate of advection using the flux form of 
* the advection equation.
*
* @param step - fraction of total time step dt to advance
* @param il,ih,jl,jh - beginning and ending indices
**********************************************************************/
void advect_uvw_velocity(double step,int il,int ih,int jl,int jh){

	// Calculate the zonal flux on the eastern side of innermost column of boundary points.
	// This will become the flux on the western side of the leftmost interior (non-boundary) cell
	compute_west(il-1,jl,jh);

	//-------------------------------------------------------
	// For each point within the east-west range
	//-------------------------------------------------------
	for(int i=il;i<ih;i++){

		// compute all fluxes for a YZ cross section
		compute_fluxes(i,jl-1,jh);

		//-------------------------------------------------------
		// For each point within the north-south range
		//-------------------------------------------------------
		LOOP2D_JK(jl,jh,1,NZ-1, UP(i,j,k) = UM(i,j,k) + step * u_tend_cell(i,j,k) )
		LOOP2D_JK(jl,jh,1,NZ-1, VP(i,j,k) = VM(i,j,k) + step * v_tend_cell(i,j,k) )
		LOOP2D_JK(jl,jh,2,NZ-1, WP(i,j,k) = WM(i,j,k) + step * w_tend_cell(i,j,k) )	

	}
	//-------------------------------------------------------
	// Calculate buoyancy term base on the microphysics option
	//-------------------------------------------------------	
	if(MICROPHYSICS_OPTION == 0){
		LOOP3D_IJK(il,ih,jl,jh,2,NZ-1,	WP(i,j,k) = WM(i,j,k) + step * dt*grav*0.5*( buoyancy_dry(i,j,k)+buoyancy_dry(i,j,k-1))	)
	}
	if(MICROPHYSICS_OPTION != 0 && !USE_ICE){
		LOOP3D_IJK(il,ih,jl,jh,2,NZ-1,	WP(i,j,k) = WM(i,j,k) + step * dt*grav*0.5*( buoyancy_warm_micro(i,j,k)+buoyancy_warm_micro(i,j,k-1))	)
	}
	if(MICROPHYSICS_OPTION != 0 && USE_ICE && MICROPHYSICS_OPTION != 3){
		LOOP3D_IJK(il,ih,jl,jh,2,NZ-1,	WP(i,j,k) = WM(i,j,k) + step * dt*grav*0.5*( buoyancy_cold_micro(i,j,k)+buoyancy_cold_micro(i,j,k-1))	)
	}
	if(MICROPHYSICS_OPTION != 0 && USE_ICE && MICROPHYSICS_OPTION == 3){
		LOOP3D_IJK(il,ih,jl,jh,2,NZ-1,	WP(i,j,k) = WM(i,j,k) + step * dt*grav*0.5*( buoyancy_cold_micro_graupel(i,j,k)+buoyancy_cold_micro_graupel(i,j,k-1))	)
	}
	
}

/*********************************************************************
* Loop through potential temperature array, interpolate from cell
* centers to cell faces, and then calculate rate of advection using
* the flux form of the advection equation.
*
* @param step - fraction of total time step dt to advance
* @param il,ih,jl,jh - beginning and ending indices
**********************************************************************/
void advect_theta(double step,int il,int ih,int jl,int jh){

	// Calculate the zonal flux on the eastern side of innermost column of boundary points.
	// This will become the flux on the western side of the leftmost interior cell
	compute_west_scalar(il-1,jl,jh,&TH(0,0,0),&THBAR(0,0,0),&THCELL(0,0),&THBCELL(0,0));

	//-------------------------------------------------------
	// For each point within the east-west range
	//-------------------------------------------------------
	for(int i=il;i<ih;i++){

		// compute all fluxes for a YZ cross section
		compute_fluxes_scalar(i,jl-1,jh,&TH(0,0,0),&THBAR(0,0,0),&THCELL(0,0),&THBCELL(0,0));

		//-------------------------------------------------------
		// For each point within the north-south range
		//-------------------------------------------------------
		for(int j=jl;j<jh;j++){
		for(int k=1;k<NZ-1;k++){
		
			// new values equal old values plus the rate of change
			THP(i,j,k) = THM(i,j,k) + step * theta_tend_cell(i,j,k);
		}}
	}
}

/*********************************************************************
* Loop through scalar array, interpolate from cell
* centers to cell faces, and then calculate rate of advection using
* the flux form of the advection equation.
*
* @param step - fraction of total time step dt to advance
* @param il,ih,jl,jh - beginning and ending indices
**********************************************************************/
void advect_scalar(double step,int il,int ih,int jl,int jh,double *varm,double *var,double *varp,struct cell *scell){

	// Calculate the zonal flux on the eastern side of innermost column of boundary points.
	// This will become the flux on the western side of the leftmost interior cell
	compute_west_scalar(il-1,jl,jh,var,scell);

	//-------------------------------------------------------
	// For each point within the east-west range
	//-------------------------------------------------------
	for(int i=il;i<ih;i++){

		// compute all fluxes for a YZ cross section
		compute_fluxes_scalar(i,jl-1,jh,var,scell);

		//-------------------------------------------------------
		// For each point within the north-south range
		//-------------------------------------------------------
		for(int j=jl;j<jh;j++){
		for(int k=1;k<NZ-1;k++){
		
			// new values equal old values plus the rate of change
			varp[INDEX(i,j,k)] = varm[INDEX(i,j,k)] + step * scalar_tend_cell(i,j,k,scell);
		}}
	}
}

/*********************************************************************
* Loop through scalar array, interpolate from cell
* centers to cell faces, and then calculate rate of advection using
* the flux form of the advection equation.
*
* @param step - fraction of total time step dt to advance
* @param il,ih,jl,jh - beginning and ending indices
**********************************************************************/
void advect_scalar(double step,int il,int ih,int jl,int jh,double *varm,double *var,double *varp,double *varbase,struct cell *scell,struct cell *bcell){

	// Calculate the zonal flux on the eastern side of innermost column of boundary points.
	// This will become the flux on the western side of the leftmost interior cell
	compute_west_scalar(il-1,jl,jh,var,varbase,scell,bcell);

	//-------------------------------------------------------
	// For each point within the east-west range
	//-------------------------------------------------------
	for(int i=il;i<ih;i++){

		// compute all fluxes for a YZ cross section
		compute_fluxes_scalar(i,jl-1,jh,var,varbase,scell,bcell);

		//-------------------------------------------------------
		// For each point within the north-south range
		//-------------------------------------------------------
		for(int j=jl;j<jh;j++){
		for(int k=1;k<NZ-1;k++){
		
			// new values equal old values plus the rate of change
			varp[INDEX(i,j,k)] = varm[INDEX(i,j,k)] + step * scalar_tend_cell(i,j,k,scell);
		}}
	}
}

/*********************************************************************
* Loop through moisture variable arrays, interpolate from cell
* centers to cell faces, and then calculate rate of advection using
* the flux form of the advection equation.
*
* @param step - fraction of total time step dt to advance
* @param il,ih,jl,jh - beginning and ending indices
**********************************************************************/
void advect_microphysics_cell(double step,int il,int ih,int jl,int jh){

	//--------------------------------------------------------
	// If Eulerian rain fall integration, calculate fall speed
	//--------------------------------------------------------
	if(RAIN_FALLOUT==1 && MICROPHYSICS_OPTION!=3)
		calculate_eulerian_fall_speed_rain(vts, qrs, il, ih, jl, jh);

	// Calculate the zonal flux on the eastern side of innermost column of boundary points.
	// This will become the flux on the western side of the leftmost interior cell
	compute_west_moisture(il-1,jl,jh);

	//-------------------------------------------------------
	// For each point within the east-west range
	//-------------------------------------------------------
	for(int i=il;i<ih;i++){

		// compute all fluxes for a YZ cross section
		compute_fluxes_moisture(i,jl-1,jh);
#if 1
		//compute_fluxes_scalar(i,jl-1,jh,qvs,m_qbar,qvcell,qvbcell);
		LOOP2D_JK(jl,jh,1,NZ-1, QVP(i,j,k) = QVM(i,j,k) + step * qv_tend_cell(i,j,k) )
			
		//compute_fluxes_scalar(i,jl-1,jh,qcs,qccell);
		LOOP2D_JK(jl,jh,1,NZ-1, QCP(i,j,k) = QCM(i,j,k) + step * qc_tend_cell(i,j,k) )
			
		//compute_fluxes_scalar_with_fallspeed(i,jl-1,jh,qrs,vts,qrcell);
		LOOP2D_JK(jl,jh,1,NZ-1, QRP(i,j,k) = QRM(i,j,k) + step * qr_tend_cell(i,j,k) )
#else
		//-------------------------------------------------------
		// For each point within the north-south range
		//-------------------------------------------------------
		for(int j=jl;j<jh;j++){
		for(int k=1;k<NZ-1;k++){
		
			// new values equal old values plus the rate of change
			QVP(i,j,k) = QVM(i,j,k) + step * qv_tend_cell(i,j,k);	// vapor
			QCP(i,j,k) = QCM(i,j,k) + step * qc_tend_cell(i,j,k);	// cloud
			QRP(i,j,k) = QRM(i,j,k) + step * qr_tend_cell(i,j,k);	// rain

		}}
#endif
	}
}

/*********************************************************************
* Loop through ice arrays, interpolate from cell
* centers to cell faces, and then calculate rate of advection using
* the flux form of the advection equation.
*
* @param step - fraction of total time step dt to advance
* @param il,ih,jl,jh - beginning and ending indices
**********************************************************************/
void advect_ice_cell(double step,int il,int ih,int jl,int jh){
	
	if(MICROPHYSICS_OPTION!=3){
		//--------------------------------------------------------
		// If Eulerian snow/ice fall integration, calculate fall speed
		//--------------------------------------------------------
		if(RAIN_FALLOUT==1 && MICROPHYSICS_OPTION!=3)
			calculate_eulerian_fall_speed_snow_ice(sts, qss, its, pis, il, ih, jl, jh);

		//-------------------------------------------------------------
		// ADVECTION OF SNOW
		//-------------------------------------------------------------
		// Calculate the zonal flux on the eastern side of innermost column of boundary points.
		// This will become the flux on the western side of the leftmost interior cell
		compute_west_scalar(il-1,jl,jh,qss,ice_cell);

		//-------------------------------------------------------
		// For each point within the east-west range
		//-------------------------------------------------------
		for(int i=il;i<ih;i++){

			// compute all fluxes for a YZ cross section
			compute_fluxes_scalar_with_fallspeed(i,jl-1,jh,qss,sts,ice_cell);

			//-------------------------------------------------------
			// For each point within the north-south range
			//-------------------------------------------------------
			for(int j=jl;j<jh;j++){
			for(int k=1;k<NZ-1;k++){
			
				// new values equal old values plus the rate of change
				QSP(i,j,k) = QSM(i,j,k) + step * scalar_tend_cell(i,j,k,ice_cell);	// snow
			}}
		}
		
		//-------------------------------------------------------------
		// ADVECTION OF CLOUD ICE
		//-------------------------------------------------------------
		// Calculate the zonal flux on the eastern side of innermost column of boundary points.
		// This will become the flux on the western side of the leftmost interior cell
		compute_west_scalar(il-1,jl,jh,qis,ice_cell);

		//-------------------------------------------------------
		// For each point within the east-west range
		//-------------------------------------------------------
		for(int i=il;i<ih;i++){

			// compute all fluxes for a YZ cross section
			compute_fluxes_scalar_with_fallspeed(i,jl-1,jh,qis,its,ice_cell);

			//-------------------------------------------------------
			// For each point within the north-south range
			//-------------------------------------------------------
			for(int j=jl;j<jh;j++){
			for(int k=1;k<NZ-1;k++){
			
				// new values equal old values plus the rate of change
				QIP(i,j,k) = QIM(i,j,k) + step * scalar_tend_cell(i,j,k,ice_cell);	// snow
			}}
		}

	} else {
		advect_scalar(step,il,ih,jl,jh,qsms,qss,qsps,ice_cell);
		advect_scalar(step,il,ih,jl,jh,qims,qis,qips,ice_cell);
		advect_scalar(step,il,ih,jl,jh,qgms,qgs,qgps,ice_cell);
		advect_scalar(step,il,ih,jl,jh,nrms,nrs,nrps,ice_cell);
		advect_scalar(step,il,ih,jl,jh,nims,nis,nips,ice_cell);
	}

}

/*********************************************************************
* Loop through water vapor array, interpolate from cell
* centers to cell faces, and then calculate rate of advection using
* the flux form of the advection equation.
*
* @param step - fraction of total time step dt to advance
* @param il,ih,jl,jh - beginning and ending indices
**********************************************************************/
void advect_qv(double step,int il,int ih,int jl,int jh){

	// Calculate the zonal flux on the eastern side of innermost column of boundary points.
	// This will become the flux on the western side of the leftmost interior cell
	compute_west_scalar(il-1,jl,jh,&QV(0,0,0),&QBAR(0,0,0),&QVCELL(0,0),&QVBCELL(0,0));

	//-------------------------------------------------------
	// For each point within the east-west range
	//-------------------------------------------------------
	for(int i=il;i<ih;i++){

		// compute all fluxes for a YZ cross section
		compute_fluxes_scalar(i,jl-1,jh,&QV(0,0,0),&QBAR(0,0,0),&QVCELL(0,0),&QVBCELL(0,0));

		//-------------------------------------------------------
		// For each point within the north-south range
		//-------------------------------------------------------
		for(int j=jl;j<jh;j++){
		for(int k=1;k<NZ-1;k++){
		
			// new values equal old values plus the rate of change
			QVP(i,j,k) = QVM(i,j,k) + step * qv_tend_cell(i,j,k);
		}}
	}
}
