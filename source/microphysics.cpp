#include "stdafx.h"

/*******************************************************************************
* 
********************************************************************************/
#define SAT_VAP_WAT(t) ( 611.2 * exp(17.67 * (t-273.15) / (t - 29.65) ) )
#define SAT_VAP_ICE(t) ( 611.2 * exp(21.8745584 * (t-273.15) / (t - 7.66) ) )
#define SAT_MIX_RATIO(e,p) ( 0.62197 * e / (p-e) )

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
			
			
			#if 0
			double asum = 0;
			double dsum = rainfall;
			
			for(int k=1;k<NZ;k++){
				
				asum += QA_mass[k]*(ZA_edges[k+1]-ZA_edges[k]);
				dsum += QD_mass[k]*(ZW(k+1)-ZW(k));
			}
			if(fabs(asum-dsum)>0.0001){
				
				printf("Rain mass compare: %d %d %f %f %f %f\n",big_i[i],big_j[j],asum,dsum,rainfall,dsum+rainfall);
			}
			#endif
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
		//printf("zero slope\n");
	} else{
		mid_slope = 0.5*(edge_slope0 + edge_slope1);
	}

	double qp_edge = 0.5 * mid_slope * (zin[p+1]-zin[p]) + qin[p];
	double qn_edge = 2.0 * qin[p] - qp_edge;
	
	if(qp_edge < 0 || qn_edge < 0){
		
		qp_edge = qin[p];
		qn_edge = qin[p];
		//printf("zero edge\n");
	}

	double t0 = (zl-zin[p]) / (zin[p+1]-zin[p]);
	double t1 = (zh-zin[p]) / (zin[p+1]-zin[p]);
	
	//printf("integrate over %f to %f\n",t0,t1);

	//if(debug){
		//printf("midpoint = %f %f %f\n",1000*((qp_edge - qn_edge) * 0.5 + qn_edge),1000*qp_edge,1000*qn_edge);
	//}

	return (qp_edge - qn_edge) * 0.5*t1*t1 + qn_edge * t1 - ((qp_edge - qn_edge) * 0.5*t0*t0 + qn_edge * t0);
}
#if 1
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

	int test_i = -199;
	int test_j = -59;


	double rainfall = 0;

	for(int i=0;i<zlevs;i++){ qout[i] = 0;}

	//if(big_i[ic]==test_i && big_j[jc]==test_j){
		//int i = 0;
		//int j = 0;
		
		//printf("e %d %d %d %d %f %f %f %f %f\n",big_i[ic],big_j[jc],i,j,zin[j],zin[j+1],zout[i],zout[i+1],qin[j]*1000);
		//}

	int highest_negative_level = 1;

	for(int j=1;j<zlevs-2;j++){
		
		int i = 0;
		
		if(zin[j] < 0 && zin[j+1] < 0){
			
			rainfall += qin[j] * (zin[j+1]-zin[j]);
			
			if(big_i[ic]==test_i && big_j[jc]==test_j){
					printf("e %d %d %d %d %f %f %f %f %f\n",big_i[ic],big_j[jc],i,j,zin[j],zin[j+1],zout[i],zout[i+1],qin[j]*1000);
				}
			
			highest_negative_level = j+1;
			
		}
	}
	
	if(zin[highest_negative_level]<0){
	
		rainfall += integrate_linear(zin,qin,highest_negative_level,zin[highest_negative_level],0,false) * 
			-zin[highest_negative_level] *  ((zin[highest_negative_level+1] - zin[highest_negative_level]) /
				 -zin[highest_negative_level])
				;
	}
	
	if(big_i[ic]==test_i && big_j[jc]==test_j){
			printf("rainfall = %f %f\n",rainfall*1000,-zin[highest_negative_level]);
		}
	
		if(big_i[ic]==test_i && big_j[jc]==test_j){
				printf("lowest level = %d\n",highest_negative_level);
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
			if(big_i[ic]==test_i && big_j[jc]==test_j){
						printf("e %d %d %d %d %f %f %f %f %f\n",big_i[ic],big_j[jc],i,j,zin[j],zin[j+1],zout[i],zout[i+1],qin[j]*1000);
				}
			qout[i] += integrate_linear(zin,qin,j,zin[j],zout[i+1],false) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}
		//-----------------------------------------------------------
		// Overlap (mA4 -> m3)
		//-----------------------------------------------------------
		if(zout[i] < zin[j+1] && zout[i+1] > zin[j+1] && zout[i] > zin[j]){
			if(big_i[ic]==test_i && big_j[jc]==test_j){
						printf("c %d %d %d %d %f %f %f %f %f\n",big_i[ic],big_j[jc],i,j,zin[j],zin[j+1],zout[i],zout[i+1],qin[j]*1000);
							}
			qout[i] += integrate_linear(zin,qin,j,zout[i],zin[j+1],false) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}
		//-----------------------------------------------------------
		// Output cell falls completely within input cell (mA5 -> m4)
		//-----------------------------------------------------------
		if(zout[i] >= zin[j] && zout[i+1] <= zin[j+1]){
			if(big_i[ic]==test_i && big_j[jc]==test_j){
					printf("a %d %d %d %d %f %f %f %f %f\n",big_i[ic],big_j[jc],i,j,zin[j],zin[j+1],zout[i],zout[i+1],qin[j]*1000);
				}
			qout[i] += integrate_linear(zin,qin,j,zout[i],zout[i+1],false) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}
		//-----------------------------------------------------------
		// Input cell falls completely within output cell
		//-----------------------------------------------------------		
		if(zout[i] < zin[j] && zout[i+1] > zin[j+1]){
			if(big_i[ic]==test_i && big_j[jc]==test_j){
					printf("d %d %d %d %d %f %f %f %f %f\n",big_i[ic],big_j[jc],i,j,zin[j],zin[j+1],zout[i],zout[i+1],qin[j]*1000);
				}
			qout[i] += integrate_linear(zin,qin,j,zin[j],zin[j+1],false) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}	
	}}
	if(big_i[ic]==test_i && big_j[jc]==test_j){
			printf("output rainfall %f %d %d\n",rainfall,big_i[ic],big_j[jc]);
		}

	if(rainfall < 0){
		//printf("output rainfall %f %d %d %f %f\n",rainfall,big_i[ic],big_j[jc],outLons[big_i[ic]],outLats[big_j[jc]]);
	}

	return rainfall;

}
#endif

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

#if 0
/*********************************************************************
* Piecewise constant interpolation which preserves total mass
*
* zin - input height levels
* zout - output height levels
* qin - input field
* qout - output field
* zlevs number of levels (i.e. array length for all input arrays)
**********************************************************************/
double piecewise_interp(double *zin,double *zout,double *qin,double *qout,int zlevs,int ic,int jc){

	double factor;

	int test_i = 13;
	int test_j = 24;

	double test[43];
	double qin_sum[43];
	double qin_sum_a[43];
	double qin_sum_b[43];
	double qin_sum_c[43];
	double zt0[43];
	double zb0[43];
	double zt1[43];
	double zb1[43];

	double rainfall = 0;

	for(int i=0;i<zlevs;i++){ qout[i] = 0; test[i] = 0; qin_sum[i] = 0; qin_sum_a[i] = 0; qin_sum_b[i] = 0; qin_sum_c[i] = 0;
		zt0[i] = 0;zb0[i] = 0;zt1[i] = 0;zb1[i] = 0;}

	for(int i=1;i<zlevs-2;i++){
	for(int j=1;j<zlevs-2;j++){
		

		if(zout[i] >= zin[j] && zout[i+1] <= zin[j+1]){
			
			factor = 1.0;//(zout[i+1] - zin[j]) / (zout[i+1]-zout[i]);
			
			test[i] += factor;
			
			qout[i] += integrate_linear(zin,qin,j,zout[i],zout[i+1],false) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
			qin_sum[j] += integrate_linear(zin,qin,j,zout[i],zout[i+1],false);
			qin_sum_a[j] = integrate_linear(zin,qin,j,zout[i],zout[i+1],false);
			if(ic==test_i && jc==test_j){
				//printf("a %f %f %f %f\n",zout[i],1000*qin[j],factor,1000*integrate_linear(zin,qin,j,zout[i],zout[i+1],false));
			}
			//qout[i] += factor*qin[j];
		}


		if(zout[i] < zin[j] && zout[i+1] < zin[j+1] && zout[i+1] > zin[j]){
			
			factor = (zout[i+1] - zin[j]) / (zout[i+1]-zout[i]);
			
			test[i] += factor;
			qin_sum[j] += integrate_linear(zin,qin,j,zin[j],zout[i+1],false);
			qin_sum_b[j] = integrate_linear(zin,qin,j,zin[j],zout[i+1],false);
			zb0[j] = zin[j];
			zt0[j] = zout[i+1];
			//qout[i] += factor*qin[j];
			if(ic==test_i && jc==test_j){
				//printf("b %f %f %f %f\n",zout[i],1000*qin[j],factor,1000*integrate_linear(zin,qin,j,zin[j],zout[i+1],false));
				
				
			}
			qout[i] += integrate_linear(zin,qin,j,zin[j],zout[i+1],false) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}
		
		if(zout[i] < zin[j+1] && zout[i+1] > zin[j+1] && zout[i] > zin[j]){
			
			factor = (zin[j+1] - zout[i]) / (zout[i+1]-zout[i]);
			
			test[i] += factor;
			qin_sum[j] += integrate_linear(zin,qin,j,zout[i],zin[j+1],false);
			qin_sum_c[j] = integrate_linear(zin,qin,j,zout[i],zin[j+1],false);
			zb1[j] = zout[i];
			zt1[j] = zin[j+1];
			//qout[i] += factor*qin[j];
			if(ic==test_i && jc==test_j){
				//printf("c %f %f %f %f\n",zout[i],1000*qin[j],factor,1000*integrate_linear(zin,qin,j,zout[i],zin[j+1],false));
			}
			qout[i] += integrate_linear(zin,qin,j,zout[i],zin[j+1],false) * ((zin[j+1] - zin[j]) / (zout[i+1]-zout[i]));
		}
	}}

	//if(ic==test_i && jc==test_j){
		
		//for(int i=0;i<zlevs;i++){
			
			//if(zin[i] > 0 && fabs(qin[i]-qin_sum[i]) > 0.0000001){
			//if(qin[10]*1000>1.0)
				//printf("%d %d %d %f %f| %f %f %f\n",i,ic,jc,qin[i]*1000,qin_sum[i]*1000,qin_sum_a[i]*1000,qin_sum_b[i]*1000,qin_sum_c[i]*1000);
				//}
		
			//}
		//}

		
	if(false && ic==test_i && jc==test_j){
		for(int i=0;i<zlevs;i++){ printf("%d %f %f %f\n",i,test[i],ZW(i),qout[i]*1000);}
		for(int i=0;i<zlevs;i++){ printf("%d %f %f %f %f %f %f %f %f %f\n",i,
		qin_sum[i]*1000,qin_sum_b[i]*1000,qin_sum_c[i]*1000,qin[i]*1000,
		zb0[i],zt0[i],zb1[i],zt1[i],1000*integrate_linear(zin,qin,i,zin[i],zin[i+1],false));}
		
		int i = 17;
		printf("integration test %f %f\n",
		1000*integrate_linear(zin,qin,i,zin[i],0.5*(zin[i+1]+zin[i]),false)+
		1000*integrate_linear(zin,qin,i,0.5*(zin[i+1]+zin[i]),zin[i+1],false),
		1000*integrate_linear(zin,qin,i,zin[i],zin[i+1],false)
		);
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
double piecewise_interp2(double *zin,double *zout,double *qin,double *qout,int zlevs,int ic,int jc){

	double factor;

	int test_i = 13;
	int test_j = 24;

	double test[43];

	double rainfall = 0;

	for(int i=0;i<zlevs;i++){ qout[i] = 0; test[i] = 0;}

	for(int i=1;i<zlevs-1;i++){
	for(int j=1;j<zlevs-1;j++){
		
		//if(zin[i] < zout[0]){
			
			//factor = (zout[i+1] - zout[i]) / (zin[j+1]-zin[j]);
			
			//rainfall += factor*qin[j];
			//}
		
		if(zout[i] > zin[j]){
			
			if(zout[i+1] <= zin[j+1] ){ // mA3->m1
				
				factor = (zout[i+1] - zout[i]) / (zout[i+1]-zout[i]);
				
				//qout[i] += integrate_linear(zin,qin,j,zout[i],zout[i+1]);
				qout[i] += factor*qin[j];
				test[i] += factor;
				if(ic==test_i && jc==test_j){ printf("a %d %f %f %f %f %f %f\n",i,zin[j],zin[j+1],zout[i],zout[i+1],qin[j]*1000,factor);}
			}
			
			if(zout[i+1] == zin[j+1] ){
			
				factor = (zout[i] - zin[j]) / (zout[i+1]-zout[i]);
				test[i] += factor;
				//qout[i-1] += factor * qin[j];
				if(ic==test_i && jc==test_j){ printf("b %d %f %f %f %f %f %f\n",i,zin[j],zin[j+1],zout[i],zout[i+1],qin[j]*1000,factor);}
			}
		}
		
		if(zout[i] > zin[j] && zout[i] < zin[j+1]){
			
			 if(zout[i-1] < zin[j] ){ // mA4-> m2 (lower portion of mA4)
			
				factor = ( zout[i] - zin[j]) / (zout[i+1]-zout[i]);
				test[i-1] += factor;
				//qout[i-1] += integrate_linear(zin,qin,j,zin[j],zout[i]);
				qout[i-1] += factor * qin[j];
				
				if(ic==test_i && jc==test_j){ printf("c %d %f %f %f %f %f %f\n",i-1,zin[j],zin[j+1],zout[i-1],zout[i],qin[j]*1000,factor);}
			}
		
			 if(zout[i+1] > zin[j+1] ){ // mA4 -> m3 (upper portion of mA4)
				
				factor = (zin[j+1] - zout[i]) / (zout[i+1]-zout[i]);
				test[i] += factor;
				//qout[i] += integrate_linear(zin,qin,j,zout[i],zin[j+1]);//factor * qin[j];
				qout[i] += factor * qin[j];
				
				if(ic==test_i && jc==test_j){ printf("d %d %f %f %f %f %f %f\n",i,zin[j],zin[j+1],zout[i],zout[i+1],qin[j]*1000,factor);}
			} 
		}
#if 0
		if(zout[i] < zin[j] && zout[i+1] > zin[j+1]){
			
			factor = 1.0;
			test[j] += factor;
			qout[i] += factor * qin[j];
			
			if(ic==test_i && jc==test_j){ printf("g %d %f %f %f %f %f %f\n",j,zin[j],zin[j+1],zout[i],zout[i+1],qin[j]*1000,factor);}
		}
#endif			
		if(zout[i] == zin[j]){
		
			 if(zout[i+1] > zin[j+1] ){
			
				factor = (zin[j+1] - zout[i]) / (zout[i+1]-zout[i]);
				test[i] += factor;
				qout[i] += factor * qin[j];
				
				if(ic==test_i && jc==test_j){ printf("e %d %f %f %f %f %f %f\n",i,zin[j],zin[j+1],zout[i],zout[i+1],qin[j]*1000,factor);}
				
			} else if(zout[i+1] == zin[j+1] ){
				
				qout[i] += qin[j];
				test[i] += 1.0;
				if(ic==test_i && jc==test_j){ printf("f %d %f %f %f %f %f %f\n",i,zin[j],zin[j+1],zout[i],zout[i+1],qin[j]*1000,1.0);}
			}
		}

	}}
	if(ic==test_i && jc==test_j){
		for(int i=0;i<zlevs;i++){ printf("%d %f\n",i,test[i]);}
	}

	return rainfall;

}
#endif

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