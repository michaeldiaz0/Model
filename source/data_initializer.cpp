#include "stdafx.h"
#include "interpolate.h"
#include "surface.h"
#include "laplacian.h"
#include "boundaries.h"
#include "initializer.h"
#include "data_initializer.h"

/******************************************************************************
* Handles to details of putting data onto the model grids. Focuses on
* initialization from code and initialization from reanalysis/processed datasets
*
*******************************************************************************/

#define LOOP3D_IJK(il,ih,jl,jh,kl,kh,equation) 	for(int i=il;i<ih;i++){ \
												for(int j=jl;j<jh;j++){ \
												for(int k=kl;k<kh;k++){ \
													equation;	\
												}}}

#define d(i,j,k) (xdim*ydim*(k) + xdim*(j) + i)
#define ALLOC(v,size) v = (double *)calloc(size,sizeof(double))
#define IVORT(i,j,k)     ivort[FULL_ARRAY_INDEX(i,j,k)]
#define ISTREAM(i,j,k) istream[FULL_ARRAY_INDEX(i,j,k)]

int rindex(int r,int k,int NR);
void rinterp(int xc,int yc,int x,int y,int *ind, double *frac, double offset);
void initialize_vortex_perturbation();
void initialize_vortex(double, double,double,double,double,double,double,double,double,double,int,int,double *u,double *v,double *th,double *q);
void initialize_basic_state_idealized();
void initialize_basic_state_idealized2();
void initialize_basic_state_ITCZ_shear();
void initialize_basic_state_ITCZ_shear2(double,double,double);
void get_base_state_from_20_88(double *,double *);

/******************************************************************************
* This subroutine is called from the main initializer to initialize the basic
* state.
*
*******************************************************************************/
void initialize_basic_state_from_subroutine(){

	if(inputs.has_shear){
	
		initialize_basic_state_ITCZ_shear2(inputs.hor_shear,inputs.vert_shear0,inputs.vert_shear1);
	
		//LOOP3D_IJK(0,NX,0,NY,0,NZ,	IQBAR(i,j,k) = -(cp * ITHBAR(i,j,k)) / Lv	)
	
	} else {
	
		initialize_basic_state_ITCZ_shear2(1.0e-4,-2.0e-3,-1.0e-3);
	}
	
	//initialize_basic_state_idealized2();
	//initialize_basic_state_idealized();
}

/******************************************************************************
* This subroutine is called from the main initializer to initialize the
* perturbation.
*******************************************************************************/
void initialize_perturbation_from_subroutine(){
		
	if(inputs.vortex_initialize == 1){
	
		initialize_vortex_perturbation();
	}
}

/******************************************************************************
* Create a localized strip of vorticity using cosine functions
*
* my_hor_shear - peak vorticity value at center of vortex strip
* lat0 - southern end of vortex strip in degrees latitude
* lat1 - northern end of vortex strip in degrees latitude
* lon0 - western end of vortex strip in degrees longitude (use -1 for left edge of domain)
* lon1 - eastern end of vortex strip in degrees longitude (use -1 for right edge of domain)
* top_height - top boundary of vortex strip in meters
* buffer_vort_degrees - length of transition to zero on left and right ends of vortex strip (degrees)
* vort - pointer to output array of vorticity
*******************************************************************************/
void create_horizontal_shear_zone(double my_hor_shear,
									double lat0,double lat1,
									double lon0,double lon1,
									double top_height,double buffer_vort_degrees,
									double *vort){
	
	int low_lat,high_lat,low_lon,high_lon;
	
	if(lon0==-1){
		low_lon = 0;
	} else {
		low_lon = get_point_from_lon(lon0);
	}
	
	if(lon1==-1){
		high_lon = NX;
	} else {
		high_lon = get_point_from_lon(lon1);
	}
	
	low_lat = get_point_from_lat(lat0);
	high_lat = get_point_from_lat(lat1);
	
	int top_vort = get_point_from_height(top_height);
	int buffer_vort = buffer_vort_degrees * 111000 / dx;
	
	double mult;

	//----------------------------------------------------------------
	// Set up the vorticity distribution
	//----------------------------------------------------------------
	for(int i=low_lon-buffer_vort;i<high_lon+buffer_vort;i++){
	for(int j=low_lat;j<high_lat;j++){
		
		for(int k=0;k<top_vort;k++){
			
			mult = 0.5*cos_profile(-0.1,0.5,frac_distance(ZU(k),0,top_height),1.0) * 

				   0.5*cos_profile(-0.5,0.5,frac_distance(j,low_lat,high_lat),1);

	   	 	if(i<low_lon){
				mult *= 0.5*cos_profile(-0.5,0.5,frac_distance(j,low_lon-buffer_vort,low_lon),1);
	   	 	}
			
			if(i > high_lon){
				mult *= 0.5*cos_profile(-0.5,0.5,frac_distance(j,high_lon,high_lon+buffer_vort),1);
			}
	   	 
	   
			vort[FULL_ARRAY_INDEX(i,j,k)] = my_hor_shear * mult;
					   
		}
	}}
	
}

/******************************************************************************
* Given a distribution of vorticity, calculate stream function
*
* vort - pointer to vorticity arrau
* stream - pointer to output stream function array
*******************************************************************************/
void get_stream_function(double *vort,double *stream){
	
	double *vortXY = (double*) calloc(NX*NY,sizeof(double));

	for(int k=0;k<NZ;k++){

		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
	
			vortXY[i*NY+j]= vort[FULL_ARRAY_INDEX(i,j,k)];
		}}

		run_laplacian_solver(vortXY);

		for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
	
			stream[FULL_ARRAY_INDEX(i,j,k)] = vortXY[i*NY+j];
			
		}}
	}
	
	free(vortXY);
}

/******************************************************************************
* Get Moist Static Energy
*******************************************************************************/
double get_MSE(double theta,double pressure,double qv,double z){

	return cp * theta*pressure + Lv*qv +  grav*z;
}

/******************************************************************************
* Flatten MSE gradient
*******************************************************************************/
void flatten_MSE_gradient(){

	LOOP3D_IJK(0,NX,0,NY,0,NZ,	QBAR(i,j,k) = -(cp * THBAR(i,j,k)) / Lv	)
		
}

/******************************************************************************
* Create vertical shear profile
*
* my_vert_shear0 - vertical shear at bottom
* my_vert_shear1 - vertical shear at top
* midlev - level of zero wind (meters)
* ktrop - vertical level above which wind shear is constant (grid point)
*******************************************************************************/
void create_vertical_shear_profile(double my_vert_shear0,double my_vert_shear1,double midlev,int ktrop,double *wind_profile){

	int btop = NZ;
	double dUdzTop = my_vert_shear0;
	double dUdzBot = my_vert_shear1;
	
	int kmidlev = get_point_from_height(midlev);

	double *dUdz = (double*) calloc(NZ,sizeof(double));

	for(int k=0;k<NZ;k++){ dUdz[k] = 0;}
	
	for(int k=0;k<=ktrop;k++){
		
		dUdz[k] = dUdzBot + (dUdzTop - dUdzBot) / (ZU(ktrop)-ZU(0)) * ZU(k);
		//printf("%d %f %f\n",k,ZU(k)*0.001,dUdz[k]*1000);
	}
	
	for(int k=ktrop;k<NZ;k++){ dUdz[k] = dUdzTop;}

	wind_profile[kmidlev] = 0;

	for(int k=kmidlev+1;k<btop;k++){ wind_profile[k] = wind_profile[k-1] + dUdz[k] * DZW(k);}
	
	for(int k=kmidlev-1;k>=0;k--){ wind_profile[k] = wind_profile[k+1] - dUdz[k] * DZW(k+1);}
	
	free(dUdz);
}

/******************************************************************************
* Get dimensional pressure from winds assuming geostrophic balance
*
*
*******************************************************************************/
void get_pressure_from_winds(double *pres,double *u){
	
	int mid_lat = 0;
	double umid;
	
	//----------------------------------------------------------------
	// Calculate the pressure profile assuming geostrophic balance
	//----------------------------------------------------------------	
	for(int i=0;i<NX-1;i++){	
	for(int k=0;k<NZ;k++){
	for(int j=1;j<NY;j++){
		
		umid = 0.25*( u[FULL_INDEX(i,j-1,k)] +
					  u[FULL_INDEX(i+1,j-1,k)] +
					  u[FULL_INDEX(i,j,k)] +
					  u[FULL_INDEX(i+1,j,k)]
					);

		
		pres[FULL_INDEX(i,j,k)] = pres[FULL_INDEX(i,j-1,k)] -
			( umid * f[j] + 
			  0*umid*umid*sin(outLats[j]*trigpi/180.0) / (R_earth*cos(outLats[j]*trigpi/180.0)) 
			) * dy;
		
		//if(i==NX/2 && k==16){ printf("%d %f %f %f %f\n",j,outLats[j],pres[FULL_ARRAY_INDEX(i,j-1,k)],umid,f[j]);}
	}}}
	
}

/******************************************************************************
* Subtract the base state from a basic state variable
*
* var - pointer to 3d variables
* base_i - x position within the array of base state
* base_j - y position within the array of base state
*******************************************************************************/
void subtract_base_state(double *var,int base_i,int base_j,int il,int ih,int jl,int jh){
	
	double *store = (double*) calloc(NZ,sizeof(double));
	
	for(int k=0;k<NZ;k++){ store[k] = var[FULL_INDEX(base_i,base_j,k)];}
	
	LOOP3D_IJK(il,ih,jl,jh,0,NZ,	var[FULL_INDEX(i,j,k)] -= store[k]	);

	free(store);

}

/******************************************************************************
* Create a basic state temperature field in hydrostatic balance with the
* pressure field 
*
* temp - pointer to output potential temperature field
* pres - pointer to input pressure field (non-dimensional units)
*******************************************************************************/
void hydrostatic_balance(double *pres,double *temp,int il,int ih,int jl,int jh){
	
	LOOP3D_IJK(il,ih,jl,jh,1,NZ,
		temp[FULL_INDEX(i,j,k)] = tbv[k]*(pres[FULL_INDEX(i,j,k)] - pres[FULL_INDEX(i,j,k-1)]) / (grav * (ZU(k)-ZU(k-1)) )
	)
}


/******************************************************************************
* Initialize a monsoon-like basic state with both horizontal and vertical shear
*
* hor_shear - value of vorticity at center of shear zone
* vert_shear0 - value of vertical shear in the uppermost troposphere
* vert_shear1 - value of vertical shear in the lowermost troposphere
*******************************************************************************/
void initialize_basic_state_ITCZ_shear2(double hor_shear,double vert_shear0,double vert_shear1){
	
	setup_memory_allocation();

	setup_vertical_height_levels();

	//initialize_from_era(&ZU(0));	// will use the basic state stratification

	setup_coordinates(0.000000,-89.462997);

	init_topography_to_zero();
	
	init_friction();

	initialize_coriolis_parameter(&outLats[0]);
	
	init_laplacian_solver_periodic_EW_zerograd_NS(NX,NY);
	
	double control_vort = 1.0e-4;	// value of vorticity used to construct humidity field
	
	//----------------------------------------------------------------
	// Get the base state values from reanalysis data
	//----------------------------------------------------------------
	int base_i = get_point_from_lon(88);
	int base_j = get_point_from_lat(20);
	
	//initialize_vertical_basic_state(base_i,base_j);
	
	get_base_state_from_20_88(tb,qb);
	
	initialize_vertical_basic_state2(tb,qb);
		
	if(VERBOSE){ print_vertical_basic_state();}
	
	//----------------------------------------------------------------
	// Initialize basic state arrays to zero
	//----------------------------------------------------------------
	int full_size = NX*NY*NZ;
	
	memset(iubar, 0,full_size*sizeof(double));
	memset(ivbar, 0,full_size*sizeof(double));
	memset(iwbar, 0,full_size*sizeof(double));
	memset(ithbar,0,full_size*sizeof(double));
	memset(iqbar, 0,full_size*sizeof(double));
	memset(ipbar, 0,full_size*sizeof(double));
	
	//----------------------------------------------------------------
	// Create vorticity strip
	//----------------------------------------------------------------
	double * ivort   = (double *)calloc(NX*NY*NZ,sizeof(double));
	double * istream = (double *)calloc(NX*NY*NZ,sizeof(double));

	create_horizontal_shear_zone(control_vort,18,21,-1,-1,17000,0,ivort);

	get_stream_function(ivort,istream);
	
	//----------------------------------------------------------------
	// Set up the vertical shear
	//----------------------------------------------------------------
	double *u_profile = (double*) calloc(NZ,sizeof(double));
	
	create_vertical_shear_profile(vert_shear0,vert_shear1,4000,37,&u_profile[0]);
	
	LOOP3D_IJK(0,NX,0,NY,0,NZ,	IUBAR(i,j,k) = u_profile[k] + 3	)
	
	//----------------------------------------------------------------
	// Add in the horizontal shear
	//----------------------------------------------------------------
	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=0;k<NZ-1;k++){
	
		IUBAR(i,j,k) += -0.25 * ( (ISTREAM(i-1,j+1,k) + ISTREAM(i  ,j+1,k)) - (ISTREAM(i-1,j-1,k) + ISTREAM(i  ,j-1,k)) ) * one_d_dy;
		IVBAR(i,j,k) += 0.25 * (  (ISTREAM(i+1,j-1,k) + ISTREAM(i+1,j+1,k)) - (ISTREAM(i-1,j-1,k) + ISTREAM(i-1,j+1,k)) ) * one_d_dx;

	}}}

	mirror_boundaries(&IUBAR(0,0,0));

	//----------------------------------------------------------------
	// Create the pressure field
	//----------------------------------------------------------------
	get_pressure_from_winds(ipbar,iubar);

	subtract_base_state(ipbar,base_i,base_j,0,NX,0,NY);
	
	mirror_boundaries(&IPBAR(0,0,0));

	// upper and lower boundaries
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		IPBAR(i,j,0) = IPBAR(i,j,1) - (ZU(1)-ZU(0)) * (IPBAR(i,j,2)-IPBAR(i,j,1)) / ( ZU(2)-ZU(1) );
		IPBAR(i,j,NZ-1) = IPBAR(i,j,NZ-2);
	}}

	//----------------------------------------------------------------
	// Construct hydrostatically balanced temperature field
	//----------------------------------------------------------------	
	hydrostatic_balance(ipbar,ithbar,0,NX,0,NY);

	//----------------------------------------------------------------
	// Convert to non-dimensional pressure
	//----------------------------------------------------------------
	LOOP3D_IJK(0,NX,0,NY,0,NZ,	IPBAR(i,j,k) = IPBAR(i,j,k) / (cp*tbv[k]) + pib[k]	)
	
#if 1
	//-----------------------------------------------------------------------
	// Relative humidity changes
	//-----------------------------------------------------------------------
	double smr,temp,pres,rh;
	double *base_relh = (double*) calloc(NZ,sizeof(double));
	
	int rh_y_min = get_point_from_lat(5);
	int rh_y_max = get_point_from_lat(15);
	int rh_y_min2 = get_point_from_lat(26);
	int rh_y_max2 = get_point_from_lat(30);
	int rh_y_min3 = get_point_from_lat(20);

	double rh_decrease = 0.1;//0.2;//0.25;
	double rh_decrease2 = 0.1;//0.1;
	double rh_decrease3 = 0;//0;
	
	double rh_slope = rh_decrease / (double)(rh_y_max-rh_y_min);
	double rh_slope2 = rh_decrease2 / (double)(rh_y_max2-rh_y_min2);
	double rh_slope3 = rh_decrease3 / (double)(rh_y_min2-rh_y_min3);

	double rh_drop[NY];
	
	for(int j=0;j<NY;j++){ rh_drop[j] = 0;} // initialize

	//---------------------------------------------------------------
	// First
	//---------------------------------------------------------------	
	for(int j=rh_y_max;j>rh_y_min;j--){		rh_drop[j] = rh_drop[j+1] + rh_slope;	}
	
	for(int j=rh_y_min;j>=0;j--){			rh_drop[j] = rh_decrease;				}
	//---------------------------------------------------------------
	// Third
	//---------------------------------------------------------------
	for(int j=rh_y_min3;j<rh_y_min2;j++){	rh_drop[j] = rh_drop[j-1] + rh_slope3;	}
	//---------------------------------------------------------------
	// Second
	//---------------------------------------------------------------
	for(int j=rh_y_min2;j<rh_y_max2;j++){	rh_drop[j] = rh_drop[j-1] + rh_slope2;	}
	
	for(int j=rh_y_max2;j<NY;j++){			rh_drop[j] = rh_drop[j-1];				}

	//for(int j=0;j<NY;j++){ printf("%f %f\n",outLats[j],rh_drop[j]);}
	//-----------------------------------------------------------------------
	// Calculate basic state relative humidity
	//-----------------------------------------------------------------------
	for(int k=0;k<NZ;k++){
		
		temp = tb[k]*pib[k];
		pres = p0*pow(pib[k],(cp/Rd));
		
		base_relh[k] = qb[k] / get_QV_Sat(temp,pres);
		//printf("%f %f %f %f\n",ZU(k),tb[k],qb[k],base_relh[k]);
		//printf("RH %d %f %f %f %f\n",k,zu[k],base_pres[k],base_relh[k],0.15*cos_profile(-0.25,0.25,(zu[k]-2000.0)/5000.0,0));
	}
	//-----------------------------------------------------------------------
	// Set mixing ratio field based on the relative humidity field
	//-----------------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=1;k<NZ-1;k++){

		temp = (ITHBAR(i,j,k)+tb[k]) * IPBAR(i,j,k);	// full, actual temperature
		pres = p0*pow(IPBAR(i,j,k),(cp/Rd));			// full, dimensional pressure

		smr = get_QV_Sat(temp,pres);					// calculate saturation mixing ratio

		//if(outLats[j]>15){ IQBAR(i,j,k) = 0;}
		//else{
		IQBAR(i,j,k) = (base_relh[k] - rh_drop[j]) * smr - qb[k];
			//}

		
		//if(IQBAR(i,j,k)>0.002){ IQBAR(i,j,k) = 0.002;}
		//if(i==NX/2 && k==9){ printf("%d %f %f\n",j,outLats[j],IQBAR(i,j,k)*1000.0);}
		//IQBAR(i,j,k) = 0;
	}}}
#endif


#if 1
	//-------------------------------------------------------------------------------------
	// Recalculate wind, pressure, and temperature fields for final value of vorticity strip
	//-------------------------------------------------------------------------------------
	create_horizontal_shear_zone(hor_shear - control_vort, 18,21,-1,-1,17000,0,ivort);

	get_stream_function(ivort,istream);

	//----------------------------------------------------------------
	// Add in the horizontal shear
	//----------------------------------------------------------------
	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=0;k<NZ-1;k++){
	
		IUBAR(i,j,k) += -0.25 * ( (ISTREAM(i-1,j+1,k) + ISTREAM(i  ,j+1,k)) - (ISTREAM(i-1,j-1,k) + ISTREAM(i  ,j-1,k)) ) * one_d_dy;
		IVBAR(i,j,k) += 0.25 * (  (ISTREAM(i+1,j-1,k) + ISTREAM(i+1,j+1,k)) - (ISTREAM(i-1,j-1,k) + ISTREAM(i-1,j+1,k)) ) * one_d_dx;

	}}}
	
	mirror_boundaries(&IUBAR(0,0,0));

	//----------------------------------------------------------------
	// Create the pressure field
	//----------------------------------------------------------------
	get_pressure_from_winds(ipbar,iubar);

	subtract_base_state(ipbar,base_i,base_j,0,NX,0,NY);
	
	mirror_boundaries(&IPBAR(0,0,0));

	// upper and lower boundaries
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		IPBAR(i,j,0) = IPBAR(i,j,1) - (ZU(1)-ZU(0)) * (IPBAR(i,j,2)-IPBAR(i,j,1)) / ( ZU(2)-ZU(1) );
		IPBAR(i,j,NZ-1) = IPBAR(i,j,NZ-2);
	}}

	//----------------------------------------------------------------
	// Construct hydrostatically balanced temperature field
	//----------------------------------------------------------------	
	hydrostatic_balance(ipbar,ithbar,0,NX,0,NY);

	//----------------------------------------------------------------
	// Convert to non-dimensional pressure
	//----------------------------------------------------------------
	LOOP3D_IJK(0,NX,0,NY,0,NZ,	IPBAR(i,j,k) = IPBAR(i,j,k) / (cp*tbv[k]) + pib[k]	)
	
#endif

	//----------------------------------------------------------------
	// Upper and lower boundary conditions for base state array
	//----------------------------------------------------------------
	upper_lower_boundaries(&ITHBAR(0,0,0));
	upper_lower_boundaries(&IQBAR(0,0,0));
	upper_lower_boundaries(&IUBAR(0,0,0));	
	upper_lower_boundaries(&IVBAR(0,0,0));	
		
	//----------------------------------------------------------------
	// Initialize friction array
	//----------------------------------------------------------------
	init_friction();

	//----------------------------------------------------------------
	// Calculate base state vertical velocity
	//----------------------------------------------------------------
	init_basic_state_vertical_velocity();
	
	free(ivort); free(istream); free(u_profile); free(base_relh);
}



/******************************************************************************
* Initiliazes a baroclinic jet basic state
*
*******************************************************************************/
void initialize_basic_state_idealized(){

	double *zgrid;
	
	if(STRETCHED_GRID){ zgrid = &zsu[0];}
	else { 				zgrid = &zu[0]; }

	double f0 = 2*7.292e-5*sin(22*3.1416/180.0);

	setup_memory_allocation();

	setup_coordinates(0,-90);

	setup_vertical_height_levels();

	initialize_coriolis_parameter(&outLats[0]);

	init_friction();


	//----------------------------------------------------------------
	// Initialize topographic array
	//----------------------------------------------------------------
	init_topography_to_zero();

	double jet_height[NY];
	double rh_profile[NY];

	int lowlat = get_point_from_lat(35);
	int higlat = get_point_from_lat(55);
	int midlat = (higlat+lowlat)/2;
	
	
	int jet_height0 = get_point_from_lat(35);
	int jet_height1 = get_point_from_lat(55);
	double jet_drop = 3000;	// meters
	
	int tropopause_height = 31;//25;
	
	double jet_max = 50.0;
	double N02T = 1.3e-4;
	double N02S = 4.0e-4;//4 * N02T;
	
	//printf("tropopause height is %f \n",zgrid[tropopause_height]);
	
	for(int j=0;j<jet_height0;j++){ jet_height[j] = zgrid[tropopause_height];}
	for(int j=jet_height1;j<NY;j++){ jet_height[j] = zgrid[tropopause_height]-jet_drop;}
	for(int j=jet_height0;j<jet_height1;j++){ jet_height[j] = zgrid[tropopause_height] - jet_drop*(double)(j-jet_height0)/(double)(jet_height1-jet_height0);}
	
	//for(int j=0;j<NY;j++){ printf("%f\n",jet_height[j]);}
	
	for(int i=0;i<NX;i++){
	for(int j=lowlat;j<higlat;j++){
		
		for(int k=0;k<NZ;k++){
	
			if(zgrid[k]<jet_height[j]){
	
				IUBAR(i,j,k) = (zgrid[k]/jet_height[j]) * jet_max*0.5*cos_profile(-0.5,0.5,(double)(j-lowlat)/(double)(higlat-lowlat),1);
			} else {
				IUBAR(i,j,k) = (jet_height[j]/zgrid[k]) * jet_max*0.5*cos_profile(-0.5,0.5,(double)(j-lowlat)/(double)(higlat-lowlat),1);
		
			}
		
			//if(k==10 && i==NX/2){ printf("%d %f %f\n",j,outLats[j],IUBAR(i,j,k)); }
		}
		
		//for(int k=jet_height[j];k<NZ;k++){
	
			
			//if(k==tropopause_height && i==NX/2){ printf("%d %f %f\n",j,outLats[j],IUBAR(i,j,k)); }
			//}
		
	}}

	tb[0] = 280;//273;
	tb[1] = tb[0];

	for(int k=2;k<tropopause_height;k++){
		
		tb[k] = N02T * tb[k-1] / grav * (zgrid[k]-zgrid[k-1]) + tb[k-1];
		//printf("%f %f\n",zgrid[k],tb[k]);
	}
	for(int k=tropopause_height;k<NZ-1;k++){
		
		tb[k] = N02S * tb[k-1] / grav * (zgrid[k]-zgrid[k-1]) + tb[k-1];
		//printf("%f %f\n",zgrid[k],tb[k]);
	}

	tb[NZ-1] = tb[NZ-2];

	double max_humid = 0.75;
	double min_humid = 0.20;

	for(int k=0;k<NZ;k++){

		if(zgrid[k]<1500){ rh_profile[k] = max_humid;}
		if(zgrid[k]<4000 && zgrid[k]>=1500){ rh_profile[k] = max_humid - (max_humid-min_humid)*(zgrid[k]-1500)/(4000-1500);}		
		if(zgrid[k]>=4000){ rh_profile[k] = min_humid;}
		//printf("%f %f\n",zgrid[k],rh_profile[k]);
	}

	for(int k=0;k<NZ;k++){ qb[k] = 0;}

	initialize_vertical_basic_state2(tb,qb);

	double base_temp,base_pres;

	for(int k=0;k<NZ;k++){
		
		base_temp = tb[k]*pib[k];
		base_pres = p0*pow(pib[k],(cp/Rd));
		
		qb[k] = get_QV_Sat(base_temp,base_pres) * rh_profile[k];
	}

	initialize_vertical_basic_state2(tb,qb);

	//for(int k=0;k<NZ;k++){ printf("%d %f %f %f %f\n",k,zgrid[k],pib[k],tb[k]*pib[k],qb[k]*1000);}


	for(int i=0;i<NX-1;i++){	
	for(int k=0;k<NZ;k++){	
	for(int j=1;j<NY;j++){
		
		IPBAR(i,j,k) = IPBAR(i,j-1,k) - 0.25*(IUBAR(i,j-1,k)+IUBAR(i+1,j-1,k)+IUBAR(i,j,k)+IUBAR(i+1,j,k)) * dy * FC(j);
		
	}}}


	for(int i=0;i<NX;i++){			
	for(int j=0;j<NY;j++){
	for(int k=1;k<NZ;k++){

		ITHBAR(i,j,k) = tbv[k]*(IPBAR(i,j,k) - IPBAR(i,j,k-1)) / (grav * (zgrid[k]-zgrid[k-1]) );// + tb[k];
	}}}
	
	
	for(int i=0;i<NX;i++){			
	for(int j=0;j<NY;j++){
		ITHBAR(i,j,0) = ITHBAR(i,j,1);
		ITHBAR(i,j,NZ-1) = ITHBAR(i,j,NZ-2);
	}}


	//----------------------------------------------------------------
	// Convert to non-dimensional pressure
	//----------------------------------------------------------------
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){

		IPBAR(i,j,k) = IPBAR(i,j,k) / (cp*tbv[k]) + pib[k];

	}}}


	double temp,pres,smr;

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=1;k<NZ-1;k++){

		temp = (ITHBAR(i,j,k)+tb[k]) * IPBAR(i,j,k);	// full, actual temperature
		pres = p0*pow(IPBAR(i,j,k),(cp/Rd));			// full, dimensional pressure

		smr = get_QV_Sat(temp,pres);					// calculate saturation mixing ratio

		IQBAR(i,j,k) = rh_profile[k] * smr - qb[k];

	}}}

	for(int i=0;i<NX;i++){			
	for(int j=0;j<NY;j++){
		IQBAR(i,j,0) = IQBAR(i,j,1);
		IQBAR(i,j,NZ-1) = IQBAR(i,j,NZ-2);
	}}

	mirror_boundaries(&IUBAR(0,0,0));
	mirror_boundaries(&IPBAR(0,0,0));
	mirror_boundaries(&ITHBAR(0,0,0));
	mirror_boundaries(&IVBAR(0,0,0));
	mirror_boundaries(&IQBAR(0,0,0));

	print_vertical_basic_state();

	//exit(0);

}

/*********************************************************************
* Initialize the perturbation fields
*
**********************************************************************/
void initialize_vortex_perturbation(){


	//----------------------------------------------------------------
	// Actions to perform on only the root process
	//----------------------------------------------------------------
	if(!PARALLEL || (PARALLEL && rank==0)){
	
		//----------------------------------------------------------------
		// Set the initial perturbation values to zero
		//----------------------------------------------------------------
		if(PARALLEL){
			memset( &IUBAR(0,0,0),0,NX*NY*NZ*sizeof(double));
			memset( &IVBAR(0,0,0),0,NX*NY*NZ*sizeof(double));
			memset(&ITHBAR(0,0,0),0,NX*NY*NZ*sizeof(double));
			memset( &IQBAR(0,0,0),0,NX*NY*NZ*sizeof(double));
			memset( &IPBAR(0,0,0),0,NX*NY*NZ*sizeof(double));
		}
		
		//----------------------------------------------------------------
		// Do the initialization
		//----------------------------------------------------------------
		int lat = get_point_from_lat(inputs.vortex_latitude);
		int lon = get_point_from_lon(inputs.vortex_longitude);


		printf("Initializing vortex at lat = %f lon = %f\n",outLats[lat],outLons[lon]);

		//initialize_vortex(1100000,3000, 11000, -2.05, 1.5,NX/2,lat);		
		//initialize_vortex(1500000,3000, 11000, -4.5, 4.05,NX/3+4*111000.0/dx,lat);
		//initialize_vortex(1500000,3000, 11000, -4.5, 4.05,NX/3,lat);
		//initialize_vortex(1500000,6000, 14000, -4.5, 4.05,NX/3,lat);
		//initialize_vortex(1000000,6000, 16000, 0.5, -0.5,lon,lat);
		//initialize_vortex(1500000,6000, 16000, -2.0, 2.0,NX/2-500.0/15.0,get_point_from_lat(18));

		if(PARALLEL){

			initialize_vortex(			
				inputs.vortex_radius,
				inputs.vortex_height_wind_max,
				inputs.vortex_vertical_extent,
				inputs.vortex_upper_warm_anomaly,
				inputs.vortex_lower_cold_anomaly,				
				inputs.vortex_rh_bottom,
				inputs.vortex_rh_top,
				inputs.vortex_rh_prime,
				inputs.vortex_rh_radius,
				inputs.vortex_rh_max,
				lon,lat,iubar,ivbar,ithbar,iqbar
			);	
			//initialize_vortex(500000*2,3000, 11000, 2.5, -2.05,lon,lat,iubar,ivbar,ithbar,iqbar);
		} else {
			
			initialize_vortex(
				inputs.vortex_radius,
				inputs.vortex_height_wind_max,
				inputs.vortex_vertical_extent,
				inputs.vortex_upper_warm_anomaly,
				inputs.vortex_lower_cold_anomaly,			
				inputs.vortex_rh_bottom,
				inputs.vortex_rh_top,
				inputs.vortex_rh_prime,
				inputs.vortex_rh_radius,
				inputs.vortex_rh_max,
				lon,lat,us,vs,ths,qvs
			);
			//initialize_vortex(500000*3,3000, 11000, 2.5, -2.05,lon,lat,ums,vms,thms,qvms);
		}
	}

	//----------------------------------------------------------------
	// Distribute the perturbation arrays from the root process to
	// the other processes
	//----------------------------------------------------------------
	if(PARALLEL){

		distributeArray(ths,&ITHBAR(0,0,0));	distributeArray(thms,&ITHBAR(0,0,0));
		distributeArray(us,  &IUBAR(0,0,0));	distributeArray(ums,  &IUBAR(0,0,0));
		distributeArray(vs,  &IVBAR(0,0,0));	distributeArray(vms,  &IVBAR(0,0,0));

		exchange(ths); exchange(thms);
		exchange(us); exchange(ums);
		exchange(vs); exchange(vms);
					
		if(USE_MICROPHYSICS){
			
			distributeArray(qvs,&IQBAR(0,0,0));	distributeArray(qvms,&IQBAR(0,0,0));
			exchange(qvs); exchange(qvms);
		}
	
	} else {
		memcpy(ums,us,NX*NY*NZ*sizeof(double));
		memcpy(vms,vs,NX*NY*NZ*sizeof(double));
		memcpy(thms,ths,NX*NY*NZ*sizeof(double));
		//memcpy(qvs,qvms,NX*NY*NZ*sizeof(double));
	}
}


/*********************************************************************
*
*
**********************************************************************/
void initialize_vortex(double r_end, double zm, double zt, double tmax, double tmin,double zb_rh,double zt_rh,
					   double rh_prime,double r_end_rh,double rh_max,int xpos, int ypos,double *u,double *v,double *th,double *q){
	
	int NR = (int)(r_end / dx);
	
	double *t_z = (double*) calloc(NZ,sizeof(double));
	double *t_r = (double*) calloc(NR,sizeof(double));
	double *phi = (double*) calloc(NR*NZ,sizeof(double));
	double *vr  = (double*) calloc(NR*NZ,sizeof(double));
	double *rh  = (double*) calloc(NR*NZ,sizeof(double));

	double zmid = 0.5*(zm+zt);
	double R,tmid;
	
	//-------------------------------------------------------------------------------------
	// Calculate temperature field
	//-------------------------------------------------------------------------------------	
	for(int k=0;k<NZ;k++){
		
		if(zsu[k]<=zm){ t_z[k] = tmin * exp( -(zsu[k]*zsu[k])/(zm*zm/8.0) );		
		} else { t_z[k] = tmax * exp( -( (zsu[k]-zmid)*(zsu[k]-zmid) ) / ( (zt-zm)*(zt-zm) / 32.0) );}
	}

	for(int r=0;r<NR;r++){ t_r[r] = exp( -dx*r*dx*r / (r_end*r_end/8.0) );}

	//-------------------------------------------------------------------------------------
	// Calculate moisture field
	//-------------------------------------------------------------------------------------
	//double zt_rh = 11000;
	//double zb_rh = -2000;
	//double zt_rh = 14000;
	//double zb_rh = 1000;
	//double rh_prime = 0.20;
	//double r_end_rh = 1500000;//900000;
	//double rh_max = 0.95;
	
	zmid = 0.5*(zt_rh+zb_rh);
	
	//-------------------------------------------------------------------------------------
	// Anomalous relative humidity
	//-------------------------------------------------------------------------------------
	for(int k=0;k<NZ;k++){
	for(int r=0;r<NR;r++){
				
		rh[rindex(r,k,NR)] =
			rh_prime * exp( -( (zsu[k]-zmid)*(zsu[k]-zmid) ) / ( (zt_rh-zb_rh)*(zt_rh-zb_rh) / 32.0) )
			* exp( -dx*r*dx*r / (r_end_rh*r_end_rh/8.0) );
	}}
	//-------------------------------------------------------------------------------------
	// Convert to mixing ratio anomaly
	//-------------------------------------------------------------------------------------
	double temp,pres,smr,rhp;

	for(int k=0;k<NZ;k++){
		
		temp = tb[k]*pib[k];
		pres = p0*pow(pib[k],(cp/Rd));
		smr = get_QV_Sat(temp,pres);
		
		for(int r=0;r<NR;r++){ 
			
			if( rh[rindex(r,k,NR)] + qb[k]/smr < rh_max ){ rh[rindex(r,k,NR)] *= smr;} 
			else { rh[rindex(r,k,NR)] = (rh_max - qb[k]/smr) * smr;}
		
			//if(r==0){ printf("%f rh = %f %f %f %f\n",zsu[k],rh[rindex(r,k,NR)]*1000.0,qb[k]*1000.0,smr*1000.0,100*(rh[rindex(r,k,NR)]+qb[k])/smr);}
		}
	}
	//-------------------------------------------------------------------------------------
	// Calculate pressure field
	//-------------------------------------------------------------------------------------
	for(int k=NZ-2;k>0;k--){
	for(int r=0;r<NR;r++){
	
		tmid = 0.5* (t_z[k]*t_r[r]+t_z[k+1]*t_r[r]) / ((tb[k]+tb[k+1])*0.5) + 0.61*rh[rindex(r,k,NR)];
	
		phi[rindex(r,k,NR)] = phi[rindex(r,k+1,NR)] - grav * tmid * (zsu[k+1]-zsu[k]);
	

		//if(r==0){
			//printf("%d %f %f\n",k,t_z[k],phi[rindex(r,k,NR)]);
		//}

	}}
	
	//-------------------------------------------------------------------------------------
	// Caluculate gradient wind balance
	//-------------------------------------------------------------------------------------
	double fcor = 2*7.292e-5*sin(outLats[ypos]*trigpi/180.0);
	double radicand;
	
	for(int k=0;k<NZ;k++){
		
		vr[rindex(0,k,NR)] = 0;
		
		for(int r=1;r<NR-1;r++){
		
			R = r*dx;
	
			if(phi[rindex(r,k,NR)]<0){
	
				vr[rindex(r,k,NR)] = -0.5*fcor*R + 
					0.5*sqrt( fcor*fcor*R*R + 4*R * 0.5*
						fabs( phi[rindex(r+1,k,NR)]-phi[rindex(r-1,k,NR)] )
							 * one_d_dx );
			} else {
				
				radicand = fcor*fcor*R*R - 4*R * 0.5* fabs( phi[rindex(r+1,k,NR)]-phi[rindex(r-1,k,NR)] ) * one_d_dx ;
				
				if(radicand<0){ radicand = 0;}
				
				vr[rindex(r,k,NR)] = 0.5*fcor*R - 0.5*sqrt(radicand) ;
				
				vr[rindex(r,k,NR)] *= -1.0;
			}
		}
	}
	
	double iproj,jproj,windmag,frac;
	int ind;
	//-------------------------------------------------------------------
	// Loop over box containing the vortex
	//-------------------------------------------------------------------
	for(int i=xpos-NR;i<xpos+NR;i++){
	for(int j=ypos-NR;j<ypos+NR;j++){
	for(int k=1;k<NZ;k++){
		
		if(i>0 && j>0 && i<NX && j<NY){
			//---------------------------------------------------------------
			// Interpolate Cartesian to radial coordinates
			//---------------------------------------------------------------	
			rinterp(xpos,ypos,i,j,&ind,&frac,0);
			//---------------------------------------------------------------
			// Don't exceed array bounds
			//---------------------------------------------------------------		
			if(ind < NR-2){
				//-----------------------------------------------------------
				// Potential temperature
				//-----------------------------------------------------------
				th[FULL_ARRAY_INDEX(i,j,k)] += t_z[k] * ( (1.0-frac)*t_r[ind] + frac*t_r[ind+1] );
				//-----------------------------------------------------------
				// Water vapor
				//-----------------------------------------------------------
				//IQBAR(i-(int)(0*meters_per_degree/dx),j+(int)(0*meters_per_degree/dx),k) 
				//			= ( (1.0-frac)*rh[rindex(ind,k,NR)] + frac*rh[rindex(ind+1,k,NR)] );	// POSSIBLE ARRAY OUT OF BOUNDS!!!
				if(PARALLEL && USE_MICROPHYSICS && rh_prime > 1.0e-10){
					q[FULL_ARRAY_INDEX(i+(int)(0.0*meters_per_degree/dx),j-(int)(0.0*meters_per_degree/dx),k)] 
								= ( (1.0-frac)*rh[rindex(ind,k,NR)] + frac*rh[rindex(ind+1,k,NR)] );
				}
				//-----------------------------------------------------------
				// Wind components (non-staggered)
				//-----------------------------------------------------------					
				if(ind != 0){

					iproj = fabs( (double)(j-ypos)) / ((double)ind + frac);
					jproj = fabs( (double)(i-xpos)) / ((double)ind + frac);
				
					windmag = (1.0-frac)*vr[rindex(ind,k,NR)] + frac*vr[rindex(ind+1,k,NR)];

					if(j<ypos){ u[FULL_ARRAY_INDEX(i,j,k)] += +windmag * iproj; }
					else { 		u[FULL_ARRAY_INDEX(i,j,k)] += -windmag * iproj; }
		
					if(i<xpos){ v[FULL_ARRAY_INDEX(i,j,k)] += -windmag * jproj; }
					else { 		v[FULL_ARRAY_INDEX(i,j,k)] += +windmag * jproj; }
				}
			}
		}
	}}}
	
	double *temp_array = (double *)calloc(NX*NY*NZ,sizeof(double));
	
	//---------------------------------------------------------------
	// Stagger wind components
	//---------------------------------------------------------------	
	for(int i=1;i<NX;i++){
	for(int j=1;j<NY;j++){
	for(int k=0;k<NZ;k++){
	
		temp_array[FULL_ARRAY_INDEX(i,j,k)] = 0.5*(u[FULL_ARRAY_INDEX(i,j,k)] + u[FULL_ARRAY_INDEX(i-1,j,k)]);
	}}}
	
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
	
		u[FULL_ARRAY_INDEX(i,j,k)] = temp_array[FULL_ARRAY_INDEX(i,j,k)];
	}}}
	
	for(int i=1;i<NX;i++){
	for(int j=1;j<NY;j++){
	for(int k=0;k<NZ;k++){
	
		temp_array[FULL_ARRAY_INDEX(i,j,k)] = 0.5*(v[FULL_ARRAY_INDEX(i,j-1,k)] + v[FULL_ARRAY_INDEX(i,j,k)]);
	}}}
	
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
	
		v[FULL_ARRAY_INDEX(i,j,k)] = temp_array[FULL_ARRAY_INDEX(i,j,k)];
	}}}	
		
	free(t_z); free(t_r); free(vr); free(phi); free(temp_array);
}

/*********************************************************************
* 
**********************************************************************/
int rindex(int r,int k,int NR){ return r+k*NR;}

/*********************************************************************
* 
**********************************************************************/
void rinterp(int xc,int yc,int x,int y,int *ind, double *frac, double offset){
	
	double dist = sqrt( (double)(x-xc)*(x-xc) + (double)(y-yc)*(y-yc) ) + offset;
	
	//printf("d = %f",dist);
	
	*ind = (int)dist;
	
	*frac = dist - (double)((int)dist);
}

/********************************************************
*
* 
*
*********************************************************/
void initialize_from_ncar(){

	/********************************************************
	* Get dimensions of data in files to determine how 
	* much memory to allocate
	*********************************************************/
	size_t dims[3];

	get_dims(geoheight_file,"lon","lat","level",dims);

	size_t xdim = dims[0], ydim = dims[1], zdim =  dims[2];

	/********************************************************
	* Create arrays large enough to hold the variable
	*********************************************************/
	int size = xdim*ydim*zdim;

	double * hgt2 = (double *)malloc(size*sizeof(double));
	double * uz2 = (double *)malloc(size*sizeof(double));
	double * vz2 = (double *)malloc(size*sizeof(double));
	double * tz2 = (double *)malloc(size*sizeof(double));

	double * hgt = (double *)malloc(size*sizeof(double));
	double * uz = (double *)malloc(size*sizeof(double));
	double * vz = (double *)malloc(size*sizeof(double));
	double * tz = (double *)malloc(size*sizeof(double));
	
	double * lat2 = (double *)malloc(ydim*sizeof(double));
	double * lon = (double *)malloc(xdim*sizeof(double));
	double * lat = (double *)malloc(ydim*sizeof(double));
	double * levs = (double *)malloc(zdim*sizeof(double));

	double * topo2 = (double *)malloc(xdim*ydim*sizeof(double));
	double * topozinterp = (double *)malloc(xdim*ydim*NZ*sizeof(double));

	/********************************************************
	* Get data from file
	*********************************************************/
	get_data(geoheight_file,"hgt", size, hgt2);
	get_data(uwind_file,"uwnd", size, uz2);
	get_data(vwind_file,"vwnd", size, vz2);
	get_data(temp_file,"air", size, tz2);
	get_data(topo_file,"hgt", xdim*ydim, topo2);

	/********************************************************
	* Reverse the y-coordinate
	*********************************************************/
	for(size_t i=0;i<xdim;i++){
	for(size_t j=0;j<ydim;j++){

		for(size_t k=0;k<zdim;k++){

			uz[d(i,j,k)] = uz2[d(i,ydim-1-j,k)];
			vz[d(i,j,k)] = vz2[d(i,ydim-1-j,k)];
			hgt[d(i,j,k)] = hgt2[d(i,ydim-1-j,k)];
			tz[d(i,j,k)] = tz2[d(i,ydim-1-j,k)];
		}

		for(int k=0;k<NZ;k++)
			topozinterp[d(i,j,k)] = topo2[d2(i,ydim-1-j)];
	}}

	/********************************************************
	* Get and process coordinates from input data
	*********************************************************/
	get_data(geoheight_file,"lat", ydim, lat2);
	get_data(geoheight_file,"lon", xdim, lon);
	get_data(geoheight_file,"level", zdim, levs);

	for(size_t i=0;i<xdim;i++)
		if(lon[i] >= 290)
			lon[i] = lon[i] - 360;

	for(size_t i=0;i<ydim;i++){ lat[i] = lat2[ydim-1-i]; }	// reverse y-coordinate

	for(size_t k=0;k<zdim;k++){ levs[k] = levs[k]*100;}

	/********************************************************
	* Convert temperature to potential temperature
	*********************************************************/
	for(size_t k=0;k<zdim;k++)
		for(size_t j=0;j<ydim;j++)
			for(size_t i=0;i<xdim;i++)
				tz[d(i,j,k)] = tz[d(i,j,k)]*pow((100000./levs[k]),(287./1004.));

	/********************************************************
	* Vertical interpolation
	*********************************************************/
	double * z 		  = (double *)malloc(NZ*sizeof(double));

	z[0] = 0;
	for(int k=1;k<NZ;k++){ z[k] = ((double)k-0.5)*dz;}

	double * uzinterp = vert_interpolate(z,hgt,uz,xdim,ydim,zdim,NZ);
	double * vzinterp = vert_interpolate(z,hgt,vz,xdim,ydim,zdim,NZ);
	double * tzinterp = vert_interpolate(z,hgt,tz,xdim,ydim,zdim,NZ);

	/********************************************************
	* Horizontal interpolation
	*********************************************************/
	double lonoffset = 90;
	double latoffset = 1;
	double dlat = 2.5;
	double dlon = 2.5;

	horz_interpolate(uzinterp,&IUBAR(0,0,0),xdim,ydim,NZ,!false,!true,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(vzinterp,&IVBAR(0,0,0),xdim,ydim,NZ,!true,!false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(tzinterp,&ITHBAR(0,0,0),xdim,ydim,NZ,!true,!true,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(topozinterp,&ITOPO(0,0,0),xdim,ydim,NZ,!true,!true,dlat,dlon,lonoffset,latoffset);

	/********************************************************
	* Construct coordinate arrays for output data
	*********************************************************/
	outLons[0] = lon[0] + lonoffset + dx/meters_per_degree;
	outLats[0] = lat[0] + dy/meters_per_degree;	

	for(int i=1;i<NX;i++){ outLons[i] = outLons[i-1] + dx/meters_per_degree;}
	for(int i=1;i<NY;i++){ outLats[i] = outLats[i-1] + dy/meters_per_degree;}

	/********************************************************
	* Free up allocated memory
	*********************************************************/
	free(uzinterp); free(vzinterp); free(tzinterp); free(z);
	free(hgt2); free(uz2); free(vz2); free(tz2);
	free(hgt); free(uz); free(vz); free(tz);
	free(levs); free(lat); free(lon); free(lat2);
	free(topozinterp); free(topo2);
}

/********************************************************
*
* 
*
*********************************************************/
void initialize_from_era(double * z){

	/********************************************************
	* Get dimensions of data in files to determine how 
	* much memory to allocate
	*********************************************************/
	size_t dims[3];

	get_dims(datafile,"lon","lat","level",dims);

	size_t xdim = dims[0]; 
	size_t ydim = dims[1];
	size_t zdim = dims[2];

	int size = xdim*ydim*zdim;

	double *piIn = (double *)malloc(zdim*sizeof(double));

	//printf("Dimensions of input data are %zu %zu %zu\n",xdim,ydim,zdim);

	/********************************************************
	* Get data from file
	*********************************************************/
	double * hgt = get_data2(datafile,"hgt", size);
	double * uz  = get_data2(datafile,"uwnd", size);
	double * vz  = get_data2(datafile,"vwnd", size);
	double * tz  = get_data2(datafile,"temp", size);
	double * qz  = get_data2(datafile2,"qv", size);
	double * topo2 = get_data2(datafile,"zsfc", xdim*ydim);
	
	/********************************************************
	* Reverse the y- and z-coordinate
	*********************************************************/
	reverse_yz_coord(uz,xdim,ydim,zdim);
	reverse_yz_coord(vz,xdim,ydim,zdim);
	reverse_yz_coord(hgt,xdim,ydim,zdim);
	reverse_yz_coord(tz,xdim,ydim,zdim);
	reverse_yz_coord(qz,xdim,ydim,zdim);

	double * topozinterp = (double *)malloc(xdim*ydim*NZ*sizeof(double));

	for(size_t i=0;i<xdim;i++){
	for(size_t j=0;j<ydim;j++){
	for(size_t k=0;k < NZ;k++){
		topozinterp[d(i,j,k)] = topo2[d2(i,ydim-1-j)];
	}}}

	if(SHIFT_PRIME_MERIDIAN){

		flip_array(hgt,xdim,ydim,zdim);
		flip_array(uz,xdim,ydim,zdim);
		flip_array(vz,xdim,ydim,zdim);
		flip_array(tz,xdim,ydim,zdim);
		flip_array(qz,xdim,ydim,zdim);
		flip_array(topozinterp,xdim,ydim,NZ);
	}

	/********************************************************
	* Get and process coordinates from input data
	*********************************************************/
	double * lat2  = get_data2(datafile,"latitude", ydim);
	double * lon   = get_data2(datafile,"longitude", xdim);
	double * levs2 = get_data2(datafile,"level", zdim);

	double * lat  = (double *)malloc(ydim*sizeof(double));
	double * levs = (double *)malloc(zdim*sizeof(double));
	
	for(size_t i=0;i<ydim;i++){ lat[i] = lat2[ydim-1-i];}			// reverse y-coordinate

	for(size_t k=0;k<zdim;k++){ levs[k] = levs2[zdim-1-k]*100.0;}// printf("%d pres %f\n",k,levs[k]);}	// reverse z-coordinate

	for(size_t k=0;k<zdim;k++){ piIn[k] = pow((levs[k]/p0),Rd/cp);}

#if 0
	for(int i=0;i<xdim;i++){
	for(int j=0;j<ydim;j++){
		
		if(lat[j] > 19 && lat[j] < 21 && lon[i] > 88 && lon[i] < 90){
		
			for(int k=0;k<zdim;k++){
		
				printf("%d %f %f %f %f\n",k,levs[k],lat[j],lon[i],qz[d(i,j,k)]*1000);
			}	
		}
	}}
#endif
	/********************************************************
	* Convert temperature to potential temperature
	*********************************************************/
	for(size_t k=0;k<zdim;k++)
		for(size_t j=0;j<ydim;j++)
			for(size_t i=0;i<xdim;i++)
				tz[d(i,j,k)] = tz[d(i,j,k)]*pow((100000.0/levs[k]),(287.0/1004.0));

	/********************************************************
	* Vertical interpolation
	*********************************************************/
	//double * z = (double *)malloc(NZ*sizeof(double));

	// construct height coordinate
	//for(int k=0;k<NZ;k++){ z[k] = ((double)k-0.5)*dz;}

	double * uzinterp = vert_interpolate(z,hgt,uz,xdim,ydim,zdim,NZ);
	double * vzinterp = vert_interpolate(z,hgt,vz,xdim,ydim,zdim,NZ);
	double * tzinterp = vert_interpolate(z,hgt,tz,xdim,ydim,zdim,NZ);
	double * qzinterp = vert_interpolate(z,hgt,qz,xdim,ydim,zdim,NZ);
	double * pzinterp = vert_interpolate_1d(z,hgt,piIn,xdim,ydim,zdim,NZ);

	/********************************************************
	* Horizontal interpolation
	*********************************************************/
	double dlat = lat[3]-lat[2];
	double dlon = lon[3]-lon[2];
	//printf("dlat = %f dlon = %f\n",dlat,dlon);
	horz_interpolate(uzinterp,&IUBAR(0,0,0),xdim,ydim,NZ,true,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(vzinterp,&IVBAR(0,0,0),xdim,ydim,NZ,false,true,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(tzinterp,&ITHBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(qzinterp,&IQBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(pzinterp,&IPBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(topozinterp,&ITOPO(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);

	/********************************************************
	* Construct coordinate arrays for output data
	*********************************************************/
	setup_coordinates(lon[0],lat[0]);

	//printf("%f %f\n",lon[0],lat[0]);

	/********************************************************
	* Free up allocated memory
	*********************************************************/
	free(uzinterp); free(vzinterp); free(tzinterp); free(qzinterp); free(pzinterp);
	free(hgt); free(uz); free(vz); free(tz); free(qz);
	free(levs); free(lat); free(lon); free(lat2); free(levs2); 
	free(topozinterp); free(topo2);	free(piIn);
}

/********************************************************
*
* 
*
*********************************************************/
void initialize_from_ERA5(double * z){

	/********************************************************
	* Get dimensions of data in files to determine how 
	* much memory to allocate
	*********************************************************/
	size_t dims[3];

	
	//const char datafile1[100] = "/Users/michaeldiaz/Desktop/output.nc";
	const char datafile1[100] = "/Users/michaeldiaz/Desktop/output.nc";
	//const char datafile3[100] = "/Users/michaeldiaz/Desktop/era5jan2018g.nc";
	const char datafile2[100] = "/Users/michaeldiaz/Desktop/era5surface.nc";

	get_dims(datafile1,"longitude","latitude","level",dims);

	size_t xdim = dims[0]; 
	size_t ydim = dims[1];
	size_t zdim = dims[2];

	int size = xdim*ydim*zdim;

	double *piIn = (double *)malloc(zdim*sizeof(double));

	//printf("Dimensions of input data are %zu %zu %zu\n",xdim,ydim,zdim);

	/********************************************************
	* Get data from file
	*********************************************************/
	double * hgt = get_data2(datafile1,"z", size);
	double * uz  = get_data2(datafile1,"u", size);
	double * vz  = get_data2(datafile1,"v", size);
	double * tz  = get_data2(datafile1,"t", size);
	double * qz  = get_data2(datafile1,"q", size);
	double * topo2 = get_data2(datafile2,"z", xdim*ydim);
	
	for(int i=0;i<size;i++){ hgt[i] /= 9.80665;}
	for(int i=0;i<xdim*ydim;i++){ topo2[i] /= 9.80665;}
	
	/********************************************************
	* Reverse the y- and z-coordinate
	*********************************************************/
	reverse_yz_coord(uz,xdim,ydim,zdim);
	reverse_yz_coord(vz,xdim,ydim,zdim);
	reverse_yz_coord(hgt,xdim,ydim,zdim);
	reverse_yz_coord(tz,xdim,ydim,zdim);
	reverse_yz_coord(qz,xdim,ydim,zdim);

	double * topozinterp = (double *)malloc(xdim*ydim*NZ*sizeof(double));

	for(size_t i=0;i<xdim;i++){
	for(size_t j=0;j<ydim;j++){
	for(size_t k=0;k < NZ;k++){
		topozinterp[d(i,j,k)] = topo2[d2(i,ydim-1-j)];
	}}}

	if(SHIFT_PRIME_MERIDIAN){

		flip_array(hgt,xdim,ydim,zdim);
		flip_array(uz,xdim,ydim,zdim);
		flip_array(vz,xdim,ydim,zdim);
		flip_array(tz,xdim,ydim,zdim);
		flip_array(qz,xdim,ydim,zdim);
		flip_array(topozinterp,xdim,ydim,NZ);
	}

	/********************************************************
	* Get and process coordinates from input data
	*********************************************************/
	double * lat2  = get_data2(datafile1,"latitude", ydim);
	double * lon   = get_data2(datafile1,"longitude", xdim);
	double * levs2 = get_data2(datafile1,"level", zdim);

	double * lat  = (double *)malloc(ydim*sizeof(double));
	double * levs = (double *)malloc(zdim*sizeof(double));
	
	for(size_t i=0;i<ydim;i++){ lat[i] = lat2[ydim-1-i];}			// reverse y-coordinate

	for(size_t k=0;k<zdim;k++){ levs[k] = levs2[zdim-1-k]*100.0;}// printf("%d pres %f\n",k,levs[k]);}	// reverse z-coordinate

	for(size_t k=0;k<zdim;k++){ piIn[k] = pow((levs[k]/p0),Rd/cp);}

	/********************************************************
	* Convert temperature to potential temperature
	*********************************************************/
	for(size_t k=0;k<zdim;k++)
		for(size_t j=0;j<ydim;j++)
			for(size_t i=0;i<xdim;i++)
				tz[d(i,j,k)] = tz[d(i,j,k)]*pow((100000.0/levs[k]),(287.0/1004.0));

	/********************************************************
	* Vertical interpolation
	*********************************************************/
	//double * z = (double *)malloc(NZ*sizeof(double));

	// construct height coordinate
	//for(int k=0;k<NZ;k++){ z[k] = ((double)k-0.5)*dz;}

	double * uzinterp = vert_interpolate(z,hgt,uz,xdim,ydim,zdim,NZ);
	double * vzinterp = vert_interpolate(z,hgt,vz,xdim,ydim,zdim,NZ);
	double * tzinterp = vert_interpolate(z,hgt,tz,xdim,ydim,zdim,NZ);
	double * qzinterp = vert_interpolate(z,hgt,qz,xdim,ydim,zdim,NZ);
	double * pzinterp = vert_interpolate_1d(z,hgt,piIn,xdim,ydim,zdim,NZ);

	/********************************************************
	* Horizontal interpolation
	*********************************************************/
	double dlat = lat[3]-lat[2];
	double dlon = lon[3]-lon[2];
	//printf("dlat = %f dlon = %f\n",dlat,dlon);
	horz_interpolate(uzinterp,&IUBAR(0,0,0),xdim,ydim,NZ,true,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(vzinterp,&IVBAR(0,0,0),xdim,ydim,NZ,false,true,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(tzinterp,&ITHBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(qzinterp,&IQBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(pzinterp,&IPBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);
	horz_interpolate(topozinterp,&ITOPO(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset,latoffset);

	/********************************************************
	* Construct coordinate arrays for output data
	*********************************************************/
	double lon_shift;
	
	if(SHIFT_PRIME_MERIDIAN){ lon_shift = 180;} else { lon_shift = 0;}

	outLons[0] = lon[0] + lonoffset + dx/meters_per_degree - lon_shift;
	outLats[0] = lat[0] + latoffset + dy/meters_per_degree;

	for(int i=1;i<NX;i++){ outLons[i] = outLons[i-1] + dx/meters_per_degree;}
	for(int i=1;i<NY;i++){ outLats[i] = outLats[i-1] + dy/meters_per_degree;}

	/********************************************************
	* Free up allocated memory
	*********************************************************/
	free(uzinterp); free(vzinterp); free(tzinterp); free(qzinterp); free(pzinterp);
	free(hgt); free(uz); free(vz); free(tz); free(qz);
	free(levs); free(lat); free(lon); free(lat2); free(levs2); 
	free(topozinterp); free(topo2);	free(piIn);
}


/********************************************************
*
* 
*
*********************************************************/
void initialize_from_ERA5_processed(double * z){

	/********************************************************
	* Get dimensions of data in files to determine how 
	* much memory to allocate
	*********************************************************/
	size_t dims[3];

	get_dims(datafile,"lon","lat","lev",dims);

	size_t xdim = dims[0]; 
	size_t ydim = dims[1];
	size_t zdim = dims[2];

	int size = xdim*ydim*zdim;

	double *piIn = (double *)malloc(zdim*sizeof(double));

	printf("Dimensions of input data are %zu %zu %zu\n",xdim,ydim,zdim);
	
	/********************************************************
	* Get data from file
	*********************************************************/
	double * hgt = get_data2(datafile,"hgt", size);
	double * uz  = get_data2(datafile,"uwnd", size);
	double * vz  = get_data2(datafile,"vwnd", size);
	double * tz  = get_data2(datafile,"temp", size);
	double * qz  = get_data2(datafile,"qv", size);
	//double * topo2 = get_data2(datafile2,"zsfc", xdim*ydim);

	for(int i=0;i<size;i++){ hgt[i] /= 9.80665;}
	
	/********************************************************
	* Reverse the y- and z-coordinate
	*********************************************************/
	reverse_yz_coord(uz,xdim,ydim,zdim);
	reverse_yz_coord(vz,xdim,ydim,zdim);
	reverse_yz_coord(hgt,xdim,ydim,zdim);
	reverse_yz_coord(tz,xdim,ydim,zdim);
	reverse_yz_coord(qz,xdim,ydim,zdim);

	//double * topozinterp = (double *)malloc(xdim*ydim*NZ*sizeof(double));
	/*
	for(size_t i=0;i<xdim;i++){
	for(size_t j=0;j<ydim;j++){
	for(size_t k=0;k < NZ;k++){
		topozinterp[d(i,j,k)] = topo2[d2(i,ydim-1-j)];
	}}}
*/
	if(SHIFT_PRIME_MERIDIAN){

		flip_array(hgt,xdim,ydim,zdim);
		flip_array(uz,xdim,ydim,zdim);
		flip_array(vz,xdim,ydim,zdim);
		flip_array(tz,xdim,ydim,zdim);
		flip_array(qz,xdim,ydim,zdim);
		//flip_array(topozinterp,xdim,ydim,NZ);
	}

	/********************************************************
	* Get and process coordinates from input data
	*********************************************************/
	double * lat2  = get_data2(datafile,"latitude", ydim);
	double * lon   = get_data2(datafile,"longitude", xdim);
	double * levs2 = get_data2(datafile,"level", zdim);

	double * lat  = (double *)malloc(ydim*sizeof(double));
	double * levs = (double *)malloc(zdim*sizeof(double));
	
	for(size_t i=0;i<ydim;i++){ lat[i] = lat2[ydim-1-i];}			// reverse y-coordinate

	for(size_t k=0;k<zdim;k++){ levs[k] = levs2[zdim-1-k]*100.0;}// printf("%d pres %f\n",k,levs[k]);}	// reverse z-coordinate

	for(size_t k=0;k<zdim;k++){ piIn[k] = pow((levs[k]/p0),Rd/cp);}

#if 0
	for(int i=0;i<xdim;i++){
	for(int j=0;j<ydim;j++){
		
		if(lat[j] > 19 && lat[j] < 21 && lon[i] > 88 && lon[i] < 90){
		
			for(int k=0;k<zdim;k++){
		
				printf("%d %f %f %f %f %f\n",k,levs[k],lat[j],lon[i],qz[d(i,j,k)]*1000,hgt[d(i,j,k)]);
			}	
		}
	}}
#endif
	//exit(0);
	/********************************************************
	* Convert temperature to potential temperature
	*********************************************************/
	for(size_t k=0;k<zdim;k++)
		for(size_t j=0;j<ydim;j++)
			for(size_t i=0;i<xdim;i++)
				tz[d(i,j,k)] = tz[d(i,j,k)]*pow((100000.0/levs[k]),(287.0/1004.0));

	/********************************************************
	* Vertical interpolation
	*********************************************************/
	//double * z = (double *)malloc(NZ*sizeof(double));

	// construct height coordinate
	//for(int k=0;k<NZ;k++){ z[k] = ((double)k-0.5)*dz;}



	double * uzinterp = vert_interpolate(z,hgt,uz,xdim,ydim,zdim,NZ);
	double * vzinterp = vert_interpolate(z,hgt,vz,xdim,ydim,zdim,NZ);
	double * tzinterp = vert_interpolate(z,hgt,tz,xdim,ydim,zdim,NZ);
	double * qzinterp = vert_interpolate(z,hgt,qz,xdim,ydim,zdim,NZ);
	double * pzinterp = vert_interpolate_1d(z,hgt,piIn,xdim,ydim,zdim,NZ);

#if 0
	for(int i=0;i<xdim;i++){
	for(int j=0;j<ydim;j++){
		
		if(lat[j] > 19 && lat[j] < 21 && lon[i] > 88 && lon[i] < 90){
		
			for(int k=0;k<NZ;k++){
		
				printf("%d %f %f %f %f %f\n",k,z[k],lat[j],lon[i],qzinterp[d(i,j,k)]*1000,pzinterp[d(i,j,k)]);
			}	
		}
	}}
#endif

	/********************************************************
	* Horizontal interpolation
	*********************************************************/
	double dlat = lat[3]-lat[2];
	double dlon = lon[3]-lon[2];
	printf("dlat = %f dlon = %f\n",dlat,dlon);
	horz_interpolate(uzinterp,&IUBAR(0,0,0),xdim,ydim,NZ,true,false,dlat,dlon,lonoffset-lon[0],latoffset-lat[0]-90);
	horz_interpolate(vzinterp,&IVBAR(0,0,0),xdim,ydim,NZ,false,true,dlat,dlon,lonoffset-lon[0],latoffset-lat[0]-90);
	horz_interpolate(tzinterp,&ITHBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset-lon[0],latoffset-lat[0]-90);
	horz_interpolate(qzinterp,&IQBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset-lon[0],latoffset-lat[0]-90);
	horz_interpolate(pzinterp,&IPBAR(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset-lon[0],latoffset-lat[0]-90);
	//horz_interpolate(topozinterp,&ITOPO(0,0,0),xdim,ydim,NZ,false,false,dlat,dlon,lonoffset-lon[0],latoffset-lat[0]-90);


	/********************************************************
	* Construct coordinate arrays for output data
	*********************************************************/
	double lon_shift;
	
	if(SHIFT_PRIME_MERIDIAN){ lon_shift = 180;} else { lon_shift = 0;}

	outLons[0] = lonoffset + dx/meters_per_degree - lon_shift;
	outLats[0] = latoffset - 90 + dy/meters_per_degree;

	for(int i=1;i<NX;i++){ outLons[i] = outLons[i-1] + dx/meters_per_degree;}
	for(int i=1;i<NY;i++){ outLats[i] = outLats[i-1] + dy/meters_per_degree;}


#if 0
	for(int i=0;i<NX;i++){
		printf("%f\n",outLons[i]);
	for(int j=0;j<NY;j++){
		
		if(outLats[j] > 19 && outLats[j] < 21 && outLons[i] > 88 && outLons[i] < 90){
		
			for(int k=0;k<NZ;k++){
		
				printf("%d %f %f %f %f %f\n",k,z[k],outLats[j],outLons[i],IQBAR(i,j,k)*1000,ITHBAR(i,j,k));
			}	
		}
	}}
#endif

	interpolate_terrain(datafile2);

	/********************************************************
	* Free up allocated memory
	*********************************************************/
	free(uzinterp); free(vzinterp); free(tzinterp); free(qzinterp); free(pzinterp);
	free(hgt); free(uz); free(vz); free(tz); free(qz);
	free(levs); free(lat); free(lon); free(lat2); free(levs2); 
	//free(topozinterp); free(topo2);
	free(piIn);
}


/********************************************************
* Values from the DiazBoos2019 basic state at 20N, 88E
* 
*
*********************************************************/
void get_base_state_from_20_88(double * tb_out,double *qb_out){

	const int size = 37;

	double tb_input[size] =
		{301.3688257116,301.8163590723,302.3725230915,303.5396454461,304.4779570865,305.4938954838,
		306.6227667893,307.8911249353,309.2454008827,310.6936441680,312.2381766112,315.4918989755,
		318.7967845618,321.8109652832,325.3755268822,330.0211783507,334.5667331740,339.1455453915,
		343.9030848237,348.2139891673,352.0680543027,353.6331293068,355.0080137497,356.5539533635,
		358.7368337136,362.1445309838,373.4121328070,429.4978633148,492.4724739660,593.3140366746,
		681.8870810434,857.2924158974,967.5530849371,1092.1733790195,1319.2834107809,1518.7488915123,
		1882.5817225223};
	
	double qb_input[size] =
		{0.0199446990,0.0190561418,0.0181378907,0.0171884018,0.0164599732,0.0157072981,
		0.0149655603,0.0142159289,0.0134605601,0.0126884429,0.0119014775,0.0103284225,
		0.0088135030,0.0074222058,0.0062207775,0.0051422818,0.0039585675,0.0027820404,
		0.0017543257,0.0009854162,0.0004371278,0.0002502617,0.0001286210,0.0000598493,
		0.0000244848,0.0000085451,0.0000032968,0.0000025111,0.0000024027,0.0000024293,
		0.0000025173,0.0000026986,0.0000028049,0.0000029210,0.0000031283,0.0000033081,
		0.0000036400};

	double z_input[size] =
		{-4.9309824095,220.0563275839,449.7083182255,684.3637424527,924.2080450767,1169.5747824156,
		1420.7582630933,1678.1441041335,1942.1632221482,2213.2407713351,2491.8242119209,3073.6699679147,
		3691.7716591048,4350.8611958576,5056.6653728972,5819.2277164608,6649.4016027568,7560.0781719656,
		8569.3345586917,9701.5591376725,10992.8973225474,11714.3712557337,12499.8515755171,13360.3068109264,
		14318.2241373647,15405.7570490178,16685.1606511117,18733.0547480514,20752.9602378743,23947.3420385707,
		26561.1613501294,31151.8561242708,33573.3309380886,35905.4984772108,39575.7896964869,42588.1612632421,
		47859.5785629266};
		
		
	vert_interpolate_1d(size,NZ,&z_input[0],&ZU(0),&tb_input[0],tb_out);
	vert_interpolate_1d(size,NZ,&z_input[0],&ZU(0),&qb_input[0],qb_out);
	
	tb_out[0]=tb_out[1];
	tb_out[NZ-1]=tb_out[NZ-2];

	qb_out[0]=qb_out[1];
	qb_out[NZ-1]=qb_out[NZ-2];
}