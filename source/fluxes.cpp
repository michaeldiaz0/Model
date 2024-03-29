#include "stdafx.h"
#include "fluxes.h"

/*******************************************************************************************
* A SET OF SUBROUTINES TO INTERPOLATE METEOROLOGICAL FIELDS FROM THEIR ORIGINAL LOCATIONS 
* ON THE MODEL'S STAGGERED GRID TO LOCATIONS TO THE FACES OF EACH CONTROL VOLUME AND TO USE THESE 
* VALUES TO CALCULATE THE FLUXES ON THE FACES OF THE CONTROL VOLUMES.
*
*
*
*
********************************************************************************************/

/***********************************************************************************************
* -------------------------------------- MACROS ------------------------------------------------
************************************************************************************************/
#if !PARALLEL
	#define SP(i,j,k) s[NZ*(NY*(i)+(j))+(k)]
	#define SB(i,j,k) sb[NZ*(NY*(i)+(j))+(k)]
#else
	#define SP(i,j,k) s[fNZ*(fNY*(i)+(j))+(k)]
	#define SB(i,j,k) sb[fNZ*(fNY*(i)+(j))+(k)]
#endif

/***********************************************************************************************
* CHOOSE A HORIZONTAL INTERPOLATION ORDER AT COMPILE TIME
************************************************************************************************/
#if HOR_ADVECTION_ORDER == 6

	#define INTERP_EAST(  VAR,SIGN,i) INTERP_6TH_EAST( VAR,i)
	#define INTERP_WEST(  VAR,SIGN,i) INTERP_6TH_WEST( VAR,i)
	#define INTERP_NORTH( VAR,SIGN,j) INTERP_6TH_NORTH(VAR,j)
	#define INTERP_SOUTH( VAR,SIGN,j) INTERP_6TH_SOUTH(VAR,j)

#elif HOR_ADVECTION_ORDER == 5

	#define INTERP_EAST(  VAR,SIGN,i) INTERP_5TH_EAST( VAR,SIGN,i)
	#define INTERP_WEST(  VAR,SIGN,i) INTERP_5TH_WEST( VAR,SIGN,i)
	#define INTERP_NORTH( VAR,SIGN,j) INTERP_5TH_NORTH(VAR,SIGN,j)
	#define INTERP_SOUTH( VAR,SIGN,j) INTERP_5TH_SOUTH(VAR,SIGN,j)

#elif HOR_ADVECTION_ORDER == 4

	#define INTERP_EAST(  VAR,SIGN,i) INTERP_4TH_EAST( VAR,i)
	#define INTERP_WEST(  VAR,SIGN,i) INTERP_4TH_WEST( VAR,i)
	#define INTERP_NORTH( VAR,SIGN,j) INTERP_4TH_NORTH(VAR,j)
	#define INTERP_SOUTH( VAR,SIGN,j) INTERP_4TH_SOUTH(VAR,j)

#elif HOR_ADVECTION_ORDER == 3

	#define INTERP_EAST(  VAR,SIGN,i) INTERP_3RD_EAST( VAR,SIGN,i)
	#define INTERP_WEST(  VAR,SIGN,i) INTERP_3RD_WEST( VAR,SIGN,i)
	#define INTERP_NORTH( VAR,SIGN,j) INTERP_3RD_NORTH(VAR,SIGN,j)
	#define INTERP_SOUTH( VAR,SIGN,j) INTERP_3RD_SOUTH(VAR,SIGN,j)

#elif HOR_ADVECTION_ORDER == 2

	#define INTERP_EAST(  VAR,SIGN,i) INTERP_2ND_EAST( VAR,i)
	#define INTERP_WEST(  VAR,SIGN,i) INTERP_2ND_WEST( VAR,i)
	#define INTERP_NORTH( VAR,SIGN,j) INTERP_2ND_NORTH(VAR,j)
	#define INTERP_SOUTH( VAR,SIGN,j) INTERP_2ND_SOUTH(VAR,j)

//#elif HOR_ADVECTION_ORDER == 1

//	#define INTERP_EAST(  VAR,SIGN,i) INTERP_1ST_EAST( VAR,SIGN,i)
//	#define INTERP_WEST(  VAR,SIGN,i) INTERP_1ST_WEST( VAR,SIGN,i)
//	#define INTERP_NORTH( VAR,SIGN,j) INTERP_1ST_NORTH(VAR,SIGN,j)
//	#define INTERP_SOUTH( VAR,SIGN,j) INTERP_1ST_SOUTH(VAR,SIGN,j)

#endif

/***********************************************************************************************
* CHOOSE A VERTICAL INTERPOLATION ORDER AT COMPILE TIME
************************************************************************************************/
#if VER_ADVECTION_ORDER == 6

	#define INTERP_TOP(   VAR,SIGN,k) INTERP_6TH_TOP(   VAR,k)
	#define INTERP_BOTTOM(VAR,SIGN,k) INTERP_6TH_BOTTOM(VAR,k)

#elif VER_ADVECTION_ORDER == 5

	#define INTERP_TOP(   VAR,SIGN,k) INTERP_5TH_TOP(   VAR,SIGN,k)
	#define INTERP_BOTTOM(VAR,SIGN,k) INTERP_5TH_BOTTOM(VAR,SIGN,k)

#elif VER_ADVECTION_ORDER == 4

	#define INTERP_TOP(   VAR,SIGN,k) INTERP_4TH_TOP(   VAR,k)
	#define INTERP_BOTTOM(VAR,SIGN,k) INTERP_4TH_BOTTOM(VAR,k)

#elif VER_ADVECTION_ORDER == 3

	#define INTERP_TOP(   VAR,SIGN,k) INTERP_3RD_TOP(   VAR,SIGN,k)
	#define INTERP_BOTTOM(VAR,SIGN,k) INTERP_3RD_BOTTOM(VAR,SIGN,k)

#elif VER_ADVECTION_ORDER == 2

	#define INTERP_TOP(   VAR,SIGN,k) INTERP_2ND_TOP(   VAR,k)
	#define INTERP_BOTTOM(VAR,SIGN,k) INTERP_2ND_BOTTOM(VAR,k)

//#elif VER_ADVECTION_ORDER == 1

//	#define INTERP_TOP(   VAR,SIGN,k) INTERP_1ST_TOP(   VAR,SIGN,k)
//	#define INTERP_BOTTOM(VAR,SIGN,k) INTERP_1ST_BOTTOM(VAR,SIGN,k)

#endif

/***********************************************************************************************
* INTERPOLATION OPERATORS
************************************************************************************************/
#define INTERP_6TH_EAST(  VAR,i) i_interp6th(VAR,i+1,i,i+2,i-1,i+3,i-2)
#define INTERP_6TH_WEST(  VAR,i) i_interp6th(VAR,i-1,i,i-2,i+1,i-3,i+2)
#define INTERP_6TH_NORTH( VAR,j) j_interp6th(VAR,j+1,j,j+2,j-1,j+3,j-2)
#define INTERP_6TH_SOUTH( VAR,j) j_interp6th(VAR,j-1,j,j-2,j+1,j-3,j+2)
#define INTERP_6TH_TOP(   VAR,k) k_interp6th(VAR,k+1,k,k+2,k-1,k+3,k-2)
#define INTERP_6TH_BOTTOM(VAR,k) k_interp6th(VAR,k-1,k,k-2,k+1,k-3,k+2)

#define INTERP_5TH_EAST(  VAR,SIGN,i) (INTERP_6TH_EAST(  VAR,i) + one60ths * (SIGN) * i_interp5th(VAR,i+3,i-2,i+2,i-1,i+1,i))
#define INTERP_5TH_WEST(  VAR,SIGN,i) (INTERP_6TH_WEST(  VAR,i) + one60ths * (SIGN) * i_interp5th(VAR,i+2,i-3,i+1,i-2,i,i-1))
#define INTERP_5TH_NORTH( VAR,SIGN,j) (INTERP_6TH_NORTH( VAR,j) + one60ths * (SIGN) * j_interp5th(VAR,j+3,j-2,j+2,j-1,j+1,j))
#define INTERP_5TH_SOUTH( VAR,SIGN,j) (INTERP_6TH_SOUTH( VAR,j) + one60ths * (SIGN) * j_interp5th(VAR,j+2,j-3,j+1,j-2,j,j-1))
#define INTERP_5TH_TOP(   VAR,SIGN,k) (INTERP_6TH_TOP(   VAR,k) + one60ths * (SIGN) * k_interp5th(VAR,k+3,k-2,k+2,k-1,k+1,k))
#define INTERP_5TH_BOTTOM(VAR,SIGN,k) (INTERP_6TH_BOTTOM(VAR,k) + one60ths * (SIGN) * k_interp5th(VAR,k+2,k-3,k+1,k-2,k,k-1))

#define INTERP_4TH_EAST(  VAR,i) i_interp4th(VAR,i+1,i,i+2,i-1)
#define INTERP_4TH_WEST(  VAR,i) i_interp4th(VAR,i-1,i,i-2,i+1)
#define INTERP_4TH_NORTH( VAR,j) j_interp4th(VAR,j+1,j,j+2,j-1)
#define INTERP_4TH_SOUTH( VAR,j) j_interp4th(VAR,j-1,j,j-2,j+1)
#define INTERP_4TH_TOP(   VAR,k) k_interp4th(VAR,k+1,k,k+2,k-1)
#define INTERP_4TH_BOTTOM(VAR,k) k_interp4th(VAR,k-1,k,k-2,k+1)

#define INTERP_3RD_EAST(  VAR,SIGN,i) (INTERP_4TH_EAST(  VAR,i) + onetwelfth * (SIGN) * i_interp3rd(VAR,i+2,i-1,i+1,i))
#define INTERP_3RD_WEST(  VAR,SIGN,i) (INTERP_4TH_WEST(  VAR,i) + onetwelfth * (SIGN) * i_interp3rd(VAR,i+1,i-2,i,i-1))
#define INTERP_3RD_NORTH( VAR,SIGN,j) (INTERP_4TH_NORTH( VAR,j) + onetwelfth * (SIGN) * j_interp3rd(VAR,j+2,j-1,j+1,j))
#define INTERP_3RD_SOUTH( VAR,SIGN,j) (INTERP_4TH_SOUTH( VAR,j) + onetwelfth * (SIGN) * j_interp3rd(VAR,j+1,j-2,j,j-1))
#define INTERP_3RD_TOP(   VAR,SIGN,k) (INTERP_4TH_TOP(   VAR,k) + onetwelfth * (SIGN) * k_interp3rd(VAR,k+2,k-1,k+1,k))
#define INTERP_3RD_BOTTOM(VAR,SIGN,k) (INTERP_4TH_BOTTOM(VAR,k) + onetwelfth * (SIGN) * k_interp3rd(VAR,k+1,k-2,k,k-1))

#define INTERP_2ND_EAST(  VAR,i) i_interp2nd(VAR,i+1,i)
#define INTERP_2ND_WEST(  VAR,i) i_interp2nd(VAR,i-1,i)
#define INTERP_2ND_NORTH( VAR,j) j_interp2nd(VAR,j+1,j)
#define INTERP_2ND_SOUTH( VAR,j) j_interp2nd(VAR,j-1,j)
#define INTERP_2ND_TOP(   VAR,k) k_interp2nd(VAR,k+1,k)
#define INTERP_2ND_BOTTOM(VAR,k) k_interp2nd(VAR,k-1,k)

#define INTERP_1ST_EAST(  VAR,SIGN,i) adv_1st(SIGN,VAR(i  ,j  ,k  ),VAR(i+1,j,  k  ))
#define INTERP_1ST_WEST(  VAR,SIGN,i) adv_1st(SIGN,VAR(i-1,j  ,k  ),VAR(i  ,j,  k  ))
#define INTERP_1ST_NORTH( VAR,SIGN,j) adv_1st(SIGN,VAR(i  ,j  ,k  ),VAR(i  ,j+1,k  ))
#define INTERP_1ST_SOUTH( VAR,SIGN,j) adv_1st(SIGN,VAR(i  ,j-1,k  ),VAR(i  ,j  ,k  ))
#define INTERP_1ST_TOP(   VAR,SIGN,k) adv_1st(SIGN,VAR(i  ,j  ,k  ),VAR(i  ,j  ,k+1))
#define INTERP_1ST_BOTTOM(VAR,SIGN,k) adv_1st(SIGN,VAR(i  ,j  ,k-1),VAR(i  ,j  ,k  ))

/***********************************************************************************************
* SIXTH ORDER INTERPOLATION FOR ADVECTION
************************************************************************************************/
#define i_interp6th(s,i0,i1,i2,i3,i4,i5) (  thirtyseven60ths * ( s(i0,j,k)+s(i1,j,k) ) - \
													two15ths * ( s(i2,j,k)+s(i3,j,k) ) + \
													one60ths * ( s(i4,j,k)+s(i5,j,k) )   )

#define j_interp6th(s,j0,j1,j2,j3,j4,j5) (  thirtyseven60ths * ( s(i,j0,k)+s(i,j1,k) ) - \
													two15ths * ( s(i,j2,k)+s(i,j3,k) ) + \
													one60ths * ( s(i,j4,k)+s(i,j5,k) )   )

#define k_interp6th(s,k0,k1,k2,k3,k4,k5) (  thirtyseven60ths * ( s(i,j,k0)+s(i,j,k1) ) - \
													two15ths * ( s(i,j,k2)+s(i,j,k3) ) + \
													one60ths * ( s(i,j,k4)+s(i,j,k5) )   )

/***********************************************************************************************
* FIFTH ORDER INTERPOLATION FOR ADVECTION (JUST THE UPWIND BIASED PART)
************************************************************************************************/
#define i_interp5th(s,i0,i1,i2,i3,i4,i5) (  -( s(i0,j,k)-s(i1,j,k) ) + \
										 5.0*( s(i2,j,k)-s(i3,j,k) ) - \
										10.0*( s(i4,j,k)-s(i5,j,k) )  )

#define j_interp5th(s,j0,j1,j2,j3,j4,j5) (  -( s(i,j0,k)-s(i,j1,k) ) + \
										 5.0*( s(i,j2,k)-s(i,j3,k) ) - \
										10.0*( s(i,j4,k)-s(i,j5,k) )  )

#define k_interp5th(s,k0,k1,k2,k3,k4,k5) (  -( s(i,j,k0)-s(i,j,k1) ) + \
										 5.0*( s(i,j,k2)-s(i,j,k3) ) - \
										10.0*( s(i,j,k4)-s(i,j,k5) )  )

/***********************************************************************************************
* FOURTH ORDER INTERPOLATION FOR ADVECTION
************************************************************************************************/
#define i_interp4th(s,i0,i1,i2,i3) (  seventwelfth * ( s(i0,j,k) + s(i1,j,k) ) - \
										onetwelfth * ( s(i2,j,k) + s(i3,j,k) )  )

#define j_interp4th(s,j0,j1,j2,j3) (  seventwelfth * ( s(i,j0,k) + s(i,j1,k) ) - \
										onetwelfth * ( s(i,j2,k) + s(i,j3,k) )  )

#define k_interp4th(s,k0,k1,k2,k3) (  seventwelfth * ( s(i,j,k0) + s(i,j,k1) ) - \
										onetwelfth * ( s(i,j,k2) + s(i,j,k3) )  )

/***********************************************************************************************
* THIRD ORDER INTERPOLATION FOR ADVECTION (JUST THE UPWIND BIASED PART)
************************************************************************************************/
#define i_interp3rd(s,i0,i1,i2,i3) (  ( s(i0,j,k)-s(i1,j,k) ) - 3.0*( s(i2,j,k)-s(i3,j,k) ) )
#define j_interp3rd(s,j0,j1,j2,j3) (  ( s(i,j0,k)-s(i,j1,k) ) - 3.0*( s(i,j2,k)-s(i,j3,k) ) )
#define k_interp3rd(s,k0,k1,k2,k3) (  ( s(i,j,k0)-s(i,j,k1) ) - 3.0*( s(i,j,k2)-s(i,j,k3) ) )

/***********************************************************************************************
* SECOND ORDER INTERPOLATION FOR ADVECTION
************************************************************************************************/
#define i_interp2nd(s,i0,i1) (  0.5* (s(i0,j,k) + s(i1,j,k) ) )
#define j_interp2nd(s,j0,j1) (  0.5* (s(i,j0,k) + s(i,j1,k) ) )
#define k_interp2nd(s,k0,k1) (  0.5* (s(i,j,k0) + s(i,j,k1) ) )

/***********************************************************************************************
* FIRST ORDER INTERPOLATION FOR ADVECTION
************************************************************************************************/
//#define i_interp1st(s,i0) (s(i0,j,k)) 
//#define j_interp1st(s,j0) (s(i,j0,k))
//#define k_interp1st(s,k0) (s(i,j,k0))

/***************************************************************************
* -------------------------- VARIABLES -------------------------------------
****************************************************************************/
const double onetwelfth = 1.0/12.0;
const double seventwelfth = 7.0/12.0;
const double thirtyseven60ths = 37.0/60.0;
const double two15ths = 2.0/15.0;
const double one60ths = 1.0/60.0;
														 
struct vel_cell *ucell,*vcell,*wcell;
struct vel_cell *ubcell,*vbcell,*wbcell;
struct cell *thcell,*thbcell;
struct cell *qvcell,*qvbcell,*qccell,*qrcell,*ice_cell;
struct sign_vel *sign_p_vel,*sign_b_vel;

/***************************************************************************
* ---------------------- FUNCTION PROTOTYPES--------------------------------
****************************************************************************/
void interpolate_velocity(int i,int jl,int jh);
void interpolate_scalar(int i,int jl,int jh,double *u,double *v,double *w,double *s,struct cell *scell);
void interpolate_scalar(int i,int jl,int jh,double *s,double *sb,struct cell *scell,struct cell *bcell);
void interpolate_scalar(int i,int jl,int jh,double *s,struct cell *scell);
void interpolate_scalar_with_fallspeed(int i,int jl,int jh,double *s,double *fall,struct cell *scell);
void interpolate_moisture(int i,int jl,int jh);
inline static char signof(double x);

/*********************************************************************
* 
**********************************************************************/
void initialize_flux_cells(int ny,int nz){

	ucell = (vel_cell*)calloc(ny*nz,sizeof(vel_cell));
	vcell = (vel_cell*)calloc(ny*nz,sizeof(vel_cell));
	wcell = (vel_cell*)calloc(ny*nz,sizeof(vel_cell));

	ubcell = (vel_cell*)calloc(ny*nz,sizeof(vel_cell));
	vbcell = (vel_cell*)calloc(ny*nz,sizeof(vel_cell));
	wbcell = (vel_cell*)calloc(ny*nz,sizeof(vel_cell));

	thcell  = (cell*)calloc(ny*nz,sizeof(cell));
	thbcell = (cell*)calloc(ny*nz,sizeof(cell));	
}

/*********************************************************************
* 
**********************************************************************/
void initialize_microphysics_cells(int ny,int nz){

	qvcell  = (cell*)malloc(ny*nz*sizeof(cell));
	qvbcell = (cell*)malloc(ny*nz*sizeof(cell));
	qccell  = (cell*)malloc(ny*nz*sizeof(cell));
	qrcell  = (cell*)malloc(ny*nz*sizeof(cell));
	ice_cell  = (cell*)calloc(ny*nz,sizeof(cell));
}

/*********************************************************************
* 
**********************************************************************/
void initialize_sign_cells(int nx,int ny,int nz){
	
	sign_p_vel = (sign_vel*)malloc(nx*ny*nz*sizeof(sign_vel));
	sign_b_vel = (sign_vel*)malloc(nx*ny*nz*sizeof(sign_vel));
	
}

/*********************************************************************
* Return the sign of x multiplied by one
*
**********************************************************************/
/*inline double adv_1st(double sign,double y0,double y1){
	
	if(sign>0){ return y1;}
	
	return y0;
}*/

/*********************************************************************
* Return the sign of x multiplied by one
*
**********************************************************************/
inline static char signof(double x){
	
	if(x>0){ return 1;}
	
	return -1;
}

/*********************************************************************
* Return the sign of x multiplied by one
*
**********************************************************************/
inline static double signof_d(double x){
	
	if(x>0){ return 1;}
	
	return -1;
}

/*********************************************************************
* Calculate sign of advecting velocity
**********************************************************************/
void compute_sign_cells(int il,int ih,int jl,int jh){

#if HOR_ADVECTION_ORDER % 2 != 0

	for(int i=il;i<ih;i++){	
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		SIGN_P_VEL(i,j,k).u = signof(U(i,j,k));
		SIGN_P_VEL(i,j,k).v = signof(V(i,j,k));
		SIGN_P_VEL(i,j,k).w = signof(W(i,j,k));

	#if !ISLINEAR
		SIGN_B_VEL(i,j,k).u = signof(UBAR(i,j,k)+U(i,j,k));
		SIGN_B_VEL(i,j,k).v = signof(VBAR(i,j,k)+V(i,j,k));
		SIGN_B_VEL(i,j,k).w = signof(WBAR(i,j,k)+W(i,j,k));
	#else
		SIGN_B_VEL(i,j,k).u = signof(UBAR(i,j,k));
		SIGN_B_VEL(i,j,k).v = signof(VBAR(i,j,k));
		SIGN_B_VEL(i,j,k).w = signof(WBAR(i,j,k));
	#endif
	
	}}}
#endif
	
}

/*********************************************************************
* Initialize flux values on the eastern side of the model domain.
*
* @param i - the x-coordinate
* @param jl,jh - the high and low index bounds for the y-coordinate
**********************************************************************/
void compute_west(int i,int jl,int jh){

	double ut,vt,wt;

	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		/***************************************************
		* Calculate advecting velocity
		****************************************************/
		UCELL(j,k).ua = 0.5*(U(i+1,j,k)+U(i,j,k));
		VCELL(j,k).ua = 0.5*(U(i+1,j,k)+U(i+1,j-1,k));
		WCELL(j,k).ua = 0.5*(U(i+1,j,k-1)+U(i+1,j,k));

		UBCELL(j,k).ua = 0.5*(UBAR(i+1,j,k)+UBAR(i,j,k));
		VBCELL(j,k).ua = 0.5*(UBAR(i+1,j,k)+UBAR(i+1,j-1,k));
		WBCELL(j,k).ua = 0.5*(UBAR(i+1,j,k-1)+UBAR(i+1,j,k));
		
		// sign of the advecting velocity
	#if !ISLINEAR
		ut = signof(UBCELL(j,k).ua+UCELL(j,k).ua);
		vt = signof(VBCELL(j,k).ua+VCELL(j,k).ua);
		wt = signof(WBCELL(j,k).ua+WCELL(j,k).ua);
	#else
		ut = signof(UBCELL(j,k).ua);
		vt = signof(VBCELL(j,k).ua);
		wt = signof(WBCELL(j,k).ua);
	#endif
		
		/***************************************************
		* Fifth order interpolation for horizontal derivatives
		* for advected velocity
		****************************************************/
		UCELL(j,k).east  = INTERP_EAST(U,ut,i);
		VCELL(j,k).east  = INTERP_EAST(V,vt,i);
		WCELL(j,k).east  = INTERP_EAST(W,wt,i);
		
		UBCELL(j,k).east  = INTERP_EAST(UBAR,signof(UCELL(j,k).ua),i);
		VBCELL(j,k).east  = INTERP_EAST(VBAR,signof(VCELL(j,k).ua),i);
		WBCELL(j,k).east  = INTERP_EAST(WBAR,signof(WCELL(j,k).ua),i);

		// advecting velocity for perturbation
	#if !ISLINEAR
		ut = UBCELL(j,k).ua+UCELL(j,k).ua;
		vt = VBCELL(j,k).ua+VCELL(j,k).ua;
		wt = WBCELL(j,k).ua+WCELL(j,k).ua;
	#else
		ut = UBCELL(j,k).ua;
		vt = VBCELL(j,k).ua;
		wt = WBCELL(j,k).ua;
	#endif

		/***************************************************
		* Calculate velocity fluxes
		****************************************************/
		UCELL(j,k).east  = UCELL(j,k).ua * UBCELL(j,k).east  + ut * UCELL(j,k).east;
		VCELL(j,k).east  = VCELL(j,k).ua * VBCELL(j,k).east  + vt * VCELL(j,k).east;
		WCELL(j,k).east  = WCELL(j,k).ua * WBCELL(j,k).east  + wt * WCELL(j,k).east;
	}}
}

/*********************************************************************
* Initialize flux values on the eastern side
*
* @param i - the x-coordinate
* @param jl,jh - the high and low index bounds for the y-coordinate
**********************************************************************/
void compute_west_uv(int i,int jl,int jh){

	double ut,vt;

	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		/***************************************************
		* Calculate advecting velocity
		****************************************************/
		UCELL(j,k).ua = 0.5*(U(i+1,j,k)+U(i,j,k));
		VCELL(j,k).ua = 0.5*(U(i+1,j,k)+U(i+1,j-1,k));

		UBCELL(j,k).ua = 0.5*(UBAR(i+1,j,k)+UBAR(i,j,k));
		VBCELL(j,k).ua = 0.5*(UBAR(i+1,j,k)+UBAR(i+1,j-1,k));
		
		// sign of the advecting velocity
	#if !ISLINEAR
		ut = signof(UBCELL(j,k).ua+UCELL(j,k).ua);
		vt = signof(VBCELL(j,k).ua+VCELL(j,k).ua);
	#else
		ut = signof(UBCELL(j,k).ua);
		vt = signof(VBCELL(j,k).ua);
	#endif
		
		/***************************************************
		* Fifth order interpolation for horizontal derivatives
		* for advected velocity
		****************************************************/
		UCELL(j,k).east  = INTERP_EAST(U,ut,i);
		VCELL(j,k).east  = INTERP_EAST(V,vt,i);
		
		UBCELL(j,k).east  = INTERP_EAST(UBAR,signof(UCELL(j,k).ua),i);
		VBCELL(j,k).east  = INTERP_EAST(VBAR,signof(VCELL(j,k).ua),i);

		// advecting velocity for perturbation
	#if !ISLINEAR
		ut = UBCELL(j,k).ua+UCELL(j,k).ua;
		vt = VBCELL(j,k).ua+VCELL(j,k).ua;
	#else
		ut = UBCELL(j,k).ua;
		vt = VBCELL(j,k).ua;
	#endif

		/***************************************************
		* Calculate velocity fluxes
		****************************************************/
		UCELL(j,k).east  = UCELL(j,k).ua * UBCELL(j,k).east  + ut * UCELL(j,k).east;
		VCELL(j,k).east  = VCELL(j,k).ua * VBCELL(j,k).east  + vt * VCELL(j,k).east;
	}}
}

/*********************************************************************
* Initialize flux values on the eastern side of a YZ cross section
* of control volumes.
*
* @param i - the x-coordinate
* @param jl,jh - the high and low index bounds for the y-coordinate
* @param s - intput perturbation scalar field
* @param sb - input background scalar field
* @param scell - output interpolated YZ field for perturbations
* @param bcell - output interpolated YZ field for background
**********************************************************************/
void compute_west_scalar(int i,int jl,int jh,double *s,double *sb,struct cell *scell,struct cell *bcell){

	double sign;

	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

	#if !ISLINEAR
		sign = signof(UBAR(i+1,j,k)+U(i+1,j,k));
	#else
		sign = signof(UBAR(i+1,j,k));
	#endif

		SCELL(j,k).east = INTERP_EAST(SP,sign,i);
		BCELL(j,k).east = INTERP_EAST(SB,signof(U(i+1,j,k)),i);
		
	#if !ISLINEAR
		SCELL(j,k).east  = U(i+1,j,k) * BCELL(j,k).east  + (UBAR(i+1,j,k)+U(i+1,j,k)) * SCELL(j,k).east;
	#else
		SCELL(j,k).east  = U(i+1,j,k) * BCELL(j,k).east  + UBAR(i+1,j,k) * SCELL(j,k).east;
	#endif
	}}
}

/*********************************************************************
* Initialize flux values on the eastern side of a YZ cross section
* of control volumes. Without basic state.
*
* @param i - the x-coordinate
* @param jl,jh - the high and low index bounds for the y-coordinate
* @param s - intput perturbation scalar field
* @param scell - output interpolated YZ field for perturbations
**********************************************************************/
void compute_west_scalar(int i,int jl,int jh,double *s,struct cell *scell){

	double sign;

	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		sign = signof(UBAR(i+1,j,k)+U(i+1,j,k));

		SCELL(j,k).east = INTERP_EAST(SP,sign,i);
		
		SCELL(j,k).east  = (UBAR(i+1,j,k)+U(i+1,j,k)) * SCELL(j,k).east;
	}}
}

/*********************************************************************
* Initialize flux values on the eastern side of a YZ cross section
* of control volumes. Without basic state.
*
* @param i - the x-coordinate
* @param jl,jh - the high and low index bounds for the y-coordinate
* @param s - intput perturbation scalar field
* @param scell - output interpolated YZ field for perturbations
**********************************************************************/
void compute_west_scalar(int i,int jl,int jh,double *u,double *s,struct cell *scell){

	double sign;

	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		sign = signof(u[INDEX(i,j,k)]);

		SCELL(j,k).east = INTERP_EAST(SP,sign,i);
		
		SCELL(j,k).east  = u[INDEX(i,j,k)] * SCELL(j,k).east;
	}}
}

/*********************************************************************
* Initialize flux values on the eastern side
**********************************************************************/
void compute_west_moisture(int i,int jl,int jh){

	double ub,usign,ubsign;

	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){

		ub = UBAR(i+1,j,k)+U(i+1,j,k);

		usign = signof(U(i+1,j,k));
		ubsign = signof(UBAR(i+1,j,k)+U(i+1,j,k));
		
		// interpolate moisture variables to the eastern face of the control volume
		QVCELL(j,k).east  = INTERP_EAST(QV,ubsign,i);
		QCCELL(j,k).east  = INTERP_EAST(QC,ubsign,i);
		QRCELL(j,k).east  = INTERP_EAST(QR,ubsign,i);
		QVBCELL(j,k).east = INTERP_EAST(QBAR,usign,i);
	
		// calculate fluxes
		QVCELL(j,k).east  = ub * QVCELL(j,k).east +  U(i+1,j,k) * QVBCELL(j,k).east;
		QCCELL(j,k).east  = ub * QCCELL(j,k).east;
		QRCELL(j,k).east  = ub * QRCELL(j,k).east;
	}}
}

/*********************************************************************
* Compute velocity fluxes
**********************************************************************/
void compute_fluxes(int i,int jl,int jh){

	interpolate_velocity(i,jl,jh);

	for(int j=jl;j<jh;j++){

		UCELL(j,0).top = 0;
		VCELL(j,0).top = 0;
		WCELL(j,0).top = 0;

		for(int k=1;k<NZ-1;k++){

	#if !ISLINEAR
			UCELL(j,k).east  = UCELL(j,k).ua * UBCELL(j,k).east  + (UBCELL(j,k).ua + UCELL(j,k).ua) * UCELL(j,k).east;
			UCELL(j,k).north = UCELL(j,k).va * UBCELL(j,k).north + (UBCELL(j,k).va + UCELL(j,k).va) * UCELL(j,k).north;
			UCELL(j,k).top   = UCELL(j,k).wa * UBCELL(j,k).top   + (UBCELL(j,k).wa + UCELL(j,k).wa) * UCELL(j,k).top;

			VCELL(j,k).east  = VCELL(j,k).ua * VBCELL(j,k).east  + (VBCELL(j,k).ua + VCELL(j,k).ua) * VCELL(j,k).east;
			VCELL(j,k).north = VCELL(j,k).va * VBCELL(j,k).north + (VBCELL(j,k).va + VCELL(j,k).va) * VCELL(j,k).north;
			VCELL(j,k).top   = VCELL(j,k).wa * VBCELL(j,k).top   + (VBCELL(j,k).wa + VCELL(j,k).wa) * VCELL(j,k).top;
		#if !HYDROSTATIC
			WCELL(j,k).east  = WCELL(j,k).ua * WBCELL(j,k).east  + (WBCELL(j,k).ua + WCELL(j,k).ua) * WCELL(j,k).east;
			WCELL(j,k).north = WCELL(j,k).va * WBCELL(j,k).north + (WBCELL(j,k).va + WCELL(j,k).va) * WCELL(j,k).north;
			WCELL(j,k).top   = WCELL(j,k).wa * WBCELL(j,k).top   + (WBCELL(j,k).wa + WCELL(j,k).wa) * WCELL(j,k).top;
		#endif
			
	#else
			UCELL(j,k).east  = UCELL(j,k).ua * UBCELL(j,k).east  + UBCELL(j,k).ua * UCELL(j,k).east;
			UCELL(j,k).north = UCELL(j,k).va * UBCELL(j,k).north + UBCELL(j,k).va * UCELL(j,k).north;
			UCELL(j,k).top   = UCELL(j,k).wa * UBCELL(j,k).top   + UBCELL(j,k).wa * UCELL(j,k).top;

			VCELL(j,k).east  = VCELL(j,k).ua * VBCELL(j,k).east  + VBCELL(j,k).ua * VCELL(j,k).east;
			VCELL(j,k).north = VCELL(j,k).va * VBCELL(j,k).north + VBCELL(j,k).va * VCELL(j,k).north;
			VCELL(j,k).top   = VCELL(j,k).wa * VBCELL(j,k).top   + VBCELL(j,k).wa * VCELL(j,k).top;
		#if !HYDROSTATIC
			WCELL(j,k).east  = WCELL(j,k).ua * WBCELL(j,k).east  + WBCELL(j,k).ua * WCELL(j,k).east;
			WCELL(j,k).north = WCELL(j,k).va * WBCELL(j,k).north + WBCELL(j,k).va * WCELL(j,k).north;
			WCELL(j,k).top   = WCELL(j,k).wa * WBCELL(j,k).top   + WBCELL(j,k).wa * WCELL(j,k).top;
		#endif
			
	#endif
			
		}
	}

}

/*********************************************************************
* Compute a scalar flux
**********************************************************************/
void compute_fluxes_scalar(int i,int jl,int jh,double *s,double *sb,struct cell *scell,struct cell *bcell){

	interpolate_scalar(i,jl,jh,s,sb,scell,bcell);

	for(int j=jl;j<jh;j++){

		SCELL(j,0).top = 0;

		for(int k=1;k<NZ-1;k++){

		#if !ISLINEAR
			SCELL(j,k).east  = U(i+1,j,k) * BCELL(j,k).east  + (UBAR(i+1,j,k)+U(i+1,j,k)) * SCELL(j,k).east;
			SCELL(j,k).north = V(i,j+1,k) * BCELL(j,k).north + (VBAR(i,j+1,k)+V(i,j+1,k)) * SCELL(j,k).north;
			SCELL(j,k).top   = W(i,j,k+1) * BCELL(j,k).top   + (WBAR(i,j,k+1)+W(i,j,k+1)) * SCELL(j,k).top;
		#else
			SCELL(j,k).east  = U(i+1,j,k) * BCELL(j,k).east  + UBAR(i+1,j,k) * SCELL(j,k).east;
			SCELL(j,k).north = V(i,j+1,k) * BCELL(j,k).north + VBAR(i,j+1,k) * SCELL(j,k).north;
			SCELL(j,k).top   = W(i,j,k+1) * BCELL(j,k).top   + WBAR(i,j,k+1) * SCELL(j,k).top;
		#endif
			
		}
	}
}

/*********************************************************************
* Compute a scalar flux (without basic state terms)
**********************************************************************/
void compute_fluxes_scalar(int i,int jl,int jh,double *s,struct cell *scell){

	interpolate_scalar(i,jl,jh,s,scell);

	for(int j=jl;j<jh;j++){

		SCELL(j,0).top = 0;

		for(int k=1;k<NZ-1;k++){

			SCELL(j,k).east  = (UBAR(i+1,j,k)+U(i+1,j,k)) * SCELL(j,k).east;
			SCELL(j,k).north = (VBAR(i,j+1,k)+V(i,j+1,k)) * SCELL(j,k).north;
			SCELL(j,k).top   = (WBAR(i,j,k+1)+W(i,j,k+1)) * SCELL(j,k).top;
		}
	}
}

/*********************************************************************
* Compute a scalar flux (without basic state terms)
**********************************************************************/
void compute_fluxes_scalar(int i,int jl,int jh,double *u,double *v,double *w,double *s,struct cell *scell){

	interpolate_scalar(i,jl,jh,u,v,w,s,scell);

	for(int j=jl;j<jh;j++){

		SCELL(j,0).top = 0;

		for(int k=1;k<NZ-1;k++){

			SCELL(j,k).east  = u[INDEX(i,j,k)] * SCELL(j,k).east;
			SCELL(j,k).north = v[INDEX(i,j,k)] * SCELL(j,k).north;
			SCELL(j,k).top   = w[INDEX(i,j,k)] * SCELL(j,k).top;
		}
	}
}

/*********************************************************************
* Compute a scalar flux with a fall speed (without basic state terms)
**********************************************************************/
void compute_fluxes_scalar_with_fallspeed(int i,int jl,int jh,double *s,double *fall,struct cell *scell){

	interpolate_scalar_with_fallspeed(i,jl,jh,s,fall,scell);

	for(int j=jl;j<jh;j++){

		SCELL(j,0).top = -0.5*(fall[INDEX(i,j,1)]+fall[INDEX(i,j,0)]) * SCELL(j,0).top;

		for(int k=1;k<NZ-1;k++){

			SCELL(j,k).east  = (UBAR(i+1,j,k)+U(i+1,j,k)) * SCELL(j,k).east;
			SCELL(j,k).north = (VBAR(i,j+1,k)+V(i,j+1,k)) * SCELL(j,k).north;
			SCELL(j,k).top   = (WBAR(i,j,k+1)+W(i,j,k+1)-0.5*(fall[INDEX(i,j,k)]+fall[INDEX(i,j,k+1)]) ) * SCELL(j,k).top;
		}
	}
}

/*********************************************************************
* Compute fluxes for moisture variables for a YZ cross section
**********************************************************************/
void compute_fluxes_moisture(int i,int jl,int jh){

	interpolate_scalar(i, jl, jh, qcs, qccell);
	interpolate_scalar_with_fallspeed(i, jl, jh, qrs, vts, qrcell);
	interpolate_scalar(i, jl, jh, qvs, m_qbar, qvcell, qvbcell);

	double ub,vb,wb;

	//---------------------------------------------------------
	// Zonal fluxes
	//---------------------------------------------------------
	for(int j=jl;j<jh;j++){

		QVCELL(j,0).top = 0;
		QCCELL(j,0).top = 0;
		QRCELL(j,0).top = -0.5*(VT(i,j,1)+VT(i,j,0)) * QRCELL(j,0).top;	// allows rain to fall out through the ground
		
		for(int k=1;k<NZ-1;k++){

			// advecting velocities
			ub = UBAR(i+1,j,k)+U(i+1,j,k);
			// vapor fluxes
			QVCELL(j,k).east  = ub * QVCELL(j,k).east  + U(i+1,j,k) * QVBCELL(j,k).east;
			// cloud water fluxes
			QCCELL(j,k).east  = ub * QCCELL(j,k).east;
			// rain water fluxes
			QRCELL(j,k).east  = ub * QRCELL(j,k).east;
		}
	}
	
	//---------------------------------------------------------
	// Meridional fluxes
	//---------------------------------------------------------
	for(int j=jl;j<jh;j++){

		for(int k=1;k<NZ-1;k++){
			// advecting velocities
			vb = VBAR(i,j+1,k)+V(i,j+1,k);
			// vapor fluxes
			QVCELL(j,k).north = vb * QVCELL(j,k).north + V(i,j+1,k) * QVBCELL(j,k).north;	
			// cloud water fluxes
			QCCELL(j,k).north = vb * QCCELL(j,k).north;
			// rain water fluxes
			QRCELL(j,k).north = vb * QRCELL(j,k).north;
		}
	}
	
	//---------------------------------------------------------
	// Vertical fluxes
	//---------------------------------------------------------
	for(int j=jl;j<jh;j++){

		for(int k=1;k<NZ-1;k++){

			// advecting velocities
			wb = WBAR(i,j,k+1)+W(i,j,k+1);
			// vapor fluxes
			QVCELL(j,k).top   = wb * QVCELL(j,k).top   + W(i,j,k+1) * QVBCELL(j,k).top;		
			// cloud water fluxes
			QCCELL(j,k).top   = wb * QCCELL(j,k).top;
			// rain water fluxes
			QRCELL(j,k).top = (wb-0.5*(VT(i,j,k+1)+VT(i,j,k))) * QRCELL(j,k).top;
		}
	}
	
}

/*********************************************************************
* Interpolate a scalar field to the faces of each control volume for
* a YZ cross section with a basic state.
*
* @param i - the x-coordinate
* @param jl,jh - the high and low index bounds for the y-coordinate
* @param s - intput perturbation scalar field
* @param sb - input background scalar field
* @param scell - output interpolated YZ field for perturbations
* @param bcell - output interpolated YZ field for background
**********************************************************************/
void interpolate_scalar(int i,int jl,int jh,double *s,double *sb,struct cell *scell,struct cell *bcell){

	int k,kmin,kmax;

	kmin = (VER_ADVECTION_ORDER+1) / 2;
	kmax = NZ - (VER_ADVECTION_ORDER+1) / 2;

	// interpolate perturbation scalar field
	interpolate_scalar(i, jl, jh, s, scell);

	//---------------------------------------------------
	// Horizontal part
	//---------------------------------------------------
	for(int j=jl;j<jh;j++){
	for(k=1;k<NZ-1;k++){

		//-------------------------------------------
		// Base state interpolation
		BCELL(j,k).west  = BCELL(j,k).east;
		BCELL(j,k).east  = INTERP_EAST( SB, SIGN_P_VEL(i+1,j,k).u, i);
		BCELL(j,k).north = INTERP_NORTH(SB, SIGN_P_VEL(i,j+1,k).v, j);
	}}

	//---------------------------------------------------
	// Upper and lower boundary points
	// Use lower order interpolations in vertical if 
	// the interpolating stencil extends beyond model domain
	//---------------------------------------------------
	#if VER_ADVECTION_ORDER > 2
	
	for(int j=jl;j<jh;j++){
		
		k = 1;
		
		BCELL(j,k).top = INTERP_2ND_TOP(SB,k);
		
		k = NZ-2;

		BCELL(j,k).top = INTERP_2ND_TOP(SB,k);
	}
	
	#endif
	
	#if VER_ADVECTION_ORDER > 4
	
	for(int j=jl;j<jh;j++){

		k = 2;

		BCELL(j,k).top = INTERP_3RD_TOP(SB, SIGN_P_VEL(i,j,k+1).w, k);
		
		k = NZ-3;

		BCELL(j,k).top = INTERP_3RD_TOP(SB, SIGN_P_VEL(i,j,k+1).w, k);
		
	}
		
	#endif

	//---------------------------------------------------
	// Interior points
	// Use requested interpolation order if the stencil
	// lies completely within the model domain
	//---------------------------------------------------
	for(int j=jl;j<jh;j++){
	for(k=kmin;k<kmax;k++){

		BCELL(j,k).top = INTERP_TOP(  SB, SIGN_P_VEL(i,j,k+1).w, k);
	}}

}

/*********************************************************************
* Interpolate a scalar field to the faces of each control volume for
* a YZ cross section without a basic state
*
* @param i - the x-coordinate
* @param jl,jh - the high and low index bounds for the y-coordinate
* @param s - intput perturbation scalar field
* @param scell - output interpolated YZ field for perturbations
**********************************************************************/
void interpolate_scalar(int i,int jl,int jh,double *s,struct cell *scell){

	int k,kmin,kmax;
	
	kmin = (VER_ADVECTION_ORDER+1) / 2;
	kmax = NZ - (VER_ADVECTION_ORDER+1) / 2;

	//---------------------------------------------------
	// Horizontal part
	//---------------------------------------------------
	for(int j=jl;j<jh;j++){
	for(k=1;k<NZ-1;k++){

		SCELL(j,k).west  = SCELL(j,k).east;
		SCELL(j,k).east  = INTERP_EAST( SP, SIGN_B_VEL(i+1,j,k).u, i);
		SCELL(j,k).north = INTERP_NORTH(SP, SIGN_B_VEL(i,j+1,k).v, j);
	}}

	//---------------------------------------------------
	// Upper and lower boundary points
	// Use lower order interpolations in vertical if 
	// the interpolating stencil extends beyond model domain
	//---------------------------------------------------
	#if VER_ADVECTION_ORDER > 2
	
	for(int j=jl;j<jh;j++){

		k = 1;

		SCELL(j,k).top = INTERP_2ND_TOP(SP,k);
		
		k = NZ-2;

		SCELL(j,k).top = INTERP_2ND_TOP(SP,k);
	}
		
	#endif
		
	#if VER_ADVECTION_ORDER > 4
	
	for(int j=jl;j<jh;j++){
		
		k = 2;
		
		SCELL(j,k).top = INTERP_3RD_TOP(SP, SIGN_B_VEL(i,j,k+1).w, k);
		
		k = NZ-3;

		SCELL(j,k).top = INTERP_3RD_TOP(SP, SIGN_B_VEL(i,j,k+1).w, k);
	}
		
	#endif

	//---------------------------------------------------
	// Interior points
	// Use requested interpolation order if the stencil
	// lies completely within the model domain
	//---------------------------------------------------
	for(int j=jl;j<jh;j++){
	for(k=kmin;k<kmax;k++){

		SCELL(j,k).top = INTERP_TOP(SP, SIGN_B_VEL(i,j,k+1).w, k);
	}}

}

/*********************************************************************
* Interpolate a scalar field to the faces of each control volume for
* a YZ cross section. Without basic state
*
* @param i - the x-coordinate
* @param jl,jh - the high and low index bounds for the y-coordinate
* @param s - intput perturbation scalar field
* @param scell - output interpolated YZ field for perturbations
**********************************************************************/
void interpolate_scalar_with_fallspeed(int i,int jl,int jh,double *s,double *fall,struct cell *scell){

	int k,kmin,kmax;
	double ut,vt,wt;
	
	kmin = (VER_ADVECTION_ORDER+1) / 2;
	kmax = NZ - (VER_ADVECTION_ORDER+1) / 2;

	//---------------------------------------------------
	// Horizontal part
	//---------------------------------------------------
	for(int j=jl;j<jh;j++){
		
		SCELL(j,0).top = 0.5*(SP(i,j,0)+SP(i,j,1));		// lower boundary values
		
		for(k=1;k<NZ-1;k++){

			SCELL(j,k).west  = SCELL(j,k).east;
			SCELL(j,k).east  = INTERP_EAST( SP,SIGN_B_VEL(i+1,j,k).u,i);
			SCELL(j,k).north = INTERP_NORTH(SP,SIGN_B_VEL(i,j+1,k).v,j);
		}
	}

	//---------------------------------------------------
	// Upper and lower boundary points
	// Use lower order interpolations in vertical if 
	// the interpolating stencil extends beyond model domain
	//---------------------------------------------------
	#if VER_ADVECTION_ORDER > 2
		
	for(int j=jl;j<jh;j++){
		
		SCELL(j,1).top   = INTERP_2ND_TOP(SP,1);
		
		SCELL(j,NZ-2).top   = INTERP_2ND_TOP(SP,NZ-2);
	}
	
	#endif
	
	#if VER_ADVECTION_ORDER > 4
	
	for(int j=jl;j<jh;j++){
		
		k = 2;
		//-------------------------------------------
		// sign of advecting velocity
		wt = signof(WBAR(i,j,k+1)+W(i,j,k+1)-0.5*(fall[INDEX(i,j,k)]+fall[INDEX(i,j,k+1)]) );
		//-------------------------------------------
		// Perturbation interpolation
		SCELL(j,k).top   = INTERP_3RD_TOP(SP,wt,k);
		
		k = NZ-3;
		//-------------------------------------------
		// sign of advecting velocity
		wt = signof(WBAR(i,j,k+1)+W(i,j,k+1)-0.5*(fall[INDEX(i,j,k)]+fall[INDEX(i,j,k+1)]) );
		//-------------------------------------------
		// Perturbation interpolation
		SCELL(j,k).top   = INTERP_3RD_TOP(SP,wt,k);
	}
	
	#endif
	
	//---------------------------------------------------
	// Interior points
	// Use requested interpolation order if the stencil
	// lies completely within the model domain
	//---------------------------------------------------
	for(int j=jl;j<jh;j++){
	for(k=kmin;k<kmax;k++){

		//-------------------------------------------
		// sign of advecting velocity
		wt = signof(WBAR(i,j,k+1)+W(i,j,k+1)-0.5*(fall[INDEX(i,j,k)]+fall[INDEX(i,j,k+1)]) );

		//-------------------------------------------
		// Perturbation interpolation
		SCELL(j,k).top   = INTERP_TOP(SP,wt,k);
	}}	

}

/*********************************************************************
* Interpolate a scalar field to the faces of each control volume for
* a YZ cross section. Without basic state
*
* @param i - the x-coordinate
* @param jl,jh - the high and low index bounds for the y-coordinate
* @param s - intput perturbation scalar field
* @param scell - output interpolated YZ field for perturbations
**********************************************************************/
void interpolate_scalar(int i,int jl,int jh,double *u,double *v,double *w,double *s,struct cell *scell){

	int k;
	double ut,vt,wt;

	for(int j=jl;j<jh;j++){
		
		k = 1;

		//-------------------------------------------
		// sign of advecting velocity
		ut = signof(u[INDEX(i,j,k)]);
		vt = signof(v[INDEX(i,j,k)]);

		//-------------------------------------------
		// Perturbation interpolation
		SCELL(j,k).west  = SCELL(j,k).east;
		SCELL(j,k).east  = INTERP_EAST( SP,ut,i);
		SCELL(j,k).north = INTERP_NORTH(SP,vt,j);
		SCELL(j,k).top   = k_interp2nd(SP,k+1,k);
		
		for(k=2;k<NZ-2;k++){

			//-------------------------------------------
			// sign of advecting velocity
			ut = signof(u[INDEX(i,j,k)]);
			vt = signof(v[INDEX(i,j,k)]);
			wt = signof(w[INDEX(i,j,k)]);

			//-------------------------------------------
			// Perturbation interpolation
			SCELL(j,k).west  = SCELL(j,k).east;
			SCELL(j,k).east  = INTERP_EAST( SP,ut,i);
			SCELL(j,k).north = INTERP_NORTH(SP,vt,j);
			SCELL(j,k).top   = INTERP_TOP(  SP,wt,k);
		}
		
		/***************************************************
		* For the uppermost physical point (NZ-1 is a "ghost"
		* point), we can't calculate the required 3rd order
		* interpolation, so can either set it to zero or
		* use a second order interpolation.
		****************************************************/
		k = NZ-2;

		//-------------------------------------------
		// sign of advecting velocity
		ut = signof(u[INDEX(i,j,k)]);
		vt = signof(v[INDEX(i,j,k)]);

		//-------------------------------------------
		// Perturbation interpolation
		SCELL(j,k).west  = SCELL(j,k).east;
		SCELL(j,k).east  = INTERP_EAST( SP,ut,i);
		SCELL(j,k).north = INTERP_NORTH(SP,vt,j);
		SCELL(j,k).top   = k_interp2nd(SP,k+1,k);
	}
}

/*********************************************************************
* Calculate advecting velocity field for faces of control volumes for
* a YZ cross section of the velocity field.
*
* @param i - the x-coordinate
* @param jl,jh - the high and low index bounds for the y-coordinate
**********************************************************************/
inline void advecting_velocity(int i,int j,int k){

	/***************************************************
	* Calculate advecting velocity for all perturbation
	* values on each cell face.
	****************************************************/
	// zonal wind
	UCELL(j,k).ua = 0.5*(U(i+1,j,k)+U(i,j,k));
	VCELL(j,k).ua = 0.5*(U(i+1,j,k)+U(i+1,j-1,k));
#if !HYDROSTATIC
	WCELL(j,k).ua = 0.5*(U(i+1,j,k-1)+U(i+1,j,k));
#endif

	VCELL(j,k).other_vel = 0.25 * (U(i,j,k)+U(i+1,j,k)+U(i,j-1,k)+U(i+1,j-1,k));

	// meridional wind
	UCELL(j,k).va = 0.5*(V(i,j+1,k)+V(i-1,j+1,k));
	VCELL(j,k).va = 0.5*(V(i,j+1,k)+V(i,j,k));
#if !HYDROSTATIC
	WCELL(j,k).va = 0.5*(V(i,j+1,k-1)+V(i,j+1,k));
#endif

	UCELL(j,k).other_vel = 0.25 * (V(i,j,k)+V(i,j+1,k)+V(i-1,j,k)+V(i-1,j+1,k));

	// vertical wind
	UCELL(j,k).wa = 0.5*(W(i,j,k+1)+W(i-1,j,k+1));
	VCELL(j,k).wa = 0.5*(W(i,j,k+1)+W(i,j-1,k+1));
#if !HYDROSTATIC
	WCELL(j,k).wa = 0.5*(W(i,j,k+1)+W(i,j,k));
#endif
	
	/***************************************************
	* Calculate advecting velocity for all background flow
	* values on each cell face.
	****************************************************/
	// zonal wind
	UBCELL(j,k).ua = 0.5*(UBAR(i+1,j,k)+UBAR(i,j,k));
	VBCELL(j,k).ua = 0.5*(UBAR(i+1,j,k)+UBAR(i+1,j-1,k));
#if !HYDROSTATIC
	WBCELL(j,k).ua = 0.5*(UBAR(i+1,j,k-1)+UBAR(i+1,j,k));
#endif

	// meridional wind
	UBCELL(j,k).va = 0.5*(VBAR(i,j+1,k)+VBAR(i-1,j+1,k));
	VBCELL(j,k).va = 0.5*(VBAR(i,j+1,k)+VBAR(i,j,k));
#if !HYDROSTATIC
	WBCELL(j,k).va = 0.5*(VBAR(i,j+1,k-1)+VBAR(i,j+1,k));
#endif

	// vertical wind
	UBCELL(j,k).wa = 0.5*(WBAR(i,j,k+1)+WBAR(i-1,j,k+1));
	VBCELL(j,k).wa = 0.5*(WBAR(i,j,k+1)+WBAR(i,j-1,k+1));
#if !HYDROSTATIC
	WBCELL(j,k).wa = 0.5*(WBAR(i,j,k+1)+WBAR(i,j,k));
#endif
}

/*********************************************************************
* Vertical interpolation using 2nd order advection
*
**********************************************************************/
void interpolate_velocity_vertical_2nd(int i,int j,int k){	
	//---------------------------------------------------------
	// Zonal velocity control volumes
	//---------------------------------------------------------
	UCELL(j,k).top = k_interp2nd(U,k+1,k);
	UBCELL(j,k).top = k_interp2nd(UBAR,k+1,k);
	//---------------------------------------------------------
	// Meridional velocity control volumes
	//---------------------------------------------------------
	VCELL(j,k).top = k_interp2nd(V,k+1,k);
	VBCELL(j,k).top = k_interp2nd(VBAR,k+1,k);
	//---------------------------------------------------------
	// Vertical velocity control volumes
	//---------------------------------------------------------
#if !HYDROSTATIC
	WCELL(j,k).top = k_interp2nd(W,k+1,k);
	WBCELL(j,k).top = k_interp2nd(WBAR,k+1,k);	
#endif
}

/*********************************************************************
*  Vertical interpolation using 3rd order advection
*
**********************************************************************/
void interpolate_velocity_vertical_3rd(int i,int j,int k){
	// the sign of the advecting velocity 
	// is needed for odd order interpolation
	double wsign;
	double wbsign;
	//---------------------------------------------------------
	// Zonal velocity control volumes
	//---------------------------------------------------------
#if HOR_ADVECTION_ORDER % 2 != 0
	wsign = signof_d(UCELL(j,k).wa);
#if !ISLINEAR
	wbsign = signof_d(UBCELL(j,k).wa+UCELL(j,k).wa);
#else
	wbsign = signof_d(UBCELL(j,k).wa);
#endif
#endif
	UCELL(j,k).top  = INTERP_3RD_TOP(U,wbsign,k);
	UBCELL(j,k).top = INTERP_3RD_TOP(UBAR,wsign,k);
	//---------------------------------------------------------
	// Meridional velocity control volumes
	//---------------------------------------------------------
#if HOR_ADVECTION_ORDER % 2 != 0
	wsign = signof_d(VCELL(j,k).wa);
#if !ISLINEAR
	wbsign = signof_d(VBCELL(j,k).wa+VCELL(j,k).wa);
#else
	wbsign = signof_d(VBCELL(j,k).wa);
#endif
#endif
	VCELL(j,k).top  = INTERP_3RD_TOP(V,wbsign,k);
	VBCELL(j,k).top = INTERP_3RD_TOP(VBAR,wsign,k);
	//---------------------------------------------------------
	// Vertical velocity control volumes
	//---------------------------------------------------------
#if !HYDROSTATIC
#if HOR_ADVECTION_ORDER % 2 != 0
	wsign = signof_d(WCELL(j,k).wa);
	#if !ISLINEAR
	wbsign = signof_d(WBCELL(j,k).wa+WCELL(j,k).wa);
	#else
	wbsign = signof_d(WBCELL(j,k).wa);
	#endif
#endif
	WCELL(j,k).top  = INTERP_3RD_TOP(W,wbsign,k);
	WBCELL(j,k).top = INTERP_3RD_TOP(WBAR,wsign,k);
#endif
}

/*********************************************************************
* Vertical interpolation of velocity using requested interpolation order
*
**********************************************************************/
void interpolate_velocity_vertical(int i,int j,int k){
	// the sign of the advecting velocity 
	// is needed for odd order interpolation
	double wsign = 1;
	double wbsign = 1;
	//---------------------------------------------------------
	// Zonal velocity control volumes
	//---------------------------------------------------------
#if HOR_ADVECTION_ORDER % 2 != 0
	wsign = signof_d(UCELL(j,k).wa);
#if !ISLINEAR
	wbsign = signof_d(UBCELL(j,k).wa+UCELL(j,k).wa);
#else
	wbsign = signof_d(UBCELL(j,k).wa);
#endif
#endif
	UCELL(j,k).top  = INTERP_TOP(U,wbsign,k);
	UBCELL(j,k).top = INTERP_TOP(UBAR,wsign,k);
	//---------------------------------------------------------
	// Meridional velocity control volumes
	//---------------------------------------------------------
#if HOR_ADVECTION_ORDER % 2 != 0
	wsign = signof_d(VCELL(j,k).wa);
#if !ISLINEAR
	wbsign = signof_d(VBCELL(j,k).wa+VCELL(j,k).wa);
#else
	wbsign = signof_d(VBCELL(j,k).wa);
#endif
#endif
	VCELL(j,k).top  = INTERP_TOP(V,wbsign,k);
	VBCELL(j,k).top = INTERP_TOP(VBAR,wsign,k);
	//---------------------------------------------------------
	// Vertical velocity control volumes
	//---------------------------------------------------------
#if !HYDROSTATIC
	#if HOR_ADVECTION_ORDER % 2 != 0
		wsign = signof_d(WCELL(j,k).wa);
		#if !ISLINEAR
			wbsign = signof_d(WBCELL(j,k).wa+WCELL(j,k).wa);
		#else
			wbsign = signof_d(WBCELL(j,k).wa);
		#endif
	#endif
	WCELL(j,k).top  = INTERP_TOP(W,wbsign,k);
	WBCELL(j,k).top = INTERP_TOP(WBAR,wsign,k);
#endif
}

/*********************************************************************
* 
* Interpolate the advected velocity onto each cell face
*
* The flux on the eastern side from the previous
* "i" interation becomes the flux on the western side.
* Note that this is a flux, not a velocity.
*
**********************************************************************/
void interpolate_velocity_horizontal(int i,int j,int k){
	// the sign of the advecting velocity 
	// is needed for odd order interpolation
	double usign,vsign,wsign;
	double ubsign,vbsign,wbsign;
	//---------------------------------------------------------
	// Zonal velocity points
	//---------------------------------------------------------
#if HOR_ADVECTION_ORDER % 2 != 0
	usign = signof_d(UCELL(j,k).ua);
	vsign = signof_d(UCELL(j,k).va);
#if !ISLINEAR
	ubsign = signof_d(UBCELL(j,k).ua+UCELL(j,k).ua);
	vbsign = signof_d(UBCELL(j,k).va+UCELL(j,k).va);
#else
	ubsign = signof_d(UBCELL(j,k).ua);
	vbsign = signof_d(UBCELL(j,k).va);
#endif
#endif	
	UCELL(j,k).west = UCELL(j,k).east;
	UCELL(j,k).east  = INTERP_EAST( U,ubsign,i);
	UCELL(j,k).north = INTERP_NORTH(U,vbsign,j);
	
	UBCELL(j,k).east  = INTERP_EAST( UBAR,usign,i);
	UBCELL(j,k).north = INTERP_NORTH(UBAR,vsign,j);

	//---------------------------------------------------------
	// Meridional velocity points
	//---------------------------------------------------------
#if HOR_ADVECTION_ORDER % 2 != 0
	usign = signof_d(VCELL(j,k).ua);
	vsign = signof_d(VCELL(j,k).va);	
#if !ISLINEAR
	ubsign = signof_d(VBCELL(j,k).ua+VCELL(j,k).ua);
	vbsign = signof_d(VBCELL(j,k).va+VCELL(j,k).va);
#else
	ubsign = signof_d(VBCELL(j,k).ua);
	vbsign = signof_d(VBCELL(j,k).va);
#endif
#endif
	VCELL(j,k).west = VCELL(j,k).east;
	VCELL(j,k).east  = INTERP_EAST( V,ubsign,i);
	VCELL(j,k).north = INTERP_NORTH(V,vbsign,j);
	
	VBCELL(j,k).east  = INTERP_EAST( VBAR,usign,i);
	VBCELL(j,k).north = INTERP_NORTH(VBAR,vsign,j);

	//---------------------------------------------------------
	// Vertical velocity points
	//---------------------------------------------------------
#if !HYDROSTATIC
	#if HOR_ADVECTION_ORDER % 2 != 0	
		usign = signof_d(WCELL(j,k).ua);
		vsign = signof_d(WCELL(j,k).va);
		#if !ISLINEAR
			ubsign = signof_d(WBCELL(j,k).ua+WCELL(j,k).ua);
			vbsign = signof_d(WBCELL(j,k).va+WCELL(j,k).va);
		#else
			ubsign = signof_d(WBCELL(j,k).ua);
			vbsign = signof_d(WBCELL(j,k).va);
		#endif
	#endif	
	WCELL(j,k).west = WCELL(j,k).east;
	WCELL(j,k).east  = INTERP_EAST( W,ubsign,i);
	WCELL(j,k).north = INTERP_NORTH(W,vbsign,j);
	
	WBCELL(j,k).east  = INTERP_EAST( WBAR,usign,i);
	WBCELL(j,k).north = INTERP_NORTH(WBAR,vsign,j);
#endif

}

/*********************************************************************
* Interpolate velocity field to the faces of each control volume for
* a YZ cross section. Also sets the fluxes on the eastern side of the 
* previous (i-1) YZ cross section equal to the fluxes on the western 
* side of the new YZ cross section.
*
* @param i - the x-coordinate
* @param jl,jh - the high and low index bounds for the y-coordinate
**********************************************************************/
void interpolate_velocity(int i,int jl,int jh){

	int k,kmin,kmax;

	kmin = (VER_ADVECTION_ORDER+1) / 2;
	kmax = NZ - (VER_ADVECTION_ORDER+1) / 2;

	double usign,vsign,wsign;
	double ubsign,vbsign,wbsign;

	for(int j=jl;j<jh;j++){

		//---------------------------------------------------
		// Lower boundary points
		// Use lower order interpolations if the stencil
		// extends below the grounnd
		//---------------------------------------------------
	
		#if VER_ADVECTION_ORDER > 2
	
		k = 1;

			advecting_velocity(i,j,k);
			interpolate_velocity_horizontal(i,j,k);
			interpolate_velocity_vertical_2nd(i,j,k);

		#endif
	
		#if VER_ADVECTION_ORDER > 4
	
		k = 2;

			advecting_velocity(i,j,k);
			interpolate_velocity_horizontal(i,j,k);
			interpolate_velocity_vertical_3rd(i,j,k);

		#endif

		for(k=kmin;k<kmax;k++){

			// advecting velocity is calculated using 2nd order interpolation
			/***************************************************
			* Calculate advecting velocity for all perturbation
			* values on each cell face.
			****************************************************/
			// zonal wind
			UCELL(j,k).ua = 0.5*(U(i+1,j,k)+U(i,j,k));
			VCELL(j,k).ua = 0.5*(U(i+1,j,k)+U(i+1,j-1,k));
		#if !HYDROSTATIC
			WCELL(j,k).ua = 0.5*(U(i+1,j,k-1)+U(i+1,j,k));
		#endif

			VCELL(j,k).other_vel = 0.25 * (U(i,j,k)+U(i+1,j,k)+U(i,j-1,k)+U(i+1,j-1,k));

			// meridional wind
			UCELL(j,k).va = 0.5*(V(i,j+1,k)+V(i-1,j+1,k));
			VCELL(j,k).va = 0.5*(V(i,j+1,k)+V(i,j,k));
		#if !HYDROSTATIC
			WCELL(j,k).va = 0.5*(V(i,j+1,k-1)+V(i,j+1,k));
		#endif

			UCELL(j,k).other_vel = 0.25 * (V(i,j,k)+V(i,j+1,k)+V(i-1,j,k)+V(i-1,j+1,k));

			// vertical wind
			UCELL(j,k).wa = 0.5*(W(i,j,k+1)+W(i-1,j,k+1));
			VCELL(j,k).wa = 0.5*(W(i,j,k+1)+W(i,j-1,k+1));
		#if !HYDROSTATIC
			WCELL(j,k).wa = 0.5*(W(i,j,k+1)+W(i,j,k));
		#endif
	
			/***************************************************
			* Calculate advecting velocity for all background flow
			* values on each cell face.
			****************************************************/
			// zonal wind
			UBCELL(j,k).ua = 0.5*(UBAR(i+1,j,k)+UBAR(i,j,k));
			VBCELL(j,k).ua = 0.5*(UBAR(i+1,j,k)+UBAR(i+1,j-1,k));
		#if !HYDROSTATIC
			WBCELL(j,k).ua = 0.5*(UBAR(i+1,j,k-1)+UBAR(i+1,j,k));
		#endif

			// meridional wind
			UBCELL(j,k).va = 0.5*(VBAR(i,j+1,k)+VBAR(i-1,j+1,k));
			VBCELL(j,k).va = 0.5*(VBAR(i,j+1,k)+VBAR(i,j,k));
		#if !HYDROSTATIC
			WBCELL(j,k).va = 0.5*(VBAR(i,j+1,k-1)+VBAR(i,j+1,k));
		#endif

			// vertical wind
			UBCELL(j,k).wa = 0.5*(WBAR(i,j,k+1)+WBAR(i-1,j,k+1));
			VBCELL(j,k).wa = 0.5*(WBAR(i,j,k+1)+WBAR(i,j-1,k+1));
		#if !HYDROSTATIC
			WBCELL(j,k).wa = 0.5*(WBAR(i,j,k+1)+WBAR(i,j,k));
		#endif

			/**********************************************************
			* Interpolate the advected velocity onto each cell face
			*
			* The flux on the eastern side from the previous
			* "i" interation becomes the flux on the western side.
			* Note that this is a flux, not a velocity.
			***********************************************************/

			//---------------------------------------------------------
			// Zonal velocity control volumes
			//---------------------------------------------------------
			// the sign of the advecting velocity 
			// is needed for 5th order interpolation
			usign = signof_d(UCELL(j,k).ua);
			vsign = signof_d(UCELL(j,k).va);
			wsign = signof_d(UCELL(j,k).wa);
		#if !ISLINEAR
			ubsign = signof_d(UBCELL(j,k).ua+UCELL(j,k).ua);
			vbsign = signof_d(UBCELL(j,k).va+UCELL(j,k).va);
			wbsign = signof_d(UBCELL(j,k).wa+UCELL(j,k).wa);
		#else
			ubsign = signof_d(UBCELL(j,k).ua);
			vbsign = signof_d(UBCELL(j,k).va);
			wbsign = signof_d(UBCELL(j,k).wa);
		#endif
			UCELL(j,k).west = UCELL(j,k).east;
			UCELL(j,k).east  = INTERP_5TH_EAST( U,ubsign,i);
			UCELL(j,k).north = INTERP_5TH_NORTH(U,vbsign,j);
			UCELL(j,k).top   = INTERP_3RD_TOP(  U,wbsign,k);
			
			UBCELL(j,k).east  = INTERP_5TH_EAST( UBAR,usign,i);
			UBCELL(j,k).north = INTERP_5TH_NORTH(UBAR,vsign,j);
			UBCELL(j,k).top   = INTERP_3RD_TOP(  UBAR,wsign,k);
		
			//---------------------------------------------------------
			// Meridional velocity control volumes
			//---------------------------------------------------------
			usign = signof_d(VCELL(j,k).ua);
			vsign = signof_d(VCELL(j,k).va);
			wsign = signof_d(VCELL(j,k).wa);
		#if !ISLINEAR
			ubsign = signof_d(VBCELL(j,k).ua+VCELL(j,k).ua);
			vbsign = signof_d(VBCELL(j,k).va+VCELL(j,k).va);
			wbsign = signof_d(VBCELL(j,k).wa+VCELL(j,k).wa);
		#else
			ubsign = signof_d(VBCELL(j,k).ua);
			vbsign = signof_d(VBCELL(j,k).va);
			wbsign = signof_d(VBCELL(j,k).wa);
		#endif
			VCELL(j,k).west  = VCELL(j,k).east;
			VCELL(j,k).east  = INTERP_5TH_EAST( V,ubsign,i);
			VCELL(j,k).north = INTERP_5TH_NORTH(V,vbsign,j);
			VCELL(j,k).top   = INTERP_3RD_TOP(  V,wbsign,k);
			
			VBCELL(j,k).east  = INTERP_5TH_EAST( VBAR,usign,i);
			VBCELL(j,k).north = INTERP_5TH_NORTH(VBAR,vsign,j);
			VBCELL(j,k).top   = INTERP_3RD_TOP(  VBAR,wsign,k);

			//---------------------------------------------------------
			// Vertical velocity control volumes
			//---------------------------------------------------------
	#if !HYDROSTATIC
			usign = signof_d(WCELL(j,k).ua);
			vsign = signof_d(WCELL(j,k).va);
			wsign = signof_d(WCELL(j,k).wa);
		#if !ISLINEAR
			ubsign = signof_d(WBCELL(j,k).ua+WCELL(j,k).ua);
			vbsign = signof_d(WBCELL(j,k).va+WCELL(j,k).va);
			wbsign = signof_d(WBCELL(j,k).wa+WCELL(j,k).wa);
		#else
			ubsign = signof_d(WBCELL(j,k).ua);
			vbsign = signof_d(WBCELL(j,k).va);
			wbsign = signof_d(WBCELL(j,k).wa);
		#endif
			WCELL(j,k).west  = WCELL(j,k).east;
			WCELL(j,k).east  = INTERP_5TH_EAST( W,ubsign,i);
			WCELL(j,k).north = INTERP_5TH_NORTH(W,vbsign,j);
			WCELL(j,k).top   = INTERP_3RD_TOP(  W,wbsign,k);
			
			WBCELL(j,k).east  = INTERP_5TH_EAST( WBAR,usign,i);
			WBCELL(j,k).north = INTERP_5TH_NORTH(WBAR,vsign,j);
			WBCELL(j,k).top   = INTERP_3RD_TOP(  WBAR,wsign,k);
	#endif
		}

		//---------------------------------------------------
		// Upper boundary points
		// Use lower order interpolations if the stencil
		// extends above the model domain
		//---------------------------------------------------

		#if VER_ADVECTION_ORDER > 4
	
		k = NZ-3;

			advecting_velocity(i,j,k);
			interpolate_velocity_horizontal(i,j,k);
			interpolate_velocity_vertical_3rd(i,j,k);

		#endif

		#if VER_ADVECTION_ORDER > 2

		k = NZ-2;

			advecting_velocity(i,j,k);
			interpolate_velocity_horizontal(i,j,k);
			interpolate_velocity_vertical_2nd(i,j,k);
	
		#endif
		
	}

}
