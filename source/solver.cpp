#include "stdafx.h"
#include "advection.h"
#include "turbulence.h"
#include "Heating.h"
#include "energy.h"
#include "damping.h"
#include "boundaries.h"
#include "initializer.h"
#include "pressure.h"
#include "budgets.h"
#include "pcomm.h"

/*******************************************************************************************
* 
* TOP LEVEL SUBROUTINES FOR SOLVING MODEL EQUATIONS
*
********************************************************************************************/

size_t file_time_counter = 0;

Heating heat;

void diff_terrain_vars(int size);

/*********************************************************************
* Integrate hydrostatic model equations. Advects all wind components,
* calculates anelastics pressure field, and then applies pressure
* gradient force to momentum field.
*
* step - fractional size of time step
* itr - which iteration within Runge-Kutta loop (used to deterime 
*       first guess for SOR algorithm
* il,ih,jl,jh - starting and ending grid indices
**********************************************************************/
void integrate_hydro(double step,int itr,int il,int ih,int jl, int jh){

	int buffer,size;

	if(PARALLEL){
		buffer = 3;
		size = fNX*fNY*fNZ;
	} else {
		buffer = 1;
		size = NX*NY*NZ;
	}
	/*****************************************************************
	* Calculate new momentum from everything but pressure
	******************************************************************/
	advect_uv_velocity(step,il+3,ih-3,jl+3,jh-3);

	if(USE_TURBULENT_STRESS){ apply_velocity_diffusion(step,il+3,ih-3,jl+3,jh-3);}

	/*****************************************************************
	* Zero out values at or below the terrain for velocity field
	******************************************************************/
	if(USE_TERRAIN){
		set_terrain(ups,&UISTOPO(0,0,0),size);
		set_terrain(vps,&VISTOPO(0,0,0),size);
	}

	/*****************************************************************
	* Solve anelastic pressure equation
	******************************************************************/
	solve_pressure(step);

	/*****************************************************************
	* Apply pressure gradient force to momentum equations
	******************************************************************/
	for(int i=il+buffer;i<ih-buffer;i++){
	for(int j=jl+buffer;j<jh-buffer;j++){
	for(int k=1;k<NZ-1;k++){
		
		UP(i,j,k) -= step*(dtx*cp*tbv[k]*(PI(i,j,k)-PI(i-1,j,k)));
		VP(i,j,k) -= step*(dty*cp*tbv[k]*(PI(i,j,k)-PI(i,j-1,k)));
		
	}}}
}

/*********************************************************************
* Integrate non-hydrostatic model equations. Advects all wind components,
* calculates anelastic pressure field, and then applies pressure
* gradient force to velocity field.
*
* step - fractional size of time step
* il,ih,jl,jh - starting and ending grid indices
**********************************************************************/
void integrate_non_hydro(double step,int il,int ih,int jl, int jh){

	int buffer,size;

	if(PARALLEL){
		buffer = 3;
		size = fNX*fNY*fNZ;
	} else {
		buffer = 1;
		size = NX*NY*NZ;
	}
	//-----------------------------------------------------------------
	// Calculate new velocity from everything but pressure
	//-----------------------------------------------------------------
	advect_uvw_velocity(step,il+3,ih-3,jl+3,jh-3);
	
	if(USE_TURBULENT_STRESS){ apply_velocity_diffusion(step,il+3,ih-3,jl+3,jh-3);}
	if(EXTRA_DIFFUSION && !STRETCHED_GRID){ diffusion_2nd_all(step,3,NX-3,3,NY-3);}//diffusion_6th_all(step,3,NX-3,3,NY-3);}
	
	//-----------------------------------------------------------------
	// Zero out values at or below the terrain for velocity field
	//-----------------------------------------------------------------
	if(USE_TERRAIN){
		set_terrain(ups,&UISTOPO(0,0,0),size);
		set_terrain(vps,&VISTOPO(0,0,0),size);
		set_terrain(wps, &ISTOPO(0,0,0),size);	// need a WISTOPO array!
	}

	//-----------------------------------------------------------------
	// Solve anelastic pressure equation
	//-----------------------------------------------------------------
	solve_pressure(step);
	
	//-----------------------------------------------------------------
	// Apply pressure gradient force to velocity field
	//-----------------------------------------------------------------
	for(int i=il+buffer;i<ih-buffer;i++){
	for(int j=jl+buffer;j<jh-buffer;j++){

		UP(i,j,1) -= step*(dtx*(PI(i,j,1)-PI(i-1,j,1)));
		VP(i,j,1) -= step*(dty*(PI(i,j,1)-PI(i,j-1,1)));

		for(int k=2;k<NZ-1;k++){
		
			UP(i,j,k) -= step*(    dtx*(PI(i,j,k)-PI(i-1,j,k)));
			VP(i,j,k) -= step*(    dty*(PI(i,j,k)-PI(i,j-1,k)));
			WP(i,j,k) -= step*(DTZW(k)*(PI(i,j,k)-PI(i,j,k-1)));
		}
	}}
	
	//if(PV_TRACER && step > 0.99){ diff_terrain_vars(size);}
}

/*********************************************************************
* Calcuate vertical velocity using anelastic continuity equation
* from Lipps and Hemler (1982)
**********************************************************************/
void w_velocity_LH(int il,int ih,int jl,int jh){
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		W(i,j,HTOPO(i,j)) = 0;
		W(i,j,HTOPO(i,j)+1) = 0;
		W(i,j,NZ-1) = 0;

		for(int k=2;k<NZ-1;k++){
			
			W(i,j,k) = (rhow[k-1]/rhow[k])*W(i,j,k-1)
				+ (rhou[k-1]/rhow[k])*
					( 
						-(UP(i+1,j,k-1) - UP(i,j,k-1) ) * one_d_dx	
		      			-(VP(i,j+1,k-1) - VP(i,j,k-1) ) * one_d_dy
					)*dz/mu[k-1];
		}
	}}
}

/*********************************************************************
* Calcuate vertical velocity using anelastic continuity equation
* from Durran (1989)
**********************************************************************/
void w_velocity(){

	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){

		W(i,j,HTOPOFULL(i,j)) = 0;	// boundary condition - zero vertical 
		W(i,j,HTOPOFULL(i,j)+1) = 0;	// velocity at top and bottom
		W(i,j,NZ-1) = 0;

		for(int k=HTOPOFULL(i,j)+2;k<NZ-1;k++){
	
			W(i,j,k) = ( (tbw[k-1]*rhow[k-1]) / (tbw[k-1]*rhow[k])  ) * W(i,j,k-1) 
				+ ( (tb[k-1]*rhou[k-1]) / (tbw[k]*rhow[k] ))*
					(
						-(UP(i+1,j,k-1) - UP(i,j,k-1)) * one_d_dx  	
		      			-(VP(i,j+1,k-1) - VP(i,j,k-1)) * one_d_dy
					)*dz;
		}
	}}

}

/*****************************************************************
* Zero out values at or below the terrain
******************************************************************/
void set_terrain(double * var,double *topo,int size){

	for(int i=0;i<size;i++){ var[i] *= topo[i];}
}

#if PARALLEL
/*****************************************************************
* 
******************************************************************/
void diff_terrain(double * varIn,double * varOut,double *topo,int size){

	for(int i=0;i<size;i++){ 
		
		if(topo[i]<0.1){ varOut[i] = -varIn[i];} 
		else 		   { varOut[i] = 0; 		 }
	}
}

/*****************************************************************
* 
******************************************************************/
void diff_terrain_vars(int size){

	diff_terrain(ups,u_bound,uistopos,size);
	diff_terrain(vps,v_bound,vistopos,size);
	diff_terrain(wps,w_bound,istopos,size);
	diff_terrain(thps,w_bound,istopos,size);

}
#endif

/*********************************************************************
* Substeps within RK3 loop
*
**********************************************************************/
void switch_array(double **var1,double **var2){
	
	double *temp = *var1;
	
	*var1 = *var2;
	*var2 = temp;
}

/*********************************************************************
* Substeps within RK3 loop
*
**********************************************************************/
void advance_inner(size_t num_bytes){
	
	switch_array(&us,&ups);
	switch_array(&vs,&vps);
	switch_array(&ths,&thps);

	if(!HYDROSTATIC){ switch_array(&ws,&wps);}	
}

/*********************************************************************
* Final step for RK3
*
**********************************************************************/
void advance_outer(size_t num_bytes){

	memcpy(ums,ups,num_bytes);
	memcpy(vms,vps,num_bytes);
	memcpy(thms,thps,num_bytes);

	if(!HYDROSTATIC){ memcpy(wms,wps,num_bytes);}

	advance_inner(num_bytes);
}

