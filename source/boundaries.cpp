#include "stdafx.h"
#include "boundaries.h"
#include "pcomm.h"

#define VERBOSE_BOUNDARIES false

#define DAMPVARIABLE(v) v(i,j,k) = v(i,j,k) - coef*v(i,j,k)

#define UPPERLOWERBOUND(v) v(i,j,0)    = v(i,j,1   );	\
						   v(i,j,NZ-1) = v(i,j,NZ-2);

#define PERIODIC_EW(v) v(0,j,k) = v(NX-6,j,k);	\
					   v(1,j,k) = v(NX-5,j,k);	\
					   v(2,j,k) = v(NX-4,j,k);	\
					   v(NX-3,j,k) = v(3,j,k);  \
					   v(NX-2,j,k) = v(4,j,k);  \
					   v(NX-1,j,k) = v(5,j,k);


double *east_coeff,*west_coeff,*north_coeff,*south_coeff;
int eb,ee,wb,we,nb,ne,sb,se;

int iwbuffer = 5;
int iebuffer = 5;
int jnbuffer = 5;
int jsbuffer = 5;

void apply_sponge_boundaries(double *var);
void upper_lower_boundaries();
void sponge_boundaries(double *,int,int);
void p_sponge_boundaries(double*,int,int,int,int,int,int,int,int);
void p_mp_sponge_boundaries(int,int,int,int,int,int,int,int);
void apply_ew_sponge(double * var,double * coeff,int il,int ih,int jl,int jh,int kl,int kh);
void apply_ns_sponge(double * var,double * coeff,int il,int ih,int jl,int jh,int kl,int kh);
void upper_lower_boundaries(double *var,int il,int ih,int jl,int jh);
void upper_lower_boundaries_zero(double *var,int il,int ih,int jl,int jh);
void periodic_EW_boundaries(double *var);
void sponge_boundaries_north_south(double *var,int ih,int jh);
void periodic_boundaries_east_west(double *var);

/*********************************************************************
* Initialize boundary stuff
**********************************************************************/
void init_boundaries(int eastLength,int westLength,int northLength,int southLength,int halo){

	//-----------------------------------------------------------
	// Arrays for damping coefficients
	//-----------------------------------------------------------	
	east_coeff  = (double*) calloc(eastLength ,sizeof(double));
	west_coeff  = (double*) calloc(westLength ,sizeof(double));
	north_coeff = (double*) calloc(northLength,sizeof(double));
	south_coeff = (double*) calloc(southLength,sizeof(double));

	//-----------------------------------------------------------
	// Beginning and ending coordinates of boundaries
	//-----------------------------------------------------------
	if(halo>0){
		eb = fNX - halo - eastLength;
		ee = fNX - halo;
		wb = halo;
		we = westLength + halo;
		nb = fNY - halo - northLength;
		ne = fNY - halo;
		sb = halo;
		se = southLength+halo;
	} else {
		eb = NX - eastLength;
		ee = NX;
		wb = 0;
		we = westLength;
		nb = NY - northLength;
		ne = NY;
		sb = 0;
		se = southLength;
	}
	
	//-----------------------------------------------------------
	// Coefficients for boundary damping
	//-----------------------------------------------------------	
	for(int i=0;i<eastLength;i++){ east_coeff[i] = 0.5*cos_profile(-0.5,0,(double)(i+1)/(double)eastLength,1);}
	
	for(int i=0;i<westLength;i++){ west_coeff[i] = 0.5*cos_profile(0,0.5,(double)i/(double)westLength,1);}
	
	for(int i=0;i<northLength;i++){ north_coeff[i] = 0.5*cos_profile(-0.5,0,(double)(i+1)/(double)northLength,1);}
	
	for(int i=0;i<southLength;i++){ south_coeff[i] = 0.5*cos_profile(0,0.5,(double)i/(double)southLength,1);}

	//-----------------------------------------------------------
	// Debugging information
	//-----------------------------------------------------------
	#if VERBOSE_BOUNDARIES
	printf("Eastern boundary indices from %d to %d\n",eb,ee);
	printf("Western boundary indices from %d to %d\n",wb,we);
	printf("Northern boundary indices from %d to %d\n",nb,ne);
	printf("Southern boundary indices from %d to %d\n",sb,se);

	printf("The eastern boundary coefficients are ");
	for(int i=0;i<eastLength;i++){ printf("%f ",east_coeff[i]);}
	printf("\n");
	
	printf("The western boundary coefficients are ");
	for(int i=0;i<westLength;i++){ printf("%f ",west_coeff[i]);}
	printf("\n");
	
	printf("The northern boundary coefficients are ");
	for(int i=0;i<northLength;i++){ printf("%f ",north_coeff[i]);}
	printf("\n");
	
	printf("The southern boundary coefficients are ");
	for(int i=0;i<southLength;i++){ printf("%f ",south_coeff[i]);}
	printf("\n");
	#endif

}

/*********************************************************************
* Apply top and bottom boundary conditions
**********************************************************************/
void apply_boundary_condition(int mpi_proc_null){

	#if PARALLEL
	//-----------------------------------------------------------
	// PARALLEL VERSION
	//-----------------------------------------------------------
		p_sponge_boundaries(ups,east,west,north,south,halo_buffer,fNX,fNY,mpi_proc_null);
		p_sponge_boundaries(vps,east,west,north,south,halo_buffer,fNX,fNY,mpi_proc_null);
		p_sponge_boundaries(wps,east,west,north,south,halo_buffer,fNX,fNY,mpi_proc_null);
		p_sponge_boundaries(thps,east,west,north,south,halo_buffer,fNX,fNY,mpi_proc_null);
		
		upper_lower_boundaries(ups,3,fNX-3,3,fNY-3);
		upper_lower_boundaries(vps,3,fNX-3,3,fNY-3);
		upper_lower_boundaries(thps,3,fNX-3,3,fNY-3);
		
		upper_lower_boundaries_zero(wps,3,fNX-3,3,fNY-3);
	#else
	//-----------------------------------------------------------
	// SERIAL VERSION
	//-----------------------------------------------------------
		//-------------------------------------------------------
		// SPONGE BOUNDARIES / SERIAL VERSION
		//-------------------------------------------------------
		if(!PERIODIC_BOUNDARIES){
			
			sponge_boundaries(ups,NX,NY);
			sponge_boundaries(vps,NX,NY);
			sponge_boundaries(thps,NX,NY);
			
			upper_lower_boundaries(ups,1,NX,1,NY);
			upper_lower_boundaries(vps,1,NX,1,NY);
			upper_lower_boundaries(thps,1,NX,1,NY);
			
			if(HYDROSTATIC){
				sponge_boundaries(ws,NX,NY);
				upper_lower_boundaries_zero(ws,1,NX,1,NY);
			} else {
				sponge_boundaries(wps,NX,NY);
				upper_lower_boundaries_zero(wps,1,NX,1,NY);				
			}
		//-------------------------------------------------------
		// PERIODIC BOUNDARIES / SERIAL VERSION
		//-------------------------------------------------------
		} else {

			sponge_boundaries_north_south(ups,NX,NY);
			sponge_boundaries_north_south(vps,NX,NY);		
			sponge_boundaries_north_south(thps,NX,NY);
			
			periodic_boundaries_east_west(ups);
			periodic_boundaries_east_west(vps);		
			periodic_boundaries_east_west(pis);
			periodic_boundaries_east_west(thps);
			
			upper_lower_boundaries(ups,1,NX,1,NY);
			upper_lower_boundaries(vps,1,NX,1,NY);
			upper_lower_boundaries(thps,1,NX,1,NY);
			
			if(HYDROSTATIC){
				sponge_boundaries_north_south(ws,NX,NY);
				periodic_boundaries_east_west(ws);
				upper_lower_boundaries_zero(ws,1,NX,1,NY);
				
			} else {
				sponge_boundaries_north_south(wps,NX,NY);
				periodic_boundaries_east_west(wps);
				upper_lower_boundaries_zero(wps,1,NX,1,NY);	
			}
			
			
			
			//periodic_ew_sponge_ns_boundaries();
		}
	#endif
}

/*********************************************************************
* Boundary conditions for microphysics variables
*
* MPI_PROC_NULL the value of a null process defined by mpi.h
**********************************************************************/
void apply_boundary_condition_microphysics(int mpi_proc_null){
	
	#if PARALLEL
	//-----------------------------------------------------------
	// PARALLEL VERSION
	//-----------------------------------------------------------		
		p_sponge_boundaries(qvps,east,west,north,south,halo_buffer,fNX,fNY,mpi_proc_null);
		p_sponge_boundaries(qrps,east,west,north,south,halo_buffer,fNX,fNY,mpi_proc_null);
		p_sponge_boundaries(qcps,east,west,north,south,halo_buffer,fNX,fNY,mpi_proc_null);

		upper_lower_boundaries(qvps,3,fNX-3,3,fNY-3);
		upper_lower_boundaries(qrps,3,fNX-3,3,fNY-3);
		upper_lower_boundaries(qcps,3,fNX-3,3,fNY-3);
				
		if(USE_ICE){
			
			p_sponge_boundaries(qips,east,west,north,south,halo_buffer,fNX,fNY,mpi_proc_null);			
			p_sponge_boundaries(qsps,east,west,north,south,halo_buffer,fNX,fNY,mpi_proc_null);

			upper_lower_boundaries(qips,3,fNX-3,3,fNY-3);
			upper_lower_boundaries(qsps,3,fNX-3,3,fNY-3);			
		}

		if(MICROPHYSICS_OPTION==3){
			p_sponge_boundaries(qgps,east,west,north,south,halo_buffer,fNX,fNY,mpi_proc_null);			
			p_sponge_boundaries(nrps,east,west,north,south,halo_buffer,fNX,fNY,mpi_proc_null);
			p_sponge_boundaries(nips,east,west,north,south,halo_buffer,fNX,fNY,mpi_proc_null);

			upper_lower_boundaries(qgps,3,fNX-3,3,fNY-3);
			upper_lower_boundaries(nrps,3,fNX-3,3,fNY-3);
			upper_lower_boundaries(nips,3,fNX-3,3,fNY-3);	

		}
	//-----------------------------------------------------------
	// SERIAL VERSION
	//-----------------------------------------------------------			
	#else
		
		if(!PERIODIC_BOUNDARIES){
			//-------------------------------------------------------
			// SPONGE BOUNDARIES / SERIAL VERSION
			//-------------------------------------------------------	
			sponge_boundaries(qvps,NX,NY);
			sponge_boundaries(qcps,NX,NY);
			sponge_boundaries(qrps,NX,NY);
		
			upper_lower_boundaries(qvps,1,NX,1,NY);
			upper_lower_boundaries(qcps,1,NX,1,NY);
			upper_lower_boundaries(qrps,1,NX,1,NY);
		
			if(USE_ICE){
				//---------------------------------------------------
				// ICE MICROPHYSICS	/ SERIAL VERSION
				//---------------------------------------------------
				sponge_boundaries(qips,NX,NY);
				sponge_boundaries(qsps,NX,NY);
			
				upper_lower_boundaries(qips,1,NX,1,NY);
				upper_lower_boundaries(qsps,1,NX,1,NY);		
			}
			
		} else {
			//-------------------------------------------------------
			// PERIODIC BOUNDARIES / SERIAL VERSION
			//-------------------------------------------------------
			sponge_boundaries_north_south(qvps,NX,NY);
			sponge_boundaries_north_south(qcps,NX,NY);
			sponge_boundaries_north_south(qrps,NX,NY);
				
			periodic_boundaries_east_west(qvps);
			periodic_boundaries_east_west(qcps);		
			periodic_boundaries_east_west(qrps);
		
			upper_lower_boundaries(qvps,1,NX,1,NY);
			upper_lower_boundaries(qcps,1,NX,1,NY);
			upper_lower_boundaries(qrps,1,NX,1,NY);	
				
			if(USE_ICE){
				//---------------------------------------------------
				// ICE MICROPHYSICS	/ SERIAL VERSION
				//---------------------------------------------------
				sponge_boundaries_north_south(qips,NX,NY);
				sponge_boundaries_north_south(qsps,NX,NY);	
		
				periodic_boundaries_east_west(qips);
				periodic_boundaries_east_west(qsps);
			
				upper_lower_boundaries(qips,1,NX,1,NY);
				upper_lower_boundaries(qsps,1,NX,1,NY);	
			}
		}
		

	#endif

}

/*********************************************************************
* Apply top and bottom boundary conditions
**********************************************************************/
void upper_lower_boundaries(double *s){
	
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		s[FULL_ARRAY_INDEX(i,j,0)] = s[FULL_ARRAY_INDEX(i,j,1)];
		s[FULL_ARRAY_INDEX(i,j,NZ-1)] = s[FULL_ARRAY_INDEX(i,j,NZ-2)];
	
	}}
}


/*********************************************************************
* Apply top and bottom boundary conditions
**********************************************************************/
void upper_lower_boundaries(double *var,int il,int ih,int jl,int jh){
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		var[INDEX(i,j,   0)] = var[INDEX(i,j,   1)];
		var[INDEX(i,j,NZ-1)] = var[INDEX(i,j,NZ-2)];
	}}
}

/*********************************************************************
* Set top and bottom boundary values to zero
**********************************************************************/
void upper_lower_boundaries_zero(double *var,int il,int ih,int jl,int jh){
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		var[INDEX(i,j,   0)] = 0;
		var[INDEX(i,j,   1)] = 0;
		var[INDEX(i,j,NZ-1)] = 0;
	}}
}

/*********************************************************************
* Mirror boundaries on the edges of the domain for a 2D array
**********************************************************************/
void mirror_boundaries_2d(double * s){

	for(int j=0;j<NY;j++){

		s[INDEX2D(   0,j)] = s[INDEX2D(   1,j)];
		s[INDEX2D(NX-1,j)] = s[INDEX2D(NX-2,j)];
	}

	for(int i=0;i<NX;i++){

		s[INDEX2D(i,   0)] = s[INDEX2D(i,1   )];
		s[INDEX2D(i,NY-1)] = s[INDEX2D(i,NY-2)];
	}

}

/*********************************************************************
* Mirror boundaries on the north and south edges of the domain
* for a 2D array
**********************************************************************/
void mirror_boundaries_ns_2d(double * s,int il,int ih){

	for(int i=il;i<ih;i++){

		s[INDEX2D(i,   0)] = s[INDEX2D(i,1   )];
		s[INDEX2D(i,NY-1)] = s[INDEX2D(i,NY-2)];
	}
}

/*********************************************************************
* Mirror boundaries on the edges of the domain for a 3D array
**********************************************************************/
void mirror_boundaries(double * s){

	for(int k=1;k<NZ-1;k++){

		for(int j=0;j<NY;j++){

			s[FULL_ARRAY_INDEX(0,j,k)] = s[FULL_ARRAY_INDEX(1,j,k)];
			s[FULL_ARRAY_INDEX(NX-1,j,k)] = s[FULL_ARRAY_INDEX(NX-2,j,k)];
		}
	
		for(int i=0;i<NX;i++){

			s[FULL_ARRAY_INDEX(i,0,k)] = s[FULL_ARRAY_INDEX(i,1,k)];
			s[FULL_ARRAY_INDEX(i,NY-1,k)] = s[FULL_ARRAY_INDEX(i,NY-2,k)];
		}
	}
}

/*********************************************************************
*
*
**********************************************************************/
void p_mirror_boundaries(double *s,int east,int west,int north,int south,int halo,int ih,int jh){

	if(west==-1){

		for(int i=0;i<halo;i++){
		for(int j=halo;j<jh-halo;j++){
		for(int k=0;k<NZ;k++){

			s[INDEX(i,j,k)] = s[INDEX(halo,j,k)];

		}}}
	}

	if(east==-1){

		for(int i=ih-halo;i<ih;i++){
		for(int j=halo;j<jh-halo;j++){
		for(int k=0;k<NZ;k++){

			s[INDEX(i,j,k)] = s[INDEX(ih-1-halo,j,k)];

		}}}
	}

	if(north==-1){

		for(int i=halo;i<ih-halo;i++){
		for(int j=jh-halo;j<jh;j++){
		for(int k=0;k<NZ;k++){

			s[INDEX(i,j,k)] = s[INDEX(i,jh-1-halo,k)];
		}}}
	}

	if(south==-1){

		for(int i=halo;i<ih-halo;i++){
		for(int j=0;j<halo;j++){
		for(int k=0;k<NZ;k++){

			s[INDEX(i,j,k)] = s[INDEX(i,halo,k)];
		}}}
	
	}

}

/*********************************************************************
*
*
**********************************************************************/
void apply_ns_sponge(double * var,double * coeff,int il,int ih,int jl,int jh,int kl,int kh){

	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=kl;k<kh;k++){

		var[INDEX(i,j,k)] -= coeff[j-jl]*var[INDEX(i,j,k)];

	}}}
}

/*********************************************************************
*
*
**********************************************************************/
void apply_ew_sponge(double * var,double * coeff,int il,int ih,int jl,int jh,int kl,int kh){

	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=kl;k<kh;k++){

		var[INDEX(i,j,k)] -= coeff[i-il]*var[INDEX(i,j,k)];

	}}}
}

/*********************************************************************
*
*
**********************************************************************/
void sponge_boundaries_north_south(double *var,int ih,int jh){
	//-----------------------------------------
	// Southern boundary of domain
	//-----------------------------------------
	apply_ns_sponge(var,south_coeff,1,ih,sb,se,1,NZ);
	//-----------------------------------------
	// Northern boundary of domain
	//-----------------------------------------
	apply_ns_sponge(var,north_coeff,1,ih,nb,ne,1,NZ);
}

/*********************************************************************
*
*
**********************************************************************/
void sponge_boundaries(double *var,int ih,int jh){
	//-----------------------------------------
	// Eastern boundary of domain
	//-----------------------------------------
	apply_ew_sponge(var,east_coeff,eb,ee,1,jh,1,NZ);
	//-----------------------------------------
	// Western boundary of domain
	//-----------------------------------------
	apply_ew_sponge(var,west_coeff,wb,we,1,jh,1,NZ);
	//-----------------------------------------
	// Southern boundary of domain
	//-----------------------------------------
	apply_ns_sponge(var,south_coeff,1,ih,sb,se,1,NZ);
	//-----------------------------------------
	// Northern boundary of domain
	//-----------------------------------------
	apply_ns_sponge(var,north_coeff,1,ih,nb,ne,1,NZ);
}

/*********************************************************************
* Apply sponge boundaries for parallel model. Uses the nullprocess ID
* and the IDs of neighboring processes to determine if the process
* is at one of the lateral boundaries
*
* var - variable to apply sponge boundaries
* east,west,north,south - ID of neighboring process
* ih,jh - maximum array index to apply boundaries
* nullprocess - the integer given to a null process
**********************************************************************/
void p_sponge_boundaries(double *var,int east,int west,int north,int south,int halo,int ih,int jh,int nullprocess){
	//-----------------------------------------
	// Eastern boundary of domain
	//-----------------------------------------
	if(east==nullprocess){ apply_ew_sponge(var,east_coeff,eb,ee,1,jh,1,NZ);}
	//-----------------------------------------
	// Western boundary of domain
	//-----------------------------------------
	if(west==nullprocess){ apply_ew_sponge(var,west_coeff,wb,we,1,jh,1,NZ);}
	//-----------------------------------------
	// Southern boundary of domain
	//-----------------------------------------
	if(south==nullprocess){ apply_ns_sponge(var,south_coeff,1,ih,sb,se,1,NZ);}
	//-----------------------------------------
	// Northern boundary of domain
	//-----------------------------------------
	if(north==nullprocess){ apply_ns_sponge(var,north_coeff,1,ih,nb,ne,1,NZ);}
}

/*********************************************************************
*
*
**********************************************************************/
void periodic_pressure_boundaries(){

	for(int k=0;k<NZ;k++){
	for(int j=0;j<NY;j++){

		PERIODIC_EW(PI)
	}}
}

/*********************************************************************
*
*
**********************************************************************/
void periodic_uvw_boundaries(){

	for(int k=0;k<NZ;k++){
	for(int j=0;j<NY;j++){

		PERIODIC_EW(UP)
		PERIODIC_EW(VP)
		PERIODIC_EW(WP)
	}}
}

/*********************************************************************
* Period east-west boundaries for serial version. Assumes three grid
* points are repeated.
*
**********************************************************************/
void periodic_boundaries_east_west(double *var){

	for(int k=0;k<NZ;k++){
	for(int j=0;j<NY;j++){

		var[FULL_ARRAY_INDEX(0,j,k)] 	= var[FULL_ARRAY_INDEX(NX-6,j,k)];
		var[FULL_ARRAY_INDEX(1,j,k)] 	= var[FULL_ARRAY_INDEX(NX-5,j,k)];
		var[FULL_ARRAY_INDEX(2,j,k)] 	= var[FULL_ARRAY_INDEX(NX-4,j,k)];
		var[FULL_ARRAY_INDEX(NX-3,j,k)] = var[FULL_ARRAY_INDEX(3,j,k)];
		var[FULL_ARRAY_INDEX(NX-2,j,k)] = var[FULL_ARRAY_INDEX(4,j,k)]; 
		var[FULL_ARRAY_INDEX(NX-1,j,k)] = var[FULL_ARRAY_INDEX(5,j,k)];
	}}
}

/*********************************************************************
* Period east-west boundaries for serial version. Assumes three grid
* points are repeated. For 2D arrays.
*
**********************************************************************/
void periodic_boundaries_east_west_2d(double *var){

	for(int j=0;j<NY;j++){

		var[INDEX2D(0,j)] 	= var[INDEX2D(NX-6,j)];
		var[INDEX2D(1,j)] 	= var[INDEX2D(NX-5,j)];
		var[INDEX2D(2,j)]   = var[INDEX2D(NX-4,j)];
		var[INDEX2D(NX-3,j)] = var[INDEX2D(3,j)];
		var[INDEX2D(NX-2,j)] = var[INDEX2D(4,j)]; 
		var[INDEX2D(NX-1,j)] = var[INDEX2D(5,j)];
	}
}
