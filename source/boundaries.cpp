#include "stdafx.h"
#include "boundaries.h"

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
			sponge_boundaries(wps,NX,NY);
			sponge_boundaries(thps,NX,NY);
			
			upper_lower_boundaries(ups,1,NX,1,NY);
			upper_lower_boundaries(vps,1,NX,1,NY);
			upper_lower_boundaries(thps,1,NX,1,NY);
			
			upper_lower_boundaries_zero(wps,1,NX,1,NY);
		//-------------------------------------------------------
		// PERIODIC BOUNDARIES / SERIAL VERSION
		//-------------------------------------------------------			
		} else {
			//printf("here\n");
			//upper_lower_boundaries();
			
			sponge_boundaries_north_south(ups,NX,NY);
			sponge_boundaries_north_south(vps,NX,NY);
			sponge_boundaries_north_south(wps,NX,NY);
			sponge_boundaries_north_south(thps,NX,NY);
			
			periodic_boundaries_east_west(ups);
			periodic_boundaries_east_west(vps);		
			periodic_boundaries_east_west(wps);	
			periodic_boundaries_east_west(pis);
			periodic_boundaries_east_west(thps);
			
			upper_lower_boundaries(ups,1,NX,1,NY);
			upper_lower_boundaries(vps,1,NX,1,NY);
			upper_lower_boundaries(thps,1,NX,1,NY);
			
			upper_lower_boundaries_zero(wps,1,NX,1,NY);
			
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
* Apply top and bottom boundary conditions
**********************************************************************/
void upper_lower_boundaries_zero(double *var,int il,int ih,int jl,int jh){
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){

		var[INDEX(i,j,   0)] = 0;
		var[INDEX(i,j,   1)] = 0;
		var[INDEX(i,j,NZ-1)] = 0;
	}}
}
#if 0
/*********************************************************************
*
*
**********************************************************************/
void mirror_boundaries(double s[NX][NY][NZ]){

	for(int k=1;k<NZ-1;k++){

		for(int j=0;j<NY;j++){

			s[0][j][k] = s[1][j][k];
			s[NX-1][j][k] = s[NX-2][j][k];
		}
	
		for(int i=0;i<NX;i++){

			s[i][0][k] = s[i][1][k];
			s[i][NY-1][k] = s[i][NY-2][k];
		}
	}
}
#endif
/*********************************************************************
*
*
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


#if 0
/*********************************************************************
*
*
**********************************************************************/
void mirror_boundaries2d(double s[NX][NY]){

	for(int j=0;j<NY;j++){

		s[0][j] = s[1][j];
		s[NX-1][j] = s[NX-2][j];
	}

	for(int i=0;i<NX;i++){

		s[i][0] = s[i][1];
		s[i][NY-1] = s[i][NY-2];
	}

}
#endif
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


#if 0

/*********************************************************************
*
*
**********************************************************************/
void zero_boundaries(){

	/*********************************
	* Eastern boundary of domain
	*********************************/
	for(int i=NX-1-iebuffer;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=1;k<NZ-1;k++){
	
		UM(i,j,k) = 0;	U(i,j,k) = 0;
		
		VM(i,j,k) = 0;	V(i,j,k) = 0;
	#if HYDROSTATIC
		W(i,j,k) = 0;
	#else
		WM(i,j,k) = 0;	W(i,j,k) = 0;
	#endif
		THM(i,j,k) = 0;	TH(i,j,k) = 0;
	}}}
 
	/*********************************
	* Western boundary of domain
	**********************************/
	for(int i=1;i<iwbuffer+1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=1;k<NZ-1;k++){
	
		UM(i,j,k) = 0;	U(i,j,k) = 0;
		
		VM(i,j,k) = 0;	V(i,j,k) = 0;
	#if HYDROSTATIC
		W(i,j,k) = 0;
	#else
		WM(i,j,k) = 0;	W(i,j,k) = 0;
	#endif
		THM(i,j,k) = 0;	TH(i,j,k) = 0;
	}}}

	/*********************************
	* Southern boundary of domain
	**********************************/
	for(int j=1;j<jsbuffer+1;j++){
	for(int i=1;i<NX-1;i++){
	for(int k=1;k<NZ-1;k++){
	
		UM(i,j,k) = 0;	U(i,j,k) = 0;
		
		VM(i,j,k) = 0;	V(i,j,k) = 0;
	#if HYDROSTATIC
		W(i,j,k) = 0;
	#else
		WM(i,j,k) = 0;	W(i,j,k) = 0;
	#endif
		THM(i,j,k) = 0;	TH(i,j,k) = 0;
	}}}

	/*********************************
	* Northern boundary of domain
	**********************************/
	for(int j=NY-1-jnbuffer;j<NY-1;j++){
	for(int i=1;i<NX-1;i++){
	for(int k=1;k<NZ-1;k++){

		UM(i,j,k) = 0;	U(i,j,k) = 0;
		
		VM(i,j,k) = 0;	V(i,j,k) = 0;
	#if HYDROSTATIC
		W(i,j,k) = 0;
	#else
		WM(i,j,k) = 0;	W(i,j,k) = 0;
	#endif
		THM(i,j,k) = 0;	TH(i,j,k) = 0;
	}}}
}
#endif

#if 0

/*********************************************************************
*
*
**********************************************************************/
void periodic_ew_sponge_ns_boundaries(){

	double coef;

	/*********************************
	* east-west period boundaries
	**********************************/
	for(int k=0;k<NZ;k++){
	for(int j=0;j<NY;j++){

		PERIODIC_EW(UP)
		PERIODIC_EW(VP)
		PERIODIC_EW(WP)
		PERIODIC_EW(PI)
		PERIODIC_EW(THP)
	}}
	
	/*********************************
	* Sponge southern boundary
	**********************************/
	for(int j=1;j<jsbuffer+1;j++){
 
		coef = 0.5 * (1.-cos(trigpi*(double)(jsbuffer-j+1)/(double)(jsbuffer)));
 
		for(int i=0;i<NX;i++){
		for(int k=1;k<NZ-1;k++){
		
			DAMPVARIABLE(UP);
			DAMPVARIABLE(VP);
			DAMPVARIABLE(WP);
			DAMPVARIABLE(THP);
		}}
	}

	/*********************************
	* Sponge northern boundary
	**********************************/
	for(int j=NY-1-jnbuffer;j<NY-1;j++){
 	
		coef = 0.5 * (1.-cos(trigpi*(double)(j+jnbuffer-NY+2)/(double)jnbuffer));
		//printf("%d coef = %f\n",j,coef);
		for(int i=0;i<NX;i++){
		for(int k=1;k<NZ-1;k++){
 
			DAMPVARIABLE(UP);
			DAMPVARIABLE(VP);
			DAMPVARIABLE(WP);
			DAMPVARIABLE(THP);
		}}
 	}
}
#endif

#if 0
/*********************************************************************
* Sponge boundary condition for microphysics variables
*
**********************************************************************/
void sp_bound_microphysics(){

	double coef;

	/*********************************
	* Eastern boundary of domain
	*********************************/
	for(int i=NX-1-iebuffer;i<NX-1;i++){
 	
		// calculate coefficients
 		coef = 0.5 * (1.-cos(trigpi*(double)(i+iebuffer-NX+2)/(double)iebuffer));

		for(int j=1;j<NY-1;j++){
		for(int k=1;k<NZ-1;k++){
 		
			DAMPVARIABLE(QVP);
			DAMPVARIABLE(QCP);
			DAMPVARIABLE(QRP);
		#if USE_ICE
			DAMPVARIABLE(QSP);
			DAMPVARIABLE(QIP);		
		#endif
 		}}
 	}

	/*********************************
	* Western boundary of domain
	**********************************/
	for(int i=1;i<iwbuffer+1;i++){
 
		// calculate coefficients
 		coef = .5 * (1.-cos(trigpi*(double)(iwbuffer-i+1)/(double)iwbuffer));

		for(int j=1;j<NY-1;j++){
		for(int k=1;k<NZ-1;k++){
		
			DAMPVARIABLE(QVP);
			DAMPVARIABLE(QCP);
			DAMPVARIABLE(QRP);
		#if USE_ICE
			DAMPVARIABLE(QSP);
			DAMPVARIABLE(QIP);		
		#endif
 		}}
	}

	/*********************************
	* Southern boundary of domain
	**********************************/
	for(int j=1;j<jsbuffer+1;j++){
 
		// calculate coefficients
		coef = .5 * (1.-cos(trigpi*(double)(jsbuffer-j+1)/(double)(jsbuffer)));

		for(int i=1;i<NX-1;i++){
		for(int k=1;k<NZ-1;k++){
		
			DAMPVARIABLE(QVP);
			DAMPVARIABLE(QCP);
			DAMPVARIABLE(QRP);
		#if USE_ICE
			DAMPVARIABLE(QSP);
			DAMPVARIABLE(QIP);		
		#endif
		}}
	}

	/*********************************
	* Northern boundary of domain
	**********************************/
	for(int j=NY-1-jnbuffer;j<NY-1;j++){
 	
		// calculate proper coefficients
		coef = .5 * (1.-cos(trigpi*(double)(j+jnbuffer-NY+2)/(double)jnbuffer));

		for(int i=1;i<NX-1;i++){
		for(int k=1;k<NZ-1;k++){
		
			DAMPVARIABLE(QVP);
			DAMPVARIABLE(QCP);
			DAMPVARIABLE(QRP);
		#if USE_ICE
			DAMPVARIABLE(QSP);
			DAMPVARIABLE(QIP);		
		#endif
		}}
 	}

	/*********************************
	* Upper and lower boundary of domain
	**********************************/
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		QVP(i,j,0) = QVP(i,j,1);
		QVP(i,j,NZ-1) = QVP(i,j,NZ-2);
		QCP(i,j,0) = QCP(i,j,1);
		QCP(i,j,NZ-1) = QCP(i,j,NZ-2);
		QRP(i,j,0) = QRP(i,j,1);
		QRP(i,j,NZ-1) = QRP(i,j,NZ-2);
	#if USE_ICE
		QSP(i,j,0) = QSP(i,j,1);
		QSP(i,j,NZ-1) = QSP(i,j,NZ-2);
		QIP(i,j,0) = QIP(i,j,1);
		QIP(i,j,NZ-1) = QIP(i,j,NZ-2);		
	#endif
	}}
}
#endif
#if 0

/*********************************************************************
*
*
**********************************************************************/
void apply_xz_sponge(double * var,double c,int j,int il,int ih){

	//printf("j = %d\n",j);

	for(int i=il;i<ih;i++){
	for(int k=1;k<NZ-1;k++){
		
		var[INDEX(i,j,k)] -= c*var[INDEX(i,j,k)];
	}}
}

/*********************************************************************
*
*
**********************************************************************/
void apply_yz_sponge(double * var,double c,int i,int jl,int jh){
	
	//printf("i = %d\n",i);
	
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		
		var[INDEX(i,j,k)] -= c*var[INDEX(i,j,k)];
	}}
}

/*********************************************************************
*
*
**********************************************************************/
double cos_profile(int myIndex,int width){
	
	//printf("c = %f\n",0.5 * (1.-cos(trigpi*(double)( myIndex )/(double)width)));
	
	return 0.5 * (1.-cos(trigpi*(double)( myIndex )/(double)width));
}

/*********************************************************************
*
*			//upper_lower_boundaries();			
			
			apply_sponge_boundaries(ups);
			apply_sponge_boundaries(vps);
			apply_sponge_boundaries(wps);
			apply_sponge_boundaries(thps);
			
**********************************************************************/
void apply_sponge_boundaries(double *var){
		
	//-----------------------------------------
	// Eastern boundary of domain
	//-----------------------------------------
	for(int i=NX-1-iebuffer;i<NX-1;i++){

		apply_yz_sponge(var,cos_profile(i+iebuffer-NX+2,iebuffer),i,1,NY-1);
 	}
	//-----------------------------------------
	// Western boundary of domain
	//-----------------------------------------
	for(int i=1;i<iwbuffer+1;i++){

		apply_yz_sponge(var,cos_profile(iwbuffer-i+1,iwbuffer),i,1,NY-1);
 	}
	//-----------------------------------------
	// Southern boundary of domain
	//-----------------------------------------	
	for(int j=1;j<jsbuffer+1;j++){
 
		apply_xz_sponge(var,cos_profile(jsbuffer-j+1,jsbuffer),j,1,NX-1);
	}
	//-----------------------------------------
	// Northern boundary of domain
	//-----------------------------------------	
	for(int j=NY-1-jnbuffer;j<NY-1;j++){
 	
		apply_xz_sponge(var,cos_profile(j+jnbuffer-NY+2,jnbuffer),j,1,NX-1);
	}
}

#endif

#if 0
/*********************************************************************
* Apply top and bottom boundary conditions
**********************************************************************/
void upper_lower_boundaries(){
	
	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){

		UPPERLOWERBOUND(UP)
		UPPERLOWERBOUND(VP)
		UPPERLOWERBOUND(THP)

	#if HYDROSTATIC
		W(i,j,0) = 0;
		W(i,j,1) = 0;
		W(i,j,NZ-1) = 0;
	#else
		WP(i,j,0) = 0;
		WP(i,j,1) = 0;
		WP(i,j,NZ-1) = 0;
	#endif

	}}
}
#endif

