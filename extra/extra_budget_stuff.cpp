
#if 0
//-------------------------------------------------------------
// Divergence
//-------------------------------------------------------------

/*
fdivg[ind] += 
	
	(
	-one_d_dx*one_d_dx*(PI(i+1,j,k) + PI(i-1,j,k) + PI(i,j+1,k) + PI(i,j-1,k)- 4.0*PI(i,j,k)) 
		
	+ FC(j) * vort_at_scalar
		
	+ 2 * one_d_dx * one_d_dx * (
		(U(i+1,j,k)-U(i,j,k)) * (V(i,j+1,k)-V(i,j,k)) - 
		(U_AT_V(U,i,j+1,k) - U_AT_V(U,i,j,k)) * ( V_AT_U(V,i+1,j,k) - V_AT_U(V,i,j,k))
		
		- (U_AT_V(UBAR,i,j+1,k) - U_AT_V(UBAR,i,j,k)) * ( V_AT_U(V,i+1,j,k) - V_AT_U(V,i,j,k))
	)
												
		- 0.5 * (U(i+1,j,k)+U(i,j,k)) * DFDY(j)
											
	) * dt
						
	- ((U(i+1,j,k)*(DIVG(i+1,j,k)+DIVG(i,j,k)) - U(i,j,k)*(DIVG(i,j,k)+DIVG(i-1,j,k)))
	+  (V(i,j+1,k)*(DIVG(i,j+1,k)+DIVG(i,j,k)) - V(i,j,k)*(DIVG(i,j,k)+DIVG(i,j-1,k))) )*0.5*dtx

	- ((UBAR(i+1,j,k)*(DIVG(i+1,j,k)+DIVG(i,j,k)) - UBAR(i,j,k)*(DIVG(i,j,k)+DIVG(i-1,j,k)))
	+  (VBAR(i,j+1,k)*(DIVG(i,j+1,k)+DIVG(i,j,k)) - VBAR(i,j,k)*(DIVG(i,j,k)+DIVG(i,j-1,k))) )*0.5*dtx
*/
	/*			
	- 0.50*( rhow[k+1]*W(i,j,k+1)*(DIVG(i,j,k+1)-DIVG(i,j,k)) + rhow[k]*W(i,j,k)*(DIVG(i,j,k)-DIVG(i,j,k-1)) ) * one_d_rhou[k] * DTZ(k)
		
	- dtx * 0.25*(
				(W(i+1,j,k) - W(i-1,j,k) + W(i+1,j,k+1) - W(i-1,j,k+1) ) * ONE_D_DZ(k) * 
				(rhow[k+1]*(U_AT_W(UBAR,i,j,k+1)+U_AT_W(U,i,j,k+1)) - rhow[k]*(U_AT_W(UBAR,i,j,k) - U_AT_W(U,i,j,k)))*one_d_rhou[k]
				)
	- dtx * 0.25*(
				(W(i,j+1,k) - W(i,j-1,k) + W(i,j+1,k+1) - W(i,j-1,k+1) ) * ONE_D_DZ(k) * 
				(rhow[k]+1*V_AT_W(V,i,j,k+1) - rhow[k]*V_AT_W(V,i,j,k))*one_d_rhou[k]
				)
			*/
	;

#endif

#if 0

/*********************************************************************
* Calculate terms in temperature budget
*
* il,ih,jl,jh - array indices
**********************************************************************/
void calculate_heat_budget(int il,int ih,int jl,int jh){
	

		

	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		//--------------------------------------------------------------
		// Hoizontal advection of perturbation temperature by perturbation wind,
		// i.e. the non-linear term
		//--------------------------------------------------------------
		pp_hor_adv[INDEX(i,j,k)] += - (
									   (U(i+1,j,k)*(TH(i+1,j,k)+TH(i,j,k)) - U(i,j,k)*(TH(i,j,k)+TH(i-1,j,k)))
									+  (V(i,j+1,k)*(TH(i,j+1,k)+TH(i,j,k)) - V(i,j,k)*(TH(i,j,k)+TH(i,j-1,k)))
										) * 0.5*dtx
									 + TH(i,j,k)*HCONV(U,V,i,j,k);

		//--------------------------------------------------------------
		// Horizontal advection of perturbation temperature by basic state wind
		//--------------------------------------------------------------		
		bp_hor_adv[INDEX(i,j,k)] += - ((UBAR(i+1,j,k)*(TH(i+1,j,k)+TH(i,j,k)) - UBAR(i,j,k)*(TH(i,j,k)+TH(i-1,j,k)))
									 +  (VBAR(i,j+1,k)*(TH(i,j+1,k)+TH(i,j,k)) - VBAR(i,j,k)*(TH(i,j,k)+TH(i,j-1,k))))*0.5*dtx
									 + TH(i,j,k)*HCONV(UBAR,VBAR,i,j,k);

		//--------------------------------------------------------------
		// Vertical advection of perturbation temperature by basic state wind
		//--------------------------------------------------------------		
		bp_hor_adv[INDEX(i,j,k)] += -0.5*( rhow[k+1]*WBAR(i,j,k+1)*(TH(i,j,k+1)+TH(i,j,k)) 
										 - rhow[k  ]*WBAR(i,j,k)*(TH(i,j,k)+TH(i,j,k-1)) ) * one_d_rhou[k] * DTZ(k);

		//--------------------------------------------------------------
		// Horizontal advection of basic state temperature by perturbation wind
		//--------------------------------------------------------------		
		pb_hor_adv[INDEX(i,j,k)] += - ((U(i+1,j,k)*(THBAR(i+1,j,k)+THBAR(i,j,k)) - U(i,j,k)*(THBAR(i,j,k)+THBAR(i-1,j,k)))
									 +  (V(i,j+1,k)*(THBAR(i,j+1,k)+THBAR(i,j,k)) - V(i,j,k)*(THBAR(i,j,k)+THBAR(i,j-1,k))))*0.5*dtx
									 + THBAR(i,j,k)*HCONV(U,V,i,j,k);

		//--------------------------------------------------------------
		// Vertical advection of basic state temperature by perturbation wind
		//--------------------------------------------------------------	
		pb_vert_adv[INDEX(i,j,k)] += -0.5*( rhow[k+1]*W(i,j,k+1)*(tb[k+1]-tb[k  ]) 
										  + rhow[k  ]*W(i,j,k  )*(tb[k  ]-tb[k-1])) 
											 * one_d_rhou[k] * DTZ(k);

		pb_vert_adv[INDEX(i,j,k)] += -0.5*( rhow[k+1]*W(i,j,k+1)*(THBAR(i,j,k+1)+THBAR(i,j,k))
				 						  - rhow[k  ]*W(i,j,k  )*(THBAR(i,j,k  )+THBAR(i,j,k-1))
										  ) * one_d_rhou[k] * DTZ(k);


		//--------------------------------------------------------------
		//
		//--------------------------------------------------------------		
		pp_vert_adv[INDEX(i,j,k)] += -0.5*( rhow[k+1]*W(i,j,k+1)*(TH(i,j,k+1)+TH(i,j,k)) - rhow[k]*W(i,j,k)*(TH(i,j,k)+TH(i,j,k-1)) ) * one_d_rhou[k] * DTZ(k);

		//--------------------------------------------------------------
		//
		//--------------------------------------------------------------		
		pp_conv[INDEX(i,j,k)] += -HCONV(U,V,i,j,k) * TH(i,j,k);
	
	}}}
	
}
#endif



#if 0
		//----------------------------------------------------------
		// Diabatic heating gradient vector
		//----------------------------------------------------------		
		dthetaDot.i = 0.5 * ( DIABATIC(i+1,j,k) - DIABATIC(i-1,j,k)) * one_d_dx;
		dthetaDot.j = 0.5 * ( DIABATIC(i,j+1,k) - DIABATIC(i,j-1,k)) * one_d_dy;
		dthetaDot.k = 0.5 * ( DIABATIC(i,j,k+1) - DIABATIC(i,j,k-1)) * ONE_D_DZ(k);

		//----------------------------------------------------------
		// Potential temperature gradient vector
		//----------------------------------------------------------		
		dtheta.i = 0.5 * (TH(i+1,j,k) - TH(i-1,j,k)) * one_d_dx;
		dtheta.j = 0.5 * (TH(i,j+1,k) - TH(i,j-1,k)) * one_d_dy;
		dtheta.k = 0.5 * (TH(i,j,k+1) - TH(i,j,k-1)) * ONE_D_DZ(k);

		dTheta.i = 0.5 * (THBAR(i+1,j,k) - THBAR(i-1,j,k)) * one_d_dx;
		dTheta.j = 0.5 * (THBAR(i,j+1,k) - THBAR(i,j-1,k)) * one_d_dy;				   
		dTheta.k = 0.5 * (THBAR(i,j,k+1) - THBAR(i,j,k-1) + tb[k+1] - tb[k-1] ) * ONE_D_DZ(k);

		//----------------------------------------------------------
		// Absolute vorticity vector
		//----------------------------------------------------------
		zeta.i = (W_AT_V(W,i,j+1,k) - W_AT_V(W,i,j,k)) * one_d_dy 	 - 	(V_AT_W(V,i,j,k+1) - V_AT_W(V,i,j,k)) * ONE_D_DZ(k);
		zeta.j = (U_AT_W(U,i,j,k+1) - U_AT_W(U,i,j,k)) * ONE_D_DZ(k) -  (W_AT_U(W,i+1,j,k) - W_AT_U(W,i,j,k)) * one_d_dx;
		zeta.k = (V_AT_U(V,i+1,j,k) - V_AT_U(V,i,j,k)) * one_d_dx 	 - 	(U_AT_V(U,i,j+1,k) - U_AT_V(U,i,j,k)) * one_d_dy;

		Zeta.i = (W_AT_V(WBAR,i,j+1,k) - W_AT_V(WBAR,i,j,k)) * one_d_dy 	- (V_AT_W(VBAR,i,j,k+1) - V_AT_W(VBAR,i,j,k)) * ONE_D_DZ(k);
		Zeta.j = (U_AT_W(UBAR,i,j,k+1) - U_AT_W(UBAR,i,j,k)) * ONE_D_DZ(k)  - (W_AT_U(WBAR,i+1,j,k) - W_AT_U(WBAR,i,j,k)) * one_d_dx;
		Zeta.k = (V_AT_U(VBAR,i+1,j,k) - V_AT_U(VBAR,i,j,k)) * one_d_dx 	- (U_AT_V(UBAR,i,j+1,k) - U_AT_V(UBAR,i,j,k)) * one_d_dy + FC(j);

		//----------------------------------------------------------
		//  Diffusional vorticity vector
		//----------------------------------------------------------	
		zetaF.i = (W_AT_V(WDIFF,i,j+1,k) - W_AT_V(WDIFF,i,j,k)) * one_d_dy 	 
			  -   (V_AT_W(VDIFF,i,j,k+1) - V_AT_W(VDIFF,i,j,k)) * ONE_D_DZ(k);
		
		zetaF.j = (U_AT_W(UDIFF,i,j,k+1) - U_AT_W(UDIFF,i,j,k)) * ONE_D_DZ(k) 
			  -   (W_AT_U(WDIFF,i+1,j,k) - W_AT_U(WDIFF,i,j,k)) * one_d_dx;
		
		zetaF.k = (V_AT_U(VDIFF,i+1,j,k) - V_AT_U(VDIFF,i,j,k)) * one_d_dx 	 
			  -   (U_AT_V(UDIFF,i,j+1,k) - U_AT_V(UDIFF,i,j,k)) * one_d_dy;

#endif


#if 0
		vort_hor_adv[INDEXT(i,j,k)] = 0.125* (vort_tilt[INDEXT(i,j+1,k  )]+vort_tilt[INDEXT(i+1,j+1,k  )]+
											  vort_tilt[INDEXT(i,j  ,k  )]+vort_tilt[INDEXT(i+1,j  ,k  )]+
											  vort_tilt[INDEXT(i,j+1,k+1)]+vort_tilt[INDEXT(i+1,j+1,k+1)]+
											  vort_tilt[INDEXT(i,j  ,k+1)]+vort_tilt[INDEXT(i+1,j  ,k+1)]);
		
		VORT_hor_adv[INDEXT(i,j,k)] = 0.125* (VORT_tilt[INDEXT(i,j+1,k  )]+VORT_tilt[INDEXT(i+1,j+1,k  )]+
											  VORT_tilt[INDEXT(i,j  ,k  )]+VORT_tilt[INDEXT(i+1,j  ,k  )]+
											  VORT_tilt[INDEXT(i,j+1,k+1)]+VORT_tilt[INDEXT(i+1,j+1,k+1)]+
											  VORT_tilt[INDEXT(i,j  ,k+1)]+VORT_tilt[INDEXT(i+1,j  ,k+1)]);

	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=1;k<NZ-1;k++){
		
		DIFFUSE_VARIABLE(vort_HOR_ADV);
		DIFFUSE_VARIABLE(VORT_hor_adv);
		DIFFUSE_VARIABLE(vort_ver_adv);
		DIFFUSE_VARIABLE(VORT_ver_adv);
		DIFFUSE_VARIABLE(vort_VER_ADV);
		DIFFUSE_VARIABLE(vort_hor_con);
		DIFFUSE_VARIABLE(cor_hor_con);
		DIFFUSE_VARIABLE(VORT_hor_con);
		DIFFUSE_VARIABLE(vort_tilt);
		DIFFUSE_VARIABLE(VORT_tilt);
		DIFFUSE_VARIABLE(vort_TILT);
		DIFFUSE_VARIABLE(cor_hor_adv);
		DIFFUSE_VARIABLE(vort_ufric);
		DIFFUSE_VARIABLE(vort_vfric);
		DIFFUSE_VARIABLE(vort_hor_adv);
	}}}

const double dx2 = 1.0/(dx*dx);
const double dy2 = 1.0/(dy*dy);
const double dz2 = 1.0/(dz*dz);

#define DIFFUSE_VARIABLE(v)	v[INDEXT(i,j,k)] = v[INDEXT(i,j,k)] + dt* ( 	 \
			+ 0.003 * (v[INDEXT(i+1,j,k)] - 2.0*v[INDEXT(i,j,k)]+v[INDEXT(i-1,j,k)])*dx2 \
			+ 0.003 * (v[INDEXT(i,j+1,k)] - 2.0*v[INDEXT(i,j,k)]+v[INDEXT(i,j-1,k)])*dy2 \
			+ 0.0003 * mu[k]*(mw[k+1]*v[INDEXT(i,j,k+1)] - mw[k+1]*v[INDEXT(i,j,k)] - mw[k]*v[INDEXT(i,j,k)] + mw[k]*v[INDEXT(i,j,k-1)])*dz2 \
			)	


-0.25* (
											(
											 (VBAR(i-1,j,k)+VBAR(i,j,k)) - (VBAR(i-1,j,k-1)+VBAR(i,j,k-1)) 	// dv/dz
											) * DTZ(k)
											*
						 			        (
											 (W(i,j-1,k)+W(i,j,k)) - (W(i-1,j-1,k)+W(i-1,j,k)) 	// dw/dx
											) * one_d_dx
										-
			   								(
			   								 (UBAR(i,j-1,k)+UBAR(i,j,k)) - (UBAR(i,j-1,k-1)+UBAR(i,j,k-1)) 	// du/dz
			   								) * DTZ(k)
											*
						 			        (
											 (W(i-1,j,k)+W(i,j,k)) - (W(i-1,j-1,k)+W(i,j-1,k)) 	// dw/dy
											) * one_d_dy
										   );
#endif