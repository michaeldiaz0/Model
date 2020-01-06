extern double sdot[NX][NY][NZ];
extern double sdotbar[NX][NY][NZ];

extern double tb3d[NX][NY][NZ];
extern double tbv3d[NX][NY][NZ];
extern double rhow3d[NX][NY][NZ];	// density at w-points
extern double qb3d[NX][NY][NZ];

extern double grhoc3d[NX][NY][NZ];	// density/Gz at cell corners
extern double grhou3d[NX][NY][NZ];	// density/Gz at u-points
extern double grhov3d[NX][NY][NZ];	// density/Gz at v-points
extern double grhos3d[NX][NY][NZ];	// density/Gz at scalar points

extern double Gx[NX][NY][NZ];	// at w-points
extern double Gy[NX][NY][NZ];	// at w-points

extern double Gz[NX][NY];
extern double Gzu[NX][NY];
extern double Gzv[NX][NY];
extern double one_d_Gz[NX][NY];
extern double one_d_Gzu[NX][NY];
extern double one_d_Gzv[NX][NY];

extern double grhoavg2d[NX][NY];

extern double sdot_at_uf[2];
extern double sdotb_at_uf[2];
extern double sdot_at_vf[2];
extern double sdotb_at_vf[2];

void sdot_interp(int,int,int);
void compute_sigma_dot(double [NX][NY][NZ],double[NX][NY][NZ]);
void compute_velocity_tend_t(double);
void advect_theta_t(double);
void compute_metric_terms();
void compute_density_terms();
void run_model_t(int,FILE *infile=NULL);

/****************************************************************************
* Kinematic fluxes of zonal velocity
*****************************************************************************/
#define XG_UU_FLUX_H ( X_UU_FLUX_H * grhos3d[i  ][j][k] )
#define XG_UU_FLUX_L ( X_UU_FLUX_L * grhos3d[i-1][j][k] )

#define YG_UV_FLUX_H ( Y_UV_FLUX_H * grhoc3d[i][j+1][k] )
#define YG_UV_FLUX_L ( Y_UV_FLUX_L * grhoc3d[i][j  ][k] )

#define ZG_UW_FLUX_H ( ( sdot_at_uf[1] * Uf[5] + (sdotb_at_uf[1]+sdot_at_uf[1]) * uf[5] ) * 0.5*( rhow3d[i-1][j][k+1]+rhow3d[i][j][k+1] ) )
#define ZG_UW_FLUX_L ( ( sdot_at_uf[0] * Uf[4] + (sdotb_at_uf[0]+sdot_at_uf[0]) * uf[4] ) * 0.5*( rhow3d[i-1][j][k  ]+rhow3d[i][j][k  ] ) )

/****************************************************************************
* Kinematic fluxes of meridional velocity
*****************************************************************************/
#define XG_VU_FLUX_H ( X_VU_FLUX_H * grhoc3d[i+1][j][k] )
#define XG_VU_FLUX_L ( X_VU_FLUX_L * grhoc3d[i  ][j][k] )

#define YG_VV_FLUX_H ( Y_VV_FLUX_H * grhos3d[i][j  ][k] )
#define YG_VV_FLUX_L ( Y_VV_FLUX_L * grhos3d[i][j-1][k] )

#define ZG_VW_FLUX_H ( ( sdot_at_vf[1] * Vf[5] + (sdotb_at_vf[1]+sdot_at_vf[1]) * vf[5] ) * 0.5*( rhow3d[i][j-1][k+1]+rhow3d[i][j][k+1] ) )
#define ZG_VW_FLUX_L ( ( sdot_at_vf[0] * Vf[4] + (sdotb_at_vf[0]+sdot_at_vf[0]) * vf[4] ) * 0.5*( rhow3d[i][j-1][k  ]+rhow3d[i][j][k  ] ) )

/****************************************************************************
* Kinematic fluxes of scalar quantities
*****************************************************************************/
#define XG_SU_FLUX_H ( X_SU_FLUX_H * grhou3d[i+1][j][k] )
#define XG_SU_FLUX_L ( X_SU_FLUX_L * grhou3d[i  ][j][k] )

#define YG_SV_FLUX_H ( Y_SV_FLUX_H * grhov3d[i][j+1][k] )
#define YG_SV_FLUX_L ( Y_SV_FLUX_L * grhov3d[i][j  ][k] )

#define ZG_SW_FLUX_H ( ( (sdot[i][j][k+1] + sdotbar[i][j][k+1]) * sf[5] + sdot[i][j][k+1] * Sf[5] ) * rhow3d[i][j][k+1] )
#define ZG_SW_FLUX_L ( ( (sdot[i][j][k  ] + sdotbar[i][j][k  ]) * sf[4] + sdot[i][j][k  ] * Sf[4] ) * rhow3d[i][j][k  ] )

/****************************************************************************
* Divergence
*****************************************************************************/
#define DRU_DX(var) ( var[i+1][j][k] * grhou3d[i+1][j][k] - var[i][j][k] * grhou3d[i][j][k] )
#define DRV_DY(var) ( var[i][j+1][k] * grhov3d[i][j+1][k] - var[i][j][k] * grhov3d[i][j][k] )
//#define DRW_DS(var) var[i][j][k+1] * grhow3d[i][j][k+1] - var[i][j][k] * grhow3d[i][j][k]

/****************************************************************************
* Derivatives
*****************************************************************************/
#define I_UPPERSUM(var) (var[i-1][j][k] + var[i][j][k] + var[i-1][j][k+1] + var[i][j][k+1])
#define I_LOWERSUM(var) (var[i-1][j][k] + var[i][j][k] + var[i-1][j][k-1] + var[i][j][k-1])
#define J_UPPERSUM(var) (var[i][j-1][k] + var[i][j][k] + var[i][j-1][k+1] + var[i][j][k+1])
#define J_LOWERSUM(var) (var[i][j-1][k] + var[i][j][k] + var[i][j-1][k-1] + var[i][j][k-1])

/*#define DS_DX(var) ( Gz[i][j] * ( var[i][j][k]*Gz[i][j] - var[i-1][j][k]*Gz[i-1][j] )/dx \*/
/*						+ 0.0625*(I_UPPERSUM(var)*I_UPPERSUM(Gx) - I_LOWERSUM(var)*I_LOWERSUM(Gx) )/dz )*/

/*#define DS_DY(var) ( Gz[i][j] * ( var[i][j][k]*Gz[i][j] - var[i][j-1][k]*Gz[i][j-1] )/dy \*/
/*						+ 0.0625*(J_UPPERSUM(var)*J_UPPERSUM(Gy) - J_LOWERSUM(var)*J_LOWERSUM(Gy) )/dz )*/

#define DS_DX(var) ( Gz[i][j] * ( var[i][j][k]*one_d_Gz[i][j] - var[i-1][j][k]*one_d_Gz[i-1][j] )/dx \
						+ 0.25*(I_UPPERSUM(var)*Gx[i][j][k+1] - I_LOWERSUM(var)*Gx[i][j][k] )/dz )

#define DS_DY(var) ( Gz[i][j] * ( var[i][j][k]*one_d_Gz[i][j] - var[i][j-1][k]*one_d_Gz[i][j-1] )/dy \
						+ 0.25*(J_UPPERSUM(var)*Gy[i][j][k+1] - J_LOWERSUM(var)*Gy[i][j][k] )/dz )


/*#define ZG_ADVECT_BASE(s) ( 0.50*( rhow3d[i][j][k+1]*sdot[i][j][k+1]*(s[i][j][k+1]-s[i][j][k  ])*Gz[i][j] \*/
/*							     + rhow3d[i][j][k  ]*sdot[i][j][k  ]*(s[i][j][k  ]-s[i][j][k-1])*Gz[i][j]	\*/
/*							    )  )*/

#define ZG_ADVECT_BASE(s) ( 0.50*( rhow3d[i][j][k+1]*w[i][j][k+1]*(s[i][j][k+1]-s[i][j][k  ])	\
							     + rhow3d[i][j][k  ]*w[i][j][k  ]*(s[i][j][k  ]-s[i][j][k-1])	\
							    ) / grhos3d[i][j][k])

/****************************************************************************
* Topography
*****************************************************************************/
#define TOPO_U ( 0.5*( topo[i-1][j][1]+topo[i][j][1] ) )
#define TOPO_V ( 0.5*( topo[i][j-1][1]+topo[i][j][1] ) )
#define TOPO_C ( 0.25 * ( topo[i-1][j][1]+topo[i-1][j-1][1]+topo[i][j-1][1]+topo[i][j][1]) )

