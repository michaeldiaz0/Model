#include "stdafx.h"
#include "budgets.h"

/************************************************************************
 * External Fortran subroutines
 * 
 * 
 * **********************************************************************/
extern"C" {

    void thompson_init(
        int *is_aerosol_aware_in,
        int *mpicomm,
        int *mpirank,
        int *mpiroot,
        int *threads,
        const char *errmsg,
        int *errflg);

    void mp_thompson(
        double *qv1d,double *qc1d,double *qi1d,double *qr1d,double *qs1d,double *qg1d,
        double *ni1d,double *nr1d,double *nc1d,double *nwfa1d,double *nifa1d,double *t1d,double *p1d,
        double *w1d,double *dzq,double *pptrain,double *pptsnow,double *pptgraul,double *pptice,
        double *rand1, double *rand2, double *rand3,int *kts,int *kte,double *dt,int *ii,int *jj,bool *,bool *,int *,
        float *,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *,
        float *,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *,
        float *,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *
    );

    void calc_refl10cm(
        double *qv1d, double *qc1d, double *qr1d, double *nr1d,double *qs1d, double *qg1d,
        double *t1d, double *p1d, double *dBZ, double *rand1, int *kts, int *kte, int *ii, int *jj, bool *melti,
        double *vt_dBZ, bool *first_time_step
    );

}

/*********************************************************************
 * 
 * 
 * *******************************************************************/
void init_thompson_microphysics(int mpirank){

    int is_aerosol_aware_in = 0;
    int mpicomm = 0;
    int mpiroot = 0;
    int threads = 1;
    const char *errmsg = "E";
    int errflg = 0;

    thompson_init(&is_aerosol_aware_in,&mpicomm,&mpirank,&mpiroot,&threads,errmsg,&errflg);
}

/*********************************************************************
 * 
 * 
 * *******************************************************************/
void run_thompson_microphysics(
    int il,int ih, int jl, int jh, int kl, int kh,int nx, int ny, int nz,
    double *qv3d, double *qc, double *qi, double *qr, double *qs, double *qg, double *ni, double *nr,
    double *t3d, double *p3d, double *w3d, double *pptrain, double *pptsnow, bool do_radar
){

    double rainfall,snowfall,graupfall,icefall;
    double theta,pressure;

    int ind = 0;
    double rand0 = 0.06,rand1=1.2,rand2=6.0;
    int one = 1;
    bool re = false;
    bool melti = true;
    bool first_time_step = false;
    double *vt_dBZ = (double *)calloc(nz,sizeof(double));


    double *nc1d = (double *)calloc(nz,sizeof(double));
    double *nwfa = (double *)calloc(nz,sizeof(double));
    double *nifa = (double *)calloc(nz,sizeof(double));
    double *t1d =  (double *)calloc(nz,sizeof(double));
    double *p1d =  (double *)calloc(nz,sizeof(double));
    float  *a    = (float  *)calloc(nz,sizeof(float));
    double *dzq =  (double *)calloc(nz,sizeof(double));
    double *qv1d = (double *)calloc(nz,sizeof(double));

    double Nt_c = 100.E6;
    double naIN1 = 0.5E6;
    const double cpRd = cp/Rd;
    double qv_sat;


    for(int i=il;i<ih;i++){
        for(int j=jl;j<jh;j++){

            for(int k=0;k<nz-1;k++){ 
                dzq[k] = ZW(k+1) - ZW(k);
            }
            dzq[nz-1] = dzq[nz-2];

            //--------------------------------------------------------------
            // Processed fields for use in microphysics
            //--------------------------------------------------------------
            for(int k=0;k<nz;k++){

                theta = t3d[INDEX(i,j,k)] + THBAR(i,j,k) + tb[k];   // full potential temperature is base state plus perturbation		
                pressure = CONVERT_PRESSURE(i,j,k) + PBAR(i,j,k);	// full pressure is base state plus perturbation
                qv1d[k] = qv3d[INDEX(i,j,k)] + QBAR(i,j,k) + qb[k];	// full vapor is base state plus perturbation
                t1d[k] = theta*pressure;					        // actual temperature
                p1d[k] = p0*pow(pressure,cpRd);						// dimensional pressure

                nc1d[k] = Nt_c/rhou[k];
                nwfa[k] = 11.1E6;
                nifa[k] = naIN1*0.01;
            }

            ind = i*ny*nz + j*nz; // index into starting locations at k = 0

            //--------------------------------------------------------------
            // Call external subroutine for Thompson microphysics
            //--------------------------------------------------------------
            mp_thompson(
                qv1d,&qc[ind],&qi[ind],&qr[ind],&qs[ind],&qg[ind],
                &ni[ind],&nr[ind],nc1d,
                nwfa,nifa,t1d,p1d,&w3d[ind],dzq,
                &rainfall,&snowfall,&graupfall,&icefall,
                &rand0,&rand1,&rand2,&kl,&kh,&dt,&i,&j,
                &re,&re,&one,
                a,a,a,a,a,a,a,a,a,a,a,a,
                a,a,a,a,a,a,a,a,a,a,a,a,
                a,a,a,a,a,a,a,a,a,a,a,a,a
            );

            //--------------------------------------------------------------
            // Convert back to the variables we need
            //--------------------------------------------------------------
            for(int k=0;k<nz;k++){

                pressure = CONVERT_PRESSURE(i,j,k) + PBAR(i,j,k);
                t3d[ind+k] = t1d[k]/pressure - THBAR(i,j,k) - tb[k];

                qv3d[ind+k] = qv1d[k] - QBAR(i,j,k) - qb[k];

                qv_sat = get_qvsat_mixed(t1d[k],p1d[k]);

                if(qv1d[k]>=qv_sat){
                    isSaturated[ind+k] = true;
                } else {
                    isSaturated[ind+k] = false;
                }
            }

            accRain[i*ny+j] += rainfall;
            accSnow[i*ny+j] += snowfall + graupfall + icefall;

            if(do_radar){
                calc_refl10cm(
                        qv1d, &qc[ind], &qr[ind], &nr[ind], &qs[ind], &qg[ind],
                        t1d, p1d, &dBZ[ind], &rand1, &kl, &kh, &i, &j, &melti,
                        vt_dBZ, &first_time_step);
            }
    

            //--------------------------------------------------------------
            // If budgets requested, store temperature and moisture tendencies
            //--------------------------------------------------------------
            if(HEAT_BUDGET || PE_BUDGET || PV_BUDGET){

                for(int k=0;k<nz;k++){
			        m_diabatic[ind+k] += t3d[ind+k] - (t1d[k]/pressure - THBAR(i,j,k) - tb[k]);
                }
            }
            if(MOISTURE_BUDGET){
                for(int k=0;k<nz;k++){
			        q_diabatic[ind+k] += qv3d[ind+k] - (qv1d[k] - QBAR(i,j,k) - qb[k]);
                }
            }
		}
    }

    free(nc1d);
    free(nwfa);
    free(nifa);
    free(t1d);
    free(qv1d);
    free(p1d);
    free(dzq);
    free(a);
    free(vt_dBZ);
}
