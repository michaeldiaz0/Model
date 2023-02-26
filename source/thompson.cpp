#include "stdafx.h"

#if HYDROSTATIC
    #define CONVERT_PRESSURE(i,j,k) PI(i,j,k)
#else
    #define CONVERT_PRESSURE(i,j,k) PI(i,j,k)/(cp*tbv[k])
#endif

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

}

/*********************************************************************
 * 
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
 * 
 * *******************************************************************/
void run_thompson_microphysics(int il,int ih, int jl, int jh, int kl, int kh,
    int nx, int ny, int nz,
    double *qv3d,double *qc,double *qi,double *qr,double *qs,double *qg,
    double *ni,double *nr,
    double *t3d,double *p3d,double *w3d,
    double *pptrain,double *pptsnow){

    double rainfall,snowfall,graupfall,icefall;
    double theta,pressure;

    int ind = 0;
    double rand0 = 0.06,rand1=1.2,rand2=6.0;
    int one = 1;
    bool re = false;


    double *nc1d = (double *)calloc(nz,sizeof(double));
    double *nwfa = (double *)calloc(nz,sizeof(double));
    double *nifa = (double *)calloc(nz,sizeof(double));
    double *t1d = (double *)calloc(nz,sizeof(double));
    double *p1d = (double *)calloc(nz,sizeof(double));
    float *a    = (float *)calloc(nz,sizeof(float));
    double *dzq = (double *)calloc(nz,sizeof(double));
    double *qv1d = (double *)calloc(nz,sizeof(double));

    double Nt_c = 100.E6;
    double naIN1 = 0.5E6;
    const double cpRd = cp/Rd;


    for(int i=il;i<ih;i++){
        for(int j=jl;j<jh;j++){

            for(int k=0;k<nz-1;k++){ dzq[k] = ZW(k+1) - ZW(k);}
            dzq[nz-1] = dzq[nz-2];

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

            ind = i*ny*nz + j*nz;

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

            for(int k=0;k<nz;k++){

                pressure = CONVERT_PRESSURE(i,j,k) + PBAR(i,j,k);
                t3d[ind+k] = t1d[k]/pressure - THBAR(i,j,k) - tb[k];

                qv3d[ind+k] = qv1d[k] - QBAR(i,j,k) - qb[k];             
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
}