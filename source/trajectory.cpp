#include "stdafx.h"
#include "interpolate.h"
#include "microphysics.h"

#define ALLOC(v,size) v = (double *)calloc(size,sizeof(double))

#define INTERP_PARCEL_TO_U(x,y,z) ( (frac1)*array_interpolate(u1, x+0.5, y, z) + (frac2)*array_interpolate(u2, x+0.5, y, z) )
#define INTERP_PARCEL_TO_V(x,y,z) ( (frac1)*array_interpolate(v1, x, y+0.5, z) + (frac2)*array_interpolate(v2, x, y+0.5, z) )
#define INTERP_PARCEL_TO_W(x,y,z) ( (frac1)*array_interpolate(w1, x, y, z+0.5) + (frac2)*array_interpolate(w2, x, y, z+0.5) )

#define SAT_VAP_WAT(t) ( 611.2 * exp(17.67 * (t-273.15) / (t - 29.65) ) )
#define SAT_VAP_ICE(t) ( 611.2 * exp(21.8745584 * (t-273.15) / (t - 7.66) ) )
#define SAT_MIX_RATIO(e,p) ( 0.62197 * e / (p-e) )

const int fileDT = 1*60*60;
const int startTime = 13;
const int endTime = 0;

double *q_zgrad;
double *t_zgrad;
double *q_ygrad;
double *t_ygrad;

void print_output(struct parcel p,int i,int t);

/*************************
* Runge-Kutta approximations
* for each coordinate
**************************/
struct k_values {

	double x;
	double y;
	double z;
};

/*************************
* Single parcel
**************************/
struct parcel {

	// initial location
	double is;
	double js;
	double ks;
	
	// current location
	double ip;
	double jp;
	double kp;
	
	// is the parcel still within
	// the model domain?
	bool isActive;
};

/*********************************************************************
* 
**********************************************************************/
struct parcel initialize_parcel(double is,double js,double ks){
	
	parcel p;
	
	p.is = is;
	p.js = js;
	p.ks = ks;
	
	p.ip = is;
	p.jp = js;
	p.kp = ks;

	p.isActive = true;
	
	return p;
}

/*********************************************************************
* 
**********************************************************************/
void add_arrays(double *array1,double *array2,int size){

	for(int i=0;i<size;i++){ array1[i] += array2[i];}

}

/*********************************************************************
* 
**********************************************************************/
void add_arrays2(double *array1,double *array2,double *array3){

	int ind;

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
		
		ind = INDEX(i,j,k);
		
		array1[ind] += (array2[ind] + array3[k]);
	}}}

}

/*********************************************************************
* 
**********************************************************************/
void add_arrays3(double *array1,double *array3){

	int ind;

	for(int i=0;i<NX;i++){
	for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
		
		ind = INDEX(i,j,k);
		
		array1[ind] += array3[k];
	}}}

}

/*********************************************************************
* 
**********************************************************************/
double array_interpolate(double *var,double is,double js,double ks){

	double var_matrix[4][4];
	double var_z[4];
	
	int imid = (int)is;
	int jmid = (int)js;
	int kmid = (int)ks;
	
	
	double imu = is - (double)imid;
	double jmu = js - (double)jmid;
	double kmu = ks - (double)kmid;
	
	for(int k=-1;k<3;k++){
		
		for(int i=-1;i<3;i++){
		for(int j=-1;j<3;j++){
	
			var_matrix[i+1][j+1] = var[INDEX(i+imid,j+jmid,k+kmid)];
		}}
		
		var_z[k+1] = BicubicInterpolate2(var_matrix,imu,jmu);
	}
	
	return CubicInterpolate(var_z[0],var_z[1],var_z[2],var_z[3],zu[kmid-1],zu[kmid],zu[kmid+1],zu[kmid+2],kmu);
}


/*********************************************************************
* 
**********************************************************************/
bool isSafeFromBoundaries(double i,double j,double k){
	
	int buffer = 10;
	int lower = 1;
	int upper = 5;
	
	if(i >= (double)(NX-1-buffer) || i < (double)buffer){ return false;}
	if(j >= (double)(NY-1-buffer) || j < (double)buffer){ return false;}	
	if(k >= (double)(NZ-1-upper)  || k < (double)lower ){ return false;}

	return true;
}

/*********************************************************************
* 
* u1,v1,w1 - velocity field at beginning of interval
* u2,v2,w2 - velocity field at end of interval
* dt - time between interval
**********************************************************************/
bool RK4_Step(struct parcel *p,double *u1,double *v1,double *w1,double *u2,double *v2,double *w2,double dt,int intervals,double r){
	
	k_values k1,k2,k3,k4;
	const double one_d_six = 1.0/6.0;

	static int counter = 0;

	double zh; 
	//double kh = convert_z_to_k(zh,50,1,dz);

	double frac1,frac2,fracInterval;

	double x,y,z;

	dt /= (double)intervals;

	double interval = 1.0/intervals;

	double dtx = dt/dx;
	double dty = dt/dy;

	if( !isSafeFromBoundaries(p->ip,p->jp,p->kp) ){
		
		p->isActive = false;
		return false;
	}

	//--------------------------------------------------------------------------------------
	// Run RK4 integration on each subinterval
	//--------------------------------------------------------------------------------------
	for(int i=0;i<intervals;i++){

		fracInterval = (double)i / (double)intervals;

		zh = convert_k_to_z(p->kp,height_lowest_level,index_lowest_level,dz);

		//--------------------------------------------------------------------------------------
		// First step
		//--------------------------------------------------------------------------------------
		x = p->ip;
		y = p->jp;
		z = p->kp;
		
		frac1 = 1.0 - fracInterval;
		frac2 = 1.0 - frac1;

		k1.x = dtx * INTERP_PARCEL_TO_U(x,y,z);
		k1.y = dty * INTERP_PARCEL_TO_V(x,y,z);
		k1.z = dt  * INTERP_PARCEL_TO_W(x,y,z);

		//--------------------------------------------------------------------------------------
		// Second step
		//--------------------------------------------------------------------------------------
		x = p->ip + r*0.5*k1.x;
		y = p->jp + r*0.5*k1.y;
		z = convert_z_to_k(zh + r*0.5*k1.z,height_lowest_level,index_lowest_level,dz);
		
		frac1 -= interval / 2.0;
		frac2 = 1.0 - frac1; 
		
		k2.x = dtx * INTERP_PARCEL_TO_U(x,y,z);
		k2.y = dty * INTERP_PARCEL_TO_V(x,y,z);
		k2.z = dt  * INTERP_PARCEL_TO_W(x,y,z);
		
		//--------------------------------------------------------------------------------------
		// Third step
		//--------------------------------------------------------------------------------------	
		x = p->ip + r*0.5*k2.x;
		y = p->jp + r*0.5*k2.y;
		z = convert_z_to_k(zh + r*0.5*k2.z,height_lowest_level,index_lowest_level,dz);

		k3.x = dtx * INTERP_PARCEL_TO_U(x,y,z);
		k3.y = dty * INTERP_PARCEL_TO_V(x,y,z);
		k3.z = dt *  INTERP_PARCEL_TO_W(x,y,z);

		//--------------------------------------------------------------------------------------
		// Fourth step
		//--------------------------------------------------------------------------------------
		x = p->ip + r*k3.x;
		y = p->jp + r*k3.y;
		z = convert_z_to_k(zh + r*k3.z,height_lowest_level,index_lowest_level,dz);
	
		frac1 -= interval / 2.0;
		frac2 = 1.0 - frac1;
	
		k4.x = dtx * INTERP_PARCEL_TO_U(x,y,z);
		k4.y = dty * INTERP_PARCEL_TO_V(x,y,z);
		k4.z = dt  * INTERP_PARCEL_TO_W(x,y,z);
	
		//--------------------------------------------------------------------------------------
		// Final value
		//--------------------------------------------------------------------------------------	
		p->ip = p->ip 			  + r*one_d_six * (k1.x + 2.0*k2.x + 2.0*k3.x + k4.x);
		p->jp = p->jp 			  + r*one_d_six * (k1.y + 2.0*k2.y + 2.0*k3.y + k4.y);
		p->kp = convert_z_to_k(zh + r*one_d_six * (k1.z + 2.0*k2.z + 2.0*k3.z + k4.z),height_lowest_level,index_lowest_level,dz);
		
		print_output(p[0],counter,0);
		counter++;
	}

	return true;
}

/*********************************************************************
* 
**********************************************************************/
void get_thermo_data(const char * myfilename,int t){

	get_data_at_time(myfilename,"theta",t,ths);
	get_data_at_time(myfilename,"pi",t,pis);
	get_data_at_time(myfilename,"qc",t,qcs);
	get_data_at_time(myfilename,"qv" ,t,qvs);
}

/*********************************************************************
* 
**********************************************************************/
void print_output(struct parcel p,int i,int t){

	//----------------------------------------------------
	// Get output variables
	//----------------------------------------------------
	double p_u,p_v,p_w,p_theta,p_qv,p_qc,b_qv,b_theta,b_pi,p_pi,zheight;
	double pressure,qvsat,esl,pd,temperature,vapor,theta;
	double cpRd = cp / Rd;
	double ti_ygrad,ti_zgrad,qi_ygrad,qi_zgrad;

	zheight = convert_k_to_z(p.kp,100,1,dz);

	p_u = array_interpolate(us,p.ip+0.5,p.jp,p.kp);
	p_v = array_interpolate(vs,p.ip,p.jp+0.5,p.kp);
	p_w = array_interpolate(ws,p.ip,p.jp,p.kp+0.5);
	p_theta = array_interpolate(ths,p.ip,p.jp,p.kp);
	p_qv = array_interpolate(qvs,p.ip,p.jp,p.kp);
	p_qc = array_interpolate(qcs,p.ip,p.jp,p.kp);
	p_pi = array_interpolate(pis,p.ip,p.jp,p.kp);

	b_pi = array_interpolate(m_pbar,p.ip,p.jp,p.kp);
	b_theta = array_interpolate(m_thbar,p.ip,p.jp,p.kp);
	b_qv = array_interpolate(m_qbar,p.ip,p.jp,p.kp);

	ti_ygrad = array_interpolate(t_ygrad,p.ip,p.jp+0.5,p.kp); 
	ti_zgrad = array_interpolate(t_zgrad,p.ip,p.jp,p.kp+0.5);
	qi_ygrad = array_interpolate(q_ygrad,p.ip,p.jp+0.5,p.kp);
	qi_zgrad = array_interpolate(q_zgrad,p.ip,p.jp,p.kp+0.5);

	theta = p_theta + b_theta;						// full potential temperature is base state plus perturbation		
	pressure = p_pi/(cp*b_theta) + b_pi;			// full pressure is base state plus perturbation
	vapor = p_qv + b_qv;		// full vapor is base state plus perturbation
	temperature = theta*pressure;					// actual temperature
	pd = p0*pow(pressure,cpRd);						// dimensional pressure

	esl = SAT_VAP_WAT(temperature);
	qvsat = SAT_MIX_RATIO(esl,pd);

	//double thetaE = theta * (1+Lv*qvsat / (cp*temperature));

	double MSE = cp*temperature + Lv*vapor + grav*convert_k_to_z(p.kp,100,1,dz);

	double cape = get_CAPE(round(p.ip),round(p.jp),round(p.kp),vapor,theta);

	printf("%d\t%d\t%.2f\t%.2f\t%.2f\t% .2f\t% .2f\t% .2f\t% .2f\t% .2f\t% .2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.1f\t%.1f\t%.0f\t%.0f\t% e\t% e\t% e\t% e\n",i,t,
	p.ip*dx/meters_per_degree + outLons[0],
	p.jp*dy/meters_per_degree + outLats[0],
	zheight,p_u,p_v,p_w*100.0,p_theta,p_qv*1000,p_qc*1000,
	b_theta,b_qv*1000,p_theta+b_theta,(p_qv+b_qv)*1000,100*(p_qv+b_qv)/qvsat,0.01*pd,cape,MSE,ti_ygrad,ti_zgrad,qi_ygrad*1000,qi_zgrad*1000
		);
}

/*********************************************************************
* 
**********************************************************************/
void trajectory_from_file(const char * myfilename){

	double llat = 21;
	double hlat = 21.2;
	double llon = 77.9;
	double hlon = 78.0;

	int size = NX*NY*NZ;

	get_data(myfilename,"zu",NZ,&zu[0]);
	get_data(myfilename,"zu",NZ,&zsu[0]);
	get_data(myfilename,"tb",NZ,&tb[0]);
	get_data(myfilename,"qb",NZ,&qb[0]);

	get_data(myfilename,"lat",NY,&outLats[0]);
	get_data(myfilename,"lon",NX,&outLons[0]);

	ALLOC(us,size); ALLOC(vs,size); ALLOC(ws,size);
	ALLOC(ums,size); ALLOC(vms,size); ALLOC(wms,size);
	ALLOC(m_ubar,size); ALLOC(m_vbar,size); ALLOC(m_wbar,size);
	
	ALLOC(ths,size); ALLOC(pis,size); ALLOC(qvs,size); ALLOC(qcs,size);
	ALLOC(m_thbar,size); ALLOC(m_qbar,size); ALLOC(m_pbar,size);

	ALLOC(q_zgrad,size); ALLOC(t_zgrad,size); ALLOC(q_ygrad,size); ALLOC(t_ygrad,size);

	get_data(myfilename,"ubar",size,m_ubar);
	get_data(myfilename,"vbar",size,m_vbar);
	get_data(myfilename,"wbar",size,m_wbar);
	get_data(myfilename,"thbar",size,m_thbar);
	get_data(myfilename,"qbar",size,m_qbar);
	get_data(myfilename,"pbar",size,m_pbar);

	get_data_at_time(myfilename,"u-wind",startTime,ums);
	get_data_at_time(myfilename,"v-wind",startTime,vms);
	get_data_at_time(myfilename,"w-wind",startTime,wms);

	add_arrays(ums,m_ubar,size);
	add_arrays(vms,m_vbar,size);
	add_arrays(wms,m_wbar,size);
	
	//parcel p = initialize_parcel(NX/2,NX/2,10);

	int il = get_point_from_lon(llon);
	int ih = get_point_from_lon(hlon);
	int jl = get_point_from_lat(llat);
	int jh = get_point_from_lat(hlat);
	int kl = 8;
	int kh = 9;

	int length = (ih-il)*(jh-jl)*(kh-kl);
	int counter = 0;
	
	//----------------------------------------------------
	// Initialize parcels
	//----------------------------------------------------
	parcel *p = (parcel *)calloc(length,sizeof(parcel));
	
	for(int i=il;i<ih;i++){
	for(int j=jl;j<jh;j++){
	for(int k=kl;k<kh;k++){
	
		p[counter] = initialize_parcel(i,j,k);
		
		counter++;
	}}}

 	get_thermo_data(myfilename,startTime);
	
	add_arrays3(m_thbar,&tb[0]);
	add_arrays3(m_qbar,&qb[0]);

	//----------------------------------------------------
	// Calculate gradients
	//----------------------------------------------------
	for(int i=1;i<NX-1;i++){
	for(int j=1;j<NY-1;j++){
	for(int k=1;k<NZ-1;k++){

		t_zgrad[INDEX(i,j,k)] = //0.5*(m_thbar[INDEX(i,j,k+1)] - m_thbar[INDEX(i,j,k  )]) / (ZU(k+1)-ZU(k  )) +
								//0.5*
									(m_thbar[INDEX(i,j,k  )] - m_thbar[INDEX(i,j,k-1)]) / (ZU(k  )-ZU(k-1));
									
		q_zgrad[INDEX(i,j,k)] = //0.5*(m_qbar[INDEX(i,j,k+1)] - m_qbar[INDEX(i,j,k  )]) / (ZU(k+1)-ZU(k  )) +
								//0.5*
									(m_qbar[INDEX(i,j,k  )] - m_qbar[INDEX(i,j,k-1)]) / (ZU(k  )-ZU(k-1));

		t_ygrad[INDEX(i,j,k)] = (m_thbar[INDEX(i,j,k)] - m_thbar[INDEX(i,j-1,k)]) / dy;																																						
		q_ygrad[INDEX(i,j,k)] = (m_qbar[INDEX(i,j,k)] - m_qbar[INDEX(i,j-1,k)]) / dy;
															
	}}}

	//----------------------------------------------------
	// Print initial values
	//----------------------------------------------------
	get_data_at_time(myfilename,"u-wind",startTime,us);
	get_data_at_time(myfilename,"v-wind",startTime,vs);
	get_data_at_time(myfilename,"w-wind",startTime,ws);

	add_arrays(us,m_ubar,size);
	add_arrays(vs,m_vbar,size);
	add_arrays(ws,m_wbar,size);
	
	for(int i=0;i<length;i++){ print_output(p[i],i,startTime);}

	//----------------------------------------------------
	// Calculate forward trajectories
	//----------------------------------------------------
	if(startTime < endTime){

		for(int t=startTime+1;t<endTime;t++){

			get_data_at_time(myfilename,"u-wind",t,us );
			get_data_at_time(myfilename,"v-wind",t,vs );
			get_data_at_time(myfilename,"w-wind",t,ws );

			get_thermo_data(myfilename,t);

			add_arrays(us,m_ubar,size);
			add_arrays(vs,m_vbar,size);
			add_arrays(ws,m_wbar,size);
		
			for(int i=0;i<length;i++){
			
				if(p[i].isActive){

					RK4_Step(&p[i],ums,vms,wms,us,vs,ws,fileDT,100,1.0);

					print_output(p[i],i,t);
				}
			}
			switch_array(&us,&ums);
			switch_array(&vs,&vms);
			switch_array(&ws,&wms);	
		}
	//----------------------------------------------------
	// Calculate backward trajectories
	//----------------------------------------------------		
	} else {

		for(int t=startTime-1;t>=endTime;t--){

			get_data_at_time(myfilename,"u-wind",t,us );
			get_data_at_time(myfilename,"v-wind",t,vs );
			get_data_at_time(myfilename,"w-wind",t,ws );

			get_thermo_data(myfilename,t);

			add_arrays(us,m_ubar,size);
			add_arrays(vs,m_vbar,size);
			add_arrays(ws,m_wbar,size);

			for(int i=0;i<length;i++){
			
				if(p[i].isActive){

					RK4_Step(&p[i],ums,vms,wms,us,vs,ws,fileDT,100,-1.0);

					//print_output(p[i],i,t);
				}
			}
			switch_array(&us,&ums);
			switch_array(&vs,&vms);
			switch_array(&ws,&wms);	
		}
	}

}