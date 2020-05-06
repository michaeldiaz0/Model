#pragma once

const int stopheating = -576;

class Heating 
{
	int * arr;
	public:
		double *heating;//[NZ];
		double z_heating[35];
		int i_eastlon,i_eastlat,i_westlon,i_westlat,i_width;
		double d_eastlon,d_eastlat,d_westlon,d_westlat,d_width;
		double equiv_precip_rate;
		double slope;
		Heating();
		//Heating(double,double,double,double,double,double);
		void initialize(double,double,double,double,double,double);
		void setLocationGrid(int,int,int,int);
		void setLocation(double,double,double,double);
		void scaleHeating(double);
		int changeSize(int,int,int);
		void applyHeating();
		void applyHeating_random();
		void p_applyHeating();
		void p_applyHeating_random();
		void heating_oval(int,int,double,double);
		void setHeatingRate(double);
		void shift(double,double);
		void shiftGrid(int,int);
		void printInfo(FILE * infile=NULL);
		int getSizeInGridPoints();
};

extern Heating heat;
