#if 0
const int outfilefreq = 60;
const int stopheating = -576;

/********************************************************
* Model grid parameters
*********************************************************/
const int NX = 760;
const int NY = 520;
const int NZ = 42;

const double dx = 5000;
const double dy = 5000;
const double dz = 500;
const double dt = 30;

const double lonoffset = 70;
const double latoffset = 96;


#endif

#if 0
const int outfilefreq = 60;
const int stopheating = -576;

/********************************************************
* Model grid parameters
*********************************************************/
const int NX = 250;
const int NY = 130;
const int NZ = 42;

const double dx = 5000;
const double dy = 5000;
const double dz = 500;
const double dt = 30;

const double lonoffset = 84;
const double latoffset = 106;


#endif

#if 0
const int outfilefreq = 50;
const int stopheating = -1;

/********************************************************
* Model grid parameters
*********************************************************/
const int NX = 188;
const int NY = 188;
const int NZ = 43;

const double dx = 4000;
const double dy = 4000;
const double dz = 500;
const double dt = 24;

const double lonoffset = 86;
const double latoffset = 105;
#endif

#if 0
//Africa domain
const int outfilefreq = 180;//48;
const int stopheating = -384;//576;

const int NX = 290;//150;
const int NY = 160;//67;
const int NZ = 42;

const double dx = 25000;//75000;
const double dy = 25000;//75000;
const double dz = 500;
const double dt = 40;//450/3;

const double lonoffset = 120+10;
const double latoffset = 80+5;

const int raydampheight = NZ - 5;

#endif

#if 1
const int outfilefreq = 180;//144;//180;//36*2;//3*36;//108;//36*2;//*4;//36;//360;//72;//120;//;225;//36;//144;//72;//36;
const int stopheating = -576*2;//3*576;//1728;//576*2;//576;//-1152;//-5;//1440;//2304;//1152;//576;

/********************************************************
* Model grid parameters
*********************************************************/
const int NX = 250;//40+150;//140;//710;//70;//600;//140;//600;//280;//140;//560;//250;
const int NY = 183;//183;//110;//65;//390;//40;//365;//65;//365;//130;//65;//360;//120;
const int NZ = 40;//35;//40;//0;

const double dx = 15000;//25000;
const double dy = 15000;//25000;
const double dz = 500;
const double dt = 40;//80;//50;//40;//75/2;

const double lonoffset = 70;//66;//70;//50;//86;//84;//70;//70;//50;//70;//78;//75;//50;
const double latoffset = 96;//96;//105.5;//90;//96;//96;//102;//102;//90;//104;//104;//96;

const int raydampheight = NZ - 5;
#endif

#if 0
const int outfilefreq = 180;//144;//180;//36*2;//3*36;//108;//36*2;//*4;//36;//360;//72;//120;//;225;//36;//144;//72;//36;
const int stopheating = -576*2;//3*576;//1728;//576*2;//576;//-1152;//-5;//1440;//2304;//1152;//576;

/********************************************************
* Model grid parameters
*********************************************************/
const int NX = 95*3-5*5-5;
const int NY = 46*3-5*1-5;
const int NZ = 42;

const double dx = 25000;
const double dy = 25000;
const double dz = 500;
const double dt = 40;

const double lonoffset = 60;
const double latoffset = 92;

const int raydampheight = NZ - 5;
#endif

#if 0
const int outfilefreq = 144;
const int stopheating = 100;

const int NX = 225;//70;//109;//60;//87;//109;//82;//82;
const int NY = 40;//53;//40;
const int NZ = 40;

const double dx = 100000;//75000;//100000;
const double dy = 100000;//75000;//100000;
const double dz = 500;
const double dt = 600;

const double lonoffset = 50;//60;//50;
const double latoffset = 90;

const int raydampheight = NZ - 20;
#endif

#if 0
const int outfilefreq = 48;//36;//3*36;//108;//36*2;//*4;//36;//360;//72;//120;//;225;//36;//144;//72;//36;
const int stopheating = 100;//576;//3*576;//1728;//576*2;//576;//-1152;//-5;//1440;//2304;//1152;//576;

/********************************************************
* Model grid parameters
*********************************************************/
const int NX = 95;//113;//75;//120;//710;//70;//600;//140;//600;//280;//140;//560;//250;
const int NY = 45+1;//45;//68;//45;//390;//40;//365;//65;//365;//130;//65;//360;//120;
const int NZ = 42;//35;//40;//0;

const double dx = 75000;
const double dy = 75000;
const double dz = 500;
const double dt = 450;

const double lonoffset = 55;//65;//50;//86;//84;//70;//70;//50;//70;//78;//75;//50;
const double latoffset = 91;//105.5;//90;//96;//96;//102;//102;//90;//104;//104;//96;

const int raydampheight = NZ - 5;
#endif
