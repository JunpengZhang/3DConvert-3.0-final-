#include "stdafx.h"
#include "project_f.h"


//int  CoordSysParam = 54;        // 选择54坐标系统或84坐标系统
int  CoordSysParam = 84;        // 选择54坐标系统或84坐标系统

#define  PI     3.14159265358979323846
#define  PSEC   206264.80624709635515647

static double   Co=PI/180.0;    // 度转弧度的系数(单位：rad/deg)


static double  	C0=6367558.49687;
static double  	C1=  32005.7801;
static double  	C2=    133.9213;
static double  	C3=      0.7032;

static double  	K0=1.57046064172E-7;
static double  	K1=5.051773759E-3;
static double  	K2=2.9837302E-5;
static double  	K3=2.38189E-7;


static double  	a54=6378245.00000000;       // a
static double  	b54=6356863.01877305;       // b
static double   f54=1/298.3;                // f
static double  	e254=6.693421622966E-3;     // e^2
static double  	e_254=6.738525414683E-3;    // e'^2

static double  	a84=6378137.00000000;       // a
static double  	b84=6356752.31424518;       // b
static double   f84=1/298.257223563;        // f
static double  	e284=6.69437999E-3;         // e^2
static double  	e_284=6.739496742E-3;       // e'^2



/* ■  MeridianABCDE  子午线弧长公式系数计算
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
用  法: void MeridianABCDE(double e2,double &A,double &B,double &C,
					double &D,double &E,double &F,double &G);
参  数: e2   椭球第一偏心率的平方;
		A-G  子午线弧长公式系数
返回值: 无.                                                              */
void MeridianABCDE(double e2,double &A,double &B,double &C,
  double &D,double &E,double &F,double &G)
{
	A = 1+0.75*e2+45.0/64.0*e2*e2+175.0/256.0*e2*e2*e2+11025.0/16384.0*e2*e2*e2*e2;
	A += 43659.0/65536.0*pow(e2,5)+693693.0/1048576.0*pow(e2,6);
	B = 0.375*e2+15.0/32.0*e2*e2+525.0/1024.0*e2*e2*e2+2205.0/4096.0*e2*e2*e2*e2;
	B += 72765.0/131072.0*pow(e2,5)+297297.0/524288.0*pow(e2,6);
	C = 15.0/256.0*e2*e2+105.0/1024.0*e2*e2*e2+2205.0/16384.0*e2*e2*e2*e2;
	C += 10395.0/65536.0*pow(e2,5)+1486485.0/8388608.0*pow(e2,6);
	D = 35.0/3072.0*e2*e2*e2+105.0/4096.0*e2*e2*e2*e2;
	D += 10395.0/262144.0*pow(e2,5)+55055.0/1048576.0*pow(e2,6);
	E = 315.0/131072.0*pow(e2,4)+3465.0/524288.0*pow(e2,5)+99099.0/8388608.0*pow(e2,6);
	F = 639.0/1310720.0*pow(e2,5)+9009.0/5242880.0*pow(e2,6);
	G = 1001.0/8388608.0*pow(e2,6);
}



/* ■  dMeridian_X  子午线弧长计算
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
用  法: double dMeridian_X(double a,double e2,double dLat);
参  数: a      椭球长半轴;
		e2     椭球第一偏心率的平方;
		dLat   纬度(弧度).
返回值: 子午线弧长.                                                      */
double dMeridian_X(double a,double e2,double dLat)
{
	double A,B,C,D,E,F,G,X;

	MeridianABCDE(e2,A,B,C,D,E,F,G);
	X = a*(1-e2)*(A*dLat-B*sin(2*dLat)+C*sin(4*dLat)-D*sin(6*dLat)
		+E*sin(8*dLat)-F*sin(10*dLat)+G*sin(12*dLat));

	return X;
}



/* ■  dMeridian_Bf  子午线弧长反算垂足纬度Bf
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
用  法: double dMeridian_Bf(double a,double e2,double X);
参  数: a   椭球长半轴;
		e2  椭球第一偏心率的平方;
		X   子午线弧长.
返回值: 返回垂足纬度Bf(单位: 弧度).                                       */
double dMeridian_Bf(double a,double e2,double X)
{
	double A, B, C, D, E, F, G, B0, Bf;

	MeridianABCDE(e2,A,B,C,D,E,F,G);
	Bf=X/(a*(1-e2))/A;
	do{
		B0=Bf;
		Bf=X/a/(1-e2)+B*sin(2*Bf)-C*sin(4*Bf)+D*sin(6*Bf)
			-E*sin(8*Bf)+F*sin(10*Bf)-G*sin(12*Bf);
		Bf=Bf/A;
	}while(fabs(Bf-B0)>1e-13);

	return Bf;
}



/* ■  Gaussxy  高斯投影正算
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
用  法: void Gaussxy(double a,double e2,double B,double L,double L0,
					double &x,double &y);
参  数: a   椭球长半轴;
		e2  椭球第一偏心率的平方;
		B   纬度(单位: 弧度);
		L   经度(单位: 弧度);
		L0  投影带中央子午线经度(单位: 度);
		x y 计算结果, y为自然值.
返回值: 无.                                                              */
void Gaussxy(double a,double e2,double B,double L,double L0,
		 double &x,double &y)
{
	double N, X, l, tB, nB;

	L=(L<0.0)? L+2.0*PI:L;
	l=L-L0*PI/180;          //计算经差
	X=dMeridian_X(a,e2,B);  //计算子午线弧长
	N=a/sqrt(1-e2*sin(B)*sin(B));
	tB=tan(B);
	nB=sqrt(e2/(1.0-e2))*cos(B);

	x = X+N*sin(B)*cos(B)*l*l/2
		+N*sin(B)*pow(cos(B),3)*(5.0-tB*tB+9.0*nB*nB+4.0*pow(nB,4))*pow(l,4)/24.0
		+N*sin(B)*pow(cos(B),5)*(61.0-58.0*tB*tB+pow(tB,4))*pow(l,6)/720.0;
	y = N*cos(B)*l+N*pow(cos(B),3)*(1.0-tB*tB+nB*nB)*pow(l,3)/6.0
		+N*pow(cos(B),5)*(5.0-18.0*tB*tB+pow(tB,4)+14.0*nB*nB-58.0*tB*tB*nB*nB)
		*pow(l,5)/120.0;

}



/* ■  GaussBL  高斯投影反算
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
用  法: void GaussBL(double a,double e2,double x,double y,double L0,
					double &B,double &L)
参  数: a   椭球长半轴;
		e2  椭球第一偏心率的平方;
		x y 高斯平面直角坐标, y为自然值;
		L0  投影带中央子午线经度;
		B L 计算结果(单位: 弧度, 其中: L以东经为准);
返回值: 无.                                                              */
void GaussBL(double a,double e2,double x,double y,double L0,
             double &B,double &L)
{
	double t,n,b,l,M,N,Bf;

	Bf=dMeridian_Bf(a,e2,x);  //计算垂足纬度
	L0=L0*PI/180;
	M=a*(1-e2)/pow(sqrt(1-e2*sin(Bf)*sin(Bf)),3);
	N=a/sqrt(1-e2*sin(Bf)*sin(Bf));
	t=tan(Bf);
	n=sqrt(e2/(1.0-e2))*cos(Bf);

	b=-t*y*y/2/M/N
		+t*pow(y,4)/24/M/pow(N,3)*(5+3*t*t+n*n-9*t*t*n*n)
		-t*pow(y,6)/720/M/pow(N,5)*(61+90*t*t+45*pow(t,4));
	l=y/N/cos(Bf)
		-pow(y,3)/6/pow(N,3)/cos(Bf)*(1+2*t*t+n*n)
		+pow(y,5)/120/pow(N,5)/cos(Bf)*(5+28*t*t+24*pow(t,4)+6*n*n+8*t*t*n*n);

	B=Bf+b;
	L=L0+l;
	if(fabs(B*PSEC) < 5e-5)  B=fabs(B);
	if(fabs(L*PSEC) < 5e-5)  L=fabs(L);

}


///////////////////////////////////////////////////////

//  Gauss Projection & Inverse-Projection (54/84)

//      B   --  latitude,  unit: degree
//      L   --  longitude,  unit: degree
//      L0  --  central meridian,  unit: degree
//      X   --  ordinate,  unit: meter
//      Y   --  abscissa,  unit: meter

///////////////////////////////////////////////////////

int BLtoXY(double B, double L, double L0, double& X, double& Y)
{
    // 度转弧度
    B = B * Co;
    L = L * Co;

    // 根据坐标系统准备参数(54/84)
    double  a = a54;
    double  e2 = e254;

    if( CoordSysParam == 84 )
    {
        a = a84;  e2 = e284;
    }

    // 调用内部的高斯正算公式
    Gaussxy(a, e2, B, L, L0, X, Y);

    // 修正纵轴偏移量500公里
    Y += 500000.0;

    return 1;
}

int XYtoBL(double X, double Y, double L0, double& B, double& L)
{
    // 修正纵轴偏移量500公里
    Y -= 500000.0;

    // 根据坐标系统准备参数(54/84)
    double  a = a54;
    double  e2 = e254;

    if( CoordSysParam == 84 )
    {
        a = a84;  e2 = e284;
    }

    // 调用内部的高斯反算公式
    GaussBL(a, e2, X, Y, L0, B, L);

    // 弧度转度
    B = B / Co;
    L = L / Co;

    return 1;
}

///////////////////////////////////////////////////////

//  Gauss Projection & Inverse-Projection (54)

//      B   --  latitude,  unit: degree
//      L   --  longitude,  unit: degree
//      L0  --  central meridian,  unit: degree
//      X   --  ordinate,  unit: meter
//      Y   --  abscissa,  unit: meter

///////////////////////////////////////////////////////

int BLtoXY_V1(double B, double L, double L0, double& X, double& Y)
{
   double Xb,l,cosB,sinB,t,m0,N,et2;

   l=L-L0;
   B=B*Co;
   l=l*Co;

   cosB=cos(B);
   sinB=sin(B);
   t=tan(B);
   N=a54/sqrt(1.0-e254*sinB*sinB);
   m0=l*cosB;
   et2=e_254*cosB*cosB;

   Xb=C0*B-cosB*(C1*sinB+C2*sinB*sinB*sinB+C3*sinB*sinB*sinB*sinB*sinB);

   X=Xb+0.5*N*t*m0*m0+(double)(1.0/24.0)*(5.0-t*t+9.0*et2+4.0*et2*et2)*N*t*m0*m0*m0*m0
      +(double)(1.0/720.0)*(61.0-58*t*t+t*t*t*t)*N*t*m0*m0*m0*m0*m0*m0;

   Y=N*m0+(double)(1.0/6.0)*(1.0-t*t+et2)*N*m0*m0*m0+(double)(1.0/120.0)*
      (5.0-18.0*t*t+t*t*t*t+14.0*et2-58.0*et2*t*t)*N*m0*m0*m0*m0*m0;

   Y+=5.0e5;
   
   return 1;
}

int XYtoBL_V1(double X, double Y, double L0, double& B, double& L)
{
   double t,V2,Bf,u,et2,l,y,N,secB;
   double cosBf, sinBf;

   y=Y-5.0e5;
   Bf=K0*X;
   sinBf=sin(Bf);
   Bf+=cos(Bf)*(K1*sinBf-K2*sinBf*sinBf*sinBf+K3*sinBf*sinBf*sinBf*sinBf*sinBf);

   cosBf=cos(Bf);
   sinBf=sin(Bf);
   t=tan(Bf);
   et2=e_254*cosBf*cosBf;
   V2=1.0+et2;
   N=a54/sqrt(1.0-e254*sinBf*sinBf);
   u=y/N;
   secB=1.0/cosBf;

   B=Bf-0.5*V2*t*u*u+(1.0/24.0)*(5.0+3.0*t*t+et2-9.0*et2*t*t)*V2*t*u*u*u*u
      -(1.0/720.0)*(61.0+90.0*t*t+45.0*t*t*t*t)*V2*t*u*u*u*u*u*u;

   l=secB*u-(1.0/6.0)*(1.0+2.0*t*t+et2)*secB*u*u*u
     +(1.0/120.0)*(5.0+28.0*t*t+24.0*t*t*t*t+6.0*et2+8.0*et2*t*t)*secB*u*u*u*u*u;

   B=B/Co;
   L=l/Co+L0;
   
   return 1;
}

///////////////////////////////////////////////////////

//  Previous Gauss Projection & Inverse-Projection

//      B   --  latitude,  unit: DMS
//      L   --  longitude,  unit: DMS
//      L0  --  central meridian,  unit: degree
//      X   --  ordinate,  unit: meter
//      Y   --  abscissa,  unit: meter

///////////////////////////////////////////////////////

int BLtoXY_V0(double B,double L,double L0,double *X,double *Y)
{
   double Xb,l,t,m0,N,et2;

   B=DMStoDD(B);
   L=DMStoDD(L);

   l=L-L0;
   B=B*Co;
   l=l*Co;

   t=tan(B);
   N=a54/sqrt(1.0-e254*sin(B)*sin(B));
   m0=l*cos(B);
   et2=e_254*cos(B)*cos(B);

   Xb=C0*B-cos(B)*(C1*sin(B)+C2*sin(B)*sin(B)*sin(B)+C3*sin(B)*sin(B)*sin(B)*sin(B)*sin(B));

   *X=Xb+0.5*N*t*m0*m0+(double)(1.0/24.0)*(5.0-t*t+9.0*et2+4.0*et2*et2)*N*t*m0*m0*m0*m0
      +(double)(1.0/720.0)*(61.0-58*t*t+t*t*t*t)*N*t*m0*m0*m0*m0*m0*m0;

   *Y=N*m0+(double)(1.0/6.0)*(1.0-t*t+et2)*N*m0*m0*m0+(double)(1.0/120.0)*
      (5.0-18.0*t*t+t*t*t*t+14.0*et2-58.0*et2*t*t)*N*m0*m0*m0*m0*m0;

   *Y+=5.0e5;

   return 1;
}

int XYtoBL_V0(double X, double Y, double L0, double *B, double *L)
{
   double t,V2,Bf,u,et2,l,y,N,secB;

   if(Y<0)  return -1;

   y=Y-5.0e5;
   Bf=K0*X;
   Bf+=cos(Bf)*(K1*sin(Bf)-K2*sin(Bf)*sin(Bf)*sin(Bf)
      +K3*sin(Bf)*sin(Bf)*sin(Bf)*sin(Bf)*sin(Bf));

   t=tan(Bf);
   et2=e_254*cos(Bf)*cos(Bf);
   V2=1.0+et2;
   N=a54/sqrt(1.0-e254*sin(Bf)*sin(Bf));
   u=y/N;
   secB=1.0/cos(Bf);

   *B=Bf-0.5*V2*t*u*u+(1.0/24.0)*(5.0+3.0*t*t+et2-9.0*et2*t*t)*V2*t*u*u*u*u
      -(1.0/720.0)*(61.0+90.0*t*t+45.0*t*t*t*t)*V2*t*u*u*u*u*u*u;

   l=secB*u-(1.0/6.0)*(1.0+2.0*t*t+et2)*secB*u*u*u
     +(1.0/120.0)*(5.0+28.0*t*t+24.0*t*t*t*t+6.0*et2+8.0*et2*t*t)*secB*u*u*u*u*u;

   *L=DDtoDMS(l/Co+L0);
   *B=DDtoDMS(*B/Co);

   return 1;
}



///////////////////////////////////////////////// 
 
//  Convert utilities

///////////////////////////////////////////////// 

double DMStoDD(double dms)
{
   int nDeg,nMin;
   double dTmp,dSec;

   nDeg=(int)dms;
   dTmp=(double)((dms-(double)nDeg)*100.0+1.0e-12);

   nMin=(int)dTmp;
   dSec=(double)((dTmp-(double)nMin)*100.0+1.0e-12);

   return((double)nDeg+(double)nMin/60.0+dSec/3600.0);
}

double DDtoDMS(double dd)
{
   int nDeg,nMin;
   double dTmp,dSec;

   nDeg=(int)dd;
   dTmp=(dd-(double)nDeg)*60.0;
   nMin=(int)dTmp;
   dSec=(dTmp-(double)nMin)*60.0;

   return((double)nDeg+nMin/100.0+dSec/10000.0);
}

double NtoL(int N)
{
    int  l = (N <= 30) ? N*6-3 : (N-60)*6-3;
    return (double)l;
}

int LtoN(double L)
{
    if(L < 0)  L += 360;
    return ((int)(L/6) + 1);
}

double LOCM(double L)
{
    return NtoL(LtoN(L));
}

char *LB_tile_name(int L, int B)
{
    static char  name[10];

    char  E_W = 'E', N_S = 'N';
    if(L < 0)
    {
	L = -L;  E_W = 'W';
    }
    if(B < 0)
    {
	B = -B;  N_S = 'S';
    }

    sprintf(name, "%c%03d%c%02d", E_W, L, N_S, B);

    return name;
}



///////////////////////////////////////////////

// 《椭球大地测量学》教材中的高斯投影计算公式

///////////////////////////////////////////////

int BLtoXY_0(double B, double L, int L0, double& X, double& Y)
{
    double  a = 6378140.000;
    double  alpha = 298.257;
    double  e2 =  0.006694384999588;
    double  e_2 = 0.006739501819473;
    double  C = 6399596.651988;
    double  C0 = 6367452.13279;
    double  C1 = 32009.8575;
    double  C2 = 133.9602;
    double  C3 = 0.6976;

    double  rho_d = 57.295779513;

    double  b = B / rho_d;
    double  l = L - L0;
    double  t = tan(b);  
    double  t2 = t*t;
    double  cosB = cos(b);
    double  sinB = sin(b);
    double  eta2 = e_2 * cosB * cosB;
    double  N = a / sqrt(1 - e2 * sinB * sinB);
    double  m0 = l * cosB / rho_d;  
    double  m2 = m0*m0;
    double  v7 = N * m0;
    double  v8 = C0 * B / rho_d;
    double  v9 = cosB * (C1*sinB + C2*sinB*sinB*sinB + C3*sinB*sinB*sinB*sinB*sinB);
    double  v10 = N * t * m2 / 2.0;
    double  v11 = (5 - t2 + 9*eta2 + 4*eta2*eta2) * N*t*m2*m2 / 24.0;
    double  v12 = (61 - 58*t2 + t2*t2) * N*t*m2*m2*m2 / 720.0;
    double  v13 = (1 - t2 + eta2) * N*m0*m2 / 6.0;
    double  v14 = (5 - 18*t2 + t2*t2 +14*eta2*eta2 - 58*eta2*t2) * N*m0*m2*m2 / 120.0;
    X = v8 - v9 + v10 + v11 + v12;
    Y = v7 + v13 + v14;

    Y += 5.0e5;
    
    return 1;
}

int XYtoBL_0(double X, double Y, int L0, double& B, double& L)
{
    Y -= 5.0e5;

    double  a = 6378140.000;
    double  alpha = 298.257;
    double  e2 =  0.006694384999588;
    double  e_2 = 0.006739501819473;
    double  K0 = 0.15704868747276e-6;
    double  K1 = 0.005052505593;
    double  K2 = 0.000029847335;
    double  K3 = 0.0000002416;

    double  rho_d = 57.295779513;

    double  Bf0 = K0 * X;
    double  cosBf0 = cos(Bf0);
    double  sinBf0 = sin(Bf0);
    double  v5 = cosBf0 * sinBf0 * (K1 - K2*sinBf0*sinBf0 + K3*sinBf0*sinBf0*sinBf0*sinBf0);
    double  Bf = Bf0 + v5;
    double  t = tan(Bf);    
    double  t2 = t*t;
    double  cosBf = cos(Bf);
    double  sinBf = sin(Bf);
    double  eta2 = e_2 * cosBf * cosBf;
    double  V2 = 1.0 + eta2;
    double  N = a / sqrt(1 - e2 * sinBf * sinBf);
    double  y_N = Y / N;    
    double  y_N2 = y_N * y_N;
    double  y_N4 = y_N2 * y_N2;
    double  v12 = -0.5 * V2 * t * y_N2;
    double  v13 = (5 + 3*t2 + eta2 - 9*eta2*t2) * V2*t*y_N4 / 24.0;
    double  v14 = -(61 + 90*t2 + 45*t2*t2) * V2*t*y_N2*y_N4 / 720.0;
    double  v16 = y_N / cosBf;
    double  v17 = -(1 + 2*t2 + eta2) / 6.0 / cosBf * y_N*y_N2;
    double  v18 = (5 + 28*t2 + 24*t2*t2 + 6*eta2 + 8*eta2*t2) / 120.0 / cosBf * y_N*y_N4;
    double  l = (v16 + v17 + v18) * rho_d;
    B = (Bf + v12 + v13 + v14) * rho_d;
    
    L = L0 + l;
    
    return 1;
}

int BLtoXY_75(double B, double L, int L0, double& X, double& Y)
{
    double  a = 6378245.000;
    double  alpha = 298.3;
    double  e2 =  0.00669342162296594;
    double  e_2 = 0.00673852541468349;
    double  C = 6399698.90178271;
    double  C0 = 6367558.49686;
    double  C1 = 32005.79642;
    double  C2 = 133.86115;
    double  C3 = 0.7031;

    double  rho_d = 57.295779513;

    double  b = B / rho_d;
    double  l = L - L0;
    double  t = tan(b);  
    double  t2 = t*t;
    double  cosB = cos(b);
    double  sinB = sin(b);
    double  eta2 = e_2 * cosB * cosB;
    double  N = a / sqrt(1 - e2 * sinB * sinB);
    double  m0 = l * cosB / rho_d;  
    double  m2 = m0*m0;
    double  v7 = N * m0;
    double  v8 = C0 * B / rho_d;
    double  v9 = cosB * (C1*sinB + C2*sinB*sinB*sinB + C3*sinB*sinB*sinB*sinB*sinB);
    double  v10 = N * t * m2 / 2.0;
    double  v11 = (5 - t2 + 9*eta2 + 4*eta2*eta2) * N*t*m2*m2 / 24.0;
    double  v12 = (61 - 58*t2 + t2*t2) * N*t*m2*m2*m2 / 720.0;
    double  v13 = (1 - t2 + eta2) * N*m0*m2 / 6.0;
    double  v14 = (5 - 18*t2 + t2*t2 +14*eta2*eta2 - 58*eta2*t2) * N*m0*m2*m2 / 120.0;
    X = v8 - v9 + v10 + v11 + v12;
    Y = v7 + v13 + v14;

    Y += 5.0e5;
    
    return 1;
}

int XYtoBL_75(double X, double Y, int L0, double& B, double& L)
{
    Y -= 5.0e5;

    double  a = 6378245.000;
    double  alpha = 298.3;
    double  e2 =  0.00669342162296594;
    double  e_2 = 0.00673852541468349;
    double  K0 = 0.15704606412188e-6;
    double  K1 = 0.005051773759;
    double  K2 = 0.000029837302;
    double  K3 = 0.000000238189;

    double  rho_d = 57.295779513;

    double  Bf0 = K0 * X;
    double  cosBf0 = cos(Bf0);
    double  sinBf0 = sin(Bf0);
    double  v5 = cosBf0 * sinBf0 * (K1 - K2*sinBf0*sinBf0 + K3*sinBf0*sinBf0*sinBf0*sinBf0);
    double  Bf = Bf0 + v5;
    double  t = tan(Bf);    
    double  t2 = t*t;
    double  cosBf = cos(Bf);
    double  sinBf = sin(Bf);
    double  eta2 = e_2 * cosBf * cosBf;
    double  V2 = 1.0 + eta2;
    double  N = a / sqrt(1 - e2 * sinBf * sinBf);
    double  y_N = Y / N;    
    double  y_N2 = y_N * y_N;
    double  y_N4 = y_N2 * y_N2;
    double  v12 = -0.5 * V2 * t * y_N2;
    double  v13 = (5 + 3*t2 + eta2 - 9*eta2*t2) * V2*t*y_N4 / 24.0;
    double  v14 = -(61 + 90*t2 + 45*t2*t2) * V2*t*y_N2*y_N4 / 720.0;
    double  v16 = y_N / cosBf;
    double  v17 = -(1 + 2*t2 + eta2) / 6.0 / cosBf * y_N*y_N2;
    double  v18 = (5 + 28*t2 + 24*t2*t2 + 6*eta2 + 8*eta2*t2) / 120.0 / cosBf * y_N*y_N4;
    double  l = (v16 + v17 + v18) * rho_d;
    B = (Bf + v12 + v13 + v14) * rho_d;
    
    L = L0 + l;
    
    return 1;
}


/*
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

   Vincenty反解公式，已知两点的大地坐标（B1，L1）、（B2，L2），
                     求大地距离S和正反方位角A12、A21。

   参数说明：
		L1      [in]		:	第一点的经度(弧度)（-PI ~PI）
		B1      [in]		:	第一点的纬度(弧度)（-PI/2 ~PI/2）
		L2      [in]		:	第二点的经度(弧度)（-PI ~PI）
		B2 	    [in]		:	第二点的纬度(弧度)（-PI/2 ~PI/2）
		s  	    [in/out]	:	两点间距离（米）
		angle12 [in/out]	:	第一点到第二点的方位角(弧度)，以正北向为零，顺时针为正[0~2Pi]
		angle21 [in/out]	:	第二点到第一点的方位角(弧度)，以正北向为零，顺时针为正[0~2Pi]
*/

void calc_two_pnt_dist_raw(double L1, double B1, double L2, double B2, double& s, double& angle12, double& angle21)
{
	double u1, u2, sinu1, cosu1, sinu2, cosu2;
	double lamda, sinlamda, coslamda;
	double tri_derta;
	double x1, x2, x3, x4, x5;
	double C, x6 ;
	double E;
    double K1, B, A, D;
	double tanangle, temp1, temp2;
	double lamda1;

    double  LONG_RAD = (CoordSysParam == 84 ? a84 : a54);
    double  SHORT_RAD = (CoordSysParam == 84 ? b84 : b54);
    double  E_2 = (CoordSysParam == 84 ? e_284 : e_254);
    double  EARTH_f = (CoordSysParam == 84 ? f84 : f54);

	if ((fabs(B1)<=1e-20) && (fabs(B2)<=1e-20))     // 赤道上两点
	{
        s = LONG_RAD*fabs(L2-L1);
		if( L2>=L1 )
        {
			angle12 = PI/2.0;
			angle21 = 3.0 * PI/2.0;
        }
		else
        {
			angle12 = 3.0 * PI/2.0;
			angle21 = PI/2.0;
        }
	}
    else if((PI/2.0-fabs(B1)<=1e-10) || (PI/2.0-fabs(B2)<=1e-10))   // 南北极点
	{
        s = SHORT_RAD*fabs(B2-B1);
		if( B2>=B1 )
        {
			angle12 = 0.0;
			angle21 = PI;
        }
		else
        {
			angle12 = PI;
			angle21 = 0.0;
        }
	}
//    {
//        s = 0.0;
//		angle12 = 0.0;
//		angle21 = PI;
//    }
    else if( fabs(B1-B2) < (1e-10)  &&  fabs(L1-L2) < (1e-10))      // 两点距离非常近
    {
        s = 0.0;
		angle12 = 0.0;
		angle21 = PI;
    }
	else
	{
		//formua (1)
		u1 = atan((1.0 - EARTH_f) * tan(B1)); 
		sinu1 = sin(u1);
		cosu1 = cos(u1);

		u2 = atan((1.0 - EARTH_f) * tan(B2));  
		sinu2 = sin(u2);
		cosu2 = cos(u2);

		//formula (2) 
		lamda = L2 - L1;

        int  k = 0;
		do{
            if(++k > 1000)  break;
            
			//fomula (3)
			lamda1 = lamda;
			sinlamda = sin(lamda);
			coslamda = cos(lamda);

			//sin(derta)=x1  formula (4)         
			x1 = pow(pow((cosu2 * sinlamda), 2) + pow((cosu1 * sinu2 - sinu1 * cosu2 * coslamda), 2), 0.5);
		
			//cos(derta)=x2  formula (5)
			x2 = sinu1 * sinu2 + cosu1 * cosu2 * coslamda;

			//tan(derta)=x3  formula (6)
			x3 = x1 / x2;

			//sin(m)=x4   formular (7)
			x4 = cosu1 * cosu2 * sinlamda / x1;

			//(cos(m))^2=x5   formula  (8)
			x5 = 1.0 - pow(x4, 2);

			//formula (9)
			C = EARTH_f * x5 * (4.0 + EARTH_f * (4.0 - 3.0 * x5)) / 16.0;

			//cos(2*derta_m)=x6  formula(10)
			x6 = x2 - 2.0 * sinu1 * sinu2 / x5;


			//formula (11)
			E = 2.0 * pow(x6, 2) - 1.0;

			//formula (12)  derta=atan(x3)
			lamda = (L2-L1) + (1.0 - C) * EARTH_f * x4 * (atan(x3) + C * x1 * (x6 + E * C * x2));

		}while(fabs(lamda-lamda1) >(0.3e-11));

        sinlamda = sin(lamda);
        coslamda = cos(lamda);


		//formula (13)
		K1 = (pow((1.0 + E_2 * x5), 0.5) - 1.0)/(pow((1.0 + E_2 * x5), 0.5) + 1.0);

		//formula (14)
		B = K1 * (1.0 - 3.0 * pow(K1, 2) / 8.0);

		//formula (15)
		A = (1.0 + 0.25 * pow(K1, 2)) / (1.0 - K1);

		//formula (16)
		D = B * x6 * (4.0 * pow(x1, 2) - 3.0) * (2.0 * E - 1.0) / 6.0;

		//formula (17)
		tri_derta = B * x1 * (x6 + 0.25 * B * (E * x2 - D));

		//result s formula (18)
		s =((atan(x3) - tri_derta) * SHORT_RAD * A);

		//计算角度

		//formular (19)
		temp1 = cosu2 * sinlamda;

		temp2 = (cosu1 * sinu2 - sinu1 * cosu2 * coslamda);

		tanangle = temp1 / temp2;

		if(temp2 < 0){
		   angle12 = atan(tanangle) + PI;
		}
		else if(temp2 > 0 && temp1 < 0){
		   angle12 = atan(tanangle) + 2.0 * PI;
		}
		else{
			angle12 = atan(tanangle);
		}

		//formular (20)
		temp1 = -cosu1 * sinlamda;

		temp2 = (cosu2 * sinu1 - sinu2 * cosu1 * coslamda);

		tanangle = temp1 / temp2;

		if(temp2 < 0){
		   angle21 = atan(tanangle) + PI;
		}
		else if(temp2 > 0 && temp1 < 0){
		   angle21 = atan(tanangle) + 2.0 * PI;
		}
		else{
			angle21 = atan(tanangle);
		}

    }
}

/*
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

   Vincenty正解公式，已知第一点大地坐标（B1，L1）和到第二点的大地距离S、方位角A12，
                     求第二点的大地坐标（B2，L2）及反方位角A21。

   参数说明：
		L1      [in]		:	第一点的经度(弧度)
		B1      [in]		:	第一点的纬度(弧度)
		s       [in]		:	两点间距离（米）
		angle   [in/out]	:	第一点到第二点的方位角(弧度)，以正北向为零，顺时针为正[0~2Pi]
		L2 	    [in/out]	:	第二点的经度(弧度)
		B2      [in/out]	:	第二点的纬度(弧度)
		angle21 [in/out]	:	第二点到第一点的方位角(弧度)，以正北向为零，顺时针为正[0~2Pi]
*/

void calc_coord_by_distance_raw(double L1, double B1, double s, double angle, double& L2, double& B2, double& angle21)
{
	double cosderta,sinderta;
	double tri_derta;
	double A12;   
	double u1,tanu1, tanderta1,sinu1,cosu1;
	double cosA12,sinA12;
	double derta1;
	double sinm,x5;  //x5=cosm^2
	double K1,A,B;
	double derta, derta2;
	double x6, E, D;
	double temp1,temp2,lamda;
	double C;

    double  LONG_RAD = (CoordSysParam == 84 ? a84 : a54);
    double  SHORT_RAD = (CoordSysParam == 84 ? b84 : b54);
    double  E_2 = (CoordSysParam == 84 ? e_284 : e_254);
    double  EARTH_f = (CoordSysParam == 84 ? f84 : f54);

    if (fabs(s)<=1e-20)        //两点距离为零，认为为同一点
    {
        L2=L1;
        B2=B1;
        angle21=(angle<PI)? angle+PI:angle-PI ;
    }
    else if ((fabs(B1)<=1e-20) && (fabs(angle-PI/2.0)<=1e-10)) //第一点在赤道上且方位角为90 °
    {
		L2=s/LONG_RAD+L1;
		B2=0.0;
        angle21 = 1.5 * PI;
	}
	else if ((fabs(B1)<=1e-20) && (fabs(angle-PI*1.5)<=1e-10))//第一点在赤道上且方位角为270 °
    {
		L2=L1-s/LONG_RAD;
		B2=0.0;
        angle21 = 0.5 * PI;
	}
    else if(PI/2.0-fabs(B1)<=1e-10) //第一点在南（北）极点
    {
        L2=L1;
        if(fabs(PI/2.0-B1)<=1e-10)//第一点在北极点
        {
            B2=B1-s/SHORT_RAD;
            angle21=0;
        }
        else                        //第一点在南极点
        {
            B2=B1+s/SHORT_RAD;
            angle21=PI;
        }
    }
	else
    {
        A12=angle;
		sinA12=sin(A12);
		cosA12=cos(A12);

		//formula(1)
		tanu1=(1.0-EARTH_f)*tan(B1);
		u1=atan(tanu1);       
		sinu1=sin(u1);
		cosu1=cos(u1);

        if (fabs(PI/2.0-A12)<=1e-10)//如果A12＝90°
        {
            derta1=(tanu1>=0.0)?PI/2.0:-PI/2.0;
        }
        else if(fabs(1.5*PI-A12)<=1e-10)//如果A12＝270°
        {
            derta1=(tanu1>=0.0)?-PI/2.0:PI/2.0;
        }
        else
        {
		    //formula(2)
		    tanderta1=tanu1/cosA12;

		    derta1=atan(tanderta1);
		    if ((cosA12<0.0)&&(tanu1>0.0))
		    {
			    derta1+=PI;
		    }
		    else if ((cosA12<0.0)&&(tanu1<0.0))
		    {
			    derta1+=-PI;
		    }
        }

		//formula(3)
		sinm=cosu1*sinA12;
		x5=1.0-pow(sinm,2);

		//formula(4)
		K1 = (pow((1.0 + E_2 * x5), 0.5) - 1.0)/(pow((1.0 + E_2 * x5), 0.5) + 1.0);

		//formula(5)
		A = (1.0 + 0.25 * pow(K1, 2)) / (1.0 - K1);

		//formula(6)
		B = K1 * (1.0 - 3.0 * pow(K1, 2) / 8.0);

		//formula(7)
		derta=s/(SHORT_RAD*A);

        int  k = 0;
		do{
            if(++k > 1000)  break;

			derta2=derta;

			cosderta=cos(derta);
			sinderta=sin(derta);

			//formular(8)
			//x6=cos(2dertam)
			//2dertam=2derta1+derta;

			//formular(9)
			x6=cos(2.0*derta1+derta);
			E=2.0*pow(x6,2)-1.0;

			//formular(10)
			D = B * x6 * (4.0 * pow(sinderta, 2) - 3.0) * (2.0 * E - 1.0) / 6.0;

			//formular(11)
			tri_derta = B * sinderta * (x6 + 0.25 * B * (E * cosderta - D));

			//formular(12)
			derta=s/(SHORT_RAD*A)+tri_derta;


		} while (fabs(derta - derta2)>(0.3e-11));

        cosderta=cos(derta);
        sinderta=sin(derta);

		//formular(13)
		B2=atan((sinu1*cosderta+cosu1*sinderta*cosA12)/((1-EARTH_f)*pow(((1-x5)+pow((sinu1*sinderta-cosu1*cosderta*cosA12),2)),0.5)));

		//formular(14)
		temp1=sinderta*sinA12;
		temp2=cosu1*cosderta-sinu1*sinderta*cosA12;
		lamda=atan(temp1/temp2);

		if ((temp2<0)&&(temp1>0))
		{
			lamda+=PI;
		}
		if ((temp2<0)&&(temp1<0))
		{
			lamda+=-PI;
		}

		//formular(15)
		C=EARTH_f*x5*(4.0+EARTH_f*(4.0-3.0*x5))/16.0;

		//formular(16)
		L2=L1+lamda-(1-C)*EARTH_f*sinm*(derta+C*sinderta*(x6+E*C*cosderta));

		if (fabs(L2)>PI)
		{
			if (L2>0)
				L2-=2*PI;
			else if (L2<0)
				L2+=2*PI;
		}

		//formular(17)
        double temp1;
        double temp2;
        double tanangle;

		temp1 = -sinm;
		temp2 = sinu1*sinderta-cosu1*cosderta*cosA12;
		tanangle = temp1 / temp2;

        if(temp2 < 0)
        {
            angle21 = atan(tanangle) + PI;
        }
		else
        {
            if(temp1 > 0)
                angle21 = atan(tanangle) + 2 * PI;
		    else
			    angle21 = atan(tanangle);
		}
	}
}

/////////////////////////////////////////////////////////////////////////

// 功能：已知两点大地坐标求其距离和方位角
// 参数定义：
//		L1    [in]		:	第一点的经度(度）（-PI ~PI）
//		B1    [in]		:	第一点的纬度(度）（-PI/2 ~PI/2）
//		L2    [in]		:	第二点的经度(度）（-PI ~PI）
//		B2 	  [in]		:	第二点的纬度(度）（-PI/2 ~PI/2）
//		s  	  [in/out]	:	两点间距离（Km）
//		angle [in/out]	:	第一点到第二点的方位角（度），以东向为零，逆时针为正[0~360]

/////////////////////////////////////////////////////////////////////////

void calc_two_pnt_dist(double L1, double B1, double L2, double B2, double& s, double& angle, double& A21)
{
	double ratio = PI / 180.0;

    L1 *= ratio;
    B1 *= ratio;
    L2 *= ratio;
    B2 *= ratio;

    calc_two_pnt_dist_raw(L1, B1, L2, B2, s, angle, A21);

    s /= 1000.0;                        // 转换到公里

    angle /= ratio;                     // 转换到度
    angle = 90.0 - angle;               // 东向基准方向角
    if(angle < 0)  angle += 360.0;      // 0 ~ 360

    A21 /= ratio;                       // 转换到度
    A21 = 90.0 - A21;                   // 东向基准方向角
    if(A21 < 0)  A21 += 360.0;          // 0 ~ 360
}

/////////////////////////////////////////////////////////////////////////

// 功能：已知一点大地坐标和到第二点的距离和方位角，求第二点的大地坐标
// 参数定义：
//		L1    [in]		:	第一点的经度（度）
//		B1    [in]		:	第一点的纬度（度）
//		s     [in]		:	两点间距离（Km）
//		angle [in]		:	第一点到第二点的方位角（度），以东向为零，逆时针为正[0~360]
//		L2 	  [in/out]	:	第二点的经度（度）
//		B2    [in/out]	:	第二点的纬度（度）

/////////////////////////////////////////////////////////////////////////

void calc_coord_by_distance(double L1, double B1, double s, double angle, double& L2, double& B2, double& A21)
{
	double ratio = PI / 180.0;

    angle = 90.0 - angle;           // 大地方位角
    while(angle < 0)  angle += 360.0;  // 0 ~ 360
	while(angle > 360)  angle -= 360.0;  // 0 ~ 360
    angle *= ratio;                 // 转换为弧度

    L1 *= ratio;                    // 转换为弧度
    B1 *= ratio;                    // 转换为弧度

    s *= 1000.0;                    // 转换到米

    calc_coord_by_distance_raw(L1, B1, s, angle, L2, B2, A21);

    L2 /= ratio;                    // 转换到度
    B2 /= ratio;                    // 转换到度

    A21 /= ratio;                   // 转换到度
    A21 = 90.0 - A21;               // 东向基准方向角
    if(A21 < 0)  A21 += 360.0;      // 0 ~ 360
}
