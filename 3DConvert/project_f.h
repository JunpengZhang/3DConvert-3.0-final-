#ifndef  _Project_Include
#define  _Project_Include

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>

///////////////////////////////////////////////////////

// ѡ������ϵͳ(54 �� ��ʾ54����ϵͳ��84 �� ��ʾ84����ϵͳ)

///////////////////////////////////////////////////////

extern int  CoordSysParam;


///////////////////////////////////////////////////////

// Vincenty���⡢������㹫ʽ

// �������壺
//		L1    	:	��һ��ľ���(�ȣ���-pi ~pi��
//		B1    	:	��һ���γ��(�ȣ���-pi/2 ~pi/2��
//		L2    	:	�ڶ���ľ���(�ȣ���-pi ~pi��
//		B2 	:	�ڶ����γ��(�ȣ���-pi/2 ~pi/2��
//		s  	:	�������루Km��
//		A12   	:	��һ�㵽�ڶ���ķ�λ�ǣ��ȣ����Զ���Ϊ�㣬��ʱ��Ϊ��[0~360]
//		A21   	:	�ڶ��㵽��һ��ķ�λ�ǣ��ȣ����Զ���Ϊ�㣬��ʱ��Ϊ��[0~360]

///////////////////////////////////////////////////////

// Vincenty���⹫ʽ����֪��һ�������꣨B1��L1���͵��ڶ���Ĵ�ؾ���S����λ��A12��
//                   ��ڶ���Ĵ�����꣨B2��L2��������λ��A21��
void  calc_coord_by_distance(double L1, double B1, double s, double A12, double& L2, double& B2, double& A21);

// Vincenty���⹫ʽ����֪����Ĵ�����꣨B1��L1������B2��L2����
//                   ���ؾ���S��������λ��A12��A21��
void  calc_two_pnt_dist(double L1, double B1, double L2, double B2, double& s, double& A12, double& A21);


//�����������������굥λ���룬����������project.cpp
void calc_coord_by_distance_raw(double L1, double B1, double s, double angle, double& L2, double& B2, double& angle21);
void calc_two_pnt_dist_raw(double L1, double B1, double L2, double B2, double& s, double& angle12, double& angle21);


///////////////////////////////////////////////////////

//  Gauss Projection & Inverse-Projection (54/84)

//      B   --  latitude, unit: degree
//      L   --  longitude,  unit: degree
//      L0  --  central meridian,  unit: degree
//      X   --  ordinate,  unit: meter
//      Y   --  abscissa,  unit: meter

///////////////////////////////////////////////////////

int   BLtoXY(double B, double L, double L0, double& X, double& Y);
int   XYtoBL(double X, double Y, double L0, double& B, double& L);


///////////////////////////////////////////////////////

//  Gauss Projection & Inverse-Projection (54)

//      B   --  latitude, unit: degree
//      L   --  longitude,  unit: degree
//      L0  --  central meridian,  unit: degree
//      X   --  ordinate,  unit: meter
//      Y   --  abscissa,  unit: meter

///////////////////////////////////////////////////////

int   BLtoXY_V1(double B, double L, double L0, double& X, double& Y);
int   XYtoBL_V1(double X, double Y, double L0, double& B, double& L);


///////////////////////////////////////////////////////

//  Gauss Projection & Inverse-Projection (54)

//      B   --  latitude, unit: DMS
//      L   --  longitude,  unit: DMS
//      L0  --  central meridian,  unit: degree
//      X   --  ordinate,  unit: meter
//      Y   --  abscissa,  unit: meter

///////////////////////////////////////////////////////

int   BLtoXY_V0(double B, double L, double L0, double* X, double* Y);
int   XYtoBL_V0(double X, double Y, double L0, double* B, double* L);


///////////////////////////////////////////////////////

//  ���ȷ��롱�롰�ȡ���﷽ʽ����

//      DMS  --  ��+��/100+��/10000
//      DD   --  ˫������, ��λ: ��

///////////////////////////////////////////////////////

double  DMStoDD(double dms);
double  DDtoDMS(double dd);


///////////////////////////////////////////////////////

//  ��˹����6�ȷִ�����ͶӰ�����������߾���

//      L  --  ���������߾���, ��λ: ��
//      N  --  ��˹����6�ȷִ���

//      l  --  ���⾭��ֵ, ��λ: ��

///////////////////////////////////////////////////////

int     LtoN(double L);
double  NtoL(int N);
double  LOCM(double l);


///////////////////////////////////////////////////////

//  ��ͼ�������½ǵ���������, ��ʽ: (E|W)xxx(N|S)xx

//      B  --  ���½ǵ�γ��, ��λ: ��
//      L  --  ���½ǵ㾭��, ��λ: ��

//  ע��: ���ص��ַ���Ϊ��̬�洢, ��Ҫ�ͷ�

///////////////////////////////////////////////////////

char  *LB_tile_name(int L, int B);


/**************************�������*************************************/
/* ��  GaussBL  ��˹ͶӰ����
������������������������������������������������������������������������������
��  ��: void GaussBL(double a,double e2,double x,double y,double L0,
double &B,double &L)
��  ��: a   ���򳤰���;
e2  �����һƫ���ʵ�ƽ��;
x y ��˹ƽ��ֱ������, yΪ��Ȼֵ;
L0  ͶӰ�����������߾���;
B L ������(��λ: ����, ����: L�Զ���Ϊ׼);
����ֵ: ��.                                                              */
void GaussBL(double a,double e2,double x,double y,double L0,
			 double &B,double &L);
/***********************************************************************/

#endif
