#ifndef  _Project_Include
#define  _Project_Include

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>

///////////////////////////////////////////////////////

// 选择坐标系统(54 — 表示54坐标系统；84 — 表示84坐标系统)

///////////////////////////////////////////////////////

extern int  CoordSysParam;


///////////////////////////////////////////////////////

// Vincenty正解、反解计算公式

// 参数定义：
//		L1    	:	第一点的经度(度）（-pi ~pi）
//		B1    	:	第一点的纬度(度）（-pi/2 ~pi/2）
//		L2    	:	第二点的经度(度）（-pi ~pi）
//		B2 	:	第二点的纬度(度）（-pi/2 ~pi/2）
//		s  	:	两点间距离（Km）
//		A12   	:	第一点到第二点的方位角（度），以东向为零，逆时针为正[0~360]
//		A21   	:	第二点到第一点的方位角（度），以东向为零，逆时针为正[0~360]

///////////////////////////////////////////////////////

// Vincenty正解公式，已知第一点大地坐标（B1，L1）和到第二点的大地距离S、方位角A12，
//                   求第二点的大地坐标（B2，L2）及反方位角A21。
void  calc_coord_by_distance(double L1, double B1, double s, double A12, double& L2, double& B2, double& A21);

// Vincenty反解公式，已知两点的大地坐标（B1，L1）、（B2，L2），
//                   求大地距离S和正反方位角A12、A21。
void  calc_two_pnt_dist(double L1, double B1, double L2, double B2, double& s, double& A12, double& A21);


//以下两函数两点坐标单位是秒，其余参数详见project.cpp
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

//  “度分秒”与“度”表达方式互换

//      DMS  --  度+分/100+秒/10000
//      DD   --  双精度数, 单位: 度

///////////////////////////////////////////////////////

double  DMStoDD(double dms);
double  DDtoDMS(double dd);


///////////////////////////////////////////////////////

//  高斯坐标6度分带号与投影的中央子午线经度

//      L  --  中央子午线经度, 单位: 度
//      N  --  高斯坐标6度分带号

//      l  --  任意经度值, 单位: 度

///////////////////////////////////////////////////////

int     LtoN(double L);
double  NtoL(int N);
double  LOCM(double l);


///////////////////////////////////////////////////////

//  按图幅的左下角点坐标命名, 格式: (E|W)xxx(N|S)xx

//      B  --  左下角点纬度, 单位: 度
//      L  --  左下角点经度, 单位: 度

//  注意: 返回的字符串为静态存储, 不要释放

///////////////////////////////////////////////////////

char  *LB_tile_name(int L, int B);


/**************************测试添加*************************************/
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
			 double &B,double &L);
/***********************************************************************/

#endif
