#include "StdAfx.h"
#include "MpGlobalfunc.h"
#include <math.h>
#include <direct.h>
#include "project_f.h"

MpGlobalFunc::MpGlobalFunc(void)
{
}

MpGlobalFunc::~MpGlobalFunc(void)
{
}

/// 线性插值计算
///
/// 已知(a0,b0)、(a1,b1)和a，求与a对应的值
double MpGlobalFunc::linear_interpolation(double a0,double b0,double a1,double b1,double a)
{
	double value;

	if( fabs(a0-a1) < 1E-7 )
		value = (b0 + b1) / 2.0;
	else
		value = b1 + (b0 - b1) * (a - a1) / (a0 - a1);

	return value;
}

void ** MpGlobalFunc::NewSpace2D(int Row,int Col,short Size)
{
	int i;
	char **pData;
	pData=new char*[Row];
	if(pData==NULL) 
		return NULL;
	pData[0] = new char[Row*Col*Size];
	memset(pData[0],0,Row*Col*Size);
	if(pData[0]==NULL)
	{
		delete [] pData;
		return NULL;
	}
	for(i=1;i<Row;i++)
		pData[i]=pData[0]+Col*Size*i;
	return (void **)pData;
}

//分配3维内存空间
void *** MpGlobalFunc::NewSpace3D(int Len,int Row,int Col,short Size)
{
	int i,j;
	char ***pData;
	pData=(char ***)NewSpace2D(Len,Row,sizeof(char *));
	if(pData==NULL) return NULL;
	pData[0][0]=new char[Len*Row*Col*Size];
	if(pData[0][0]==NULL)
	{
		MpGlobalFunc::DeleteSpace2D(pData,Len);
		return NULL;
	}
	for(i=0;i<Len;i++)
	for(j=0;j<Row;j++)
		pData[i][j] = pData[0][0]+Row*Col*Size*i+Col*Size*j;
	return (void ***)pData;
}

void ** MpGlobalFunc::NewSpace2D_1(int Row,int Col,short Size)
{
	int i,j;
	char **pData;

	pData=new char*[Row];
	if(pData==NULL) 
		return NULL;

	/// 对每一行分配空间
	for( i=0; i<Row; i++ )
	{
		pData[i] = NULL;
		pData[i] = new char[Col*Size];
		if( pData[i]==NULL )
		{
			for( j=i-1; j>=0; j-- )
			{
				delete [] pData[j];
				pData[j] = NULL;
			}
			break;
		}		
	}

	/// 如果分配失败
	if(pData[0]==NULL)
	{
		delete [] pData;
		return NULL;
	}

	return (void **)pData;
}


//两点画线，给出两端点浮点坐标，左下角点的浮点坐标和分辨率，返回两点间线段的坐标
int MpGlobalFunc::PointToLine( double x1,double y1,double x2,double y2,
				double xo,double yo,double cellsize,
				int &NumPoi,int * (&pX),int * (&pY) )
{
	double pa_a,pa_b,pa_c;
	int i,delta,curX,curY;
	double curx,cury;
	int X1,Y1,X2,Y2;

	if(pX!=NULL||pY!=NULL)
		return -1;

	//求出两个端点的象素坐标
	X1 = int((x1-xo)/cellsize+0.5);
	Y1 = int((y1-yo)/cellsize+0.5);
	X2 = int((x2-xo)/cellsize+0.5);
	Y2 = int((y2-yo)/cellsize+0.5);

	//如果是一个点
	if(X1==X2&&Y1==Y2)
	{
		pX = new int[1];
		pY = new int[1];
		pX[0] = X1;
		pY[0] = Y1;
		NumPoi = 1;
		return 1;
	}
		
	//求出直线方程
	PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c);

	//线扫描
	if( fabs(pa_b-1.0)<qxwLimMinValue )
	{
		NumPoi = abs(X1-X2)+1;
		pX = new int[NumPoi];
		pY = new int[NumPoi];
		pX[0] = X1;
		pY[0] = Y1;
		pX[NumPoi-1] = X2;
		pY[NumPoi-1] = Y2;
		
		delta = (X1<X2) ? 1 : -1 ;
		for( i=1,curX=X1+delta ; i<NumPoi-1 ; i++,curX+=delta )
		{
			pX[i] = curX;
			curx = xo + curX*cellsize;
			cury = -pa_a*curx-pa_c;
			curY = int((cury-yo)/cellsize+0.5);
			pY[i] = curY;
		}
	}
	else
	{
		NumPoi = abs(Y1-Y2)+1;
		pX = new int[NumPoi];
		pY = new int[NumPoi];
		pX[0] = X1;
		pY[0] = Y1;
		pX[NumPoi-1] = X2;
		pY[NumPoi-1] = Y2;

		delta = (Y1<Y2) ? 1 : -1 ;
		for( i=1,curY=Y1+delta ; i<NumPoi-1 ; i++,curY+=delta )
		{
			pY[i] = curY;
			cury = yo + curY*cellsize;
			curx = -pa_b*cury-pa_c;
			curX = int((curx-xo)/cellsize+0.5);
			pX[i] = curX;
		}
	}
	
	return 1;
}

//两点画线，给出两端点浮点坐标，左下角点的浮点坐标和两个方向的分辨率，返回两点间线段的坐标
// 注意，两点画线限定了x和y坐标的范围，如果超出部分按最大或最小值取
// 该函数仍有改进余地，由于限幅的影响，可能出现两个点相同坐标，目前没有排除相同坐标
int MpGlobalFunc::PointToLineWithRgn( double x1,double y1,double x2,double y2,
				double xo,double yo,double csX,double csY, int minX,int maxX,int minY,int maxY,
				int &NumPoi,int * (&pX),int * (&pY) )
{
	int i,RV;

	/// 首先调用不限制范围的两点画线算法
	RV = PointToLine( x1, y1, x2, y2, xo, yo, csX, csY, NumPoi, pX, pY );
	if(RV!=1)
		return RV;

	/// 对输出的x y坐标进行限幅
	for( i=0; i<NumPoi; i++ )
	{
		pX[i] = (pX[i]>maxX)? maxX : pX[i];
		pY[i] = (pY[i]>maxY)? maxY : pY[i];

		pX[i] = (pX[i]<minX)? minX : pX[i];
		pY[i] = (pY[i]<minY)? minY : pY[i];
	}

	return 1;
}

//两点画线，给出两端点浮点坐标，左下角点的浮点坐标和两个方向的分辨率，返回两点间线段的坐标
int MpGlobalFunc::PointToLine( double x1,double y1,double x2,double y2,
				double xo,double yo,double csX,double csY,
				int &NumPoi,int * (&pX),int * (&pY) )
{
	double pa_a,pa_b,pa_c;
	int i,delta,curX,curY;
	double curx,cury;
	int X1,Y1,X2,Y2;

	if(pX!=NULL||pY!=NULL)
		return -1;

	//求出两个端点的象素坐标
	X1 = (int)floor((x1-xo)/csX+0.5);
	Y1 = (int)floor((y1-yo)/csY+0.5);
	X2 = (int)floor((x2-xo)/csX+0.5);
	Y2 = (int)floor((y2-yo)/csY+0.5);

	//如果是一个点
	if(X1==X2&&Y1==Y2)
	{
		pX = new int[1];
		pY = new int[1];
		pX[0] = X1;
		pY[0] = Y1;
		NumPoi = 1;
		return 1;
	}
		
	//求出直线方程
	PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c);

	//线扫描
	if( fabs(pa_b-1.0)<qxwLimMinValue )
	{
		NumPoi = abs(X1-X2)+1;
		pX = new int[NumPoi];
		pY = new int[NumPoi];
		pX[0] = X1;
		pY[0] = Y1;
		pX[NumPoi-1] = X2;
		pY[NumPoi-1] = Y2;
		
		delta = (X1<X2) ? 1 : -1 ;
		for( i=1,curX=X1+delta ; i<NumPoi-1 ; i++,curX+=delta )
		{
			pX[i] = curX;
			curx = xo + curX*csX;
			cury = -pa_a*curx-pa_c;
			curY = (int)floor((cury-yo)/csY+0.5);
			pY[i] = curY;
		}
	}
	else
	{
		NumPoi = abs(Y1-Y2)+1;
		pX = new int[NumPoi];
		pY = new int[NumPoi];
		pX[0] = X1;
		pY[0] = Y1;
		pX[NumPoi-1] = X2;
		pY[NumPoi-1] = Y2;

		delta = (Y1<Y2) ? 1 : -1 ;
		for( i=1,curY=Y1+delta ; i<NumPoi-1 ; i++,curY+=delta )
		{
			pY[i] = curY;
			cury = yo + curY*csY;
			curx = -pa_b*cury-pa_c;
			curX = (int)floor((curx-xo)/csX+0.5);
			pX[i] = curX;
		}
	}
	
	return 1;
}



//输入两个角度，求角度差，单位度
int MpGlobalFunc::DiffOfTwoDir(double dir1,double dir2,bool bDu,double &res)
{
	double dif1,dif2,dif3;
	double PI2;

	if(bDu)
		PI2 = 360;
	else
		PI2 = 2*qxwPI;

	dir1 -= floor(dir1/PI2)*PI2;
	dir2 -= floor(dir2/PI2)*PI2;
	dif1 = fabs(dir1-dir2);
	dif2 = fabs(dir1-dir2+PI2);
	dif3 = fabs(dir1-dir2-PI2);
	res = (dif1<dif2)?dif1:dif2;
	res = (res<dif3)?res:dif3;
	return 1;
}

/******************************************************
函数名：InterLinear1D
描述：
    一维线性插值函数。给出两点坐标（x1,v1）、（x2,v2），
线性插值得X=x时Y的值v。注意x需要在x1与x2之间。

调用本函数的函数清单：
	
输入参数：
	double x1,  点1的X坐标
	double v1,	点1的Y坐标
	double x2,  点2的X坐标
	double v2,  点2的Y坐标
	double x, 	待插值的X    
输出参数：
	double & v, 插值结果

返回值：
	1，成功；
	-1，失败。
********************************************************/
int MpGlobalFunc::InterLinear1D( double x1, double v1 ,double x2, double v2, double x, double & v )
{
	double minX,maxX;
	
	if(x1<x2)
	{
		minX = x1;
		maxX = x2;
	}
	else
	{
		minX = x2;
		maxX = x1;
	}

	//如果x1与x2太接近，认为是一个值，则直接将v1作为结果返回
	if( fabs(x1-x2) < qxwLimMinValue )
	{
		v = v1;
		return 1;
	}

	//如果x不在x1与x2之间，错误返回
	//if( x<minX || x>maxX )
	//	return -1;
		
	//线性插值
	v = v1 + (v2-v1) / (x2-x1) * (x-x1) ;
	return 1;
}

int MpGlobalFunc::InterLinear2D( double x1, double x2, double y1, double y2, double v11, double v12, 
				  double v21, double v22, double x, double y, double &v )
{
	double v1,v2;
	int RV;
	
	//首先对x进行插值
	RV = InterLinear1D(x1,v11,x2,v21, x, v1);
	if(RV==-1)	return -1;

	RV = InterLinear1D(x1,v12,x2,v22, x, v2);
	if(RV==-1)	return -1;

	//再对y进行查值
	RV = InterLinear1D(y1,v1,y2,v2, y, v);
	if(RV==-1)	return -1;

	return 1;
}


/************************************************************************
	函数名称:	PointToLineEquation   
	功能描述:	给出两点坐标，求直线方程。
	算法思想:	
	参数列表:   
				x1，线段端点1的X坐标
				y1，线段端点1的Y坐标
				x2，线段端点2的X坐标
				y2，线段端点2的Y坐标
				para_a，直线方程系数
				para_b，直线方程系数
				para_c，直线方程系数
				minV，精度
				
	返回结果：  
				-1，两点太近，不进行计算
				1，成功
************************************************************************/
int MpGlobalFunc::PointToLineEquation(double x1,double y1,double x2,double y2,double &para_a,
						double &para_b,double &para_c,double minV)
{
	//如果两点太近，不进行计算直接返回
	if( sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) ) < minV )
		return ERR_TOONEAR;
		
	//直线接近水平
	if( fabs(x1-x2)>fabs(y1-y2) )
	{
		para_a = -(y1-y2)/(x1-x2);
		para_c = -y1-para_a*x1;
		para_b = 1.0;
	}

	//直线接近垂直
	else
	{
		para_b = -(x1-x2)/(y1-y2);
		para_c = -x1-para_b*y1;
		para_a = 1.0;
	}
	
	return 1;
}

/************************************************************************
	函数名称:	Cal_TwoLine_Cross   
	功能描述:	计算两直线交点坐标。
	算法思想:	
	参数列表:   
				a1，直线1标准方程参数
				b1，直线1标准方程参数
				c1，直线1标准方程参数
				a2，直线2标准方程参数
				b2，直线2标准方程参数
				c2，直线2标准方程参数
				x，交点X坐标
				y，交点Y坐标
				
	返回值果：  
				-1，无法计算
				1，成功
************************************************************************/
int MpGlobalFunc::Cal_TwoLine_Cross( double a1, double b1, double c1, double a2, 
					  double b2, double c2, double &x, double &y )
{
	double Up,Down;
	
	Up = (b1*c2-b2*c1);
	Down = (a1*b2-a2*b1);
	if( fabs(Up)>fabs(Down) && fabs(Down)/fabs(Up)<1e-300 )	//判断是否可能溢出
		return ERR_CANNOTPROC;
	x = Up/Down;	

	Up = (c2*a1-c1*a2);
	Down = (b1*a2-b2*a1);
	if( fabs(Up)>fabs(Down) && fabs(Down)/fabs(Up)<1e-300 )	//判断是否可能溢出
		return ERR_CANNOTPROC;
	y = Up/Down;
	
	return 1;
}

int MpGlobalFunc::Cal_TwoLine_Cross( double x1, double y1, double x2, double y2, 
					  double x3, double y3, double x4, double y4, 
					  double &x, double &y, double minV )
{
	double a1,b1,c1,a2,b2,c2;
	int RV;

	RV = PointToLineEquation( x1, y1, x2, y2, a1, b1, c1, minV );
	if(RV!=1)
		return RV;

	RV = PointToLineEquation( x3, y3, x4, y4, a2, b2, c2, minV );
	if(RV!=1)
		return RV;

	RV = Cal_TwoLine_Cross( a1, b1, c1, a2, b2, c2, x, y );
	if(RV!=1)
		return RV;	
	
	return 1;
}

/************************************************************************
	函数名称:	Cal_PoiLine_Relation   
	功能描述:	判断点与线段的位置关系。
	算法思想:	
	参数列表:   
				poiX，点的X坐标
				poiY，点的Y坐标
				x1，线段端点1的X坐标
				y1，线段端点1的Y坐标
				x2，线段端点2的X坐标
				y2，线段端点2的Y坐标
				minV，计算精度，可不进行输入，默认为0.00001
				
	返回结果：  
				-1，线段两端点太近，不进行计算
				1，点与线段端点重合
				2，点为线段内的点
				3，点在线段所在的直线上，但在线段外，
				4，点不在线段所在的直线上
************************************************************************/
int MpGlobalFunc::Cal_PoiLine_Relation( double poiX, double poiY, double x1, double y1, 
					   double x2, double y2, double minV )
{
	double pa_a,pa_b,pa_c;
	double dis;
	int bOnLine,bPort,bInLine;

	//如果两点太近，不进行计算直接返回
	if( sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) ) < minV )
		return ERR_TOONEAR;

	bOnLine = bPort = bInLine = 0;		//初始化标志符
	PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c,minV);	//求直线方程
	
	//利用点与直线的距离判断点是否在直线上
	dis = fabs(pa_a*poiX+pa_b*poiY+pa_c) / sqrt(pa_a*pa_a+pa_b*pa_b) ;
	if( dis<minV )
		bOnLine = 1 ;

	//判断点是否是线段端点
	if( 
		( sqrt( (poiX-x1)*(poiX-x1)+(poiY-y1)*(poiY-y1) ) < minV ) ||
		( sqrt( (poiX-x2)*(poiX-x2)+(poiY-y2)*(poiY-y2) ) < minV )
		)
		bPort = 1;
	
	//判断点是否在线段两端点坐标范围内
	if( fabs(x1-x2) >= fabs(y1-y2) )	//斜率小于1，横着的线段
	{
		if( poiX>min(x1,x2) && poiX<max(x1,x2) )
			bInLine = 1;
	}
	else							//斜率小于1，竖着的线段
	{
		if( poiY>min(y1,y2) && poiY<max(y1,y2) )
			bInLine = 1;
	}

	if(bPort==1)
		return 1;
	else if( bOnLine==1 && bInLine==1 )
		return 2;
	else if( bOnLine==1 )
		return 3;
	else
		return 4;	
}

/************************************************************************
	函数名称:	Cal_TwoLine_Relation   
	功能描述:	判断线段与线段的位置关系。
	算法思想:	
	参数列表:   
				x1，线段1端点1的X坐标
				y1，线段1端点1的Y坐标
				x2，线段1端点2的X坐标
				y2，线段1端点2的Y坐标
				x3，线段2端点1的X坐标
				y3，线段2端点1的Y坐标
				x4，线段2端点2的X坐标
				y4，线段2端点2的Y坐标
				Type，线段关系类型
				Lin1Type，交点与线段1的关系
				Lin2Type，交点与线段2的关系
				minV，计算精度，可不进行输入，默认为0.00001
				
	结果：
			Type =	0，线段与线段所在直线重合，Line1Type与Line2Type无效
					1，线段与线段所在直线平行，Line1Type与Line2Type无效
					2，线段与线段所在直线相交
			Line1Type =	0，交点在线段1外
						1，交点为线段1端点
						2，交点在线段1内
			Line2Type =	0，交点在线段2外
						1，交点为线段2端点
						2，交点在线段2内
	  
	返回值：  
				-1，点距离太近，不进行计算
				1，成功
************************************************************************/
int MpGlobalFunc::Cal_TwoLine_Relation( double x1, double y1, double x2, double y2,
						 double x3, double y3, double x4, double y4, 
						 int &Type, int &Line1Type, int &Line2Type, double minV )
{
	double a1,b1,c1,a2,b2,c2;
	double dis1,dis2,dis3,croX,croY;
	double minX,maxX,minY,maxY;
	double x0,y0;
	int RV;

	//分别判断两条线段的两个端点是否太近，太近则函数直接返回
	dis1 = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
	dis2 = sqrt( (x3-x4)*(x3-x4) + (y3-y4)*(y3-y4) );
	if( dis1<=minV || dis2<=minV )
		return ERR_TOONEAR;

	//计算两线段直线方程
	PointToLineEquation(x1,y1,x2,y2,a1,b1,c1,minV);
	PointToLineEquation(x3,y3,x4,y4,a2,b2,c2,minV);

	//取线段1的中点
	x0 = (x1+x2)/2;
	y0 = (y1+y2)/2;

	//计算线段1三个点到线段2所在直线的距离
	dis1 = fabs(a2*x1+b2*y1+c2) / sqrt(a2*a2+b2*b2) ;
	dis2 = fabs(a2*x2+b2*y2+c2) / sqrt(a2*a2+b2*b2) ;
	dis3 = fabs(a2*x0+b2*y0+c2) / sqrt(a2*a2+b2*b2) ;

	//判断射线与直线的关系
	if( fabs(dis1)<=minV && fabs(dis2)<=minV && fabs(dis3)<=minV )
		Type = 0;
	else if( fabs(dis1-dis2)<=minV && fabs(dis2-dis3)<=minV )
		Type = 1;
	else
		Type = 2;

	//如果重合或者平行，直接返回，不再进行判断
	if(Type!=2)
		return 1;

	RV = Cal_TwoLine_Cross(a1,b1,c1,a2,b2,c2,croX,croY);	//计算两条直线交点坐标

	//判线段1与交点的情况
	maxX = max(x1,x2);
	minX = min(x1,x2);
	maxY = max(y1,y2);
	minY = min(y1,y2);
	dis1 = sqrt( (croX-x1)*(croX-x1) + (croY-y1)*(croY-y1) );
	dis2 = sqrt( (croX-x2)*(croX-x2) + (croY-y2)*(croY-y2) );
	if( dis1<=minV || dis2<=minV )
		Line1Type = 1;
	else if( 
		( fabs(x1-x2)>=fabs(y1-y2) && croX>minX && croX<maxX ) ||
		( fabs(x1-x2)<fabs(y1-y2) && croY>minY && croY<maxY )
		)
		Line1Type = 2;
	else
		Line1Type = 0;

	//判线段2与交点的情况
	maxX = max(x3,x4);
	minX = min(x3,x4);
	maxY = max(y3,y4);
	minY = min(y3,y4);
	dis1 = sqrt( (croX-x3)*(croX-x3) + (croY-y3)*(croY-y3) );
	dis2 = sqrt( (croX-x4)*(croX-x4) + (croY-y4)*(croY-y4) );
	if( dis1<=minV || dis2<=minV )
		Line2Type = 1;
	else if( 
		( fabs(x3-x4)>=fabs(y3-y4) && croX>minX && croX<maxX ) ||
		( fabs(x3-x4)<fabs(y3-y4) && croY>minY && croY<maxY )
		)
		Line2Type = 2;
	else
		Line2Type = 0;

	return 1;
}


/************************************************************************
	函数名称:	Cal_RadialLine_Relation   
	功能描述:	判断射线与线段的关系。
	算法思想:	
	参数列表:   
				r_x，射线端点X坐标
				r_y，射线端点Y坐标
				r_dir，射线方向，x正向为零，逆时针正，单位度
				x1，线段端点1的X坐标
				y1，线段端点1的Y坐标
				x2，线段端点2的X坐标
				y2，线段端点2的Y坐标
				Type，射线与线段的关系类型
				RadType，交点与射线的关系
				LineType，交点与线段的关系
				minV，计算精度，可不进行输入，默认为0.00001
				
	结果：
			type =	0，射线与线段所在直线重合，RadType与LineeType无效
					1，射线与线段所在直线平行，RadType与LineeType无效
					2，射线与线段所在直线所在直线相交
			RadType =	0，交点不在射线上
						1，交点为射线端点
						2，交点在射线上
			LineType =	0，交点在线段外
						1，交点为线段端点
						2，交点在线段内		
	  
	返回值：  
		-1，点距离太近，不进行计算
		1，成功
************************************************************************/
int MpGlobalFunc::Cal_RadialLine_Relation( double r_x, double r_y, double r_dir, 
							double x1, double y1, double x2, double y2, 
							int &Type, int &RadType,int &LineType, double minV)
{
	double rx1,ry1,rx2,ry2;
	double dis1,dis2,dis3,pa_a,pa_b,pa_c;
	double r_a,r_b,r_c,maxX,minX,maxY,minY;
	double crossX,crossY;
	double ang;
	double sinA,cosA;
	double len=1000;
	int RV;

	ang = r_dir/180*qxwPI;	//转化为弧度

	//在射线上以一定间隔距离（len）取两点
	rx1 = r_x + len*cos( ang );
	ry1 = r_y + len*sin( ang );
	rx2 = rx1 + len*cos( ang );
	ry2 = ry1 + len*sin( ang );

	//计算线段直线方程
	RV = PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c);
	if(RV!=1)
		return RV;

	//计算射线三点与直线的距离
	dis1 = fabs( r_x*pa_a + r_y*pa_b + pa_c ) / sqrt( pa_a*pa_a + pa_b*pa_b );
	dis2 = fabs( rx1*pa_a + ry1*pa_b + pa_c ) / sqrt( pa_a*pa_a + pa_b*pa_b );
	dis3 = fabs( rx2*pa_a + ry2*pa_b + pa_c ) / sqrt( pa_a*pa_a + pa_b*pa_b );

	//判断射线与直线的关系
	if( fabs(dis1)<=minV && fabs(dis2)<=minV && fabs(dis3)<=minV )
		Type = 0;
	else if( fabs(dis1-dis2)<=minV && fabs(dis2-dis3)<=minV )
		Type = 1;
	else
		Type = 2;

	//如果射线与线段所在直线不相交，直接返回
	if(Type!=2)
		return 1;

	//求射线直线方程
	PointToLineEquation(r_x,r_y,rx1,ry1,r_a,r_b,r_c);

	//计算两条直线交点坐标
	//crossX = ( pa_b*r_c - r_b*pa_c ) / ( pa_a*r_b - r_a*pa_b );
	//crossY = ( r_c*pa_a - pa_c*r_a ) / ( pa_b*r_a - r_b*pa_a );
	RV = Cal_TwoLine_Cross(pa_a,pa_b,pa_c,r_a,r_b,r_c,crossX,crossY);
	if(RV!=1)
		return RV;
	
	//判断交点与射线的关系
	dis1 = sqrt( (crossX-r_x)*(crossX-r_x) + (crossY-r_y)*(crossY-r_y) );
	if( dis1<minV )
		RadType = 1 ;
	else
	{
		sinA = (crossY-r_y)/dis1;
		cosA = (crossX-r_x)/dis1;	
		if( sin(ang)*sinA>=0 && cos(ang)*cosA>=0 )
			RadType = 2; 
		else
			RadType = 0;
	}
	
	//判断交点与线段的关系
	maxX = max(x1,x2);
	minX = min(x1,x2);
	maxY = max(y1,y2);
	minY = min(y1,y2);
	dis1 = sqrt( (crossX-x1)*(crossX-x1) + (crossY-y1)*(crossY-y1) );
	dis2 = sqrt( (crossX-x2)*(crossX-x2) + (crossY-y2)*(crossY-y2) );
	if( dis1<minV || dis2<minV )
		LineType = 1;
	else if( 
		( fabs(x1-x2)>=fabs(y1-y2) && crossX>minX && crossX<maxX ) ||
		( fabs(x1-x2)<fabs(y1-y2) && crossY>minY && crossY<maxY )
		)
		LineType = 2;
	else
		LineType = 0;
	
	return 1;
}

/************************************************************************
	函数名称:	Cal_PoiPolygon_Relation   
	功能描述:	判断点与多边形的位置关系。
	算法思想:	
	参数列表:   
				NumPoi，多边形角点个数
				pX，多边形角点X坐标数组指针，顺时针排列
				pY，多边形角点Y坐标数组指针，顺时针排列
				PoiX，点X坐标
				PoiY，点Y坐标				
				minV，计算精度，可不进行输入，默认为0.00001
		  
	返回值：  
				0，内部
				1，外部
				2，点在多边形边上
				-1，计算出错或者本函数未能有效判断出二者关系
************************************************************************/
int MpGlobalFunc::Cal_PoiPolygon_Relation( int NumPoi, double *pX, double *pY, 
						  double PoiX, double PoiY, double minV )
{
	int CroDotNumL,CroDotNumR;	//左右交点个数
	double Dir0=0;				//起始方向，0度
	double DirDlt=20;			//方向步长，20度
	double curDir;
	int i;
	int EviNumLim=3;			//证据数目门限
	double x1,y1,x2,y2;			//线段两端点
	double x0,y0,dir;			//射线起点和方向
	int RV;
	int Type,RadType,LineType;	//射线与线段关系函数返回值
	int evidenceOut,evidenceIn;	//证据变量
	double maxX,minX,maxY,minY,maxDis,minDis,dis,LenAxis;
	double minLen;
	double angView;		//点对多边形的视场角
	int NumLoop;		//总循环次数

	//首先判断点是否在多边形的边上
	//并求多边形最短边长
	minLen = sqrt( (pX[0]-pX[1])*(pX[0]-pX[1]) + (pY[0]-pY[1])*(pY[0]-pY[1]) );
	for( i=0; i<NumPoi; i++ )
	{
		//确定线段两端点
		x1 = pX[i];
		y1 = pY[i];
		if(i==NumPoi-1)
		{
			x2 = pX[0];
			y2 = pY[0];
		}
		else
		{
			x2 = pX[i+1];
			y2 = pY[i+1];
		}

		dis = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
		minLen = min(minLen,dis);

		//判断点是否在当前的边上，如果在，直接返回
		RV = Cal_PoiLine_Relation( PoiX, PoiY, x1, y1, x2, y2, minV );
		if( RV==1 || RV==2 )
			return 2;
	}

	//求多边形外接矩形范围
	//并且求点到多边形顶点的最远和最近距离
	maxX = minX = pX[0];
	maxY = minY = pY[0];
	maxDis = minDis = sqrt( (PoiX-pX[0])*(PoiX-pX[0]) + (PoiY-pY[0])*(PoiY-pY[0]) );
	for( i=1; i<NumPoi; i++ )
	{
		maxX = max(maxX,pX[i]);
		minX = min(minX,pX[i]);
		maxY = max(maxY,pY[i]);
		minY = min(minY,pY[i]);
		dis = sqrt( (PoiX-pX[i])*(PoiX-pX[i]) + (PoiY-pY[i])*(PoiY-pY[i]) );
		minDis = min(minDis,dis);
		maxDis = max(maxDis,dis);
	}

	//求外接矩形对角线长，即此多边形覆盖区域的轴长
	maxX = maxX-minX;
	maxY = maxY-minY;
	LenAxis = sqrt(maxX*maxX+maxY*maxY);

	//如果点到多边形顶点的最近距离比轴长还长，说明点在多边形外
	if( minDis>LenAxis+minV )
		return 1;	

	//求点对多边形的最小视场角
	angView = atan2( minLen/2, maxDis ) *2 ;
	angView = angView/qxwPI*180;

	//点不在多边形的边上，通过射线与边相交法判断二者关系
	//对角度作循环，从0度开始，步长为30度渐增，每个角度判断射线与边的相交情况
	//如果存在两个角度证实点在多边形内部，则认为是内部点，函数返回
	//如果存在两个角度证实点在多边形外部，则认为是外部点，函数返回
	evidenceOut = evidenceIn = 0;
	NumLoop = 1;
	for( curDir=Dir0; ; curDir+=DirDlt )
	{
		//如果总循环次数大于5，退出
		if(NumLoop>5)
			break;
		
		//如果大于等于Dir0+180
		if( curDir>=Dir0+180 )
		{
			Dir0 += angView/5;
			curDir = Dir0;
			NumLoop ++;
		}
		
		//对每个边进行循环计算，计算左、右射线与边的交点
		CroDotNumR = CroDotNumL = 0;
		for( i=0; i<NumPoi; i++ )
		{
			//确定线段两端点
			x1 = pX[i];
			y1 = pY[i];
			if(i==NumPoi-1)
			{
				x2 = pX[0];
				y2 = pY[0];
			}
			else
			{
				x2 = pX[i+1];
				y2 = pY[i+1];
			}

			//确定右射线起点和方向
			x0 = PoiX;
			y0 = PoiY;
			dir = curDir/180*qxwPI;

			//判断射线和线段的相交情况
			Cal_RadialLine_Relation( x0, y0, curDir, x1, y1, x2, y2, Type, RadType, LineType );
			if( Type==2 && LineType==2 && RadType==2 )
				CroDotNumR++;

			//确定左射线起点和方向
			x0 = PoiX;
			y0 = PoiY;
			dir = (curDir+180)/180*qxwPI;

			//判断射线和线段的相交情况
			Cal_RadialLine_Relation( x0, y0, curDir+180, x1, y1, x2, y2, Type, RadType, LineType );
			if( Type==2 && LineType==2 && RadType==2 )
				CroDotNumL++;				
		}

		//如果左右交点个数和相加为偶数，此方向有效
		if( (CroDotNumL+CroDotNumR)%2==0 )	//&& (CroDotNumL+CroDotNumR)>0
		{
			if( CroDotNumL%2==0 && CroDotNumR%2==0 )
				evidenceOut ++;
			else
				evidenceIn ++;
		}

		//如果内部或外部的证据次数大于等于3，退出循环
		if( evidenceOut>=EviNumLim || evidenceIn>=EviNumLim )
			break;
	}

	if(evidenceIn>=EviNumLim)
		return 0;
	else if(evidenceOut>=EviNumLim)
		return 1;
	else
		return ERR_CANNOTPROC;
}


/************************************************************************
	函数名称:	Cal_LinePolygon_Relation   
	功能描述:	判断线段与多边形的位置关系。
	算法思想:	
	参数列表:   
				NP，多边形角点个数
				pX，多边形角点X坐标数组指针，顺时针排列
				pY，多边形角点Y坐标数组指针，顺时针排列
				X1，线段端点1X坐标
				Y1，线段端点1Y坐标
				X2，线段端点2X坐标
				Y2，线段端点2Y坐标
				minV，计算精度，可不进行输入，默认为0.00001
		  
	返回值：  
			-1，计算出错
			0，线段在多边形内
			1，线段在多边形外
			2，线段与多边形相交
			3，线段在多边形的边上
			4，点相交
************************************************************************/
int MpGlobalFunc::Cal_LinePolygon_Relation( int NP, double *pX, double *pY, double X1, double Y1, 
							 double X2, double Y2, double minV )
{
	double crox,croy;
	int i,j,num;
	int evIn;	//点在多边形内部的个数
	int evOut;	//点在多边形外部的个数
	int evOn;	//点在多边形边上的个数
	int RV,RV1,RV2;
	double x1,x2,y1,y2;
	int Type,Line1Type,Line2Type;
	double *px,*py;
	
	px = py = NULL;

	//判断线段两个端点是否位于多边形的一条边内
	//如果在，说明线段在多边形边上
	for( i=0; i<NP; i++ )
	{
		x1 = pX[i];
		y1 = pY[i];
		if(i==NP-1)
		{
			x2 = pX[0];
			y2 = pY[0];
		}else{
			x2 = pX[i+1];
			y2 = pY[i+1];
		}

		//两端点1与边的关系
		RV1 = Cal_PoiLine_Relation( X1, Y1, x1, y1, x2, y2, minV );
		RV2 = Cal_PoiLine_Relation( X2, Y2, x1, y1, x2, y2, minV );
		if( (RV1==1||RV1==2) && (RV2==1||RV2==2) )
			return 3;
	}

	//判断线段两个端点与多边形关系
	RV1 = Cal_PoiPolygon_Relation( NP, pX, pY, X1, Y1, minV );
	RV2 = Cal_PoiPolygon_Relation( NP, pX, pY, X2, Y2, minV );
	if( RV1 ==-1 ) 
		return RV1;
	if( RV2 ==-1 )
		return RV2;

	//如果一个在内部一个在外部说明相交
	if( ( RV1==0 && RV2==1 ) || ( RV1==1 && RV2==0 ) )
		return 2;

	px = new double[NP*2];
	py = new double[NP*2];
	num=0;

	if( px==NULL || py==NULL )
	{
		if(px)
			delete [] px;
		if(py)
			delete [] py;
		px = py = NULL;
		return ERR_NEWSPACE;
	}

	for( i=0; i<NP; i++ )
	{
		x1 = pX[i];
		y1 = pY[i];
		
		if(i==NP-1)
		{
			x2 = pX[0];
			y2 = pY[0];
		}else{
			x2 = pX[i+1];
			y2 = pY[i+1];
		}

		//当前边与线段的关系
		RV = Cal_TwoLine_Relation( x1, y1, x2, y2, X1, Y1, X2, Y2, 
			Type, Line1Type, Line2Type, minV );

		//如果线段和边相交，交点加入交点数组中
		if( 
			Type==2 && 
			( Line1Type==1 || Line1Type==2 ) && 
			( Line2Type==1 || Line2Type==2 )
			)
		{
			Cal_TwoLine_Cross( x1, y1, x2, y2, X1, Y1, X2, Y2, crox, croy, minV );
			px[num] = crox;
			py[num] = croy;
			num++;
		}		
	}

	//将线段端点也加入数组
	px[num] = X1;
	py[num] = Y1;
	num++;
	px[num] = X2;
	py[num] = Y2;
	num++;
	
	//对交点数组进行排序
	if( fabs(X1-X2)>fabs(Y1-Y2) ) //横向，按X排序
	{
		for( i=num-1; i>0; i-- )
		{
			for(j=0;j<i;j++)
			{
				if(px[j]>px[i])
				{
					x1 = px[j];
					px[j] = px[i];
					px[i] = x1;
					y1 = py[j];
					py[j] = py[i];
					py[i] = y1;
				}
			}
		}

	}
	else						//纵向，按Y排序
	{
		for( i=num-1; i>0; i-- )
		{
			for(j=0;j<i;j++)
			{
				if(py[j]>py[i])
				{
					x1 = px[j];
					px[j] = px[i];
					px[i] = x1;
					y1 = py[j];
					py[j] = py[i];
					py[i] = y1;
				}
			}
		}
	}

	//对两个相邻交点的中点进行判断
	evIn = evOut = evOn = 0;
	for ( i=0; i<num-1; i++ )
	{
		x1 = (px[i]+px[i+1]) /2 ;
		y1 = (py[i]+py[i+1]) /2 ;

		RV = Cal_PoiPolygon_Relation( NP, pX, pY, x1, y1, minV );
		if(RV==0)
			evIn++;
		if(RV==1)
			evOut++;
		if(RV==2)
			evOn++;
	}

	delete [] px;
	delete [] py;
	px = py = NULL;

	//只有外部，没有内部也没有边上的，外部
	if( evIn==0 && evOn==0 && evOut>0 )	
		return 1;

	//有外部，没有内部，有边上的，点相交
	if( evIn==0 && evOn>0 && evOut>0 )	
		return 4;

	//只有内部，内部
	if( evIn>=0 && evOn==0 && evOut==0 )
		return 0;

	//有内部，有边上，没有外部，内部
	if( evIn>=0 && evOn>0 && evOut==0 )
		return 0;

	//有内部，无边上，有外部，相交
	if( evIn>0 && evOn==0 && evOut>0 )
		return 2;

	//有内部，有边上，有外部，相交
	if( evIn>0 && evOn>0 && evOut>0 )
		return 2;

	return ERR_CANNOTPROC;
}



/************************************************************************
	函数名称:	Cal_TwoPolygon_Relation   
	功能描述:	判断多边形1与相对于多边形2的位置关系。
	算法思想:	
	参数列表:   
				NP1，多边形1角点个数
				pX1，多边形1角点X坐标数组指针，顺时针排列
				pY1，多边形1角点Y坐标数组指针，顺时针排列
				NP2，多边形2角点个数
				pX2，多边形2角点X坐标数组指针，顺时针排列
				pY2，多边形2角点Y坐标数组指针，顺时针排列
				minV，计算精度，可不进行输入，默认为0.00001
		  
	返回值：  
			-1，计算出错
			0，1在2内部
			1，1与2相交
			2，1与2相邻，有边重合或有顶点在另外多边形的边上
			3，1与2相离，互为外部
			4，2在1内部
			5，重合
************************************************************************/
int MpGlobalFunc::Cal_TwoPolygon_Relation( int NP1, double *pX1, double *pY1,
							int NP2, double *pX2, double *pY2, double minV )
{
	int i;
	int RV;
	double x1,x2,y1,y2;
	int Num10;	//多边形1的边在多边形2内部的个数
	int Num11;	//多边形1的边在多边形2外部的个数
	int Num12;	//多边形1的边与多边形2相交的个数
	int Num13;	//多边形1的边在多边形2边上的个数
	int Num14;	//多边形1的与多边形2点相交的个数
	int Num20;	//多边形2的边在多边形1内部的个数
	int Num21;	//多边形2的边在多边形1外部的个数
	int Num22;	//多边形2的边与多边形1相交的个数
	int Num23;	//多边形2的边在多边形1边上的个数
	int Num24;	//多边形2的与多边形1点相交的个数
	
	//先对多边形1的每个边与多边形2关系进行判断
	Num10 = Num11 = Num12 = Num13 = Num14 = 0;
	for( i=0; i<NP1; i++ )
	{
		x1 = pX1[i];
		y1 = pY1[i];
		
		if(i==NP1-1)
		{
			x2 = pX1[0];
			y2 = pY1[0];
		}else{
			x2 = pX1[i+1];
			y2 = pY1[i+1];
		}

		RV = Cal_LinePolygon_Relation( NP2, pX2, pY2, x1, y1, x2, y2, minV );
		if(RV==0)
			Num10++;
		if(RV==1)
			Num11++;
		if(RV==2)
			Num12++;
		if(RV==3)
			Num13++;	
		if(RV==4)
			Num14++;	
		if(RV<0)
			return RV;
	}

	//再对多边形2的每个边与多边形1关系进行判断
	Num20 = Num21 = Num22 = Num23 = Num24 = 0;
	for( i=0; i<NP2; i++ )
	{
		x1 = pX2[i];
		y1 = pY2[i];
		if(i==NP2-1)
		{
			x2 = pX2[0];
			y2 = pY2[0];
		}else{
			x2 = pX2[i+1];
			y2 = pY2[i+1];
		}

		RV = Cal_LinePolygon_Relation( NP1, pX1, pY1, x1, y1, x2, y2, minV );
		if(RV==0)
			Num20++;
		if(RV==1)
			Num21++;
		if(RV==2)
			Num22++;
		if(RV==3)
			Num23++;	
		if(RV==4)
			Num24++;	
		if(RV<0)
			return RV;
	}

	//重合
	if( Num13==NP1 && Num23==NP2 )
		return 5;

	//多边形1的每条边均在多边形2内部或者在多边形2的边上，且必有内部
	//则多边形1在多边形2内部
	if( Num10+Num13==NP1 && Num10>0 )	
		return 0;

	//多边形2的每条边均在多边形1内部或者在多边形1的边上，且必有内部
	//则多边形2在多边形1内部
	if( Num20+Num23==NP2 && Num20>0 )	
		return 4;

	//1、2相交只有两种情况：（1）1的边与2相交 （2）1的边与2点交且有内部
	if( 
		Num12>0 || 
		( Num10>0 && Num14>0 )
		)
		return 1;

	//多边形1的边均在多边形2外部，并且多边形2不在多边形1内部
	//二者相离
	if(Num11==NP1)
		return 3;

	//多边形1的边与多边形2的关系点交、外部或位于边上重合，且不能只是外部
	if( Num11+Num13+Num14==NP1 && Num11!=NP1 )
		return 2;

	return ERR_CANNOTPROC;
}


unsigned int MpGlobalFunc::GetFileLen(const char *pFN)
{
	FILE *pF=NULL;
	unsigned int Len;
	
	pF = fopen(pFN,"rb");
	if(pF==NULL)
		return 0;

	fseek(pF,0,SEEK_END);
	Len = ftell(pF);
	fclose(pF);
	return Len;
}

int MpGlobalFunc::ReadFileToBuf(const char * FN,char *(&pDat),unsigned int &Len )
{
	FILE *pF = NULL;

	Len = 0;
	pDat = NULL;

	pF = fopen( FN, "rb" );
	if(pF==NULL)
		return 1;

	fseek(pF,0,SEEK_END);
	Len = ftell(pF);

	if( Len==0 )
	{
		fclose(pF);
		return 1;
	}

	pDat = new char[Len+1];
	fseek( pF, 0, SEEK_SET );
	fread( pDat, 1, Len, pF );
	fclose(pF);

	pDat[Len] = '\0';
	return 1;
}

int MpGlobalFunc::CopyFileToFile(const char *FNOld, const char *FNNew, int bOverWrite )
{
	char * pDat;
	int bEx,RV;
	FILE *pF;
	unsigned int Len;

	//检测文件是否存在
	pF = fopen(FNNew,"rb");
	if( pF )
	{
		bEx = 1;
		fclose(pF);
	}
	else
		bEx = 0;
	

	//如果文件已经存在且覆盖模式为否，直接返回0
	if( bEx==1 && bOverWrite==0 )
		return 0;

	//读入文件，如果读出0个字节，直接返回0
	RV = ReadFileToBuf( FNOld, pDat, Len );
	//20161223lls add start 304静态分析问题修改
	if(pDat == NULL)
	{
		return 0;
	}
	//20161223lls add end 304静态分析问题修改
	if( Len==0 )
	{
		return 0;
	}

	//写入新文件
	pF = fopen(FNNew,"wb+");
	//20161223lls add start 304静态分析问题修改
	if(pF == NULL)
	{
		return 0;
	}
	//20161223lls add end 304静态分析问题修改
	fwrite( pDat, 1, Len, pF );
	fclose(pF);

	delete [] pDat;
	return 1;
}


int MpGlobalFunc::stringMakeUp(char * pDat)
{
	int i,len;

	if( pDat==NULL )
		return 1;

	len = (int)strlen(pDat);

	for( i=0; i<len; i++ )
	{
		if( pDat[i]>=97 && pDat[i]<=122 )
			pDat[i] = pDat[i]-32;
	}

	return 1;
}

///////////////////////////////////////////////
/// 字符串解析函数							///
/// 编写：祁加强							///
/// 输入参数：								///
/// Src			源字符串					///
/// Splitter	分割字符串列表				///
/// 输出参数：								///
/// vecStrings	解析得到的字符串列表		///
///////////////////////////////////////////////
int MpGlobalFunc::SplitString( CString Src, std::vector<CString> Splitter, std::vector<CString> &vecStrings )
{

	for(UINT i=1; i<Splitter.size(); i++)
	{
		Src.Replace(Splitter[i], Splitter[0]);
	}

	vecStrings.clear();

	Src.TrimLeft(Splitter[0]);
	int pos = Src.FindOneOf(Splitter[0]);
	while ( pos >= 0)
	{
		vecStrings.push_back( Src.Left(pos) );
		Src.TrimLeft(vecStrings[vecStrings.size()-1]);
		Src.TrimLeft(Splitter[0]);
		pos = Src.FindOneOf(Splitter[0]);
	}
	vecStrings.push_back( Src );

	return 0;
}




int MpGlobalFunc::JieCheng(int a)
{
	int b;

	if(a<2)
		return a;
	else
	{
		b = JieCheng(a-1);
		return a*b;
	}
}


int MpGlobalFunc::zuhe(int NumAll,int NumSel,int &NumRes, int **(&ppRes) )
{
	int Hei,Wid;
	int i,j,k;
	int **ppR=NULL;
	int NumR;
	int maxV;
	int SNBeg;
	int RV;

	if(NumSel<1)
		return -1;

	Hei = JieCheng(NumAll);
	Wid = NumAll;
	ppRes = (int **)NewSpace2D(Hei,Wid,sizeof(int));

	for(i=0;i<Hei;i++)
		for(j=0;j<Wid;j++)
			ppRes[i][j]=-1;

	NumRes = 0;

	/// 如果是一个元素的组合，直接输出每个元素
	if(NumSel==1)
	{
		for(i=0;i<NumAll;i++)
		{
			ppRes[NumRes][0] = i;
			NumRes++;
		}
		return 1;
	}

	/// 如果不是一个的组合，递归调用数目减1的组合
	RV = zuhe(NumAll,NumSel-1,NumR, ppR );
	if(RV!=1)
	{
		if (ppR) //20161225 304pingce
		{
			delete []ppR;
		}
		return -1;
	}
		

	/// 对减1组合的所有情况增加1个元素
	for( i=0; i<NumR; i++ )
	{
		/// 当前组合结果的最大值
		maxV = -1;
		for(j=0;j<NumAll;j++)
		{
			if(ppR[i][j]==-1)
				break;
			maxV = max(maxV,ppR[i][j]);
		}
		SNBeg = j;

		/// 对当前NumSel-1元素组合结果追加一个元素
		/// 形成多个NumSel元素组合
		for( j=maxV+1; j<NumAll; j++ )
		{
			/// 将前面NumSel-1个元素复制
			for( k=0; k<SNBeg; k++ )
				ppRes[NumRes][k] = ppR[i][k];

			/// 追加一个元素
			ppRes[NumRes][SNBeg] = j;
			NumRes++;
		}
	}

	MpGlobalFunc::DeleteSpace2D(ppR);
	return 1;		
}

/// 两个规则矩形是否相交
/// 1，相交，0，不相交
int MpGlobalFunc::Chk2RegularRect(double xl1,double yb1,double xr1,double yt1,
						   double xl2,double yb2,double xr2,double yt2)
{
	/// 判断两个规则矩形是否相交，等于判断两个规则矩形是否不相交
	/// 如果不相交，则至少需要满足四个条件之一
	/// 矩形1在矩形2左侧，1在2右侧，1在2上上，1在2下方
	/// 如果四个条件都不成立，说明两个矩形相交

	/// 矩形1在矩形2左侧
	if( xr1<xl2 )
		return 0;

	/// 矩形1在矩形2右侧
	if( xl1>xr2 )
		return 0;

	/// 矩形1在矩形2下方
	if( yt1<yb2 )
		return 0;

	/// 矩形1在矩形2上方
	if( yb1>yt2 )
		return 0;

	/// 以上四种情况都不成立，则两矩形相交
	return 1;

	///// 以下是原来的另一种思路，对相交不相交的情况总结不全，算法错误
	//
	//// 判断第一个外接矩形四个角点是否位于第2个矩形内
	//// 再判断第二个外接矩形四个角点是否位于第1个矩形内
	//// 如果两个都否，说明两外接矩形不相干

	//// 第一个矩形左下角点是否在第二个矩形内
	//if( xl1>=xl2 && xl1<=xr2 && yb1>=yb2 && yb1<=yt2 )
	//	return 1;

	//// 第一个矩形左上角点是否在第二个矩形内
	//if( xl1>=xl2 && xl1<=xr2 && yt1>=yb2 && yt1<=yt2 )
	//	return 1;

	//// 第一个矩形右上角点是否在第二个矩形内
	//if( xr1>=xl2 && xr1<=xr2 && yt1>=yb2 && yt1<=yt2 )
	//	return 1;

	//// 第一个矩形右下角点是否在第二个矩形内
	//if( xr1>=xl2 && xr1<=xr2 && yb1>=yb2 && yb1<=yt2 )
	//	return 1;

	//// 第二个矩形左下角点是否在第一个矩形内
	//if( xl2>=xl1 && xl2<=xr1 && yb2>=yb1 && yb2<=yt1 )
	//	return 1;

	//// 第二个矩形左上角点是否在第一个矩形内
	//if( xl2>=xl1 && xl2<=xr1 && yt2>=yb1 && yt2<=yt1 )
	//	return 1;

	//// 第二个矩形右上角点是否在第一个矩形内
	//if( xr2>=xl1 && xr2<=xr1 && yt2>=yb1 && yt2<=yt1 )
	//	return 1;

	//// 第二个矩形右下角点是否在第一个矩形内
	//if( xr2>=xl1 && xr2<=xr1 && yb2>=yb1 && yb2<=yt1 )
	//	return 1;

	//return 0;
}

int MpGlobalFunc::SuperCopyDir(const char *DirS,const char *DirO,int NumTry,int InterTime)
{
	CString Dir1,Dir2,curName;
	CString D1,D2;
	CFileFind myFind;
	int i;
	int RV;
	int bOK;
	BOOL bSuc;
	char BakCurPath[500];
	
	/// 检查输入的源目录和目标目录字符串是否正常
	if( DirS==NULL || strlen(DirS)==0 )
		return -1;

	if( DirO==NULL || strlen(DirO)==0 )
		return -1;

	if( NumTry<1 || InterTime<1 )
		return -1;
	
	/// 检查目标文件夹是否存在，如果存在，先删除
	Dir2 = DirO;
	Dir2.TrimRight("\\");
	if( myFind.FindFile(Dir2) )
	{
		/// 连续删除NumTry次，成功直接跳出
		for( bOK=0,i=0; i<NumTry; i++ )
		{
			RV = RMDirAndFiles(Dir2);
			if(RV==1)
			{
				bOK=1;
				break;
			}

			/// 两次重试间隔InterTime毫秒
			Sleep(InterTime);
		}

		/// 判断是否删除成功，不成功报错返回
		if(bOK==0)
		{
			return -1;
		}		
	}

	/// 创建目标目录，连续多次创建，成功直接跳出
	for( bOK=0,i=0; i<NumTry; i++ )
	{
		RV = _mkdir(Dir2);
		if(RV==0)
		{
			bOK = 1;
			break;
		}

		/// 两次重试间隔InterTime毫秒
		Sleep(InterTime);
	}

	/// 如果未创建成功，报错返回
	if(bOK==0)
		return -1;

	//保存进程当前目录
	GetCurrentDirectory(499,BakCurPath);

	/// 对源目录中所有文件和目录进行处理
	Dir1 = DirS;
	Dir1.TrimRight("\\");
	SetCurrentDirectory(Dir1);
	BOOL bWorking = myFind.FindFile("*.*");
	while (bWorking)
	{
		bWorking = myFind.FindNextFile();
		
		curName = myFind.GetFileName();
		if(curName=="."||curName=="..")
			continue;

		D1 = Dir1 + "\\" + curName;
		D2 = Dir2 + "\\" + curName;

		/// 如果是目录，递归调用本函数
		/// 注意，递归调用本函数中已经包含多次尝试，因此这里不再多次尝试
		if( myFind.IsDirectory() )
		{			
			RV = SuperCopyDir(D1,D2,NumTry,InterTime);
			if(RV!=1)
			{
				SetCurrentDirectory(BakCurPath);
				return -1;				
			}

			continue;
		}

		/// 走到这里，说明是文件，需要进行拷贝
		for( bOK=0,i=0; i<NumTry; i++ )
		{
			bSuc = CopyFile(D1,D2,FALSE);
			if(bSuc)
			{
				bOK=1;
				break;
			}

			/// 两次重试间隔InterTime毫秒
			Sleep(InterTime);
		}
		if(bOK==0)
		{
			SetCurrentDirectory(BakCurPath);
			return -1;
		}
	}

	/// 到这里，说明目录中所有子目录和文件都拷贝成功了
	SetCurrentDirectory(BakCurPath);
	return 1;
}

/***************************************
功能 ：删除目录以及目录中的文件和子目录。
算法：递归调用
***************************************/
int MpGlobalFunc::RMDirAndFiles(LPCSTR pDirName)
{
	CFileFind myFind;
	CString DirName,curName;
	char BakCurPath[500];
	int i,RV=1;
	
	//保存进程当前目录
	GetCurrentDirectory(499,BakCurPath);	

	DirName.Format("%s",pDirName);
	DirName.TrimRight("\\");
	
	if(!myFind.FindFile(DirName))
	{
		return 1;
	}

	SetCurrentDirectory(DirName);
	BOOL bWorking = myFind.FindFile("*.*");
	while (bWorking)
	{
		bWorking = myFind.FindNextFile();
		
		curName = myFind.GetFileName();
		if(curName=="."||curName=="..")
			continue;

		//如果是目录，则递归调用
		curName = myFind.GetFilePath();
		if(myFind.IsDirectory())
		{
			i = RMDirAndFiles(LPCSTR(curName));
			RV = (i==-1) ? -1 : RV ;
		}
		
		//如果是文件，则删除文件
		else
		{
			i=0;
			while( !DeleteFile( LPCSTR(curName) ) && ++i<=500 )
			{
				;
			}
			RV = (i>=500) ? -1 : RV ;
		}
		
		if(RV==-1)
			break;
	}
	myFind.Close();

	//恢复进程当前目录
	bWorking = SetCurrentDirectory(BakCurPath);

	if(RV==-1)
		return RV;

	//删除目标目录
	i=0;	
	while( _rmdir(DirName)!=0 && ++i<=500 )
	{
		;
	}
	RV = (i>=500) ? -1 : RV ;
	return RV;
}

// 移动一个窗口到另一个窗口的特定位置
int MpGlobalFunc::PositionWndAroundAnother(CWnd & WndStill, CWnd & Wnd2Move, int iPosFlag)	//未测试！！！！！！
{
	CRect rectWnd2Move;
	Wnd2Move.GetWindowRect(&rectWnd2Move);

	CRect rectWndStill;
	WndStill.GetWindowRect(&rectWndStill);

	switch(iPosFlag)
	{
	case 1: // 底边左对齐
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left,
				rectWndStill.bottom,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 2: // 左下角 (外)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.bottom,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 3: // 左边下对齐
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 4: // 左边上对齐
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.top,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 5: // 左上角 (外)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break; 

	case 6: // 顶边左对齐
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left,
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 7: // 顶边右对齐
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;
	case 8: // 右上角 (外)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right,
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 9: // 右边上对齐
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right,
				rectWndStill.top,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 10: // 右边下对齐
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right,
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 11: // 右下角 (外)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right,
				rectWndStill.bottom,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 12: // 底边右对齐
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.bottom,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 13: // 左下角 (内)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left,
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 14: // 左上角 (内)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left,
				rectWndStill.top,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break; 
	case 15: // 右上角 (内)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.top,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 16: // 右下角 (内)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 17: // 中央
		{
			Wnd2Move.MoveWindow(
				(rectWndStill.right - rectWndStill.left)/2 - (rectWnd2Move.right - rectWnd2Move.left)/2,
				(rectWndStill.bottom - rectWndStill.top)/2 - (rectWnd2Move.bottom - rectWnd2Move.top)/2,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	default:
		break;
	}

	return 0;
}

int MpGlobalFunc::RMObjFileInPath(const char * Path,const char *ObjStr)
{
	CFileFind myFind;
	CString DirName,curName,curName0;
	char BakCurPath[500];
	int i,RV=1;
	
	//保存进程当前目录
	GetCurrentDirectory(499,BakCurPath);	

	DirName.Format("%s",Path);
	DirName.TrimRight("\\");
	
	if(!myFind.FindFile(DirName))
		return 1;

	SetCurrentDirectory(DirName);
	BOOL bWorking = myFind.FindFile("*.*");
	while (bWorking)
	{
		bWorking = myFind.FindNextFile();
		
		curName0 = myFind.GetFileName();
		if(curName=="."||curName=="..")
			continue;

		//如果是目录，跳下一个
		curName = myFind.GetFilePath();
		if(myFind.IsDirectory())
		{
			continue;
		}
		
		//如果是文件，需要判断文件是否需要删除
		else
		{
			/// 如果文件名中不包含目标字符串，跳下一个
			if( curName0.Find(ObjStr)==-1 )
				continue;

			/// 删除当前文件
			i=0;
			while( !DeleteFile( LPCSTR(curName) ) && ++i<=500 )
			{
				;
			}
			RV = (i>=500) ? -1 : RV ;
		}
		
		if(RV==-1)
			break;
	}
	myFind.Close();

	//恢复进程当前目录
	bWorking = SetCurrentDirectory(BakCurPath);

	if(RV==-1)
		return RV;
	else
		return RV;
}

//整型数字字符串合法性检查函数
int MpGlobalFunc::StringIsInt(const char * pObj)
{
	int RV;
	const char *pChar;

	/// 先检查是否数字
	RV = StringIsNumber(pObj);
	if(RV==0)
		return 0;

	/// 再检查有无小数点
	pChar = strstr(pObj,".");
	if(pChar)
		return 0;

	return 1;
}

//数字字符串合法性检查函数
int MpGlobalFunc::StringIsNumber(const char * pObj)
{
	int i;
	int len;
	int posDot;		//小数点的位置
	int numDot;		//小数点的数目
	int posMinus;	//负号的位置
	int numMinus;	//负号的数目
	char maxChar='9';
	char minChar='0';

	len = strlen(pObj);
	if(len==0)
		return 0;
	
	//统计负号与小数点的位置和数目
	posMinus = posDot = -1;
	numMinus = numDot = 0;
	for(i=0;i<len;i++)
	{
		if( pObj[i]=='-' )
		{
			posMinus = i;
			numMinus++;
		}
		if( pObj[i]=='.' )
		{
			posDot = i;
			numDot++;
		}
	}

	//如果有负号，只能出现一次且必须在第一个字符位置，否则非法
	if(numMinus>1)
		return 0;
	else if( (numMinus==1) && (posMinus!=0) )
		return 0;

	//如果有小数点，只能出现一次且其前后都必须有数字
	if(numDot>1)
		return 0;
	else if( (numDot==1) && ( (posDot==0)||(posDot==len-1) ) )
		return 0;

	//统计除了负号和小数点之外，是否每个字符都在0到9之间
	for(i=0;i<len;i++)
	{
		if( (pObj[i]=='-') || (pObj[i]=='.') )
			continue;
		if( pObj[i]<minChar || pObj[i]>maxChar )
			return 0;	
	}

	return 1;
}

int MpGlobalFunc::ChkRegionCover(double xl,double yb,double xr,double yt,int NumReg,
						  double *pXl,double *pYb,double *pXr,double *pYt,double cs)
{
	char **ppF;
	int Hei,Wid,i,j,SN;
	double xo,yo,x,y;
	double x1,y1,x2,y2;
	int i1,i2,j1,j2;
	
	/// 判断输入数据范围是否正常
	if(cs<=0)
		return -1;

	if( xl>=xr || yb>=yt )
		return -1;

	if(NumReg<0)
		return -1;

	if( NumReg==0 )
		return 0;

	xo = xl + 0.5*cs;
	yo = yb + 0.5*cs;

	/// 计算标志矩阵象素尺寸
	Hei = (int)ceil( (yt-yb)/cs );
	Wid = (int)ceil( (xr-xl)/cs );

	/// 为标志矩阵分配内存空间
	/// 注意标志矩阵行号从下向上递增
	ppF = (char **)NewSpace2D(Hei,Wid,sizeof(char));
	if(ppF==NULL)
		return -1;

	/// 初始化标志矩阵
	for( i=0; i<Hei; i++ )
		memset( ppF[i], 0, sizeof(char)*Wid );

	/// 可以使用两种方法进行标记，

	/// 方法一
	/// 一种是对目标区域每个象素循环，
	/// 判断每个象素是否位于其中一个矩形内，位于则标记1

	/// 方法二
	/// 对矩形区域进行循环，找到每个矩形区域与目标区域的重合部分，
	/// 将重合部分的象素标记为1

	/// 方法一代码 未完成
	/*if(0)
	{
		/// 对每个象素逐个循环，判断是否被矩形数组覆盖
		for( i=0; i<Hei; i++ )
		{
			for( j=0; j<Wid; j++ )
			{
				x = xo + j*cs;
				y = yo + i*cs;


			}
		}
	}*/

	/// 方法二代码
	if(1)
	{
		/// 对每个矩形区域逐个处理
		for( SN=0; SN<NumReg; SN++ )
		{
			/// 找到重合部分
			x1 = max( xl, pXl[SN] );
			x2 = min( xr, pXr[SN] );
			y1 = max( yb, pYb[SN] );
			y2 = min( yt, pYt[SN] );

			/// 判断重合部分是否有效，无效则直接跳下一循环
			if( x1>=x2 || y1>=y2 )
				continue;

			/// 有效则计算该有效区域对应的坐标象素范围
			j1 = int((x1-xo)/cs + 0.5);
			j2 = int((x2-xo)/cs + 0.5);
			i1 = int((y1-yo)/cs + 0.5);
			i2 = int((y2-yo)/cs + 0.5);

			/// 将重合范围标记为1
			for( i=i1; i<=i2; i++ )
			{
				for( j=j1; j<=j2; j++ )
				{
					ppF[i][j] = 1;
				}
			}
		}
	}

	/// 全部标记完，判断是否被覆盖
    for( i=0; i<Hei; i++ )
	{
		for( j=0; j<Wid; j++ )
		{
			/// 一旦发现有象素未被标记，说明区域未覆盖，返回
			if(ppF[i][j]==0)
			{
				MpGlobalFunc::DeleteSpace2D(ppF);
				return 0;
			}
		}
	}

	/// 走到这里，说明全部标记为1，目标区域被覆盖
	MpGlobalFunc::DeleteSpace2D(ppF);
	return 1;
}

unsigned char ** MpGlobalFunc::LineToPic24bit(int Num,double *pX,double *pY,double cs,int &Hei,int &Wid)
{
	unsigned char **ppTmp;
	unsigned char **ppD;
	int i,j;

	ppTmp = LineToPic8bit(Num,pX,pY,cs,Hei,Wid);
	if(ppTmp==NULL)
		return NULL;

	ppD = (unsigned char **)NewSpace2D(Hei,Wid*3,sizeof(char));
	for( i=0; i<Hei; i++ )
	for( j=0; j<Wid; j++ )
	{
		ppD[i][j*3] = ppD[i][j*3+1] = ppD[i][j*3+2] = ppTmp[i][j];
	}

	MpGlobalFunc::DeleteSpace2D(ppTmp);
	return ppD;
}

unsigned char ** MpGlobalFunc::LineToPic8bit(int Num,double *pX,double *pY,double cs,int &Hei,int &Wid)
{
	double minX,maxX,minY,maxY;
	int i,j;
	unsigned char ** ppD=NULL;
	int NumPoi;
	int *pXc,*pYc;

	/// 统计折线坐标范围
	minX = maxX = pX[0];
	minY = maxY = pY[0];
	for( i=0; i<Num; i++ )
	{
		minX = min(minX,pX[i]);
		maxX = max(maxX,pX[i]);
		minY = min(minY,pY[i]);
		maxY = max(maxY,pY[i]);
	}

	/// 计算矩阵尺寸并分配空间
	Hei = (int)ceil( (maxY-minY)/cs ) + 1;
	Wid = (int)ceil( (maxX-minX)/cs ) + 1;
	ppD = (unsigned char **)NewSpace2D(Hei,Wid,sizeof(char));
	if(ppD==NULL)
		return NULL;
	for( i=0; i<Hei; i++ )
		memset( ppD[i], 0, sizeof(char)*Wid );

	/// 对折线逐个处理，点画线
	for( i=0; i<Num-1; i++ )
	{

		/// 两点画线
		pXc = pYc = NULL;
		PointToLine( pX[i],pY[i],pX[i+1],pY[i+1], minX,minY,cs, NumPoi,pXc,pYc );

		/// 将当前线段覆盖的象素置为1
		for( j=0; j<NumPoi; j++ )
		{
			ppD[ Hei-1-pYc[j] ][ pXc[j] ] = 1;
		}

		delete [] pXc;
		delete [] pYc;
	}

	/// 背景从0改为255
	for( i=0; i<Hei; i++ )
	for( j=0; j<Wid; j++ )
	{
		if(ppD[i][j]==0)
			ppD[i][j] = 255;
		else
			ppD[i][j] = 0;
	}

	return ppD;
}

int MpGlobalFunc::GetObjFileNameInPath(const char *Path,const char *ObjStr,CString &FN,int SN)
{
	CFileFind myFind;
	CString DirName,curName;
	char BakCurPath[500];
	int RV=1;
	int Num=0;

	FN = "";
		
	//保存进程当前目录
	GetCurrentDirectory(499,BakCurPath);	

	DirName.Format("%s",Path);
	DirName.TrimRight("\\");
	
	if(!myFind.FindFile(DirName))
		return 1;

	SetCurrentDirectory(DirName);
	BOOL bWorking = myFind.FindFile("*.*");
	while (bWorking)
	{
		bWorking = myFind.FindNextFile();
		
		curName = myFind.GetFileName();
		if(curName=="."||curName=="..")
			continue;

		//如果是目录，眺至下一个
		//curName = myFind.GetFilePath();
		if(myFind.IsDirectory())
		{
			continue;
		}
		
		//如果是文件，需要判断文件是否需要删除
		else
		{
			/// 如果文件名中不包含目标字符串，跳下一个
			if( curName.Find(ObjStr)==-1 )
				continue;

			/// 数目加1
			Num++;

			//如果数目与序号一致，输出
			if(Num==SN)
			{
				FN = curName;
				return 1;
			}
		}
	}
	myFind.Close();

	//恢复进程当前目录
	bWorking = SetCurrentDirectory(BakCurPath);

	return 0;
}

int MpGlobalFunc::GetObjFileNumInPath(const char *Path,const char *ObjStr,int &Num)
{
	CFileFind myFind;
	CString DirName,curName;
	char BakCurPath[500];
	int RV=1;

	Num = 0;
	
	//保存进程当前目录
	GetCurrentDirectory(499,BakCurPath);	

	DirName.Format("%s",Path);
	DirName.TrimRight("\\");
	
	if(!myFind.FindFile(DirName))
		return 1;

	SetCurrentDirectory(DirName);
	BOOL bWorking = myFind.FindFile("*.*");
	while (bWorking)
	{
		bWorking = myFind.FindNextFile();
		
		curName = myFind.GetFileName();
		if(curName=="."||curName=="..")
			continue;

		//如果是目录，跳下一个
		//curName = myFind.GetFilePath();
		if(myFind.IsDirectory())
		{
			continue;
		}
		
		//如果是文件，需要判断文件是否需要删除
		else
		{
			/// 如果文件名中不包含目标字符串，跳下一个
			if( curName.Find(ObjStr)==-1 )
				continue;

			/// 数目加1
			Num++;
		}
	}
	myFind.Close();

	//恢复进程当前目录
	bWorking = SetCurrentDirectory(BakCurPath);

	return 1;
}

/// 检查线段是否在圆形外面，是返回1，不是返回0，错误返回-1
int MpGlobalFunc::IsLineOutCircle(double x1,double y1,double x2,double y2,
								  double xo,double yo,double Rad,double minV)
{
	double pa,pb,pc;	//< 直线方程参数
	double dis;
	double xc1,yc1,xc2,yc2;
	int RV;
	double v1,v2;

	/// 求直线方程
	RV = PointToLineEquation(x1, y1,x2, y2,pa, pb,pc,minV);
	if(RV!=1)
		return -1;

	/// 求圆心到直线距离
	dis = fabs( pa*xo + pb*yo + pc ) / sqrt( pa*pa + pb*pb );

	/// 如果圆心到直线距离大于半径，线段一定在圆外面
	if(dis>Rad)
		return 1;

	/// 判断线段所在直线与圆的交点，交点一般为两个，
	/// 如果直线与圆相切，两个角点合为一个
	RV = Cal_LineCircle_Cross( pa, pb, pc, xo, yo, Rad, xc1, yc1, xc2, yc2,minV);
	if(RV==-1)
		return -1;

	/// 判断两个交点是否在线段外面
	/// 如果线段两点x差值大，对x进行判断，否则对y进行判断
	if( fabs(y2-y1)>fabs(x2-x1) )
	{
		v1 = min(y1,y2);
		v2 = max(y1,y2);
		if( (yc1<v1||yc1>v2) && (yc2<v1||yc2>v2) )
			return 1;
		else
			return 0;
	}
	else
	{
		v1 = min(x1,x2);
		v2 = max(x1,x2);
		if( (xc1<v1||xc1>v2) && (xc2<v1||xc2>v2) )
			return 1;
		else
			return 0;
	}
}

/// 将字符串首尾的空格去掉，注意直接在字符串指针中改写
int MpGlobalFunc::StrDelSpace(char * pObj)
{
	int len,i;
//	char *pPos;
	char *buf;

	if( pObj==NULL || strlen(pObj)==0 )
		return 1;

	len = strlen(pObj);

	buf = new char[len+1];
	if(buf==NULL)
		return -1;

	/// 去除后边空格
	if(pObj[len-1]==' ')
	{
		/// 从后向前找到第一个不是空格的字符，直到开始
		for(i=len-2;i>=0;i--)
		{
			if(pObj[i]!=' ')
			{			
				break;
			}
		}
		pObj[i+1] = '\0';
	}

	/// 去除前边空格
	if(pObj[0]==' ')
	{
		/// 从前向后找到第一个不是空格的字符，直到结束符
		for(i=1;i<len;i++)
		{
			if(pObj[i]!=' ')
			{			
				break;
			}
		}

		/// 将i指示位置之后的内容向前拷贝，用一个内存过度
		strcpy(buf,pObj);
		strcpy(pObj,buf+i);
	}

	delete [] buf;
	return 1;
}

/// 计算直线与圆的交点，返回0，无交点，返回1，两个交点，返回2，一个交点，返回－1，出错
int MpGlobalFunc::Cal_LineCircle_Cross(double pa,double pb,double pc,double xo,double yo,
						double R,double &xc1,double &yc1,double &xc2,double &yc2,double minV)
{
	double m,n;
	double a,b,c;
	double delta,dis;

	/// 根据直线方程的参数判断对x进行计算还是对y进行计算
	if( fabs(pb)>fabs(pa) )	//<  对x进行计算
	{
		/// 将ax+by+c=0转换为y=mx+n的形式
		m = -1*pa/pb;
		n = -1*pc/pb;

		/// 将y=mx+n代入圆方程，形成关于x的二次方程
		a = 1 + m*m;
		b = -2*xo + 2*m*(n-yo);
		c = xo*xo + (n-yo)*(n-yo) - R*R;

		/// 根据二次方程求解
		delta = b*b - 4*a*c;

		/// 无交点
		if(delta<0)
			return 0;

		/// 计算两个根
		xc1 = ( -1*b + sqrt(delta) ) / (2*a);
		yc1 = m*xc1 + n;
		xc2 = ( -1*b - sqrt(delta) ) / (2*a);
		yc2 = m*xc2 + n;

		/// 判断根的个数
		dis = sqrt( (xc1-xc2)*(xc1-xc2)+(yc1-yc2)*(yc1-yc2) );
        if( dis<minV)
		{
			return 2;
		}
		else
		{
			return 1;
		}
	}

	else	//<  对y进行计算
	{
		/// 将ax+by+c=0转换为x=my+n的形式
		m = -1*pb/pa;
		n = -1*pc/pa;

		/// 将x=my+n代入圆方程，形成关于y的二次方程
		a = 1 + m*m;
		b = -2*yo + 2*m*(n-xo);
		c = yo*yo + (n-xo)*(n-xo) - R*R;

		/// 根据二次方程求解
		delta = b*b - 4*a*c;

		/// 无交点
		if(delta<0)
			return 0;

		/// 计算两个根
		yc1 = ( -1*b + sqrt(delta) ) / (2*a);
		xc1 = m*yc1 + n;
		yc2 = ( -1*b - sqrt(delta) ) / (2*a);
		xc2 = m*yc2 + n;

		/// 判断根的个数
		dis = sqrt( (xc1-xc2)*(xc1-xc2)+(yc1-yc2)*(yc1-yc2) );
        if( dis<minV)
		{
			return 2;
		}
		else
		{
			return 1;
		}
	}
}
//函数功能：给定相机位置、观察点位置、视场角度、求取所生成的梯形区域大小
//输入：szy20141023 add 
//eye 相机位置(有带号的高斯坐标)
//view 观察位置(有带号的高斯坐标)
//fov_x 视场纵角度 弧度
//fov_y 视场纵角度 弧度
//输出：result 梯形区域大小(无带号高斯坐标)
////调试所用数据（为了方便可以在函数入口处进行检测）
//szypolygen_GS eye1,view1;
//eye1.p_x = 16528629.71;
//eye1.p_y = 4484293.03;
//eye1.p_z = 3284.715;
//vector<szypolygen_GS>Eye;
//
//vector<szypolygen_GS>Result;
//Eye.push_back(eye1);
//
//view1.p_x = 16524888.59;
//view1.p_y = 4480170.31;
//view1.p_z = 929.046;
//vector<szypolygen_GS>sss;
//sss.push_back(view1);
//MpGlobalFunc a;

//int f = a.DL_calc_JXQY(Eye,sss, 4.6*qxwPI/180, 3.7*qxwPI/180, Result);


int MpGlobalFunc::DL_calc_JXQY(vector<szypolygen_GS> eye, vector<szypolygen_GS> view, double fov_x, double fov_y, vector<szypolygen_GS> &result)
{
	int nDH;
	double angle21;
	if (eye[0].p_x >= 1000000.0)//建筑物坐标含带号
	{
		nDH = int(eye[0].p_x / 1000000);//获取带号
		eye[0].p_x = eye[0].p_x - nDH * 1000000;//去掉带号
	}
	else
	{
		AfxMessageBox("所提供的位置信息无带号，请添加");
		return 0;
	}
	if (view[0].p_x >= 1000000.0)//建筑物坐标含带号
	{
		nDH = int(view[0].p_x / 1000000);//获取带号
		view[0].p_x = view[0].p_x - nDH * 1000000;//去掉带号
	}
	else
	{
		AfxMessageBox("所提供的位置信息无带号，请添加");
		return 0;
	}
	//根据由AB、AH线段长度求取HB的长度，其中H点是A点投影到地面上与地面垂直的点
	double AB = sqrt(pow(eye[0].p_x-view[0].p_x, 2)
		+pow(eye[0].p_y-view[0].p_y, 2)
		+pow(eye[0].p_z-view[0].p_z, 2)
		);
	double AH = eye[0].p_z-view[0].p_z;
	double HB = sqrt(AB*AB-AH*AH);
	//获取角HAB的角度，并根据角HAB的角度和已知的角CAB的角度求取角HAC的角度
	double jiao_HAB = atan(HB/AH);
	double jiao_HAC = jiao_HAB-fov_x/2;
	double HC = AH * tan(jiao_HAC);
	//得到BC、BD长度
	double BC = HB - HC;
	double HD = AH * tan(jiao_HAC+fov_x);
	double BD = HD -HB;

	szypolygen_GS C,E,F,D,G,N;
	szypolygen_DD c,e,f,d,g,n;
	//求取C点的高斯坐标
	C.p_x = HC/HB*(view[0].p_x-eye[0].p_x)+eye[0].p_x;
	C.p_y = HC/HB*(view[0].p_y-eye[0].p_y)+eye[0].p_y;
	C.p_z =view[0].p_z;
	//求取AC的线段长度
	double AC = sqrt(pow(eye[0].p_x-C.p_x, 2)
		+pow(eye[0].p_y-C.p_y, 2)
		+pow(eye[0].p_z-C.p_z, 2)
		);
	//根据C点的高斯坐标求取E、F点的经纬度坐标
	double EC = AC*tan(fov_y/2);
	double jiao_CB = atan((view[0].p_y-C.p_y)/(view[0].p_x-C.p_x));
	//获取C点到E点、F点的角度
	double jiao_CE = jiao_CB*180/qxwPI+90;
	double jiao_CF = jiao_CB*180/qxwPI-90;
	//转换前先将C点的高斯坐标转换为经纬度坐标
	XYtoBL(C.p_y,C.p_x,nDH*6-3,c.DD_y,c.DD_x);
	calc_coord_by_distance(c.DD_x, c.DD_y, EC/1000, jiao_CE, e.DD_x, e.DD_y, angle21);
	calc_coord_by_distance(c.DD_x, c.DD_y, EC/1000, jiao_CF, f.DD_x, f.DD_y, angle21);

	//求取D点的高斯坐标
	D.p_x = HD/HB*(view[0].p_x-eye[0].p_x)+eye[0].p_x;
	D.p_y = HD/HB*(view[0].p_y-eye[0].p_y)+eye[0].p_y;
	D.p_z =view[0].p_z;
	//求取AD的线段长度
	double AD = sqrt(pow(eye[0].p_x-D.p_x, 2)
		+pow(eye[0].p_y-D.p_y, 2)
		+pow(eye[0].p_z-D.p_z, 2)
		);
	//根据D点的高斯坐标求取G、N点的经纬度坐标
	double DG = AD*tan(fov_y/2);
	
	//获取D点到G点、N点的角度
	double jiao_DG = jiao_CB*180/qxwPI+90;
	double jiao_DN = jiao_CB*180/qxwPI-90;
	//转换前先将D点的高斯坐标转换为经纬度坐标
	XYtoBL(D.p_y,D.p_x,nDH*6-3,d.DD_y,d.DD_x);
	calc_coord_by_distance(d.DD_x, d.DD_y, DG/1000, jiao_DG, g.DD_x, g.DD_y, angle21);
	calc_coord_by_distance(d.DD_x, d.DD_y, DG/1000, jiao_DN, n.DD_x, n.DD_y, angle21);
	BLtoXY(e.DD_y,e.DD_x,nDH*6-3,E.p_y,E.p_x);
	BLtoXY(f.DD_y,f.DD_x,nDH*6-3,F.p_y,F.p_x);
	BLtoXY(n.DD_y,n.DD_x,nDH*6-3,N.p_y,N.p_x);
	BLtoXY(g.DD_y,g.DD_x,nDH*6-3,G.p_y,G.p_x);
	E.p_z = F.p_z = N.p_z = G.p_z = view[0].p_z;
	result.push_back(E);
	result.push_back(F);
	result.push_back(N);
	result.push_back(G);
	return 1;
}