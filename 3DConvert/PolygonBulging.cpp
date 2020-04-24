#include "stdafx.h"
#include "MpGlobalFunc.h"
#include "PolygonBulging.h"
#include <iostream>
#include <fstream>

//分配2维内存空间
void **NewSpace2D(int Row,int Col,short Size)
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
void ***NewSpace3D(int Len,int Row,int Col,short Size)
{
	int i,j;
	char ***pData;
	pData=(char ***)NewSpace2D(Len,Row,sizeof(char *));
	if(pData==NULL) return NULL;
	pData[0][0]=new char[Len*Row*Col*Size];
	if(pData[0][0]==NULL)
	{
		DeleteSpace2D(pData,Len);
		return NULL;
	}
	for(i=0;i<Len;i++)
	for(j=0;j<Row;j++)
		pData[i][j] = pData[0][0]+Row*Col*Size*i+Col*Size*j;
	return (void ***)pData;
}

//分配4维内存空间
void ****NewSpace4D(int Len1,int Len2,int Row,int Col,short Size)
{
	int i,j,k,Page,Vol;
	char ****pData;
	pData=(char ****)NewSpace3D(Len1,Len2,Row,sizeof(char *));
	if(pData==NULL) return NULL;
	pData[0][0][0]=new char[Len1*Len2*Row*Col*Size];
	if(pData[0][0][0]==NULL)
	{
		DeleteSpace3D(pData,Len1,Len2);
		return NULL;
	}
	Page = Row*Col*Size;
	Vol = Page*Len2;
	for(i=0;i<Len1;i++)
	for(j=0;j<Len2;j++)
	for(k=0;k<Row;k++)
		pData[i][j][k]=pData[0][0][0]+Vol*i+Page*j+Col*Size*k;
	return (void ****)pData;
}

//删除特定文件夹中所有的文件，不删除目录


/***************************************
功能 ：删除目录以及目录中的文件和子目录。
算法：递归调用
***************************************/

//两点画线，给出两端点浮点坐标，左下角点的浮点坐标和分辨率，返回两点间线段的坐标
int PointToLine( double x1,double y1,double x2,double y2,
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

//两点画线，给出两端点坐标，返回两点间线段的坐标
int PointToLine(int x1,int y1,int x2,int y2,int &NumPoi,int * (&pX),int *(&pY))
{
	double pa_a,pa_b,pa_c;
	int i,delta,curX,curY;
	
	if(pX!=NULL||pY!=NULL)
		return -1;
	
	//如果是一个点
	if(x1==x2&&y1==y2)
	{
		pX = new int[1];
		pY = new int[1];
		pX[0] = x1;
		pY[0] = y1;
		NumPoi = 1;
		return 1;
	}
	
	//求出直线方程
	PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c);

	//线扫描
	if( fabs(pa_b-1.0)<qxwLimMinValue )
	{
		NumPoi = abs(x1-x2)+1;
		pX = new int[NumPoi];
		pY = new int[NumPoi];
		delta = (x1<x2)?1:-1;
		for( i=0,curX=x1; i<NumPoi; i++,curX+=delta )
		{
			pX[i] = curX;
			pY[i] = (int)( -pa_a*pX[i]-pa_c +0.5 );
		}
	}
	else
	{
		NumPoi = abs(y1-y2)+1;
		pX = new int[NumPoi];
		pY = new int[NumPoi];
		delta = (y1<y2)?1:-1;
		for( i=0,curY=y1; i<NumPoi; i++,curY+=delta )
		{
			pY[i] = curY;
			pX[i] = (int)( -pa_b*pY[i]-pa_c +0.5) ;
		}
	}
	
	return 1;
}

//根据端点坐标求出直线标准方程
// para_a * x + para_b * y + c = 0
int PointToLineEquation(double x1,double y1,double x2,double y2,
						double &para_a,double &para_b,double &para_c)
{
	if( fabs(x1-x2)<qxwLimMinValue && fabs(y1-y2)<qxwLimMinValue )
		return -1;
	
	if( fabs(x1-x2)>fabs(y1-y2) )
	{
		para_a = -(y1-y2)/(x1-x2);
		para_c = -y1-para_a*x1;
		para_b = 1.0;
	}
	else
	{
		para_b = -(x1-x2)/(y1-y2);
		para_c = -x1-para_b*y1;
		para_a = 1.0;
	}
	
	return 1;
}


int PointToLineEquation(int X1,int Y1,int X2,int Y2,
						double &para_a,double &para_b,double &para_c)
{
	double x1,x2,y1,y2;

	if( abs(X1-X2)<qxwLimMinValue && abs(Y1-Y2)<qxwLimMinValue )
		return -1;
	
	x1 = (double)X1;
	x2 = (double)X2;
	y1 = (double)Y1;
	y2 = (double)Y2;
	if( fabs(x1-x2)>fabs(y1-y2) )
	{
		para_a = -(y1-y2)/(x1-x2);
		para_c = -y1-para_a*x1;
		para_b = 1.0;
	}
	else
	{
		para_b = -(x1-x2)/(y1-y2);
		para_c = -x1-para_b*y1;
		para_a = 1.0;
	}
	
	return 1;
}

//给出两点（点1，点2）坐标，计算点1到点2的方向，x为正东，y为正北
int TwoPoint2Angle(double x1,double y1,double x2,double y2,int bHudu,double &AngleRes)
{
	double x2subx1,y2suby1;

	x2subx1 = x2-x1;
	y2suby1 = y2-y1;
	
	if( fabs(x2subx1)<qxwLimMinValue )
		AngleRes = (y2>y1) ? qxwPI_DIV_2: -1*qxwPI_DIV_2 ;
	else if( fabs(y2suby1)<qxwLimMinValue )
		AngleRes = (x2>x1) ? 0.0 : qxwPI ;
	else if( x2subx1>0 && y2suby1>0 )		//第一象限
		AngleRes = atan(y2suby1/x2subx1);
	else if( x2subx1<0 && y2suby1>0 )		//第二象限
		AngleRes = qxwPI + atan(y2suby1/x2subx1);
	else if( x2subx1<0 && y2suby1<0 )		//第三象限
		AngleRes = atan(y2suby1/x2subx1)-qxwPI;
	else if( x2subx1>0 && y2suby1<0 )		//第四象限
		AngleRes = atan(y2suby1/x2subx1);

	if(bHudu==0)
		AngleRes = AngleRes*180/qxwPI;

	return 1;
}

/**** 求射线和线段的角点数目
给出射线端点坐标和方向角(正东为零，逆时针为正)，给出线段端点坐标
结果信息Res说明：
	1，射线与直线不重合，且与线段不相交
	2，射线与直线不重合，交点是普通点
	3，射线与直线不重合，交点是射线的端点
	4，射线与直线不重合，交点是线段的端点
	5，射线与直线不重合，交点是线段且射线的端点
	6，射线与直线重合，但线段不在射线上
	7，射线与直线重合，且线段在射线上
***************/
int RadialMeetLine2Poi(
					   double r_x0,		//射线横坐标
					   double r_y0,		//射线纵坐标
					   double r_dir,	//射线方向，正东为零，逆时针为正，弧度
					   double l_x1,		//线段端点坐标
					   double l_y1,
					   double l_x2,
					   double l_y2,
					   int &Res			//结果，详细信息见说明
					   )
{
	double pa_a,pa_b,pa_c;					//直线方程参数
	double l_xmin,l_xmax,l_ymin,l_ymax;		//线段坐标范围
	double r_xt,r_yt,r_l;					//临时点坐标
	double xPoi,yPoi;						//交点坐标
	double dis1,dis2,dismax,disPoi;
	int bPoiR;								//为射线端点
	int bPoiL;								//为线段端点
	
	bPoiR = bPoiL = 0;
	
	//求线段坐标范围
	l_xmin = (l_x1<l_x2) ? l_x1 : l_x2 ;
	l_xmax = (l_x1<l_x2) ? l_x2 : l_x1 ;
	l_ymin = (l_y1<l_y2) ? l_y1 : l_y2 ;
	l_ymax = (l_y1<l_y2) ? l_y2 : l_y1 ;
	
	//求直线方程
	PointToLineEquation( l_x1, l_y1, l_x2, l_y2, pa_a, pa_b, pa_c );

	//取射线上距顶点10处的临时点
	r_l = 10.0;
	r_xt = r_x0 + r_l*cos(r_dir);
	r_yt = r_y0 + r_l*sin(r_dir);

	//判断是否重合，如果两个点都在直线上，重合
	if( fabs(pa_a*r_x0+pa_b*r_y0+pa_c) < qxwLimMinValue &&
		fabs(pa_a*r_xt+pa_b*r_yt+pa_c) < qxwLimMinValue )
	{
		//求线段两端点到射线起点的距离
		if( fabs(pa_a-1)<qxwLimMinValue )
		{
			dis1 = (l_y1-r_y0)/sin(r_dir);
			dis2 = (l_y2-r_y0)/sin(r_dir);
		}
		else
		{
			dis1 = (l_x1-r_x0)/cos(r_dir);
			dis2 = (l_x2-r_x0)/cos(r_dir);
		}
		dismax = (dis1>dis2) ? dis1 : dis2 ;

		//判断线段是否在射线上，如果最大距离比零小，不在射线上
		if(dismax<0)
			Res = 6;
		else
			Res = 7;
		return 1;		
	}

	//此时射线与线段肯定不重合，求交点
	disPoi = pa_a*cos(r_dir) + pa_b*sin(r_dir) ;
	if( fabs(disPoi)<qxwLimMinValue )	//此时平行
	{
		Res = 1;	return 1;
	}

	//此时射线一定和直线相交，求交点
	disPoi = ( 0 - pa_c - pa_a*r_x0 - pa_b*r_y0 )/disPoi ;
	xPoi = r_x0+disPoi*cos(r_dir);
	yPoi = r_y0+disPoi*sin(r_dir);

	//如果交点是线段端点之一
	if( ( fabs(xPoi-l_x1)<qxwLimMinValue && fabs(yPoi-l_y1)<qxwLimMinValue ) ||
		( fabs(xPoi-l_x2)<qxwLimMinValue && fabs(yPoi-l_y2)<qxwLimMinValue ) )
	{
		bPoiL = 1;	
	}
	//如果交点在线段内
	else if( xPoi>=l_xmin && xPoi<=l_xmax && yPoi>=l_ymin && yPoi<=l_ymax )
	{
		bPoiL = 0;
	}
	//如果不是线段端点且交点坐标不在线段范围内，不相交
	else
	{
		Res = 1;	return 1;
	}
	
	//如果交点是射线起点
	if( fabs(xPoi-r_x0)<qxwLimMinValue && fabs(yPoi-r_y0)<qxwLimMinValue  )
	{
		bPoiR = 1;
	}

	//运行到这里，说明肯定相交，判断交点情况
	if( bPoiL==0 && bPoiR==0 )	//普通点
		return 2;
	else if( bPoiR==1 && bPoiL==0 )
		return 3;
	else if( bPoiR==0 && bPoiL==1 )
		return 4;
	else
		return 5;	
}

//数字字符串合法性检查函数
int StringIsNumber(const char * pObj)
{
	int i;
	int len;
	int posDot;		//小数点的位置
	int numDot;		//小数点的数目
	int posMinus;	//负号的位置
	int numMinus;	//负号的数目
	char maxChar='9';
	char minChar='0';

	len = (int)strlen(pObj);
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
		return -1;
	else if( (numMinus==1) && (posMinus!=0) )
		return -1;

	//如果有小数点，只能出现一次且其前后都必须有数字
	if(numDot>1)
		return -1;
	else if( (numDot==1) && ( (posDot==0)||(posDot==len-1) ) )
		return -1;

	//统计除了负号和小数点之外，是否每个字符都在0到9之间
	for(i=0;i<len;i++)
	{
		if( (pObj[i]=='-') || (pObj[i]=='.') )
			continue;
		if( pObj[i]<minChar || pObj[i]>maxChar )
			return -1;	
	}

	return 1;
}


//输入两个角度，求角度差，单位度
int DiffOfTwoDir(double dir1,double dir2,bool bDu,double &res)
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

//获取当前时间，并转换为字符串
int CreateCurTimeString(char * pObjStr)
{
	time_t curTime;
	struct tm *pCurTime;

	if(pObjStr==NULL)
		return -1;

	time(&curTime);
	pCurTime = localtime(&curTime);
	if (pCurTime == NULL)//20161212 304 静态分析问题更改
	{
		cout << "获取时间失败." << endl;
		return 0;
	}
	sprintf(pObjStr,"%s",asctime(pCurTime));
	pObjStr[strlen(pObjStr)-1]='\0';

	return 1;
}

//获取当前时间，并转换为字符串，同时纪录当前时间
int CreateCurTimeString(time_t &curTime, char * pObjStr)
{
	struct tm *pCurTime;

	if(pObjStr==NULL)
		return -1;

	time(&curTime);
	pCurTime = localtime(&curTime);
	if (pCurTime == NULL)//20161212 304 静态分析问题更改
	{
		cout << "获取时间失败." << endl;
		return 0;
	}
	sprintf(pObjStr,"%d-%02d-%02d %02d:%02d:%02d",
		pCurTime->tm_year+1900,pCurTime->tm_mon+1,pCurTime->tm_mday,
		pCurTime->tm_hour,pCurTime->tm_min,pCurTime->tm_sec);
	
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
int InterLinear1D( double x1, double v1 ,double x2, double v2, double x, double & v )
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

int InterLinear2D( double x1, double x2, double y1, double y2, double v11, double v12, 
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

int InterLinear3D( double x1, double x2, double y1, double y2, double z1, double z2, 
				  double v111, double v121, double v211, double v221, 
				  double v112, double v122, double v212, double v222, 
				  double x, double y, double z, double &v )
{
	double v1,v2;
	int RV;
	
	//首先对x、y进行插值
	RV = InterLinear2D( x1, x2, y1, y2, v111, v121, v211, v221 , x, y, v1 );
	if(RV==-1)	return -1;

	RV = InterLinear2D( x1, x2, y1, y2, v112, v122, v212, v222 , x, y, v2 );
	if(RV==-1)	return -1;

	//再对z进行插值
	RV = InterLinear1D( z1, v1, z2, v2, z, v );
	if(RV==-1)	return -1;

	return 1;
}


/*****************************
判断点与线段的关系。

  -1,两点太近，不进行计算
  1，点为线段的端点
  2，点为线段内的点
  3，点在直线上，在线段外，
  4，点不在直线上
******************************/
int Cal_PoiLine_Relation( double poiX, double poiY, double x1, double y1, 
					   double x2, double y2, double minV=1e-5 )
{
	double pa_a,pa_b,pa_c;
	double dis;
	int bOnLine,bPort,bInLine;

	//如果两点太近，不进行计算直接返回
	if( sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) ) < minV )
		return -1;

	bOnLine = bPort = bInLine = 0;		//初始化标志符
	PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c);	//求直线方程
	dis = pa_a*poiX + pa_b*poiY + pa_c;

	//判断点是否在直线上
	if( fabs(dis)<minV )
		bOnLine = 1 ;

	//判断点是否是线段端点
	if( 
		( fabs(poiX-x1)<minV && fabs(poiY-y1)<minV ) ||
		( fabs(poiX-x2)<minV && fabs(poiY-y2)<minV )
		)
		bPort = 1;
	
	//判断点是否在线段两端点坐标范围内
	if( fabs(pa_a-1)>=fabs(pa_b-1) )	//斜率小于1，横着的线段
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

/*****************************
判断线段与线段的关系。

 Type = 0，线段与线段重合
		1，线段与线段平行
		2，线段与线段所在直线相交
 Line1Type =	0，交点在线段1外
				1，交点为线段1端点
				2，交点在线段1内
 Line2Type =	0，交点在线段2外
				1，交点为线段2端点
				2，交点在线段2内
  返回值 
	1，成功
	-1，无法判断关系，可能程序出错或者两点距离太近
******************************/
int Cal_TwoLine_Relation( double x1, double y1, double x2, double y2,
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
		return -1;

	//计算两线段直线方程
	PointToLineEquation(x1,y1,x2,y2,a1,b1,c1);
	PointToLineEquation(x3,y3,x4,y4,a2,b2,c2);

	//取线段1的中点
	x0 = (x1+x2)/2;
	y0 = (y1+y2)/2;

	//计算线段1三个点到线段2所在直线的距离
	dis1 = fabs(a2*x1+b2*y1+c2) / sqrt(a2*a2+b2*b2) ;
	dis2 = fabs(a2*x2+b2*y2+c2) / sqrt(a2*a2+b2*b2) ;
	dis3 = fabs(a2*x0+b2*y0+c2) / sqrt(a2*a2+b2*b2) ;

	//判断射线与直线的关系
	if( fabs(dis1-dis2)<=minV && fabs(dis2-dis3)<=minV && fabs(dis1)<=minV )
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
		( fabs(a1-1)>=fabs(b1-1) && croX>minX && croX<maxX ) ||
		( fabs(a1-1)<fabs(b1-1) && croY>minY && croY<maxY )
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
		( fabs(a2-1)>=fabs(b2-1) && croX>minX && croX<maxX ) ||
		( fabs(a2-1)<fabs(b2-1) && croY>minY && croY<maxY )
		)
		Line2Type = 2;
	else
		Line2Type = 0;

	return 1;
}


int Cal_TwoLine_Cross( double a1, double b1, double c1, double a2, 
					  double b2, double c2, double &x, double &y )
{
	double Up,Down;

	Up = (b1*c2-b2*c1);
	Down = (a1*b2-a2*b1);

	if( fabs(Up)>fabs(Down) && fabs(Down)/fabs(Up)<1e-300 )
		return -1;
	x = Up/Down;	

	Up = (c2*a1-c1*a2);
	Down = (b1*a2-b2*a1);
	if( fabs(Up)>fabs(Down) && fabs(Down)/fabs(Up)<1e-300 )
		return -1;
	y = Up/Down;
	
	return 1;
}



/*****************************
判断射线与线段的关系。射线用端点与方向表示，
(r_x,r_y)为端点坐标，r_dir为方向，东零，逆时针正，单位度
线段用两端点表示，(x1,y1)(x2,y2)为端点坐标

 type = 0，射线与线段重合
		1，射线与线段平行
		2，射线与线段所在直线相交
 RadType =	0，交点不在射线上
			1，交点为射线端点
			2，交点在射线上
 LineType = 0，交点在线段外
			1，交点为线段端点
			2，交点在线段内
******************************/
int Cal_RadialLine_Relation( double r_x, double r_y, double r_dir, 
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

	ang = r_dir/180*qxwPI;	//转化为弧度

	//在射线上以一定间隔距离（len）取两点
	rx1 = r_x + len*cos( ang );
	ry1 = r_y + len*sin( ang );
	rx2 = rx1 + len*cos( ang );
	ry2 = ry1 + len*sin( ang );

	//计算线段直线方程
	PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c);

	//计算射线两点与直线的距离
	dis1 = fabs( r_x*pa_a + r_y*pa_b + pa_c ) / sqrt( pa_a*pa_a + pa_b*pa_b );
	dis2 = fabs( rx1*pa_a + ry1*pa_b + pa_c ) / sqrt( pa_a*pa_a + pa_b*pa_b );
	dis3 = fabs( rx2*pa_a + ry2*pa_b + pa_c ) / sqrt( pa_a*pa_a + pa_b*pa_b );

	//判断射线与直线的关系
	if( fabs(dis1-dis2)<=minV && fabs(dis2-dis3)<=minV && fabs(dis1)<=minV )
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

	//计算交点坐标
	crossX = ( pa_b*r_c - r_b*pa_c ) / ( pa_a*r_b - r_a*pa_b );
	crossY = ( r_c*pa_a - pa_c*r_a ) / ( pa_b*r_a - r_b*pa_a );

	//判断交点与射线的关系
	dis1 = sqrt( (crossX-r_x)*(crossX-r_x) + (crossY-r_y)*(crossY-r_y) );
	sinA = (crossY-r_y)/dis1;
	cosA = (crossX-r_x)/dis1;
	if( dis1<minV )
		RadType = 1 ;
	else if( sin(ang)*sinA>=0 && cos(ang)*cosA>=0 )
		RadType = 2 ; 
	else
		RadType = 0;

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
		( fabs(pa_a-1)>=fabs(pa_b-1) && crossX>minX && crossX<maxX ) ||
		( fabs(pa_a-1)<fabs(pa_b-1) && crossY>minY && crossY<maxY )
		)
		LineType = 2;
	else
		LineType = 0;
	
	return 1;
}



/**********
判断点与多边形的关系，内部还是外部。

  0，内部
  1，外部
  -1，计算出错或者本函数未能有效判断出二者关系
******************************/
int Cal_PoiPolygon_Relation( int NumPoi, double *pX, double *pY, 
						  double PoiX, double PoiY, double minV )
{
	int CroDotNumL,CroDotNumR;	//左右交点个数
	double Dir0=0;				//起始方向，0度
	double DirDlt=20;			//方向步长，45度
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

	//首先判断点是否在多边形的边上，如果点在任一一边上，
	//说明点在多边形内部
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

		//判断点是否在当前的边上，如果在，认为点在多边形内部，直接返回
		RV = Cal_PoiLine_Relation( PoiX, PoiY, x1, y1, x2, y2, minV );
		if( RV==1 || RV==2 )
			return 0;		
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
	angView = atan2( minLen/2, maxDis );
	angView *= 2;
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
		if( (CroDotNumL+CroDotNumR)%2==0 && (CroDotNumL+CroDotNumR)>0 )
		{
			if( CroDotNumL%2==0 && CroDotNumR%2==0 )
				evidenceOut ++;
			else
				evidenceIn ++;
		}

		//如果内部或外部的证据次数大于等于2，退出循环
		if( evidenceOut>=EviNumLim || evidenceIn>=EviNumLim )
			break;
	}

	if(evidenceIn>=EviNumLim)
		return 0;
	else if(evidenceOut>=EviNumLim)
		return 1;
	else
		return -1;
}



/************************
判断两个多边形之间的关系

  -1，计算出错
  0，分离
  1，相交
  2，1在2内部
  3，2在1内部
  4，近邻
  5，重合
**************************/
int Cal_TwoPolygon_Relation( int NP1, double *pX1, double *pY1,
							int NP2, double *pX2, double *pY2, double minV )
{
	double curX,curY;
	int i,j;
	int RV;
	int eviOut,eviIn;
	int relation12,relation21;		//1，内部，2，外部，3内外都有
	double x1,x2,x3,x4,y1,y2,y3,y4;
	int Type,Line1Type,Line2Type;

	//首先对第一个多边形的每个顶点与第二个多边形的关系进行判断
	eviOut = eviIn = 0;
	for( i=0; i<NP1; i++ )
	{
		curX = pX1[i];
		curY = pY1[i];
		RV = Cal_PoiPolygon_Relation( NP2, pX2, pY2, curX, curY, minV );
		if( RV==0 )
			eviIn++;
		else if(RV==1)
			eviOut++;
	}
	if( eviIn>0 && eviOut==0 )
		relation12 = 1;
	else if( eviIn==0 && eviOut>0 )
		relation12 = 2;
	else
		relation12 = 3;

	//对第二个多边形的每个顶点与第一个多边形的关系进行判断
	eviOut = eviIn = 0;
	for( i=0; i<NP2; i++ )
	{
		curX = pX2[i];
		curY = pY2[i];
		RV = Cal_PoiPolygon_Relation( NP1, pX1, pY1, curX, curY, minV );
		if( RV==0 )
			eviIn++;
		else if(RV==1)
			eviOut++;
	}
	if( eviIn>0 && eviOut==0 )
		relation21 = 1;
	else if( eviIn==0 && eviOut>0 )
		relation21 = 2;
	else
		relation21 = 3;

	//对边两个多边形边相交情况进行判断，
	//如果多边形相交，那么两多边形一定分别至少有一条边与对方的边线段内交
	//反过来也成立
	for( i=0; i<NP1; i++ )
	{
		x1 = pX1[i];
		y1 = pY1[i];
		if(i==NP1-1)
		{
			x2 = pX1[0];
			y2 = pY1[0];
		}
		else
		{
			x2 = pX1[i+1];
			y2 = pY1[i+1];
		}

		for( j=0; j<NP2; j++ )
		{
			x3 = pX2[j];
			y3 = pY2[j];
			if(j==NP2-1)
			{
				x4 = pX2[0];
				y4 = pY2[0];
			}
			else
			{
				x4 = pX2[j+1];
				y4 = pY2[j+1];
			}

			Cal_TwoLine_Relation( x1,y1,x2,y2,x3,y3,x4,y4, Type, Line1Type, Line2Type );

			//如果有边相交，则多边形一定相交
			if( Type==2 && Line1Type==2 && Line2Type==2 )
			{
				return 1;
			}
		}
	}

	//根据上述两个关系判断两多边形的关系
	if( relation12==1 && relation21==1  )
		return 5;
	else if( relation12==2 && relation21==2 )
		return 0;
	else if( relation12==1 )
		return 2;
	else if( relation21==1 )
		return 3;
	else if( relation12==2 && relation21==3 )
		return 4;
	else if( relation12==3 && relation21==2 )
		return 4;
	else		//此时relation12==3 && relation21==3，可能相交和近邻，本算法暂时无法判断
		return -1;	
}


int Cal_Normal_Value_By_Line(double para,double &Val)
{
	static double g_ValTab[33][10] = {
{ 0.5000,  0.5040,  0.5080,  0.5120,  0.5160,  0.5199,  0.5239,  0.5279,  0.5319,  0.5359 },
{ 0.5398,  0.5438,  0.5478,  0.5517,  0.5557,  0.5596,  0.5636,  0.5675,  0.5714,  0.5753 },
{ 0.5793,  0.5832,  0.5817,  0.5910,  0.5948,  0.5987,  0.6026,  0.6064,  0.6103,  0.6141 },
{ 0.6179,  0.6217,  0.6255,  0.6293,  0.6331,  0.6368,  0.6406,  0.6443,  0.6480,  0.6517 },
{ 0.6554,  0.6591,  0.6628,  0.6664,  0.6700,  0.6736,  0.6772,  0.6808,  0.6844,  0.6879 },

{ 0.6915,  0.6950,  0.6985,  0.7019,  0.7054,  0.7088,  0.7123,  0.7157,  0.7190,  0.7224 },
{ 0.7257,  0.7291,  0.7324,  0.7357,  0.7389,  0.7422,  0.7454,  0.7486,  0.7517,  0.7549 },
{ 0.7580,  0.7611,  0.7642,  0.7673,  0.7703,  0.7734,  0.7764,  0.7794,  0.7823,  0.7852 },
{ 0.7881,  0.7910,  0.7939,  0.7967,  0.7995,  0.8023,  0.8051,  0.8078,  0.8106,  0.8133 },
{ 0.8159,  0.8186,  0.8212,  0.8238,  0.8264,  0.8289,  0.8315,  0.8340,  0.8365,  0.8389 },

{ 0.8413,  0.8438,  0.8461,  0.8485,  0.8508,  0.8531,  0.8554,  0.8577,  0.8599,  0.8621 },
{ 0.8643,  0.8665,  0.8686,  0.8708,  0.8729,  0.8749,  0.8770,  0.8790,  0.8810,  0.8830 },
{ 0.8849,  0.8869,  0.8888,  0.8907,  0.8925,  0.8944,  0.8962,  0.8980,  0.8997,  0.9015 },
{ 0.9032,  0.9049,  0.9066,  0.9082,  0.9099,  0.9115,  0.9131,  0.9147,  0.9162,  0.9177 },
{ 0.9192,  0.9207,  0.9222,  0.9236,  0.9251,  0.9265,  0.9279,  0.9292,  0.9306,  0.9319 },

{ 0.9332,  0.9345,  0.9357,  0.9370,  0.9382,  0.9394,  0.9406,  0.9418,  0.9430,  0.9441 },
{ 0.9452,  0.9463,  0.9474,  0.9484,  0.9495,  0.9505,  0.9515,  0.9525,  0.9535,  0.9545 },
{ 0.9554,  0.9564,  0.9573,  0.9582,  0.9591,  0.9599,  0.9608,  0.9616,  0.9625,  0.9633 },
{ 0.9641,  0.9649,  0.9656,  0.9664,  0.9671,  0.9678,  0.9686,  0.9693,  0.9700,  0.9706 },
{ 0.9713,  0.9719,  0.9726,  0.9732,  0.9738,  0.9744,  0.9750,  0.9756,  0.9762,  0.9767 },

{ 0.9772,  0.9778,  0.9783,  0.9788,  0.9793,  0.9798,  0.9803,  0.9808,  0.9812,  0.9817 },
{ 0.9821,  0.9826,  0.9830,  0.9834,  0.9838,  0.9842,  0.9846,  0.9850,  0.9854,  0.9857 },
{ 0.9861,  0.9864,  0.9868,  0.9871,  0.9875,  0.9878,  0.9881,  0.9884,  0.9887,  0.9890 },
{ 0.9893,  0.9896,  0.9898,  0.9901,  0.9904,  0.9906,  0.9909,  0.9911,  0.9913,  0.9916 },
{ 0.9918,  0.9920,  0.9922,  0.9925,  0.9927,  0.9929,  0.9931,  0.9932,  0.9934,  0.9936 },

{ 0.9938,  0.9940,  0.9941,  0.9943,  0.9945,  0.9946,  0.9948,  0.9949,  0.9951,  0.9952 },
{ 0.9953,  0.9955,  0.9956,  0.9957,  0.9959,  0.9960,  0.9961,  0.9962,  0.9963,  0.9964 },
{ 0.9965,  0.9966,  0.9967,  0.9968,  0.9969,  0.9970,  0.9971,  0.9972,  0.9973,  0.9974 },
{ 0.9974,  0.9975,  0.9976,  0.9977,  0.9977,  0.9978,  0.9979,  0.9979,  0.9980,  0.9981 },
{ 0.9981,  0.9982,  0.9983,  0.9983,  0.9984,  0.9984,  0.9985,  0.9985,  0.9986,  0.9986 },

{ 0.9987,  0.9987,  0.9987,  0.9988,  0.9988,  0.9989,  0.9989,  0.9989,  0.9990,  0.9990 },
{ 0.9990,  0.9991,  0.9991,  0.9991,  0.9992,  0.9992,  0.9992,  0.9992,  0.9993,  0.9993 },
{ 0.9993,  0.9993,  0.9994,  0.9994,  0.9994,  0.9994,  0.9994,  0.9995,  0.9995,  0.9995 }
	};

	int bNeg=0;
	int SNRow,SNCol;
	double v1,v2,x1,x2;

	//本函数不能计算超过3.29倍sigma的值
	if( para>3.29 )
	{
		Val = 0.9996;
		return -1;
	}
	if( para<-3.29 )
	{
		Val = 0.0004;
		return -1;
	}

	//对于负数，转化为正数进行计算
	if(para<0)
	{
		para = fabs(para);
		bNeg = 1;
	}

	//求行序号与列序号
	SNRow = (int)(para*10);
	SNCol = (int)(para*100-SNRow*10);

	//处理正好等于3.29的情况，注意浮点数不能直接判断相等
	if( fabs(para-3.29)<0.000001)
	{
		Val = 0.9995;
		if(bNeg==1)
			Val = 1-Val;
		return 1;
	}

	x1 = SNRow / 10.0  + SNCol / 100.0;
	x2 = x1+0.01;

	v1 = g_ValTab[SNRow][SNCol];
	SNCol++;
	if(SNCol==10)
	{
		SNCol = 0;
		SNRow++;
	}
	v2 = g_ValTab[SNRow][SNCol];

	//线性插值
	Val = v1 + (para-x1)/(x2-x1)*(v2-v1);

	if(bNeg==1)
	{
		Val = 1-Val;
	}

	return 1;
}


///////////////////////////////////////////////////////////////////
////////////////  以下是矩阵操作的相关函数/////////////////////////

//矩阵相乘
int MatMul(double * mat1,int h1,int w1,double *mat2,int h2,int w2,double *matOut)
{
	int i,j,k,ho,wo;
	double val;

	double *mat11 = new double[h1*w1];
	double *mat21 = new double[h2*w2];

	memcpy(mat11,mat1,sizeof(double)*h1*w1);
	memcpy(mat21,mat2,sizeof(double)*h2*w2);
	
	//如果前者矩阵宽不等于后者矩阵高，报错返回
	if(w1!=h2)
		return -1;

	ho = h1;
	wo = w2;

	for(i=0;i<ho;i++)
	for(j=0;j<wo;j++)
	{
		val = 0.0;
		for(k=0;k<w1;k++)
			val += mat11[i*w1+k]*mat21[k*w2+j];
		matOut[i*wo+j] = val;
	}

	delete [] mat11;
	delete [] mat21;

	return 1;
}

//判断点是否在线的右侧
int Point_OnRight_Line(double x1,double y1,double x2, double y2,double xo, double yo )
{
	double dir1,dir2;

	dir1 = atan2( y2-y1, x2-x1 );
	dir2 = atan2( yo-y2, xo-x2 );
	dir1 = dir1/qxwPI*180;
	dir2 = dir2/qxwPI*180;

	if( 
		( dir2<dir1 && dir2>dir1-180 ) ||
		( dir2+360<dir1 && dir2+360>dir1-180 ) ||
		( dir2-360<dir1 && dir2-360>dir1-180 )
		)
		return 1;
	else
		return 0;
}


int ConcaveToBulge_iteration( int Num, double *pX, double *pY, vector<qxwPolygon> &vecPoly )
{
	int i,j,SNCon,SNObj;
	int SN,SN1,SN2,SN3,SN4,bRight;
	double dir1,dir2;
	int bConCave;
	int type,Line1Type,Line2Type;
	int bValid,bFind;
	double *pNew1X,*pNew1Y,*pNew2X,*pNew2Y;
	int Num1,Num2;
	int RV;

	/// 调试输出信息用
	static int NumCal=0;
	int Num00;
//	FILE *pF;
	char FN[100]; 
	int bTs = 1;
	double xMin,xMax,yMin,yMax;
	double Wid,tmp;
	int b1,b2;
	
	if(bTs)
	{
		NumCal ++;
		Num00 = NumCal;		
		b1 = b2 = 0;

		/// 统计图形范围，确定显示图形的范围
		xMin = xMax = pX[0];
		yMin = yMax = pY[0];
		for( i=0; i<Num; i++ )
		{
			xMin = min(pX[i],xMin);
			xMax = max(pX[i],xMax);
			yMin = min(pY[i],yMin);
			yMax = max(pY[i],yMax);
		}

		Wid = max(xMax-xMin,yMax-yMin);
		Wid *= 1.2;
		
		tmp = (Wid-(xMax-xMin))/2;
		xMin -= tmp;
		xMax += tmp;
		tmp = (Wid-(yMax-yMin))/2;
		yMin -= tmp;
		yMax += tmp;

		xMin = floor(xMin);
		yMin = floor(yMin);
		xMax = ceil(xMax);
		yMax = ceil(yMax);

		//pF = fopen(FN,"wt+");
		//if (pF == NULL)//20161212 304 静态分析问题更改
		//{
		//	cout << "文件打开失败." << endl;
		//	return 0;
		//}
		///// 输出本身多边形
		//fprintf(pF,"x_%d = [\t",Num00);
		//for( i=0; i<Num; i++ )
		//{
		//	fprintf(pF,"%.0f\t",pX[i]);
		//}
		//fprintf(pF,"%.0f\t];\n",pX[0]);

		//fprintf(pF,"y_%d = [\t",Num00);
		//for( i=0; i<Num; i++ )
		//{
		//	fprintf(pF,"%.0f\t",pY[i]);
		//}
		//fprintf(pF,"%.0f\t];\n",pY[0]);

		//fprintf(pF,"\n");
	}
		
	//首先判断此多边形是否凸多边形，如果不是，找到凹点
	for( bConCave=0,i=0; i<Num; i++ )
	{
		if(i==0)
			dir1 = atan2( pY[i]-pY[Num-1], pX[i]-pX[Num-1] ); 
		else
			dir1 = atan2( pY[i]-pY[i-1], pX[i]-pX[i-1] );
	
		if(i==Num-1)
			dir2 = atan2( pY[0]-pY[i], pX[0]-pX[i] );
		else
			dir2 = atan2( pY[i+1]-pY[i], pX[i+1]-pX[i] );

		dir1 = dir1/qxwPI*180;
		dir2 = dir2/qxwPI*180;

		if( 
			( dir2>dir1 && dir2<dir1+180 ) ||
			( dir2+360>dir1 && dir2+360<dir1+180 ) ||
			( dir2-360>dir1 && dir2-360<dir1+180 )
			)
		{
			bConCave = 1;
			SNCon = i;
			break;
		}
	}

	//如果不是凹多边形，直接将此多边形加入输出多边形列表返回
	if( bConCave==0 )
	{
		qxwPolygon poly;
		qxwPoint2D poi;

		for( i=0; i<Num; i++ )
		{
			poi.x = pX[i];
			poi.y = pY[i];
			
			poly.push_back(poi);
		}
		vecPoly.push_back(poly);

		if(bTs)
		{
			/// 输出自身
			//fprintf(pF,"x_%d_0 = [\t",Num00);
			//for( i=0; i<Num; i++ )
			//{
			//	fprintf(pF,"%.0f\t",pX[i]);
			//}
			//fprintf(pF,"%.0f\t];\n",pX[0]);

			//fprintf(pF,"y_%d_0 = [\t",Num00);
			//for( i=0; i<Num; i++ )
			//{
			//	fprintf(pF,"%.0f\t",pY[i]);
			//}
			//fprintf(pF,"%.0f\t];\n",pY[0]);

			//fprintf(pF,"\n");

			//char buf[1000];
			//sprintf(buf, "figure(%d), plot(x_%d_0,y_%d_0), axis([%f %f %f %f]); \n", Num00,Num00,Num00,xMin,xMax,yMin,yMax);
			//fprintf(pF,"%s",buf);

			//fclose(pF);
		}
		
		return 1;
	}

	//确定当前边的起止序号
	SN1 = SNCon;
	SN2 = SNCon+1;
	if(SNCon==Num-1)
		SN2 = 0;
	
	//分割出一个三角形，作为一个输出多边形
	for( bFind=0,i=0; i<Num; i++ )
	{
		if( i==SN1 || i==SN2 )
			continue;

		bRight = Point_OnRight_Line( pX[SN1], pY[SN1], pX[SN2], pY[SN2], pX[i], pY[i] );
		if( bRight!=1 )
			continue;

		//选定点i作为三角形的一个顶点，SN1、SN2、i构成一个三角形
		//判断多边形的边是否与三角形的两个边相交
		for( bValid=1,j=0; j<Num; j++ )
		{
			SN3 = j;
			SN4 = j+1;
			if(j==Num-1)
				SN4 = 0;

			/// 判断当前边是否与三角形两个边重合
			/// 如果重合不进行判断

			/// 边(SN1,i) 是否与边(SN3,SN4)重合
			if( ( SN1==SN3 && i==SN4 ) || ( SN1==SN4 && i==SN3 ) )
				continue;

			/// 边(SN2,i) 是否与边(SN3,SN4)重合
			if( ( SN2==SN3 && i==SN4 ) || ( SN2==SN4 && i==SN3 ) )
				continue;

			
//			//如果当前边顶点与三角形任一一个顶点重合，不进行计算
//			if( 
//				SN3==SN1 || SN3==SN2 || SN3==i || 
//				SN4==SN1 || SN4==SN2 || SN4==i  
//				)
//				continue;

			/// 判断相交情况，如果相交且交点在线段中间，说明该三角形不是多边形内部
			Cal_TwoLine_Relation( pX[SN1], pY[SN1], pX[i], pY[i], 
				pX[SN3], pY[SN3], pX[SN4], pY[SN4], type, Line1Type, Line2Type, 1e-8 );
			if( type==2 && Line2Type==2 && Line1Type==2 )
			{
				bValid = 0;
				break;
			}

			Cal_TwoLine_Relation( pX[SN2], pY[SN2], pX[i], pY[i], 
				pX[SN3], pY[SN3], pX[SN4], pY[SN4], type, Line1Type, Line2Type, 1e-8 );
			if( type==2 && Line2Type==2 && Line1Type==2 )
			{
				bValid = 0;
				break;
			}
		}

		//如果有效，说明找到了目标三角形
		if(bValid==1)
		{
			bFind = 1;
			SNObj = i;
			break;
		}
	}

	//应该能够找到目标三角形，否则有问题
	if(bFind==0)
		return -1;

	//将目标三角形(SN1、SN2、SBObj构成)加入输出队列
//	pPoly = new QBFPolygon;
//	pPoly->pNext = NULL;
//	pPoly->Num = 3;
//	pPoly->pX = new double[3];
//	pPoly->pY = new double[3];
//	pPoly->pX[0] = pX[SN1];
//	pPoly->pY[0] = pY[SN1];
//	pPoly->pX[1] = pX[SN2];
//	pPoly->pY[1] = pY[SN2];
//	pPoly->pX[2] = pX[SNObj];
//	pPoly->pY[2] = pY[SNObj];
//	pOut->pEnd->pNext = pPoly;
//	pOut->pEnd = pPoly;
//	pOut->Num ++;

	qxwPolygon poly;
	qxwPoint2D poi;

	poi.x = pX[SN1];
	poi.y = pY[SN1];
	poly.push_back(poi);
	poi.x = pX[SN2];
	poi.y = pY[SN2];
	poly.push_back(poi);
	poi.x = pX[SNObj];
	poi.y = pY[SNObj];
	poly.push_back(poi);

	vecPoly.push_back(poly);

	/// 输出中间三角形
	if(bTs)
	{
		//fprintf(pF,"x_%d_0 = [\t",NumCal);

		//fprintf(pF,"%.0f\t",pX[SN1]);
		//fprintf(pF,"%.0f\t",pX[SN2]);
		//fprintf(pF,"%.0f\t",pX[SNObj]);
		//
		//fprintf(pF,"%.0f\t];\n",pX[SN1]);

		//fprintf(pF,"y_%d_0 = [\t",NumCal);

		//fprintf(pF,"%.0f\t",pY[SN1]);
		//fprintf(pF,"%.0f\t",pY[SN2]);
		//fprintf(pF,"%.0f\t",pY[SNObj]);

		//fprintf(pF,"%.0f\t];\n",pY[SN1]);

		//fprintf(pF,"\n");
	}
	
	
	//将剩余的多边形为输入，重复调度本函数，最多可能产生两个新多边形
	
	//处理第一个多边形
	if( SNObj<SN2 )
		Num1 = SNObj+Num-SN2+1;
	else
		Num1 = SNObj-SN2+1;
	if( Num1>2 )
	{
		//生成新的多边形
		pNew1X = new double[Num1];
		pNew1Y = new double[Num1];
		for( SN=0,i=SN2; SN<Num1; SN++,i++ )
		{
			if(i==Num)
				i=0;
			pNew1X[SN] = pX[i];
			pNew1Y[SN] = pY[i];
		}

		if(bTs)
		{
			b1 = 1;

		/*	/// 输出第一个多边形
			fprintf(pF,"x_%d_1 = [\t",Num00);
			for( i=0; i<Num1; i++ )
			{
				fprintf(pF,"%.0f\t",pNew1X[i]);
			}
			fprintf(pF,"%.0f\t];\n",pNew1X[0]);

			fprintf(pF,"y_%d_1 = [\t",Num00);
			for( i=0; i<Num1; i++ )
			{
				fprintf(pF,"%.0f\t",pNew1Y[i]);
			}
			fprintf(pF,"%.0f\t];\n",pNew1Y[0]);

			fprintf(pF,"\n");*/
		}
		

		//递归调用本函数进行多边形拆分
		RV = ConcaveToBulge_iteration( Num1, pNew1X, pNew1Y, vecPoly );
		delete [] pNew1X;
		delete [] pNew1Y;

		if(RV!=1)
			return -1;	
	}

	//处理第二个多边形
	if( SN1<SNObj )
		Num2 = SN1+Num-SNObj+1;
	else
		Num2 = SN1-SNObj+1;
	pNew2X = new double[Num2];
	pNew2Y = new double[Num2];
	for( SN=0,i=SNObj; SN<Num2; SN++,i++ )
	{
		if(i==Num)
			i=0;
		pNew2X[SN] = pX[i];
		pNew2Y[SN] = pY[i];
	}
	
	if(bTs)
	{
		if(Num2>0)
			b2 = 1;
		
		/// 输出第二个多边形
		//fprintf(pF,"x_%d_2 = [\t",Num00);
		//for( i=0; i<Num2; i++ )
		//{
		//	fprintf(pF,"%.0f\t",pNew2X[i]);
		//}
		//fprintf(pF,"%.0f\t];\n",pNew2X[0]);

		//fprintf(pF,"y_%d_2 = [\t",Num00);
		//for( i=0; i<Num2; i++ )
		//{
		//	fprintf(pF,"%.0f\t",pNew2Y[i]);
		//}
		//fprintf(pF,"%.0f\t];\n",pNew2Y[0]);

		//fprintf(pF,"\n");

		/// 绘图命令
		char buf1[100];
		char buf2[100];

		if(b1)
			sprintf(buf1,",x_%d_1,y_%d_1,'r'",Num00,Num00);
		else
			sprintf(buf1," ");

		if(b2)
			sprintf(buf2,",x_%d_2,y_%d_2,'g'",Num00,Num00);
		else
			sprintf(buf2," ");

		char buf[1000];
		sprintf( buf, "figure(%d), plot(x_%d_0, y_%d_0, 'b' %s %s , 'LineWidth', 2 ), axis([%f %f %f %f]); \n ",
			Num00,Num00,Num00,buf1,buf2,xMin,xMax,yMin,yMax);
//		fprintf(pF,"%s",buf);

//		fclose(pF);
	}
	
	
	//递归调用本函数进行多边形拆分
	RV = ConcaveToBulge_iteration( Num2, pNew2X, pNew2Y, vecPoly );
	delete [] pNew2X;
	delete [] pNew2Y;

	if(RV!=1)
		return -1;
	
	return 1;
}

int ConcaveToBulge_iteration( int Num, double *pX, double *pY, QBFPolygonList *pOut )
{
	int i,j,SNCon,SNObj;
	int SN,SN1,SN2,SN3,SN4,bRight;
	double dir1,dir2;
	int bConCave;
	int type,Line1Type,Line2Type;
	int bValid,bFind;
	QBFPolygon *pPoly;
	double *pNew1X,*pNew1Y,*pNew2X,*pNew2Y;
	int Num1,Num2;
	int RV;
	
	//首先判断此多边形是否凸多边形，如果不是，找到凹点
	for( bConCave=0,i=0; i<Num; i++ )
	{
		if(i==0)
			dir1 = atan2( pY[i]-pY[Num-1], pX[i]-pX[Num-1] ); 
		else
			dir1 = atan2( pY[i]-pY[i-1], pX[i]-pX[i-1] );
	
		if(i==Num-1)
			dir2 = atan2( pY[0]-pY[i], pX[0]-pX[i] );
		else
			dir2 = atan2( pY[i+1]-pY[i], pX[i+1]-pX[i] );

		dir1 = dir1/qxwPI*180;
		dir2 = dir2/qxwPI*180;

		if( 
			( dir2>dir1 && dir2<dir1+180 ) ||
			( dir2+360>dir1 && dir2+360<dir1+180 ) ||
			( dir2-360>dir1 && dir2-360<dir1+180 )
			)
		{
			bConCave = 1;
			SNCon = i;
			break;
		}
	}

	//如果不是凹多边形，直接将此多边形加入输出多边形列表返回
	if( bConCave==0 )
	{
		pPoly = new QBFPolygon;
		pPoly->pNext = NULL;
		pPoly->Num = Num;
		pPoly->pX = new double[Num];
		pPoly->pY = new double[Num];
		memcpy( pPoly->pX, pX, sizeof(double)*Num );
		memcpy( pPoly->pY, pY, sizeof(double)*Num );
		pOut->pEnd->pNext = pPoly;
		pOut->pEnd = pPoly;
		pOut->Num ++;
		return 1;
	}

	//确定当前边的起止序号
	SN1 = SNCon;
	SN2 = SNCon+1;
	if(SNCon==Num-1)
		SN2 = 0;
	
	//分割出一个三角形，作为一个输出多边形
	for( bFind=0,i=0; i<Num; i++ )
	{
		if( i==SN1 || i==SN2 )
			continue;

		bRight = Point_OnRight_Line( pX[SN1], pY[SN1], pX[SN2], pY[SN2], pX[i], pY[i] );
		if( bRight!=1 )
			continue;

		//选定点i作为三角形的一个顶点，SN1、SN2、i构成一个三角形
		//判断多边形的边是否与三角形的两个边相交
		for( bValid=1,j=0; j<Num; j++ )
		{
			SN3 = j;
			SN4 = j+1;
			if(j==Num-1)
				SN4 = 0;

			//如果当前边顶点与三角形任一一个顶点重合，不进行计算
			if( 
				SN3==SN1 || SN3==SN2 || SN3==i || 
				SN4==SN1 || SN4==SN2 || SN4==i  
				)
				continue;

			//判断相交情况
			Cal_TwoLine_Relation( pX[SN1], pY[SN1], pX[i], pY[i], 
				pX[SN3], pY[SN3], pX[SN4], pY[SN4], type, Line1Type, Line2Type, 1e-8 );
			if( type==2 && Line2Type==2 && Line1Type==2 )
			{
				bValid = 0;
				break;
			}

			Cal_TwoLine_Relation( pX[SN2], pY[SN2], pX[i], pY[i], 
				pX[SN3], pY[SN3], pX[SN4], pY[SN4], type, Line1Type, Line2Type, 1e-8 );
			if( type==2 && Line2Type==2 && Line1Type==2 )
			{
				bValid = 0;
				break;
			}
		}

		//如果有效，说明找到了目标三角形
		if(bValid==1)
		{
			bFind = 1;
			SNObj = i;
			break;
		}
	}

	//应该能够找到目标三角形，否则有问题
	if(bFind==0)
		return -1;

	//将目标三角形(SN1、SN2、SBObj构成)加入输出队列
	pPoly = new QBFPolygon;
	pPoly->pNext = NULL;
	pPoly->Num = 3;
	pPoly->pX = new double[3];
	pPoly->pY = new double[3];
	pPoly->pX[0] = pX[SN1];
	pPoly->pY[0] = pY[SN1];
	pPoly->pX[1] = pX[SN2];
	pPoly->pY[1] = pY[SN2];
	pPoly->pX[2] = pX[SNObj];
	pPoly->pY[2] = pY[SNObj];
	pOut->pEnd->pNext = pPoly;
	pOut->pEnd = pPoly;
	pOut->Num ++;
	
	//将剩余的多边形为输入，重复调度本函数，最多可能产生两个新多边形
	
	//处理第一个多边形
	if( SNObj<SN2 )
		Num1 = SNObj+Num-SN2+1;
	else
		Num1 = SNObj-SN2+1;
	if( Num1>2 )
	{
		//生成新的多边形
		pNew1X = new double[Num1];
		pNew1Y = new double[Num1];
		for( SN=0,i=SN2; SN<Num1; SN++,i++ )
		{
			if(i==Num)
				i=0;
			pNew1X[SN] = pX[i];
			pNew1Y[SN] = pY[i];
		}

		//递归调用本函数进行多边形拆分
		RV = ConcaveToBulge_iteration( Num1, pNew1X, pNew1Y, pOut );
		delete [] pNew1X;
		delete [] pNew1Y;

		if(RV!=1)
			return -1;		
	}

	//处理第二个多边形
	if( SN1<SNObj )
		Num2 = SN1+Num-SNObj+1;
	else
		Num2 = SN1-SNObj+1;
	pNew2X = new double[Num2];
	pNew2Y = new double[Num2];
	for( SN=0,i=SNObj; SN<Num2; SN++,i++ )
	{
		if(i==Num)
			i=0;
		pNew2X[SN] = pX[i];
		pNew2Y[SN] = pY[i];
	}
	
	//递归调用本函数进行多边形拆分
	RV = ConcaveToBulge_iteration( Num2, pNew2X, pNew2Y, pOut );
	delete [] pNew2X;
	delete [] pNew2Y;

	if(RV!=1)
		return -1;
	
	return 1;
}

/****************************************************
功能：
	输入任意多边形，将其拆分为数个凸多边形。

输入：
    文件名，文件中以文本存放多边形角点个数并按顺时针顺序存放
每个角点的x、y坐标。

输出：
    将拆分结果写入要求的文本文件。
	依次存放多边形个数n
	多边形1 的顶点个数，顺时针存放顶点x、y坐标
	多边形2 的顶点个数，顺时针存放顶点x、y坐标
	。。。。。。
	多边形n 的顶点个数，顺时针存放顶点x、y坐标

返回值：
	1，成功
	－1，失败
***********************************************************/
int ConcaveToBulgePolygon( char *FN1, char *FN2 )
{
	QBFPolygonList ListPoly;
	QBFPolygon *pHead,*pCur;
	double *pX,*pY;
	int Num,i,j,RV;
	FILE * pF;

	//初始化多边形链表
	pHead = new QBFPolygon;
	pHead->Num = 0;
	pHead->pNext = NULL;
	pHead->pX = NULL;
	pHead->pY = NULL;
	ListPoly.Num = 1;
	ListPoly.pHead = ListPoly.pEnd = pHead;
	
	//打开输入文件
	pF = fopen( FN1, "rt" );
	if( pF==NULL )
	{
		while(pHead)//20161212 304 静态分析问题更改
		{
			pCur = pHead->pNext;
			if(pHead->Num>0)
			{
				delete [] pHead->pX;
				delete [] pHead->pY;
			}
			delete pHead;
			pHead = pCur;
		}
		return -1;
	}

	fscanf( pF, "%d", &Num );	//读入角点个数

	//分配空间
	pX = new double[Num];
	pY = new double[Num];

	//依次读入角点坐标
	for( i=0; i<Num; i++ )
	{
		fscanf( pF, "%lf", &(pX[i]) );
		fscanf( pF, "%lf", &(pY[i]) );
	}
	fclose(pF);

	//调用递归函数拆分多边形
	RV = ConcaveToBulge_iteration( Num, pX, pY, &ListPoly );
	if(RV!=1)
	{
		//释放空间
		delete [] pX;
		delete [] pY;
		ListPoly.pHead = ListPoly.pEnd = NULL;
		ListPoly.Num = 0;
		while(pHead)
		{
			pCur = pHead->pNext;
			if(pHead->Num>0)
			{
				delete [] pHead->pX;
				delete [] pHead->pY;
			}
			delete pHead;
			pHead = pCur;
		}
		return -1;
	}

	pCur = ListPoly.pHead->pNext;
	pF = fopen( FN2, "wt+" );
	if(pF == NULL) //20161225 304pingce
	{ 
		delete [] pX;
		delete [] pY;
		return -1;
	}
	fprintf( pF, "%d\n\n", ListPoly.Num-1 );
	for( i=0; i<ListPoly.Num-1; i++ )
	{
		fprintf( pF, "%d\n", pCur->Num );
		for( j=0; j<pCur->Num; j++ )
			fprintf( pF, "%f\t%f\n", pCur->pX[j], pCur->pY[j] );
		fprintf(pF,"\n");
		pCur = pCur->pNext;
	}
	fclose(pF);

	//释放空间
	delete [] pX;
	delete [] pY;
	ListPoly.pHead = NULL;
	ListPoly.Num = 0;
	while(pHead)
	{
		pCur = pHead->pNext;
		if(pHead->Num>0)
		{
			delete [] pHead->pX;
			delete [] pHead->pY;
		}
		delete pHead;
		pHead = pCur;
	}

	return 1;
}


//测试相关函数

/*************************************************
 测试凹多边形拆分为凸多边形递归函数的函数 

从文本文件读入多边形角点个数以及每个角点的x、y坐标，然后进行拆分，并写入文件
*******************************************************************/
int ts_ConcaveToBulge_iteration()
{
	FILE *pF;
	double pX[100],pY[100];
	int Num,i,j;
	QBFPolygonList ListPoly;
	QBFPolygon *pHead,*pCur;
	int RV;

	RV =  ConcaveToBulgePolygon( "c:\\JD.txt", "c:\\res.txt" );
	return 1;

	//pHead = new QBFPolygon;
	//pHead->Num = 0;
	//pHead->pNext = NULL;
	//pHead->pX = NULL;
	//pHead->pY = NULL;
	//ListPoly.Num=1;
	//ListPoly.pHead = ListPoly.pEnd = pHead;
	//
	//pF = fopen("c:\\JD.txt","rt");
	//fscanf( pF, "%d", &Num );
	//for( i=0; i<Num; i++ )
	//{
	//	fscanf( pF, "%lf", &(pX[i]) );
	//	fscanf( pF, "%lf", &(pY[i]) );
	//}
	//fclose(pF);

	//RV = ConcaveToBulge_iteration(Num,pX,pY,&ListPoly);

	//pCur = ListPoly.pHead->pNext;
	//pF = fopen("c:\\res.txt","wt+");
	//fprintf( pF, "%d\n\n", ListPoly.Num-1 );
	//for( i=0; i<ListPoly.Num-1; i++ )
	//{
	//	for( j=0; j<pCur->Num; j++ )
	//		fprintf( pF, "%f %f\n", pCur->pX[j], pCur->pY[j] );
	//	fprintf(pF,"\n");
	//	pCur = pCur->pNext;
	//}
	//fclose(pF);

	////释放空间
	//ListPoly.pHead = NULL;
	//ListPoly.Num = 0;
	//while(pHead)
	//{
	//	pCur = pHead->pNext;
	//	if(pHead->Num>0)
	//	{
	//		delete [] pHead->pX;
	//		delete [] pHead->pY;
	//	}
	//	delete pHead;
	//	pHead = pCur;
	//}

	//return 1;
}
