// qxwBaseFunc.h

#include <math.h>
#include <direct.h>
#include <STRING.H>
#include <Time.h>
#include <vector>
#include <fstream>
using namespace std;

using std::vector;

#ifndef __qxwBaseFunc__
#define __qxwBaseFunc__

#ifndef NULL
#define NULL 0
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

//const double qxwPI = 3.1415926535897932384626433832795;
//const double qxwLimMinValue = 1e-30;
//const double qxwPI_DIV_2 = qxwPI/2.0;
//const double qxwPI_DIV_4 = qxwPI/4.0;
//const double qxwPI_DIV_3 = qxwPI/3.0;
//const double qxwPI_DIV_6 = qxwPI/6.0;
//const double qxwPI_MUL_2 = qxwPI*2.0;
//const double qxwPI_MUL_4 = qxwPI*4.0;

//const double DoubleMax = 1e+308;

//判断一个数是否为偶数
//不能适用于浮点数
template <class PT>
bool IsEven(PT data)
{
	if( data%2==0 )
		return TRUE;
	else
		return FALSE;
}

//分配2维内存空间
void **NewSpace2D( int Row, int Col, short Size );

//释放2维内存空间
template <class PT>
void DeleteSpace2D( PT ** (&pData), int Row=0 )
{
	if( pData==NULL ) return;
	if( pData[0]!=NULL ) delete [] pData[0];
	delete [] pData;
	pData = NULL;
}

//分配3维内存空间
void ***NewSpace3D( int Len, int Row, int Col, short Size );

//释放3维内存空间
template <class PT>
void DeleteSpace3D( PT *** (&pData), int Len=0, int Row=0 )
{
	if( pData==NULL ) 
		return;
	if( pData[0][0]!=NULL ) 
		delete [] pData[0][0];
	DeleteSpace2D(pData,Len);
	pData =NULL;
}

//分配4维内存空间
void ****NewSpace4D( int Len1, int Len2, int Row, int Col, short Size );

//释放4维内存空间
template <class PT>
void DeleteSpace4D(PT **** (&pData),int Len1,int Len2,int Row=0)
{
	if(pData==NULL)
		return;
	if(pData[0][0][0]!=NULL) 
		delete [] pData[0][0][0];
	DeleteSpace3D(pData,Len1,Len2);
	pData =NULL;
}

/*****************  函数：ArrayOverturn()  ***********************
功能：一维数组翻转，目的指针和源指针可以指向同一个空间
*****************************************************************/
template <class PT>
int ArrayOverturn(PT *pOrData,PT *pDesData,int Len)
{
	PT *pTemp;
	int i;

	if( pOrData==NULL || pDesData==NULL )
		return -1;
	pTemp = new PT[Len];
	if(pTemp==NULL)
		return -1;
	memcpy(pTemp,pOrData,sizeof(PT)*Len);
	for( i=0; i<Len; i++)
		pDesData[i] = pTemp[Len-1-i];
	delete [] pTemp;
	return 1;
}

/*****************  函数：ArrayOverturn_2D()  ********************
功能：二维数组翻转，UpDown表示上下翻转，LeftRight表示左右翻转，目
的指针和源指针可以指向同一空间
*****************************************************************/
template <class PT>
int ArrayOverturn2D(PT **pOrData,PT **pDesData,int Hei,int Wid,
					 bool UpDown,bool LeftRight)
{
	PT *pTemp=NULL;
	int i;

	//如果内存空间异常，什么都不做返回
	if( pOrData==NULL || pDesData==NULL )
		return -1;
	for( i=0; i<Hei; i++)
	{
		if( pOrData[i]==NULL || pDesData[i]==NULL )
			return -1;
	}

	//申请临时内存空间
	pTemp = new PT[Wid];
	if(pTemp==NULL)	return -1;
	
	//上下左右都不翻
	if((!UpDown)&&(!LeftRight))
	{
		for(i=0;i<Hei;i++)
		{
			memcpy(pTemp,pOrData[i],sizeof(PT)*Wid);
			memcpy(pDesData[i],pTemp,sizeof(PT)*Wid);
		}
		delete [] pTemp;
		return 1;
	}

	//开始翻转
	if(UpDown)
	{
		for( i=0; i<Hei/2; i++)
		{
			memcpy(pTemp,pOrData[i],sizeof(PT)*Wid);
			memcpy(pDesData[i],pOrData[Hei-1-i],sizeof(PT)*Wid);
			memcpy(pDesData[Hei-1-i],pTemp,sizeof(PT)*Wid);
		}
		if(LeftRight)
		{
			for( i=0; i<Hei; i++)
				ArrayOverturn(pDesData[i],pDesData[i],Wid);
		}
	}
	else
	{
		for( i=0; i<Hei; i++)
			ArrayOverturn(pOrData[i],pDesData[i],Wid);
	}
			
	delete [] pTemp;
	return 1;
}

/************** 函数：EnlargeArray2D_Near()  ******************
功能：
    数组插值放大，最近插值
**************************************************************/
template <class PT>
int EnlargeArray2D_Near(PT ** pSource,int SmallHei,int SmallWid,
						PT ** pDest,int BigHei,int BigWid)
{
	int i,j;
	double te1,te2;

	//检查像素尺寸是否正常
	if( SmallHei<1 || SmallWid<1 || BigHei<2 || BigWid<2 )
		return -2;
	//检查内存空间是否正常
	if( pSource==NULL || pDest==NULL )
		return -1;
	for( i=0; i<SmallHei; i++)
	{
		if( pSource[i]==NULL )
			return -1;
	}
	for( i=0; i<BigHei; i++)
	{
		if( pDest[i]==NULL )
			return -1;
	}

	te1 = (SmallHei-1.0)/(BigHei-1);
	te2 = (SmallWid-1.0)/(BigWid-1);
	for( i=0; i<BigHei; i++ )
	for( j=0; j<BigWid; j++ )
		pDest[i][j]=pSource[(int)(i*te1+0.5)][(int)(j*te2+0.5)];
	return 1;
}

//计算一副图象的相关长度
//注意相关长度到底是半长度还是全长度
template <class PT>
int CorLen(PT **pDataIn,int Hei,int Wid,double th,int &Lv,int &Lh)
{
	double Var,sumSqGray,sumGray,EGray,Num,Co;
	int i,j,l,smaHei,smaWid;
	double *pLv,*pLh;

	pLv = new double[Hei];
	pLh = new double[Wid];
	//计算方差
	sumSqGray = sumGray = 0.0;
	Num = Hei*Wid-1;
	for(i=0;i<Hei;i++)
	for(j=0;j<Wid;j++)
	{
		sumGray += pDataIn[i][j];
		sumSqGray += 1.0*pDataIn[i][j]*pDataIn[i][j];
	}
	EGray = sumGray/(Num+1);
	Var = sumSqGray/Num-(sumGray/Num)*(sumGray/Num);

	//计算横向自相关序列
	for(l=0;l<Wid;l++)
	{
		smaWid = Wid-l;
		Co = 0.0;
		for(i=0;i<Hei;i++)
		for(j=0;j<smaWid;j++)
			Co += (pDataIn[i][j]-EGray)*(pDataIn[i][j+l]-EGray);
		pLh[l] = Co/Hei/smaWid/Var;
		if(pLh[l]<=th)
		{
			Lh = l;
			break;
		}
	}
	Lh = Lh+1;

	//计算纵向相关长度
	for(l=0;l<Hei;l++)
	{
		smaHei = Hei-l;
		Co = 0.0;
		for(i=0;i<smaHei;i++)
		for(j=0;j<Wid;j++)
			Co += (pDataIn[i][j]-EGray)*(pDataIn[i+l][j]-EGray);
		pLv[l] = Co/Wid/smaHei/Var;
		if(pLv[l]<=th)
		{
			Lv = l;
			break;
		}
	}
	Lv = Lv+1;
	
	delete [] pLh;
	delete [] pLv;
	return 1;
}

//删除特定文件夹中所有的文件，不删除目录
int RMFilesFromDir(const char * DirName);

/***************************************
功能 ：删除目录以及目录中的文件和子目录。
算法：递归调用
***************************************/
int RMDirAndFiles(const char* pDirName);

int PointToLineEquation(double x1,double y1,double x2,double y2,
						double &para_a,double &para_b,double &para_c);

int PointToLineEquation(int x1,int y1,int x2,int y2,
						double &para_a,double &para_b,double &para_c);


//给出两点坐标，计算线段偏移坐标序列
int PointToLine(int x1,int y1,int x2,int y2,int &NumPoi,
				int * (&pX),int *(&pY));

int PointToLine( double x1,double y1,double x2,double y2,
				double xo,double yo,double cellsize,
				int &NumPoi,int * (&pX),int * (&pY) );

int TwoPoint2Angle(double x1,double y1,double x2,double y2,int bHudu,double &AngleRes);

/**** 求射线和线段的交点数目
给出射线端点坐标和方向角(正东为零，逆时针为正)，给出线段端点坐标
结果信息Res说明：
	1，射线与直线不重合，且与线段不相交
	2，射线与直线不重合，且与线段相交于普通点
	3，射线与直线不重合，且与线段相交于线段端点
	4，射线与直线重合，但线段不在射线上
	5，射线与直线重合，且线段在射线上
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
					   );


//数字字符串合法性检查函数
int StringIsNumber(const char * pObj);
int DiffOfTwoDir(double dir1,double dir2,bool bDu,double &res);

int CreateCurTimeString(char * pObjStr);
int CreateCurTimeString(time_t &curTime, char * pObjStr);

int InterLinear1D( double x1, double v1 ,double x2, double v2, double x, double & v );

int InterLinear2D( double x1, double x2, double y1, double y2, double v11, double v12, 
				  double v21, double v22, double x, double y, double &v );

int InterLinear3D( double x1, double x2, double y1, double y2, double z1, double z2, 
				  double v111, double v121, double v211, double v221, 
				  double v112, double v122, double v212, double v222, 
				  double x, double y, double z, double &v );

long GetFileLen(char *pFN);
int CombinFile(char *pSouFN,char *pObjFN,int Num);
int SplitFile(char *pSouFN,char *pObjFN);

int Cal_PoiLine_Relation( double poiX, double poiY, double x1, double y1, 
					   double x2, double y2, double minV );
int Cal_RadialLine_Relation( double r_x, double r_y, double r_dir, 
							double x1, double y1, double x2, double y2, 
							int &Type, int &RadType,int &LineType, double minV=1e-5 );
int Cal_PoiPolygon_Relation( int NumPoi, double *pX, double *pY, 
						  double PoiX, double PoiY, double minV=1e-5 );

int Cal_TwoPolygon_Relation( int NP1, double *pX1, double *pY1,
							int NP2, double *pX2, double *pY2, double minV=1e-5 );

int Cal_TwoLine_Relation( double x1, double y1, double x2, double y2,
						 double x3, double y3, double x4, double y4, 
						 int &Type, int &Line1Type, int &Line2Type, double minV=1e-5 );

int Cal_TwoLine_Cross( double a1, double b1, double c1, double a2, 
					  double b2, double c2, double &x, double &y );

int Cal_Normal_Value_By_Line(double para,double &Val);


int MatMul(double * mat1,int h1,int w1,double *mat2,int h2,int w2,double *matOut);

//数据结构

//多边形
typedef struct tagQBFPolygon{
	struct tagQBFPolygon *pNext;
	int Num;
	double *pX;
	double *pY;		
}QBFPolygon;

//多边形链表
typedef struct tagQBFPolygonList{
	struct tagQBFPolygon * pHead;
	struct tagQBFPolygon * pEnd;
	int Num;
}QBFPolygonList;

typedef struct tagQxwPoint2D{
	double x;
	double y;
}qxwPoint2D;

typedef vector<qxwPoint2D> qxwPolygon;


int ConcaveToBulgePolygon( char *FN1, char *FN2 );
int ConcaveToBulge_iteration( int Num, double *pX, double *pY, QBFPolygonList *pOut );
int ConcaveToBulge_iteration( int Num, double *pX, double *pY, vector<qxwPolygon> &vecPoly );
int Point_OnRight_Line(double x1,double y1,double x2, double y2,double xo, double yo );

//测试相关函数
int ts_ConcaveToBulge_iteration();

#endif
