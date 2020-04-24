#pragma once

#include <vector>

using namespace std;

#ifndef NULL
#define NULL 0
#endif

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#ifndef NULL
#define NULL            0
#endif

const double qxwLimMinValue = 1e-30;
const double qxwPI = 3.1415926535897932384626433832795;
struct  szypolygen_GS//szy20141023 add
{
	double p_x;
	double p_y;
	double p_z;
};
struct  szypolygen_DD//szy20141023 add
{
	double DD_x;
	double DD_y;
	double DD_z;
};
//�������Ը��ô���

const int ERR_NEWSPACE = -1;			//�����ڴ����
const int ERR_CANNOTPROC = -2;			//���������������������޷�����
const int ERR_TOONEAR = -3;				//����̫�����޷�����

//const double qxwPI = 3.1415926535897932384626433832795;
//const double qxwLimMinValue = 1e-30;
const double qxwPI_DIV_2 = qxwPI/2.0;
const double qxwPI_DIV_4 = qxwPI/4.0;
const double qxwPI_DIV_3 = qxwPI/3.0;
const double qxwPI_DIV_6 = qxwPI/6.0;
const double qxwPI_MUL_2 = qxwPI*2.0;
const double qxwPI_MUL_4 = qxwPI*4.0;
const double DoubleMax = 1e+308;

class MpGlobalFunc
{
public:
	MpGlobalFunc(void);
	~MpGlobalFunc(void);

	/// ���Բ�ֵ����
	///
	/// ��֪(a0,b0)��(a1,b1)��a������a��Ӧ��ֵ
	static double linear_interpolation(double a0,double b0,double a1,double b1,double a);

	//����2ά�ڴ�ռ�
static  void **NewSpace2D( int Row, int Col, short Size );
static  void **NewSpace2D_1( int Row, int Col, short Size );

static  void *** NewSpace3D(int Len,int Row,int Col,short Size);

//�ͷ�2ά�ڴ�ռ�
template <class PT>
static void DeleteSpace2D( PT ** (&pData), int Row=0 )
{
	if( pData==NULL ) return;
	if( pData[0]!=NULL ) delete [] pData[0];
	delete [] pData;
	pData = NULL;
}

//�ͷ�3ά�ڴ�ռ�
template <class PT>
static void DeleteSpace3D( PT *** (&pData), int Len=0, int Row=0 )
{
	if( pData==NULL ) 
		return;
	if( pData[0][0]!=NULL ) 
		delete [] pData[0][0];
	DeleteSpace2D(pData,Len);
	pData =NULL;
}

template <class PT>
static void DeleteSpace2D_1( PT ** (&pData), int Row )
{
	int i;

	if( pData==NULL ) return;

	for( i=0; i<Row; i++ )
	{
		if( pData[i] != NULL )
			delete [] pData[i];
		pData[i] = NULL;
	}
	
	delete [] pData;
	pData = NULL;
}

template <class PT>
static void DeleteSpace1D( PT * (&pData) )
{
	if(pData==NULL)
		return;
	
	delete [] pData;
	pData = NULL;
}

//���㻭�ߣ��������˵㸡�����꣬���½ǵ�ĸ�������ͷֱ��ʣ�����������߶ε�����
static int PointToLine( double x1,double y1,double x2,double y2,
				double xo,double yo,double cellsize,
				int &NumPoi,int * (&pX),int * (&pY) );

//���㻭�ߣ��������˵㸡�����꣬���½ǵ�ĸ����������������ķֱ��ʣ�����������߶ε�����
static int PointToLine( double x1,double y1,double x2,double y2,
				double xo,double yo,double csX,double csY,
				int &NumPoi,int * (&pX),int * (&pY) );

//���㻭�ߣ��������˵㸡�����꣬���½ǵ�ĸ����������������ķֱ��ʣ�����������߶ε�����
// ע�⣬���㻭���޶���x��y����ķ�Χ������������ְ�������Сֵȡ
static int PointToLineWithRgn( double x1,double y1,double x2,double y2,
				double xo,double yo,double csX,double csY, int minX,int maxX,int minY,int maxY,
				int &NumPoi,int * (&pX),int * (&pY) );

//���ݶ˵��������ֱ�߱�׼����
// para_a * x + para_b * y + c = 0
//int PointToLineEquation(double x1,double y1,double x2,double y2,
//						double &para_a,double &para_b,double &para_c);

//���������Ƕȣ���ǶȲ��λ��
static  int DiffOfTwoDir(double dir1,double dir2,bool bDu,double &res);

static  int InterLinear1D( double x1, double v1 ,double x2, double v2, double x, double & v );

static  int InterLinear2D( double x1, double x2, double y1, double y2, double v11, double v12, 
				  double v21, double v22, double x, double y, double &v );






static  int PointToLineEquation(double x1,double y1,double x2,double y2,double &para_a,
						double &para_b,double &para_c,double minV=1e-5);

 static int Cal_RadialLine_Relation( double r_x, double r_y, double r_dir, 
							double x1, double y1, double x2, double y2, 
							int &Type, int &RadType,int &LineType, double minV=1e-5);

static  int Cal_TwoLine_Cross( double x1, double y1, double x2, double y2, 
					  double x3, double y3, double x4, double y4, 
					  double &x, double &y, double minV=1e-5 );
static  int Cal_TwoLine_Cross( double a1, double b1, double c1, double a2, 
					  double b2, double c2, double &x, double &y );

//�����߶�
static  int Cal_PoiLine_Relation( double poiX, double poiY, double x1, double y1, 
					   double x2, double y2, double minV=1e-5 );

//�߶����߶�
static  int Cal_TwoLine_Relation( double x1, double y1, double x2, double y2,
						 double x3, double y3, double x4, double y4, 
						 int &Type, int &Line1Type, int &Line2Type, double minV=1e-5 );

//��������
static  int Cal_PoiPolygon_Relation( int NumPoi, double *pX, double *pY, 
						  double PoiX, double PoiY, double minV=1e-5 );

//�߶�������
static  int Cal_LinePolygon_Relation( int NP, double *pX, double *pY, double X1, double Y1, 
							 double X2, double Y2, double minV=1e-5 );

//�����������
static  int Cal_TwoPolygon_Relation( int NP1, double *pX1, double *pY1,
							int NP2, double *pX2, double *pY2, double minV=1e-5 );

static  unsigned int GetFileLen(const char *pFN);

static  int ReadFileToBuf(const char * FN,char *(&pDat),unsigned int &Len );

static  int CopyFileToFile(const char *FNOld, const char *FNNew, int bOverWite );

static  int stringMakeUp(char * pDat);

static int SplitString( CString Src, std::vector<CString> Splitter, std::vector<CString> &vecStrings );

static int JieCheng(int a);

static int zuhe(int NumALl,int NumSel,int &NumRes, int **(&pRes) ); 

static int Chk2RegularRect(double xl1,double yb1,double xr1,double yt1,
						   double xl2,double yb2,double xr2,double yt2);

static int SuperCopyDir(const char *DirS,const char *DirO,int NumTry=500,int InterTime=1000);

static int RMDirAndFiles(LPCSTR pDirName);

static int RMObjFileInPath(const char * Path,const char *ObjStr);

static int GetObjFileNumInPath(const char *Path,const char *ObjStr,int &Num);

static int GetObjFileNameInPath(const char *Path,const char *ObjStr,CString &FN,int SN=1);

static int IsLineOutCircle(double x1,double y1,double x2,double y2,
						   double xo,double yo,double Rad,double minV=1e-5);

static int StringIsInt(const char * pObj);

static int StringIsNumber(const char * pObj);

static int Cal_LineCircle_Cross(double pa,double pb,double pc,double xo,double yo,double R,
								double &xc1,double &yc1,double &xc2,double &yc2,double minV=1e-5);

static int StrDelSpace(char * pObj);


static unsigned char ** LineToPic8bit(int Num,double *pX,double *pY,double cs,int &Hei,int &Wid);
static unsigned char ** LineToPic24bit(int Num,double *pX,double *pY,double cs,int &Hei,int &Wid);

/// ���Ŀ�������Ƿ񱻾������鸲�ǣ�Ŀ���������귶ΧΪxl yb xr yt 
/// ������������º����Ͻǵ�����ͨ��ָ�����������դ�񻯵ļ�鷽�����ֱ���Ϊcs
/// ����1�����ǣ�0��������
static int ChkRegionCover(double xl,double yb,double xr,double yt,int NumReg,
						  double *pXl,double *pYb,double *pXr,double *pYt,double cs);

// �ƶ�һ�����ڵ���һ�����ڵ��ض�λ��
static int PositionWndAroundAnother(CWnd & WndStill, CWnd & Wnd2Move, int iPosFlag);
//===�������λ�á��۲��λ�á��ӳ��Ƕȡ���ȡ�����ɵ����������С
static int DL_calc_JXQY(vector<szypolygen_GS> eye, vector<szypolygen_GS> view, double fov_x, double fov_y, vector<szypolygen_GS> &result);
};

