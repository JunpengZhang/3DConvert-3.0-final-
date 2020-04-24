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

//�ж�һ�����Ƿ�Ϊż��
//���������ڸ�����
template <class PT>
bool IsEven(PT data)
{
	if( data%2==0 )
		return TRUE;
	else
		return FALSE;
}

//����2ά�ڴ�ռ�
void **NewSpace2D( int Row, int Col, short Size );

//�ͷ�2ά�ڴ�ռ�
template <class PT>
void DeleteSpace2D( PT ** (&pData), int Row=0 )
{
	if( pData==NULL ) return;
	if( pData[0]!=NULL ) delete [] pData[0];
	delete [] pData;
	pData = NULL;
}

//����3ά�ڴ�ռ�
void ***NewSpace3D( int Len, int Row, int Col, short Size );

//�ͷ�3ά�ڴ�ռ�
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

//����4ά�ڴ�ռ�
void ****NewSpace4D( int Len1, int Len2, int Row, int Col, short Size );

//�ͷ�4ά�ڴ�ռ�
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

/*****************  ������ArrayOverturn()  ***********************
���ܣ�һά���鷭ת��Ŀ��ָ���Դָ�����ָ��ͬһ���ռ�
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

/*****************  ������ArrayOverturn_2D()  ********************
���ܣ���ά���鷭ת��UpDown��ʾ���·�ת��LeftRight��ʾ���ҷ�ת��Ŀ
��ָ���Դָ�����ָ��ͬһ�ռ�
*****************************************************************/
template <class PT>
int ArrayOverturn2D(PT **pOrData,PT **pDesData,int Hei,int Wid,
					 bool UpDown,bool LeftRight)
{
	PT *pTemp=NULL;
	int i;

	//����ڴ�ռ��쳣��ʲô����������
	if( pOrData==NULL || pDesData==NULL )
		return -1;
	for( i=0; i<Hei; i++)
	{
		if( pOrData[i]==NULL || pDesData[i]==NULL )
			return -1;
	}

	//������ʱ�ڴ�ռ�
	pTemp = new PT[Wid];
	if(pTemp==NULL)	return -1;
	
	//�������Ҷ�����
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

	//��ʼ��ת
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

/************** ������EnlargeArray2D_Near()  ******************
���ܣ�
    �����ֵ�Ŵ������ֵ
**************************************************************/
template <class PT>
int EnlargeArray2D_Near(PT ** pSource,int SmallHei,int SmallWid,
						PT ** pDest,int BigHei,int BigWid)
{
	int i,j;
	double te1,te2;

	//������سߴ��Ƿ�����
	if( SmallHei<1 || SmallWid<1 || BigHei<2 || BigWid<2 )
		return -2;
	//����ڴ�ռ��Ƿ�����
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

//����һ��ͼ�����س���
//ע����س��ȵ����ǰ볤�Ȼ���ȫ����
template <class PT>
int CorLen(PT **pDataIn,int Hei,int Wid,double th,int &Lv,int &Lh)
{
	double Var,sumSqGray,sumGray,EGray,Num,Co;
	int i,j,l,smaHei,smaWid;
	double *pLv,*pLh;

	pLv = new double[Hei];
	pLh = new double[Wid];
	//���㷽��
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

	//����������������
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

	//����������س���
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

//ɾ���ض��ļ��������е��ļ�����ɾ��Ŀ¼
int RMFilesFromDir(const char * DirName);

/***************************************
���� ��ɾ��Ŀ¼�Լ�Ŀ¼�е��ļ�����Ŀ¼��
�㷨���ݹ����
***************************************/
int RMDirAndFiles(const char* pDirName);

int PointToLineEquation(double x1,double y1,double x2,double y2,
						double &para_a,double &para_b,double &para_c);

int PointToLineEquation(int x1,int y1,int x2,int y2,
						double &para_a,double &para_b,double &para_c);


//�����������꣬�����߶�ƫ����������
int PointToLine(int x1,int y1,int x2,int y2,int &NumPoi,
				int * (&pX),int *(&pY));

int PointToLine( double x1,double y1,double x2,double y2,
				double xo,double yo,double cellsize,
				int &NumPoi,int * (&pX),int * (&pY) );

int TwoPoint2Angle(double x1,double y1,double x2,double y2,int bHudu,double &AngleRes);

/**** �����ߺ��߶εĽ�����Ŀ
�������߶˵�����ͷ����(����Ϊ�㣬��ʱ��Ϊ��)�������߶ζ˵�����
�����ϢRes˵����
	1��������ֱ�߲��غϣ������߶β��ཻ
	2��������ֱ�߲��غϣ������߶��ཻ����ͨ��
	3��������ֱ�߲��غϣ������߶��ཻ���߶ζ˵�
	4��������ֱ���غϣ����߶β���������
	5��������ֱ���غϣ����߶���������
***************/
int RadialMeetLine2Poi(
					   double r_x0,		//���ߺ�����
					   double r_y0,		//����������
					   double r_dir,	//���߷�������Ϊ�㣬��ʱ��Ϊ��������
					   double l_x1,		//�߶ζ˵�����
					   double l_y1,
					   double l_x2,
					   double l_y2,
					   int &Res			//�������ϸ��Ϣ��˵��
					   );


//�����ַ����Ϸ��Լ�麯��
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

//���ݽṹ

//�����
typedef struct tagQBFPolygon{
	struct tagQBFPolygon *pNext;
	int Num;
	double *pX;
	double *pY;		
}QBFPolygon;

//���������
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

//������غ���
int ts_ConcaveToBulge_iteration();

#endif
