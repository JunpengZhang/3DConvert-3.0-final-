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

/// ���Բ�ֵ����
///
/// ��֪(a0,b0)��(a1,b1)��a������a��Ӧ��ֵ
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

//����3ά�ڴ�ռ�
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

	/// ��ÿһ�з���ռ�
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

	/// �������ʧ��
	if(pData[0]==NULL)
	{
		delete [] pData;
		return NULL;
	}

	return (void **)pData;
}


//���㻭�ߣ��������˵㸡�����꣬���½ǵ�ĸ�������ͷֱ��ʣ�����������߶ε�����
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

	//��������˵����������
	X1 = int((x1-xo)/cellsize+0.5);
	Y1 = int((y1-yo)/cellsize+0.5);
	X2 = int((x2-xo)/cellsize+0.5);
	Y2 = int((y2-yo)/cellsize+0.5);

	//�����һ����
	if(X1==X2&&Y1==Y2)
	{
		pX = new int[1];
		pY = new int[1];
		pX[0] = X1;
		pY[0] = Y1;
		NumPoi = 1;
		return 1;
	}
		
	//���ֱ�߷���
	PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c);

	//��ɨ��
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

//���㻭�ߣ��������˵㸡�����꣬���½ǵ�ĸ����������������ķֱ��ʣ�����������߶ε�����
// ע�⣬���㻭���޶���x��y����ķ�Χ������������ְ�������Сֵȡ
// �ú������иĽ���أ������޷���Ӱ�죬���ܳ�����������ͬ���꣬Ŀǰû���ų���ͬ����
int MpGlobalFunc::PointToLineWithRgn( double x1,double y1,double x2,double y2,
				double xo,double yo,double csX,double csY, int minX,int maxX,int minY,int maxY,
				int &NumPoi,int * (&pX),int * (&pY) )
{
	int i,RV;

	/// ���ȵ��ò����Ʒ�Χ�����㻭���㷨
	RV = PointToLine( x1, y1, x2, y2, xo, yo, csX, csY, NumPoi, pX, pY );
	if(RV!=1)
		return RV;

	/// �������x y��������޷�
	for( i=0; i<NumPoi; i++ )
	{
		pX[i] = (pX[i]>maxX)? maxX : pX[i];
		pY[i] = (pY[i]>maxY)? maxY : pY[i];

		pX[i] = (pX[i]<minX)? minX : pX[i];
		pY[i] = (pY[i]<minY)? minY : pY[i];
	}

	return 1;
}

//���㻭�ߣ��������˵㸡�����꣬���½ǵ�ĸ����������������ķֱ��ʣ�����������߶ε�����
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

	//��������˵����������
	X1 = (int)floor((x1-xo)/csX+0.5);
	Y1 = (int)floor((y1-yo)/csY+0.5);
	X2 = (int)floor((x2-xo)/csX+0.5);
	Y2 = (int)floor((y2-yo)/csY+0.5);

	//�����һ����
	if(X1==X2&&Y1==Y2)
	{
		pX = new int[1];
		pY = new int[1];
		pX[0] = X1;
		pY[0] = Y1;
		NumPoi = 1;
		return 1;
	}
		
	//���ֱ�߷���
	PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c);

	//��ɨ��
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



//���������Ƕȣ���ǶȲ��λ��
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
��������InterLinear1D
������
    һά���Բ�ֵ�����������������꣨x1,v1������x2,v2����
���Բ�ֵ��X=xʱY��ֵv��ע��x��Ҫ��x1��x2֮�䡣

���ñ������ĺ����嵥��
	
���������
	double x1,  ��1��X����
	double v1,	��1��Y����
	double x2,  ��2��X����
	double v2,  ��2��Y����
	double x, 	����ֵ��X    
���������
	double & v, ��ֵ���

����ֵ��
	1���ɹ���
	-1��ʧ�ܡ�
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

	//���x1��x2̫�ӽ�����Ϊ��һ��ֵ����ֱ�ӽ�v1��Ϊ�������
	if( fabs(x1-x2) < qxwLimMinValue )
	{
		v = v1;
		return 1;
	}

	//���x����x1��x2֮�䣬���󷵻�
	//if( x<minX || x>maxX )
	//	return -1;
		
	//���Բ�ֵ
	v = v1 + (v2-v1) / (x2-x1) * (x-x1) ;
	return 1;
}

int MpGlobalFunc::InterLinear2D( double x1, double x2, double y1, double y2, double v11, double v12, 
				  double v21, double v22, double x, double y, double &v )
{
	double v1,v2;
	int RV;
	
	//���ȶ�x���в�ֵ
	RV = InterLinear1D(x1,v11,x2,v21, x, v1);
	if(RV==-1)	return -1;

	RV = InterLinear1D(x1,v12,x2,v22, x, v2);
	if(RV==-1)	return -1;

	//�ٶ�y���в�ֵ
	RV = InterLinear1D(y1,v1,y2,v2, y, v);
	if(RV==-1)	return -1;

	return 1;
}


/************************************************************************
	��������:	PointToLineEquation   
	��������:	�����������꣬��ֱ�߷��̡�
	�㷨˼��:	
	�����б�:   
				x1���߶ζ˵�1��X����
				y1���߶ζ˵�1��Y����
				x2���߶ζ˵�2��X����
				y2���߶ζ˵�2��Y����
				para_a��ֱ�߷���ϵ��
				para_b��ֱ�߷���ϵ��
				para_c��ֱ�߷���ϵ��
				minV������
				
	���ؽ����  
				-1������̫���������м���
				1���ɹ�
************************************************************************/
int MpGlobalFunc::PointToLineEquation(double x1,double y1,double x2,double y2,double &para_a,
						double &para_b,double &para_c,double minV)
{
	//�������̫���������м���ֱ�ӷ���
	if( sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) ) < minV )
		return ERR_TOONEAR;
		
	//ֱ�߽ӽ�ˮƽ
	if( fabs(x1-x2)>fabs(y1-y2) )
	{
		para_a = -(y1-y2)/(x1-x2);
		para_c = -y1-para_a*x1;
		para_b = 1.0;
	}

	//ֱ�߽ӽ���ֱ
	else
	{
		para_b = -(x1-x2)/(y1-y2);
		para_c = -x1-para_b*y1;
		para_a = 1.0;
	}
	
	return 1;
}

/************************************************************************
	��������:	Cal_TwoLine_Cross   
	��������:	������ֱ�߽������ꡣ
	�㷨˼��:	
	�����б�:   
				a1��ֱ��1��׼���̲���
				b1��ֱ��1��׼���̲���
				c1��ֱ��1��׼���̲���
				a2��ֱ��2��׼���̲���
				b2��ֱ��2��׼���̲���
				c2��ֱ��2��׼���̲���
				x������X����
				y������Y����
				
	����ֵ����  
				-1���޷�����
				1���ɹ�
************************************************************************/
int MpGlobalFunc::Cal_TwoLine_Cross( double a1, double b1, double c1, double a2, 
					  double b2, double c2, double &x, double &y )
{
	double Up,Down;
	
	Up = (b1*c2-b2*c1);
	Down = (a1*b2-a2*b1);
	if( fabs(Up)>fabs(Down) && fabs(Down)/fabs(Up)<1e-300 )	//�ж��Ƿ�������
		return ERR_CANNOTPROC;
	x = Up/Down;	

	Up = (c2*a1-c1*a2);
	Down = (b1*a2-b2*a1);
	if( fabs(Up)>fabs(Down) && fabs(Down)/fabs(Up)<1e-300 )	//�ж��Ƿ�������
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
	��������:	Cal_PoiLine_Relation   
	��������:	�жϵ����߶ε�λ�ù�ϵ��
	�㷨˼��:	
	�����б�:   
				poiX�����X����
				poiY�����Y����
				x1���߶ζ˵�1��X����
				y1���߶ζ˵�1��Y����
				x2���߶ζ˵�2��X����
				y2���߶ζ˵�2��Y����
				minV�����㾫�ȣ��ɲ��������룬Ĭ��Ϊ0.00001
				
	���ؽ����  
				-1���߶����˵�̫���������м���
				1�������߶ζ˵��غ�
				2����Ϊ�߶��ڵĵ�
				3�������߶����ڵ�ֱ���ϣ������߶��⣬
				4���㲻���߶����ڵ�ֱ����
************************************************************************/
int MpGlobalFunc::Cal_PoiLine_Relation( double poiX, double poiY, double x1, double y1, 
					   double x2, double y2, double minV )
{
	double pa_a,pa_b,pa_c;
	double dis;
	int bOnLine,bPort,bInLine;

	//�������̫���������м���ֱ�ӷ���
	if( sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) ) < minV )
		return ERR_TOONEAR;

	bOnLine = bPort = bInLine = 0;		//��ʼ����־��
	PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c,minV);	//��ֱ�߷���
	
	//���õ���ֱ�ߵľ����жϵ��Ƿ���ֱ����
	dis = fabs(pa_a*poiX+pa_b*poiY+pa_c) / sqrt(pa_a*pa_a+pa_b*pa_b) ;
	if( dis<minV )
		bOnLine = 1 ;

	//�жϵ��Ƿ����߶ζ˵�
	if( 
		( sqrt( (poiX-x1)*(poiX-x1)+(poiY-y1)*(poiY-y1) ) < minV ) ||
		( sqrt( (poiX-x2)*(poiX-x2)+(poiY-y2)*(poiY-y2) ) < minV )
		)
		bPort = 1;
	
	//�жϵ��Ƿ����߶����˵����귶Χ��
	if( fabs(x1-x2) >= fabs(y1-y2) )	//б��С��1�����ŵ��߶�
	{
		if( poiX>min(x1,x2) && poiX<max(x1,x2) )
			bInLine = 1;
	}
	else							//б��С��1�����ŵ��߶�
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
	��������:	Cal_TwoLine_Relation   
	��������:	�ж��߶����߶ε�λ�ù�ϵ��
	�㷨˼��:	
	�����б�:   
				x1���߶�1�˵�1��X����
				y1���߶�1�˵�1��Y����
				x2���߶�1�˵�2��X����
				y2���߶�1�˵�2��Y����
				x3���߶�2�˵�1��X����
				y3���߶�2�˵�1��Y����
				x4���߶�2�˵�2��X����
				y4���߶�2�˵�2��Y����
				Type���߶ι�ϵ����
				Lin1Type���������߶�1�Ĺ�ϵ
				Lin2Type���������߶�2�Ĺ�ϵ
				minV�����㾫�ȣ��ɲ��������룬Ĭ��Ϊ0.00001
				
	�����
			Type =	0���߶����߶�����ֱ���غϣ�Line1Type��Line2Type��Ч
					1���߶����߶�����ֱ��ƽ�У�Line1Type��Line2Type��Ч
					2���߶����߶�����ֱ���ཻ
			Line1Type =	0���������߶�1��
						1������Ϊ�߶�1�˵�
						2���������߶�1��
			Line2Type =	0���������߶�2��
						1������Ϊ�߶�2�˵�
						2���������߶�2��
	  
	����ֵ��  
				-1�������̫���������м���
				1���ɹ�
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

	//�ֱ��ж������߶ε������˵��Ƿ�̫����̫������ֱ�ӷ���
	dis1 = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
	dis2 = sqrt( (x3-x4)*(x3-x4) + (y3-y4)*(y3-y4) );
	if( dis1<=minV || dis2<=minV )
		return ERR_TOONEAR;

	//�������߶�ֱ�߷���
	PointToLineEquation(x1,y1,x2,y2,a1,b1,c1,minV);
	PointToLineEquation(x3,y3,x4,y4,a2,b2,c2,minV);

	//ȡ�߶�1���е�
	x0 = (x1+x2)/2;
	y0 = (y1+y2)/2;

	//�����߶�1�����㵽�߶�2����ֱ�ߵľ���
	dis1 = fabs(a2*x1+b2*y1+c2) / sqrt(a2*a2+b2*b2) ;
	dis2 = fabs(a2*x2+b2*y2+c2) / sqrt(a2*a2+b2*b2) ;
	dis3 = fabs(a2*x0+b2*y0+c2) / sqrt(a2*a2+b2*b2) ;

	//�ж�������ֱ�ߵĹ�ϵ
	if( fabs(dis1)<=minV && fabs(dis2)<=minV && fabs(dis3)<=minV )
		Type = 0;
	else if( fabs(dis1-dis2)<=minV && fabs(dis2-dis3)<=minV )
		Type = 1;
	else
		Type = 2;

	//����غϻ���ƽ�У�ֱ�ӷ��أ����ٽ����ж�
	if(Type!=2)
		return 1;

	RV = Cal_TwoLine_Cross(a1,b1,c1,a2,b2,c2,croX,croY);	//��������ֱ�߽�������

	//���߶�1�뽻������
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

	//���߶�2�뽻������
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
	��������:	Cal_RadialLine_Relation   
	��������:	�ж��������߶εĹ�ϵ��
	�㷨˼��:	
	�����б�:   
				r_x�����߶˵�X����
				r_y�����߶˵�Y����
				r_dir�����߷���x����Ϊ�㣬��ʱ��������λ��
				x1���߶ζ˵�1��X����
				y1���߶ζ˵�1��Y����
				x2���߶ζ˵�2��X����
				y2���߶ζ˵�2��Y����
				Type���������߶εĹ�ϵ����
				RadType�����������ߵĹ�ϵ
				LineType���������߶εĹ�ϵ
				minV�����㾫�ȣ��ɲ��������룬Ĭ��Ϊ0.00001
				
	�����
			type =	0���������߶�����ֱ���غϣ�RadType��LineeType��Ч
					1���������߶�����ֱ��ƽ�У�RadType��LineeType��Ч
					2���������߶�����ֱ������ֱ���ཻ
			RadType =	0�����㲻��������
						1������Ϊ���߶˵�
						2��������������
			LineType =	0���������߶���
						1������Ϊ�߶ζ˵�
						2���������߶���		
	  
	����ֵ��  
		-1�������̫���������м���
		1���ɹ�
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

	ang = r_dir/180*qxwPI;	//ת��Ϊ����

	//����������һ��������루len��ȡ����
	rx1 = r_x + len*cos( ang );
	ry1 = r_y + len*sin( ang );
	rx2 = rx1 + len*cos( ang );
	ry2 = ry1 + len*sin( ang );

	//�����߶�ֱ�߷���
	RV = PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c);
	if(RV!=1)
		return RV;

	//��������������ֱ�ߵľ���
	dis1 = fabs( r_x*pa_a + r_y*pa_b + pa_c ) / sqrt( pa_a*pa_a + pa_b*pa_b );
	dis2 = fabs( rx1*pa_a + ry1*pa_b + pa_c ) / sqrt( pa_a*pa_a + pa_b*pa_b );
	dis3 = fabs( rx2*pa_a + ry2*pa_b + pa_c ) / sqrt( pa_a*pa_a + pa_b*pa_b );

	//�ж�������ֱ�ߵĹ�ϵ
	if( fabs(dis1)<=minV && fabs(dis2)<=minV && fabs(dis3)<=minV )
		Type = 0;
	else if( fabs(dis1-dis2)<=minV && fabs(dis2-dis3)<=minV )
		Type = 1;
	else
		Type = 2;

	//����������߶�����ֱ�߲��ཻ��ֱ�ӷ���
	if(Type!=2)
		return 1;

	//������ֱ�߷���
	PointToLineEquation(r_x,r_y,rx1,ry1,r_a,r_b,r_c);

	//��������ֱ�߽�������
	//crossX = ( pa_b*r_c - r_b*pa_c ) / ( pa_a*r_b - r_a*pa_b );
	//crossY = ( r_c*pa_a - pa_c*r_a ) / ( pa_b*r_a - r_b*pa_a );
	RV = Cal_TwoLine_Cross(pa_a,pa_b,pa_c,r_a,r_b,r_c,crossX,crossY);
	if(RV!=1)
		return RV;
	
	//�жϽ��������ߵĹ�ϵ
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
	
	//�жϽ������߶εĹ�ϵ
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
	��������:	Cal_PoiPolygon_Relation   
	��������:	�жϵ������ε�λ�ù�ϵ��
	�㷨˼��:	
	�����б�:   
				NumPoi������νǵ����
				pX������νǵ�X��������ָ�룬˳ʱ������
				pY������νǵ�Y��������ָ�룬˳ʱ������
				PoiX����X����
				PoiY����Y����				
				minV�����㾫�ȣ��ɲ��������룬Ĭ��Ϊ0.00001
		  
	����ֵ��  
				0���ڲ�
				1���ⲿ
				2�����ڶ���α���
				-1�����������߱�����δ����Ч�жϳ����߹�ϵ
************************************************************************/
int MpGlobalFunc::Cal_PoiPolygon_Relation( int NumPoi, double *pX, double *pY, 
						  double PoiX, double PoiY, double minV )
{
	int CroDotNumL,CroDotNumR;	//���ҽ������
	double Dir0=0;				//��ʼ����0��
	double DirDlt=20;			//���򲽳���20��
	double curDir;
	int i;
	int EviNumLim=3;			//֤����Ŀ����
	double x1,y1,x2,y2;			//�߶����˵�
	double x0,y0,dir;			//�������ͷ���
	int RV;
	int Type,RadType,LineType;	//�������߶ι�ϵ��������ֵ
	int evidenceOut,evidenceIn;	//֤�ݱ���
	double maxX,minX,maxY,minY,maxDis,minDis,dis,LenAxis;
	double minLen;
	double angView;		//��Զ���ε��ӳ���
	int NumLoop;		//��ѭ������

	//�����жϵ��Ƿ��ڶ���εı���
	//����������̱߳�
	minLen = sqrt( (pX[0]-pX[1])*(pX[0]-pX[1]) + (pY[0]-pY[1])*(pY[0]-pY[1]) );
	for( i=0; i<NumPoi; i++ )
	{
		//ȷ���߶����˵�
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

		//�жϵ��Ƿ��ڵ�ǰ�ı��ϣ�����ڣ�ֱ�ӷ���
		RV = Cal_PoiLine_Relation( PoiX, PoiY, x1, y1, x2, y2, minV );
		if( RV==1 || RV==2 )
			return 2;
	}

	//��������Ӿ��η�Χ
	//������㵽����ζ������Զ���������
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

	//����Ӿ��ζԽ��߳������˶���θ���������᳤
	maxX = maxX-minX;
	maxY = maxY-minY;
	LenAxis = sqrt(maxX*maxX+maxY*maxY);

	//����㵽����ζ�������������᳤������˵�����ڶ������
	if( minDis>LenAxis+minV )
		return 1;	

	//���Զ���ε���С�ӳ���
	angView = atan2( minLen/2, maxDis ) *2 ;
	angView = angView/qxwPI*180;

	//�㲻�ڶ���εı��ϣ�ͨ����������ཻ���ж϶��߹�ϵ
	//�ԽǶ���ѭ������0�ȿ�ʼ������Ϊ30�Ƚ�����ÿ���Ƕ��ж�������ߵ��ཻ���
	//������������Ƕ�֤ʵ���ڶ�����ڲ�������Ϊ���ڲ��㣬��������
	//������������Ƕ�֤ʵ���ڶ�����ⲿ������Ϊ���ⲿ�㣬��������
	evidenceOut = evidenceIn = 0;
	NumLoop = 1;
	for( curDir=Dir0; ; curDir+=DirDlt )
	{
		//�����ѭ����������5���˳�
		if(NumLoop>5)
			break;
		
		//������ڵ���Dir0+180
		if( curDir>=Dir0+180 )
		{
			Dir0 += angView/5;
			curDir = Dir0;
			NumLoop ++;
		}
		
		//��ÿ���߽���ѭ�����㣬��������������ߵĽ���
		CroDotNumR = CroDotNumL = 0;
		for( i=0; i<NumPoi; i++ )
		{
			//ȷ���߶����˵�
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

			//ȷ�����������ͷ���
			x0 = PoiX;
			y0 = PoiY;
			dir = curDir/180*qxwPI;

			//�ж����ߺ��߶ε��ཻ���
			Cal_RadialLine_Relation( x0, y0, curDir, x1, y1, x2, y2, Type, RadType, LineType );
			if( Type==2 && LineType==2 && RadType==2 )
				CroDotNumR++;

			//ȷ�����������ͷ���
			x0 = PoiX;
			y0 = PoiY;
			dir = (curDir+180)/180*qxwPI;

			//�ж����ߺ��߶ε��ཻ���
			Cal_RadialLine_Relation( x0, y0, curDir+180, x1, y1, x2, y2, Type, RadType, LineType );
			if( Type==2 && LineType==2 && RadType==2 )
				CroDotNumL++;				
		}

		//������ҽ�����������Ϊż�����˷�����Ч
		if( (CroDotNumL+CroDotNumR)%2==0 )	//&& (CroDotNumL+CroDotNumR)>0
		{
			if( CroDotNumL%2==0 && CroDotNumR%2==0 )
				evidenceOut ++;
			else
				evidenceIn ++;
		}

		//����ڲ����ⲿ��֤�ݴ������ڵ���3���˳�ѭ��
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
	��������:	Cal_LinePolygon_Relation   
	��������:	�ж��߶������ε�λ�ù�ϵ��
	�㷨˼��:	
	�����б�:   
				NP������νǵ����
				pX������νǵ�X��������ָ�룬˳ʱ������
				pY������νǵ�Y��������ָ�룬˳ʱ������
				X1���߶ζ˵�1X����
				Y1���߶ζ˵�1Y����
				X2���߶ζ˵�2X����
				Y2���߶ζ˵�2Y����
				minV�����㾫�ȣ��ɲ��������룬Ĭ��Ϊ0.00001
		  
	����ֵ��  
			-1���������
			0���߶��ڶ������
			1���߶��ڶ������
			2���߶��������ཻ
			3���߶��ڶ���εı���
			4�����ཻ
************************************************************************/
int MpGlobalFunc::Cal_LinePolygon_Relation( int NP, double *pX, double *pY, double X1, double Y1, 
							 double X2, double Y2, double minV )
{
	double crox,croy;
	int i,j,num;
	int evIn;	//���ڶ�����ڲ��ĸ���
	int evOut;	//���ڶ�����ⲿ�ĸ���
	int evOn;	//���ڶ���α��ϵĸ���
	int RV,RV1,RV2;
	double x1,x2,y1,y2;
	int Type,Line1Type,Line2Type;
	double *px,*py;
	
	px = py = NULL;

	//�ж��߶������˵��Ƿ�λ�ڶ���ε�һ������
	//����ڣ�˵���߶��ڶ���α���
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

		//���˵�1��ߵĹ�ϵ
		RV1 = Cal_PoiLine_Relation( X1, Y1, x1, y1, x2, y2, minV );
		RV2 = Cal_PoiLine_Relation( X2, Y2, x1, y1, x2, y2, minV );
		if( (RV1==1||RV1==2) && (RV2==1||RV2==2) )
			return 3;
	}

	//�ж��߶������˵������ι�ϵ
	RV1 = Cal_PoiPolygon_Relation( NP, pX, pY, X1, Y1, minV );
	RV2 = Cal_PoiPolygon_Relation( NP, pX, pY, X2, Y2, minV );
	if( RV1 ==-1 ) 
		return RV1;
	if( RV2 ==-1 )
		return RV2;

	//���һ�����ڲ�һ�����ⲿ˵���ཻ
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

		//��ǰ�����߶εĹ�ϵ
		RV = Cal_TwoLine_Relation( x1, y1, x2, y2, X1, Y1, X2, Y2, 
			Type, Line1Type, Line2Type, minV );

		//����߶κͱ��ཻ��������뽻��������
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

	//���߶ζ˵�Ҳ��������
	px[num] = X1;
	py[num] = Y1;
	num++;
	px[num] = X2;
	py[num] = Y2;
	num++;
	
	//�Խ��������������
	if( fabs(X1-X2)>fabs(Y1-Y2) ) //���򣬰�X����
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
	else						//���򣬰�Y����
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

	//���������ڽ�����е�����ж�
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

	//ֻ���ⲿ��û���ڲ�Ҳû�б��ϵģ��ⲿ
	if( evIn==0 && evOn==0 && evOut>0 )	
		return 1;

	//���ⲿ��û���ڲ����б��ϵģ����ཻ
	if( evIn==0 && evOn>0 && evOut>0 )	
		return 4;

	//ֻ���ڲ����ڲ�
	if( evIn>=0 && evOn==0 && evOut==0 )
		return 0;

	//���ڲ����б��ϣ�û���ⲿ���ڲ�
	if( evIn>=0 && evOn>0 && evOut==0 )
		return 0;

	//���ڲ����ޱ��ϣ����ⲿ���ཻ
	if( evIn>0 && evOn==0 && evOut>0 )
		return 2;

	//���ڲ����б��ϣ����ⲿ���ཻ
	if( evIn>0 && evOn>0 && evOut>0 )
		return 2;

	return ERR_CANNOTPROC;
}



/************************************************************************
	��������:	Cal_TwoPolygon_Relation   
	��������:	�ж϶����1������ڶ����2��λ�ù�ϵ��
	�㷨˼��:	
	�����б�:   
				NP1�������1�ǵ����
				pX1�������1�ǵ�X��������ָ�룬˳ʱ������
				pY1�������1�ǵ�Y��������ָ�룬˳ʱ������
				NP2�������2�ǵ����
				pX2�������2�ǵ�X��������ָ�룬˳ʱ������
				pY2�������2�ǵ�Y��������ָ�룬˳ʱ������
				minV�����㾫�ȣ��ɲ��������룬Ĭ��Ϊ0.00001
		  
	����ֵ��  
			-1���������
			0��1��2�ڲ�
			1��1��2�ཻ
			2��1��2���ڣ��б��غϻ��ж������������εı���
			3��1��2���룬��Ϊ�ⲿ
			4��2��1�ڲ�
			5���غ�
************************************************************************/
int MpGlobalFunc::Cal_TwoPolygon_Relation( int NP1, double *pX1, double *pY1,
							int NP2, double *pX2, double *pY2, double minV )
{
	int i;
	int RV;
	double x1,x2,y1,y2;
	int Num10;	//�����1�ı��ڶ����2�ڲ��ĸ���
	int Num11;	//�����1�ı��ڶ����2�ⲿ�ĸ���
	int Num12;	//�����1�ı�������2�ཻ�ĸ���
	int Num13;	//�����1�ı��ڶ����2���ϵĸ���
	int Num14;	//�����1��������2���ཻ�ĸ���
	int Num20;	//�����2�ı��ڶ����1�ڲ��ĸ���
	int Num21;	//�����2�ı��ڶ����1�ⲿ�ĸ���
	int Num22;	//�����2�ı�������1�ཻ�ĸ���
	int Num23;	//�����2�ı��ڶ����1���ϵĸ���
	int Num24;	//�����2��������1���ཻ�ĸ���
	
	//�ȶԶ����1��ÿ����������2��ϵ�����ж�
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

	//�ٶԶ����2��ÿ����������1��ϵ�����ж�
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

	//�غ�
	if( Num13==NP1 && Num23==NP2 )
		return 5;

	//�����1��ÿ���߾��ڶ����2�ڲ������ڶ����2�ı��ϣ��ұ����ڲ�
	//������1�ڶ����2�ڲ�
	if( Num10+Num13==NP1 && Num10>0 )	
		return 0;

	//�����2��ÿ���߾��ڶ����1�ڲ������ڶ����1�ı��ϣ��ұ����ڲ�
	//������2�ڶ����1�ڲ�
	if( Num20+Num23==NP2 && Num20>0 )	
		return 4;

	//1��2�ֻཻ�������������1��1�ı���2�ཻ ��2��1�ı���2�㽻�����ڲ�
	if( 
		Num12>0 || 
		( Num10>0 && Num14>0 )
		)
		return 1;

	//�����1�ı߾��ڶ����2�ⲿ�����Ҷ����2���ڶ����1�ڲ�
	//��������
	if(Num11==NP1)
		return 3;

	//�����1�ı�������2�Ĺ�ϵ�㽻���ⲿ��λ�ڱ����غϣ��Ҳ���ֻ���ⲿ
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

	//����ļ��Ƿ����
	pF = fopen(FNNew,"rb");
	if( pF )
	{
		bEx = 1;
		fclose(pF);
	}
	else
		bEx = 0;
	

	//����ļ��Ѿ������Ҹ���ģʽΪ��ֱ�ӷ���0
	if( bEx==1 && bOverWrite==0 )
		return 0;

	//�����ļ����������0���ֽڣ�ֱ�ӷ���0
	RV = ReadFileToBuf( FNOld, pDat, Len );
	//20161223lls add start 304��̬���������޸�
	if(pDat == NULL)
	{
		return 0;
	}
	//20161223lls add end 304��̬���������޸�
	if( Len==0 )
	{
		return 0;
	}

	//д�����ļ�
	pF = fopen(FNNew,"wb+");
	//20161223lls add start 304��̬���������޸�
	if(pF == NULL)
	{
		return 0;
	}
	//20161223lls add end 304��̬���������޸�
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
/// �ַ�����������							///
/// ��д�����ǿ							///
/// ���������								///
/// Src			Դ�ַ���					///
/// Splitter	�ָ��ַ����б�				///
/// ���������								///
/// vecStrings	�����õ����ַ����б�		///
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

	/// �����һ��Ԫ�ص���ϣ�ֱ�����ÿ��Ԫ��
	if(NumSel==1)
	{
		for(i=0;i<NumAll;i++)
		{
			ppRes[NumRes][0] = i;
			NumRes++;
		}
		return 1;
	}

	/// �������һ������ϣ��ݹ������Ŀ��1�����
	RV = zuhe(NumAll,NumSel-1,NumR, ppR );
	if(RV!=1)
	{
		if (ppR) //20161225 304pingce
		{
			delete []ppR;
		}
		return -1;
	}
		

	/// �Լ�1��ϵ������������1��Ԫ��
	for( i=0; i<NumR; i++ )
	{
		/// ��ǰ��Ͻ�������ֵ
		maxV = -1;
		for(j=0;j<NumAll;j++)
		{
			if(ppR[i][j]==-1)
				break;
			maxV = max(maxV,ppR[i][j]);
		}
		SNBeg = j;

		/// �Ե�ǰNumSel-1Ԫ����Ͻ��׷��һ��Ԫ��
		/// �γɶ��NumSelԪ�����
		for( j=maxV+1; j<NumAll; j++ )
		{
			/// ��ǰ��NumSel-1��Ԫ�ظ���
			for( k=0; k<SNBeg; k++ )
				ppRes[NumRes][k] = ppR[i][k];

			/// ׷��һ��Ԫ��
			ppRes[NumRes][SNBeg] = j;
			NumRes++;
		}
	}

	MpGlobalFunc::DeleteSpace2D(ppR);
	return 1;		
}

/// ������������Ƿ��ཻ
/// 1���ཻ��0�����ཻ
int MpGlobalFunc::Chk2RegularRect(double xl1,double yb1,double xr1,double yt1,
						   double xl2,double yb2,double xr2,double yt2)
{
	/// �ж�������������Ƿ��ཻ�������ж�������������Ƿ��ཻ
	/// ������ཻ����������Ҫ�����ĸ�����֮һ
	/// ����1�ھ���2��࣬1��2�Ҳ࣬1��2���ϣ�1��2�·�
	/// ����ĸ���������������˵�����������ཻ

	/// ����1�ھ���2���
	if( xr1<xl2 )
		return 0;

	/// ����1�ھ���2�Ҳ�
	if( xl1>xr2 )
		return 0;

	/// ����1�ھ���2�·�
	if( yt1<yb2 )
		return 0;

	/// ����1�ھ���2�Ϸ�
	if( yb1>yt2 )
		return 0;

	/// ������������������������������ཻ
	return 1;

	///// ������ԭ������һ��˼·�����ཻ���ཻ������ܽ᲻ȫ���㷨����
	//
	//// �жϵ�һ����Ӿ����ĸ��ǵ��Ƿ�λ�ڵ�2��������
	//// ���жϵڶ�����Ӿ����ĸ��ǵ��Ƿ�λ�ڵ�1��������
	//// �����������˵������Ӿ��β����

	//// ��һ���������½ǵ��Ƿ��ڵڶ���������
	//if( xl1>=xl2 && xl1<=xr2 && yb1>=yb2 && yb1<=yt2 )
	//	return 1;

	//// ��һ���������Ͻǵ��Ƿ��ڵڶ���������
	//if( xl1>=xl2 && xl1<=xr2 && yt1>=yb2 && yt1<=yt2 )
	//	return 1;

	//// ��һ���������Ͻǵ��Ƿ��ڵڶ���������
	//if( xr1>=xl2 && xr1<=xr2 && yt1>=yb2 && yt1<=yt2 )
	//	return 1;

	//// ��һ���������½ǵ��Ƿ��ڵڶ���������
	//if( xr1>=xl2 && xr1<=xr2 && yb1>=yb2 && yb1<=yt2 )
	//	return 1;

	//// �ڶ����������½ǵ��Ƿ��ڵ�һ��������
	//if( xl2>=xl1 && xl2<=xr1 && yb2>=yb1 && yb2<=yt1 )
	//	return 1;

	//// �ڶ����������Ͻǵ��Ƿ��ڵ�һ��������
	//if( xl2>=xl1 && xl2<=xr1 && yt2>=yb1 && yt2<=yt1 )
	//	return 1;

	//// �ڶ����������Ͻǵ��Ƿ��ڵ�һ��������
	//if( xr2>=xl1 && xr2<=xr1 && yt2>=yb1 && yt2<=yt1 )
	//	return 1;

	//// �ڶ����������½ǵ��Ƿ��ڵ�һ��������
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
	
	/// ��������ԴĿ¼��Ŀ��Ŀ¼�ַ����Ƿ�����
	if( DirS==NULL || strlen(DirS)==0 )
		return -1;

	if( DirO==NULL || strlen(DirO)==0 )
		return -1;

	if( NumTry<1 || InterTime<1 )
		return -1;
	
	/// ���Ŀ���ļ����Ƿ���ڣ�������ڣ���ɾ��
	Dir2 = DirO;
	Dir2.TrimRight("\\");
	if( myFind.FindFile(Dir2) )
	{
		/// ����ɾ��NumTry�Σ��ɹ�ֱ������
		for( bOK=0,i=0; i<NumTry; i++ )
		{
			RV = RMDirAndFiles(Dir2);
			if(RV==1)
			{
				bOK=1;
				break;
			}

			/// �������Լ��InterTime����
			Sleep(InterTime);
		}

		/// �ж��Ƿ�ɾ���ɹ������ɹ�������
		if(bOK==0)
		{
			return -1;
		}		
	}

	/// ����Ŀ��Ŀ¼��������δ������ɹ�ֱ������
	for( bOK=0,i=0; i<NumTry; i++ )
	{
		RV = _mkdir(Dir2);
		if(RV==0)
		{
			bOK = 1;
			break;
		}

		/// �������Լ��InterTime����
		Sleep(InterTime);
	}

	/// ���δ�����ɹ���������
	if(bOK==0)
		return -1;

	//������̵�ǰĿ¼
	GetCurrentDirectory(499,BakCurPath);

	/// ��ԴĿ¼�������ļ���Ŀ¼���д���
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

		/// �����Ŀ¼���ݹ���ñ�����
		/// ע�⣬�ݹ���ñ��������Ѿ�������γ��ԣ�������ﲻ�ٶ�γ���
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

		/// �ߵ����˵�����ļ�����Ҫ���п���
		for( bOK=0,i=0; i<NumTry; i++ )
		{
			bSuc = CopyFile(D1,D2,FALSE);
			if(bSuc)
			{
				bOK=1;
				break;
			}

			/// �������Լ��InterTime����
			Sleep(InterTime);
		}
		if(bOK==0)
		{
			SetCurrentDirectory(BakCurPath);
			return -1;
		}
	}

	/// �����˵��Ŀ¼��������Ŀ¼���ļ��������ɹ���
	SetCurrentDirectory(BakCurPath);
	return 1;
}

/***************************************
���� ��ɾ��Ŀ¼�Լ�Ŀ¼�е��ļ�����Ŀ¼��
�㷨���ݹ����
***************************************/
int MpGlobalFunc::RMDirAndFiles(LPCSTR pDirName)
{
	CFileFind myFind;
	CString DirName,curName;
	char BakCurPath[500];
	int i,RV=1;
	
	//������̵�ǰĿ¼
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

		//�����Ŀ¼����ݹ����
		curName = myFind.GetFilePath();
		if(myFind.IsDirectory())
		{
			i = RMDirAndFiles(LPCSTR(curName));
			RV = (i==-1) ? -1 : RV ;
		}
		
		//������ļ�����ɾ���ļ�
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

	//�ָ����̵�ǰĿ¼
	bWorking = SetCurrentDirectory(BakCurPath);

	if(RV==-1)
		return RV;

	//ɾ��Ŀ��Ŀ¼
	i=0;	
	while( _rmdir(DirName)!=0 && ++i<=500 )
	{
		;
	}
	RV = (i>=500) ? -1 : RV ;
	return RV;
}

// �ƶ�һ�����ڵ���һ�����ڵ��ض�λ��
int MpGlobalFunc::PositionWndAroundAnother(CWnd & WndStill, CWnd & Wnd2Move, int iPosFlag)	//δ���ԣ�����������
{
	CRect rectWnd2Move;
	Wnd2Move.GetWindowRect(&rectWnd2Move);

	CRect rectWndStill;
	WndStill.GetWindowRect(&rectWndStill);

	switch(iPosFlag)
	{
	case 1: // �ױ������
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left,
				rectWndStill.bottom,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 2: // ���½� (��)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.bottom,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 3: // ����¶���
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 4: // ����϶���
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.top,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 5: // ���Ͻ� (��)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break; 

	case 6: // ���������
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left,
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 7: // �����Ҷ���
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;
	case 8: // ���Ͻ� (��)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right,
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 9: // �ұ��϶���
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right,
				rectWndStill.top,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 10: // �ұ��¶���
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right,
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 11: // ���½� (��)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right,
				rectWndStill.bottom,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 12: // �ױ��Ҷ���
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.bottom,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 13: // ���½� (��)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left,
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 14: // ���Ͻ� (��)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.left,
				rectWndStill.top,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break; 
	case 15: // ���Ͻ� (��)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.top,
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 16: // ���½� (��)
		{
			Wnd2Move.MoveWindow(
				rectWndStill.right - (rectWnd2Move.right - rectWnd2Move.left),
				rectWndStill.bottom - (rectWnd2Move.bottom - rectWnd2Move.top),
				rectWnd2Move.right - rectWnd2Move.left,
				rectWnd2Move.bottom - rectWnd2Move.top 
				);
		}
		break;

	case 17: // ����
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
	
	//������̵�ǰĿ¼
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

		//�����Ŀ¼������һ��
		curName = myFind.GetFilePath();
		if(myFind.IsDirectory())
		{
			continue;
		}
		
		//������ļ�����Ҫ�ж��ļ��Ƿ���Ҫɾ��
		else
		{
			/// ����ļ����в�����Ŀ���ַ���������һ��
			if( curName0.Find(ObjStr)==-1 )
				continue;

			/// ɾ����ǰ�ļ�
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

	//�ָ����̵�ǰĿ¼
	bWorking = SetCurrentDirectory(BakCurPath);

	if(RV==-1)
		return RV;
	else
		return RV;
}

//���������ַ����Ϸ��Լ�麯��
int MpGlobalFunc::StringIsInt(const char * pObj)
{
	int RV;
	const char *pChar;

	/// �ȼ���Ƿ�����
	RV = StringIsNumber(pObj);
	if(RV==0)
		return 0;

	/// �ټ������С����
	pChar = strstr(pObj,".");
	if(pChar)
		return 0;

	return 1;
}

//�����ַ����Ϸ��Լ�麯��
int MpGlobalFunc::StringIsNumber(const char * pObj)
{
	int i;
	int len;
	int posDot;		//С�����λ��
	int numDot;		//С�������Ŀ
	int posMinus;	//���ŵ�λ��
	int numMinus;	//���ŵ���Ŀ
	char maxChar='9';
	char minChar='0';

	len = strlen(pObj);
	if(len==0)
		return 0;
	
	//ͳ�Ƹ�����С�����λ�ú���Ŀ
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

	//����и��ţ�ֻ�ܳ���һ���ұ����ڵ�һ���ַ�λ�ã�����Ƿ�
	if(numMinus>1)
		return 0;
	else if( (numMinus==1) && (posMinus!=0) )
		return 0;

	//�����С���㣬ֻ�ܳ���һ������ǰ�󶼱���������
	if(numDot>1)
		return 0;
	else if( (numDot==1) && ( (posDot==0)||(posDot==len-1) ) )
		return 0;

	//ͳ�Ƴ��˸��ź�С����֮�⣬�Ƿ�ÿ���ַ�����0��9֮��
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
	
	/// �ж��������ݷ�Χ�Ƿ�����
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

	/// �����־�������سߴ�
	Hei = (int)ceil( (yt-yb)/cs );
	Wid = (int)ceil( (xr-xl)/cs );

	/// Ϊ��־��������ڴ�ռ�
	/// ע���־�����кŴ������ϵ���
	ppF = (char **)NewSpace2D(Hei,Wid,sizeof(char));
	if(ppF==NULL)
		return -1;

	/// ��ʼ����־����
	for( i=0; i<Hei; i++ )
		memset( ppF[i], 0, sizeof(char)*Wid );

	/// ����ʹ�����ַ������б�ǣ�

	/// ����һ
	/// һ���Ƕ�Ŀ������ÿ������ѭ����
	/// �ж�ÿ�������Ƿ�λ������һ�������ڣ�λ������1

	/// ������
	/// �Ծ����������ѭ�����ҵ�ÿ������������Ŀ��������غϲ��֣�
	/// ���غϲ��ֵ����ر��Ϊ1

	/// ����һ���� δ���
	/*if(0)
	{
		/// ��ÿ���������ѭ�����ж��Ƿ񱻾������鸲��
		for( i=0; i<Hei; i++ )
		{
			for( j=0; j<Wid; j++ )
			{
				x = xo + j*cs;
				y = yo + i*cs;


			}
		}
	}*/

	/// ����������
	if(1)
	{
		/// ��ÿ�����������������
		for( SN=0; SN<NumReg; SN++ )
		{
			/// �ҵ��غϲ���
			x1 = max( xl, pXl[SN] );
			x2 = min( xr, pXr[SN] );
			y1 = max( yb, pYb[SN] );
			y2 = min( yt, pYt[SN] );

			/// �ж��غϲ����Ƿ���Ч����Ч��ֱ������һѭ��
			if( x1>=x2 || y1>=y2 )
				continue;

			/// ��Ч��������Ч�����Ӧ���������ط�Χ
			j1 = int((x1-xo)/cs + 0.5);
			j2 = int((x2-xo)/cs + 0.5);
			i1 = int((y1-yo)/cs + 0.5);
			i2 = int((y2-yo)/cs + 0.5);

			/// ���غϷ�Χ���Ϊ1
			for( i=i1; i<=i2; i++ )
			{
				for( j=j1; j<=j2; j++ )
				{
					ppF[i][j] = 1;
				}
			}
		}
	}

	/// ȫ������꣬�ж��Ƿ񱻸���
    for( i=0; i<Hei; i++ )
	{
		for( j=0; j<Wid; j++ )
		{
			/// һ������������δ����ǣ�˵������δ���ǣ�����
			if(ppF[i][j]==0)
			{
				MpGlobalFunc::DeleteSpace2D(ppF);
				return 0;
			}
		}
	}

	/// �ߵ����˵��ȫ�����Ϊ1��Ŀ�����򱻸���
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

	/// ͳ���������귶Χ
	minX = maxX = pX[0];
	minY = maxY = pY[0];
	for( i=0; i<Num; i++ )
	{
		minX = min(minX,pX[i]);
		maxX = max(maxX,pX[i]);
		minY = min(minY,pY[i]);
		maxY = max(maxY,pY[i]);
	}

	/// �������ߴ粢����ռ�
	Hei = (int)ceil( (maxY-minY)/cs ) + 1;
	Wid = (int)ceil( (maxX-minX)/cs ) + 1;
	ppD = (unsigned char **)NewSpace2D(Hei,Wid,sizeof(char));
	if(ppD==NULL)
		return NULL;
	for( i=0; i<Hei; i++ )
		memset( ppD[i], 0, sizeof(char)*Wid );

	/// ��������������㻭��
	for( i=0; i<Num-1; i++ )
	{

		/// ���㻭��
		pXc = pYc = NULL;
		PointToLine( pX[i],pY[i],pX[i+1],pY[i+1], minX,minY,cs, NumPoi,pXc,pYc );

		/// ����ǰ�߶θ��ǵ�������Ϊ1
		for( j=0; j<NumPoi; j++ )
		{
			ppD[ Hei-1-pYc[j] ][ pXc[j] ] = 1;
		}

		delete [] pXc;
		delete [] pYc;
	}

	/// ������0��Ϊ255
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
		
	//������̵�ǰĿ¼
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

		//�����Ŀ¼��������һ��
		//curName = myFind.GetFilePath();
		if(myFind.IsDirectory())
		{
			continue;
		}
		
		//������ļ�����Ҫ�ж��ļ��Ƿ���Ҫɾ��
		else
		{
			/// ����ļ����в�����Ŀ���ַ���������һ��
			if( curName.Find(ObjStr)==-1 )
				continue;

			/// ��Ŀ��1
			Num++;

			//�����Ŀ�����һ�£����
			if(Num==SN)
			{
				FN = curName;
				return 1;
			}
		}
	}
	myFind.Close();

	//�ָ����̵�ǰĿ¼
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
	
	//������̵�ǰĿ¼
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

		//�����Ŀ¼������һ��
		//curName = myFind.GetFilePath();
		if(myFind.IsDirectory())
		{
			continue;
		}
		
		//������ļ�����Ҫ�ж��ļ��Ƿ���Ҫɾ��
		else
		{
			/// ����ļ����в�����Ŀ���ַ���������һ��
			if( curName.Find(ObjStr)==-1 )
				continue;

			/// ��Ŀ��1
			Num++;
		}
	}
	myFind.Close();

	//�ָ����̵�ǰĿ¼
	bWorking = SetCurrentDirectory(BakCurPath);

	return 1;
}

/// ����߶��Ƿ���Բ�����棬�Ƿ���1�����Ƿ���0�����󷵻�-1
int MpGlobalFunc::IsLineOutCircle(double x1,double y1,double x2,double y2,
								  double xo,double yo,double Rad,double minV)
{
	double pa,pb,pc;	//< ֱ�߷��̲���
	double dis;
	double xc1,yc1,xc2,yc2;
	int RV;
	double v1,v2;

	/// ��ֱ�߷���
	RV = PointToLineEquation(x1, y1,x2, y2,pa, pb,pc,minV);
	if(RV!=1)
		return -1;

	/// ��Բ�ĵ�ֱ�߾���
	dis = fabs( pa*xo + pb*yo + pc ) / sqrt( pa*pa + pb*pb );

	/// ���Բ�ĵ�ֱ�߾�����ڰ뾶���߶�һ����Բ����
	if(dis>Rad)
		return 1;

	/// �ж��߶�����ֱ����Բ�Ľ��㣬����һ��Ϊ������
	/// ���ֱ����Բ���У������ǵ��Ϊһ��
	RV = Cal_LineCircle_Cross( pa, pb, pc, xo, yo, Rad, xc1, yc1, xc2, yc2,minV);
	if(RV==-1)
		return -1;

	/// �ж����������Ƿ����߶�����
	/// ����߶�����x��ֵ�󣬶�x�����жϣ������y�����ж�
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

/// ���ַ�����β�Ŀո�ȥ����ע��ֱ�����ַ���ָ���и�д
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

	/// ȥ����߿ո�
	if(pObj[len-1]==' ')
	{
		/// �Ӻ���ǰ�ҵ���һ�����ǿո���ַ���ֱ����ʼ
		for(i=len-2;i>=0;i--)
		{
			if(pObj[i]!=' ')
			{			
				break;
			}
		}
		pObj[i+1] = '\0';
	}

	/// ȥ��ǰ�߿ո�
	if(pObj[0]==' ')
	{
		/// ��ǰ����ҵ���һ�����ǿո���ַ���ֱ��������
		for(i=1;i<len;i++)
		{
			if(pObj[i]!=' ')
			{			
				break;
			}
		}

		/// ��iָʾλ��֮���������ǰ��������һ���ڴ����
		strcpy(buf,pObj);
		strcpy(pObj,buf+i);
	}

	delete [] buf;
	return 1;
}

/// ����ֱ����Բ�Ľ��㣬����0���޽��㣬����1���������㣬����2��һ�����㣬���أ�1������
int MpGlobalFunc::Cal_LineCircle_Cross(double pa,double pb,double pc,double xo,double yo,
						double R,double &xc1,double &yc1,double &xc2,double &yc2,double minV)
{
	double m,n;
	double a,b,c;
	double delta,dis;

	/// ����ֱ�߷��̵Ĳ����ж϶�x���м��㻹�Ƕ�y���м���
	if( fabs(pb)>fabs(pa) )	//<  ��x���м���
	{
		/// ��ax+by+c=0ת��Ϊy=mx+n����ʽ
		m = -1*pa/pb;
		n = -1*pc/pb;

		/// ��y=mx+n����Բ���̣��γɹ���x�Ķ��η���
		a = 1 + m*m;
		b = -2*xo + 2*m*(n-yo);
		c = xo*xo + (n-yo)*(n-yo) - R*R;

		/// ���ݶ��η������
		delta = b*b - 4*a*c;

		/// �޽���
		if(delta<0)
			return 0;

		/// ����������
		xc1 = ( -1*b + sqrt(delta) ) / (2*a);
		yc1 = m*xc1 + n;
		xc2 = ( -1*b - sqrt(delta) ) / (2*a);
		yc2 = m*xc2 + n;

		/// �жϸ��ĸ���
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

	else	//<  ��y���м���
	{
		/// ��ax+by+c=0ת��Ϊx=my+n����ʽ
		m = -1*pb/pa;
		n = -1*pc/pa;

		/// ��x=my+n����Բ���̣��γɹ���y�Ķ��η���
		a = 1 + m*m;
		b = -2*yo + 2*m*(n-xo);
		c = yo*yo + (n-xo)*(n-xo) - R*R;

		/// ���ݶ��η������
		delta = b*b - 4*a*c;

		/// �޽���
		if(delta<0)
			return 0;

		/// ����������
		yc1 = ( -1*b + sqrt(delta) ) / (2*a);
		xc1 = m*yc1 + n;
		yc2 = ( -1*b - sqrt(delta) ) / (2*a);
		xc2 = m*yc2 + n;

		/// �жϸ��ĸ���
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
//�������ܣ��������λ�á��۲��λ�á��ӳ��Ƕȡ���ȡ�����ɵ����������С
//���룺szy20141023 add 
//eye ���λ��(�д��ŵĸ�˹����)
//view �۲�λ��(�д��ŵĸ�˹����)
//fov_x �ӳ��ݽǶ� ����
//fov_y �ӳ��ݽǶ� ����
//�����result ���������С(�޴��Ÿ�˹����)
////�����������ݣ�Ϊ�˷�������ں�����ڴ����м�⣩
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
	if (eye[0].p_x >= 1000000.0)//���������꺬����
	{
		nDH = int(eye[0].p_x / 1000000);//��ȡ����
		eye[0].p_x = eye[0].p_x - nDH * 1000000;//ȥ������
	}
	else
	{
		AfxMessageBox("���ṩ��λ����Ϣ�޴��ţ������");
		return 0;
	}
	if (view[0].p_x >= 1000000.0)//���������꺬����
	{
		nDH = int(view[0].p_x / 1000000);//��ȡ����
		view[0].p_x = view[0].p_x - nDH * 1000000;//ȥ������
	}
	else
	{
		AfxMessageBox("���ṩ��λ����Ϣ�޴��ţ������");
		return 0;
	}
	//������AB��AH�߶γ�����ȡHB�ĳ��ȣ�����H����A��ͶӰ������������洹ֱ�ĵ�
	double AB = sqrt(pow(eye[0].p_x-view[0].p_x, 2)
		+pow(eye[0].p_y-view[0].p_y, 2)
		+pow(eye[0].p_z-view[0].p_z, 2)
		);
	double AH = eye[0].p_z-view[0].p_z;
	double HB = sqrt(AB*AB-AH*AH);
	//��ȡ��HAB�ĽǶȣ������ݽ�HAB�ĽǶȺ���֪�Ľ�CAB�ĽǶ���ȡ��HAC�ĽǶ�
	double jiao_HAB = atan(HB/AH);
	double jiao_HAC = jiao_HAB-fov_x/2;
	double HC = AH * tan(jiao_HAC);
	//�õ�BC��BD����
	double BC = HB - HC;
	double HD = AH * tan(jiao_HAC+fov_x);
	double BD = HD -HB;

	szypolygen_GS C,E,F,D,G,N;
	szypolygen_DD c,e,f,d,g,n;
	//��ȡC��ĸ�˹����
	C.p_x = HC/HB*(view[0].p_x-eye[0].p_x)+eye[0].p_x;
	C.p_y = HC/HB*(view[0].p_y-eye[0].p_y)+eye[0].p_y;
	C.p_z =view[0].p_z;
	//��ȡAC���߶γ���
	double AC = sqrt(pow(eye[0].p_x-C.p_x, 2)
		+pow(eye[0].p_y-C.p_y, 2)
		+pow(eye[0].p_z-C.p_z, 2)
		);
	//����C��ĸ�˹������ȡE��F��ľ�γ������
	double EC = AC*tan(fov_y/2);
	double jiao_CB = atan((view[0].p_y-C.p_y)/(view[0].p_x-C.p_x));
	//��ȡC�㵽E�㡢F��ĽǶ�
	double jiao_CE = jiao_CB*180/qxwPI+90;
	double jiao_CF = jiao_CB*180/qxwPI-90;
	//ת��ǰ�Ƚ�C��ĸ�˹����ת��Ϊ��γ������
	XYtoBL(C.p_y,C.p_x,nDH*6-3,c.DD_y,c.DD_x);
	calc_coord_by_distance(c.DD_x, c.DD_y, EC/1000, jiao_CE, e.DD_x, e.DD_y, angle21);
	calc_coord_by_distance(c.DD_x, c.DD_y, EC/1000, jiao_CF, f.DD_x, f.DD_y, angle21);

	//��ȡD��ĸ�˹����
	D.p_x = HD/HB*(view[0].p_x-eye[0].p_x)+eye[0].p_x;
	D.p_y = HD/HB*(view[0].p_y-eye[0].p_y)+eye[0].p_y;
	D.p_z =view[0].p_z;
	//��ȡAD���߶γ���
	double AD = sqrt(pow(eye[0].p_x-D.p_x, 2)
		+pow(eye[0].p_y-D.p_y, 2)
		+pow(eye[0].p_z-D.p_z, 2)
		);
	//����D��ĸ�˹������ȡG��N��ľ�γ������
	double DG = AD*tan(fov_y/2);
	
	//��ȡD�㵽G�㡢N��ĽǶ�
	double jiao_DG = jiao_CB*180/qxwPI+90;
	double jiao_DN = jiao_CB*180/qxwPI-90;
	//ת��ǰ�Ƚ�D��ĸ�˹����ת��Ϊ��γ������
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