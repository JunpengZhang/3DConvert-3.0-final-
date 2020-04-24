#include "stdafx.h"
#include "MpGlobalFunc.h"
#include "PolygonBulging.h"
#include <iostream>
#include <fstream>

//����2ά�ڴ�ռ�
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

//����3ά�ڴ�ռ�
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

//����4ά�ڴ�ռ�
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

//ɾ���ض��ļ��������е��ļ�����ɾ��Ŀ¼


/***************************************
���� ��ɾ��Ŀ¼�Լ�Ŀ¼�е��ļ�����Ŀ¼��
�㷨���ݹ����
***************************************/

//���㻭�ߣ��������˵㸡�����꣬���½ǵ�ĸ�������ͷֱ��ʣ�����������߶ε�����
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

//���㻭�ߣ��������˵����꣬����������߶ε�����
int PointToLine(int x1,int y1,int x2,int y2,int &NumPoi,int * (&pX),int *(&pY))
{
	double pa_a,pa_b,pa_c;
	int i,delta,curX,curY;
	
	if(pX!=NULL||pY!=NULL)
		return -1;
	
	//�����һ����
	if(x1==x2&&y1==y2)
	{
		pX = new int[1];
		pY = new int[1];
		pX[0] = x1;
		pY[0] = y1;
		NumPoi = 1;
		return 1;
	}
	
	//���ֱ�߷���
	PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c);

	//��ɨ��
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

//���ݶ˵��������ֱ�߱�׼����
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

//�������㣨��1����2�����꣬�����1����2�ķ���xΪ������yΪ����
int TwoPoint2Angle(double x1,double y1,double x2,double y2,int bHudu,double &AngleRes)
{
	double x2subx1,y2suby1;

	x2subx1 = x2-x1;
	y2suby1 = y2-y1;
	
	if( fabs(x2subx1)<qxwLimMinValue )
		AngleRes = (y2>y1) ? qxwPI_DIV_2: -1*qxwPI_DIV_2 ;
	else if( fabs(y2suby1)<qxwLimMinValue )
		AngleRes = (x2>x1) ? 0.0 : qxwPI ;
	else if( x2subx1>0 && y2suby1>0 )		//��һ����
		AngleRes = atan(y2suby1/x2subx1);
	else if( x2subx1<0 && y2suby1>0 )		//�ڶ�����
		AngleRes = qxwPI + atan(y2suby1/x2subx1);
	else if( x2subx1<0 && y2suby1<0 )		//��������
		AngleRes = atan(y2suby1/x2subx1)-qxwPI;
	else if( x2subx1>0 && y2suby1<0 )		//��������
		AngleRes = atan(y2suby1/x2subx1);

	if(bHudu==0)
		AngleRes = AngleRes*180/qxwPI;

	return 1;
}

/**** �����ߺ��߶εĽǵ���Ŀ
�������߶˵�����ͷ����(����Ϊ�㣬��ʱ��Ϊ��)�������߶ζ˵�����
�����ϢRes˵����
	1��������ֱ�߲��غϣ������߶β��ཻ
	2��������ֱ�߲��غϣ���������ͨ��
	3��������ֱ�߲��غϣ����������ߵĶ˵�
	4��������ֱ�߲��غϣ��������߶εĶ˵�
	5��������ֱ�߲��غϣ��������߶������ߵĶ˵�
	6��������ֱ���غϣ����߶β���������
	7��������ֱ���غϣ����߶���������
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
					   )
{
	double pa_a,pa_b,pa_c;					//ֱ�߷��̲���
	double l_xmin,l_xmax,l_ymin,l_ymax;		//�߶����귶Χ
	double r_xt,r_yt,r_l;					//��ʱ������
	double xPoi,yPoi;						//��������
	double dis1,dis2,dismax,disPoi;
	int bPoiR;								//Ϊ���߶˵�
	int bPoiL;								//Ϊ�߶ζ˵�
	
	bPoiR = bPoiL = 0;
	
	//���߶����귶Χ
	l_xmin = (l_x1<l_x2) ? l_x1 : l_x2 ;
	l_xmax = (l_x1<l_x2) ? l_x2 : l_x1 ;
	l_ymin = (l_y1<l_y2) ? l_y1 : l_y2 ;
	l_ymax = (l_y1<l_y2) ? l_y2 : l_y1 ;
	
	//��ֱ�߷���
	PointToLineEquation( l_x1, l_y1, l_x2, l_y2, pa_a, pa_b, pa_c );

	//ȡ�����Ͼඥ��10������ʱ��
	r_l = 10.0;
	r_xt = r_x0 + r_l*cos(r_dir);
	r_yt = r_y0 + r_l*sin(r_dir);

	//�ж��Ƿ��غϣ���������㶼��ֱ���ϣ��غ�
	if( fabs(pa_a*r_x0+pa_b*r_y0+pa_c) < qxwLimMinValue &&
		fabs(pa_a*r_xt+pa_b*r_yt+pa_c) < qxwLimMinValue )
	{
		//���߶����˵㵽�������ľ���
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

		//�ж��߶��Ƿ��������ϣ�������������С������������
		if(dismax<0)
			Res = 6;
		else
			Res = 7;
		return 1;		
	}

	//��ʱ�������߶ο϶����غϣ��󽻵�
	disPoi = pa_a*cos(r_dir) + pa_b*sin(r_dir) ;
	if( fabs(disPoi)<qxwLimMinValue )	//��ʱƽ��
	{
		Res = 1;	return 1;
	}

	//��ʱ����һ����ֱ���ཻ���󽻵�
	disPoi = ( 0 - pa_c - pa_a*r_x0 - pa_b*r_y0 )/disPoi ;
	xPoi = r_x0+disPoi*cos(r_dir);
	yPoi = r_y0+disPoi*sin(r_dir);

	//����������߶ζ˵�֮һ
	if( ( fabs(xPoi-l_x1)<qxwLimMinValue && fabs(yPoi-l_y1)<qxwLimMinValue ) ||
		( fabs(xPoi-l_x2)<qxwLimMinValue && fabs(yPoi-l_y2)<qxwLimMinValue ) )
	{
		bPoiL = 1;	
	}
	//����������߶���
	else if( xPoi>=l_xmin && xPoi<=l_xmax && yPoi>=l_ymin && yPoi<=l_ymax )
	{
		bPoiL = 0;
	}
	//��������߶ζ˵��ҽ������겻���߶η�Χ�ڣ����ཻ
	else
	{
		Res = 1;	return 1;
	}
	
	//����������������
	if( fabs(xPoi-r_x0)<qxwLimMinValue && fabs(yPoi-r_y0)<qxwLimMinValue  )
	{
		bPoiR = 1;
	}

	//���е����˵���϶��ཻ���жϽ������
	if( bPoiL==0 && bPoiR==0 )	//��ͨ��
		return 2;
	else if( bPoiR==1 && bPoiL==0 )
		return 3;
	else if( bPoiR==0 && bPoiL==1 )
		return 4;
	else
		return 5;	
}

//�����ַ����Ϸ��Լ�麯��
int StringIsNumber(const char * pObj)
{
	int i;
	int len;
	int posDot;		//С�����λ��
	int numDot;		//С�������Ŀ
	int posMinus;	//���ŵ�λ��
	int numMinus;	//���ŵ���Ŀ
	char maxChar='9';
	char minChar='0';

	len = (int)strlen(pObj);
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
		return -1;
	else if( (numMinus==1) && (posMinus!=0) )
		return -1;

	//�����С���㣬ֻ�ܳ���һ������ǰ�󶼱���������
	if(numDot>1)
		return -1;
	else if( (numDot==1) && ( (posDot==0)||(posDot==len-1) ) )
		return -1;

	//ͳ�Ƴ��˸��ź�С����֮�⣬�Ƿ�ÿ���ַ�����0��9֮��
	for(i=0;i<len;i++)
	{
		if( (pObj[i]=='-') || (pObj[i]=='.') )
			continue;
		if( pObj[i]<minChar || pObj[i]>maxChar )
			return -1;	
	}

	return 1;
}


//���������Ƕȣ���ǶȲ��λ��
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

//��ȡ��ǰʱ�䣬��ת��Ϊ�ַ���
int CreateCurTimeString(char * pObjStr)
{
	time_t curTime;
	struct tm *pCurTime;

	if(pObjStr==NULL)
		return -1;

	time(&curTime);
	pCurTime = localtime(&curTime);
	if (pCurTime == NULL)//20161212 304 ��̬�����������
	{
		cout << "��ȡʱ��ʧ��." << endl;
		return 0;
	}
	sprintf(pObjStr,"%s",asctime(pCurTime));
	pObjStr[strlen(pObjStr)-1]='\0';

	return 1;
}

//��ȡ��ǰʱ�䣬��ת��Ϊ�ַ�����ͬʱ��¼��ǰʱ��
int CreateCurTimeString(time_t &curTime, char * pObjStr)
{
	struct tm *pCurTime;

	if(pObjStr==NULL)
		return -1;

	time(&curTime);
	pCurTime = localtime(&curTime);
	if (pCurTime == NULL)//20161212 304 ��̬�����������
	{
		cout << "��ȡʱ��ʧ��." << endl;
		return 0;
	}
	sprintf(pObjStr,"%d-%02d-%02d %02d:%02d:%02d",
		pCurTime->tm_year+1900,pCurTime->tm_mon+1,pCurTime->tm_mday,
		pCurTime->tm_hour,pCurTime->tm_min,pCurTime->tm_sec);
	
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

int InterLinear2D( double x1, double x2, double y1, double y2, double v11, double v12, 
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

int InterLinear3D( double x1, double x2, double y1, double y2, double z1, double z2, 
				  double v111, double v121, double v211, double v221, 
				  double v112, double v122, double v212, double v222, 
				  double x, double y, double z, double &v )
{
	double v1,v2;
	int RV;
	
	//���ȶ�x��y���в�ֵ
	RV = InterLinear2D( x1, x2, y1, y2, v111, v121, v211, v221 , x, y, v1 );
	if(RV==-1)	return -1;

	RV = InterLinear2D( x1, x2, y1, y2, v112, v122, v212, v222 , x, y, v2 );
	if(RV==-1)	return -1;

	//�ٶ�z���в�ֵ
	RV = InterLinear1D( z1, v1, z2, v2, z, v );
	if(RV==-1)	return -1;

	return 1;
}


/*****************************
�жϵ����߶εĹ�ϵ��

  -1,����̫���������м���
  1����Ϊ�߶εĶ˵�
  2����Ϊ�߶��ڵĵ�
  3������ֱ���ϣ����߶��⣬
  4���㲻��ֱ����
******************************/
int Cal_PoiLine_Relation( double poiX, double poiY, double x1, double y1, 
					   double x2, double y2, double minV=1e-5 )
{
	double pa_a,pa_b,pa_c;
	double dis;
	int bOnLine,bPort,bInLine;

	//�������̫���������м���ֱ�ӷ���
	if( sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) ) < minV )
		return -1;

	bOnLine = bPort = bInLine = 0;		//��ʼ����־��
	PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c);	//��ֱ�߷���
	dis = pa_a*poiX + pa_b*poiY + pa_c;

	//�жϵ��Ƿ���ֱ����
	if( fabs(dis)<minV )
		bOnLine = 1 ;

	//�жϵ��Ƿ����߶ζ˵�
	if( 
		( fabs(poiX-x1)<minV && fabs(poiY-y1)<minV ) ||
		( fabs(poiX-x2)<minV && fabs(poiY-y2)<minV )
		)
		bPort = 1;
	
	//�жϵ��Ƿ����߶����˵����귶Χ��
	if( fabs(pa_a-1)>=fabs(pa_b-1) )	//б��С��1�����ŵ��߶�
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

/*****************************
�ж��߶����߶εĹ�ϵ��

 Type = 0���߶����߶��غ�
		1���߶����߶�ƽ��
		2���߶����߶�����ֱ���ཻ
 Line1Type =	0���������߶�1��
				1������Ϊ�߶�1�˵�
				2���������߶�1��
 Line2Type =	0���������߶�2��
				1������Ϊ�߶�2�˵�
				2���������߶�2��
  ����ֵ 
	1���ɹ�
	-1���޷��жϹ�ϵ�����ܳ����������������̫��
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

	//�ֱ��ж������߶ε������˵��Ƿ�̫����̫������ֱ�ӷ���
	dis1 = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
	dis2 = sqrt( (x3-x4)*(x3-x4) + (y3-y4)*(y3-y4) );
	if( dis1<=minV || dis2<=minV )
		return -1;

	//�������߶�ֱ�߷���
	PointToLineEquation(x1,y1,x2,y2,a1,b1,c1);
	PointToLineEquation(x3,y3,x4,y4,a2,b2,c2);

	//ȡ�߶�1���е�
	x0 = (x1+x2)/2;
	y0 = (y1+y2)/2;

	//�����߶�1�����㵽�߶�2����ֱ�ߵľ���
	dis1 = fabs(a2*x1+b2*y1+c2) / sqrt(a2*a2+b2*b2) ;
	dis2 = fabs(a2*x2+b2*y2+c2) / sqrt(a2*a2+b2*b2) ;
	dis3 = fabs(a2*x0+b2*y0+c2) / sqrt(a2*a2+b2*b2) ;

	//�ж�������ֱ�ߵĹ�ϵ
	if( fabs(dis1-dis2)<=minV && fabs(dis2-dis3)<=minV && fabs(dis1)<=minV )
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
		( fabs(a1-1)>=fabs(b1-1) && croX>minX && croX<maxX ) ||
		( fabs(a1-1)<fabs(b1-1) && croY>minY && croY<maxY )
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
�ж��������߶εĹ�ϵ�������ö˵��뷽���ʾ��
(r_x,r_y)Ϊ�˵����꣬r_dirΪ���򣬶��㣬��ʱ��������λ��
�߶������˵��ʾ��(x1,y1)(x2,y2)Ϊ�˵�����

 type = 0���������߶��غ�
		1���������߶�ƽ��
		2���������߶�����ֱ���ཻ
 RadType =	0�����㲻��������
			1������Ϊ���߶˵�
			2��������������
 LineType = 0���������߶���
			1������Ϊ�߶ζ˵�
			2���������߶���
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

	ang = r_dir/180*qxwPI;	//ת��Ϊ����

	//����������һ��������루len��ȡ����
	rx1 = r_x + len*cos( ang );
	ry1 = r_y + len*sin( ang );
	rx2 = rx1 + len*cos( ang );
	ry2 = ry1 + len*sin( ang );

	//�����߶�ֱ�߷���
	PointToLineEquation(x1,y1,x2,y2,pa_a,pa_b,pa_c);

	//��������������ֱ�ߵľ���
	dis1 = fabs( r_x*pa_a + r_y*pa_b + pa_c ) / sqrt( pa_a*pa_a + pa_b*pa_b );
	dis2 = fabs( rx1*pa_a + ry1*pa_b + pa_c ) / sqrt( pa_a*pa_a + pa_b*pa_b );
	dis3 = fabs( rx2*pa_a + ry2*pa_b + pa_c ) / sqrt( pa_a*pa_a + pa_b*pa_b );

	//�ж�������ֱ�ߵĹ�ϵ
	if( fabs(dis1-dis2)<=minV && fabs(dis2-dis3)<=minV && fabs(dis1)<=minV )
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

	//���㽻������
	crossX = ( pa_b*r_c - r_b*pa_c ) / ( pa_a*r_b - r_a*pa_b );
	crossY = ( r_c*pa_a - pa_c*r_a ) / ( pa_b*r_a - r_b*pa_a );

	//�жϽ��������ߵĹ�ϵ
	dis1 = sqrt( (crossX-r_x)*(crossX-r_x) + (crossY-r_y)*(crossY-r_y) );
	sinA = (crossY-r_y)/dis1;
	cosA = (crossX-r_x)/dis1;
	if( dis1<minV )
		RadType = 1 ;
	else if( sin(ang)*sinA>=0 && cos(ang)*cosA>=0 )
		RadType = 2 ; 
	else
		RadType = 0;

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
		( fabs(pa_a-1)>=fabs(pa_b-1) && crossX>minX && crossX<maxX ) ||
		( fabs(pa_a-1)<fabs(pa_b-1) && crossY>minY && crossY<maxY )
		)
		LineType = 2;
	else
		LineType = 0;
	
	return 1;
}



/**********
�жϵ������εĹ�ϵ���ڲ������ⲿ��

  0���ڲ�
  1���ⲿ
  -1�����������߱�����δ����Ч�жϳ����߹�ϵ
******************************/
int Cal_PoiPolygon_Relation( int NumPoi, double *pX, double *pY, 
						  double PoiX, double PoiY, double minV )
{
	int CroDotNumL,CroDotNumR;	//���ҽ������
	double Dir0=0;				//��ʼ����0��
	double DirDlt=20;			//���򲽳���45��
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

	//�����жϵ��Ƿ��ڶ���εı��ϣ����������һһ���ϣ�
	//˵�����ڶ�����ڲ�
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

		//�жϵ��Ƿ��ڵ�ǰ�ı��ϣ�����ڣ���Ϊ���ڶ�����ڲ���ֱ�ӷ���
		RV = Cal_PoiLine_Relation( PoiX, PoiY, x1, y1, x2, y2, minV );
		if( RV==1 || RV==2 )
			return 0;		
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
	angView = atan2( minLen/2, maxDis );
	angView *= 2;
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
		if( (CroDotNumL+CroDotNumR)%2==0 && (CroDotNumL+CroDotNumR)>0 )
		{
			if( CroDotNumL%2==0 && CroDotNumR%2==0 )
				evidenceOut ++;
			else
				evidenceIn ++;
		}

		//����ڲ����ⲿ��֤�ݴ������ڵ���2���˳�ѭ��
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
�ж����������֮��Ĺ�ϵ

  -1���������
  0������
  1���ཻ
  2��1��2�ڲ�
  3��2��1�ڲ�
  4������
  5���غ�
**************************/
int Cal_TwoPolygon_Relation( int NP1, double *pX1, double *pY1,
							int NP2, double *pX2, double *pY2, double minV )
{
	double curX,curY;
	int i,j;
	int RV;
	int eviOut,eviIn;
	int relation12,relation21;		//1���ڲ���2���ⲿ��3���ⶼ��
	double x1,x2,x3,x4,y1,y2,y3,y4;
	int Type,Line1Type,Line2Type;

	//���ȶԵ�һ������ε�ÿ��������ڶ�������εĹ�ϵ�����ж�
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

	//�Եڶ�������ε�ÿ���������һ������εĹ�ϵ�����ж�
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

	//�Ա���������α��ཻ��������жϣ�
	//���������ཻ����ô�������һ���ֱ�������һ������Է��ı��߶��ڽ�
	//������Ҳ����
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

			//����б��ཻ��������һ���ཻ
			if( Type==2 && Line1Type==2 && Line2Type==2 )
			{
				return 1;
			}
		}
	}

	//��������������ϵ�ж�������εĹ�ϵ
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
	else		//��ʱrelation12==3 && relation21==3�������ཻ�ͽ��ڣ����㷨��ʱ�޷��ж�
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

	//���������ܼ��㳬��3.29��sigma��ֵ
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

	//���ڸ�����ת��Ϊ�������м���
	if(para<0)
	{
		para = fabs(para);
		bNeg = 1;
	}

	//��������������
	SNRow = (int)(para*10);
	SNCol = (int)(para*100-SNRow*10);

	//�������õ���3.29�������ע�⸡��������ֱ���ж����
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

	//���Բ�ֵ
	Val = v1 + (para-x1)/(x2-x1)*(v2-v1);

	if(bNeg==1)
	{
		Val = 1-Val;
	}

	return 1;
}


///////////////////////////////////////////////////////////////////
////////////////  �����Ǿ����������غ���/////////////////////////

//�������
int MatMul(double * mat1,int h1,int w1,double *mat2,int h2,int w2,double *matOut)
{
	int i,j,k,ho,wo;
	double val;

	double *mat11 = new double[h1*w1];
	double *mat21 = new double[h2*w2];

	memcpy(mat11,mat1,sizeof(double)*h1*w1);
	memcpy(mat21,mat2,sizeof(double)*h2*w2);
	
	//���ǰ�߾�������ں��߾���ߣ�������
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

//�жϵ��Ƿ����ߵ��Ҳ�
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

	/// ���������Ϣ��
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

		/// ͳ��ͼ�η�Χ��ȷ����ʾͼ�εķ�Χ
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
		//if (pF == NULL)//20161212 304 ��̬�����������
		//{
		//	cout << "�ļ���ʧ��." << endl;
		//	return 0;
		//}
		///// �����������
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
		
	//�����жϴ˶�����Ƿ�͹����Σ�������ǣ��ҵ�����
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

	//������ǰ�����Σ�ֱ�ӽ��˶���μ������������б���
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
			/// �������
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

	//ȷ����ǰ�ߵ���ֹ���
	SN1 = SNCon;
	SN2 = SNCon+1;
	if(SNCon==Num-1)
		SN2 = 0;
	
	//�ָ��һ�������Σ���Ϊһ����������
	for( bFind=0,i=0; i<Num; i++ )
	{
		if( i==SN1 || i==SN2 )
			continue;

		bRight = Point_OnRight_Line( pX[SN1], pY[SN1], pX[SN2], pY[SN2], pX[i], pY[i] );
		if( bRight!=1 )
			continue;

		//ѡ����i��Ϊ�����ε�һ�����㣬SN1��SN2��i����һ��������
		//�ж϶���εı��Ƿ��������ε��������ཻ
		for( bValid=1,j=0; j<Num; j++ )
		{
			SN3 = j;
			SN4 = j+1;
			if(j==Num-1)
				SN4 = 0;

			/// �жϵ�ǰ���Ƿ����������������غ�
			/// ����غϲ������ж�

			/// ��(SN1,i) �Ƿ����(SN3,SN4)�غ�
			if( ( SN1==SN3 && i==SN4 ) || ( SN1==SN4 && i==SN3 ) )
				continue;

			/// ��(SN2,i) �Ƿ����(SN3,SN4)�غ�
			if( ( SN2==SN3 && i==SN4 ) || ( SN2==SN4 && i==SN3 ) )
				continue;

			
//			//�����ǰ�߶�������������һһ�������غϣ������м���
//			if( 
//				SN3==SN1 || SN3==SN2 || SN3==i || 
//				SN4==SN1 || SN4==SN2 || SN4==i  
//				)
//				continue;

			/// �ж��ཻ���������ཻ�ҽ������߶��м䣬˵���������β��Ƕ�����ڲ�
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

		//�����Ч��˵���ҵ���Ŀ��������
		if(bValid==1)
		{
			bFind = 1;
			SNObj = i;
			break;
		}
	}

	//Ӧ���ܹ��ҵ�Ŀ�������Σ�����������
	if(bFind==0)
		return -1;

	//��Ŀ��������(SN1��SN2��SBObj����)�����������
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

	/// ����м�������
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
	
	
	//��ʣ��Ķ����Ϊ���룬�ظ����ȱ������������ܲ��������¶����
	
	//�����һ�������
	if( SNObj<SN2 )
		Num1 = SNObj+Num-SN2+1;
	else
		Num1 = SNObj-SN2+1;
	if( Num1>2 )
	{
		//�����µĶ����
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

		/*	/// �����һ�������
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
		

		//�ݹ���ñ��������ж���β��
		RV = ConcaveToBulge_iteration( Num1, pNew1X, pNew1Y, vecPoly );
		delete [] pNew1X;
		delete [] pNew1Y;

		if(RV!=1)
			return -1;	
	}

	//����ڶ��������
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
		
		/// ����ڶ��������
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

		/// ��ͼ����
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
	
	
	//�ݹ���ñ��������ж���β��
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
	
	//�����жϴ˶�����Ƿ�͹����Σ�������ǣ��ҵ�����
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

	//������ǰ�����Σ�ֱ�ӽ��˶���μ������������б���
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

	//ȷ����ǰ�ߵ���ֹ���
	SN1 = SNCon;
	SN2 = SNCon+1;
	if(SNCon==Num-1)
		SN2 = 0;
	
	//�ָ��һ�������Σ���Ϊһ����������
	for( bFind=0,i=0; i<Num; i++ )
	{
		if( i==SN1 || i==SN2 )
			continue;

		bRight = Point_OnRight_Line( pX[SN1], pY[SN1], pX[SN2], pY[SN2], pX[i], pY[i] );
		if( bRight!=1 )
			continue;

		//ѡ����i��Ϊ�����ε�һ�����㣬SN1��SN2��i����һ��������
		//�ж϶���εı��Ƿ��������ε��������ཻ
		for( bValid=1,j=0; j<Num; j++ )
		{
			SN3 = j;
			SN4 = j+1;
			if(j==Num-1)
				SN4 = 0;

			//�����ǰ�߶�������������һһ�������غϣ������м���
			if( 
				SN3==SN1 || SN3==SN2 || SN3==i || 
				SN4==SN1 || SN4==SN2 || SN4==i  
				)
				continue;

			//�ж��ཻ���
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

		//�����Ч��˵���ҵ���Ŀ��������
		if(bValid==1)
		{
			bFind = 1;
			SNObj = i;
			break;
		}
	}

	//Ӧ���ܹ��ҵ�Ŀ�������Σ�����������
	if(bFind==0)
		return -1;

	//��Ŀ��������(SN1��SN2��SBObj����)�����������
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
	
	//��ʣ��Ķ����Ϊ���룬�ظ����ȱ������������ܲ��������¶����
	
	//�����һ�������
	if( SNObj<SN2 )
		Num1 = SNObj+Num-SN2+1;
	else
		Num1 = SNObj-SN2+1;
	if( Num1>2 )
	{
		//�����µĶ����
		pNew1X = new double[Num1];
		pNew1Y = new double[Num1];
		for( SN=0,i=SN2; SN<Num1; SN++,i++ )
		{
			if(i==Num)
				i=0;
			pNew1X[SN] = pX[i];
			pNew1Y[SN] = pY[i];
		}

		//�ݹ���ñ��������ж���β��
		RV = ConcaveToBulge_iteration( Num1, pNew1X, pNew1Y, pOut );
		delete [] pNew1X;
		delete [] pNew1Y;

		if(RV!=1)
			return -1;		
	}

	//����ڶ��������
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
	
	//�ݹ���ñ��������ж���β��
	RV = ConcaveToBulge_iteration( Num2, pNew2X, pNew2Y, pOut );
	delete [] pNew2X;
	delete [] pNew2Y;

	if(RV!=1)
		return -1;
	
	return 1;
}

/****************************************************
���ܣ�
	�����������Σ�������Ϊ����͹����Ρ�

���룺
    �ļ������ļ������ı���Ŷ���νǵ��������˳ʱ��˳����
ÿ���ǵ��x��y���ꡣ

�����
    ����ֽ��д��Ҫ����ı��ļ���
	���δ�Ŷ���θ���n
	�����1 �Ķ��������˳ʱ���Ŷ���x��y����
	�����2 �Ķ��������˳ʱ���Ŷ���x��y����
	������������
	�����n �Ķ��������˳ʱ���Ŷ���x��y����

����ֵ��
	1���ɹ�
	��1��ʧ��
***********************************************************/
int ConcaveToBulgePolygon( char *FN1, char *FN2 )
{
	QBFPolygonList ListPoly;
	QBFPolygon *pHead,*pCur;
	double *pX,*pY;
	int Num,i,j,RV;
	FILE * pF;

	//��ʼ�����������
	pHead = new QBFPolygon;
	pHead->Num = 0;
	pHead->pNext = NULL;
	pHead->pX = NULL;
	pHead->pY = NULL;
	ListPoly.Num = 1;
	ListPoly.pHead = ListPoly.pEnd = pHead;
	
	//�������ļ�
	pF = fopen( FN1, "rt" );
	if( pF==NULL )
	{
		while(pHead)//20161212 304 ��̬�����������
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

	fscanf( pF, "%d", &Num );	//����ǵ����

	//����ռ�
	pX = new double[Num];
	pY = new double[Num];

	//���ζ���ǵ�����
	for( i=0; i<Num; i++ )
	{
		fscanf( pF, "%lf", &(pX[i]) );
		fscanf( pF, "%lf", &(pY[i]) );
	}
	fclose(pF);

	//���õݹ麯����ֶ����
	RV = ConcaveToBulge_iteration( Num, pX, pY, &ListPoly );
	if(RV!=1)
	{
		//�ͷſռ�
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

	//�ͷſռ�
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


//������غ���

/*************************************************
 ���԰�����β��Ϊ͹����εݹ麯���ĺ��� 

���ı��ļ��������νǵ�����Լ�ÿ���ǵ��x��y���꣬Ȼ����в�֣���д���ļ�
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

	////�ͷſռ�
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
