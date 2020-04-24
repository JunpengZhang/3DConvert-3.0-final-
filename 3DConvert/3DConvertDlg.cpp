
// 3DConvertDlg.cpp : ʵ���ļ�
//

#include "stdafx.h"
#include "3DConvert.h"
#include "3DConvertDlg.h"
#include "afxdialogex.h"
#include "math.h"
#include "PolygonBulging.h"
#include "pointclockpre.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>




#ifdef _DEBUG
#define new DEBUG_NEW
#endif
#define PI 3.1415926535897932384626433832795

// ����Ӧ�ó��򡰹��ڡ��˵���� CAboutDlg �Ի���

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// �Ի�������
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

// ʵ��
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CMy3DConvertDlg �Ի���




CMy3DConvertDlg::CMy3DConvertDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(CMy3DConvertDlg::IDD, pParent)
	, In_Path(_T(""))
	, m_Radio(FALSE)
	, Out_Path(_T(""))
	, m_DH(_T(""))
	{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CMy3DConvertDlg::DoDataExchange(CDataExchange* pDX)
{
CDialogEx::DoDataExchange(pDX);
DDX_Text(pDX, IDC_EDIT1, In_Path);
DDX_Text(pDX, IDC_EDIT2, Out_Path);
DDX_Radio(pDX, IDC_RADIO3, m_Radio);
DDX_Text(pDX, IDC_EDIT3, m_DH);

}

BEGIN_MESSAGE_MAP(CMy3DConvertDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON1, &CMy3DConvertDlg::OnBnClickedButton1)
	ON_BN_CLICKED(IDOK, &CMy3DConvertDlg::OnBnClickedOk)
	ON_BN_CLICKED(IDC_BUTTON2, &CMy3DConvertDlg::OnBnClickedButton2)
END_MESSAGE_MAP()


// CMy3DConvertDlg ��Ϣ�������

BOOL CMy3DConvertDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// ��������...���˵�����ӵ�ϵͳ�˵��С�

	// IDM_ABOUTBOX ������ϵͳ���Χ�ڡ�
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// ���ô˶Ի����ͼ�ꡣ��Ӧ�ó��������ڲ��ǶԻ���ʱ����ܽ��Զ�
	//  ִ�д˲���
	SetIcon(m_hIcon, TRUE);			// ���ô�ͼ��
	SetIcon(m_hIcon, FALSE);		// ����Сͼ��

	// TODO: �ڴ���Ӷ���ĳ�ʼ������

	return TRUE;  // ���ǽ��������õ��ؼ������򷵻� TRUE
}

void CMy3DConvertDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// �����Ի��������С����ť������Ҫ����Ĵ���
//  �����Ƹ�ͼ�ꡣ����ʹ���ĵ�/��ͼģ�͵� MFC Ӧ�ó���
//  �⽫�ɿ���Զ���ɡ�

void CMy3DConvertDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // ���ڻ��Ƶ��豸������

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// ʹͼ���ڹ����������о���
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// ����ͼ��
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//���û��϶���С������ʱϵͳ���ô˺���ȡ�ù��
//��ʾ��
HCURSOR CMy3DConvertDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

int CMy3DConvertDlg::OnChangeByteOrder(int indata)
{
	char ss[9];
	char ee[8];
	unsigned long val=unsigned long(indata);
	ultoa(val,ss,16);// ��ʮ�����Ƶ���(val)ת��һ���ַ���(ss)��
	int i;
	int length=strlen(ss);
	if(length!=8)
	{
		for(i=0;i<8-length;i++)
			ee[i]='0';
		for(i=0;i<length;i++)
			ee[i+8-length]=ss[i];
		for(i=0;i<8;i++)
			ss[i]=ee[i];
	}
	////****** ���е���
	int t;
	t=ss[0];ss[0]=ss[6];ss[6]=t;
	t=ss[1];ss[1]=ss[7];ss[7]=t;
	t=ss[2];ss[2]=ss[4];ss[4]=t;
	t=ss[3];ss[3]=ss[5];ss[5]=t;
	////******
	//****** ������ʮ�������� (val) ���ַ��� (ss) �е�ʮ��������ת��ʮ������
	int value=0;
	for(i=0;i<8;i++)
	{
		int k;
		if(ss[i]=='a'||ss[i]=='b'||ss[i]=='c'||ss[i]=='d'||ss[i]=='e'||ss[i]=='f')
			k=10+ss[i]-'a';
		else
			k=ss[i]-'0';
		//printf("k=%d\n",k);
		value=value+int(k*pow(16.0,7-i));
	}
	//printf("value=%d\n",value);
	return(value);
}
//string����ת������
//int_temp ���������
//string_temp �����ַ���
void CMy3DConvertDlg::str2int(int &int_temp,const string &string_temp)  
{  
	stringstream stream(string_temp);  
	stream>>int_temp;  
}

string CMy3DConvertDlg::TrimLeft(string str)
{
	for (string::size_type i=0; i != str.size(); ++i)
	{
		if (!isspace(str[i]))
		{
			return str.substr(i);
		}
	}
	return "";
}

int CMy3DConvertDlg::ReadBuildingsFrom3d(const char * path,_3DFileInfo& _3dFile)
{
	setlocale(LC_ALL,"Chinese-Simplified");
	fstream infile;
	infile.open(path , ios::in);
	if(!infile.is_open()) 
	{
		return -1;
	}
	int numOfLayer,numOfPolygon,numOfPoint;
	string strLine;
	int lineNumber = 1;		// �кţ���һ�е��к���1�����������У�
	string Ln5;
	while (lineNumber <= 6)
	{
		if (lineNumber == 6)
		{
			infile >>numOfLayer;//���ļ�����
		}

		getline(infile,strLine);

		switch(lineNumber)
		{
		case 1:
			_3dFile.versionNum = strLine;
			break;

		case 2: _3dFile.Ln2 = atoi(strLine.c_str());
			break;

		case 3: _3dFile.Meter = strLine;
			break;

		case 5: Ln5 = strLine;
			break;

		case 6: _3dFile.layerNum = numOfLayer;
			break;
		default:
			break;
		}

		lineNumber++;
	}

	//��ȡ������ֵ
	double minX;
	double maxX;
	double minY;
	double maxY;
	double dH;
	double dH1;
	Ln5 = TrimLeft(Ln5);
	sscanf(Ln5.c_str(),"%lf %lf %lf %lf %lf %lf" ,
		&minX, &maxX, &minY, &maxY, &dH, &dH1);
	_3dFile.faceMaxMinGs.push_back(minX);
	_3dFile.faceMaxMinGs.push_back(maxX);
	_3dFile.faceMaxMinGs.push_back(minY);
	_3dFile.faceMaxMinGs.push_back(maxY);
	_3dFile.faceMaxMinGs.push_back(dH);
	_3dFile.faceMaxMinGs.push_back(dH1);

	for (int kk=0; kk<numOfLayer; kk++)
	{
		getline(infile,strLine);//����
		_3dFile.layerName = strLine;

		getline(infile,strLine);//�ò���ʾ���
		_3dFile.layerShow = atoi(strLine.c_str());

		getline(infile,strLine);//�ò�������
		_3dFile.layerLineNum = atoi(strLine.c_str());

		infile >>numOfPolygon;
		_3dFile.modNum = numOfPolygon;//�ܵ�����
		getline(infile,strLine);

		getline(infile,strLine);
		_3dFile.Ln11 = atoi(strLine.c_str());

		getline(infile,strLine);
		_3dFile.Ln12 = atoi(strLine.c_str());

		getline(infile,strLine);
		_3dFile.Ln13 = atoi(strLine.c_str());

		/////20190117zqs--���ļ�����//////
		infile.close();


		//////////////////////////////////
		for(int i = 0; i<numOfPolygon; i ++)
		{
			infile.open("path , ios::add");
			if (!infile.is_open())
			{
				return -1;
			}

			Building tmpBuilding;
			polygon tmpFace;
			getline(infile,strLine);
			strLine = TrimLeft(strLine);
			char buf[100];
			int a,b,c,d,anatomicNum,f,e,g,h,p,q,k,l,m,n,r,s;//ÿ�����һ����ֵ
			double a1,b1,c1;
			int polyid;
			sscanf(strLine.c_str()," %d %d %d %d %d %d %d %lf %d %d %d %d %d %d %d %d %d %d %d %lf %lf %[^\0]s" ,
				&numOfPoint, &polyid, &b, &c, &d,&e,&g,&a1,&anatomicNum,&h,&s,&f,&r,&p,&q,&k,&l,&m,&n,&b1,&c1,buf);

			tmpFace.nodeNum = numOfPoint;//����ܸ���
			tmpFace.id = polyid;//��ID
			tmpFace.PolygonRGBs.iB = b;
			tmpFace.PolygonRGBs.iG = c;
			tmpFace.PolygonRGBs.iR = d;

			tmpFace.faceLn1Col6 = e;
			tmpFace.faceLn1Col7 = g;
			tmpFace.faceLn1Col8 = a1;

			tmpFace.faceLn1Col9 = anatomicNum;
			tmpFace.faceLn1Col10 = h;
			tmpFace.faceLn1Col11 = s;
			tmpFace.faceLn1Col12 = f;
			tmpFace.faceLn1Col13 = r;
			tmpFace.faceLn1Col14 = p;
			tmpFace.faceLn1Col15 = q;
			tmpFace.faceLn1Col16 = k;
			tmpFace.faceLn1Col17 = l;
			tmpFace.faceLn1Col18 = m;
			tmpFace.faceLn1Col19 = n;
			tmpFace.faceLn1Col20 = b1;
			tmpFace.faceLn1Col21 = c1;

			string strTemp = buf;
			string buildId = strTemp.substr(0,strTemp.find("#\\/#"));//��������
			tmpBuilding.id = buildId;

			for (int j=0; j<anatomicNum; j++)//��ȡ�ʷ�����
			{
				infile>>a;
				tmpFace.subFace.push_back(a);
			}

			getline(infile,strLine);

			//��ȡ���������
			tmpFace.points.resize(numOfPoint);
			for(int j=0; j<numOfPoint ;j++)
			{
				infile>>tmpFace.points[j].x >>tmpFace.points[j].y >> tmpFace.points[j].z >> tmpFace.points[j].nodeCol4 >>tmpFace.points[j].nodeCol5;
				getline(infile,strLine);
			}

			//�ж��Ƿ���ڵ�ǰ���������������浽��Ӧ�������±��
			bool isExist = false;
			int index=0;
			int samNum = (int)_3dFile.buildInfo.size();
			for (index; index<samNum; index++)
			{
				if (buildId == _3dFile.buildInfo[index].id)
				{
					isExist = true;
					break;
				}
			}

			if (isExist)
			{
				_3dFile.buildInfo[index].polygons.push_back(tmpFace);
			}
			else
			{
				tmpBuilding.polygons.push_back(tmpFace);
				_3dFile.buildInfo.push_back(tmpBuilding);
			}	

		}
	infile.close();
	}

	//infile.close();
	return 1;
}

//������
int CMy3DConvertDlg::BuildAnilysize_CulFaceOrigin(_3DFileInfo tmp3DFile,CString error_outputpath)
{

	FILE* PF;
    double dMaxDis = 1;//������ֵ �ɵ� ��0.5��ʼ���ֲ����� ������Ҫ������ǰ0��1����
	double dX,dY,dZ;
	//��ȡ���������
	int allBuildNum = (int)tmp3DFile.buildInfo.size();
	PF = fopen(error_outputpath,"w+");
	for (int buildNum=0; buildNum<allBuildNum;buildNum++)
	{
		int polygon_num = tmp3DFile.buildInfo[buildNum].polygons.size();
		for (int polygonNum=0;polygonNum<polygon_num;polygonNum++)
		{
			int nPointNum = tmp3DFile.buildInfo[buildNum].polygons[polygonNum].points.size();
			double *pPointDX = new double[nPointNum];//��̬�����ڴ�
			double *pPointDY = new double[nPointNum];
			double *pPointDZ = new double[nPointNum];
			CCoorDS Origin,  CoorX,  CoorY;
			CoorX.dx = CoorX.dy = CoorX.dz = Origin.dx = Origin.dy = Origin.dz = CoorY.dx= CoorY.dy = CoorY.dz = 0.0;//���
			//��ȡÿһ��������ĵ�һ����x��y��z��Ϣ
			dX = dY = dZ = 0.0;//���
			dX = tmp3DFile.buildInfo[buildNum].polygons[0].points[0].x;
			dY = tmp3DFile.buildInfo[buildNum].polygons[0].points[0].y;
			dZ = tmp3DFile.buildInfo[buildNum].polygons[0].points[0].z;

			//�õ�ÿ����������Ե�һ�����һ��������λ��������Ϣ
			for (int point_num=0;point_num<nPointNum;point_num++)
			{
				double temp_pPointDX = tmp3DFile.buildInfo[buildNum].polygons[polygonNum].points[point_num].x;//�õ������е��X����
				double temp_pPointDY = tmp3DFile.buildInfo[buildNum].polygons[polygonNum].points[point_num].y;
				double temp_pPointDZ = tmp3DFile.buildInfo[buildNum].polygons[polygonNum].points[point_num].z;
 				
 				pPointDX[point_num] = temp_pPointDX - dX;//�õ�ÿ����������Ե�һ�����һ�������Ϣ
 				pPointDY[point_num] = temp_pPointDY - dY;
 				pPointDZ[point_num] = temp_pPointDZ - dZ;
			}
			int nRet = CulFaceOrigin(pPointDX, pPointDY, pPointDZ, nPointNum, Origin, CoorX, CoorY,dMaxDis);
			if (nRet != 1)
			{
				fprintf(PF,"������%s�ĵ�%d���治����\n",tmp3DFile.buildInfo[buildNum].id.c_str(),polygonNum);
				//return nRet;
			}
		}
	}
	fclose(PF);	
	return 1;
}

bool CMy3DConvertDlg::ComputeMFaX(double *X,double *Y,double *Z,int Num,bool IsClockWise,vector<double>&MFX)
{
	//���������С��3�����޷�����
	if (Num<3)
	{
		return false;
	}

	//����ǰ3�����ж�ƽ���Ƿ���ƴ�ֱXYƽ��,���ж��䷨���Ƿ�ֱ��Z��
	double tmpFX[3];
	tmpFX[0]=(Y[1]-Y[0])*(Z[2]-Z[1])-(Y[2]-Y[1])*(Z[1]-Z[0]);
	tmpFX[1]=(Z[1]-Z[0])*(X[2]-X[1])-(Z[2]-Z[1])*(X[1]-X[0]);
	tmpFX[2]=(X[1]-X[0])*(Y[2]-Y[1])-(X[2]-X[1])*(Y[1]-Y[0]);

	bool IsFaceVerXY=false;
	//ƽ����Z��
	if(abs(tmpFX[2])<=0.000001)
	{
		IsFaceVerXY=true;
	}

	//������ͶӰ
	double *TYX=new double[Num];
	double *TYY=new double[Num];
	if(IsFaceVerXY)
	{
		//X�ᴹֱ������
		if(tmpFX[1]<=0.000001)
		{
			//��YZƽ��ͶӰ
			for(int i=0;i<Num;i++)
			{
				TYX[i]=Y[i];
				TYY[i]=Z[i];
			}
		}


		else
		{
			//��XZƽ��ͶӰ
			for(int i=0;i<Num;i++)
			{
				TYX[i]=X[i];
				TYY[i]=Z[i];
			}
		}
	}

	else
	{
		//��XYƽ��ͶӰ
		for(int i=0;i<Num;i++)
		{
			TYX[i]=X[i];
			TYY[i]=Y[i];
		}
	}

	bool Flag=false;
	//�жϵ�2�������Ƿ񰼶���
	double Value0=(TYY[2]-TYY[1])*(TYX[1]-TYX[0])-(TYY[1]-TYY[0])*(TYX[2]-TYX[1]);

	//�ҳ�͹����
	int ConvexP;
	FindConvexP(Num,TYX,TYY,ConvexP);
	int beforeP,afterP;
	beforeP=ConvexP==0?Num-1:ConvexP-1;
	afterP=ConvexP==Num-1?0:ConvexP+1;
	double tmpValue=(TYY[afterP]-TYY[ConvexP])*(TYX[ConvexP]-TYX[beforeP])-(TYY[ConvexP]-TYY[beforeP])*(TYX[afterP]-TYX[ConvexP]);
	if(tmpValue*Value0<0)
	{
		Flag=true;
	}

	if(TYX)
		delete []TYX;
	if(TYY)
		delete []TYY;
	//���㷨��
	MFX.resize(3);
	if (Flag)
	{
		//�а����㣬����ȡ��
		MFX[0]=-tmpFX[0];
		MFX[1]=-tmpFX[1];
		MFX[2]=-tmpFX[2];

	}
	else
	{
		//�ް�����
		MFX[0]=tmpFX[0];
		MFX[1]=tmpFX[1];
		MFX[2]=tmpFX[2];
	}

	//����ǵ㰴˳ʱ�����������ȡ��
	if(IsClockWise)
	{
		MFX[0]=-MFX[0];
		MFX[1]=-MFX[1];
		MFX[2]=-MFX[2];
	}

	vector<double> tmpMFX;
	tmpMFX.resize(3);
	//�����һ��
	tmpMFX[0]=MFX[0]/sqrt(MFX[0]*MFX[0]+MFX[1]*MFX[1]+MFX[2]*MFX[2]);
	tmpMFX[1]=MFX[1]/sqrt(MFX[0]*MFX[0]+MFX[1]*MFX[1]+MFX[2]*MFX[2]);
	tmpMFX[2]=MFX[2]/sqrt(MFX[0]*MFX[0]+MFX[1]*MFX[1]+MFX[2]*MFX[2]);
	MFX=tmpMFX;
	return true;
}

bool CMy3DConvertDlg::FindConvexP( int Num,double *X,double *Y,int &ConvexP )
{
	if(Num<=0)
	{
		return false;
	}
	//�ҳ�����Y
	double MaxY=Y[0];
	for(int i=0;i<Num;i++)
	{
		MaxY=MaxY>Y[i]?MaxY:Y[i];
	}

	//������Y��ֻһ�����ҳ�����Y��Ӧ������X

	vector<int>Point;
	for (int i=0;i<Num;i++)
	{
		if(abs(Y[i]-MaxY)<0.000001)
		{
			Point.push_back(i);
		}
	}

	int Max=0;
	double MaxX=X[Point[0]];
	for(int i=0;i<int(Point.size());i++)
	{
		if(MaxX<X[Point[i]])
		{
			Max=Point[i];
		}
	}

	ConvexP=Max;
	return true;
}

int CMy3DConvertDlg::CulFaceOrigin(double *pX, double *pY, double *pZ, int nNum, CCoorDS &Origin, CCoorDS &CoorX, CCoorDS &CoorY, double dMinDis)
{
	double x, y, h, xmin, ymin, xmax, ymax, hb;
	hb = 99999;
	xmin = ymin = 9999999999;
	xmax = ymax = -999999999;
	int j, SNb;
	SNb = -1;
	for (j=0; j<nNum; j++)
	{
		x = pX[j];
		y = pY[j];
		h = pZ[j];
		xmin = min(xmin,x);
		xmax = max(xmax,x);
		ymin = min(ymin,y);
		ymax = max(ymax,y);

		//ͳ�ƽǵ�������͵����
		if( h<hb )
		{
			SNb = j;
			hb = h;
		}	
	}

	//////////////////////////////////////////////////////////////////////////
	//����Ϊ��ǰ�潨��ƽ������ϵ
	//��ȷ��ԭ�����꣬X�����Y����
	//////////////////////////////////////////////////////////////////////////

	//ȷ����Ͷ����ǰ�������
	int SN1 = (SNb==0) ? nNum-1 : SNb-1;
	//ȷ����Ͷ���ĺ󶥵����
	int SN2 = (SNb==(nNum-1)) ? 0 : SNb+1;

	//������������ʸ��
	C3DVector v1, v2, v, vX, vY, vo;
	v1.Set(pX[SNb]-pX[SN1],pY[SNb]-pY[SN1],pZ[SNb]-pZ[SN1]);
	v2.Set(pX[SN2]-pX[SNb],pY[SN2]-pY[SNb],pZ[SN2]-pZ[SNb]);

	//ʸ����˵õ��淨ʸ
	double xo, yo, ho;
	// 	v = (v1<v2)*(-1);  //ֱ�Ӳ��δ���ǰ���������
	// 	v.Normalize(1);
	// 	C3DVector VOld = v;
	//  ���ǰ�����ε��淨����㷽��
	//  true /*�����������涥�㰴��˳ʱ��˳�����ע�����϶��´�ֱ�۲�ˮƽ��ʱ������˳���ŷ���˳ʱ��˳��*/
	vector<double> MFX;
	ComputeMFaX(pX,pY,pZ,nNum,true,MFX);
	if (MFX.size() == 0)
	{
		return -1;
	}
	v.x = MFX[0];
	v.y = MFX[1];
	v.z = MFX[2];

	xo = sqrt(v.x*v.x+v.y*v.y);
	ho = fabs(v.z);

	//�������ƽ��
	//if( xo<ho && xo/ho<0.00001 )// Sam  2013-1-30 16:24:16
	if( xo<ho && xo/ho<0.005 )
	{
		//������Ӿ������½ǵ�Ϊԭ��
		xo = xmin;
		yo = ymin;
		ho = pZ[0];
		//����ΪX��������ΪY����������ϵ
		vX.Set(1,0,0);
		vY.Set(0,1,0);			
	}
	//����治��ƽ��
	else
	{
		//������Ͷ���Ϊԭ��
		xo = pX[SNb];
		yo = pY[SNb];
		ho = pZ[SNb];

		//ȷ��X����
		vX.z = 0;
		vX.x = -1*v.y;
		vX.y = v.x;
		vX.Normalize(1);

		//ȷ��Y����
		vY = v;
		vY = vY<vX;
		vY.Normalize(1);

		//ȷ������СX����
		xmin = 99999999;
		for( j=0; j<nNum; j++ )
		{
			v1.x = pX[j]-xo;
			v1.y = pY[j]-yo;
			v1.z = pZ[j]-ho;
			x = v1*vX;
			y = v1*vY;
			xmin = min(xmin,x);
		}

		//��ԭ����X������ƽ��xmin
		vo.Set(xo,yo,ho);
		vo = vo + vX*xmin;
		xo = vo.x;
		yo = vo.y;
		ho = vo.z;
	}


	//Ϊ�����ԭ�������ֵ
	Origin.dx = xo;
	Origin.dy = yo;
	Origin.dz = ho;

	CoorX.dx = vX.x;
	CoorX.dy = vX.y;
	CoorX.dz = vX.z;

	CoorY.dx = vY.x;
	CoorY.dy = vY.y;
	CoorY.dz = vY.z;


	v.Normalize(1);
	// �������ڵ㹲����
	for (int i=0; i<nNum; i++)
	{
		double dx = pX[i] - Origin.dx;
		double dy = pY[i] - Origin.dy;
		double dz = pZ[i] - Origin.dz;

		C3DVector VecFaceOriginToPT(dx,dy,dz);
		C3DVector FaShiOfFace = v;

		double dis = fabs(VecFaceOriginToPT.x* FaShiOfFace.x + VecFaceOriginToPT.y* FaShiOfFace.y + VecFaceOriginToPT.z* FaShiOfFace.z);
		if (dis > dMinDis)
		{
			return -1;
		}
	}

	return 1;
}
//�������������Сֵ����
//Num ������
//Array ��������
//MaxMin ��������Сֵ���飬MaxMin[0]Ϊ��Сֵ��MaxMin[1]Ϊ���ֵ
void CMy3DConvertDlg::findMinAndMax(int Num,double *Array,double *MaxMin)
	{
	double tempMax,tempMin;
	int i;
	tempMax = Array[0];
	tempMin = Array[0];
	for(i=1;i<Num;i++)
		{
		if(Array[i]>tempMax) tempMax = Array[i];
		if(Array[i]<tempMin) tempMin = Array[i];
		}
	MaxMin[0] = tempMin;
	MaxMin[1] = tempMax;
	return ;
	}

int CMy3DConvertDlg::ClockWiseXY(int Num,double *PointsX,double *PointsY)
{
	int flag;
	for (int a=0;a<Num;a++)//������ת��Ϊ�����������ʽ
	{
		m_vctClockX.push_back(PointsX[a]);
		m_vctClockY.push_back(PointsY[a]);
	}
	flag = ClockWisePro(Num, m_vctClockX, m_vctClockY, m_vctClockRX, m_vctClockRY);//���нǵ����򣬽��ǵ㰴˳ʱ��˳�����
	for (int b=0;b<Num;b++)//������ת����͹����λ��ֺ��������ʽ
	{
		PointsX[b] = m_vctClockRX[b];
		PointsY[b] = m_vctClockRY[b];
	}
	m_vctClockRX.clear();
	m_vctClockRY.clear();
	m_vctClockX.clear();
	m_vctClockY.clear();
	return flag;
}

int CMy3DConvertDlg::ClockWiseXZ(int Num,double *PointsX,double *PointsZ)
{
	int flag;
	for (int a=0;a<Num;a++)//������ת��Ϊ�����������ʽ
	{
		m_vctClockX.push_back(PointsX[a]);
		m_vctClockY.push_back(PointsZ[a]);
	}
	flag = ClockWisePro(Num, m_vctClockX, m_vctClockY, m_vctClockRX, m_vctClockRY);//���нǵ����򣬽��ǵ㰴˳ʱ��˳�����
	for (int b=0;b<Num;b++)//������ת����͹����λ��ֺ��������ʽ
	{
		PointsX[b] = m_vctClockRX[b];
		PointsZ[b] = m_vctClockRY[b];
	}
	m_vctClockRX.clear();
	m_vctClockRY.clear();
	m_vctClockX.clear();
	m_vctClockY.clear();
	return flag;
}

int CMy3DConvertDlg::ClockWiseYZ(int Num,double *PointsY,double *PointsZ)
{
	int flag;
	for (int a=0;a<Num;a++)//������ת��Ϊ�����������ʽ
	{
		m_vctClockX.push_back(PointsY[a]);
		m_vctClockY.push_back(PointsZ[a]);
	}
	flag = ClockWisePro(Num, m_vctClockX, m_vctClockY, m_vctClockRX, m_vctClockRY);//���нǵ����򣬽��ǵ㰴˳ʱ��˳�����
	for (int b=0;b<Num;b++)//������ת����͹����λ��ֺ��������ʽ
	{
		PointsY[b] = m_vctClockRX[b];
		PointsZ[b] = m_vctClockRY[b];
	}
	m_vctClockRX.clear();
	m_vctClockRY.clear();
	m_vctClockX.clear();
	m_vctClockY.clear();
	return flag;
}
//�����ʷֺ���������������Ƭ��д�ļ�
//Num �����
//PointsX ����x��������
//PointsY ����y��������
//PointsZ ����z��������
//transferTo3D ��������ļ���
void CMy3DConvertDlg::creatTriangulation(int Num,double *PointsX,double *PointsY,double *PointsZ,FILE * transferTo3D,int polygonNumber,int Building)
{
	vector<qxwPolygon> vecPoly;
	vector<vector<int>> vctNum;
	double MaxMinX[2],MaxMinY[2],MaxMinZ[2],MinMaxX,MinMaxY,MinMaxZ;
	//////////////
	double angle1,angle2;
	int flag1 = 0,trigger = 0;

	//���߷�����ά�潵Ϊƽ��
	int	indexNvector = -1;
	for(int i = 0;i<Num-2;i++)
	{
		double SideA,SideB,SideC,pABC,sABC;
		SideA = sqrt(pow((PointsX[i]-PointsX[i+1]),2)+pow((PointsY[i]-PointsY[i+1]),2)+pow((PointsZ[i]-PointsZ[i+1]),2));
		SideB = sqrt(pow((PointsX[i]-PointsX[i+2]),2)+pow((PointsY[i]-PointsY[i+2]),2)+pow((PointsZ[i]-PointsZ[i+2]),2));
		SideC = sqrt(pow((PointsX[i+2]-PointsX[i+1]),2)+pow((PointsY[i+2]-PointsY[i+1]),2)+pow((PointsZ[i+2]-PointsZ[i+1]),2));
		pABC = (SideA+SideB+SideC)/2;
		sABC = sqrt(pABC*(pABC-SideA)*(pABC-SideB)*(pABC-SideC));	//���׹�ʽ������ж������Ƿ���
		if(sABC <= 10e-2)												//����ж���ֵ �ɱ�
			continue;
		else														//�������������Ϊ0���������㷨����
		{
			indexNvector = i;
			break;
		}
		/*
		double cosA;
		vectorA.x = PointsX[i]-PointsX[i+1];	//�˴�Ҳ���������������ж������Ƿ���
		vectorA.y = PointsY[i]-PointsY[i+1];
		vectorA.z = PointsZ[i]-PointsZ[i+1];
		vectorB.x = PointsX[i]-PointsX[i+2];
		vectorB.y = PointsY[i]-PointsY[i+2];
		vectorB.z = PointsZ[i]-PointsZ[i+2];
		cosA = (vectorA.x * vectorB.x + vectorA.y * vectorB.y + vectorA.z * vectorB.z) / (sqrt(pow(vectorA.x,2)+pow(vectorA.y,2)+pow(vectorA.z,2))+sqrt(pow(vectorB.x,2)+pow(vectorB.y,2)+pow(vectorB.z,2)));
		if(cosA > 10e-6)												
			continue;
		else														//������������ֵΪ0���������㷨����
		{
			indexNvector = i;
			break;
		}*/

	}
	int flag = 0;
	int flagE = 0;
	if(indexNvector == -1)
	{
		CString temp;
		temp.Format(_T("��%d����ķ�����������,�뼰ʱ����"),polygonNumber);
		//AfxMessageBox(temp);
		fprintf(corrigendumFILE,"	��%d����(Building#%d)�ķ�����������\n",polygonNumber,Building);
	}
	else
	{
		vectorA.x = PointsX[indexNvector]-PointsX[indexNvector+1];
		vectorA.y = PointsY[indexNvector]-PointsY[indexNvector+1];
		vectorA.z = PointsZ[indexNvector]-PointsZ[indexNvector+1];
		vectorB.x = PointsX[indexNvector]-PointsX[indexNvector+2];
		vectorB.y = PointsY[indexNvector]-PointsY[indexNvector+2];
		vectorB.z = PointsZ[indexNvector]-PointsZ[indexNvector+2];
		normVector.x = vectorA.y*vectorB.z-vectorA.z*vectorB.y;
		normVector.y = vectorB.x*vectorA.z-vectorA.x*vectorB.z;
		normVector.z = vectorA.x*vectorB.y-vectorA.y*vectorB.x;			//������
		/*
		findMinAndMax(Num,PointsX,MaxMinX);
		findMinAndMax(Num,PointsY,MaxMinY);
		findMinAndMax(Num,PointsZ,MaxMinZ);
		MinMaxX = MaxMinX[1]-MaxMinX[0];
		MinMaxY = MaxMinY[1]-MaxMinY[0];
		MinMaxZ = MaxMinZ[1]-MaxMinZ[0];*/

		if(abs(normVector.x) <= 10e-6)
		{
			if((abs(normVector.y) <= 10e-6)||(abs(normVector.y) < abs(normVector.z)))
			{
				flag = ClockWiseXY(Num,PointsX,PointsY);
				trigger = 1;
				ConcaveToBulge_iteration(Num,PointsX,PointsY,vecPoly);
				TransferToDotNum(Num,PointsX,PointsY,vecPoly,vctNum);
			}
			if((abs(normVector.z) <= 10e-6)||(abs(normVector.y) > abs(normVector.z)))
			{
				flag = ClockWiseXZ(Num,PointsX,PointsZ);
				trigger = 2;
				ConcaveToBulge_iteration(Num,PointsX,PointsZ,vecPoly);
				TransferToDotNum(Num,PointsX,PointsZ,vecPoly,vctNum);
			}
		}
		else if(abs(normVector.y <= 10e-6))
		{
			if((abs(normVector.z)<= 10e-6)||(abs(normVector.z) < abs(normVector.x)))
			{
				flag = ClockWiseYZ(Num,PointsY,PointsZ);
				trigger = 3;
				ConcaveToBulge_iteration(Num,PointsY,PointsZ,vecPoly);
				TransferToDotNum(Num,PointsY,PointsZ,vecPoly,vctNum);
			}
			if(abs(normVector.z) > abs(normVector.x))
			{
				flag = ClockWiseXY(Num,PointsX,PointsY);
				trigger = 1;
				ConcaveToBulge_iteration(Num,PointsX,PointsY,vecPoly);
				TransferToDotNum(Num,PointsX,PointsY,vecPoly,vctNum);
			}
		}
		else if(abs(normVector.z) <= 10e-6)
		{
			if(abs(normVector.y) > abs(normVector.x))
			{
				flag = ClockWiseXZ(Num,PointsX,PointsZ);
				trigger = 2;
				ConcaveToBulge_iteration(Num,PointsX,PointsZ,vecPoly);
				TransferToDotNum(Num,PointsX,PointsZ,vecPoly,vctNum);
			}
			if(abs(normVector.y) < abs(normVector.x))
			{
				flag = ClockWiseYZ(Num,PointsY,PointsZ);
				trigger = 3;
				ConcaveToBulge_iteration(Num,PointsY,PointsZ,vecPoly);
				TransferToDotNum(Num,PointsY,PointsZ,vecPoly,vctNum);
			}
		}
		else
		{
			if((abs(normVector.x) <= abs(normVector.z))&&(abs(normVector.y) <= abs(normVector.z)))
			{
				flag = ClockWiseXY(Num,PointsX,PointsY);
				trigger = 1;
				ConcaveToBulge_iteration(Num,PointsX,PointsY,vecPoly);
				TransferToDotNum(Num,PointsX,PointsY,vecPoly,vctNum);
			}
			if((abs(normVector.x) <= abs(normVector.y))&&(abs(normVector.z) <= abs(normVector.y)))
			{
				flag = ClockWiseXZ(Num,PointsX,PointsZ);
				trigger = 2;
				ConcaveToBulge_iteration(Num,PointsX,PointsZ,vecPoly);
				TransferToDotNum(Num,PointsX,PointsZ,vecPoly,vctNum);
			}
			if((abs(normVector.y) <= abs(normVector.x))&&(abs(normVector.z) <= abs(normVector.x)))
			{
				flag = ClockWiseYZ(Num,PointsY,PointsZ);
				trigger = 3;
				ConcaveToBulge_iteration(Num,PointsY,PointsZ,vecPoly);
				TransferToDotNum(Num,PointsY,PointsZ,vecPoly,vctNum);
			}
		}
	}
	if((vecPoly.size() == 0)||(vctNum.size() == 0))
	{
		flag = ClockWiseXY(Num,PointsX,PointsY);
		trigger = 1;
		ConcaveToBulge_iteration(Num,PointsX,PointsY,vecPoly);
		TransferToDotNum(Num,PointsX,PointsY,vecPoly,vctNum);
		flagE += vctNum.size();

		flag = ClockWiseXZ(Num,PointsX,PointsZ);
		trigger = 2;
		ConcaveToBulge_iteration(Num,PointsX,PointsZ,vecPoly);
		TransferToDotNum(Num,PointsX,PointsZ,vecPoly,vctNum);
		flagE += vctNum.size();

		flag = ClockWiseYZ(Num,PointsY,PointsZ);
		trigger = 3;
		ConcaveToBulge_iteration(Num,PointsY,PointsZ,vecPoly);
		TransferToDotNum(Num,PointsY,PointsZ,vecPoly,vctNum);
		flagE += vctNum.size();
	}
	/*for (int i=0;i<Num-2;i++)
	{
		angle1 = atan2(PointsX[i+1]-PointsX[i],PointsY[i+1]-PointsY[i]);
		angle2 = atan2(PointsX[i+2]-PointsX[i+1],PointsY[i+2]-PointsY[i+1]);
		if (fabs(angle1)>10e-6 && fabs(angle2)>10e-6)
		{
			if((fabs(angle1-angle2)>10e-6) && (fabs(angle1-angle2+PI)>10e-6) && (fabs(angle1-angle2-PI)>10e-6))
			{
				flag1 = 1;//���ڲ����ߵĵ�
				break;
			}
		}

	}
	//int flag = 0;
	if (flag1==0)//XY����ȫ����
	{
		findMinAndMax(Num,PointsX,MaxMinX);
		findMinAndMax(Num,PointsY,MaxMinY);
		MinMaxX = MaxMinX[1]-MaxMinX[0];
		MinMaxY = MaxMinY[1]-MaxMinY[0];
		if(MinMaxX<=MinMaxY)
		{
			for (int a=0;a<Num;a++)//������ת��Ϊ�����������ʽ
			{
				m_vctClockX.push_back(PointsY[a]);
				m_vctClockY.push_back(PointsZ[a]);
			}
			flag = ClockWisePro(Num, m_vctClockX, m_vctClockY, m_vctClockRX, m_vctClockRY);//���нǵ����򣬽��ǵ㰴˳ʱ��˳�����
			trigger = 3;
			for (int b=0;b<Num;b++)//������ת����͹����λ��ֺ��������ʽ
			{
				PointsY[b] = m_vctClockRX[b];
				PointsZ[b] = m_vctClockRY[b];
			}
			m_vctClockRX.clear();
			m_vctClockRY.clear();
			m_vctClockX.clear();
			m_vctClockY.clear();

			ConcaveToBulge_iteration(Num,PointsY,PointsZ,vecPoly);
			TransferToDotNum(Num,PointsY,PointsZ,vecPoly,vctNum);

		}
		else
		{
			for (int a=0;a<Num;a++)//������ת��Ϊ�����������ʽ
			{
				m_vctClockX.push_back(PointsX[a]);
				m_vctClockY.push_back(PointsZ[a]);
			}
			flag = ClockWisePro(Num, m_vctClockX, m_vctClockY, m_vctClockRX, m_vctClockRY);//���нǵ����򣬽��ǵ㰴˳ʱ��˳�����
			trigger = 2;
			for (int b=0;b<Num;b++)//������ת����͹����λ��ֺ��������ʽ
			{
				PointsX[b] = m_vctClockRX[b];
				PointsZ[b] = m_vctClockRY[b];
			}
			m_vctClockRX.clear();
			m_vctClockRY.clear();
			m_vctClockX.clear();
			m_vctClockY.clear();
			ConcaveToBulge_iteration(Num,PointsX,PointsZ,vecPoly);
			TransferToDotNum(Num,PointsX,PointsZ,vecPoly,vctNum);
		}
	}
	else//XY���겻����
	{
		for (int a=0;a<Num;a++)//������ת��Ϊ�����������ʽ
			{
			m_vctClockX.push_back(PointsX[a]);
			m_vctClockY.push_back(PointsY[a]);
			}
		flag = ClockWisePro(Num, m_vctClockX, m_vctClockY, m_vctClockRX, m_vctClockRY);//���нǵ����򣬽��ǵ㰴˳ʱ��˳�����
		trigger = 1;
		for (int b=0;b<Num;b++)//������ת����͹����λ��ֺ��������ʽ
			{
			PointsX[b] = m_vctClockRX[b];
			PointsY[b] = m_vctClockRY[b];
			}
		m_vctClockRX.clear();
		m_vctClockRY.clear();
		m_vctClockX.clear();
		m_vctClockY.clear();
		ConcaveToBulge_iteration(Num,PointsX,PointsY,vecPoly);
		TransferToDotNum(Num,PointsX,PointsY,vecPoly,vctNum);
	}
	/////////////
	*/	
	if (flag == 1)		//���ת��˳�򣬽����˳��ת����
	{
		for (int a=0;a<vctNum.size();a++)
		{
			for (int b=0;b<vctNum[a].size();b++)
			{
				int x = vctNum[a][b];
				vctNum[a][b] = Num-x-1;
			}
		}
		
		if(trigger == 1)
		{
			double temp = 0.0;
			for(int i=0;i<(Num/2);i++)
			{
				temp = PointsX[i];
				PointsX[i] = PointsX[Num-1 - i];
				PointsX[Num-1 - i] = temp;
				temp = PointsY[i];
				PointsY[i] = PointsY[Num-1 - i];
				PointsY[Num-1 - i] = temp;
			}
		}
		if(trigger == 2)
		{
			double temp = 0.0;
			for(int i=0;i<(Num/2);i++)
			{
				temp = PointsX[i];
				PointsX[i] = PointsX[Num-1 - i];
				PointsX[Num-1 - i] = temp;
				temp = PointsZ[i];
				PointsZ[i] = PointsZ[Num-1 - i];
				PointsZ[Num-1 - i] = temp;
			}
		}
		if(trigger == 3)
		{
			double temp = 0.0;
			for(int i=0;i<(Num/2);i++)
			{
				temp = PointsY[i];
				PointsY[i] = PointsY[Num-1 - i];
				PointsY[Num-1 - i] = temp;
				temp = PointsZ[i];
				PointsZ[i] = PointsZ[Num-1 - i];
				PointsZ[Num-1 - i] = temp;
			}
		}
		
	}
	vector<pointsArray> vctTriangle;
	for (int i=0;i<vctNum.size();i++)
	{
		Triangulation(vctNum[i],vctTriangle);
	}
	
	if(vctTriangle.size() == 0)
	{
		CString temp;
		if(flagE == 0)
		{
			temp.Format(_T("��%d�����������ʷ��������\n���á�	3	0	1	2�����棬�뼰ʱ����"),polygonNumber);
			//AfxMessageBox(temp);
			fprintf(corrigendumFILE,"	��%d������(Building#%d)�����ʷ��������,���ܿ��ܴ����߽���\n",polygonNumber,Building);
		}
		else
		{
			temp.Format(_T("��%d�����������ʷ��������\n���á�	3	0	1	2�����棬�뼰ʱ����"),polygonNumber);
			//AfxMessageBox(temp);
			fprintf(corrigendumFILE,"	��%d������(Building#%d)�����ʷ��������\n",polygonNumber,Building);
		}
		
		for(int i = 0;i<Num-2;i++)
		{
			fprintf(transferTo3D,"           3           0           %d           %d\n",i+1,i+2);
		}
	}
	else
		for (int i=0;i<vctTriangle.size();i++)
		{
			fprintf(transferTo3D,"           %d           %d           %d           %d\n",vctTriangle[i].a,vctTriangle[i].b,vctTriangle[i].c,vctTriangle[i].d);
		}
	
}

//���͹����λ��ֺ���ݵ�����ƥ������
//Num �����
//pX ԭʼ������ǰ��x��������
//pY ԭʼ������ǰ��y��������
//vecPoly ���ֺ����������
//vctNum �洢����ŵ���������νṹΪ<�����<�����>>�����ջ��ֺ�ĵ�˳������ԭʼ�����
bool CMy3DConvertDlg::TransferToDotNum(int Num,double *pX,double *pY,vector<qxwPolygon> vecPoly,vector<vector<int>> &vctNum)
{
	int iTmp = 0;
	for(int a=0;a<vecPoly.size();a++)//ѭ��͹����θ���
		{
		qxwPolygon vctDot = vecPoly[a];
		vector<int> vctiTmp;
		for(int b=0;b<vctDot.size();b++)//ѭ��һ��͹����εĵ�
			{
				vector<int> vctTmp;
				for(int c=0;c<Num;c++)//���õ�X������pX�ȶԣ���������õ�X������ͬ��pX����ż�¼��vctTmp��
					{
					if(fabs(pX[c]-vctDot[b].x)<10e-6)
						vctTmp.push_back(c);
					}
				if(0 == vctTmp.size())
					return 0;
				for(int d=0;d<vctTmp.size();d++)//��vctTmp��ѭ��pY���꣬�ҵ���õ�Y������ͬ�ĵ���ţ���¼��vctiTmp��
					{
					if(fabs(pY[vctTmp[d]]-vctDot[b].y)<10e-6)
						vctiTmp.push_back(vctTmp[d]);
					}
			}
		vctNum.push_back(vctiTmp);
		}
	return 1;
}
//���������ʷֹ�ϵ
//vctNum ���������
//vctTriangle ��������ʷֵ�����
bool CMy3DConvertDlg::Triangulation(vector<int> vctNum,vector<pointsArray> &vctTriangle)
{
	Node nDot;
	int iSize = vctNum.size();
	int iTmp = 0;
	for(int a=0;a<iSize-2;a++)
		{
		nDot.a = 3;
		nDot.b = vctNum[0];
		nDot.c = vctNum[iTmp+1];
		nDot.d = vctNum[iTmp+2];
		vctTriangle.push_back(nDot);
		iTmp++;
		}
	return 1;
}

void CMy3DConvertDlg::coplanerJudge(int pointsNum,int polygonNum,double *PointsX, double *PointsY,double *PointsZ)
{
	int i;
	for(i=0;i<pointsNum-3;i++)					//���߷��й���
	{
		vectorA.x = PointsX[i]-PointsX[i+1];
		vectorA.y = PointsY[i]-PointsY[i+1];
		vectorA.z = PointsZ[i]-PointsZ[i+1];
		vectorB.x = PointsX[i]-PointsX[i+2];
		vectorB.y = PointsY[i]-PointsY[i+2];
		vectorB.z = PointsZ[i]-PointsZ[i+2];
		normVector.x = vectorA.y*vectorB.z-vectorA.z*vectorB.y;
		normVector.y = vectorB.x*vectorA.z-vectorA.x*vectorB.z;
		normVector.z = vectorA.x*vectorB.y-vectorA.y*vectorB.x;
		vectorA.x = PointsX[i+1]-PointsX[i+2];
		vectorA.y = PointsY[i+1]-PointsY[i+2];
		vectorA.z = PointsZ[i+1]-PointsZ[i+2];
		vectorB.x = PointsX[i+1]-PointsX[i+3];
		vectorB.y = PointsY[i+1]-PointsY[i+3];
		vectorB.z = PointsZ[i+1]-PointsZ[i+3];
		normVectorTemp.x = vectorA.y*vectorB.z-vectorA.z*vectorB.y;
		normVectorTemp.y = vectorB.x*vectorA.z-vectorA.x*vectorB.z;
		normVectorTemp.z = vectorA.x*vectorB.y-vectorA.y*vectorB.x;
		/*if(((abs(normVector.x)/abs(normVectorTemp.x)-abs(normVector.y)/abs(normVectorTemp.y))<=10e-2)&&((abs(normVector.x)/abs(normVectorTemp.x)-abs(normVector.z)/abs(normVectorTemp.z))<=10e-2)&&((abs(normVector.y)/abs(normVectorTemp.y)-abs(normVector.z)/abs(normVectorTemp.z))<=10e-2))
		{
			CString temp;
			temp.Format(_T("��%d������ڲ�����ĵ�,��ע�Ⲣ����"),polygonNum);
		    AfxMessageBox(temp);
			continue;
		}*/
	}
	for(i = 0;i<pointsNum-2;i++)
	{
		double SideA,SideB,SideC,pABC,sABC;
		SideA = sqrt(pow((PointsX[i]-PointsX[i+1]),2)+pow((PointsY[i]-PointsY[i+1]),2)+pow((PointsZ[i]-PointsZ[i+1]),2));
		SideB = sqrt(pow((PointsX[i]-PointsX[i+2]),2)+pow((PointsY[i]-PointsY[i+2]),2)+pow((PointsZ[i]-PointsZ[i+2]),2));
		SideC = sqrt(pow((PointsX[i+2]-PointsX[i+1]),2)+pow((PointsY[i+2]-PointsY[i+1]),2)+pow((PointsZ[i+2]-PointsZ[i+1]),2));
		pABC = (SideA+SideB+SideC)/2;
		sABC = sqrt(pABC*(pABC-SideA)*(pABC-SideB)*(pABC-SideC));	//���׹�ʽ������ж������Ƿ���
		if(sABC <= 10e-2)												//����жϷ�ֵ �ɱ�
			continue;
		else														//�������������Ϊ0���������㷨����
		{
			break;
		}
	}
	double a,b,c,d;
	a = (PointsY[i+1]-PointsY[i])*(PointsZ[i+2]-PointsZ[i])-(PointsZ[i+1]-PointsZ[i])*(PointsY[i+2]-PointsY[i]);
	b = (PointsZ[i+1]-PointsZ[i])*(PointsX[i+2]-PointsX[i])-(PointsX[i+1]-PointsX[i])*(PointsZ[i+2]-PointsZ[i]);
	c = (PointsX[i+1]-PointsX[i])*(PointsY[i+2]-PointsY[i])-(PointsY[i+1]-PointsY[i])*(PointsX[i+2]-PointsX[i]);
	d = 0-(a*PointsX[i]+b*PointsY[i]+c*PointsZ[i]);
	for(i = 0;i<pointsNum;i++)
	{
		double result;
		result = a*PointsX[i]+b*PointsY[i]+c*PointsZ[i]+d;
		if(abs(result)>=1)
		{
			CString temp;
			temp.Format(_T("��%d������ڲ�����ĵ�,��ע�Ⲣ����"),polygonNum);
		    AfxMessageBox(temp);
			continue;
		}
	}
	return ;
}
//��shp�ļ�д3d�ļ����ܺ���
//m_ShpFile_fp shp�ļ�ָ��
//transferTo3D 3d�ļ�ָ��
int CMy3DConvertDlg::OnReadPolygonShp(FILE * m_ShpFile_fp,DBFHandle hDBF,FILE * transferTo3D)
{

	int RecordNumber;
	int ContentLength;
	int lenthOfBuildingPointsArray = 0;
	double buildingPointsArray[100][3]={0.0,0.0,0.0};
	int polygonNumber = 0;
	int wrongCount = 0;
	while((fread(&RecordNumber,sizeof(int),1,m_ShpFile_fp)!=0))
	{
		fread(&ContentLength,sizeof(int),   1,m_ShpFile_fp);
		RecordNumber      = OnChangeByteOrder (RecordNumber);
		printf("\n\nRecordNumber=%d\n",RecordNumber);
		ContentLength     = OnChangeByteOrder (ContentLength);
		printf("ContentLength=%d\n",ContentLength);
		int shapeType;
		double Box[4];
		int NumParts;
		int NumPoints;
		int *Parts;
		int i,j;
		fread(&shapeType,    sizeof(int),   1,m_ShpFile_fp);

		if (15 != shapeType)
		{
			CString Hint;
			Hint.Format("shp�ļ��ļ�������Ϊ%d,�޷�ת���ļ�����Ҫ��������Ϊ15��PlygonZ��\n shp�ļ��������Ͳο���\n0 Null Shape\n1 Point\n3 PolyLine\n5 Polygon\n8 MultiPoint\n11 PointZ\n13 PolyLineZ\n15 PolygonZ\n18 MultiPointZ\n21 PointM\n23 PolyLineM\n25 PolygonM\n28 MultiPointM\n31 MultiPatch\n",shapeType); 
			AfxMessageBox(Hint);
		}
		//printf("shapeType=%d\n",shapeType);

		//��Box
		for(i=0;i<4;i++)
			{
			fread(Box+i,     sizeof(double),1,m_ShpFile_fp);
			printf("Box[%d]=%lf\n",i,*(Box+i));
			}
		//��NumParts��NumPoints
		fread(&NumParts,     sizeof(int),   1,m_ShpFile_fp);
		//printf("NumParts=%d\n",NumParts);
		fread(&NumPoints,    sizeof(int),   1,m_ShpFile_fp);
		//printf("NumPoints=%d\n",NumPoints);

		//��Parts��Points
		Parts       =new int[NumParts];
		for(i=0;i<NumParts;i++)
			{
			fread(Parts+i,   sizeof(int),   1,m_ShpFile_fp);
			//printf("Parts+%d=%d\n",i,*(Parts+i));
			}
		int actualNumPoints;
		actualNumPoints = NumPoints;
		int R,G,B;
		double Triangulation;
		R=220;
		G=220;
		B=220;
	
		double *PointsX;
		double *PointsY;    
		PointsX =new double[NumPoints];
		PointsY =new double[NumPoints];

		for(j=0;j<NumPoints;j++)
			{
			fread(PointsX+j, sizeof(double),1,m_ShpFile_fp);
			fread(PointsY+j, sizeof(double),1,m_ShpFile_fp);
			
			//printf("X[%d]=%lf\n",j,*(PointsX+j));
			//printf("Y[%d]=%lf\n",j,*(PointsY+j));
			
			}
		double Zrange[2];
		for(j=0;j<2;j++)
			{
			fread(Zrange+j, sizeof(double),1,m_ShpFile_fp);
			//printf("Zrange[%d]=%lf\n",j,*(Zrange+j));
			}
		double *PointsZ;
		PointsZ =new double[NumPoints];
		double *actualPointsX;
		double *actualPointsY;
		double *actualPointsZ;
		actualPointsX =new double[NumPoints];
		actualPointsY =new double[NumPoints];
		actualPointsZ =new double[NumPoints];
		int counter = 0;
		for(j=0;j<NumPoints;j++)
		{
			fread(PointsZ+j, sizeof(double),1,m_ShpFile_fp);
			if(j==0)
			{
				actualPointsX[counter] = PointsX[j];
				actualPointsY[counter] = PointsY[j];
				actualPointsZ[counter] = PointsZ[j];
				counter++;
			}
			int flagZ = 0;
			for(int k=0;k<j;k++)
			{
				if(PointsX[k]==PointsX[j]&&PointsY[k]==PointsY[j]&&PointsZ[k]==PointsZ[j])
				{
					actualNumPoints--;
					flagZ = 1;
					break;
				}
			}
			if(flagZ == 0 && j != 0)
			{
				actualPointsX[counter] = PointsX[j];
				actualPointsY[counter] = PointsY[j];
				actualPointsZ[counter] = PointsZ[j];
				counter++;
			}
		}
		double Mrange[2];
		for(j=0;j<2;j++)
			{
			fread(Mrange+j, sizeof(double),1,m_ShpFile_fp);
			//printf("Mrange[%d]=%lf\n",j,*(Mrange+j));
			}
		double *PointsM;  
		PointsM =new double[NumPoints];        
		for(j=0;j<NumPoints;j++)
			{
			fread(PointsM+j, sizeof(double),1,m_ShpFile_fp);
			//printf("M[%d]=%lf\n",j,*(PointsM+j));
			}
		delete[] PointsM;
		delete[] PointsX;
		delete[] PointsY;
		delete[] PointsZ;
		Triangulation = (actualNumPoints-2)*4;
		
		int Building = 0;
		if(actualNumPoints <= 2)	
		{
			CString temp;
			wrongCount++;
			Building = DBFReadIntegerAttribute(hDBF,polygonNumber-1,0);
			temp.Format(_T("��%d�����ݵ����Ϊ%d,��ע�Ⲣ����"),RecordNumber,actualNumPoints);
		    //AfxMessageBox(temp);
			fprintf(corrigendumFILE,"	��%d������(Building#%d)����Ϊ%d\n",RecordNumber,Building,actualNumPoints);
			continue;
		}

		/*if(actualNumPoints > 3)						//�����ж�
		{
			coplanerJudge(actualNumPoints,polygonNumber+wrongCount,actualPointsX,actualPointsY,actualPointsZ);
		}*/
		
		polygonNumber++;

		Building = DBFReadIntegerAttribute(hDBF,polygonNumber-1+wrongCount,0);
//��������Լ��
		
		/*
		int numbOfSamePoints = 0;
		if(buildingPointsArray[0][0]==0.0&&buildingPointsArray[0][1]==0.0&&buildingPointsArray[0][2]==0.0) //�������ʼֵ
		{
			for(int i=0;i<actualNumPoints;i++)
			{
				buildingPointsArray[i][0]=actualPointsX[i];
				buildingPointsArray[i][1]=actualPointsY[i];
				buildingPointsArray[i][2]=actualPointsZ[i];
				lenthOfBuildingPointsArray++;
			}
		}
		else


		{
			for(int j=0;j<actualNumPoints;j++)
			{
				int flagJK = 0;
				for(int k=0;k<lenthOfBuildingPointsArray;k++)
				{
					if(actualPointsX[j]==buildingPointsArray[k][0]&&actualPointsY[j]==buildingPointsArray[k][1]&&actualPointsZ[j]==buildingPointsArray[k][2])	
					{
						numbOfSamePoints++;
						break;
					}
					else flagJK++;
				}
				if(flagJK == lenthOfBuildingPointsArray)
				{
					buildingPointsArray[lenthOfBuildingPointsArray][0] = actualPointsX[j];
					buildingPointsArray[lenthOfBuildingPointsArray][1] = actualPointsY[j];
					buildingPointsArray[lenthOfBuildingPointsArray][2] = actualPointsZ[j];
					lenthOfBuildingPointsArray++;
				}
			}	
			if(numbOfSamePoints < 2)
			{
				BuildingNumb++;
				lenthOfBuildingPointsArray = 0;
				for(int i=0;i<actualNumPoints;i++)
				{
					buildingPointsArray[i][0]=actualPointsX[i];
					buildingPointsArray[i][1]=actualPointsY[i];
					buildingPointsArray[i][2]=actualPointsZ[i];
					lenthOfBuildingPointsArray++;
				}
			}
		}
		*/
		
		fprintf(transferTo3D,"          %d             %d     %d   %d   %d     0     0                1.000000000            %d             0       0     1   255     2     0     1     0     1     0                  1.000000000                0.000000000    BUILDING %d#\\/#\n",actualNumPoints,polygonNumber,R,G,B,int(Triangulation),Building);

		if(actualNumPoints>2) creatTriangulation(actualNumPoints,actualPointsX,actualPointsY,actualPointsZ,transferTo3D,polygonNumber,Building);

		for(j=0;j<actualNumPoints;j++)
			{
				if(actualPointsZ[j]<100)
					//fprintf(transferTo3D,"            %s%.11lf               %.10lf                    %.15lf                    -1.000000000                    -1.000000000\n",m_DH,actualPointsX[j],actualPointsY[j],actualPointsZ[j]);
					fprintf(transferTo3D,"            %s%.6lf               %.6lf                    %.6lf                    -1.000000                    -1.000000\n",m_DH,actualPointsX[j],actualPointsY[j],actualPointsZ[j]);
				else
					//fprintf(transferTo3D,"            %s%.11lf               %.10lf                    %.15lf                    -1.000000000                    -1.000000000\n",m_DH,actualPointsX[j],actualPointsY[j],actualPointsZ[j]);
					fprintf(transferTo3D,"            %s%.6lf               %.6lf                   %.6lf                    -1.000000                    -1.000000\n",m_DH,actualPointsX[j],actualPointsY[j],actualPointsZ[j]);
			}
		delete[] actualPointsX;
		delete[] actualPointsY;
		delete[] actualPointsZ;

		delete[] Parts;
		}
	return polygonNumber;
}
//��shp�ļ�д3d�ļ���������д�ļ�ͷ
//m_InFILE shp�ļ�ָ��
//DBFHandle hDBF dbf�ļ����
//m_OutFILE 3d�ļ�ָ��
void CMy3DConvertDlg::transferSHPto3D(FILE * m_InFILE,DBFHandle hDBF,FILE * m_OutFILE)
{
	int fileCode =-1;
    int fileLength=-1;
	int version=-1;
	int shapeType=-1;

	int i;
	int FileCode;
	int Unused;
	int FileLength;
	int Version;
	int ShapeType;
	double Xmin;
	double Ymin;
	double Xmax;
	double Ymax;
	double Zmin;
	double Zmax;
	double Mmin;
	double Mmax;

	fread(&FileCode,sizeof(int),1,m_InFILE);
	FileCode = OnChangeByteOrder(FileCode);
	printf("FileCode=%d\n",FileCode);

	for(i=0;i<5;i++)
		fread(&Unused,sizeof(int),1,m_InFILE);
	fread(&FileLength,sizeof(int),1,m_InFILE);
	FileLength=OnChangeByteOrder(FileLength);
	printf("FileLength=%d\n",FileLength);

	fread(&Version,sizeof(int),1,m_InFILE);
	printf("Version=%d\n",Version);

	fread(&ShapeType,sizeof(int),1,m_InFILE);
	printf("ShapeType=%d\n",ShapeType);

	fread(&Xmin,         sizeof(double),1,m_InFILE);
	printf("Xmin=%lf\n",Xmin);

	fread(&Ymin,         sizeof(double),1,m_InFILE);
	printf("Ymin=%lf\n",Ymin);

	fread(&Xmax,         sizeof(double),1,m_InFILE);
	printf("Xmax=%lf\n",Xmax);

	fread(&Ymax,         sizeof(double),1,m_InFILE);
	printf("Ymax=%lf\n",Ymax);

	fread(&Zmin,         sizeof(double),1,m_InFILE);
	printf("Zmin=%lf\n",Zmin);

	fread(&Zmax,        sizeof(double),1,m_InFILE);
	printf("Zmax=%lf\n",Zmax);

	fread(&Mmin,         sizeof(double),1,m_InFILE);
	printf("Mmin=%lf\n",Mmin);

	fread(&Mmax,         sizeof(double),1,m_InFILE);
	printf("Mmax=%lf\n",Mmax);

	// ��ȡ�����ļ�ͷ�����ݽ���
//----------------------------------------------------------------
	//�߶�ȡ�������ݱ�д��.3d�ļ�
	fprintf(m_OutFILE,"IMAGIS3D2.0\n");
	fprintf(m_OutFILE,"\t\t0\n");
	fprintf(m_OutFILE,"Meter\n");
	fprintf(m_OutFILE,"#\\/#\n");
	fprintf(m_OutFILE,"\t\t%s%.9lf             %s%.9lf            %.9f            %.9f                %.9f                %.9f\n",m_DH,Xmin,m_DH,Xmax,Ymin,Ymax,Zmin,Zmax);
	fprintf(m_OutFILE,"\t\t1\n");
	fprintf(m_OutFILE,"Free Layer\n");
	fprintf(m_OutFILE,"\t\t0\n");
	fprintf(m_OutFILE,"\t\t0\n");
	fprintf(m_OutFILE,"\t\t0\n");
	fprintf(m_OutFILE,"\t\t0\n");
	fprintf(m_OutFILE,"\t\t0\n");
	fprintf(m_OutFILE,"\t\t0\n");
	int totalPolygonNum;
	totalPolygonNum = OnReadPolygonShp(m_InFILE,hDBF,m_OutFILE);

	fseek(m_OutFILE,0,SEEK_SET);

	/*char s[300];
	int cout=0;
	int len;
	for(i=0;i<9;i++)
	{
		fgets(s,300,transferTo3D);
		len = strlen(s);
		cout+=len;
	}
	fseek(transferTo3D,cout-1,SEEK_SET);
	//fprintf(transferTo3D,"pain6dlsty");
	fprintf(transferTo3D,"         %d\n",totalPolygonNum);*/
	fprintf(m_OutFILE,"IMAGIS3D2.0\n");
	fprintf(m_OutFILE,"\t\t0\n");
	fprintf(m_OutFILE,"Meter\n");
	fprintf(m_OutFILE,"#\\/#\n");
	fprintf(m_OutFILE,"\t\t%s%.9lf             %s%.9lf            %.9f            %.9f                %.9f                %.9f\n",m_DH,Xmin,m_DH,Xmax,Ymin,Ymax,Zmin,Zmax);
	fprintf(m_OutFILE,"\t\t1\n");
	fprintf(m_OutFILE,"Free Layer\n");
	fprintf(m_OutFILE,"\t\t0\n");
	fprintf(m_OutFILE,"\t\t0\n");
	fprintf(m_OutFILE,"\t\t%d\n",totalPolygonNum);
	fprintf(m_OutFILE,"\t\t0\n");
	fprintf(m_OutFILE,"\t\t0\n");
	fprintf(m_OutFILE,"\t\t0\n");
	fclose(m_OutFILE);
	ReadBuildingsFrom3d(Out_PathX,Info3DFile);
	Out_PathX = Out_Path+"\\3d�ļ���������.txt";
	BuildAnilysize_CulFaceOrigin(Info3DFile,Out_PathX);
	DBFClose(hDBF);

	return ;
}

//��3d�ļ�ת��Ϊshp�ļ�������
//infile ����3D�ļ�ָ��
//hSHP dbf�ļ����
//hDBF 3d�ļ�ָ��
void CMy3DConvertDlg::transfer3DtoSHP(fstream &infile,SHPHandle hSHP,DBFHandle hDBF)
{
	double Xmin;
	double Ymin;
	double Xmax;
	double Ymax;
	double Zmin;
	double Zmax;
// 	fstream infile("D:\\2018.4.18\\projects\\test5.15\\test_0319.3d");
// 	SHPHandle hSHP = SHPCreate("D:\\2018.4.18\\projects\\test5.15\\0319SHP@1.shp",SHPT_POLYGONZ);
// 	DBFHandle hDBF = DBFCreate("D:\\2018.4.18\\projects\\test5.15\\0319SHP@1.dbf");
	if(!infile)
	{
		cout<<".3d file can't open"<<endl;
	}
	else
	{
		
		string a;

/////20180601 zqs add, to recognize tab//////////////
		for (int i=0;i<4;i++)
		{	
			infile >> a;
		}
		infile >> Xmin;
		infile >> Xmax;
		infile >> Ymin;
		infile >> Ymax;
		infile >> Zmin;
		infile >> Zmax;
/////20180601 zqs add end//////////////

		//for(int i=0;i<4;i++) 
		//{
		//	getline(infile, a, '\n');
		//	cout<<a<<endl;
		//}
		//std::cout.precision(16);
		//getline(infile, a, '\n');
		//int cursor=0;
		//char b[50];
		//int index=0;
		//while(a[cursor]==' ')  cursor++;
		//while(a[cursor]!=' ')  b[index++]=a[cursor++];
		//b[index]='\0';
		//Xmin = stod(b);
		//index = 0;
		//while(a[cursor]==' ')  cursor++;
		//while(a[cursor]!=' ')  b[index++]=a[cursor++];
		//b[index]='\0';
		//Xmax = stod(b);
		//index = 0;
		//while(a[cursor]==' ')  cursor++;
		//while(a[cursor]!=' ')  b[index++]=a[cursor++];
		//b[index]='\0';
		//Ymin = stod(b);
		//index = 0;
		//while(a[cursor]==' ')  cursor++;
		//while(a[cursor]!=' ')  b[index++]=a[cursor++];
		//b[index]='\0';
		//Ymax = stod(b);
		//index = 0;
		//while(a[cursor]==' ')  cursor++;
		//while(a[cursor]!=' ')  b[index++]=a[cursor++];
		//b[index]='\0';
		//Zmin = stod(b);
		//index = 0;
		//while(a[cursor]==' ')  cursor++;
		//int lenth = 0;
		//lenth = a.size();
		//while(a[cursor]!=' ')  
		//{
		//	b[index++]=a[cursor++];
		//	if(cursor == lenth) break;
		//}
		//b[index]='\0';
		//Zmax = stod(b);

		//cout<<Xmin<<"\t"<<Xmax<<"\t"<<Ymin<<"\t"<<Ymax<<"\t"<<Zmin<<"\t"<<Zmax<<endl;

		for(int i=0;i<5;i++) 
		{
			getline(infile, a, '\n');
		}

		int totlePolygonNumb;
		//getline(infile, a, '\n');
		//totlePolygonNumb = stoi(a);		
		//cout<<totlePolygonNumb<<endl;
		infile >> totlePolygonNumb;
        for(int i=0;i<3;i++) infile >> a;
		getline(infile, a, '\n');
		DBFAddField(hDBF,"SS_ID",FTInteger,8,0);
		DBFAddField(hDBF,"STORIES",FTInteger,8,0);
		DBFAddField(hDBF,"shade_flag",FTString,14,0);
		DBFAddField(hDBF,"shade_valu",FTString,255,0);
		int count=0;
		int wrongCount = 0;
		while (count<totlePolygonNumb)  
		{  
		
			int pointsNumb;
			int polygonNumb;
			int buildingNumb;
			

			/*getline(infile, a, '\n');
			int cursor=0;
			char b[50];
			int index=0;
			while(a[cursor]==' ')  cursor++;
			while(a[cursor]!=' ')  b[index++]=a[cursor++];
			b[index]='\0';
			str2int(pointsNumb,b);
			index = 0;
			while(a[cursor]==' ')  cursor++;
			while(a[cursor]!=' ')  b[index++]=a[cursor++];
			b[index]='\0';
			str2int(polygonNumb,b);
			index = 0;
			while(a[cursor]!='G') cursor++;
			cursor++;
			while(a[cursor]!='#')  b[index++]=a[cursor++];
			b[index]='\0';
			str2int(buildingNumb,b);
			cout<<pointsNumb<<"\t"<<polygonNumb<<"\t"<<buildingNumb<<endl;*/
			infile >> pointsNumb;
			infile >> polygonNumb;
			for(int i=0;i<20;i++)	infile >> a;
			infile >> buildingNumb;
			getline(infile, a, '\n');
			int lenth = 0;
			//�жϲ���ȥ�����С�ڵ���2�Ĳ����������
			
			/*if(pointsNumb <= 2)
			{
				wrongCount++;
				int wrongCounter = 0;
				while(wrongCounter<pointsNumb)
				{
					getline(infile, a, '\n');
					lenth = a.length();
					int flag = 0;
					for(int i = 0;i<lenth;i++)
					{
						if(a[i] == '.') 
						{
							flag = 1;
							break;
						}
					}
					if(flag == 1)	wrongCounter++;
				}
				count++;
				continue;
			}*/

			if(pointsNumb <= 2) 
			{
				CString temp;
				temp.Format(_T("��%d�����ݵ����Ϊ%d,��ע�Ⲣ����"),polygonNumb,pointsNumb);
				//AfxMessageBox(temp);
				fprintf(corrigendumFILE,"	��%d������(Building#%d)�����Ϊ%d\n",polygonNumb,buildingNumb,pointsNumb);
				
			}

			//��ȡ���Թ��ʷ�
			int x,y;
			infile >> x;
			if(x!=3)
				for(int i=0; i<x ;i++)	
				{
					infile >> y;
				}
			else
				for (int i=0;i<((pointsNumb-2)*4-1);i++)
					infile >> x;

			//getline(infile, a, '\n');
			//int triangulate = 0;
			//lenth = a.length();
			//for(int i=0;i<lenth;i++) 
			//	if(a[i]!=' ')	
			//	{
			//		while(i++)
			//		{
			//			if(i == lenth||a[i]==' ')	break;
			//		}
			//		triangulate++;
			//	}
			//if(pointsNumb > 2)
			//{
			//	if(triangulate == 4)
			//		for(int j=0;j<(pointsNumb-2)-1;j++)  getline(infile, a, '\n');
			//	else if(triangulate == 6)
			//		for(int j=0;j<(ceil((double(pointsNumb-2)*4/6))-1);j++)  getline(infile, a, '\n');
			//}


			pointInfo *pointsArray;    
			pointsArray = new pointInfo[pointsNumb];
			for(int j=0;j<pointsNumb;j++)	
			{
				/*std::cout.precision(15);
				getline(infile, a, '\n');
				int cursor=0;
				char b[50];
				int index=0;
				while(a[cursor]==' ')  cursor++;
				while(a[cursor]!=' ')  b[index++]=a[cursor++];
				b[index]='\0';
				pointsArray[j].x = stod(b);
				//printf("%.9lf\t",pointsArray[j].x);
				index = 0;
				while(a[cursor]==' ')  cursor++;
				while(a[cursor]!=' ')  b[index++]=a[cursor++];
				b[index]='\0';
				pointsArray[j].y = stod(b);
				//printf("%.9lf\t",pointsArray[j].y);
				index = 0;
				while(a[cursor]==' ')  cursor++;
				int lenth = 0;
				lenth = a.size();
				while(a[cursor]!=' ')  
				{
					b[index++]=a[cursor++];
					if(cursor == lenth) break;
				}
				b[index]='\0';
				pointsArray[j].z = stod(b);
				//printf("%.9lf\n",pointsArray[j].z);*/
				infile >> pointsArray[j].x;
				infile >> pointsArray[j].y;
				infile >> pointsArray[j].z;
				infile >> a;
				infile >> a;
			}
//			getline(infile, a, '\n');

			double *X,*Y,*Z,*M;
			int partStart[1]={0},partType[1]={SHPT_POLYGONZ};
			pointsNumb++;
			
				X =new double[pointsNumb];
				Y =new double[pointsNumb];
				Z =new double[pointsNumb];
				M =new double[pointsNumb];
				for(int i=0;i<pointsNumb-1;i++)
				{
					X[i] = pointsArray[i].x;
					Y[i] = pointsArray[i].y;
					Z[i] = pointsArray[i].z;
					M[i] = 0;
				}
				X[pointsNumb-1] = pointsArray[0].x;
				Y[pointsNumb-1] = pointsArray[0].y;
				Z[pointsNumb-1] = pointsArray[0].z;
				M[pointsNumb-1] = 0;
				/*if(pointsNumb == 22)
				{
					X[pointsNumb-2] = pointsArray[0].x;
					Y[pointsNumb-2] = pointsArray[0].y;
					Z[pointsNumb-2] = pointsArray[0].z;
					M[pointsNumb-2] = 0;
				}
			}*/
			
			//׷��shp��Ŀ
			SHPObject * oSHP = SHPCreateObject(SHPT_POLYGONZ,polygonNumb-wrongCount,1,partStart,partType,pointsNumb,X,Y,Z,M);
			SHPWriteObject(hSHP,-1,oSHP);
			SHPDestroyObject(oSHP);
			delete[] X;
			delete[] Y;
			delete[] Z;
			delete[] M;

			DBFWriteIntegerAttribute(hDBF,count,0,buildingNumb);
			DBFWriteIntegerAttribute(hDBF,count,1,2);
			DBFWriteStringAttribute(hDBF,count,2,"AUTOMATIC");
			delete[] pointsArray;
			//cout<<"******************************************"<<endl;
			count++;  
			//if(polygonNumb == 2000) break;
		}
		SHPClose(hSHP);
		DBFClose(hDBF);
		//��ԭ��ͶӰ��������.prj�ļ�
		FILE * m_OutPrj;
		m_OutPrj = fopen(Out_Path+"\\OUTPUT.prj","w");
		fprintf(m_OutPrj,"PROJCS[\"WGS_1984\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0000000000000000,298.2572229328696400]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0000000000000000],PARAMETER[\"False_Northing\",0.0000000000000000],PARAMETER[\"Central_Meridian\",117.0000000000002600],PARAMETER[\"Scale_Factor\",1.0000000000000000],PARAMETER[\"Latitude_Of_Origin\",0.0000000000000000],UNIT[\"Meter\",1.0]]\n");
		fclose(m_OutPrj);
	}
}


void CMy3DConvertDlg::OnBnClickedButton1()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	UpdateData(TRUE);	//������������
	if (!m_Radio)//SHP�ļ�ת3d�ļ�
	{
		CString pathName,fileName,fileTitle;
		char * filters = _T("SHP�ļ�(*.shp)|*.shp");
		CFileDialog m_MyOpenDialog(true,NULL,"*.shp",OFN_ALLOWMULTISELECT | OFN_ENABLESIZING | OFN_HIDEREADONLY,filters);
		m_MyOpenDialog.m_ofn.nMaxFile = 500 * MAX_PATH;
		char * ch = new TCHAR[m_MyOpenDialog.m_ofn.nMaxFile];
		m_MyOpenDialog.m_ofn.lpstrFile = ch;			
		ZeroMemory(m_MyOpenDialog.m_ofn.lpstrFile,sizeof(TCHAR) * m_MyOpenDialog.m_ofn.nMaxFile);
		if(m_MyOpenDialog.DoModal() == IDOK)
		{
			POSITION	pos_file;
			pos_file = m_MyOpenDialog.GetStartPosition();
			while(pos_file != NULL)
			{
				pathName = m_MyOpenDialog.GetNextPathName(pos_file);
				ary_fileName.Add(pathName);
				int length = pathName.GetLength();
				for(int i = length-1;i > 0;i--)
				{
					if(pathName.GetAt(i) == '\\')
					{
						fileName = pathName.Right(length-i-1);
						break;
					}
				}
				fileTitle = fileName.Left(fileName.GetLength()-4);
				ary_fileTitle.Add(fileTitle);
			}
		}
		/*int XX = ary_fileTitle.GetSize();
		CString BBB;
		for(int j = 0;j<XX;j++)
		{
			BBB = ary_fileName.GetAt(j);
		}*/
		delete[] ch;
		m_FileName = ary_fileName.GetAt(0);
		UpdateData(FALSE);	//������ʾ
		
	}
	else if (m_Radio)//3d�ļ�תSHP�ļ�
	{
		CString pathName,fileName,fileTitle;
		char * filters = _T("SHP�ļ�(*.3d)|*.3d");
		CFileDialog m_MyOpenDialog(true,NULL,"*.3d",OFN_ALLOWMULTISELECT | OFN_ENABLESIZING | OFN_HIDEREADONLY,filters);
		m_MyOpenDialog.m_ofn.nMaxFile = 500 * MAX_PATH;
		char * ch = new TCHAR[m_MyOpenDialog.m_ofn.nMaxFile];
		m_MyOpenDialog.m_ofn.lpstrFile = ch;			
		ZeroMemory(m_MyOpenDialog.m_ofn.lpstrFile,sizeof(TCHAR) * m_MyOpenDialog.m_ofn.nMaxFile);
		if(m_MyOpenDialog.DoModal() == IDOK)
		{
			POSITION	pos_file;
			pos_file = m_MyOpenDialog.GetStartPosition();
			while(pos_file != NULL)
			{
				pathName = m_MyOpenDialog.GetNextPathName(pos_file);
				ary_fileName.Add(pathName);
				int length = pathName.GetLength();
				for(int i = length-1;i > 0;i--)
				{
					if(pathName.GetAt(i) == '\\')
					{
						fileName = pathName.Right(length-i-1);
						break;
					}
				}
				fileTitle = fileName.Left(fileName.GetLength()-3);
				ary_fileTitle.Add(fileTitle);
			}
		}
		delete[] ch;
		m_FileName = fileName.GetAt(0);
		UpdateData(FALSE);	//������ʾ
	}
}



void CMy3DConvertDlg::OnBnClickedButton2()
	{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	UpdateData(TRUE);	//������������

	BROWSEINFO bi;		//
	memset(&bi, 0, sizeof(BROWSEINFO) );    //��ʼ��BROWSEINFO�ṹbi�ĵ�ַΪȫ"0"
	CString StrFolder;

	StrFolder.Format("D:\\");	//��ʼ·��

	if (!Out_Path.IsEmpty())
	{
		StrFolder = Out_Path;
	}

	int CALLBACK SetSelect(HWND ,  UINT ,  LPARAM ,  LPARAM );	//���ó�ʼѡ��·��
	bi.lpfn = SetSelect;
	bi.lParam =(LPARAM) StrFolder.GetBuffer(StrFolder.GetLength());
	StrFolder.ReleaseBuffer();

	bi.lpszTitle = "��ѡ��Ŀ���ļ���";

	LPITEMIDLIST idl = SHBrowseForFolder(&bi);  //�򿪶Ի���
	if (idl == NULL)
		{
		return;
		}
	SHGetPathFromIDList(idl, StrFolder.GetBuffer(MAX_PATH));
	StrFolder.ReleaseBuffer();
	Out_Path.Format("%s", StrFolder);	//�õ����������ļ���·��
	//////////////////////////////////////////////////////////////////////////

	UpdateData(FALSE);	//������ʾ
	}

int CALLBACK SetSelect(HWND hwnd,  UINT uMsg,  LPARAM lParam,  LPARAM lpData)	//callback�ع�
{
	if (uMsg == BFFM_INITIALIZED)
		{
		::SendMessage(hwnd,BFFM_SETSELECTION,TRUE,lpData);
		}
	return 0;
}

void CMy3DConvertDlg::OnBnClickedOk()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	//CDialogEx::OnOK();
	UpdateData(TRUE);

	FILE * m_InFILE;
	FILE * m_OutFILE;
	if (!m_Radio)
	{
		Out_PathX = Out_Path+"\\SHPת3d�����ļ�.txt";
		corrigendumFILE = fopen(Out_PathX,"w+");
		for(int i=0;i<ary_fileName.GetSize();i++)
		{
			In_Path = ary_fileName.GetAt(i);
			if((m_InFILE=fopen(In_Path,"rb"))==NULL)    
			{
				AfxMessageBox("�����ļ�·��Ϊ�գ�");
				//exit(0);
			}
			else
			{
				In_PathX = In_Path;
				In_PathX = In_PathX.Left(In_PathX.GetLength()-3)+"dbf";
				DBFHandle hDBF = DBFOpen(In_PathX,"rb");
				Out_PathX = Out_Path+"\\"+ary_fileTitle.GetAt(i)+".3d";
				m_OutFILE=fopen(Out_PathX,"w+");
				fprintf(corrigendumFILE,"��%s���ļ�����\n",ary_fileTitle.GetAt(i));
				transferSHPto3D(m_InFILE,hDBF,m_OutFILE);
				fclose(m_InFILE);
			}
		}
		AfxMessageBox("ת�����");
	}
	
	if (m_Radio)//3d�ļ�תSHP�ļ�
	{
		Out_PathX = Out_Path+"\\3dתSHP�����ļ�.txt";
		corrigendumFILE = fopen(Out_PathX,"w+");
		for(int i=0;i<ary_fileName.GetSize();i++)
		{
			In_Path = ary_fileName.GetAt(i);
			if(In_Path == "")   
			{
				AfxMessageBox("�����ļ�·��Ϊ�գ�");
				//exit(0);
			}
			else
			{
				fstream infile(In_Path);
				SHPHandle hSHP = SHPCreate(Out_Path + "\\"+ary_fileTitle.GetAt(i)+".shp",SHPT_POLYGONZ);
				DBFHandle hDBF = DBFCreate(Out_Path +"\\"+ary_fileTitle.GetAt(i)+".dbf");
				fprintf(corrigendumFILE,"��%s���ļ�����\n",ary_fileTitle.GetAt(i));
				transfer3DtoSHP(infile,hSHP,hDBF);	
			}
		}
		AfxMessageBox("ת�����");
	}
	fclose(corrigendumFILE);
		
}






