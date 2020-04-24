
// 3DConvertDlg.h : ͷ�ļ�
//

#pragma once

#include "resource.h"
#include "PolygonBulging.h"
#include <vector>
#include "pointclockpre.h"
#include "shapefil.h"
#include <stdio.h>
#include <tchar.h>
#include <stdlib.h>
#include <conio.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <sstream>
#include "3DVector.h"
#include "BuildingStruct.h"

using namespace std;

// CMy3DConvertDlg �Ի���

class CMy3DConvertDlg : public CDialogEx
{

// ����
public:
	CMy3DConvertDlg(CWnd* pParent = NULL);	// ��׼���캯��

// �Ի�������
	enum { IDD = IDD_MY3DCONVERT_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV ֧��


// ʵ��
protected:
	HICON m_hIcon;

	// ���ɵ���Ϣӳ�亯��

	//��ʼ���Ի���
	virtual BOOL OnInitDialog();
	//��ϵͳ������
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	//����ɫ
	afx_msg void OnPaint();
	//��ѯʱ�϶�ͼ��
	afx_msg HCURSOR OnQueryDragIcon();
	//������Ϣӳ��
	DECLARE_MESSAGE_MAP()
public:
	//�����ť�������ļ���
	afx_msg void OnBnClickedButton1();
	//ȷ��ת���İ�ť
	afx_msg void OnBnClickedOk();
	//�����ť������ļ���
	afx_msg void OnBnClickedButton2();

	CString In_Path;
	CString In_PathX;
	CString Out_Path;
	CString Out_PathX;
	CString m_FileName;
	BOOL m_Radio;
	//������
	struct pointsArray
	{
		int a;
		int b;
		int c;
		int d;
	};
	//����Ϣ
	struct pointInfo
	{	
		double x;
		double y;
		double z;
	}; 

	pointInfo vectorA,vectorB,normVector,normVectorTemp;
	double_vector m_vctClockX;
	double_vector m_vctClockY;
	double_vector m_vctClockRX;
	double_vector m_vctClockRY;

	//�����ֽ�˳��
	int OnChangeByteOrder(int indata);
	//��ȡ����ε�shp�ļ�
	int OnReadPolygonShp(FILE * m_ShpFile_fp,DBFHandle hDBF,FILE * transferTo3D);
	
	typedef pointsArray Node;
	//ת�Ƶ�����
	bool TransferToDotNum(int Num,double *pX,double *pY,vector<qxwPolygon> vecPoly,vector<vector<int>> &vctNum);
	//���Ƕ�λ
	bool Triangulation(vector<int> vctNum,vector<pointsArray> &vctTriangle);

	void str2int(int &int_temp,const string &string_temp);

	//������Сֵ�����ֵ
	void findMinAndMax(int Num,double *Array,double *MaxMin);
	//�������Ƕ�λ
	void creatTriangulation(int Num,double *PointsX,double *PointsY,double *PointsZ,FILE * transferTo3D,int polygonNumber,int Building);
	//��shpתΪ3D
	void transferSHPto3D(FILE * m_InFILE,DBFHandle hDBF,FILE * m_OutFILE);
	//��3DתΪshp
	void transfer3DtoSHP(fstream &infile,SHPHandle hSHP,DBFHandle hDBF);
	//��XY˳ʱ��
	int CMy3DConvertDlg::ClockWiseXY(int Num,double *PointsX,double *PointsY);
	//��XZ˳ʱ��
	int CMy3DConvertDlg::ClockWiseXZ(int Num,double *PointsX,double *PointsZ);
	//��YZ˳ʱ��
	int CMy3DConvertDlg::ClockWiseYZ(int Num,double *PointsY,double *PointsZ);
	//�ж�����
	void CMy3DConvertDlg::coplanerJudge(int pointsNum,int polygonNum,double *PointsX, double *PointsY,double *PointsZ);
	//���������С����ͷ
	int CMy3DConvertDlg::BuildAnilysize_CulFaceOrigin(_3DFileInfo tmp3DFile,CString error_outputpath);
	//��3D�ļ��ж�ȡ������
	int CMy3DConvertDlg::ReadBuildingsFrom3d(const char * path,_3DFileInfo& _3dFile);
	//�������
	string CMy3DConvertDlg::TrimLeft(string str);
	//��ͷ
	int CMy3DConvertDlg::CulFaceOrigin(double *pX, double *pY, double *pZ, int nNum, CCoorDS &Origin, CCoorDS &CoorX, CCoorDS &CoorY, double dMinDis);
	//����MFAX	
	bool CMy3DConvertDlg::ComputeMFaX(double *X,double *Y,double *Z,int Num,bool IsClockWise,vector<double>&MFX);
	//Ѱ��ͻ����
	bool CMy3DConvertDlg::FindConvexP( int Num,double *X,double *Y,int &ConvexP );

	FILE * corrigendumFILE;
	_3DFileInfo Info3DFile;

	CArray<CString,CString>	ary_fileName;
	CArray<CString,CString>	ary_fileTitle;

	CString m_DH;
	};