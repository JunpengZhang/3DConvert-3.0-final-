
// 3DConvertDlg.h : 头文件
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

// CMy3DConvertDlg 对话框

class CMy3DConvertDlg : public CDialogEx
{

// 构造
public:
	CMy3DConvertDlg(CWnd* pParent = NULL);	// 标准构造函数

// 对话框数据
	enum { IDD = IDD_MY3DCONVERT_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;

	// 生成的消息映射函数

	//初始化对话框
	virtual BOOL OnInitDialog();
	//在系统命令下
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	//上颜色
	afx_msg void OnPaint();
	//查询时拖动图标
	afx_msg HCURSOR OnQueryDragIcon();
	//声明消息映射
	DECLARE_MESSAGE_MAP()
public:
	//浏览按钮（输入文件）
	afx_msg void OnBnClickedButton1();
	//确认转换的按钮
	afx_msg void OnBnClickedOk();
	//浏览按钮（输出文件）
	afx_msg void OnBnClickedButton2();

	CString In_Path;
	CString In_PathX;
	CString Out_Path;
	CString Out_PathX;
	CString m_FileName;
	BOOL m_Radio;
	//点数组
	struct pointsArray
	{
		int a;
		int b;
		int c;
		int d;
	};
	//点信息
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

	//更改字节顺序
	int OnChangeByteOrder(int indata);
	//读取多边形的shp文件
	int OnReadPolygonShp(FILE * m_ShpFile_fp,DBFHandle hDBF,FILE * transferTo3D);
	
	typedef pointsArray Node;
	//转移到点编号
	bool TransferToDotNum(int Num,double *pX,double *pY,vector<qxwPolygon> vecPoly,vector<vector<int>> &vctNum);
	//三角定位
	bool Triangulation(vector<int> vctNum,vector<pointsArray> &vctTriangle);

	void str2int(int &int_temp,const string &string_temp);

	//查找最小值和最大值
	void findMinAndMax(int Num,double *Array,double *MaxMin);
	//创建三角定位
	void creatTriangulation(int Num,double *PointsX,double *PointsY,double *PointsZ,FILE * transferTo3D,int polygonNumber,int Building);
	//将shp转为3D
	void transferSHPto3D(FILE * m_InFILE,DBFHandle hDBF,FILE * m_OutFILE);
	//将3D转为shp
	void transfer3DtoSHP(fstream &infile,SHPHandle hSHP,DBFHandle hDBF);
	//以XY顺时针
	int CMy3DConvertDlg::ClockWiseXY(int Num,double *PointsX,double *PointsY);
	//沿XZ顺时针
	int CMy3DConvertDlg::ClockWiseXZ(int Num,double *PointsX,double *PointsZ);
	//沿YZ顺时针
	int CMy3DConvertDlg::ClockWiseYZ(int Num,double *PointsY,double *PointsZ);
	//判定共面
	void CMy3DConvertDlg::coplanerJudge(int pointsNum,int polygonNum,double *PointsX, double *PointsY,double *PointsZ);
	//建立任意大小的列头
	int CMy3DConvertDlg::BuildAnilysize_CulFaceOrigin(_3DFileInfo tmp3DFile,CString error_outputpath);
	//从3D文件中读取建筑物
	int CMy3DConvertDlg::ReadBuildingsFrom3d(const char * path,_3DFileInfo& _3dFile);
	//向左剪切
	string CMy3DConvertDlg::TrimLeft(string str);
	//列头
	int CMy3DConvertDlg::CulFaceOrigin(double *pX, double *pY, double *pZ, int nNum, CCoorDS &Origin, CCoorDS &CoorX, CCoorDS &CoorY, double dMinDis);
	//计算MFAX	
	bool CMy3DConvertDlg::ComputeMFaX(double *X,double *Y,double *Z,int Num,bool IsClockWise,vector<double>&MFX);
	//寻找突出点
	bool CMy3DConvertDlg::FindConvexP( int Num,double *X,double *Y,int &ConvexP );

	FILE * corrigendumFILE;
	_3DFileInfo Info3DFile;

	CArray<CString,CString>	ary_fileName;
	CArray<CString,CString>	ary_fileTitle;

	CString m_DH;
	};