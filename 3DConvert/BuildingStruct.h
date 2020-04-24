#pragma once
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

struct _RGB
{
	int iR;
	int iG;
	int iB;
	_RGB()
	{
		iR= 0;
		iG= 0;
		iB= 0;
	}

	_RGB(const _RGB &p)
	{
		iR = p.iR;
		iG = p.iG;
		iB = p.iB;
	}
	_RGB(int R, int G, int B,int buildID)
	{
		iR = R;
		iG = G;
		iB = B;
	}
	void operator=(const _RGB & p)
	{
		iR = p.iR;
		iG = p.iG;
		iB = p.iB;
	}
};
struct Point{

	double x;
	double y;
	double z;

	Point(double a,double b,double c):x(a),y(b),z(c){}

	Point()
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}

	void operator=(const Point &p)
	{
		x = p.x;
		y = p.y;
		z = p.z;
	}

	ostream& operator<<(ostream& output)
	{
		output << x<<y<<z<<endl;
		return output;
	}
};
///三维点
struct point3d
{
	double x;
	double y;
	double z;
	double nodeCol4;
	double nodeCol5;
	point3d()
	{
		nodeCol5=nodeCol4=x=y=z=.0;
	}
	point3d(const point3d &p)
	{
		x = p.x;
		y = p.y;
		z = p.z;
		nodeCol4=p.nodeCol4;
		nodeCol5=p.nodeCol5;
	}
	point3d(double X, double Y, double Z)
	{
		x = X;
		y = Y;
		z = Z;
	}
	void operator=(const point3d & p)
	{
		x = p.x;
		y = p.y;
		z = p.z;
	}
	point3d operator+(const point3d& p) const
	{
		return point3d(x+p.x, y+p.y, z+p.z);
	}

	point3d operator-(const point3d& p) const
	{
		return point3d(x-p.x, y-p.y, z-p.z);
	}
};
struct polygon
{
	int nodeNum;//面角点个数
	int id;//面ID
	_RGB PolygonRGBs;//添加RGB by zht 2015 06 02
	int faceLn1Col6;//每个面模块第一行，第六列值
	int faceLn1Col7;//每个面模块第一行，第七列值

	double faceLn1Col8;//每个面模块第一行，第8列值
	int faceLn1Col9;//每个面模块第一行，第9列值
	int faceLn1Col10;//每个面模块第一行，第10列值
	int faceLn1Col11;//每个面模块第一行，第11列值
	int faceLn1Col12;//每个面模块第一行，第12列值
	int faceLn1Col13;//每个面模块第一行，第13列值
	int faceLn1Col14;//每个面模块第一行，第14列值
	int faceLn1Col15;//每个面模块第一行，第15列值
	int faceLn1Col16;//每个面模块第一行，第16列值
	int faceLn1Col17;//每个面模块第一行，第17列值
	int faceLn1Col18;//每个面模块第一行，第18列值
	int faceLn1Col19;//每个面模块第一行，第18列值
	double faceLn1Col20;//每个面模块第一行，第18列值
	double faceLn1Col21;//每个面模块第一行，第18列值
	vector<int>subFace;//每个面第三行四行的信息

	vector<point3d> points;
	polygon()
	{
		points.clear();
	}
	/*polygon(const polygon & p)
	{
		id = p.id;
		for (vector<point3d>::size_type i=0; i<p.points.size(); ++i)
		{
			point3d point = p.points[i];
			points.push_back(point);
		}
	}*/
};

///建筑物
struct Building 
{
	string id;///建筑物名称
	vector<polygon> polygons;///每个面
	/*	vector<int> m_Index;///每个面对应的Index*/

	//添加RGB by zht 2015 06 02
	void operator=(const Building & p)
	{
		polygons.clear();
		id = p.id;
		for (vector<polygon>::size_type i=0; i<p.polygons.size(); ++i)
		{
			polygon _polygon = p.polygons[i];
			polygons.push_back(_polygon);
		}
	}
};

//3d文件结构体信息，部分参数由于不知物理含义采用行列格式命名
struct _3DFileInfo
{
	string versionNum;//版本号
	int Ln2;//第二行参数
	string Meter;//固定样式
	vector<double>faceMaxMinGs;//第五行表示最小最大高斯坐标XY
	int layerNum;//层总数
	string layerName;//层的名
	int layerShow;//该层是否显示
	int layerLineNum;//该层线条数
	int modNum;//模块数
	int Ln11;//十一行
	int Ln12;//12行
	int Ln13;//13行

	vector<Building>buildInfo;//建筑物信息
};



struct Build
{
	string id;///建筑物名称
	vector<Point> Points;
	void operator=(const Build & p)
	{
		id = p.id;
		for (vector<polygon>::size_type i=0; i<p.Points.size(); ++i)
		{
			Point pt = p.Points[i];
			Points.push_back(pt);
		}
	}
};
typedef struct 
{
	double dx;
	double dy;
	double dz;
}CCoorDS;
struct XMLPOINT 
{
	CString sX;
	CString sY;
	CString sZ;
	CString sU;
	CString sV;
	XMLPOINT()
	{
		sX = "";
		sY = "";
		sZ = "";
		sU = "";
		sV = "";
	}
};