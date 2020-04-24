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
///��ά��
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
	int nodeNum;//��ǵ����
	int id;//��ID
	_RGB PolygonRGBs;//���RGB by zht 2015 06 02
	int faceLn1Col6;//ÿ����ģ���һ�У�������ֵ
	int faceLn1Col7;//ÿ����ģ���һ�У�������ֵ

	double faceLn1Col8;//ÿ����ģ���һ�У���8��ֵ
	int faceLn1Col9;//ÿ����ģ���һ�У���9��ֵ
	int faceLn1Col10;//ÿ����ģ���һ�У���10��ֵ
	int faceLn1Col11;//ÿ����ģ���һ�У���11��ֵ
	int faceLn1Col12;//ÿ����ģ���һ�У���12��ֵ
	int faceLn1Col13;//ÿ����ģ���һ�У���13��ֵ
	int faceLn1Col14;//ÿ����ģ���һ�У���14��ֵ
	int faceLn1Col15;//ÿ����ģ���һ�У���15��ֵ
	int faceLn1Col16;//ÿ����ģ���һ�У���16��ֵ
	int faceLn1Col17;//ÿ����ģ���һ�У���17��ֵ
	int faceLn1Col18;//ÿ����ģ���һ�У���18��ֵ
	int faceLn1Col19;//ÿ����ģ���һ�У���18��ֵ
	double faceLn1Col20;//ÿ����ģ���һ�У���18��ֵ
	double faceLn1Col21;//ÿ����ģ���һ�У���18��ֵ
	vector<int>subFace;//ÿ������������е���Ϣ

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

///������
struct Building 
{
	string id;///����������
	vector<polygon> polygons;///ÿ����
	/*	vector<int> m_Index;///ÿ�����Ӧ��Index*/

	//���RGB by zht 2015 06 02
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

//3d�ļ��ṹ����Ϣ�����ֲ������ڲ�֪������������и�ʽ����
struct _3DFileInfo
{
	string versionNum;//�汾��
	int Ln2;//�ڶ��в���
	string Meter;//�̶���ʽ
	vector<double>faceMaxMinGs;//�����б�ʾ��С����˹����XY
	int layerNum;//������
	string layerName;//�����
	int layerShow;//�ò��Ƿ���ʾ
	int layerLineNum;//�ò�������
	int modNum;//ģ����
	int Ln11;//ʮһ��
	int Ln12;//12��
	int Ln13;//13��

	vector<Building>buildInfo;//��������Ϣ
};



struct Build
{
	string id;///����������
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