#ifndef  _3D_VECTOR_H_
#define _3D_VECTOR_H_

//wb

#include "math.h"

const double C3DVector_MinV = 1e-300;

class C3DVector  
{
public:
	C3DVector()	;
	virtual ~C3DVector();

	double x, y, z;

	// 用户构造函数20110907
	C3DVector(double X, double Y, double Z) 
	{ 
		x = X; y = Y; z = Z;
	}

	// 定义矢量之间的'+'法  vector
	C3DVector operator+(C3DVector vVector)
	{
		// 返回结果
		return C3DVector(vVector.x + x, vVector.y + y, vVector.z + z);
	}

	// 定义矢量之间的'-'法 
	C3DVector operator-(C3DVector vVector)
	{
		// 返回矢量相减的结果
		return C3DVector(x - vVector.x, y - vVector.y, z - vVector.z);
	}
	
	// 定义矢量与数的'*'法
	C3DVector operator*(double num)
	{
		// 返回结果
		return C3DVector(x * num, y * num, z * num);
	}

	// 定义矢量与矢量的'*'法
	double operator*(C3DVector v)
	{
		// 返回结果
		return x * v.x + y * v.y + z * v.z ;
	}

	// 定义矢量与矢量的'X'法
	C3DVector operator <(C3DVector v)
	{
		// 返回结果
		return C3DVector( y*v.z - z*v.y , z*v.x - x*v.z , x*v.y -y*v.x );
	}

	// 定义矢量与数的'/'法
	C3DVector operator/(double num)
	{
		// 返回结果
		return C3DVector(x / num, y / num, z / num);
	}

	//求模
	double Magnitude()
	{
		return sqrt( x * x + y * y + z * z );
	}

	// ljl 2013-8-7 13:15:58

	//  两矢量之间的夹角(弧度)

	double CalcAngleByVector(C3DVector v);
	

	// 已知一角度(弧度) 和一个矢量 解另一矢量
	C3DVector CalcOtherVectorByAngle(double dAngle);
	
	//规格化
	int Normalize(double Mag=1);

	int Set(double X, double Y, double Z);



};

#endif 