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

	// �û����캯��20110907
	C3DVector(double X, double Y, double Z) 
	{ 
		x = X; y = Y; z = Z;
	}

	// ����ʸ��֮���'+'��  vector
	C3DVector operator+(C3DVector vVector)
	{
		// ���ؽ��
		return C3DVector(vVector.x + x, vVector.y + y, vVector.z + z);
	}

	// ����ʸ��֮���'-'�� 
	C3DVector operator-(C3DVector vVector)
	{
		// ����ʸ������Ľ��
		return C3DVector(x - vVector.x, y - vVector.y, z - vVector.z);
	}
	
	// ����ʸ��������'*'��
	C3DVector operator*(double num)
	{
		// ���ؽ��
		return C3DVector(x * num, y * num, z * num);
	}

	// ����ʸ����ʸ����'*'��
	double operator*(C3DVector v)
	{
		// ���ؽ��
		return x * v.x + y * v.y + z * v.z ;
	}

	// ����ʸ����ʸ����'X'��
	C3DVector operator <(C3DVector v)
	{
		// ���ؽ��
		return C3DVector( y*v.z - z*v.y , z*v.x - x*v.z , x*v.y -y*v.x );
	}

	// ����ʸ��������'/'��
	C3DVector operator/(double num)
	{
		// ���ؽ��
		return C3DVector(x / num, y / num, z / num);
	}

	//��ģ
	double Magnitude()
	{
		return sqrt( x * x + y * y + z * z );
	}

	// ljl 2013-8-7 13:15:58

	//  ��ʸ��֮��ļн�(����)

	double CalcAngleByVector(C3DVector v);
	

	// ��֪һ�Ƕ�(����) ��һ��ʸ�� ����һʸ��
	C3DVector CalcOtherVectorByAngle(double dAngle);
	
	//���
	int Normalize(double Mag=1);

	int Set(double X, double Y, double Z);



};

#endif 