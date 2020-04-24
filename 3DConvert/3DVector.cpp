//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "3DVector.h"


C3DVector::C3DVector()
{
	x=y=z=0;	
}

C3DVector::~C3DVector()
{

}

int C3DVector::Set(double X, double Y, double Z)
{
	x = X;
	y = Y;
	z = Z;
	return 1;
}

//规格化
int C3DVector::Normalize(double Mag)
{
	double mag0,para;

	Mag = (Mag<0)? -1.0*Mag : Mag; 
	mag0 = Magnitude();

//	if( Mag>mag0 && mag0/Mag<0.00001 )
//		return 1;

	if(mag0>Mag)
	{
		para = Mag/mag0;
		x *= para;
		y *= para;
		z *= para;
	}
	else 
	{
		para = mag0/Mag;

		if( fabs(x)>para && para/fabs(x)<C3DVector_MinV )
			return -1;
		if( fabs(y)>para && para/fabs(y)<C3DVector_MinV )
			return -1;
		if( fabs(z)>para && para/fabs(z)<C3DVector_MinV )
			return -1;

		x /= para;
		y /= para;
		z /= para;				
	}

	return 1;
}


double C3DVector::CalcAngleByVector(C3DVector v)
{
	double d1 =0.0;
	double d2 = 0.0;
	d1 = Magnitude();
	d2 = v.Magnitude();

	double dRide =0.0;
	dRide = this->operator *(v);

	double dCos = 0.0;
	dCos = dRide/(d1*d2);

	return acos(dCos);
}

// 已知一角度(弧度) 和一个矢量 解另一矢量
C3DVector C3DVector::CalcOtherVectorByAngle(double dAngle)
{
	double d1 = 0.0;
	d1 = Magnitude();

	double dTan =0.0;
	dTan = tan(dAngle);

	double dGd = 0.0;
	dGd = dTan * d1;

	return C3DVector(0-x,0-y,dGd);
}