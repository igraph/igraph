#include "RayVector.h"
#include <math.h>

namespace igraph {

Vector::Vector()
{
	mI = mJ = mK = 0.0;
}

Vector::Vector(const Point& vStartPoint, const Point& vEndPoint)
{
	mI = vEndPoint.X() - vStartPoint.X();
	mJ = vEndPoint.Y() - vStartPoint.Y();
	mK = vEndPoint.Z() - vStartPoint.Z();
}

Vector::Vector(double vI, double vJ, double vK)
{
	mI = vI;
	mJ = vJ;
	mK = vK;
}

Vector::~Vector()
{}

// returns a unit vector of this vector
Vector Vector::Normalize() const
{
	double magnitude = Magnitude();
	return Vector(mI/magnitude, mJ/magnitude, mK/magnitude);
}

void Vector::NormalizeThis()
{
	*this = Normalize();
}

void Vector::ReverseDirection()
{
	*this = *this * -1.0;
}

bool Vector::IsSameDirection(const Vector& rVector) const
{
	return ( this->Normalize().Dot(rVector.Normalize()) > 0.0 );
}


void Vector::I(double vI)
{
	mI = vI;
}

double Vector::I() const
{
	return mI;
}

void Vector::J(double vJ)
{
	mJ = vJ;
}

double Vector::J() const
{
	return mJ;
}
void Vector::K(double vK)
{
	mK = vK;
}

double Vector::K() const
{
	return mK;
}

// returns the dot product of this and rVector
double Vector::Dot(const Vector& rVector) const
{
	return mI*rVector.I() + mJ*rVector.J() + mK*rVector.K();
}

// returns the cross product of this and vVector
Vector Vector::Cross(const Vector& rVector) const
{
	return Vector(mJ*rVector.K() - rVector.J()*mK, -1.0*(mI*rVector.K() - rVector.I()*mK), mI*rVector.J() - rVector.I()*mJ);
}

// returns the sum of this vector with another vector
Vector Vector::operator+ (Vector vRhs) const
{
	return Vector(mI + vRhs.I(), mJ + vRhs.J(), mK + vRhs.K());
}

// returns the sume of a vector and a Point
Point Vector::operator+ (Point vRhs) const
{
	return Point(mI + vRhs.X(), mJ + vRhs.Y(), mK + vRhs.Z());
}

// returns the difference of two vectors
Vector Vector::operator- (Vector vRhs) const
{
	return Vector(mI-vRhs.I(), mJ-vRhs.J(), mK-vRhs.K());
}

// returns multiplication of a scalar with this vector
Vector Vector::operator* (double vRhs) const
{
	return Vector(mI*vRhs, mJ*vRhs, mK*vRhs);
}

// converts this vector to a point
Point Vector::ToPoint() const
{
	return Point(mI,mJ,mK);
}

// returns the magnitude
double Vector::Magnitude() const
{
	return sqrt(mI*mI + mJ*mJ + mK*mK);
}

} // namespace igraph
