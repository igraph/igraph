#include "Point.h"
#include <math.h>

namespace igraph {

Point::Point()
{
	X(0.0);
	Y(0.0);
	Z(0.0);
	Name(0);
}

Point::Point(double vX, double vY, double vZ, int vName)
{
	X(vX);
	Y(vY);
	Z(vZ);
	Name(vName);
}

Point::Point(double vX, double vY, double vZ)
{
	X(vX);
	Y(vY);
	Z(vZ);
	Name(0);
}

Point::~Point()
{}

double Point::X() const
{
	return mX;
}

void Point::X(double vX)
{
	mX = vX;
}

double Point::Y() const
{
	return mY;
}

void Point::Y(double vY)
{
	mY = vY;
}

double Point::Z() const
{
	return mZ;
}

void Point::Z(double vZ)
{
	mZ = vZ;
}

int Point::Name() const
{
	return mName;
}

void Point::Name(int vName)
{
	mName = vName;
}

double Point::Distance(const Point& rPoint) const
{
	return sqrt( (rPoint.X() - mX)*(rPoint.X() - mX) + (rPoint.Y() - mY)*(rPoint.Y() - mY) + (rPoint.Z() - mZ)*(rPoint.Z() - mZ) );
}

bool Point::operator==(const Point& vRhs) const
{
	bool result = true;
/*
	if ( mX + .001 <= vRhs.X() )
		result = false;
	if ( mX - .001 >= vRhs.X() )
		result = false;
	if ( mY + .001 <= vRhs.Y() )
		result = false;
	if ( mY - .001 >= vRhs.Y() )
		result = false;
	if ( mZ + .001 <= vRhs.Z() )
		result = false;
	if ( mZ - .001 >= vRhs.Z() )
		result = false;
*/
	if ( mX != vRhs.X() )
		result = false;
	if ( mY != vRhs.Y() )
		result = false;
	if ( mZ != vRhs.Z() )
		result = false;


	return result;
}

} // namespace igraph
