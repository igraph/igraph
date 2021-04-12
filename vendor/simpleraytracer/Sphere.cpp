#include "Sphere.h"
#include <math.h>

namespace igraph {

Sphere::Sphere()
{}

Sphere::Sphere(Point vCenter, double vRadius)
{
	Type("Sphere");
	mCenter = vCenter;
	mRadius = vRadius;
}

Sphere::~Sphere()
{
}

bool Sphere::Intersect(const Ray& vRay, Point& vIntersectPoint) const
{
	double c;
	Vector V;
	Vector EO(vRay.Origin(), mCenter);
	double v;
	double disc;
	double d;
	Vector E(Point(0,0,0), vRay.Origin()); // E = vector from origin to ray origin
	Vector P;
	c = mCenter.Distance(vRay.Origin()); //c = distance from eye to center of sphere
	V = vRay.Direction();
	V.NormalizeThis();
	v = EO.Dot(V);
	double v2 = V.Dot(EO.Normalize());
	if (v2 >= 0.0)
	{
		disc = mRadius*mRadius - (EO.Dot(EO) - v*v);
		if (disc <= 0)
			return false;
		else
		{
			d = sqrt(disc);
			P = E + V*(v-d);
			vIntersectPoint = P.ToPoint();
			return true;
		}
	}
	else
		return false;
}

Vector Sphere::Normal(const Point& rSurfacePoint, const Point& rOffSurface) const
{
	// currently does not take rOffSurface point into account,
	// it should check if this point is inside the sphere, if it is
	// return a normal facing the center.

	Vector radius_vector (mCenter, rSurfacePoint);
	return (radius_vector.Normalize());
}

double Sphere::Radius() const
{
	return mRadius;
}
const Point& Sphere::Center() const
{
	return mCenter;
}

} // namespace igraph
