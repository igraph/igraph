/** Sphere.h
*/

#ifndef SPHERE_H
#define SPHERE_H

#include "Shape.h"

namespace igraph {

class Sphere : public Shape
{
public:
	Sphere();
	Sphere(Point vCenter, double vRadius);
	~Sphere();

	virtual bool Intersect(const Ray& vRay, Point& vIntersectPoint) const;
	virtual Vector Normal(const Point& rSurfacePoint, const Point& rOffSurface) const;

	double Radius() const;
	const Point& Center() const;

private:
	Point mCenter;
	double mRadius;
};

} // namespace igraph

#endif
