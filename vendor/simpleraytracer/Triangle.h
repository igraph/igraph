/** Triangle.h
*/

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Shape.h"

namespace igraph {

class Triangle : public Shape
{
public:
	Triangle();
	Triangle(const Point& rPoint1, const Point& rPoint2, const Point& rPoint3);
	~Triangle();

	virtual bool Intersect(const Ray& vRay, Point& vIntersectPoint) const;
	virtual Vector Normal(const Point& rSurfacePoint, const Point& rOffSurface) const;

private:
	Point mPoint1, mPoint2, mPoint3;
};

} // namespace igraph

#endif
