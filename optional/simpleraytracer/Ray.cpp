#include "Ray.h"

namespace igraph {

Ray::Ray()
{}

Ray::~Ray()
{}

Ray::Ray(const Point& rOrigin, const Vector& rDirection)
{
	Direction(rDirection);
	Origin(rOrigin);
	
}

Ray::Ray(const Point& rOrigin, const Point& rEndPoint)
{
	Direction(Vector(rOrigin,rEndPoint));
	Origin(rOrigin);
}

const Point& Ray::Origin() const
{
	return mOrigin;
}

void Ray::Origin(Point vOrigin)
{
	mOrigin = vOrigin;
}

const Vector& Ray::Direction() const
{
	return mDirection;
}

void Ray::Direction(Vector vDirection)
{
	mDirection = vDirection;
}

} // namespace igraph
