#include "Triangle.h"
#include <math.h>

namespace igraph {

Triangle::Triangle()
{}

Triangle::Triangle(const Point& rPoint1, const Point& rPoint2, const Point& rPoint3)
{
	Type("Triangle");
	mPoint1 = rPoint1;
	mPoint2 = rPoint2;
	mPoint3 = rPoint3;
}

Triangle::~Triangle()
{
}

bool Triangle::Intersect(const Ray& vRay, Point& rIntersectPoint) const
{

	Vector pointb_minus_pointa (mPoint1, mPoint2);
	Vector pointb_minus_pointc (mPoint1, mPoint3);
/*
	Vector plane_normal = pointb_minus_pointa.Cross(pointb_minus_pointc);	

	// get the plane normal facing the right way:
	Vector plane_normal_normalized = plane_normal.Normalize();
	Vector triangle_to_ray_origin = Vector(mPoint1, vRay.Origin() );
	triangle_to_ray_origin.NormalizeThis();
	if ( plane_normal_normalized.Dot(triangle_to_ray_origin) < 0.0 )
	{
		plane_normal = plane_normal * -1.0;
		plane_normal_normalized = plane_normal_normalized * -1.0;
	}

	// check that the ray is actually facing the triangle
	Vector ray_direction_normalized = vRay.Direction().Normalize();
	if ( plane_normal_normalized.Dot(ray_direction_normalized) > 0.0 )
		return false;
*/
	Vector plane_normal = this->Normal(mPoint1, vRay.Origin());
	Vector ray_direction_normalized = vRay.Direction().Normalize();
	if ( plane_normal.IsSameDirection(ray_direction_normalized) )
		return false;

	Vector b_minus_u (vRay.Origin(), mPoint2);

	double t = plane_normal.Dot(b_minus_u) / plane_normal.Dot(vRay.Direction());
	Point p =  (vRay.Direction() * t) + vRay.Origin();

	Vector p_minus_a (mPoint1, p);
	Vector p_minus_b (mPoint2, p);
	Vector p_minus_c (mPoint3, p);
	Vector pointc_minus_pointb (mPoint2, mPoint3);
	Vector pointa_minus_pointc (mPoint3, mPoint1);

	double test1 = (pointb_minus_pointa.Cross(p_minus_a)).Dot(plane_normal);
	double test2 = (pointc_minus_pointb.Cross(p_minus_b)).Dot(plane_normal);
	double test3 = (pointa_minus_pointc.Cross(p_minus_c)).Dot(plane_normal);

	if ((test1 > 0 && test2 > 0 && test3 > 0) || (test1 < 0 && test2 < 0 && test3 < 0))
	{
		rIntersectPoint = p;
		return true; 
	}
	else
		return false;
}

Vector Triangle::Normal(const Point& rSurfacePoint, const Point& rOffSurface) const
{
	Vector pointb_minus_pointa (mPoint1, mPoint2);
	Vector pointb_minus_pointc (mPoint1, mPoint3);
	Vector plane_normal = pointb_minus_pointa.Cross(pointb_minus_pointc).Normalize();	

	// get the plane normal facing the right way:
	Vector triangle_to_off_surface_point = Vector(mPoint1, rOffSurface );
	triangle_to_off_surface_point.NormalizeThis();
	if ( !plane_normal.IsSameDirection(triangle_to_off_surface_point) )
	{
		plane_normal.ReverseDirection();
	}

	return plane_normal;
}

} // namespace igraph
