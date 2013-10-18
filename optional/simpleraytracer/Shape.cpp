#include "Shape.h"
#include "unit_limiter.h"

namespace igraph {

Shape::Shape()
{
	mName = 0;
	mAmbientReflectivity = .6;
	mSpecularReflectivity = 0;
	mDiffuseReflectivity = 0;
	mSpecularSize = 64;
}

Shape::~Shape()
{}

int Shape::Name() const
{
	return mName;
}

void Shape::Name(int vName)
{
	mName = vName;
}

const Color& Shape::ShapeColor() const
{
	return mShapeColor;
}

void Shape::ShapeColor(const Color& rColor)
{
	mShapeColor = rColor;
}

double Shape::AmbientReflectivity() const
{
	return mAmbientReflectivity;
}
double Shape::SpecularReflectivity() const
{
	return mSpecularReflectivity;
}
double Shape::DiffuseReflectivity() const
{
	return mDiffuseReflectivity;
}

void Shape::AmbientReflectivity(double rReflectivity)
{
	mAmbientReflectivity = unit_limiter(rReflectivity);
}
void Shape::SpecularReflectivity(double rReflectivity)
{
	mSpecularReflectivity = unit_limiter(rReflectivity);
}
void Shape::DiffuseReflectivity(double rReflectivity)
{
	mDiffuseReflectivity = unit_limiter(rReflectivity);
}

Ray Shape::Reflect(const Point& rReflectFrom, const Ray& rIncidentRay) const
{
	Ray result; // the reflected ray
	Vector result_direction; // the reflected direction vector
	Vector incident_unit = rIncidentRay.Direction().Normalize();
	Vector normal = this->Normal(rReflectFrom, rIncidentRay.Origin() );
	if ( !normal.IsSameDirection(incident_unit) )
		normal.ReverseDirection(); // we want the normal in the same direction of the incident ray.

	result.Origin(rReflectFrom);
	result.Direction( normal*2.0*normal.Dot(incident_unit) - incident_unit );
/*
	if ( normal.Dot(rIncidentRay.Direction().Normalize()) < 0.0 )
		normal.ReverseDirection();

	result.Origin(rReflectFrom);
	result.Direction((normal*2.0) - rIncidentRay.Direction().Normalize());
*/

	return result;
}

const string& Shape::Type() const
{
	return mType;
}

void Shape::Type(const string& rType)
{
	mType = rType;
}

int Shape::SpecularSize() const
{
	return mSpecularSize;
}

void Shape::SpecularSize(int vSpecularSize)
{
	mSpecularSize = vSpecularSize;
}

} // namespace igraph
