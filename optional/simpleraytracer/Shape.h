/** Shape.h
 */

#ifndef SHAPE_H
#define SHAPE_H

#include <string>
#include "Color.h"
#include "Ray.h"
#include "Point.h"

#include <list>
using namespace std;

namespace igraph {

class Shape
{
public:
	Shape();
	virtual ~Shape();

	virtual bool Intersect(const Ray& rRay, Point& rIntersectPoint) const = 0;
	virtual Vector Normal(const Point& rSurfacePoint, const Point& rOffSurface) const = 0; 
	// returns a normalized vector that is the normal of this shape from the surface point
	// it also takes the rOffSurface point into account, for example:
	//  if rSurfacePoint is on top of a triangle, then the normal returned will be going up.

	Ray Reflect(const Point& rReflectFrom, const Ray& rRay) const;
	
	void Name(int vName);
	int Name() const;

	const Color& ShapeColor() const;
	void ShapeColor(const Color& rColor);

	double SpecularReflectivity() const;
	void SpecularReflectivity(double rReflectivity);
	double DiffuseReflectivity() const;
	void DiffuseReflectivity(double rReflectivity);
	double AmbientReflectivity() const;
	void AmbientReflectivity(double rReflectivity);

	int SpecularSize() const;
	void SpecularSize(int vSpecularSize);

	const string& Type() const;
	void Type(const string& rType);

private:
	int mName;
	string mType;
	Color mShapeColor;
	double mSpecularReflectivity; // from 0 to 1
	int mSpecularSize; // 1 to 64
	double mDiffuseReflectivity; // from 0 to 1
	double mAmbientReflectivity; // from 0 to 1
};

typedef list<Shape*> ShapeList;
typedef list<Shape*>::iterator ShapeListIterator;

} // namespace igraph

#endif
