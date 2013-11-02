/** RayTraceCanvas.h
 */

#ifndef RAY_TRACER_H
#define RAY_TRACER_H


#include<string>
#include "Point.h"
#include "Shape.h"
#include "Color.h"
#include "Light.h"

namespace igraph {

class Image
{
 public: 
  int width, height;
  double *red, *green, *blue, *trans;
};

class RayTracer
{

public:
	RayTracer();
	~RayTracer();

	void RayTrace(Image &result);

	void AddShape(Shape* pShape);
	void AddLight(Light* pLight);

	void BackgroundColor(const Color& rBackgroundColor);
	void EyePoint(const Point& rEyePoint);
	void AmbientColor(const Color& rAmbient);
	void AmbientIntensity(double vAmbientIntensity);

private:
	
	Color Render(const Ray& rRay, bool vIsReflecting = false, const Shape* pReflectingFrom = 0 ); // vEyeRay should be true if the ray we are tracing is a ray from the eye, otherwise it should be false
	Shape* QueryScene(const Ray& rRay, Point& rIntersectionPoint, bool vIsReflecting = false, const Shape* pReflectingFrom = 0);
	double Shade(const Shape* pShapeToShade, const Point& rPointOnShapeToShade);
	double Specular(const Shape* pShapeToShade, const Point& rPointOnShapeToShade, const Light* pLight);	
		
	Color mBackgroundColor;
	Color mAmbientColor;
	Point mEyePoint;
	Color mSpecularColor;
	double mAmbientIntensity;	

	ShapeList* mpShapes;
	LightList* mpLights;

	int mRecursions;
	int mRecursionLimit;
	int mAntiAliasDetail;
};

} // namespace igraph

#endif
