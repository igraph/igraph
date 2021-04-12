#include "Light.h"
#include "unit_limiter.h"

namespace igraph {

Light::Light() : mLightPoint(0,0,0)
{
	mIntensity = 0.1;
}

Light::Light(const Point& rLightPoint) : mLightPoint(rLightPoint)
{
	mIntensity = 0.1;
}

Light::~Light()
{}

const Point& Light::LightPoint() const
{
	return mLightPoint;
}
void Light::LightPoint(const Point& rLightPoint)
{
	mLightPoint = rLightPoint;
}
double Light::Intensity() const
{
	return mIntensity;
}
void Light::Intensity(double vIntensity)
{
	mIntensity = unit_limiter(vIntensity);
}

const Color& Light::LightColor() const
{
	return mLightColor;
}

void Light::LightColor(const Color& rLightColor)
{
	mLightColor = rLightColor;
}

} // namespace igraph
