#include "RayTracer.h"
#include "unit_limiter.h"
#include <fstream>
#include <iostream>

namespace igraph {

RayTracer::RayTracer() : mBackgroundColor(0,0,0,0), mAmbientColor(0,0,0), mEyePoint(0,0,0), mSpecularColor(1,1,1)
{
	// begin settings
	mAmbientIntensity = .7;
	mRecursionLimit = 700;
	mAntiAliasDetail = 1;
	// end settings

	mRecursions = 0;

	mpShapes = new ShapeList;
	mpLights = new LightList;
}

RayTracer::~RayTracer()
{
	ShapeListIterator iter1 = mpShapes->begin();
	while ( iter1 != mpShapes->end() )
	{
	  delete *iter1;
	  iter1++;
	}
	delete mpShapes;

	LightListIterator iter2 = mpLights->begin();
	while ( iter2 != mpLights->end() )
	{
	  delete *iter2;
	  iter2++;
	}

	delete mpLights;

}

void RayTracer::RayTrace(Image &result)
{
        int mWidth=result.width;
	int mHeight=result.height;
	Ray eye_ray(mEyePoint,Vector(0,0,1));
	Color draw_color;
	double i_inc, j_inc, anti_alias_i_inc, anti_alias_j_inc; // amount to increment the ray in each direction
	double i, j, anti_alias_i, anti_alias_j; // the i and j values of the ray
	int pixel_x, pixel_y, anti_alias_pixel_x, anti_alias_pixel_y; // the pixels being drawn

	double average_red_byte, average_green_byte, average_blue_byte, average_trans_byte;
	int anti_alias_count; // the number of anti aliases (used in averaging)
	int idx=0;

	i_inc = 2.0/(double)mWidth;
	j_inc = 2.0/(double)mHeight;
	anti_alias_i_inc = 1.0/(double)mAntiAliasDetail;
	anti_alias_j_inc = 1.0/(double)mAntiAliasDetail;

	pixel_y = 0;
	j = 1.0;
	for (; pixel_y < mHeight; j -= j_inc, pixel_y++)
	{
		pixel_x = 0;
		i = -1.0;
		for (; pixel_x < mWidth; i += i_inc, pixel_x++)
		{
			anti_alias_pixel_y = 0;
			anti_alias_j = 0.0;
 			average_red_byte = 0;
			average_green_byte = 0;
			average_blue_byte = 0;
			average_trans_byte = 0;
			anti_alias_count = 0;
			for (; anti_alias_pixel_y < mAntiAliasDetail; anti_alias_j += anti_alias_j_inc, anti_alias_pixel_y++)
			{
				anti_alias_pixel_x = 0;
				anti_alias_i = 0.0;

				for (; anti_alias_pixel_x < mAntiAliasDetail; anti_alias_i += anti_alias_i_inc, anti_alias_pixel_x++)
				{
					anti_alias_count++;
					eye_ray.Direction( Vector(i+(anti_alias_i*i_inc),j+(anti_alias_j*j_inc),1.0) );
					draw_color = Render(eye_ray);

					average_red_byte = average_red_byte + ((double)draw_color.RedByte() - average_red_byte)/(double)anti_alias_count;
					average_green_byte = average_green_byte + ((double)draw_color.GreenByte() - average_green_byte)/(double)anti_alias_count;
					average_blue_byte = average_blue_byte + ((double)draw_color.BlueByte() - average_blue_byte)/(double)anti_alias_count;
					average_trans_byte = average_trans_byte + ((double)draw_color.TransparentByte() - average_trans_byte)/(double)anti_alias_count;
				}
			}
			
			result.red  [idx] = average_red_byte/255;
			result.green[idx] = average_green_byte/255;
			result.blue [idx] = average_blue_byte/255;
			result.trans[idx] = average_trans_byte/255;
			idx++;
		}
	}
}

Color RayTracer::Render(const Ray& rRay, bool vIsReflecting, const Shape* pReflectingFrom )
{
	mRecursions++;
	Shape* closest_shape;
	Point intersect_point;
	Color result;
	if (vIsReflecting)
		closest_shape = QueryScene(rRay, intersect_point, vIsReflecting, pReflectingFrom);
	else
		closest_shape = QueryScene(rRay, intersect_point);

	if (closest_shape == NULL && !vIsReflecting)
	{
		mRecursions = 0;
		return mBackgroundColor;
	}
	if (closest_shape == NULL && vIsReflecting)
	{
		mRecursions = 0;
		return mAmbientColor*mAmbientIntensity;
	}
	if ( mRecursions > mRecursionLimit )
	{
		mRecursions = 0;
		return Color(0,0,0); // mAmbientColor*mAmbientIntensity;
	}
	result = closest_shape->ShapeColor()*Shade(closest_shape, intersect_point);

	Ray backwards_ray(intersect_point,rRay.Direction()*-1);
	if ( closest_shape->DiffuseReflectivity() > 0.0 )
		result = result + (Render( closest_shape->Reflect(intersect_point,backwards_ray), true, closest_shape )*closest_shape->DiffuseReflectivity());

	return (result + mSpecularColor);
}

double RayTracer::Shade(const Shape* pShapeToShade, const Point& rPointOnShapeToShade)
{
	double intensity = mAmbientIntensity * pShapeToShade->AmbientReflectivity(); // the ambient intensity of the scene
	Ray light_ray; // the ray that goes from the intersection point to the light sources
	double dot_product;
	Shape* closest_shape; // the shape closest from the intersection point to the light source
	Point light_intersect; // the intersection point of the ray that goes from the intersection point to the light source 
	light_ray.Origin(rPointOnShapeToShade); // lightRay. org= object. intersect;
	Ray light_ray_from_light;

	LightListIterator iter = mpLights->begin();
	mSpecularColor.Red(0);
	mSpecularColor.Green(0);
	mSpecularColor.Blue(0);

	while ( iter != mpLights->end() ) // foreach light in LightList do
	{

		light_ray.Direction(Vector(rPointOnShapeToShade,(*iter)->LightPoint())); // lightRay. dir= light. dir
		light_ray_from_light.Origin((*iter)->LightPoint());
		light_ray_from_light.Direction(Vector((*iter)->LightPoint(),rPointOnShapeToShade));

		closest_shape = QueryScene(light_ray_from_light, light_intersect);
		if ( closest_shape == NULL || (closest_shape == pShapeToShade && light_ray.Direction().Dot(pShapeToShade->Normal(rPointOnShapeToShade, light_ray_from_light.Origin() )) >= 0.0 )  ) //if (QueryScene( lightRay)= NIL)
		{
			Vector normal_vector = pShapeToShade->Normal(rPointOnShapeToShade, Point() );
			dot_product = normal_vector.Dot(light_ray.Direction().Normalize());
			dot_product *= (*iter)->Intensity();

			if (dot_product < 0.0)
			{
				if (pShapeToShade->Type() == "Triangle")
					dot_product = dot_product*-1.0;
				else
					dot_product = 0.0;
			}
			intensity = unit_limiter( intensity + dot_product );

			if ( light_ray.Direction().Dot(pShapeToShade->Normal(rPointOnShapeToShade, light_ray_from_light.Origin() )) >= 0.0 )
			{
				double specular = Specular(pShapeToShade, rPointOnShapeToShade, *iter);
				mSpecularColor = mSpecularColor + Color(specular,specular,specular);
			}
		}

		iter++;
	}

	return intensity;
}

double RayTracer::Specular(const Shape* pShapeToShade, const Point& rPointOnShapeToShade, const Light* pLight)
{
	Ray reflected = pShapeToShade->Reflect(rPointOnShapeToShade,Ray(rPointOnShapeToShade, pLight->LightPoint()));

	Vector eye_vector(rPointOnShapeToShade, mEyePoint);
	Vector reflected_vector = reflected.Direction().Normalize();
	eye_vector.NormalizeThis();
	double dot_product = eye_vector.Dot(reflected_vector);	
	
	int n = pShapeToShade->SpecularSize();
	double specular_intensity = dot_product/(n - n*dot_product+ dot_product);
	return unit_limiter(specular_intensity*pLight->Intensity());
}

Shape* RayTracer::QueryScene(const Ray& rRay, Point& rIntersectionPoint, bool vIsReflecting, const Shape* pReflectingFrom)
{
	Shape* closest_shape = NULL;
	Point intersect_point;
	double closest_distance;
	double intersect_distance;
	bool found_intersection = false;

	ShapeListIterator iter = mpShapes->begin();
	while ( iter != mpShapes->end() )
	{
		if ( (*iter)->Intersect( rRay, intersect_point ) )
		{
			intersect_distance =  intersect_point.Distance(rRay.Origin());
			if ( !found_intersection && (*iter) != pReflectingFrom)
			{
				found_intersection = true;
				rIntersectionPoint = intersect_point;
				closest_shape = *iter;
				closest_distance = intersect_distance;
			}
			else if ( intersect_distance < closest_distance && (*iter) != pReflectingFrom )
			{
				rIntersectionPoint = intersect_point;
				closest_shape = *iter;
				closest_distance = intersect_distance;
			}
		}
		iter++;
	}
	
	return closest_shape;
}

void RayTracer::AddShape(Shape* pShape)
{
	// should check if a shape with the same name already exists
	mpShapes->push_back(pShape);
}
void RayTracer::AddLight(Light* pLight)
{
	// should check if a shape with the same name already exists
	mpLights->push_back(pLight);
}

void RayTracer::BackgroundColor(const Color& rBackgroundColor)
{
	mBackgroundColor = rBackgroundColor;
}
void RayTracer::EyePoint(const Point& rEyePoint)
{
	mEyePoint = rEyePoint;
}
void RayTracer::AmbientColor(const Color& rAmbientColor)
{
	mAmbientColor = rAmbientColor;
}
void RayTracer::AmbientIntensity(double vAmbientIntensity)
{
	mAmbientIntensity = unit_limiter(vAmbientIntensity);
}

} // namespace igraph
