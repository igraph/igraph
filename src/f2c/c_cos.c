#include "f2c.h"

#ifdef KR_headers
extern double sin(), cos(), sinh(), cosh();

VOID c_cos(r, z) f2c_complex *r, *z;
#else
#undef abs
#include "math.h"
#ifdef __cplusplus
extern "C" {
#endif

void c_cos(f2c_complex *r, f2c_complex *z)
#endif
{
	double zi = z->i, zr = z->r;
	r->r =   cos(zr) * cosh(zi);
	r->i = - sin(zr) * sinh(zi);
	}
#ifdef __cplusplus
}
#endif
