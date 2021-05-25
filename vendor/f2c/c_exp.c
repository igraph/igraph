#include "f2c.h"

#ifdef KR_headers
extern double exp(), cos(), sin();

 VOID c_exp(r, z) f2c_complex *r, *z;
#else
#undef abs
#include "math.h"
#ifdef __cplusplus
extern "C" {
#endif

void c_exp(f2c_complex *r, f2c_complex *z)
#endif
{
	double expx, zi = z->i;

	expx = exp(z->r);
	r->r = expx * cos(zi);
	r->i = expx * sin(zi);
	}
#ifdef __cplusplus
}
#endif
