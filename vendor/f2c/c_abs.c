#include "f2c.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
extern double f__cabs();

double c_abs(z) f2c_complex *z;
#else
extern double f__cabs(double, double);

double c_abs(f2c_complex *z)
#endif
{
return( f__cabs( z->r, z->i ) );
}
#ifdef __cplusplus
}
#endif
