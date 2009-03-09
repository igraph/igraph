#include "config.h"
#include "f2c.h"

#ifdef KR_headers
double igraphd_sign(a,b) doublereal *a, *b;
#else
double igraphd_sign(doublereal *a, doublereal *b)
#endif
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}
