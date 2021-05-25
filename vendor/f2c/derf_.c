#include "f2c.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
double derf_(x) doublereal *x;
#else
double derf_(doublereal *x)
#endif
{
return( erf(*x) );
}
#ifdef __cplusplus
}
#endif
