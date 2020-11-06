#include "f2c.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifndef REAL
#define REAL double
#endif

#ifdef KR_headers
#ifndef HAVE_ERFC
double erfc();
#endif
REAL erfc_(x) real *x;
#else
#ifndef HAVE_ERFC
extern double erfc(double);
#endif
REAL erfc_(real *x)
#endif
{
return( erfc((double)*x) );
}
#ifdef __cplusplus
}
#endif
