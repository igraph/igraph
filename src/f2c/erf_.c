#include "f2c.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef REAL
#define REAL double
#endif

#ifdef KR_headers
REAL erf_(x) real *x;
#else
REAL erf_(real *x)
#endif
{
return( erf((double)*x) );
}
#ifdef __cplusplus
}
#endif
