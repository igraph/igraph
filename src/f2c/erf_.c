#include "f2c.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifndef REAL
#define REAL double
#endif

#ifdef KR_headers
#ifndef HAVE_ERF
double erf();
#endif
REAL erf_(x) real *x;
#else
#ifndef HAVE_ERF
extern double erf(double);
#endif
REAL erf_(real *x)
#endif
{
return( erf((double)*x) );
}
#ifdef __cplusplus
}
#endif
