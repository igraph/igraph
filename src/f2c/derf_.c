#include "f2c.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
#ifndef HAVE_ERF
double erf();
#endif
double derf_(x) doublereal *x;
#else
#ifndef HAVE_ERF
extern double erf(double);
#endif
double derf_(doublereal *x)
#endif
{
return( erf(*x) );
}
#ifdef __cplusplus
}
#endif
