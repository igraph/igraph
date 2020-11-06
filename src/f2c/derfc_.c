#include "f2c.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
#ifndef HAVE_ERFC
extern double erfc();
#endif
double derfc_(x) doublereal *x;
#else
#ifndef HAVE_ERFC
extern double erfc(double);
#endif
double derfc_(doublereal *x)
#endif
{
return( erfc(*x) );
}
#ifdef __cplusplus
}
#endif
