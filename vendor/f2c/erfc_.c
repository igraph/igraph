#include "f2c.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef REAL
#define REAL double
#endif

#ifdef KR_headers
REAL erfc_(x) real *x;
#else
REAL erfc_(real *x)
#endif
{
return( erfc((double)*x) );
}
#ifdef __cplusplus
}
#endif
