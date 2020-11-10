#include "f2c.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
double derfc_(x) doublereal *x;
#else
double derfc_(doublereal *x)
#endif
{
return( erfc(*x) );
}
#ifdef __cplusplus
}
#endif
