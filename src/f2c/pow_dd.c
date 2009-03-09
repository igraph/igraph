#include "config.h"
#include "f2c.h"

#ifdef KR_headers
double pow();
double igraphpow_dd(ap, bp) doublereal *ap, *bp;
#else
#undef abs
#include "math.h"
double igraphpow_dd(doublereal *ap, doublereal *bp)
#endif
{
return(pow(*ap, *bp) );
}
