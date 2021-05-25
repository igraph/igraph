#include "f2c.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
double r_imag(z) f2c_complex *z;
#else
double r_imag(f2c_complex *z)
#endif
{
return(z->i);
}
#ifdef __cplusplus
}
#endif
