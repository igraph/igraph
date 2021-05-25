#include "f2c.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
VOID pow_ci(p, a, b) 	/* p = a**b  */
 f2c_complex *p, *a; integer *b;
#else
extern void pow_zi(doublecomplex*, doublecomplex*, integer*);
void pow_ci(f2c_complex *p, f2c_complex *a, integer *b) 	/* p = a**b  */
#endif
{
doublecomplex p1, a1;

a1.r = a->r;
a1.i = a->i;

pow_zi(&p1, &a1, b);

p->r = p1.r;
p->i = p1.i;
}
#ifdef __cplusplus
}
#endif
