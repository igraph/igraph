#include "stdio.h"
#include "f2c.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
extern void sig_die(void);

int abort_(void)
#else
extern void sig_die(const char*,int);

int abort_(void)
#endif
{
sig_die("Fortran abort routine called", 1);
return 0;	/* not reached */
}
#ifdef __cplusplus
}
#endif
