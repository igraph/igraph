#include "stdio.h"
#include "f2c.h"
#include "igraph_error.h"

#ifdef KR_headers
extern void f_exit();
int s_stop(s, n) char *s; ftnlen n;
#else
#undef abs
#undef min
#undef max
#include "stdlib.h"
#ifdef __cplusplus
extern "C" {
#endif
#ifdef __cplusplus
extern "C" {
#endif
void f_exit(void);

int s_stop(char *s, ftnlen n)
#endif
{
	IGRAPH_FATAL("STOP statement executed from f2c code");
	/*
int i;

if(n > 0)
	{
	fprintf(stderr, "STOP ");
	for(i = 0; i<n ; ++i)
		putc(*s++, stderr);
	fprintf(stderr, " statement executed\n");
	}
#ifdef NO_ONEXIT
f_exit();
#endif
exit(0);
*/

/* We cannot avoid (useless) compiler diagnostics here:		*/
/* some compilers complain if there is no return statement,	*/
/* and others complain that this one cannot be reached.		*/

return 0; /* NOT REACHED */
}
#ifdef __cplusplus
}
#endif
#ifdef __cplusplus
}
#endif
