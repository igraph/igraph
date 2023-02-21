#include "stdio.h"
#include "signal.h"

#include "igraph_error.h"

#ifndef SIGIOT
#ifdef SIGABRT
#define SIGIOT SIGABRT
#endif
#endif

#ifdef KR_headers
void sig_die(s, kill) char *s; int kill;
#else
#include "stdlib.h"
#ifdef __cplusplus
extern "C" {
#endif
#ifdef __cplusplus
extern "C" {
#endif
 extern void f_exit(void);

void sig_die(const char *s, int kill)
#endif
{
	/* print error message, then clear buffers */
	fprintf(stderr, "%s\n", s);

	if(kill)
		{
		fflush(stderr);
		f_exit();
		fflush(stderr);
		/* now get a core */
#ifdef SIGIOT
		signal(SIGIOT, SIG_DFL);
#endif
		}
	else {
#ifdef NO_ONEXIT
		f_exit();
#endif
		}
	IGRAPH_FATAL("sig_die() called from f2c code");
	}
#ifdef __cplusplus
}
#endif
#ifdef __cplusplus
}
#endif
