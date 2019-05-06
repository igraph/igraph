
#ifndef CLIQUER_MISC_H
#define CLIQUER_MISC_H

#include "cliquerconf.h"

/*
 * We #define boolean instead of using a typedef because nauty.h uses it
 * also.  AFAIK, there is no way to check for an existing typedef, and
 * re-typedefing is illegal (even when using exactly the same datatype!).
 */
#ifndef boolean
#define boolean int
#endif


/*
 * The original cliquer source has some functions incorrectly marked as unused,
 * thus leave this undefined.
 */
#define UNUSED_FUNCTION


/*
 * Default inlining directive:  "inline"
 */
#ifndef INLINE
#define INLINE inline
#endif


#include <stdio.h>
#include <stdlib.h>

#ifndef ASSERT
#ifdef USING_R
#include <R.h>
#define ASSERT(expr) \
        if (!(expr)) { \
	        error("cliquer file %s: line %d: assertion failed: " \
			"(%s)\n",__FILE__,__LINE__,#expr); \
	}
#else
#define ASSERT(expr) \
        if (!(expr)) { \
		fprintf(stderr,"cliquer file %s: line %d: assertion failed: " \
			"(%s)\n",__FILE__,__LINE__,#expr); \
		abort(); \
	}
#endif
#endif /* !ASSERT */


#ifndef FALSE
#define FALSE (0)
#endif
#ifndef TRUE
#define TRUE (!FALSE)
#endif


#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif
#ifndef ABS
#define ABS(v)  (((v)<0)?(-(v)):(v))
#endif

#endif /* !CLIQUER_MISC_H */

