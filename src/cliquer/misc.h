
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
 * Default value for UNUSED_FUNCTION:  use "__attribute__((unused))" for
 * GCC versions that support it, otherwise leave blank.
 */
#ifndef UNUSED_FUNCTION
# if     __GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ > 4)
#  define UNUSED_FUNCTION __attribute__((unused))
# else
#  define UNUSED_FUNCTION
# endif
#endif  /* !UNUSED_FUNCTION */


/*
 * Default inlining directive:  "inline"
 */
#ifndef INLINE
#define INLINE inline
#endif


#include <stdio.h>
#include <stdlib.h>

#ifndef ASSERT
#define ASSERT(expr) \
        if (!(expr)) { \
		fprintf(stderr,"cliquer file %s: line %d: assertion failed: " \
			"(%s)\n",__FILE__,__LINE__,#expr); \
		abort(); \
	}
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

