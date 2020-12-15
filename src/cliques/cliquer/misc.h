
#ifndef CLIQUER_MISC_H
#define CLIQUER_MISC_H

#include "igraph_error.h"
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
#define ASSERT IGRAPH_ASSERT
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

