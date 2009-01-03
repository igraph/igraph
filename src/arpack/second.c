/* second.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "config.h"
#include "arpack_internal.h"

/* Subroutine */ int igraphsecond_(t)
real *t;
{
    extern doublereal igraphetime_();
    static real t1, tarray[2];



/*  -- LAPACK auxiliary routine (preliminary version) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     July 26, 1991 */

/*  Purpose */
/*  ======= */

/*  SECOND returns the user time for a process in seconds. */
/*  This version gets the time from the system function ETIME. */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    t1 = igraphetime_(tarray);
    *t = tarray[0];
    return 0;

/*     End of SECOND */

} /* igraphsecond_ */

