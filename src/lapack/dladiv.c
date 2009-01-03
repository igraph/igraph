/* igraphdladiv.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "config.h"
#include "arpack_internal.h"

/* Subroutine */ int igraphdladiv_(a, b, c__, d__, p, q)
doublereal *a, *b, *c__, *d__, *p, *q;
{
    static doublereal e, f;


/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLADIV performs complex division in  real arithmetic */

/*                        a + i*b */
/*             p + i*q = --------- */
/*                        c + i*d */

/*  The algorithm is due to Robert L. Smith and can be found */
/*  in D. Knuth, The art of Computer Programming, Vol.2, p.195 */

/*  Arguments */
/*  ========= */

/*  A       (input) DOUBLE PRECISION */
/*  B       (input) DOUBLE PRECISION */
/*  C       (input) DOUBLE PRECISION */
/*  D       (input) DOUBLE PRECISION */
/*          The scalars a, b, c, and d in the above expression. */

/*  P       (output) DOUBLE PRECISION */
/*  Q       (output) DOUBLE PRECISION */
/*          The scalars p and q in the above expression. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    if (abs(*d__) < abs(*c__)) {
	e = *d__ / *c__;
	f = *c__ + *d__ * e;
	*p = (*a + *b * e) / f;
	*q = (*b - *a * e) / f;
    } else {
	e = *c__ / *d__;
	f = *d__ + *c__ * e;
	*p = (*b + *a * e) / f;
	*q = (-(*a) + *b * e) / f;
    }

    return 0;

/*     End of DLADIV */

} /* igraphdladiv_ */

