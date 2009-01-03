/* igraphdasum.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "config.h"
#include "arpack_internal.h"

doublereal igraphdasum_(n, dx, incx)
integer *n;
doublereal *dx;
integer *incx;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static integer i__, m;
    static doublereal dtemp;
    static integer ix, mp1;


/*     takes the sum of the absolute values. */
/*     uses unrolled loops for increment equal to one. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified to correct problem with negative increment, 8/21/90. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += (d__1 = dx[ix], abs(d__1));
	ix += *incx;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for increment equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += (d__1 = dx[i__], abs(d__1));
/* L30: */
    }
    if (*n < 6) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 6) {
	dtemp = dtemp + (d__1 = dx[i__], abs(d__1)) + (d__2 = dx[i__ + 1], 
		abs(d__2)) + (d__3 = dx[i__ + 2], abs(d__3)) + (d__4 = dx[i__ 
		+ 3], abs(d__4)) + (d__5 = dx[i__ + 4], abs(d__5)) + (d__6 = 
		dx[i__ + 5], abs(d__6));
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* igraphdasum_ */

