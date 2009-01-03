/* dnrm2.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "config.h"
#include "arpack_internal.h"

doublereal igraphdnrm2_(n, dx, incx)
integer *n;
doublereal *dx;
integer *incx;
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal cutlo = 1.415686533102923e-146;
    static doublereal cuthi = 1.34078079299426e154;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal xmax;
    static integer next, i__, j, ix;
    static doublereal hitest, sum;

    /* Assigned format variables */
    static char *next_fmt;

    /* Parameter adjustments */
    --dx;

    /* Function Body */

/*     euclidean norm of the n-vector stored in dx() with storage */
/*     increment incx . */
/*     if    n .le. 0 return with result = 0. */
/*     if n .ge. 1 then incx must be .ge. 1 */

/*           c.l.lawson, 1978 jan 08 */
/*     modified to correct problem with negative increment, 8/21/90. */

/*     four phase method     using two built-in constants that are */
/*     hopefully applicable to all machines. */
/*         cutlo = maximum of  dsqrt(u/eps)  over all known machines. */
/*         cuthi = minimum of  dsqrt(v)      over all known machines. */
/*     where */
/*         eps = smallest no. such that eps + 1. .gt. 1. */
/*         u   = smallest positive no.   (underflow limit) */
/*         v   = largest  no.            (overflow  limit) */

/*     brief outline of algorithm.. */

/*     phase 1    scans zero components. */
/*     move to phase 2 when a component is nonzero and .le. cutlo */
/*     move to phase 3 when a component is .gt. cutlo */
/*     move to phase 4 when a component is .ge. cuthi/m */
/*     where m = n for x() real and m = 2*n for complex. */

/*     values for cutlo and cuthi.. */
/*     from the environmental parameters listed in the imsl converter */
/*     document the limiting values are as follows.. */
/*     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are */
/*                   univac and dec at 2**(-103) */
/*                   thus cutlo = 2**(-51) = 4.44089e-16 */
/*     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec. */
/*                   thus cuthi = 2**(63.5) = 1.30438e19 */
/*     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec. */
/*                   thus cutlo = 2**(-33.5) = 8.23181d-11 */
/*     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19 */
/*     data cutlo, cuthi / 8.232d-11,  1.304d19 / */
/*     data cutlo, cuthi / 4.441e-16,  1.304e19 / */
/*     data cutlo, cuthi / 8.232d-11,  1.304d19 / */
/* .... from Ed Anderson */

    if (*n > 0) {
	goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
    next_fmt = fmt_30;
    sum = zero;
    i__ = 1;
    if (*incx < 0) {
	i__ = (-(*n) + 1) * *incx + 1;
    }
    ix = 1;
/*                                                 begin main loop */
L20:
    switch ((int)next) {
	case 0: goto L30;
	case 1: goto L50;
	case 2: goto L70;
	case 3: goto L110;
    }
L30:
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
	goto L85;
    }
    next = 1;
    next_fmt = fmt_50;
    xmax = zero;

/*                        phase 1.  sum is zero */

L50:
    if (dx[i__] == zero) {
	goto L200;
    }
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
	goto L85;
    }

/*                                prepare for phase 2. */
    next = 2;
    next_fmt = fmt_70;
    goto L105;

/*                                prepare for phase 4. */

L100:
    next = 3;
    next_fmt = fmt_110;
    sum = sum / dx[i__] / dx[i__];
L105:
    xmax = (d__1 = dx[i__], abs(d__1));
    goto L115;

/*                   phase 2.  sum is small. */
/*                             scale to avoid destructive underflow. */

L70:
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
	goto L75;
    }

/*                     common code for phases 2 and 4. */
/*                     in phase 4 sum is large.  scale to avoid overflow. */

L110:
    if ((d__1 = dx[i__], abs(d__1)) <= xmax) {
	goto L115;
    }
/* Computing 2nd power */
    d__1 = xmax / dx[i__];
    sum = one + sum * (d__1 * d__1);
    xmax = (d__1 = dx[i__], abs(d__1));
    goto L200;

L115:
/* Computing 2nd power */
    d__1 = dx[i__] / xmax;
    sum += d__1 * d__1;
    goto L200;


/*                  prepare for phase 3. */

L75:
    sum = sum * xmax * xmax;


/*     for real or d.p. set hitest = cuthi/n */
/*     for complex      set hitest = cuthi/(2*n) */

L85:
    hitest = cuthi / (real) (*n);

/*                   phase 3.  sum is mid-range.  no scaling. */

    i__1 = *n;
    for (j = ix; j <= i__1; ++j) {
	if ((d__1 = dx[i__], abs(d__1)) >= hitest) {
	    goto L100;
	}
/* Computing 2nd power */
	d__1 = dx[i__];
	sum += d__1 * d__1;
	i__ += *incx;
/* L95: */
    }
    ret_val = sqrt(sum);
    goto L300;

L200:
    ++ix;
    i__ += *incx;
    if (ix <= *n) {
	goto L20;
    }

/*              end of main loop. */

/*              compute square root and adjust for scaling. */

    ret_val = xmax * sqrt(sum);
L300:
    return ret_val;
} /* igraphdnrm2_ */

