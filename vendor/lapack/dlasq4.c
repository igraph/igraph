/*  -- translated by f2c (version 20240504).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* > \brief \b DLASQ4 computes an approximation to the smallest eigenvalue using values of d from the previous
 transform. Used by sbdsqr.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLASQ4 + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq4.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq4.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq4.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN,   
                            DN1, DN2, TAU, TTYPE, G )   

         INTEGER            I0, N0, N0IN, PP, TTYPE   
         DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU   
         DOUBLE PRECISION   Z( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DLASQ4 computes an approximation TAU to the smallest eigenvalue   
   > using values of d from the previous transform.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] I0   
   > \verbatim   
   >          I0 is INTEGER   
   >        First index.   
   > \endverbatim   
   >   
   > \param[in] N0   
   > \verbatim   
   >          N0 is INTEGER   
   >        Last index.   
   > \endverbatim   
   >   
   > \param[in] Z   
   > \verbatim   
   >          Z is DOUBLE PRECISION array, dimension ( 4*N )   
   >        Z holds the qd array.   
   > \endverbatim   
   >   
   > \param[in] PP   
   > \verbatim   
   >          PP is INTEGER   
   >        PP=0 for ping, PP=1 for pong.   
   > \endverbatim   
   >   
   > \param[in] N0IN   
   > \verbatim   
   >          N0IN is INTEGER   
   >        The value of N0 at start of EIGTEST.   
   > \endverbatim   
   >   
   > \param[in] DMIN   
   > \verbatim   
   >          DMIN is DOUBLE PRECISION   
   >        Minimum value of d.   
   > \endverbatim   
   >   
   > \param[in] DMIN1   
   > \verbatim   
   >          DMIN1 is DOUBLE PRECISION   
   >        Minimum value of d, excluding D( N0 ).   
   > \endverbatim   
   >   
   > \param[in] DMIN2   
   > \verbatim   
   >          DMIN2 is DOUBLE PRECISION   
   >        Minimum value of d, excluding D( N0 ) and D( N0-1 ).   
   > \endverbatim   
   >   
   > \param[in] DN   
   > \verbatim   
   >          DN is DOUBLE PRECISION   
   >        d(N)   
   > \endverbatim   
   >   
   > \param[in] DN1   
   > \verbatim   
   >          DN1 is DOUBLE PRECISION   
   >        d(N-1)   
   > \endverbatim   
   >   
   > \param[in] DN2   
   > \verbatim   
   >          DN2 is DOUBLE PRECISION   
   >        d(N-2)   
   > \endverbatim   
   >   
   > \param[out] TAU   
   > \verbatim   
   >          TAU is DOUBLE PRECISION   
   >        This is the shift.   
   > \endverbatim   
   >   
   > \param[out] TTYPE   
   > \verbatim   
   >          TTYPE is INTEGER   
   >        Shift type.   
   > \endverbatim   
   >   
   > \param[in,out] G   
   > \verbatim   
   >          G is REAL   
   >        G is passed as an argument in order to save its value between   
   >        calls to DLASQ4.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup auxOTHERcomputational   

   > \par Further Details:   
    =====================   
   >   
   > \verbatim   
   >   
   >  CNST1 = 9/16   
   > \endverbatim   
   >   
    =====================================================================   
   Subroutine */ int igraphdlasq4_(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, integer *n0in, doublereal *dmin__, doublereal *dmin1, 
	doublereal *dmin2, doublereal *dn, doublereal *dn1, doublereal *dn2, 
	doublereal *tau, integer *ttype, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal s, a2, b1, b2;
    integer i4, nn, np;
    doublereal gam, gap1, gap2;


/*  -- LAPACK computational routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       A negative DMIN forces the shift to take that absolute value   
       TTYPE records the type of shift.   

       Parameter adjustments */
    --z__;

    /* Function Body */
    if (*dmin__ <= 0.) {
	*tau = -(*dmin__);
	*ttype = -1;
	return 0;
    }

    nn = (*n0 << 2) + *pp;
    if (*n0in == *n0) {

/*        No eigenvalues deflated. */

	if (*dmin__ == *dn || *dmin__ == *dn1) {

	    b1 = sqrt(z__[nn - 3]) * sqrt(z__[nn - 5]);
	    b2 = sqrt(z__[nn - 7]) * sqrt(z__[nn - 9]);
	    a2 = z__[nn - 7] + z__[nn - 5];

/*           Cases 2 and 3. */

	    if (*dmin__ == *dn && *dmin1 == *dn1) {
		gap2 = *dmin2 - a2 - *dmin2 * .25;
		if (gap2 > 0. && gap2 > b2) {
		    gap1 = a2 - *dn - b2 / gap2 * b2;
		} else {
		    gap1 = a2 - *dn - (b1 + b2);
		}
		if (gap1 > 0. && gap1 > b1) {
/* Computing MAX */
		    d__1 = *dn - b1 / gap1 * b1, d__2 = *dmin__ * .5;
		    s = max(d__1,d__2);
		    *ttype = -2;
		} else {
		    s = 0.;
		    if (*dn > b1) {
			s = *dn - b1;
		    }
		    if (a2 > b1 + b2) {
/* Computing MIN */
			d__1 = s, d__2 = a2 - (b1 + b2);
			s = min(d__1,d__2);
		    }
/* Computing MAX */
		    d__1 = s, d__2 = *dmin__ * .333;
		    s = max(d__1,d__2);
		    *ttype = -3;
		}
	    } else {

/*              Case 4. */

		*ttype = -4;
		s = *dmin__ * .25;
		if (*dmin__ == *dn) {
		    gam = *dn;
		    a2 = 0.;
		    if (z__[nn - 5] > z__[nn - 7]) {
			return 0;
		    }
		    b2 = z__[nn - 5] / z__[nn - 7];
		    np = nn - 9;
		} else {
		    np = nn - (*pp << 1);
		    b2 = z__[np - 2];
		    gam = *dn1;
		    if (z__[np - 4] > z__[np - 2]) {
			return 0;
		    }
		    a2 = z__[np - 4] / z__[np - 2];
		    if (z__[nn - 9] > z__[nn - 11]) {
			return 0;
		    }
		    b2 = z__[nn - 9] / z__[nn - 11];
		    np = nn - 13;
		}

/*              Approximate contribution to norm squared from I < NN-1. */

		a2 += b2;
		i__1 = (*i0 << 2) - 1 + *pp;
		for (i4 = np; i4 >= i__1; i4 += -4) {
		    if (b2 == 0.) {
			goto L20;
		    }
		    b1 = b2;
		    if (z__[i4] > z__[i4 - 2]) {
			return 0;
		    }
		    b2 *= z__[i4] / z__[i4 - 2];
		    a2 += b2;
		    if (max(b2,b1) * 100. < a2 || .563 < a2) {
			goto L20;
		    }
/* L10: */
		}
L20:
		a2 *= 1.05;

/*              Rayleigh quotient residual bound. */

		if (a2 < .563) {
		    s = gam * (1. - sqrt(a2)) / (a2 + 1.);
		}
	    }
	} else if (*dmin__ == *dn2) {

/*           Case 5. */

	    *ttype = -5;
	    s = *dmin__ * .25;

/*           Compute contribution to norm squared from I > NN-2. */

	    np = nn - (*pp << 1);
	    b1 = z__[np - 2];
	    b2 = z__[np - 6];
	    gam = *dn2;
	    if (z__[np - 8] > b2 || z__[np - 4] > b1) {
		return 0;
	    }
	    a2 = z__[np - 8] / b2 * (z__[np - 4] / b1 + 1.);

/*           Approximate contribution to norm squared from I < NN-2. */

	    if (*n0 - *i0 > 2) {
		b2 = z__[nn - 13] / z__[nn - 15];
		a2 += b2;
		i__1 = (*i0 << 2) - 1 + *pp;
		for (i4 = nn - 17; i4 >= i__1; i4 += -4) {
		    if (b2 == 0.) {
			goto L40;
		    }
		    b1 = b2;
		    if (z__[i4] > z__[i4 - 2]) {
			return 0;
		    }
		    b2 *= z__[i4] / z__[i4 - 2];
		    a2 += b2;
		    if (max(b2,b1) * 100. < a2 || .563 < a2) {
			goto L40;
		    }
/* L30: */
		}
L40:
		a2 *= 1.05;
	    }

	    if (a2 < .563) {
		s = gam * (1. - sqrt(a2)) / (a2 + 1.);
	    }
	} else {

/*           Case 6, no information to guide us. */

	    if (*ttype == -6) {
		*g += (1. - *g) * .333;
	    } else if (*ttype == -18) {
		*g = .083250000000000005;
	    } else {
		*g = .25;
	    }
	    s = *g * *dmin__;
	    *ttype = -6;
	}

    } else if (*n0in == *n0 + 1) {

/*        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN. */

	if (*dmin1 == *dn1 && *dmin2 == *dn2) {

/*           Cases 7 and 8. */

	    *ttype = -7;
	    s = *dmin1 * .333;
	    if (z__[nn - 5] > z__[nn - 7]) {
		return 0;
	    }
	    b1 = z__[nn - 5] / z__[nn - 7];
	    b2 = b1;
	    if (b2 == 0.) {
		goto L60;
	    }
	    i__1 = (*i0 << 2) - 1 + *pp;
	    for (i4 = (*n0 << 2) - 9 + *pp; i4 >= i__1; i4 += -4) {
		a2 = b1;
		if (z__[i4] > z__[i4 - 2]) {
		    return 0;
		}
		b1 *= z__[i4] / z__[i4 - 2];
		b2 += b1;
		if (max(b1,a2) * 100. < b2) {
		    goto L60;
		}
/* L50: */
	    }
L60:
	    b2 = sqrt(b2 * 1.05);
/* Computing 2nd power */
	    d__1 = b2;
	    a2 = *dmin1 / (d__1 * d__1 + 1.);
	    gap2 = *dmin2 * .5 - a2;
	    if (gap2 > 0. && gap2 > b2 * a2) {
/* Computing MAX */
		d__1 = s, d__2 = a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2);
		s = max(d__1,d__2);
	    } else {
/* Computing MAX */
		d__1 = s, d__2 = a2 * (1. - b2 * 1.01);
		s = max(d__1,d__2);
		*ttype = -8;
	    }
	} else {

/*           Case 9. */

	    s = *dmin1 * .25;
	    if (*dmin1 == *dn1) {
		s = *dmin1 * .5;
	    }
	    *ttype = -9;
	}

    } else if (*n0in == *n0 + 2) {

/*        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.   

          Cases 10 and 11. */

	if (*dmin2 == *dn2 && z__[nn - 5] * 2. < z__[nn - 7]) {
	    *ttype = -10;
	    s = *dmin2 * .333;
	    if (z__[nn - 5] > z__[nn - 7]) {
		return 0;
	    }
	    b1 = z__[nn - 5] / z__[nn - 7];
	    b2 = b1;
	    if (b2 == 0.) {
		goto L80;
	    }
	    i__1 = (*i0 << 2) - 1 + *pp;
	    for (i4 = (*n0 << 2) - 9 + *pp; i4 >= i__1; i4 += -4) {
		if (z__[i4] > z__[i4 - 2]) {
		    return 0;
		}
		b1 *= z__[i4] / z__[i4 - 2];
		b2 += b1;
		if (b1 * 100. < b2) {
		    goto L80;
		}
/* L70: */
	    }
L80:
	    b2 = sqrt(b2 * 1.05);
/* Computing 2nd power */
	    d__1 = b2;
	    a2 = *dmin2 / (d__1 * d__1 + 1.);
	    gap2 = z__[nn - 7] + z__[nn - 9] - sqrt(z__[nn - 11]) * sqrt(z__[
		    nn - 9]) - a2;
	    if (gap2 > 0. && gap2 > b2 * a2) {
/* Computing MAX */
		d__1 = s, d__2 = a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2);
		s = max(d__1,d__2);
	    } else {
/* Computing MAX */
		d__1 = s, d__2 = a2 * (1. - b2 * 1.01);
		s = max(d__1,d__2);
	    }
	} else {
	    s = *dmin2 * .25;
	    *ttype = -11;
	}
    } else if (*n0in > *n0 + 2) {

/*        Case 12, more than two eigenvalues deflated. No information. */

	s = 0.;
	*ttype = -12;
    }

    *tau = s;
    return 0;

/*     End of DLASQ4 */

} /* igraphdlasq4_ */

