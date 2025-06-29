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

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__10 = 10;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__11 = 11;

/* > \brief \b DLASQ2 computes all the eigenvalues of the symmetric positive definite tridiagonal matrix assoc
iated with the qd Array Z to high relative accuracy. Used by sbdsqr and sstegr.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLASQ2 + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq2.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq2.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq2.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DLASQ2( N, Z, INFO )   

         INTEGER            INFO, N   
         DOUBLE PRECISION   Z( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DLASQ2 computes all the eigenvalues of the symmetric positive   
   > definite tridiagonal matrix associated with the qd array Z to high   
   > relative accuracy are computed to high relative accuracy, in the   
   > absence of denormalization, underflow and overflow.   
   >   
   > To see the relation of Z to the tridiagonal matrix, let L be a   
   > unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and   
   > let U be an upper bidiagonal matrix with 1's above and diagonal   
   > Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the   
   > symmetric tridiagonal to which it is similar.   
   >   
   > Note : DLASQ2 defines a logical variable, IEEE, which is true   
   > on machines which follow ieee-754 floating-point standard in their   
   > handling of infinities and NaNs, and false otherwise. This variable   
   > is passed to DLASQ3.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >        The number of rows and columns in the matrix. N >= 0.   
   > \endverbatim   
   >   
   > \param[in,out] Z   
   > \verbatim   
   >          Z is DOUBLE PRECISION array, dimension ( 4*N )   
   >        On entry Z holds the qd array. On exit, entries 1 to N hold   
   >        the eigenvalues in decreasing order, Z( 2*N+1 ) holds the   
   >        trace, and Z( 2*N+2 ) holds the sum of the eigenvalues. If   
   >        N > 2, then Z( 2*N+3 ) holds the iteration count, Z( 2*N+4 )   
   >        holds NDIVS/NIN^2, and Z( 2*N+5 ) holds the percentage of   
   >        shifts that failed.   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
   >        = 0: successful exit   
   >        < 0: if the i-th argument is a scalar and had an illegal   
   >             value, then INFO = -i, if the i-th argument is an   
   >             array and the j-entry had an illegal value, then   
   >             INFO = -(i*100+j)   
   >        > 0: the algorithm failed   
   >              = 1, a split was marked by a positive value in E   
   >              = 2, current block of Z not diagonalized after 100*N   
   >                   iterations (in inner while loop).  On exit Z holds   
   >                   a qd array with the same eigenvalues as the given Z.   
   >              = 3, termination criterion of outer while loop not met   
   >                   (program created more than N unreduced blocks)   
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
   >  Local Variables: I0:N0 defines a current unreduced segment of Z.   
   >  The shifts are accumulated in SIGMA. Iteration count is in ITER.   
   >  Ping-pong is controlled by PP (alternates between 0 and 1).   
   > \endverbatim   
   >   
    =====================================================================   
   Subroutine */ int igraphdlasq2_(integer *n, doublereal *z__, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal d__, e, g;
    integer k;
    doublereal s, t;
    integer i0, i1, i4, n0, n1;
    doublereal dn;
    integer pp;
    doublereal dn1, dn2, dee, eps, tau, tol;
    integer ipn4;
    doublereal tol2;
    logical ieee;
    integer nbig;
    doublereal dmin__, emin, emax;
    integer kmin, ndiv, iter;
    doublereal qmin, temp, qmax, zmax;
    integer splt;
    doublereal dmin1, dmin2;
    integer nfail;
    doublereal desig, trace, sigma;
    integer iinfo;
    doublereal tempe, tempq;
    integer ttype;
    extern /* Subroutine */ int igraphdlasq3_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, logical *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal igraphdlamch_(char *);
    doublereal deemin;
    integer iwhila, iwhilb;
    doublereal oldemn, safmin;
    extern /* Subroutine */ int igraphxerbla_(char *, integer *, ftnlen);
    extern integer igraphilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int igraphdlasrt_(char *, integer *, doublereal *, 
	    integer *);


/*  -- LAPACK computational routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       Test the input arguments.   
       (in case DLASQ2 is not called by DLASQ1)   

       Parameter adjustments */
    --z__;

    /* Function Body */
    *info = 0;
    eps = igraphdlamch_("Precision");
    safmin = igraphdlamch_("Safe minimum");
    tol = eps * 100.;
/* Computing 2nd power */
    d__1 = tol;
    tol2 = d__1 * d__1;

    if (*n < 0) {
	*info = -1;
	igraphxerbla_("DLASQ2", &c__1, (ftnlen)6);
	return 0;
    } else if (*n == 0) {
	return 0;
    } else if (*n == 1) {

/*        1-by-1 case. */

	if (z__[1] < 0.) {
	    *info = -201;
	    igraphxerbla_("DLASQ2", &c__2, (ftnlen)6);
	}
	return 0;
    } else if (*n == 2) {

/*        2-by-2 case. */

	if (z__[2] < 0. || z__[3] < 0.) {
	    *info = -2;
	    igraphxerbla_("DLASQ2", &c__2, (ftnlen)6);
	    return 0;
	} else if (z__[3] > z__[1]) {
	    d__ = z__[3];
	    z__[3] = z__[1];
	    z__[1] = d__;
	}
	z__[5] = z__[1] + z__[2] + z__[3];
	if (z__[2] > z__[3] * tol2) {
	    t = (z__[1] - z__[3] + z__[2]) * .5;
	    s = z__[3] * (z__[2] / t);
	    if (s <= t) {
		s = z__[3] * (z__[2] / (t * (sqrt(s / t + 1.) + 1.)));
	    } else {
		s = z__[3] * (z__[2] / (t + sqrt(t) * sqrt(t + s)));
	    }
	    t = z__[1] + (s + z__[2]);
	    z__[3] *= z__[1] / t;
	    z__[1] = t;
	}
	z__[2] = z__[3];
	z__[6] = z__[2] + z__[1];
	return 0;
    }

/*     Check for negative data and compute sums of q's and e's. */

    z__[*n * 2] = 0.;
    emin = z__[2];
    qmax = 0.;
    zmax = 0.;
    d__ = 0.;
    e = 0.;

    i__1 = *n - 1 << 1;
    for (k = 1; k <= i__1; k += 2) {
	if (z__[k] < 0.) {
	    *info = -(k + 200);
	    igraphxerbla_("DLASQ2", &c__2, (ftnlen)6);
	    return 0;
	} else if (z__[k + 1] < 0.) {
	    *info = -(k + 201);
	    igraphxerbla_("DLASQ2", &c__2, (ftnlen)6);
	    return 0;
	}
	d__ += z__[k];
	e += z__[k + 1];
/* Computing MAX */
	d__1 = qmax, d__2 = z__[k];
	qmax = max(d__1,d__2);
/* Computing MIN */
	d__1 = emin, d__2 = z__[k + 1];
	emin = min(d__1,d__2);
/* Computing MAX */
	d__1 = max(qmax,zmax), d__2 = z__[k + 1];
	zmax = max(d__1,d__2);
/* L10: */
    }
    if (z__[(*n << 1) - 1] < 0.) {
	*info = -((*n << 1) + 199);
	igraphxerbla_("DLASQ2", &c__2, (ftnlen)6);
	return 0;
    }
    d__ += z__[(*n << 1) - 1];
/* Computing MAX */
    d__1 = qmax, d__2 = z__[(*n << 1) - 1];
    qmax = max(d__1,d__2);
    zmax = max(qmax,zmax);

/*     Check for diagonality. */

    if (e == 0.) {
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    z__[k] = z__[(k << 1) - 1];
/* L20: */
	}
	igraphdlasrt_("D", n, &z__[1], &iinfo);
	z__[(*n << 1) - 1] = d__;
	return 0;
    }

    trace = d__ + e;

/*     Check for zero data. */

    if (trace == 0.) {
	z__[(*n << 1) - 1] = 0.;
	return 0;
    }

/*     Check whether the machine is IEEE conformable. */

    ieee = igraphilaenv_(&c__10, "DLASQ2", "N", &c__1, &c__2, &c__3, &c__4, (ftnlen)
	    6, (ftnlen)1) == 1 && igraphilaenv_(&c__11, "DLASQ2", "N", &c__1, &c__2,
	     &c__3, &c__4, (ftnlen)6, (ftnlen)1) == 1;

/*     Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...). */

    for (k = *n << 1; k >= 2; k += -2) {
	z__[k * 2] = 0.;
	z__[(k << 1) - 1] = z__[k];
	z__[(k << 1) - 2] = 0.;
	z__[(k << 1) - 3] = z__[k - 1];
/* L30: */
    }

    i0 = 1;
    n0 = *n;

/*     Reverse the qd-array, if warranted. */

    if (z__[(i0 << 2) - 3] * 1.5 < z__[(n0 << 2) - 3]) {
	ipn4 = i0 + n0 << 2;
	i__1 = i0 + n0 - 1 << 1;
	for (i4 = i0 << 2; i4 <= i__1; i4 += 4) {
	    temp = z__[i4 - 3];
	    z__[i4 - 3] = z__[ipn4 - i4 - 3];
	    z__[ipn4 - i4 - 3] = temp;
	    temp = z__[i4 - 1];
	    z__[i4 - 1] = z__[ipn4 - i4 - 5];
	    z__[ipn4 - i4 - 5] = temp;
/* L40: */
	}
    }

/*     Initial split checking via dqd and Li's test. */

    pp = 0;

    for (k = 1; k <= 2; ++k) {

	d__ = z__[(n0 << 2) + pp - 3];
	i__1 = (i0 << 2) + pp;
	for (i4 = (n0 - 1 << 2) + pp; i4 >= i__1; i4 += -4) {
	    if (z__[i4 - 1] <= tol2 * d__) {
		z__[i4 - 1] = -0.;
		d__ = z__[i4 - 3];
	    } else {
		d__ = z__[i4 - 3] * (d__ / (d__ + z__[i4 - 1]));
	    }
/* L50: */
	}

/*        dqd maps Z to ZZ plus Li's test. */

	emin = z__[(i0 << 2) + pp + 1];
	d__ = z__[(i0 << 2) + pp - 3];
	i__1 = (n0 - 1 << 2) + pp;
	for (i4 = (i0 << 2) + pp; i4 <= i__1; i4 += 4) {
	    z__[i4 - (pp << 1) - 2] = d__ + z__[i4 - 1];
	    if (z__[i4 - 1] <= tol2 * d__) {
		z__[i4 - 1] = -0.;
		z__[i4 - (pp << 1) - 2] = d__;
		z__[i4 - (pp << 1)] = 0.;
		d__ = z__[i4 + 1];
	    } else if (safmin * z__[i4 + 1] < z__[i4 - (pp << 1) - 2] && 
		    safmin * z__[i4 - (pp << 1) - 2] < z__[i4 + 1]) {
		temp = z__[i4 + 1] / z__[i4 - (pp << 1) - 2];
		z__[i4 - (pp << 1)] = z__[i4 - 1] * temp;
		d__ *= temp;
	    } else {
		z__[i4 - (pp << 1)] = z__[i4 + 1] * (z__[i4 - 1] / z__[i4 - (
			pp << 1) - 2]);
		d__ = z__[i4 + 1] * (d__ / z__[i4 - (pp << 1) - 2]);
	    }
/* Computing MIN */
	    d__1 = emin, d__2 = z__[i4 - (pp << 1)];
	    emin = min(d__1,d__2);
/* L60: */
	}
	z__[(n0 << 2) - pp - 2] = d__;

/*        Now find qmax. */

	qmax = z__[(i0 << 2) - pp - 2];
	i__1 = (n0 << 2) - pp - 2;
	for (i4 = (i0 << 2) - pp + 2; i4 <= i__1; i4 += 4) {
/* Computing MAX */
	    d__1 = qmax, d__2 = z__[i4];
	    qmax = max(d__1,d__2);
/* L70: */
	}

/*        Prepare for the next iteration on K. */

	pp = 1 - pp;
/* L80: */
    }

/*     Initialise variables to pass to DLASQ3. */

    ttype = 0;
    dmin1 = 0.;
    dmin2 = 0.;
    dn = 0.;
    dn1 = 0.;
    dn2 = 0.;
    g = 0.;
    tau = 0.;

    iter = 2;
    nfail = 0;
    ndiv = n0 - i0 << 1;

    i__1 = *n + 1;
    for (iwhila = 1; iwhila <= i__1; ++iwhila) {
	if (n0 < 1) {
	    goto L170;
	}

/*        While array unfinished do   

          E(N0) holds the value of SIGMA when submatrix in I0:N0   
          splits from the rest of the array, but is negated. */

	desig = 0.;
	if (n0 == *n) {
	    sigma = 0.;
	} else {
	    sigma = -z__[(n0 << 2) - 1];
	}
	if (sigma < 0.) {
	    *info = 1;
	    return 0;
	}

/*        Find last unreduced submatrix's top index I0, find QMAX and   
          EMIN. Find Gershgorin-type bound if Q's much greater than E's. */

	emax = 0.;
	if (n0 > i0) {
	    emin = (d__1 = z__[(n0 << 2) - 5], abs(d__1));
	} else {
	    emin = 0.;
	}
	qmin = z__[(n0 << 2) - 3];
	qmax = qmin;
	for (i4 = n0 << 2; i4 >= 8; i4 += -4) {
	    if (z__[i4 - 5] <= 0.) {
		goto L100;
	    }
	    if (qmin >= emax * 4.) {
/* Computing MIN */
		d__1 = qmin, d__2 = z__[i4 - 3];
		qmin = min(d__1,d__2);
/* Computing MAX */
		d__1 = emax, d__2 = z__[i4 - 5];
		emax = max(d__1,d__2);
	    }
/* Computing MAX */
	    d__1 = qmax, d__2 = z__[i4 - 7] + z__[i4 - 5];
	    qmax = max(d__1,d__2);
/* Computing MIN */
	    d__1 = emin, d__2 = z__[i4 - 5];
	    emin = min(d__1,d__2);
/* L90: */
	}
	i4 = 4;

L100:
	i0 = i4 / 4;
	pp = 0;

	if (n0 - i0 > 1) {
	    dee = z__[(i0 << 2) - 3];
	    deemin = dee;
	    kmin = i0;
	    i__2 = (n0 << 2) - 3;
	    for (i4 = (i0 << 2) + 1; i4 <= i__2; i4 += 4) {
		dee = z__[i4] * (dee / (dee + z__[i4 - 2]));
		if (dee <= deemin) {
		    deemin = dee;
		    kmin = (i4 + 3) / 4;
		}
/* L110: */
	    }
	    if (kmin - i0 << 1 < n0 - kmin && deemin <= z__[(n0 << 2) - 3] * 
		    .5) {
		ipn4 = i0 + n0 << 2;
		pp = 2;
		i__2 = i0 + n0 - 1 << 1;
		for (i4 = i0 << 2; i4 <= i__2; i4 += 4) {
		    temp = z__[i4 - 3];
		    z__[i4 - 3] = z__[ipn4 - i4 - 3];
		    z__[ipn4 - i4 - 3] = temp;
		    temp = z__[i4 - 2];
		    z__[i4 - 2] = z__[ipn4 - i4 - 2];
		    z__[ipn4 - i4 - 2] = temp;
		    temp = z__[i4 - 1];
		    z__[i4 - 1] = z__[ipn4 - i4 - 5];
		    z__[ipn4 - i4 - 5] = temp;
		    temp = z__[i4];
		    z__[i4] = z__[ipn4 - i4 - 4];
		    z__[ipn4 - i4 - 4] = temp;
/* L120: */
		}
	    }
	}

/*        Put -(initial shift) into DMIN.   

   Computing MAX */
	d__1 = 0., d__2 = qmin - sqrt(qmin) * 2. * sqrt(emax);
	dmin__ = -max(d__1,d__2);

/*        Now I0:N0 is unreduced.   
          PP = 0 for ping, PP = 1 for pong.   
          PP = 2 indicates that flipping was applied to the Z array and   
                 and that the tests for deflation upon entry in DLASQ3   
                 should not be performed. */

	nbig = (n0 - i0 + 1) * 100;
	i__2 = nbig;
	for (iwhilb = 1; iwhilb <= i__2; ++iwhilb) {
	    if (i0 > n0) {
		goto L150;
	    }

/*           While submatrix unfinished take a good dqds step. */

	    igraphdlasq3_(&i0, &n0, &z__[1], &pp, &dmin__, &sigma, &desig, &qmax, &
		    nfail, &iter, &ndiv, &ieee, &ttype, &dmin1, &dmin2, &dn, &
		    dn1, &dn2, &g, &tau);

	    pp = 1 - pp;

/*           When EMIN is very small check for splits. */

	    if (pp == 0 && n0 - i0 >= 3) {
		if (z__[n0 * 4] <= tol2 * qmax || z__[(n0 << 2) - 1] <= tol2 *
			 sigma) {
		    splt = i0 - 1;
		    qmax = z__[(i0 << 2) - 3];
		    emin = z__[(i0 << 2) - 1];
		    oldemn = z__[i0 * 4];
		    i__3 = n0 - 3 << 2;
		    for (i4 = i0 << 2; i4 <= i__3; i4 += 4) {
			if (z__[i4] <= tol2 * z__[i4 - 3] || z__[i4 - 1] <= 
				tol2 * sigma) {
			    z__[i4 - 1] = -sigma;
			    splt = i4 / 4;
			    qmax = 0.;
			    emin = z__[i4 + 3];
			    oldemn = z__[i4 + 4];
			} else {
/* Computing MAX */
			    d__1 = qmax, d__2 = z__[i4 + 1];
			    qmax = max(d__1,d__2);
/* Computing MIN */
			    d__1 = emin, d__2 = z__[i4 - 1];
			    emin = min(d__1,d__2);
/* Computing MIN */
			    d__1 = oldemn, d__2 = z__[i4];
			    oldemn = min(d__1,d__2);
			}
/* L130: */
		    }
		    z__[(n0 << 2) - 1] = emin;
		    z__[n0 * 4] = oldemn;
		    i0 = splt + 1;
		}
	    }

/* L140: */
	}

	*info = 2;

/*        Maximum number of iterations exceeded, restore the shift   
          SIGMA and place the new d's and e's in a qd array.   
          This might need to be done for several blocks */

	i1 = i0;
	n1 = n0;
L145:
	tempq = z__[(i0 << 2) - 3];
	z__[(i0 << 2) - 3] += sigma;
	i__2 = n0;
	for (k = i0 + 1; k <= i__2; ++k) {
	    tempe = z__[(k << 2) - 5];
	    z__[(k << 2) - 5] *= tempq / z__[(k << 2) - 7];
	    tempq = z__[(k << 2) - 3];
	    z__[(k << 2) - 3] = z__[(k << 2) - 3] + sigma + tempe - z__[(k << 
		    2) - 5];
	}

/*        Prepare to do this on the previous block if there is one */

	if (i1 > 1) {
	    n1 = i1 - 1;
	    while(i1 >= 2 && z__[(i1 << 2) - 5] >= 0.) {
		--i1;
	    }
	    sigma = -z__[(n1 << 2) - 1];
	    goto L145;
	}
	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {
	    z__[(k << 1) - 1] = z__[(k << 2) - 3];

/*        Only the block 1..N0 is unfinished.  The rest of the e's   
          must be essentially zero, although sometimes other data   
          has been stored in them. */

	    if (k < n0) {
		z__[k * 2] = z__[(k << 2) - 1];
	    } else {
		z__[k * 2] = 0.;
	    }
	}
	return 0;

/*        end IWHILB */

L150:

/* L160: */
	;
    }

    *info = 3;
    return 0;

/*     end IWHILA */

L170:

/*     Move q's to the front. */

    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	z__[k] = z__[(k << 2) - 3];
/* L180: */
    }

/*     Sort and compute sum of eigenvalues. */

    igraphdlasrt_("D", n, &z__[1], &iinfo);

    e = 0.;
    for (k = *n; k >= 1; --k) {
	e += z__[k];
/* L190: */
    }

/*     Store trace, sum(eigenvalues) and information on performance. */

    z__[(*n << 1) + 1] = trace;
    z__[(*n << 1) + 2] = e;
    z__[(*n << 1) + 3] = (doublereal) iter;
/* Computing 2nd power */
    i__1 = *n;
    z__[(*n << 1) + 4] = (doublereal) ndiv / (doublereal) (i__1 * i__1);
    z__[(*n << 1) + 5] = nfail * 100. / (doublereal) iter;
    return 0;

/*     End of DLASQ2 */

} /* igraphdlasq2_ */

