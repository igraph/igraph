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

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b33 = 1.;

/* > \brief \b DSTERF   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DSTERF + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsterf.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsterf.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsterf.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DSTERF( N, D, E, INFO )   

         INTEGER            INFO, N   
         DOUBLE PRECISION   D( * ), E( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DSTERF computes all eigenvalues of a symmetric tridiagonal matrix   
   > using the Pal-Walker-Kahan variant of the QL or QR algorithm.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix.  N >= 0.   
   > \endverbatim   
   >   
   > \param[in,out] D   
   > \verbatim   
   >          D is DOUBLE PRECISION array, dimension (N)   
   >          On entry, the n diagonal elements of the tridiagonal matrix.   
   >          On exit, if INFO = 0, the eigenvalues in ascending order.   
   > \endverbatim   
   >   
   > \param[in,out] E   
   > \verbatim   
   >          E is DOUBLE PRECISION array, dimension (N-1)   
   >          On entry, the (n-1) subdiagonal elements of the tridiagonal   
   >          matrix.   
   >          On exit, E has been destroyed.   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
   >          = 0:  successful exit   
   >          < 0:  if INFO = -i, the i-th argument had an illegal value   
   >          > 0:  the algorithm failed to find all of the eigenvalues in   
   >                a total of 30*N iterations; if INFO = i, then i   
   >                elements of E have not converged to zero.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date November 2011   

   > \ingroup auxOTHERcomputational   

    =====================================================================   
   Subroutine */ int igraphdsterf_(integer *n, doublereal *d__, doublereal *e, 
	integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    doublereal c__;
    integer i__, l, m;
    doublereal p, r__, s;
    integer l1;
    doublereal bb, rt1, rt2, eps, rte;
    integer lsv;
    doublereal eps2, oldc;
    integer lend;
    doublereal rmax;
    integer jtot;
    extern /* Subroutine */ int igraphdlae2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    doublereal gamma, alpha, sigma, anorm;
    extern doublereal igraphdlapy2_(doublereal *, doublereal *), igraphdlamch_(char *);
    integer iscale;
    extern /* Subroutine */ int igraphdlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *);
    doublereal oldgam, safmin;
    extern /* Subroutine */ int igraphxerbla_(char *, integer *, ftnlen);
    doublereal safmax;
    extern doublereal igraphdlanst_(char *, integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int igraphdlasrt_(char *, integer *, doublereal *, 
	    integer *);
    integer lendsv;
    doublereal ssfmin;
    integer nmaxit;
    doublereal ssfmax;


/*  -- LAPACK computational routine (version 3.4.0) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       November 2011   


    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    --e;
    --d__;

    /* Function Body */
    *info = 0;

/*     Quick return if possible */

    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	igraphxerbla_("DSTERF", &i__1, (ftnlen)6);
	return 0;
    }
    if (*n <= 1) {
	return 0;
    }

/*     Determine the unit roundoff for this environment. */

    eps = igraphdlamch_("E");
/* Computing 2nd power */
    d__1 = eps;
    eps2 = d__1 * d__1;
    safmin = igraphdlamch_("S");
    safmax = 1. / safmin;
    ssfmax = sqrt(safmax) / 3.;
    ssfmin = sqrt(safmin) / eps2;
    rmax = igraphdlamch_("O");

/*     Compute the eigenvalues of the tridiagonal matrix. */

    nmaxit = *n * 30;
    sigma = 0.;
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration   
       for each block, according to whether top or bottom diagonal   
       element is smaller. */

    l1 = 1;

L10:
    if (l1 > *n) {
	goto L170;
    }
    if (l1 > 1) {
	e[l1 - 1] = 0.;
    }
    i__1 = *n - 1;
    for (m = l1; m <= i__1; ++m) {
	if ((d__3 = e[m], abs(d__3)) <= sqrt((d__1 = d__[m], abs(d__1))) * 
		sqrt((d__2 = d__[m + 1], abs(d__2))) * eps) {
	    e[m] = 0.;
	    goto L30;
	}
/* L20: */
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

/*     Scale submatrix in rows and columns L to LEND */

    i__1 = lend - l + 1;
    anorm = igraphdlanst_("M", &i__1, &d__[l], &e[l]);
    iscale = 0;
    if (anorm == 0.) {
	goto L10;
    }
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	igraphdlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	igraphdlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	igraphdlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	igraphdlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info);
    }

    i__1 = lend - 1;
    for (i__ = l; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = e[i__];
	e[i__] = d__1 * d__1;
/* L40: */
    }

/*     Choose between QL and QR iteration */

    if ((d__1 = d__[lend], abs(d__1)) < (d__2 = d__[l], abs(d__2))) {
	lend = lsv;
	l = lendsv;
    }

    if (lend >= l) {

/*        QL Iteration   

          Look for small subdiagonal element. */

L50:
	if (l != lend) {
	    i__1 = lend - 1;
	    for (m = l; m <= i__1; ++m) {
		if ((d__2 = e[m], abs(d__2)) <= eps2 * (d__1 = d__[m] * d__[m 
			+ 1], abs(d__1))) {
		    goto L70;
		}
/* L60: */
	    }
	}
	m = lend;

L70:
	if (m < lend) {
	    e[m] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L90;
	}

/*        If remaining matrix is 2 by 2, use DLAE2 to compute its   
          eigenvalues. */

	if (m == l + 1) {
	    rte = sqrt(e[l]);
	    igraphdlae2_(&d__[l], &rte, &d__[l + 1], &rt1, &rt2);
	    d__[l] = rt1;
	    d__[l + 1] = rt2;
	    e[l] = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L50;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

/*        Form shift. */

	rte = sqrt(e[l]);
	sigma = (d__[l + 1] - p) / (rte * 2.);
	r__ = igraphdlapy2_(&sigma, &c_b33);
	sigma = p - rte / (sigma + d_sign(&r__, &sigma));

	c__ = 1.;
	s = 0.;
	gamma = d__[m] - sigma;
	p = gamma * gamma;

/*        Inner loop */

	i__1 = l;
	for (i__ = m - 1; i__ >= i__1; --i__) {
	    bb = e[i__];
	    r__ = p + bb;
	    if (i__ != m - 1) {
		e[i__ + 1] = s * r__;
	    }
	    oldc = c__;
	    c__ = p / r__;
	    s = bb / r__;
	    oldgam = gamma;
	    alpha = d__[i__];
	    gamma = c__ * (alpha - sigma) - s * oldgam;
	    d__[i__ + 1] = oldgam + (alpha - gamma);
	    if (c__ != 0.) {
		p = gamma * gamma / c__;
	    } else {
		p = oldc * bb;
	    }
/* L80: */
	}

	e[l] = s * p;
	d__[l] = sigma + gamma;
	goto L50;

/*        Eigenvalue found. */

L90:
	d__[l] = p;

	++l;
	if (l <= lend) {
	    goto L50;
	}
	goto L150;

    } else {

/*        QR Iteration   

          Look for small superdiagonal element. */

L100:
	i__1 = lend + 1;
	for (m = l; m >= i__1; --m) {
	    if ((d__2 = e[m - 1], abs(d__2)) <= eps2 * (d__1 = d__[m] * d__[m 
		    - 1], abs(d__1))) {
		goto L120;
	    }
/* L110: */
	}
	m = lend;

L120:
	if (m > lend) {
	    e[m - 1] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L140;
	}

/*        If remaining matrix is 2 by 2, use DLAE2 to compute its   
          eigenvalues. */

	if (m == l - 1) {
	    rte = sqrt(e[l - 1]);
	    igraphdlae2_(&d__[l], &rte, &d__[l - 1], &rt1, &rt2);
	    d__[l] = rt1;
	    d__[l - 1] = rt2;
	    e[l - 1] = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L100;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

/*        Form shift. */

	rte = sqrt(e[l - 1]);
	sigma = (d__[l - 1] - p) / (rte * 2.);
	r__ = igraphdlapy2_(&sigma, &c_b33);
	sigma = p - rte / (sigma + d_sign(&r__, &sigma));

	c__ = 1.;
	s = 0.;
	gamma = d__[m] - sigma;
	p = gamma * gamma;

/*        Inner loop */

	i__1 = l - 1;
	for (i__ = m; i__ <= i__1; ++i__) {
	    bb = e[i__];
	    r__ = p + bb;
	    if (i__ != m) {
		e[i__ - 1] = s * r__;
	    }
	    oldc = c__;
	    c__ = p / r__;
	    s = bb / r__;
	    oldgam = gamma;
	    alpha = d__[i__ + 1];
	    gamma = c__ * (alpha - sigma) - s * oldgam;
	    d__[i__] = oldgam + (alpha - gamma);
	    if (c__ != 0.) {
		p = gamma * gamma / c__;
	    } else {
		p = oldc * bb;
	    }
/* L130: */
	}

	e[l - 1] = s * p;
	d__[l] = sigma + gamma;
	goto L100;

/*        Eigenvalue found. */

L140:
	d__[l] = p;

	--l;
	if (l >= lend) {
	    goto L100;
	}
	goto L150;

    }

/*     Undo scaling if necessary */

L150:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	igraphdlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
    }
    if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	igraphdlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
    }

/*     Check for no convergence to an eigenvalue after a total   
       of N*MAXIT iterations. */

    if (jtot < nmaxit) {
	goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (e[i__] != 0.) {
	    ++(*info);
	}
/* L160: */
    }
    goto L180;

/*     Sort eigenvalues in increasing order. */

L170:
    igraphdlasrt_("I", n, &d__[1], info);

L180:
    return 0;

/*     End of DSTERF */

} /* igraphdsterf_ */

