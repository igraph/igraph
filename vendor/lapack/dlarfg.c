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

/* > \brief \b DLARFG generates an elementary reflector (Householder matrix).   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLARFG + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfg.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfg.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfg.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )   

         INTEGER            INCX, N   
         DOUBLE PRECISION   ALPHA, TAU   
         DOUBLE PRECISION   X( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DLARFG generates a real elementary reflector H of order n, such   
   > that   
   >   
   >       H * ( alpha ) = ( beta ),   H**T * H = I.   
   >           (   x   )   (   0  )   
   >   
   > where alpha and beta are scalars, and x is an (n-1)-element real   
   > vector. H is represented in the form   
   >   
   >       H = I - tau * ( 1 ) * ( 1 v**T ) ,   
   >                     ( v )   
   >   
   > where tau is a real scalar and v is a real (n-1)-element   
   > vector.   
   >   
   > If the elements of x are all zero, then tau = 0 and H is taken to be   
   > the unit matrix.   
   >   
   > Otherwise  1 <= tau <= 2.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the elementary reflector.   
   > \endverbatim   
   >   
   > \param[in,out] ALPHA   
   > \verbatim   
   >          ALPHA is DOUBLE PRECISION   
   >          On entry, the value alpha.   
   >          On exit, it is overwritten with the value beta.   
   > \endverbatim   
   >   
   > \param[in,out] X   
   > \verbatim   
   >          X is DOUBLE PRECISION array, dimension   
   >                         (1+(N-2)*abs(INCX))   
   >          On entry, the vector x.   
   >          On exit, it is overwritten with the vector v.   
   > \endverbatim   
   >   
   > \param[in] INCX   
   > \verbatim   
   >          INCX is INTEGER   
   >          The increment between elements of X. INCX > 0.   
   > \endverbatim   
   >   
   > \param[out] TAU   
   > \verbatim   
   >          TAU is DOUBLE PRECISION   
   >          The value tau.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup doubleOTHERauxiliary   

    =====================================================================   
   Subroutine */ int igraphdlarfg_(integer *n, doublereal *alpha, doublereal *x, 
	integer *incx, doublereal *tau)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    integer j, knt;
    doublereal beta;
    extern doublereal igraphdnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int igraphdscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    doublereal xnorm;
    extern doublereal igraphdlapy2_(doublereal *, doublereal *), igraphdlamch_(char *);
    doublereal safmin, rsafmn;


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       Parameter adjustments */
    --x;

    /* Function Body */
    if (*n <= 1) {
	*tau = 0.;
	return 0;
    }

    i__1 = *n - 1;
    xnorm = igraphdnrm2_(&i__1, &x[1], incx);

    if (xnorm == 0.) {

/*        H  =  I */

	*tau = 0.;
    } else {

/*        general case */

	d__1 = igraphdlapy2_(alpha, &xnorm);
	beta = -d_sign(&d__1, alpha);
	safmin = igraphdlamch_("S") / igraphdlamch_("E");
	knt = 0;
	if (abs(beta) < safmin) {

/*           XNORM, BETA may be inaccurate; scale X and recompute them */

	    rsafmn = 1. / safmin;
L10:
	    ++knt;
	    i__1 = *n - 1;
	    igraphdscal_(&i__1, &rsafmn, &x[1], incx);
	    beta *= rsafmn;
	    *alpha *= rsafmn;
	    if (abs(beta) < safmin) {
		goto L10;
	    }

/*           New BETA is at most 1, at least SAFMIN */

	    i__1 = *n - 1;
	    xnorm = igraphdnrm2_(&i__1, &x[1], incx);
	    d__1 = igraphdlapy2_(alpha, &xnorm);
	    beta = -d_sign(&d__1, alpha);
	}
	*tau = (beta - *alpha) / beta;
	i__1 = *n - 1;
	d__1 = 1. / (*alpha - beta);
	igraphdscal_(&i__1, &d__1, &x[1], incx);

/*        If ALPHA is subnormal, it may lose relative accuracy */

	i__1 = knt;
	for (j = 1; j <= i__1; ++j) {
	    beta *= safmin;
/* L20: */
	}
	*alpha = beta;
    }

    return 0;

/*     End of DLARFG */

} /* igraphdlarfg_ */

