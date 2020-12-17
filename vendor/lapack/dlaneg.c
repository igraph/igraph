/*  -- translated by f2c (version 20191129).
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

/* > \brief \b DLANEG computes the Sturm count.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLANEG + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaneg.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaneg.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaneg.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         INTEGER FUNCTION DLANEG( N, D, LLD, SIGMA, PIVMIN, R )   

         INTEGER            N, R   
         DOUBLE PRECISION   PIVMIN, SIGMA   
         DOUBLE PRECISION   D( * ), LLD( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DLANEG computes the Sturm count, the number of negative pivots   
   > encountered while factoring tridiagonal T - sigma I = L D L^T.   
   > This implementation works directly on the factors without forming   
   > the tridiagonal matrix T.  The Sturm count is also the number of   
   > eigenvalues of T less than sigma.   
   >   
   > This routine is called from DLARRB.   
   >   
   > The current routine does not use the PIVMIN parameter but rather   
   > requires IEEE-754 propagation of Infinities and NaNs.  This   
   > routine also has no input range restrictions but does require   
   > default exception handling such that x/0 produces Inf when x is   
   > non-zero, and Inf/Inf produces NaN.  For more information, see:   
   >   
   >   Marques, Riedy, and Voemel, "Benefits of IEEE-754 Features in   
   >   Modern Symmetric Tridiagonal Eigensolvers," SIAM Journal on   
   >   Scientific Computing, v28, n5, 2006.  DOI 10.1137/050641624   
   >   (Tech report version in LAWN 172 with the same title.)   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix.   
   > \endverbatim   
   >   
   > \param[in] D   
   > \verbatim   
   >          D is DOUBLE PRECISION array, dimension (N)   
   >          The N diagonal elements of the diagonal matrix D.   
   > \endverbatim   
   >   
   > \param[in] LLD   
   > \verbatim   
   >          LLD is DOUBLE PRECISION array, dimension (N-1)   
   >          The (N-1) elements L(i)*L(i)*D(i).   
   > \endverbatim   
   >   
   > \param[in] SIGMA   
   > \verbatim   
   >          SIGMA is DOUBLE PRECISION   
   >          Shift amount in T - sigma I = L D L^T.   
   > \endverbatim   
   >   
   > \param[in] PIVMIN   
   > \verbatim   
   >          PIVMIN is DOUBLE PRECISION   
   >          The minimum pivot in the Sturm sequence.  May be used   
   >          when zero pivots are encountered on non-IEEE-754   
   >          architectures.   
   > \endverbatim   
   >   
   > \param[in] R   
   > \verbatim   
   >          R is INTEGER   
   >          The twist index for the twisted factorization that is used   
   >          for the negcount.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup auxOTHERauxiliary   

   > \par Contributors:   
    ==================   
   >   
   >     Osni Marques, LBNL/NERSC, USA \n   
   >     Christof Voemel, University of California, Berkeley, USA \n   
   >     Jason Riedy, University of California, Berkeley, USA \n   
   >   
    ===================================================================== */
integer igraphdlaneg_(integer *n, doublereal *d__, doublereal *lld, doublereal *
	sigma, doublereal *pivmin, integer *r__)
{
    /* System generated locals */
    integer ret_val, i__1, i__2, i__3, i__4;

    /* Local variables */
    integer j;
    doublereal p, t;
    integer bj;
    doublereal tmp;
    integer neg1, neg2;
    doublereal bsav, gamma, dplus;
    extern logical igraphdisnan_(doublereal *);
    integer negcnt;
    logical sawnan;
    doublereal dminus;


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   

       Some architectures propagate Infinities and NaNs very slowly, so   
       the code computes counts in BLKLEN chunks.  Then a NaN can   
       propagate at most BLKLEN columns before being detected.  This is   
       not a general tuning parameter; it needs only to be just large   
       enough that the overhead is tiny in common cases.   
       Parameter adjustments */
    --lld;
    --d__;

    /* Function Body */
    negcnt = 0;
/*     I) upper part: L D L^T - SIGMA I = L+ D+ L+^T */
    t = -(*sigma);
    i__1 = *r__ - 1;
    for (bj = 1; bj <= i__1; bj += 128) {
	neg1 = 0;
	bsav = t;
/* Computing MIN */
	i__3 = bj + 127, i__4 = *r__ - 1;
	i__2 = min(i__3,i__4);
	for (j = bj; j <= i__2; ++j) {
	    dplus = d__[j] + t;
	    if (dplus < 0.) {
		++neg1;
	    }
	    tmp = t / dplus;
	    t = tmp * lld[j] - *sigma;
/* L21: */
	}
	sawnan = igraphdisnan_(&t);
/*     Run a slower version of the above loop if a NaN is detected.   
       A NaN should occur only with a zero pivot after an infinite   
       pivot.  In that case, substituting 1 for T/DPLUS is the   
       correct limit. */
	if (sawnan) {
	    neg1 = 0;
	    t = bsav;
/* Computing MIN */
	    i__3 = bj + 127, i__4 = *r__ - 1;
	    i__2 = min(i__3,i__4);
	    for (j = bj; j <= i__2; ++j) {
		dplus = d__[j] + t;
		if (dplus < 0.) {
		    ++neg1;
		}
		tmp = t / dplus;
		if (igraphdisnan_(&tmp)) {
		    tmp = 1.;
		}
		t = tmp * lld[j] - *sigma;
/* L22: */
	    }
	}
	negcnt += neg1;
/* L210: */
    }

/*     II) lower part: L D L^T - SIGMA I = U- D- U-^T */
    p = d__[*n] - *sigma;
    i__1 = *r__;
    for (bj = *n - 1; bj >= i__1; bj += -128) {
	neg2 = 0;
	bsav = p;
/* Computing MAX */
	i__3 = bj - 127;
	i__2 = max(i__3,*r__);
	for (j = bj; j >= i__2; --j) {
	    dminus = lld[j] + p;
	    if (dminus < 0.) {
		++neg2;
	    }
	    tmp = p / dminus;
	    p = tmp * d__[j] - *sigma;
/* L23: */
	}
	sawnan = igraphdisnan_(&p);
/*     As above, run a slower version that substitutes 1 for Inf/Inf. */

	if (sawnan) {
	    neg2 = 0;
	    p = bsav;
/* Computing MAX */
	    i__3 = bj - 127;
	    i__2 = max(i__3,*r__);
	    for (j = bj; j >= i__2; --j) {
		dminus = lld[j] + p;
		if (dminus < 0.) {
		    ++neg2;
		}
		tmp = p / dminus;
		if (igraphdisnan_(&tmp)) {
		    tmp = 1.;
		}
		p = tmp * d__[j] - *sigma;
/* L24: */
	    }
	}
	negcnt += neg2;
/* L230: */
    }

/*     III) Twist index   
         T was shifted by SIGMA initially. */
    gamma = t + *sigma + p;
    if (gamma < 0.) {
	++negcnt;
    }
    ret_val = negcnt;
    return ret_val;
} /* igraphdlaneg_ */

