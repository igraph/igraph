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

/* > \brief \b DLARRJ performs refinement of the initial estimates of the eigenvalues of the matrix T.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLARRJ + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrj.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrj.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrj.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DLARRJ( N, D, E2, IFIRST, ILAST,   
                            RTOL, OFFSET, W, WERR, WORK, IWORK,   
                            PIVMIN, SPDIAM, INFO )   

         INTEGER            IFIRST, ILAST, INFO, N, OFFSET   
         DOUBLE PRECISION   PIVMIN, RTOL, SPDIAM   
         INTEGER            IWORK( * )   
         DOUBLE PRECISION   D( * ), E2( * ), W( * ),   
        $                   WERR( * ), WORK( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > Given the initial eigenvalue approximations of T, DLARRJ   
   > does  bisection to refine the eigenvalues of T,   
   > W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial   
   > guesses for these eigenvalues are input in W, the corresponding estimate   
   > of the error in these guesses in WERR. During bisection, intervals   
   > [left, right] are maintained by storing their mid-points and   
   > semi-widths in the arrays W and WERR respectively.   
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
   >          The N diagonal elements of T.   
   > \endverbatim   
   >   
   > \param[in] E2   
   > \verbatim   
   >          E2 is DOUBLE PRECISION array, dimension (N-1)   
   >          The Squares of the (N-1) subdiagonal elements of T.   
   > \endverbatim   
   >   
   > \param[in] IFIRST   
   > \verbatim   
   >          IFIRST is INTEGER   
   >          The index of the first eigenvalue to be computed.   
   > \endverbatim   
   >   
   > \param[in] ILAST   
   > \verbatim   
   >          ILAST is INTEGER   
   >          The index of the last eigenvalue to be computed.   
   > \endverbatim   
   >   
   > \param[in] RTOL   
   > \verbatim   
   >          RTOL is DOUBLE PRECISION   
   >          Tolerance for the convergence of the bisection intervals.   
   >          An interval [LEFT,RIGHT] has converged if   
   >          RIGHT-LEFT.LT.RTOL*MAX(|LEFT|,|RIGHT|).   
   > \endverbatim   
   >   
   > \param[in] OFFSET   
   > \verbatim   
   >          OFFSET is INTEGER   
   >          Offset for the arrays W and WERR, i.e., the IFIRST-OFFSET   
   >          through ILAST-OFFSET elements of these arrays are to be used.   
   > \endverbatim   
   >   
   > \param[in,out] W   
   > \verbatim   
   >          W is DOUBLE PRECISION array, dimension (N)   
   >          On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are   
   >          estimates of the eigenvalues of L D L^T indexed IFIRST through   
   >          ILAST.   
   >          On output, these estimates are refined.   
   > \endverbatim   
   >   
   > \param[in,out] WERR   
   > \verbatim   
   >          WERR is DOUBLE PRECISION array, dimension (N)   
   >          On input, WERR( IFIRST-OFFSET ) through WERR( ILAST-OFFSET ) are   
   >          the errors in the estimates of the corresponding elements in W.   
   >          On output, these errors are refined.   
   > \endverbatim   
   >   
   > \param[out] WORK   
   > \verbatim   
   >          WORK is DOUBLE PRECISION array, dimension (2*N)   
   >          Workspace.   
   > \endverbatim   
   >   
   > \param[out] IWORK   
   > \verbatim   
   >          IWORK is INTEGER array, dimension (2*N)   
   >          Workspace.   
   > \endverbatim   
   >   
   > \param[in] PIVMIN   
   > \verbatim   
   >          PIVMIN is DOUBLE PRECISION   
   >          The minimum pivot in the Sturm sequence for T.   
   > \endverbatim   
   >   
   > \param[in] SPDIAM   
   > \verbatim   
   >          SPDIAM is DOUBLE PRECISION   
   >          The spectral diameter of T.   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
   >          Error flag.   
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
   > Beresford Parlett, University of California, Berkeley, USA \n   
   > Jim Demmel, University of California, Berkeley, USA \n   
   > Inderjit Dhillon, University of Texas, Austin, USA \n   
   > Osni Marques, LBNL/NERSC, USA \n   
   > Christof Voemel, University of California, Berkeley, USA   

    =====================================================================   
   Subroutine */ int igraphdlarrj_(integer *n, doublereal *d__, doublereal *e2, 
	integer *ifirst, integer *ilast, doublereal *rtol, integer *offset, 
	doublereal *w, doublereal *werr, doublereal *work, integer *iwork, 
	doublereal *pivmin, doublereal *spdiam, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    integer i__, j, k, p;
    doublereal s;
    integer i1, i2, ii;
    doublereal fac, mid;
    integer cnt;
    doublereal tmp, left;
    integer iter, nint, prev, next, savi1;
    doublereal right, width, dplus;
    integer olnint, maxitr;


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   



       Parameter adjustments */
    --iwork;
    --work;
    --werr;
    --w;
    --e2;
    --d__;

    /* Function Body */
    *info = 0;

    maxitr = (integer) ((log(*spdiam + *pivmin) - log(*pivmin)) / log(2.)) + 
	    2;

/*     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].   
       The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while   
       Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 )   
       for an unconverged interval is set to the index of the next unconverged   
       interval, and is -1 or 0 for a converged interval. Thus a linked   
       list of unconverged intervals is set up. */

    i1 = *ifirst;
    i2 = *ilast;
/*     The number of unconverged intervals */
    nint = 0;
/*     The last unconverged interval found */
    prev = 0;
    i__1 = i2;
    for (i__ = i1; i__ <= i__1; ++i__) {
	k = i__ << 1;
	ii = i__ - *offset;
	left = w[ii] - werr[ii];
	mid = w[ii];
	right = w[ii] + werr[ii];
	width = right - mid;
/* Computing MAX */
	d__1 = abs(left), d__2 = abs(right);
	tmp = max(d__1,d__2);
/*        The following test prevents the test of converged intervals */
	if (width < *rtol * tmp) {
/*           This interval has already converged and does not need refinement.   
             (Note that the gaps might change through refining the   
              eigenvalues, however, they can only get bigger.)   
             Remove it from the list. */
	    iwork[k - 1] = -1;
/*           Make sure that I1 always points to the first unconverged interval */
	    if (i__ == i1 && i__ < i2) {
		i1 = i__ + 1;
	    }
	    if (prev >= i1 && i__ <= i2) {
		iwork[(prev << 1) - 1] = i__ + 1;
	    }
	} else {
/*           unconverged interval found */
	    prev = i__;
/*           Make sure that [LEFT,RIGHT] contains the desired eigenvalue   

             Do while( CNT(LEFT).GT.I-1 ) */

	    fac = 1.;
L20:
	    cnt = 0;
	    s = left;
	    dplus = d__[1] - s;
	    if (dplus < 0.) {
		++cnt;
	    }
	    i__2 = *n;
	    for (j = 2; j <= i__2; ++j) {
		dplus = d__[j] - s - e2[j - 1] / dplus;
		if (dplus < 0.) {
		    ++cnt;
		}
/* L30: */
	    }
	    if (cnt > i__ - 1) {
		left -= werr[ii] * fac;
		fac *= 2.;
		goto L20;
	    }

/*           Do while( CNT(RIGHT).LT.I ) */

	    fac = 1.;
L50:
	    cnt = 0;
	    s = right;
	    dplus = d__[1] - s;
	    if (dplus < 0.) {
		++cnt;
	    }
	    i__2 = *n;
	    for (j = 2; j <= i__2; ++j) {
		dplus = d__[j] - s - e2[j - 1] / dplus;
		if (dplus < 0.) {
		    ++cnt;
		}
/* L60: */
	    }
	    if (cnt < i__) {
		right += werr[ii] * fac;
		fac *= 2.;
		goto L50;
	    }
	    ++nint;
	    iwork[k - 1] = i__ + 1;
	    iwork[k] = cnt;
	}
	work[k - 1] = left;
	work[k] = right;
/* L75: */
    }
    savi1 = i1;

/*     Do while( NINT.GT.0 ), i.e. there are still unconverged intervals   
       and while (ITER.LT.MAXITR) */

    iter = 0;
L80:
    prev = i1 - 1;
    i__ = i1;
    olnint = nint;
    i__1 = olnint;
    for (p = 1; p <= i__1; ++p) {
	k = i__ << 1;
	ii = i__ - *offset;
	next = iwork[k - 1];
	left = work[k - 1];
	right = work[k];
	mid = (left + right) * .5;
/*        semiwidth of interval */
	width = right - mid;
/* Computing MAX */
	d__1 = abs(left), d__2 = abs(right);
	tmp = max(d__1,d__2);
	if (width < *rtol * tmp || iter == maxitr) {
/*           reduce number of unconverged intervals */
	    --nint;
/*           Mark interval as converged. */
	    iwork[k - 1] = 0;
	    if (i1 == i__) {
		i1 = next;
	    } else {
/*              Prev holds the last unconverged interval previously examined */
		if (prev >= i1) {
		    iwork[(prev << 1) - 1] = next;
		}
	    }
	    i__ = next;
	    goto L100;
	}
	prev = i__;

/*        Perform one bisection step */

	cnt = 0;
	s = mid;
	dplus = d__[1] - s;
	if (dplus < 0.) {
	    ++cnt;
	}
	i__2 = *n;
	for (j = 2; j <= i__2; ++j) {
	    dplus = d__[j] - s - e2[j - 1] / dplus;
	    if (dplus < 0.) {
		++cnt;
	    }
/* L90: */
	}
	if (cnt <= i__ - 1) {
	    work[k - 1] = mid;
	} else {
	    work[k] = mid;
	}
	i__ = next;
L100:
	;
    }
    ++iter;
/*     do another loop if there are still unconverged intervals   
       However, in the last iteration, all intervals are accepted   
       since this is the best we can do. */
    if (nint > 0 && iter <= maxitr) {
	goto L80;
    }


/*     At this point, all the intervals have converged */
    i__1 = *ilast;
    for (i__ = savi1; i__ <= i__1; ++i__) {
	k = i__ << 1;
	ii = i__ - *offset;
/*        All intervals marked by '0' have been refined. */
	if (iwork[k - 1] == 0) {
	    w[ii] = (work[k - 1] + work[k]) * .5;
	    werr[ii] = work[k] - w[ii];
	}
/* L110: */
    }

    return 0;

/*     End of DLARRJ */

} /* igraphdlarrj_ */

