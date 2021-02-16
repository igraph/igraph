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

/* > \brief \b DLARRC computes the number of eigenvalues of the symmetric tridiagonal matrix.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLARRC + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrc.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrc.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrc.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DLARRC( JOBT, N, VL, VU, D, E, PIVMIN,   
                                     EIGCNT, LCNT, RCNT, INFO )   

         CHARACTER          JOBT   
         INTEGER            EIGCNT, INFO, LCNT, N, RCNT   
         DOUBLE PRECISION   PIVMIN, VL, VU   
         DOUBLE PRECISION   D( * ), E( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > Find the number of eigenvalues of the symmetric tridiagonal matrix T   
   > that are in the interval (VL,VU] if JOBT = 'T', and of L D L^T   
   > if JOBT = 'L'.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] JOBT   
   > \verbatim   
   >          JOBT is CHARACTER*1   
   >          = 'T':  Compute Sturm count for matrix T.   
   >          = 'L':  Compute Sturm count for matrix L D L^T.   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix. N > 0.   
   > \endverbatim   
   >   
   > \param[in] VL   
   > \verbatim   
   >          VL is DOUBLE PRECISION   
   > \endverbatim   
   >   
   > \param[in] VU   
   > \verbatim   
   >          VU is DOUBLE PRECISION   
   >          The lower and upper bounds for the eigenvalues.   
   > \endverbatim   
   >   
   > \param[in] D   
   > \verbatim   
   >          D is DOUBLE PRECISION array, dimension (N)   
   >          JOBT = 'T': The N diagonal elements of the tridiagonal matrix T.   
   >          JOBT = 'L': The N diagonal elements of the diagonal matrix D.   
   > \endverbatim   
   >   
   > \param[in] E   
   > \verbatim   
   >          E is DOUBLE PRECISION array, dimension (N)   
   >          JOBT = 'T': The N-1 offdiagonal elements of the matrix T.   
   >          JOBT = 'L': The N-1 offdiagonal elements of the matrix L.   
   > \endverbatim   
   >   
   > \param[in] PIVMIN   
   > \verbatim   
   >          PIVMIN is DOUBLE PRECISION   
   >          The minimum pivot in the Sturm sequence for T.   
   > \endverbatim   
   >   
   > \param[out] EIGCNT   
   > \verbatim   
   >          EIGCNT is INTEGER   
   >          The number of eigenvalues of the symmetric tridiagonal matrix T   
   >          that are in the interval (VL,VU]   
   > \endverbatim   
   >   
   > \param[out] LCNT   
   > \verbatim   
   >          LCNT is INTEGER   
   > \endverbatim   
   >   
   > \param[out] RCNT   
   > \verbatim   
   >          RCNT is INTEGER   
   >          The left and right negcounts of the interval.   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
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
   Subroutine */ int igraphdlarrc_(char *jobt, integer *n, doublereal *vl, 
	doublereal *vu, doublereal *d__, doublereal *e, doublereal *pivmin, 
	integer *eigcnt, integer *lcnt, integer *rcnt, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    integer i__;
    doublereal sl, su, tmp, tmp2;
    logical matt;
    extern logical igraphlsame_(char *, char *);
    doublereal lpivot, rpivot;


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       Parameter adjustments */
    --e;
    --d__;

    /* Function Body */
    *info = 0;
    *lcnt = 0;
    *rcnt = 0;
    *eigcnt = 0;
    matt = igraphlsame_(jobt, "T");
    if (matt) {
/*        Sturm sequence count on T */
	lpivot = d__[1] - *vl;
	rpivot = d__[1] - *vu;
	if (lpivot <= 0.) {
	    ++(*lcnt);
	}
	if (rpivot <= 0.) {
	    ++(*rcnt);
	}
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = e[i__];
	    tmp = d__1 * d__1;
	    lpivot = d__[i__ + 1] - *vl - tmp / lpivot;
	    rpivot = d__[i__ + 1] - *vu - tmp / rpivot;
	    if (lpivot <= 0.) {
		++(*lcnt);
	    }
	    if (rpivot <= 0.) {
		++(*rcnt);
	    }
/* L10: */
	}
    } else {
/*        Sturm sequence count on L D L^T */
	sl = -(*vl);
	su = -(*vu);
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lpivot = d__[i__] + sl;
	    rpivot = d__[i__] + su;
	    if (lpivot <= 0.) {
		++(*lcnt);
	    }
	    if (rpivot <= 0.) {
		++(*rcnt);
	    }
	    tmp = e[i__] * d__[i__] * e[i__];

	    tmp2 = tmp / lpivot;
	    if (tmp2 == 0.) {
		sl = tmp - *vl;
	    } else {
		sl = sl * tmp2 - *vl;
	    }

	    tmp2 = tmp / rpivot;
	    if (tmp2 == 0.) {
		su = tmp - *vu;
	    } else {
		su = su * tmp2 - *vu;
	    }
/* L20: */
	}
	lpivot = d__[*n] + sl;
	rpivot = d__[*n] + su;
	if (lpivot <= 0.) {
	    ++(*lcnt);
	}
	if (rpivot <= 0.) {
	    ++(*rcnt);
	}
    }
    *eigcnt = *rcnt - *lcnt;
    return 0;

/*     end of DLARRC */

} /* igraphdlarrc_ */

