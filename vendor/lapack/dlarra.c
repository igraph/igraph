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

/* > \brief \b DLARRA computes the splitting points with the specified threshold.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLARRA + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarra.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarra.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarra.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DLARRA( N, D, E, E2, SPLTOL, TNRM,   
                             NSPLIT, ISPLIT, INFO )   

         INTEGER            INFO, N, NSPLIT   
         DOUBLE PRECISION    SPLTOL, TNRM   
         INTEGER            ISPLIT( * )   
         DOUBLE PRECISION   D( * ), E( * ), E2( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > Compute the splitting points with threshold SPLTOL.   
   > DLARRA sets any "small" off-diagonal elements to zero.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix. N > 0.   
   > \endverbatim   
   >   
   > \param[in] D   
   > \verbatim   
   >          D is DOUBLE PRECISION array, dimension (N)   
   >          On entry, the N diagonal elements of the tridiagonal   
   >          matrix T.   
   > \endverbatim   
   >   
   > \param[in,out] E   
   > \verbatim   
   >          E is DOUBLE PRECISION array, dimension (N)   
   >          On entry, the first (N-1) entries contain the subdiagonal   
   >          elements of the tridiagonal matrix T; E(N) need not be set.   
   >          On exit, the entries E( ISPLIT( I ) ), 1 <= I <= NSPLIT,   
   >          are set to zero, the other entries of E are untouched.   
   > \endverbatim   
   >   
   > \param[in,out] E2   
   > \verbatim   
   >          E2 is DOUBLE PRECISION array, dimension (N)   
   >          On entry, the first (N-1) entries contain the SQUARES of the   
   >          subdiagonal elements of the tridiagonal matrix T;   
   >          E2(N) need not be set.   
   >          On exit, the entries E2( ISPLIT( I ) ),   
   >          1 <= I <= NSPLIT, have been set to zero   
   > \endverbatim   
   >   
   > \param[in] SPLTOL   
   > \verbatim   
   >          SPLTOL is DOUBLE PRECISION   
   >          The threshold for splitting. Two criteria can be used:   
   >          SPLTOL<0 : criterion based on absolute off-diagonal value   
   >          SPLTOL>0 : criterion that preserves relative accuracy   
   > \endverbatim   
   >   
   > \param[in] TNRM   
   > \verbatim   
   >          TNRM is DOUBLE PRECISION   
   >          The norm of the matrix.   
   > \endverbatim   
   >   
   > \param[out] NSPLIT   
   > \verbatim   
   >          NSPLIT is INTEGER   
   >          The number of blocks T splits into. 1 <= NSPLIT <= N.   
   > \endverbatim   
   >   
   > \param[out] ISPLIT   
   > \verbatim   
   >          ISPLIT is INTEGER array, dimension (N)   
   >          The splitting points, at which T breaks up into blocks.   
   >          The first block consists of rows/columns 1 to ISPLIT(1),   
   >          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),   
   >          etc., and the NSPLIT-th consists of rows/columns   
   >          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
   >          = 0:  successful exit   
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
   Subroutine */ int igraphdlarra_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *e2, doublereal *spltol, doublereal *tnrm, integer *nsplit,
	 integer *isplit, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__;
    doublereal tmp1, eabs;


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       Parameter adjustments */
    --isplit;
    --e2;
    --e;
    --d__;

    /* Function Body */
    *info = 0;
/*     Compute splitting points */
    *nsplit = 1;
    if (*spltol < 0.) {
/*        Criterion based on absolute off-diagonal value */
	tmp1 = abs(*spltol) * *tnrm;
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    eabs = (d__1 = e[i__], abs(d__1));
	    if (eabs <= tmp1) {
		e[i__] = 0.;
		e2[i__] = 0.;
		isplit[*nsplit] = i__;
		++(*nsplit);
	    }
/* L9: */
	}
    } else {
/*        Criterion that guarantees relative accuracy */
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    eabs = (d__1 = e[i__], abs(d__1));
	    if (eabs <= *spltol * sqrt((d__1 = d__[i__], abs(d__1))) * sqrt((
		    d__2 = d__[i__ + 1], abs(d__2)))) {
		e[i__] = 0.;
		e2[i__] = 0.;
		isplit[*nsplit] = i__;
		++(*nsplit);
	    }
/* L10: */
	}
    }
    isplit[*nsplit] = *n;
    return 0;

/*     End of DLARRA */

} /* igraphdlarra_ */

