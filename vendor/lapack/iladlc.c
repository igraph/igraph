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

/* > \brief \b ILADLC scans a matrix for its last non-zero column.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download ILADLC + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlc.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladlc.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladlc.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         INTEGER FUNCTION ILADLC( M, N, A, LDA )   

         INTEGER            M, N, LDA   
         DOUBLE PRECISION   A( LDA, * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > ILADLC scans A for its last non-zero column.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] M   
   > \verbatim   
   >          M is INTEGER   
   >          The number of rows of the matrix A.   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The number of columns of the matrix A.   
   > \endverbatim   
   >   
   > \param[in] A   
   > \verbatim   
   >          A is DOUBLE PRECISION array, dimension (LDA,N)   
   >          The m by n matrix A.   
   > \endverbatim   
   >   
   > \param[in] LDA   
   > \verbatim   
   >          LDA is INTEGER   
   >          The leading dimension of the array A. LDA >= max(1,M).   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup auxOTHERauxiliary   

    ===================================================================== */
integer igraphiladlc_(integer *m, integer *n, doublereal *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, ret_val, i__1;

    /* Local variables */
    integer i__;


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       Quick test for the common case where one corner is non-zero.   
       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    if (*n == 0) {
	ret_val = *n;
    } else if (a[*n * a_dim1 + 1] != 0. || a[*m + *n * a_dim1] != 0.) {
	ret_val = *n;
    } else {
/*     Now scan each column from the end, returning with the first non-zero. */
	for (ret_val = *n; ret_val >= 1; --ret_val) {
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (a[i__ + ret_val * a_dim1] != 0.) {
		    return ret_val;
		}
	    }
	}
    }
    return ret_val;
} /* igraphiladlc_ */

