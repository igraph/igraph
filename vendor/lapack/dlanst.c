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

/* > \brief \b DLANST returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a real symmetric tridiagonal matrix.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLANST + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanst.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanst.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanst.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )   

         CHARACTER          NORM   
         INTEGER            N   
         DOUBLE PRECISION   D( * ), E( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DLANST  returns the value of the one norm,  or the Frobenius norm, or   
   > the  infinity norm,  or the  element of  largest absolute value  of a   
   > real symmetric tridiagonal matrix A.   
   > \endverbatim   
   >   
   > \return DLANST   
   > \verbatim   
   >   
   >    DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
   >             (   
   >             ( norm1(A),         NORM = '1', 'O' or 'o'   
   >             (   
   >             ( normI(A),         NORM = 'I' or 'i'   
   >             (   
   >             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   
   >   
   > where  norm1  denotes the  one norm of a matrix (maximum column sum),   
   > normI  denotes the  infinity norm  of a matrix  (maximum row sum) and   
   > normF  denotes the  Frobenius norm of a matrix (square root of sum of   
   > squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] NORM   
   > \verbatim   
   >          NORM is CHARACTER*1   
   >          Specifies the value to be returned in DLANST as described   
   >          above.   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix A.  N >= 0.  When N = 0, DLANST is   
   >          set to zero.   
   > \endverbatim   
   >   
   > \param[in] D   
   > \verbatim   
   >          D is DOUBLE PRECISION array, dimension (N)   
   >          The diagonal elements of A.   
   > \endverbatim   
   >   
   > \param[in] E   
   > \verbatim   
   >          E is DOUBLE PRECISION array, dimension (N-1)   
   >          The (n-1) sub-diagonal or super-diagonal elements of A.   
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
doublereal igraphdlanst_(char *norm, integer *n, doublereal *d__, doublereal *e)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__;
    doublereal sum, scale;
    extern logical igraphlsame_(char *, char *);
    doublereal anorm;
    extern logical igraphdisnan_(doublereal *);
    extern /* Subroutine */ int igraphdlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       Parameter adjustments */
    --e;
    --d__;

    /* Function Body */
    if (*n <= 0) {
	anorm = 0.;
    } else if (igraphlsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	anorm = (d__1 = d__[*n], abs(d__1));
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sum = (d__1 = d__[i__], abs(d__1));
	    if (anorm < sum || igraphdisnan_(&sum)) {
		anorm = sum;
	    }
	    sum = (d__1 = e[i__], abs(d__1));
	    if (anorm < sum || igraphdisnan_(&sum)) {
		anorm = sum;
	    }
/* L10: */
	}
    } else if (igraphlsame_(norm, "O") || *(unsigned char *)
	    norm == '1' || igraphlsame_(norm, "I")) {

/*        Find norm1(A). */

	if (*n == 1) {
	    anorm = abs(d__[1]);
	} else {
	    anorm = abs(d__[1]) + abs(e[1]);
	    sum = (d__1 = e[*n - 1], abs(d__1)) + (d__2 = d__[*n], abs(d__2));
	    if (anorm < sum || igraphdisnan_(&sum)) {
		anorm = sum;
	    }
	    i__1 = *n - 1;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		sum = (d__1 = d__[i__], abs(d__1)) + (d__2 = e[i__], abs(d__2)
			) + (d__3 = e[i__ - 1], abs(d__3));
		if (anorm < sum || igraphdisnan_(&sum)) {
		    anorm = sum;
		}
/* L20: */
	    }
	}
    } else if (igraphlsame_(norm, "F") || igraphlsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	if (*n > 1) {
	    i__1 = *n - 1;
	    igraphdlassq_(&i__1, &e[1], &c__1, &scale, &sum);
	    sum *= 2;
	}
	igraphdlassq_(n, &d__[1], &c__1, &scale, &sum);
	anorm = scale * sqrt(sum);
    }

    ret_val = anorm;
    return ret_val;

/*     End of DLANST */

} /* igraphdlanst_ */

