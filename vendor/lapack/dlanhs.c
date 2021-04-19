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

/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLANHS returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of an upper Hessenberg matrix.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLANHS + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanhs.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanhs.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanhs.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         DOUBLE PRECISION FUNCTION DLANHS( NORM, N, A, LDA, WORK )   

         CHARACTER          NORM   
         INTEGER            LDA, N   
         DOUBLE PRECISION   A( LDA, * ), WORK( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DLANHS  returns the value of the one norm,  or the Frobenius norm, or   
   > the  infinity norm,  or the  element of  largest absolute value  of a   
   > Hessenberg matrix A.   
   > \endverbatim   
   >   
   > \return DLANHS   
   > \verbatim   
   >   
   >    DLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
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
   >          Specifies the value to be returned in DLANHS as described   
   >          above.   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix A.  N >= 0.  When N = 0, DLANHS is   
   >          set to zero.   
   > \endverbatim   
   >   
   > \param[in] A   
   > \verbatim   
   >          A is DOUBLE PRECISION array, dimension (LDA,N)   
   >          The n by n upper Hessenberg matrix A; the part of A below the   
   >          first sub-diagonal is not referenced.   
   > \endverbatim   
   >   
   > \param[in] LDA   
   > \verbatim   
   >          LDA is INTEGER   
   >          The leading dimension of the array A.  LDA >= max(N,1).   
   > \endverbatim   
   >   
   > \param[out] WORK   
   > \verbatim   
   >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),   
   >          where LWORK >= N when NORM = 'I'; otherwise, WORK is not   
   >          referenced.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup doubleOTHERauxiliary   

    ===================================================================== */
doublereal igraphdlanhs_(char *norm, integer *n, doublereal *a, integer *lda, 
	doublereal *work)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j;
    doublereal sum, scale;
    extern logical igraphlsame_(char *, char *);
    doublereal value = 0.;
    extern logical igraphdisnan_(doublereal *);
    extern /* Subroutine */ int igraphdlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


   =====================================================================   


       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;

    /* Function Body */
    if (*n == 0) {
	value = 0.;
    } else if (igraphlsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	value = 0.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	    i__3 = *n, i__4 = j + 1;
	    i__2 = min(i__3,i__4);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
		if (value < sum || igraphdisnan_(&sum)) {
		    value = sum;
		}
/* L10: */
	    }
/* L20: */
	}
    } else if (igraphlsame_(norm, "O") || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

	value = 0.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = 0.;
/* Computing MIN */
	    i__3 = *n, i__4 = j + 1;
	    i__2 = min(i__3,i__4);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
/* L30: */
	    }
	    if (value < sum || igraphdisnan_(&sum)) {
		value = sum;
	    }
/* L40: */
	}
    } else if (igraphlsame_(norm, "I")) {

/*        Find normI(A). */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[i__] = 0.;
/* L50: */
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	    i__3 = *n, i__4 = j + 1;
	    i__2 = min(i__3,i__4);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
/* L60: */
	    }
/* L70: */
	}
	value = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sum = work[i__];
	    if (value < sum || igraphdisnan_(&sum)) {
		value = sum;
	    }
/* L80: */
	}
    } else if (igraphlsame_(norm, "F") || igraphlsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	    i__3 = *n, i__4 = j + 1;
	    i__2 = min(i__3,i__4);
	    igraphdlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
/* L90: */
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANHS */

} /* igraphdlanhs_ */

