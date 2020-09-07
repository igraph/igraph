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
static doublereal c_b10 = -1.;
static doublereal c_b12 = 1.;

/* > \brief \b DPOTF2 computes the Cholesky factorization of a symmetric/Hermitian positive definite matrix (u
nblocked algorithm).   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DPOTF2 + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpotf2.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpotf2.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpotf2.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO )   

         CHARACTER          UPLO   
         INTEGER            INFO, LDA, N   
         DOUBLE PRECISION   A( LDA, * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DPOTF2 computes the Cholesky factorization of a real symmetric   
   > positive definite matrix A.   
   >   
   > The factorization has the form   
   >    A = U**T * U ,  if UPLO = 'U', or   
   >    A = L  * L**T,  if UPLO = 'L',   
   > where U is an upper triangular matrix and L is lower triangular.   
   >   
   > This is the unblocked version of the algorithm, calling Level 2 BLAS.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] UPLO   
   > \verbatim   
   >          UPLO is CHARACTER*1   
   >          Specifies whether the upper or lower triangular part of the   
   >          symmetric matrix A is stored.   
   >          = 'U':  Upper triangular   
   >          = 'L':  Lower triangular   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix A.  N >= 0.   
   > \endverbatim   
   >   
   > \param[in,out] A   
   > \verbatim   
   >          A is DOUBLE PRECISION array, dimension (LDA,N)   
   >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
   >          n by n upper triangular part of A contains the upper   
   >          triangular part of the matrix A, and the strictly lower   
   >          triangular part of A is not referenced.  If UPLO = 'L', the   
   >          leading n by n lower triangular part of A contains the lower   
   >          triangular part of the matrix A, and the strictly upper   
   >          triangular part of A is not referenced.   
   >   
   >          On exit, if INFO = 0, the factor U or L from the Cholesky   
   >          factorization A = U**T *U  or A = L*L**T.   
   > \endverbatim   
   >   
   > \param[in] LDA   
   > \verbatim   
   >          LDA is INTEGER   
   >          The leading dimension of the array A.  LDA >= max(1,N).   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
   >          = 0: successful exit   
   >          < 0: if INFO = -k, the k-th argument had an illegal value   
   >          > 0: if INFO = k, the leading minor of order k is not   
   >               positive definite, and the factorization could not be   
   >               completed.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup doublePOcomputational   

    =====================================================================   
   Subroutine */ int igraphdpotf2_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer j;
    doublereal ajj;
    extern doublereal igraphddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int igraphdscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical igraphlsame_(char *, char *);
    extern /* Subroutine */ int igraphdgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    logical upper;
    extern logical igraphdisnan_(doublereal *);
    extern /* Subroutine */ int igraphxerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    upper = igraphlsame_(uplo, "U");
    if (! upper && ! igraphlsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	igraphxerbla_("DPOTF2", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Compute the Cholesky factorization A = U**T *U. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

/*           Compute U(J,J) and test for non-positive-definiteness. */

	    i__2 = j - 1;
	    ajj = a[j + j * a_dim1] - igraphddot_(&i__2, &a[j * a_dim1 + 1], &c__1, 
		    &a[j * a_dim1 + 1], &c__1);
	    if (ajj <= 0. || igraphdisnan_(&ajj)) {
		a[j + j * a_dim1] = ajj;
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    a[j + j * a_dim1] = ajj;

/*           Compute elements J+1:N of row J. */

	    if (j < *n) {
		i__2 = j - 1;
		i__3 = *n - j;
		igraphdgemv_("Transpose", &i__2, &i__3, &c_b10, &a[(j + 1) * a_dim1 
			+ 1], lda, &a[j * a_dim1 + 1], &c__1, &c_b12, &a[j + (
			j + 1) * a_dim1], lda);
		i__2 = *n - j;
		d__1 = 1. / ajj;
		igraphdscal_(&i__2, &d__1, &a[j + (j + 1) * a_dim1], lda);
	    }
/* L10: */
	}
    } else {

/*        Compute the Cholesky factorization A = L*L**T. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness. */

	    i__2 = j - 1;
	    ajj = a[j + j * a_dim1] - igraphddot_(&i__2, &a[j + a_dim1], lda, &a[j 
		    + a_dim1], lda);
	    if (ajj <= 0. || igraphdisnan_(&ajj)) {
		a[j + j * a_dim1] = ajj;
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    a[j + j * a_dim1] = ajj;

/*           Compute elements J+1:N of column J. */

	    if (j < *n) {
		i__2 = *n - j;
		i__3 = j - 1;
		igraphdgemv_("No transpose", &i__2, &i__3, &c_b10, &a[j + 1 + 
			a_dim1], lda, &a[j + a_dim1], lda, &c_b12, &a[j + 1 + 
			j * a_dim1], &c__1);
		i__2 = *n - j;
		d__1 = 1. / ajj;
		igraphdscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
	    }
/* L20: */
	}
    }
    goto L40;

L30:
    *info = j;

L40:
    return 0;

/*     End of DPOTF2 */

} /* igraphdpotf2_ */

