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
static doublereal c_b8 = 0.;
static doublereal c_b14 = -1.;

/* > \brief \b DSYTD2 reduces a symmetric matrix to real symmetric tridiagonal form by an orthogonal similarit
y transformation (unblocked algorithm).   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DSYTD2 + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytd2.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytd2.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytd2.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )   

         CHARACTER          UPLO   
         INTEGER            INFO, LDA, N   
         DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal   
   > form T by an orthogonal similarity transformation: Q**T * A * Q = T.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] UPLO   
   > \verbatim   
   >          UPLO is CHARACTER*1   
   >          Specifies whether the upper or lower triangular part of the   
   >          symmetric matrix A is stored:   
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
   >          n-by-n upper triangular part of A contains the upper   
   >          triangular part of the matrix A, and the strictly lower   
   >          triangular part of A is not referenced.  If UPLO = 'L', the   
   >          leading n-by-n lower triangular part of A contains the lower   
   >          triangular part of the matrix A, and the strictly upper   
   >          triangular part of A is not referenced.   
   >          On exit, if UPLO = 'U', the diagonal and first superdiagonal   
   >          of A are overwritten by the corresponding elements of the   
   >          tridiagonal matrix T, and the elements above the first   
   >          superdiagonal, with the array TAU, represent the orthogonal   
   >          matrix Q as a product of elementary reflectors; if UPLO   
   >          = 'L', the diagonal and first subdiagonal of A are over-   
   >          written by the corresponding elements of the tridiagonal   
   >          matrix T, and the elements below the first subdiagonal, with   
   >          the array TAU, represent the orthogonal matrix Q as a product   
   >          of elementary reflectors. See Further Details.   
   > \endverbatim   
   >   
   > \param[in] LDA   
   > \verbatim   
   >          LDA is INTEGER   
   >          The leading dimension of the array A.  LDA >= max(1,N).   
   > \endverbatim   
   >   
   > \param[out] D   
   > \verbatim   
   >          D is DOUBLE PRECISION array, dimension (N)   
   >          The diagonal elements of the tridiagonal matrix T:   
   >          D(i) = A(i,i).   
   > \endverbatim   
   >   
   > \param[out] E   
   > \verbatim   
   >          E is DOUBLE PRECISION array, dimension (N-1)   
   >          The off-diagonal elements of the tridiagonal matrix T:   
   >          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.   
   > \endverbatim   
   >   
   > \param[out] TAU   
   > \verbatim   
   >          TAU is DOUBLE PRECISION array, dimension (N-1)   
   >          The scalar factors of the elementary reflectors (see Further   
   >          Details).   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
   >          = 0:  successful exit   
   >          < 0:  if INFO = -i, the i-th argument had an illegal value.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup doubleSYcomputational   

   > \par Further Details:   
    =====================   
   >   
   > \verbatim   
   >   
   >  If UPLO = 'U', the matrix Q is represented as a product of elementary   
   >  reflectors   
   >   
   >     Q = H(n-1) . . . H(2) H(1).   
   >   
   >  Each H(i) has the form   
   >   
   >     H(i) = I - tau * v * v**T   
   >   
   >  where tau is a real scalar, and v is a real vector with   
   >  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in   
   >  A(1:i-1,i+1), and tau in TAU(i).   
   >   
   >  If UPLO = 'L', the matrix Q is represented as a product of elementary   
   >  reflectors   
   >   
   >     Q = H(1) H(2) . . . H(n-1).   
   >   
   >  Each H(i) has the form   
   >   
   >     H(i) = I - tau * v * v**T   
   >   
   >  where tau is a real scalar, and v is a real vector with   
   >  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),   
   >  and tau in TAU(i).   
   >   
   >  The contents of A on exit are illustrated by the following examples   
   >  with n = 5:   
   >   
   >  if UPLO = 'U':                       if UPLO = 'L':   
   >   
   >    (  d   e   v2  v3  v4 )              (  d                  )   
   >    (      d   e   v3  v4 )              (  e   d              )   
   >    (          d   e   v4 )              (  v1  e   d          )   
   >    (              d   e  )              (  v1  v2  e   d      )   
   >    (                  d  )              (  v1  v2  v3  e   d  )   
   >   
   >  where d and e denote diagonal and off-diagonal elements of T, and vi   
   >  denotes an element of the vector defining H(i).   
   > \endverbatim   
   >   
    =====================================================================   
   Subroutine */ int igraphdsytd2_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tau, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__;
    extern doublereal igraphddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    doublereal taui;
    extern /* Subroutine */ int igraphdsyr2_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    doublereal alpha;
    extern logical igraphlsame_(char *, char *);
    extern /* Subroutine */ int igraphdaxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    logical upper;
    extern /* Subroutine */ int igraphdsymv_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *), igraphdlarfg_(integer *, doublereal *,
	     doublereal *, integer *, doublereal *), igraphxerbla_(char *, integer *
	    , ftnlen);


/*  -- LAPACK computational routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       Test the input parameters   

       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;

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
	igraphxerbla_("DSYTD2", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return 0;
    }

    if (upper) {

/*        Reduce the upper triangle of A */

	for (i__ = *n - 1; i__ >= 1; --i__) {

/*           Generate elementary reflector H(i) = I - tau * v * v**T   
             to annihilate A(1:i-1,i+1) */

	    igraphdlarfg_(&i__, &a[i__ + (i__ + 1) * a_dim1], &a[(i__ + 1) * a_dim1 
		    + 1], &c__1, &taui);
	    e[i__] = a[i__ + (i__ + 1) * a_dim1];

	    if (taui != 0.) {

/*              Apply H(i) from both sides to A(1:i,1:i) */

		a[i__ + (i__ + 1) * a_dim1] = 1.;

/*              Compute  x := tau * A * v  storing x in TAU(1:i) */

		igraphdsymv_(uplo, &i__, &taui, &a[a_offset], lda, &a[(i__ + 1) * 
			a_dim1 + 1], &c__1, &c_b8, &tau[1], &c__1);

/*              Compute  w := x - 1/2 * tau * (x**T * v) * v */

		alpha = taui * -.5 * igraphddot_(&i__, &tau[1], &c__1, &a[(i__ + 1) 
			* a_dim1 + 1], &c__1);
		igraphdaxpy_(&i__, &alpha, &a[(i__ + 1) * a_dim1 + 1], &c__1, &tau[
			1], &c__1);

/*              Apply the transformation as a rank-2 update:   
                   A := A - v * w**T - w * v**T */

		igraphdsyr2_(uplo, &i__, &c_b14, &a[(i__ + 1) * a_dim1 + 1], &c__1, 
			&tau[1], &c__1, &a[a_offset], lda);

		a[i__ + (i__ + 1) * a_dim1] = e[i__];
	    }
	    d__[i__ + 1] = a[i__ + 1 + (i__ + 1) * a_dim1];
	    tau[i__] = taui;
/* L10: */
	}
	d__[1] = a[a_dim1 + 1];
    } else {

/*        Reduce the lower triangle of A */

	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Generate elementary reflector H(i) = I - tau * v * v**T   
             to annihilate A(i+2:n,i) */

	    i__2 = *n - i__;
/* Computing MIN */
	    i__3 = i__ + 2;
	    igraphdlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3,*n) + i__ *
		     a_dim1], &c__1, &taui);
	    e[i__] = a[i__ + 1 + i__ * a_dim1];

	    if (taui != 0.) {

/*              Apply H(i) from both sides to A(i+1:n,i+1:n) */

		a[i__ + 1 + i__ * a_dim1] = 1.;

/*              Compute  x := tau * A * v  storing y in TAU(i:n-1) */

		i__2 = *n - i__;
		igraphdsymv_(uplo, &i__2, &taui, &a[i__ + 1 + (i__ + 1) * a_dim1], 
			lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b8, &tau[
			i__], &c__1);

/*              Compute  w := x - 1/2 * tau * (x**T * v) * v */

		i__2 = *n - i__;
		alpha = taui * -.5 * igraphddot_(&i__2, &tau[i__], &c__1, &a[i__ + 
			1 + i__ * a_dim1], &c__1);
		i__2 = *n - i__;
		igraphdaxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[
			i__], &c__1);

/*              Apply the transformation as a rank-2 update:   
                   A := A - v * w**T - w * v**T */

		i__2 = *n - i__;
		igraphdsyr2_(uplo, &i__2, &c_b14, &a[i__ + 1 + i__ * a_dim1], &c__1,
			 &tau[i__], &c__1, &a[i__ + 1 + (i__ + 1) * a_dim1], 
			lda);

		a[i__ + 1 + i__ * a_dim1] = e[i__];
	    }
	    d__[i__] = a[i__ + i__ * a_dim1];
	    tau[i__] = taui;
/* L20: */
	}
	d__[*n] = a[*n + *n * a_dim1];
    }

    return 0;

/*     End of DSYTD2 */

} /* igraphdsytd2_ */

