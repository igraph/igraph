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
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static doublereal c_b22 = -1.;
static doublereal c_b23 = 1.;

/* > \brief \b DSYTRD   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DSYTRD + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrd.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrd.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrd.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )   

         CHARACTER          UPLO   
         INTEGER            INFO, LDA, LWORK, N   
         DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ),   
        $                   WORK( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DSYTRD reduces a real symmetric matrix A to real symmetric   
   > tridiagonal form T by an orthogonal similarity transformation:   
   > Q**T * A * Q = T.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] UPLO   
   > \verbatim   
   >          UPLO is CHARACTER*1   
   >          = 'U':  Upper triangle of A is stored;   
   >          = 'L':  Lower triangle of A is stored.   
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
   >          N-by-N upper triangular part of A contains the upper   
   >          triangular part of the matrix A, and the strictly lower   
   >          triangular part of A is not referenced.  If UPLO = 'L', the   
   >          leading N-by-N lower triangular part of A contains the lower   
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
   > \param[out] WORK   
   > \verbatim   
   >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))   
   >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   
   > \endverbatim   
   >   
   > \param[in] LWORK   
   > \verbatim   
   >          LWORK is INTEGER   
   >          The dimension of the array WORK.  LWORK >= 1.   
   >          For optimum performance LWORK >= N*NB, where NB is the   
   >          optimal blocksize.   
   >   
   >          If LWORK = -1, then a workspace query is assumed; the routine   
   >          only calculates the optimal size of the WORK array, returns   
   >          this value as the first entry of the WORK array, and no error   
   >          message related to LWORK is issued by XERBLA.   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
   >          = 0:  successful exit   
   >          < 0:  if INFO = -i, the i-th argument had an illegal value   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date November 2011   

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
   Subroutine */ int igraphdsytrd_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tau, doublereal *
	work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__, j, nb, kk, nx, iws;
    extern logical igraphlsame_(char *, char *);
    integer nbmin, iinfo;
    logical upper;
    extern /* Subroutine */ int igraphdsytd2_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *), igraphdsyr2k_(char *, char *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *), igraphdlatrd_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *), igraphxerbla_(char *, 
	    integer *, ftnlen);
    extern integer igraphilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    integer ldwork, lwkopt;
    logical lquery;


/*  -- LAPACK computational routine (version 3.4.0) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       November 2011   


    =====================================================================   


       Test the input parameters   

       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    upper = igraphlsame_(uplo, "U");
    lquery = *lwork == -1;
    if (! upper && ! igraphlsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*lwork < 1 && ! lquery) {
	*info = -9;
    }

    if (*info == 0) {

/*        Determine the block size. */

	nb = igraphilaenv_(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
	lwkopt = *n * nb;
	work[1] = (doublereal) lwkopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	igraphxerbla_("DSYTRD", &i__1, (ftnlen)6);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	work[1] = 1.;
	return 0;
    }

    nx = *n;
    iws = 1;
    if (nb > 1 && nb < *n) {

/*        Determine when to cross over from blocked to unblocked code   
          (last block is always handled by unblocked code).   

   Computing MAX */
	i__1 = nb, i__2 = igraphilaenv_(&c__3, "DSYTRD", uplo, n, &c_n1, &c_n1, &
		c_n1, (ftnlen)6, (ftnlen)1);
	nx = max(i__1,i__2);
	if (nx < *n) {

/*           Determine if workspace is large enough for blocked code. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  determine the   
                minimum value of NB, and reduce NB or force use of   
                unblocked code by setting NX = N.   

   Computing MAX */
		i__1 = *lwork / ldwork;
		nb = max(i__1,1);
		nbmin = igraphilaenv_(&c__2, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1,
			 (ftnlen)6, (ftnlen)1);
		if (nb < nbmin) {
		    nx = *n;
		}
	    }
	} else {
	    nx = *n;
	}
    } else {
	nb = 1;
    }

    if (upper) {

/*        Reduce the upper triangle of A.   
          Columns 1:kk are handled by the unblocked method. */

	kk = *n - (*n - nx + nb - 1) / nb * nb;
	i__1 = kk + 1;
	i__2 = -nb;
	for (i__ = *n - nb + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
		i__2) {

/*           Reduce columns i:i+nb-1 to tridiagonal form and form the   
             matrix W which is needed to update the unreduced part of   
             the matrix */

	    i__3 = i__ + nb - 1;
	    igraphdlatrd_(uplo, &i__3, &nb, &a[a_offset], lda, &e[1], &tau[1], &
		    work[1], &ldwork);

/*           Update the unreduced submatrix A(1:i-1,1:i-1), using an   
             update of the form:  A := A - V*W**T - W*V**T */

	    i__3 = i__ - 1;
	    igraphdsyr2k_(uplo, "No transpose", &i__3, &nb, &c_b22, &a[i__ * a_dim1 
		    + 1], lda, &work[1], &ldwork, &c_b23, &a[a_offset], lda);

/*           Copy superdiagonal elements back into A, and diagonal   
             elements into D */

	    i__3 = i__ + nb - 1;
	    for (j = i__; j <= i__3; ++j) {
		a[j - 1 + j * a_dim1] = e[j - 1];
		d__[j] = a[j + j * a_dim1];
/* L10: */
	    }
/* L20: */
	}

/*        Use unblocked code to reduce the last or only block */

	igraphdsytd2_(uplo, &kk, &a[a_offset], lda, &d__[1], &e[1], &tau[1], &iinfo);
    } else {

/*        Reduce the lower triangle of A */

	i__2 = *n - nx;
	i__1 = nb;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {

/*           Reduce columns i:i+nb-1 to tridiagonal form and form the   
             matrix W which is needed to update the unreduced part of   
             the matrix */

	    i__3 = *n - i__ + 1;
	    igraphdlatrd_(uplo, &i__3, &nb, &a[i__ + i__ * a_dim1], lda, &e[i__], &
		    tau[i__], &work[1], &ldwork);

/*           Update the unreduced submatrix A(i+ib:n,i+ib:n), using   
             an update of the form:  A := A - V*W**T - W*V**T */

	    i__3 = *n - i__ - nb + 1;
	    igraphdsyr2k_(uplo, "No transpose", &i__3, &nb, &c_b22, &a[i__ + nb + 
		    i__ * a_dim1], lda, &work[nb + 1], &ldwork, &c_b23, &a[
		    i__ + nb + (i__ + nb) * a_dim1], lda);

/*           Copy subdiagonal elements back into A, and diagonal   
             elements into D */

	    i__3 = i__ + nb - 1;
	    for (j = i__; j <= i__3; ++j) {
		a[j + 1 + j * a_dim1] = e[j];
		d__[j] = a[j + j * a_dim1];
/* L30: */
	    }
/* L40: */
	}

/*        Use unblocked code to reduce the last or only block */

	i__1 = *n - i__ + 1;
	igraphdsytd2_(uplo, &i__1, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[i__], 
		&tau[i__], &iinfo);
    }

    work[1] = (doublereal) lwkopt;
    return 0;

/*     End of DSYTRD */

} /* igraphdsytrd_ */

