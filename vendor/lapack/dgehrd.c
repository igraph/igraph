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
static integer c__65 = 65;
static doublereal c_b25 = -1.;
static doublereal c_b26 = 1.;

/* > \brief \b DGEHRD   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DGEHRD + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgehrd.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgehrd.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgehrd.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )   

         INTEGER            IHI, ILO, INFO, LDA, LWORK, N   
         DOUBLE PRECISION  A( LDA, * ), TAU( * ), WORK( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DGEHRD reduces a real general matrix A to upper Hessenberg form H by   
   > an orthogonal similarity transformation:  Q**T * A * Q = H .   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix A.  N >= 0.   
   > \endverbatim   
   >   
   > \param[in] ILO   
   > \verbatim   
   >          ILO is INTEGER   
   > \endverbatim   
   >   
   > \param[in] IHI   
   > \verbatim   
   >          IHI is INTEGER   
   >   
   >          It is assumed that A is already upper triangular in rows   
   >          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally   
   >          set by a previous call to DGEBAL; otherwise they should be   
   >          set to 1 and N respectively. See Further Details.   
   >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   
   > \endverbatim   
   >   
   > \param[in,out] A   
   > \verbatim   
   >          A is DOUBLE PRECISION array, dimension (LDA,N)   
   >          On entry, the N-by-N general matrix to be reduced.   
   >          On exit, the upper triangle and the first subdiagonal of A   
   >          are overwritten with the upper Hessenberg matrix H, and the   
   >          elements below the first subdiagonal, with the array TAU,   
   >          represent the orthogonal matrix Q as a product of elementary   
   >          reflectors. See Further Details.   
   > \endverbatim   
   >   
   > \param[in] LDA   
   > \verbatim   
   >          LDA is INTEGER   
   >          The leading dimension of the array A.  LDA >= max(1,N).   
   > \endverbatim   
   >   
   > \param[out] TAU   
   > \verbatim   
   >          TAU is DOUBLE PRECISION array, dimension (N-1)   
   >          The scalar factors of the elementary reflectors (see Further   
   >          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to   
   >          zero.   
   > \endverbatim   
   >   
   > \param[out] WORK   
   > \verbatim   
   >          WORK is DOUBLE PRECISION array, dimension (LWORK)   
   >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   
   > \endverbatim   
   >   
   > \param[in] LWORK   
   > \verbatim   
   >          LWORK is INTEGER   
   >          The length of the array WORK.  LWORK >= max(1,N).   
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
   >          < 0:  if INFO = -i, the i-th argument had an illegal value.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date November 2011   

   > \ingroup doubleGEcomputational   

   > \par Further Details:   
    =====================   
   >   
   > \verbatim   
   >   
   >  The matrix Q is represented as a product of (ihi-ilo) elementary   
   >  reflectors   
   >   
   >     Q = H(ilo) H(ilo+1) . . . H(ihi-1).   
   >   
   >  Each H(i) has the form   
   >   
   >     H(i) = I - tau * v * v**T   
   >   
   >  where tau is a real scalar, and v is a real vector with   
   >  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on   
   >  exit in A(i+2:ihi,i), and tau in TAU(i).   
   >   
   >  The contents of A are illustrated by the following example, with   
   >  n = 7, ilo = 2 and ihi = 6:   
   >   
   >  on entry,                        on exit,   
   >   
   >  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )   
   >  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )   
   >  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )   
   >  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )   
   >  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )   
   >  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )   
   >  (                         a )    (                          a )   
   >   
   >  where a denotes an element of the original matrix A, h denotes a   
   >  modified element of the upper Hessenberg matrix H, and vi denotes an   
   >  element of the vector defining H(i).   
   >   
   >  This file is a slight modification of LAPACK-3.0's DGEHRD   
   >  subroutine incorporating improvements proposed by Quintana-Orti and   
   >  Van de Geijn (2006). (See DLAHR2.)   
   > \endverbatim   
   >   
    =====================================================================   
   Subroutine */ int igraphdgehrd_(integer *n, integer *ilo, integer *ihi, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    integer i__, j;
    doublereal t[4160]	/* was [65][64] */;
    integer ib;
    doublereal ei;
    integer nb, nh, nx, iws;
    extern /* Subroutine */ int igraphdgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);
    integer nbmin, iinfo;
    extern /* Subroutine */ int igraphdtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *), igraphdaxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), igraphdgehd2_(integer *, integer *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *), igraphdlahr2_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *), 
	    igraphdlarfb_(char *, char *, char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *), igraphxerbla_(char *, integer *, ftnlen);
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
    --tau;
    --work;

    /* Function Body */
    *info = 0;
/* Computing MIN */
    i__1 = 64, i__2 = igraphilaenv_(&c__1, "DGEHRD", " ", n, ilo, ihi, &c_n1, (
	    ftnlen)6, (ftnlen)1);
    nb = min(i__1,i__2);
    lwkopt = *n * nb;
    work[1] = (doublereal) lwkopt;
    lquery = *lwork == -1;
    if (*n < 0) {
	*info = -1;
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
	*info = -2;
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*lwork < max(1,*n) && ! lquery) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	igraphxerbla_("DGEHRD", &i__1, (ftnlen)6);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero */

    i__1 = *ilo - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tau[i__] = 0.;
/* L10: */
    }
    i__1 = *n - 1;
    for (i__ = max(1,*ihi); i__ <= i__1; ++i__) {
	tau[i__] = 0.;
/* L20: */
    }

/*     Quick return if possible */

    nh = *ihi - *ilo + 1;
    if (nh <= 1) {
	work[1] = 1.;
	return 0;
    }

/*     Determine the block size   

   Computing MIN */
    i__1 = 64, i__2 = igraphilaenv_(&c__1, "DGEHRD", " ", n, ilo, ihi, &c_n1, (
	    ftnlen)6, (ftnlen)1);
    nb = min(i__1,i__2);
    nbmin = 2;
    iws = 1;
    if (nb > 1 && nb < nh) {

/*        Determine when to cross over from blocked to unblocked code   
          (last block is always handled by unblocked code)   

   Computing MAX */
	i__1 = nb, i__2 = igraphilaenv_(&c__3, "DGEHRD", " ", n, ilo, ihi, &c_n1, (
		ftnlen)6, (ftnlen)1);
	nx = max(i__1,i__2);
	if (nx < nh) {

/*           Determine if workspace is large enough for blocked code */

	    iws = *n * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  determine the   
                minimum value of NB, and reduce NB or force use of   
                unblocked code   

   Computing MAX */
		i__1 = 2, i__2 = igraphilaenv_(&c__2, "DGEHRD", " ", n, ilo, ihi, &
			c_n1, (ftnlen)6, (ftnlen)1);
		nbmin = max(i__1,i__2);
		if (*lwork >= *n * nbmin) {
		    nb = *lwork / *n;
		} else {
		    nb = 1;
		}
	    }
	}
    }
    ldwork = *n;

    if (nb < nbmin || nb >= nh) {

/*        Use unblocked code below */

	i__ = *ilo;

    } else {

/*        Use blocked code */

	i__1 = *ihi - 1 - nx;
	i__2 = nb;
	for (i__ = *ilo; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = nb, i__4 = *ihi - i__;
	    ib = min(i__3,i__4);

/*           Reduce columns i:i+ib-1 to Hessenberg form, returning the   
             matrices V and T of the block reflector H = I - V*T*V**T   
             which performs the reduction, and also the matrix Y = A*V*T */

	    igraphdlahr2_(ihi, &i__, &ib, &a[i__ * a_dim1 + 1], lda, &tau[i__], t, &
		    c__65, &work[1], &ldwork);

/*           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the   
             right, computing  A := A - Y * V**T. V(i+ib,ib-1) must be set   
             to 1 */

	    ei = a[i__ + ib + (i__ + ib - 1) * a_dim1];
	    a[i__ + ib + (i__ + ib - 1) * a_dim1] = 1.;
	    i__3 = *ihi - i__ - ib + 1;
	    igraphdgemm_("No transpose", "Transpose", ihi, &i__3, &ib, &c_b25, &
		    work[1], &ldwork, &a[i__ + ib + i__ * a_dim1], lda, &
		    c_b26, &a[(i__ + ib) * a_dim1 + 1], lda);
	    a[i__ + ib + (i__ + ib - 1) * a_dim1] = ei;

/*           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the   
             right */

	    i__3 = ib - 1;
	    igraphdtrmm_("Right", "Lower", "Transpose", "Unit", &i__, &i__3, &c_b26,
		     &a[i__ + 1 + i__ * a_dim1], lda, &work[1], &ldwork);
	    i__3 = ib - 2;
	    for (j = 0; j <= i__3; ++j) {
		igraphdaxpy_(&i__, &c_b25, &work[ldwork * j + 1], &c__1, &a[(i__ + 
			j + 1) * a_dim1 + 1], &c__1);
/* L30: */
	    }

/*           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the   
             left */

	    i__3 = *ihi - i__;
	    i__4 = *n - i__ - ib + 1;
	    igraphdlarfb_("Left", "Transpose", "Forward", "Columnwise", &i__3, &
		    i__4, &ib, &a[i__ + 1 + i__ * a_dim1], lda, t, &c__65, &a[
		    i__ + 1 + (i__ + ib) * a_dim1], lda, &work[1], &ldwork);
/* L40: */
	}
    }

/*     Use unblocked code to reduce the rest of the matrix */

    igraphdgehd2_(n, &i__, ihi, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
    work[1] = (doublereal) iws;

    return 0;

/*     End of DGEHRD */

} /* igraphdgehrd_ */

