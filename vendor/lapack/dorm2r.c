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

/* > \brief \b DORM2R multiplies a general matrix by the orthogonal matrix from a QR factorization determined 
by sgeqrf (unblocked algorithm).   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DORM2R + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorm2r.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorm2r.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorm2r.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,   
                            WORK, INFO )   

         CHARACTER          SIDE, TRANS   
         INTEGER            INFO, K, LDA, LDC, M, N   
         DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DORM2R overwrites the general real m by n matrix C with   
   >   
   >       Q * C  if SIDE = 'L' and TRANS = 'N', or   
   >   
   >       Q**T* C  if SIDE = 'L' and TRANS = 'T', or   
   >   
   >       C * Q  if SIDE = 'R' and TRANS = 'N', or   
   >   
   >       C * Q**T if SIDE = 'R' and TRANS = 'T',   
   >   
   > where Q is a real orthogonal matrix defined as the product of k   
   > elementary reflectors   
   >   
   >       Q = H(1) H(2) . . . H(k)   
   >   
   > as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n   
   > if SIDE = 'R'.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] SIDE   
   > \verbatim   
   >          SIDE is CHARACTER*1   
   >          = 'L': apply Q or Q**T from the Left   
   >          = 'R': apply Q or Q**T from the Right   
   > \endverbatim   
   >   
   > \param[in] TRANS   
   > \verbatim   
   >          TRANS is CHARACTER*1   
   >          = 'N': apply Q  (No transpose)   
   >          = 'T': apply Q**T (Transpose)   
   > \endverbatim   
   >   
   > \param[in] M   
   > \verbatim   
   >          M is INTEGER   
   >          The number of rows of the matrix C. M >= 0.   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The number of columns of the matrix C. N >= 0.   
   > \endverbatim   
   >   
   > \param[in] K   
   > \verbatim   
   >          K is INTEGER   
   >          The number of elementary reflectors whose product defines   
   >          the matrix Q.   
   >          If SIDE = 'L', M >= K >= 0;   
   >          if SIDE = 'R', N >= K >= 0.   
   > \endverbatim   
   >   
   > \param[in] A   
   > \verbatim   
   >          A is DOUBLE PRECISION array, dimension (LDA,K)   
   >          The i-th column must contain the vector which defines the   
   >          elementary reflector H(i), for i = 1,2,...,k, as returned by   
   >          DGEQRF in the first k columns of its array argument A.   
   >          A is modified by the routine but restored on exit.   
   > \endverbatim   
   >   
   > \param[in] LDA   
   > \verbatim   
   >          LDA is INTEGER   
   >          The leading dimension of the array A.   
   >          If SIDE = 'L', LDA >= max(1,M);   
   >          if SIDE = 'R', LDA >= max(1,N).   
   > \endverbatim   
   >   
   > \param[in] TAU   
   > \verbatim   
   >          TAU is DOUBLE PRECISION array, dimension (K)   
   >          TAU(i) must contain the scalar factor of the elementary   
   >          reflector H(i), as returned by DGEQRF.   
   > \endverbatim   
   >   
   > \param[in,out] C   
   > \verbatim   
   >          C is DOUBLE PRECISION array, dimension (LDC,N)   
   >          On entry, the m by n matrix C.   
   >          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.   
   > \endverbatim   
   >   
   > \param[in] LDC   
   > \verbatim   
   >          LDC is INTEGER   
   >          The leading dimension of the array C. LDC >= max(1,M).   
   > \endverbatim   
   >   
   > \param[out] WORK   
   > \verbatim   
   >          WORK is DOUBLE PRECISION array, dimension   
   >                                   (N) if SIDE = 'L',   
   >                                   (M) if SIDE = 'R'   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
   >          = 0: successful exit   
   >          < 0: if INFO = -i, the i-th argument had an illegal value   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup doubleOTHERcomputational   

    =====================================================================   
   Subroutine */ int igraphdorm2r_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    integer i__, i1, i2, i3, ic, jc, mi, ni, nq;
    doublereal aii;
    logical left;
    extern /* Subroutine */ int igraphdlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *);
    extern logical igraphlsame_(char *, char *);
    extern /* Subroutine */ int igraphxerbla_(char *, integer *, ftnlen);
    logical notran;


/*  -- LAPACK computational routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = igraphlsame_(side, "L");
    notran = igraphlsame_(trans, "N");

/*     NQ is the order of Q */

    if (left) {
	nq = *m;
    } else {
	nq = *n;
    }
    if (! left && ! igraphlsame_(side, "R")) {
	*info = -1;
    } else if (! notran && ! igraphlsame_(trans, "T")) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > nq) {
	*info = -5;
    } else if (*lda < max(1,nq)) {
	*info = -7;
    } else if (*ldc < max(1,*m)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	igraphxerbla_("DORM2R", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0 || *k == 0) {
	return 0;
    }

    if (left && ! notran || ! left && notran) {
	i1 = 1;
	i2 = *k;
	i3 = 1;
    } else {
	i1 = *k;
	i2 = 1;
	i3 = -1;
    }

    if (left) {
	ni = *n;
	jc = 1;
    } else {
	mi = *m;
	ic = 1;
    }

    i__1 = i2;
    i__2 = i3;
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	if (left) {

/*           H(i) is applied to C(i:m,1:n) */

	    mi = *m - i__ + 1;
	    ic = i__;
	} else {

/*           H(i) is applied to C(1:m,i:n) */

	    ni = *n - i__ + 1;
	    jc = i__;
	}

/*        Apply H(i) */

	aii = a[i__ + i__ * a_dim1];
	a[i__ + i__ * a_dim1] = 1.;
	igraphdlarf_(side, &mi, &ni, &a[i__ + i__ * a_dim1], &c__1, &tau[i__], &c__[
		ic + jc * c_dim1], ldc, &work[1]);
	a[i__ + i__ * a_dim1] = aii;
/* L10: */
    }
    return 0;

/*     End of DORM2R */

} /* igraphdorm2r_ */

