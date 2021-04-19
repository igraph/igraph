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
static integer c_n1 = -1;

/* > \brief \b DORGHR   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DORGHR + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorghr.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorghr.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorghr.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )   

         INTEGER            IHI, ILO, INFO, LDA, LWORK, N   
         DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DORGHR generates a real orthogonal matrix Q which is defined as the   
   > product of IHI-ILO elementary reflectors of order N, as returned by   
   > DGEHRD:   
   >   
   > Q = H(ilo) H(ilo+1) . . . H(ihi-1).   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix Q. N >= 0.   
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
   >          ILO and IHI must have the same values as in the previous call   
   >          of DGEHRD. Q is equal to the unit matrix except in the   
   >          submatrix Q(ilo+1:ihi,ilo+1:ihi).   
   >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   
   > \endverbatim   
   >   
   > \param[in,out] A   
   > \verbatim   
   >          A is DOUBLE PRECISION array, dimension (LDA,N)   
   >          On entry, the vectors which define the elementary reflectors,   
   >          as returned by DGEHRD.   
   >          On exit, the N-by-N orthogonal matrix Q.   
   > \endverbatim   
   >   
   > \param[in] LDA   
   > \verbatim   
   >          LDA is INTEGER   
   >          The leading dimension of the array A. LDA >= max(1,N).   
   > \endverbatim   
   >   
   > \param[in] TAU   
   > \verbatim   
   >          TAU is DOUBLE PRECISION array, dimension (N-1)   
   >          TAU(i) must contain the scalar factor of the elementary   
   >          reflector H(i), as returned by DGEHRD.   
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
   >          The dimension of the array WORK. LWORK >= IHI-ILO.   
   >          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is   
   >          the optimal blocksize.   
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

   > \ingroup doubleOTHERcomputational   

    =====================================================================   
   Subroutine */ int igraphdorghr_(integer *n, integer *ilo, integer *ihi, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, nb, nh, iinfo;
    extern /* Subroutine */ int igraphxerbla_(char *, integer *, ftnlen);
    extern integer igraphilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int igraphdorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    integer lwkopt;
    logical lquery;


/*  -- LAPACK computational routine (version 3.4.0) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       November 2011   


    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    nh = *ihi - *ilo;
    lquery = *lwork == -1;
    if (*n < 0) {
	*info = -1;
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
	*info = -2;
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*lwork < max(1,nh) && ! lquery) {
	*info = -8;
    }

    if (*info == 0) {
	nb = igraphilaenv_(&c__1, "DORGQR", " ", &nh, &nh, &nh, &c_n1, (ftnlen)6, (
		ftnlen)1);
	lwkopt = max(1,nh) * nb;
	work[1] = (doublereal) lwkopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	igraphxerbla_("DORGHR", &i__1, (ftnlen)6);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	work[1] = 1.;
	return 0;
    }

/*     Shift the vectors which define the elementary reflectors one   
       column to the right, and set the first ilo and the last n-ihi   
       rows and columns to those of the unit matrix */

    i__1 = *ilo + 1;
    for (j = *ihi; j >= i__1; --j) {
	i__2 = j - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    a[i__ + j * a_dim1] = 0.;
/* L10: */
	}
	i__2 = *ihi;
	for (i__ = j + 1; i__ <= i__2; ++i__) {
	    a[i__ + j * a_dim1] = a[i__ + (j - 1) * a_dim1];
/* L20: */
	}
	i__2 = *n;
	for (i__ = *ihi + 1; i__ <= i__2; ++i__) {
	    a[i__ + j * a_dim1] = 0.;
/* L30: */
	}
/* L40: */
    }
    i__1 = *ilo;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    a[i__ + j * a_dim1] = 0.;
/* L50: */
	}
	a[j + j * a_dim1] = 1.;
/* L60: */
    }
    i__1 = *n;
    for (j = *ihi + 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    a[i__ + j * a_dim1] = 0.;
/* L70: */
	}
	a[j + j * a_dim1] = 1.;
/* L80: */
    }

    if (nh > 0) {

/*        Generate Q(ilo+1:ihi,ilo+1:ihi) */

	igraphdorgqr_(&nh, &nh, &nh, &a[*ilo + 1 + (*ilo + 1) * a_dim1], lda, &tau[*
		ilo], &work[1], lwork, &iinfo);
    }
    work[1] = (doublereal) lwkopt;
    return 0;

/*     End of DORGHR */

} /* igraphdorghr_ */

