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
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> DGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matr
ices</b>   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DGEEV + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeev.f
">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeev.f
">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeev.f
">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,   
                           LDVR, WORK, LWORK, INFO )   

         CHARACTER          JOBVL, JOBVR   
         INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N   
         DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),   
        $                   WI( * ), WORK( * ), WR( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DGEEV computes for an N-by-N real nonsymmetric matrix A, the   
   > eigenvalues and, optionally, the left and/or right eigenvectors.   
   >   
   > The right eigenvector v(j) of A satisfies   
   >                  A * v(j) = lambda(j) * v(j)   
   > where lambda(j) is its eigenvalue.   
   > The left eigenvector u(j) of A satisfies   
   >               u(j)**H * A = lambda(j) * u(j)**H   
   > where u(j)**H denotes the conjugate-transpose of u(j).   
   >   
   > The computed eigenvectors are normalized to have Euclidean norm   
   > equal to 1 and largest component real.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] JOBVL   
   > \verbatim   
   >          JOBVL is CHARACTER*1   
   >          = 'N': left eigenvectors of A are not computed;   
   >          = 'V': left eigenvectors of A are computed.   
   > \endverbatim   
   >   
   > \param[in] JOBVR   
   > \verbatim   
   >          JOBVR is CHARACTER*1   
   >          = 'N': right eigenvectors of A are not computed;   
   >          = 'V': right eigenvectors of A are computed.   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix A. N >= 0.   
   > \endverbatim   
   >   
   > \param[in,out] A   
   > \verbatim   
   >          A is DOUBLE PRECISION array, dimension (LDA,N)   
   >          On entry, the N-by-N matrix A.   
   >          On exit, A has been overwritten.   
   > \endverbatim   
   >   
   > \param[in] LDA   
   > \verbatim   
   >          LDA is INTEGER   
   >          The leading dimension of the array A.  LDA >= max(1,N).   
   > \endverbatim   
   >   
   > \param[out] WR   
   > \verbatim   
   >          WR is DOUBLE PRECISION array, dimension (N)   
   > \endverbatim   
   >   
   > \param[out] WI   
   > \verbatim   
   >          WI is DOUBLE PRECISION array, dimension (N)   
   >          WR and WI contain the real and imaginary parts,   
   >          respectively, of the computed eigenvalues.  Complex   
   >          conjugate pairs of eigenvalues appear consecutively   
   >          with the eigenvalue having the positive imaginary part   
   >          first.   
   > \endverbatim   
   >   
   > \param[out] VL   
   > \verbatim   
   >          VL is DOUBLE PRECISION array, dimension (LDVL,N)   
   >          If JOBVL = 'V', the left eigenvectors u(j) are stored one   
   >          after another in the columns of VL, in the same order   
   >          as their eigenvalues.   
   >          If JOBVL = 'N', VL is not referenced.   
   >          If the j-th eigenvalue is real, then u(j) = VL(:,j),   
   >          the j-th column of VL.   
   >          If the j-th and (j+1)-st eigenvalues form a complex   
   >          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and   
   >          u(j+1) = VL(:,j) - i*VL(:,j+1).   
   > \endverbatim   
   >   
   > \param[in] LDVL   
   > \verbatim   
   >          LDVL is INTEGER   
   >          The leading dimension of the array VL.  LDVL >= 1; if   
   >          JOBVL = 'V', LDVL >= N.   
   > \endverbatim   
   >   
   > \param[out] VR   
   > \verbatim   
   >          VR is DOUBLE PRECISION array, dimension (LDVR,N)   
   >          If JOBVR = 'V', the right eigenvectors v(j) are stored one   
   >          after another in the columns of VR, in the same order   
   >          as their eigenvalues.   
   >          If JOBVR = 'N', VR is not referenced.   
   >          If the j-th eigenvalue is real, then v(j) = VR(:,j),   
   >          the j-th column of VR.   
   >          If the j-th and (j+1)-st eigenvalues form a complex   
   >          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and   
   >          v(j+1) = VR(:,j) - i*VR(:,j+1).   
   > \endverbatim   
   >   
   > \param[in] LDVR   
   > \verbatim   
   >          LDVR is INTEGER   
   >          The leading dimension of the array VR.  LDVR >= 1; if   
   >          JOBVR = 'V', LDVR >= N.   
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
   >          The dimension of the array WORK.  LWORK >= max(1,3*N), and   
   >          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good   
   >          performance, LWORK must generally be larger.   
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
   >          > 0:  if INFO = i, the QR algorithm failed to compute all the   
   >                eigenvalues, and no eigenvectors have been computed;   
   >                elements i+1:N of WR and WI contain eigenvalues which   
   >                have converged.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup doubleGEeigen   

    =====================================================================   
   Subroutine */ int igraphdgeev_(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl, 
	integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work, 
	integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, 
	    i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, k;
    doublereal r__, cs, sn;
    integer ihi;
    doublereal scl;
    integer ilo;
    doublereal dum[1], eps;
    integer ibal;
    char side[1];
    doublereal anrm;
    integer ierr, itau;
    extern /* Subroutine */ int igraphdrot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    integer iwrk, nout;
    extern doublereal igraphdnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int igraphdscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical igraphlsame_(char *, char *);
    extern doublereal igraphdlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int igraphdlabad_(doublereal *, doublereal *), igraphdgebak_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *), 
	    igraphdgebal_(char *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *);
    logical scalea;
    extern doublereal igraphdlamch_(char *);
    doublereal cscale;
    extern doublereal igraphdlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *);
    extern /* Subroutine */ int igraphdgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), igraphdlascl_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *);
    extern integer igraphidamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int igraphdlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *), 
	    igraphdlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), igraphxerbla_(char *, integer *, ftnlen);
    logical select[1];
    extern integer igraphilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    doublereal bignum;
    extern /* Subroutine */ int igraphdorghr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), igraphdhseqr_(char *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *), igraphdtrevc_(char *, char *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *);
    integer minwrk, maxwrk;
    logical wantvl;
    doublereal smlnum;
    integer hswork;
    logical lquery, wantvr;


/*  -- LAPACK driver routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --wr;
    --wi;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;

    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    wantvl = igraphlsame_(jobvl, "V");
    wantvr = igraphlsame_(jobvr, "V");
    if (! wantvl && ! igraphlsame_(jobvl, "N")) {
	*info = -1;
    } else if (! wantvr && ! igraphlsame_(jobvr, "N")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldvl < 1 || wantvl && *ldvl < *n) {
	*info = -9;
    } else if (*ldvr < 1 || wantvr && *ldvr < *n) {
	*info = -11;
    }

/*     Compute workspace   
        (Note: Comments in the code beginning "Workspace:" describe the   
         minimal amount of workspace needed at that point in the code,   
         as well as the preferred amount for good performance.   
         NB refers to the optimal block size for the immediately   
         following subroutine, as returned by ILAENV.   
         HSWORK refers to the workspace preferred by DHSEQR, as   
         calculated below. HSWORK is computed assuming ILO=1 and IHI=N,   
         the worst case.) */

    if (*info == 0) {
	if (*n == 0) {
	    minwrk = 1;
	    maxwrk = 1;
	} else {
	    maxwrk = (*n << 1) + *n * igraphilaenv_(&c__1, "DGEHRD", " ", n, &c__1, 
		    n, &c__0, (ftnlen)6, (ftnlen)1);
	    if (wantvl) {
		minwrk = *n << 2;
/* Computing MAX */
		i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * igraphilaenv_(&c__1, 
			"DORGHR", " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)
			1);
		maxwrk = max(i__1,i__2);
		igraphdhseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[
			1], &vl[vl_offset], ldvl, &work[1], &c_n1, info);
		hswork = (integer) work[1];
/* Computing MAX */
		i__1 = maxwrk, i__2 = *n + 1, i__1 = max(i__1,i__2), i__2 = *
			n + hswork;
		maxwrk = max(i__1,i__2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *n << 2;
		maxwrk = max(i__1,i__2);
	    } else if (wantvr) {
		minwrk = *n << 2;
/* Computing MAX */
		i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * igraphilaenv_(&c__1, 
			"DORGHR", " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)
			1);
		maxwrk = max(i__1,i__2);
		igraphdhseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[
			1], &vr[vr_offset], ldvr, &work[1], &c_n1, info);
		hswork = (integer) work[1];
/* Computing MAX */
		i__1 = maxwrk, i__2 = *n + 1, i__1 = max(i__1,i__2), i__2 = *
			n + hswork;
		maxwrk = max(i__1,i__2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *n << 2;
		maxwrk = max(i__1,i__2);
	    } else {
		minwrk = *n * 3;
		igraphdhseqr_("E", "N", n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[
			1], &vr[vr_offset], ldvr, &work[1], &c_n1, info);
		hswork = (integer) work[1];
/* Computing MAX */
		i__1 = maxwrk, i__2 = *n + 1, i__1 = max(i__1,i__2), i__2 = *
			n + hswork;
		maxwrk = max(i__1,i__2);
	    }
	    maxwrk = max(maxwrk,minwrk);
	}
	work[1] = (doublereal) maxwrk;

	if (*lwork < minwrk && ! lquery) {
	    *info = -13;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	igraphxerbla_("DGEEV ", &i__1, (ftnlen)6);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Get machine constants */

    eps = igraphdlamch_("P");
    smlnum = igraphdlamch_("S");
    bignum = 1. / smlnum;
    igraphdlabad_(&smlnum, &bignum);
    smlnum = sqrt(smlnum) / eps;
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

    anrm = igraphdlange_("M", n, n, &a[a_offset], lda, dum);
    scalea = FALSE_;
    if (anrm > 0. && anrm < smlnum) {
	scalea = TRUE_;
	cscale = smlnum;
    } else if (anrm > bignum) {
	scalea = TRUE_;
	cscale = bignum;
    }
    if (scalea) {
	igraphdlascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr);
    }

/*     Balance the matrix   
       (Workspace: need N) */

    ibal = 1;
    igraphdgebal_("B", n, &a[a_offset], lda, &ilo, &ihi, &work[ibal], &ierr);

/*     Reduce to upper Hessenberg form   
       (Workspace: need 3*N, prefer 2*N+N*NB) */

    itau = ibal + *n;
    iwrk = itau + *n;
    i__1 = *lwork - iwrk + 1;
    igraphdgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1,
	     &ierr);

    if (wantvl) {

/*        Want left eigenvectors   
          Copy Householder vectors to VL */

	*(unsigned char *)side = 'L';
	igraphdlacpy_("L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl)
		;

/*        Generate orthogonal matrix in VL   
          (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

	i__1 = *lwork - iwrk + 1;
	igraphdorghr_(n, &ilo, &ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk],
		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VL   
          (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

	iwrk = itau;
	i__1 = *lwork - iwrk + 1;
	igraphdhseqr_("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &
		vl[vl_offset], ldvl, &work[iwrk], &i__1, info);

	if (wantvr) {

/*           Want left and right eigenvectors   
             Copy Schur vectors to VR */

	    *(unsigned char *)side = 'B';
	    igraphdlacpy_("F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr);
	}

    } else if (wantvr) {

/*        Want right eigenvectors   
          Copy Householder vectors to VR */

	*(unsigned char *)side = 'R';
	igraphdlacpy_("L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr)
		;

/*        Generate orthogonal matrix in VR   
          (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

	i__1 = *lwork - iwrk + 1;
	igraphdorghr_(n, &ilo, &ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk],
		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VR   
          (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

	iwrk = itau;
	i__1 = *lwork - iwrk + 1;
	igraphdhseqr_("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &
		vr[vr_offset], ldvr, &work[iwrk], &i__1, info);

    } else {

/*        Compute eigenvalues only   
          (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

	iwrk = itau;
	i__1 = *lwork - iwrk + 1;
	igraphdhseqr_("E", "N", n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &
		vr[vr_offset], ldvr, &work[iwrk], &i__1, info);
    }

/*     If INFO > 0 from DHSEQR, then quit */

    if (*info > 0) {
	goto L50;
    }

    if (wantvl || wantvr) {

/*        Compute left and/or right eigenvectors   
          (Workspace: need 4*N) */

	igraphdtrevc_(side, "B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl,
		 &vr[vr_offset], ldvr, n, &nout, &work[iwrk], &ierr);
    }

    if (wantvl) {

/*        Undo balancing of left eigenvectors   
          (Workspace: need N) */

	igraphdgebak_("B", "L", n, &ilo, &ihi, &work[ibal], n, &vl[vl_offset], ldvl,
		 &ierr);

/*        Normalize left eigenvectors and make largest component real */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (wi[i__] == 0.) {
		scl = 1. / igraphdnrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
		igraphdscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
	    } else if (wi[i__] > 0.) {
		d__1 = igraphdnrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
		d__2 = igraphdnrm2_(n, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
		scl = 1. / igraphdlapy2_(&d__1, &d__2);
		igraphdscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
		igraphdscal_(n, &scl, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
		i__2 = *n;
		for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
		    d__1 = vl[k + i__ * vl_dim1];
/* Computing 2nd power */
		    d__2 = vl[k + (i__ + 1) * vl_dim1];
		    work[iwrk + k - 1] = d__1 * d__1 + d__2 * d__2;
/* L10: */
		}
		k = igraphidamax_(n, &work[iwrk], &c__1);
		igraphdlartg_(&vl[k + i__ * vl_dim1], &vl[k + (i__ + 1) * vl_dim1], 
			&cs, &sn, &r__);
		igraphdrot_(n, &vl[i__ * vl_dim1 + 1], &c__1, &vl[(i__ + 1) * 
			vl_dim1 + 1], &c__1, &cs, &sn);
		vl[k + (i__ + 1) * vl_dim1] = 0.;
	    }
/* L20: */
	}
    }

    if (wantvr) {

/*        Undo balancing of right eigenvectors   
          (Workspace: need N) */

	igraphdgebak_("B", "R", n, &ilo, &ihi, &work[ibal], n, &vr[vr_offset], ldvr,
		 &ierr);

/*        Normalize right eigenvectors and make largest component real */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (wi[i__] == 0.) {
		scl = 1. / igraphdnrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
		igraphdscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
	    } else if (wi[i__] > 0.) {
		d__1 = igraphdnrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
		d__2 = igraphdnrm2_(n, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
		scl = 1. / igraphdlapy2_(&d__1, &d__2);
		igraphdscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
		igraphdscal_(n, &scl, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
		i__2 = *n;
		for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
		    d__1 = vr[k + i__ * vr_dim1];
/* Computing 2nd power */
		    d__2 = vr[k + (i__ + 1) * vr_dim1];
		    work[iwrk + k - 1] = d__1 * d__1 + d__2 * d__2;
/* L30: */
		}
		k = igraphidamax_(n, &work[iwrk], &c__1);
		igraphdlartg_(&vr[k + i__ * vr_dim1], &vr[k + (i__ + 1) * vr_dim1], 
			&cs, &sn, &r__);
		igraphdrot_(n, &vr[i__ * vr_dim1 + 1], &c__1, &vr[(i__ + 1) * 
			vr_dim1 + 1], &c__1, &cs, &sn);
		vr[k + (i__ + 1) * vr_dim1] = 0.;
	    }
/* L40: */
	}
    }

/*     Undo scaling if necessary */

L50:
    if (scalea) {
	i__1 = *n - *info;
/* Computing MAX */
	i__3 = *n - *info;
	i__2 = max(i__3,1);
	igraphdlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[*info + 
		1], &i__2, &ierr);
	i__1 = *n - *info;
/* Computing MAX */
	i__3 = *n - *info;
	i__2 = max(i__3,1);
	igraphdlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[*info + 
		1], &i__2, &ierr);
	if (*info > 0) {
	    i__1 = ilo - 1;
	    igraphdlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[1], 
		    n, &ierr);
	    i__1 = ilo - 1;
	    igraphdlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[1], 
		    n, &ierr);
	}
    }

    work[1] = (doublereal) maxwrk;
    return 0;

/*     End of DGEEV */

} /* igraphdgeev_ */

