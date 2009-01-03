/* igraphdtrsen.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "config.h"
#include "arpack_internal.h"

/* Table of constant values */

static integer c_n1 = -1;

/* Subroutine */ int igraphdtrsen_(job, compq, select, n, t, ldt, q, ldq, wr, wi, m,
	 s, sep, work, lwork, iwork, liwork, info, job_len, compq_len)
char *job, *compq;
logical *select;
integer *n;
doublereal *t;
integer *ldt;
doublereal *q;
integer *ldq;
doublereal *wr, *wi;
integer *m;
doublereal *s, *sep, *work;
integer *lwork, *iwork, *liwork, *info;
ftnlen job_len;
ftnlen compq_len;
{
    /* System generated locals */
    integer q_dim1, q_offset, t_dim1, t_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer kase;
    static logical pair;
    static integer ierr;
    static logical swap;
    static integer k;
    static doublereal scale;
    extern logical igraphlsame_();
    static logical wantq, wants;
    static doublereal rnorm;
    static integer n1, n2, kk;
    extern doublereal igraphdlange_();
    static integer nn, ks;
    extern /* Subroutine */ int igraphdlacon_(), igraphdlacpy_(), igraphxerbla_();
    static logical wantbh;
    extern /* Subroutine */ int igraphdtrexc_();
    static logical wantsp;
    extern /* Subroutine */ int igraphdtrsyl_();
    static doublereal est;


/*  -- LAPACK routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTRSEN reorders the real Schur factorization of a real matrix */
/*  A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in */
/*  the leading diagonal blocks of the upper quasi-triangular matrix T, */
/*  and the leading columns of Q form an orthonormal basis of the */
/*  corresponding right invariant subspace. */

/*  Optionally the routine computes the reciprocal condition numbers of */
/*  the cluster of eigenvalues and/or the invariant subspace. */

/*  T must be in Schur canonical form (as returned by DHSEQR), that is, */
/*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each */
/*  2-by-2 diagonal block has its diagonal elemnts equal and its */
/*  off-diagonal elements of opposite sign. */

/*  Arguments */
/*  ========= */

/*  JOB     (input) CHARACTER*1 */
/*          Specifies whether condition numbers are required for the */
/*          cluster of eigenvalues (S) or the invariant subspace (SEP): */
/*          = 'N': none; */
/*          = 'E': for eigenvalues only (S); */
/*          = 'V': for invariant subspace only (SEP); */
/*          = 'B': for both eigenvalues and invariant subspace (S and */
/*                 SEP). */

/*  COMPQ   (input) CHARACTER*1 */
/*          = 'V': update the matrix Q of Schur vectors; */
/*          = 'N': do not update Q. */

/*  SELECT  (input) LOGICAL array, dimension (N) */
/*          SELECT specifies the eigenvalues in the selected cluster. To */
/*          select a real eigenvalue w(j), SELECT(j) must be set to */
/*          .TRUE.. To select a complex conjugate pair of eigenvalues */
/*          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block, */
/*          either SELECT(j) or SELECT(j+1) or both must be set to */
/*          .TRUE.; a complex conjugate pair of eigenvalues must be */
/*          either both included in the cluster or both excluded. */

/*  N       (input) INTEGER */
/*          The order of the matrix T. N >= 0. */

/*  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N) */
/*          On entry, the upper quasi-triangular matrix T, in Schur */
/*          canonical form. */
/*          On exit, T is overwritten by the reordered matrix T, again in */
/*          Schur canonical form, with the selected eigenvalues in the */
/*          leading diagonal blocks. */

/*  LDT     (input) INTEGER */
/*          The leading dimension of the array T. LDT >= max(1,N). */

/*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*          On entry, if COMPQ = 'V', the matrix Q of Schur vectors. */
/*          On exit, if COMPQ = 'V', Q has been postmultiplied by the */
/*          orthogonal transformation matrix which reorders T; the */
/*          leading M columns of Q form an orthonormal basis for the */
/*          specified invariant subspace. */
/*          If COMPQ = 'N', Q is not referenced. */

/*  LDQ     (input) INTEGER */
/*          The leading dimension of the array Q. */
/*          LDQ >= 1; and if COMPQ = 'V', LDQ >= N. */

/*  WR      (output) DOUBLE PRECISION array, dimension (N) */
/*  WI      (output) DOUBLE PRECISION array, dimension (N) */
/*          The real and imaginary parts, respectively, of the reordered */
/*          eigenvalues of T. The eigenvalues are stored in the same */
/*          order as on the diagonal of T, with WR(i) = T(i,i) and, if */
/*          T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 and */
/*          WI(i+1) = -WI(i). Note that if a complex eigenvalue is */
/*          sufficiently ill-conditioned, then its value may differ */
/*          significantly from its value before reordering. */

/*  M       (output) INTEGER */
/*          The dimension of the specified invariant subspace. */
/*          0 < = M <= N. */

/*  S       (output) DOUBLE PRECISION */
/*          If JOB = 'E' or 'B', S is a lower bound on the reciprocal */
/*          condition number for the selected cluster of eigenvalues. */
/*          S cannot underestimate the true reciprocal condition number */
/*          by more than a factor of sqrt(N). If M = 0 or N, S = 1. */
/*          If JOB = 'N' or 'V', S is not referenced. */

/*  SEP     (output) DOUBLE PRECISION */
/*          If JOB = 'V' or 'B', SEP is the estimated reciprocal */
/*          condition number of the specified invariant subspace. If */
/*          M = 0 or N, SEP = norm(T). */
/*          If JOB = 'N' or 'E', SEP is not referenced. */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK) */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK. */
/*          If JOB = 'N', LWORK >= max(1,N); */
/*          if JOB = 'E', LWORK >= M*(N-M); */
/*          if JOB = 'V' or 'B', LWORK >= 2*M*(N-M). */

/*  IWORK   (workspace) INTEGER array, dimension (LIWORK) */
/*          IF JOB = 'N' or 'E', IWORK is not referenced. */

/*  LIWORK  (input) INTEGER */
/*          The dimension of the array IWORK. */
/*          If JOB = 'N' or 'E', LIWORK >= 1; */
/*          if JOB = 'V' or 'B', LIWORK >= M*(N-M). */

/*  INFO    (output) INTEGER */
/*          = 0: successful exit */
/*          < 0: if INFO = -i, the i-th argument had an illegal value */
/*          = 1: reordering of T failed because some eigenvalues are too */
/*               close to separate (the problem is very ill-conditioned); */
/*               T may have been partially reordered, and WR and WI */
/*               contain the eigenvalues in the same order as in T; S and */
/*               SEP (if requested) are set to zero. */

/*  Further Details */
/*  =============== */

/*  DTRSEN first collects the selected eigenvalues by computing an */
/*  orthogonal transformation Z to move them to the top left corner of T. */
/*  In other words, the selected eigenvalues are the eigenvalues of T11 */
/*  in: */

/*                Z'*T*Z = ( T11 T12 ) n1 */
/*                         (  0  T22 ) n2 */
/*                            n1  n2 */

/*  where N = n1+n2 and Z' means the transpose of Z. The first n1 columns */
/*  of Z span the specified invariant subspace of T. */

/*  If T has been obtained from the real Schur factorization of a matrix */
/*  A = Q*T*Q', then the reordered real Schur factorization of A is given */
/*  by A = (Q*Z)*(Z'*T*Z)*(Q*Z)', and the first n1 columns of Q*Z span */
/*  the corresponding invariant subspace of A. */

/*  The reciprocal condition number of the average of the eigenvalues of */
/*  T11 may be returned in S. S lies between 0 (very badly conditioned) */
/*  and 1 (very well conditioned). It is computed as follows. First we */
/*  compute R so that */

/*                         P = ( I  R ) n1 */
/*                             ( 0  0 ) n2 */
/*                               n1 n2 */

/*  is the projector on the invariant subspace associated with T11. */
/*  R is the solution of the Sylvester equation: */

/*                        T11*R - R*T22 = T12. */

/*  Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote */
/*  the two-norm of M. Then S is computed as the lower bound */

/*                      (1 + F-norm(R)**2)**(-1/2) */

/*  on the reciprocal of 2-norm(P), the true reciprocal condition number. */
/*  S cannot underestimate 1 / 2-norm(P) by more than a factor of */
/*  sqrt(N). */

/*  An approximate error bound for the computed average of the */
/*  eigenvalues of T11 is */

/*                         EPS * norm(T) / S */

/*  where EPS is the machine precision. */

/*  The reciprocal condition number of the right invariant subspace */
/*  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP. */
/*  SEP is defined as the separation of T11 and T22: */

/*                     sep( T11, T22 ) = sigma-min( C ) */

/*  where sigma-min(C) is the smallest singular value of the */
/*  n1*n2-by-n1*n2 matrix */

/*     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) ) */

/*  I(m) is an m by m identity matrix, and kprod denotes the Kronecker */
/*  product. We estimate sigma-min(C) by the reciprocal of an estimate of */
/*  the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C) */
/*  cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2). */

/*  When SEP is small, small changes in T can cause large changes in */
/*  the invariant subspace. An approximate bound on the maximum angular */
/*  error in the computed right invariant subspace is */

/*                      EPS * norm(T) / SEP */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and test the input parameters */

    /* Parameter adjustments */
    --select;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1 * 1;
    t -= t_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1 * 1;
    q -= q_offset;
    --wr;
    --wi;
    --work;
    --iwork;

    /* Function Body */
    wantbh = igraphlsame_(job, "B", (ftnlen)1, (ftnlen)1);
    wants = igraphlsame_(job, "E", (ftnlen)1, (ftnlen)1) || wantbh;
    wantsp = igraphlsame_(job, "V", (ftnlen)1, (ftnlen)1) || wantbh;
    wantq = igraphlsame_(compq, "V", (ftnlen)1, (ftnlen)1);

    *info = 0;
    if (! igraphlsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! wants && ! wantsp) {
	*info = -1;
    } else if (! igraphlsame_(compq, "N", (ftnlen)1, (ftnlen)1) && ! wantq) {
	*info = -2;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ldt < max(1,*n)) {
	*info = -6;
    } else if (*ldq < 1 || wantq && *ldq < *n) {
	*info = -8;
    } else {

/*        Set M to the dimension of the specified invariant subspace, */
/*        and test LWORK and LIWORK. */

	*m = 0;
	pair = FALSE_;
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    if (pair) {
		pair = FALSE_;
	    } else {
		if (k < *n) {
		    if (t[k + 1 + k * t_dim1] == 0.) {
			if (select[k]) {
			    ++(*m);
			}
		    } else {
			pair = TRUE_;
			if (select[k] || select[k + 1]) {
			    *m += 2;
			}
		    }
		} else {
		    if (select[*n]) {
			++(*m);
		    }
		}
	    }
/* L10: */
	}

	n1 = *m;
	n2 = *n - *m;
	nn = n1 * n2;

	if (*lwork < 1 || wants && ! wantsp && *lwork < nn || wantsp && *
		lwork < nn << 1) {
	    *info = -15;
	} else if (*liwork < 1 || wantsp && *liwork < nn) {
	    *info = -17;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	igraphxerbla_("DTRSEN", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == *n || *m == 0) {
	if (wants) {
	    *s = 1.;
	}
	if (wantsp) {
	    *sep = igraphdlange_("1", n, n, &t[t_offset], ldt, &work[1], (ftnlen)1);
	}
	goto L40;
    }

/*     Collect the selected blocks at the top-left corner of T. */

    ks = 0;
    pair = FALSE_;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (pair) {
	    pair = FALSE_;
	} else {
	    swap = select[k];
	    if (k < *n) {
		if (t[k + 1 + k * t_dim1] != 0.) {
		    pair = TRUE_;
		    swap = swap || select[k + 1];
		}
	    }
	    if (swap) {
		++ks;

/*              Swap the K-th block to position KS. */

		ierr = 0;
		kk = k;
		if (k != ks) {
		    igraphdtrexc_(compq, n, &t[t_offset], ldt, &q[q_offset], ldq, &
			    kk, &ks, &work[1], &ierr, (ftnlen)1);
		}
		if (ierr == 1 || ierr == 2) {

/*                 Blocks too close to swap: exit. */

		    *info = 1;
		    if (wants) {
			*s = 0.;
		    }
		    if (wantsp) {
			*sep = 0.;
		    }
		    goto L40;
		}
		if (pair) {
		    ++ks;
		}
	    }
	}
/* L20: */
    }

    if (wants) {

/*        Solve Sylvester equation for R: */

/*           T11*R - R*T22 = scale*T12 */

	igraphdlacpy_("F", &n1, &n2, &t[(n1 + 1) * t_dim1 + 1], ldt, &work[1], &n1, 
		(ftnlen)1);
	igraphdtrsyl_("N", "N", &c_n1, &n1, &n2, &t[t_offset], ldt, &t[n1 + 1 + (n1 
		+ 1) * t_dim1], ldt, &work[1], &n1, &scale, &ierr, (ftnlen)1, 
		(ftnlen)1);

/*        Estimate the reciprocal of the condition number of the cluster */
/*        of eigenvalues. */

	rnorm = igraphdlange_("F", &n1, &n2, &work[1], &n1, &work[1], (ftnlen)1);
	if (rnorm == 0.) {
	    *s = 1.;
	} else {
	    *s = scale / (sqrt(scale * scale / rnorm + rnorm) * sqrt(rnorm));
	}
    }

    if (wantsp) {

/*        Estimate sep(T11,T22). */

	est = 0.;
	kase = 0;
L30:
	igraphdlacon_(&nn, &work[nn + 1], &work[1], &iwork[1], &est, &kase);
	if (kase != 0) {
	    if (kase == 1) {

/*              Solve  T11*R - R*T22 = scale*X. */

		igraphdtrsyl_("N", "N", &c_n1, &n1, &n2, &t[t_offset], ldt, &t[n1 + 
			1 + (n1 + 1) * t_dim1], ldt, &work[1], &n1, &scale, &
			ierr, (ftnlen)1, (ftnlen)1);
	    } else {

/*              Solve  T11'*R - R*T22' = scale*X. */

		igraphdtrsyl_("T", "T", &c_n1, &n1, &n2, &t[t_offset], ldt, &t[n1 + 
			1 + (n1 + 1) * t_dim1], ldt, &work[1], &n1, &scale, &
			ierr, (ftnlen)1, (ftnlen)1);
	    }
	    goto L30;
	}

	*sep = scale / est;
    }

L40:

/*     Store the output eigenvalues in WR and WI. */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	wr[k] = t[k + k * t_dim1];
	wi[k] = 0.;
/* L50: */
    }
    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	if (t[k + 1 + k * t_dim1] != 0.) {
	    wi[k] = sqrt((d__1 = t[k + (k + 1) * t_dim1], abs(d__1))) * sqrt((
		    d__2 = t[k + 1 + k * t_dim1], abs(d__2)));
	    wi[k + 1] = -wi[k];
	}
/* L60: */
    }
    return 0;

/*     End of DTRSEN */

} /* igraphdtrsen_ */

