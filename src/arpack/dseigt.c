/* dseigt.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "config.h"
#include "arpack_internal.h"

/* Common Block Declarations */

struct {
    integer logfil, ndigit, mgetv0, msaupd, msaup2, msaitr, mseigt, msapps, 
	    msgets, mseupd, mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, 
	    mneupd, mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd;
} debug_;

#define debug_1 debug_

struct {
    integer nopx, nbx, nrorth, nitref, nrstrt;
    real tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv, tnaupd, 
	    tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv, tcaupd, tcaup2, 
	    tcaitr, tceigh, tcgets, tcapps, tcconv, tmvopx, tmvbx, tgetv0, 
	    titref, trvec;
} timing_;

#define timing_1 timing_

/* Table of constant values */

static integer c__1 = 1;

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: dseigt */

/* \Description: */
/*  Compute the eigenvalues of the current symmetric tridiagonal matrix */
/*  and the corresponding error bounds given the current residual norm. */

/* \Usage: */
/*  call dseigt */
/*     ( RNORM, N, H, LDH, EIG, BOUNDS, WORKL, IERR ) */

/* \Arguments */
/*  RNORM   Double precision scalar.  (INPUT) */
/*          RNORM contains the residual norm corresponding to the current */
/*          symmetric tridiagonal matrix H. */

/*  N       Integer.  (INPUT) */
/*          Size of the symmetric tridiagonal matrix H. */

/*  H       Double precision N by 2 array.  (INPUT) */
/*          H contains the symmetric tridiagonal matrix with the */
/*          subdiagonal in the first column starting at H(2,1) and the */
/*          main diagonal in second column. */

/*  LDH     Integer.  (INPUT) */
/*          Leading dimension of H exactly as declared in the calling */
/*          program. */

/*  EIG     Double precision array of length N.  (OUTPUT) */
/*          On output, EIG contains the N eigenvalues of H possibly */
/*          unsorted.  The BOUNDS arrays are returned in the */
/*          same sorted order as EIG. */

/*  BOUNDS  Double precision array of length N.  (OUTPUT) */
/*          On output, BOUNDS contains the error estimates corresponding */
/*          to the eigenvalues EIG.  This is equal to RNORM times the */
/*          last components of the eigenvectors corresponding to the */
/*          eigenvalues in EIG. */

/*  WORKL   Double precision work array of length 3*N.  (WORKSPACE) */
/*          Private (replicated) array on each PE or array allocated on */
/*          the front end. */

/*  IERR    Integer.  (OUTPUT) */
/*          Error exit flag from dstqrb. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  real */

/* \Routines called: */
/*     dstqrb  ARPACK routine that computes the eigenvalues and the */
/*             last components of the eigenvectors of a symmetric */
/*             and tridiagonal matrix. */
/*     second  ARPACK utility routine for timing. */
/*     dvout   ARPACK utility routine that prints vectors. */
/*     dcopy   Level 1 BLAS that copies one vector to another. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \Revision history: */
/*     xx/xx/92: Version ' 2.4' */

/* \SCCS Information: @(#) */
/* FILE: seigt.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2 */

/* \Remarks */
/*     None */

/* \EndLib */

/* ----------------------------------------------------------------------- */

/* Subroutine */ int igraphdseigt_(rnorm, n, h__, ldh, eig, bounds, workl, ierr)
doublereal *rnorm;
integer *n;
doublereal *h__;
integer *ldh;
doublereal *eig, *bounds, *workl;
integer *ierr;
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int igraphdcopy_(), igraphdvout_();
    static real t0, t1;
    extern /* Subroutine */ int igraphsecond_();
    static integer msglvl;
    extern /* Subroutine */ int igraphdstqrb_();


/*     %----------------------------------------------------% */
/*     | Include files for debugging and timing information | */
/*     %----------------------------------------------------% */


/* \SCCS Information: @(#) */
/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

/*     %---------------------------------% */
/*     | See debug.doc for documentation | */
/*     %---------------------------------% */

/*     %------------------% */
/*     | Scalar Arguments | */
/*     %------------------% */

/*     %--------------------------------% */
/*     | See stat.doc for documentation | */
/*     %--------------------------------% */

/* \SCCS Information: @(#) */
/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */



/*     %-----------------% */
/*     | Array Arguments | */
/*     %-----------------% */


/*     %------------% */
/*     | Parameters | */
/*     %------------% */


/*     %---------------% */
/*     | Local Scalars | */
/*     %---------------% */


/*     %----------------------% */
/*     | External Subroutines | */
/*     %----------------------% */


/*     %-----------------------% */
/*     | Executable Statements | */
/*     %-----------------------% */

/*     %-------------------------------% */
/*     | Initialize timing statistics  | */
/*     | & message level for debugging | */
/*     %-------------------------------% */

    /* Parameter adjustments */
    --workl;
    --bounds;
    --eig;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1 * 1;
    h__ -= h_offset;

    /* Function Body */
    igraphsecond_(&t0);
    msglvl = debug_1.mseigt;

    if (msglvl > 0) {
	igraphdvout_(&debug_1.logfil, n, &h__[(h_dim1 << 1) + 1], &debug_1.ndigit, 
		"_seigt: main diagonal of matrix H", (ftnlen)33);
	if (*n > 1) {
	    i__1 = *n - 1;
	    igraphdvout_(&debug_1.logfil, &i__1, &h__[h_dim1 + 2], &debug_1.ndigit, 
		    "_seigt: sub diagonal of matrix H", (ftnlen)32);
	}
    }

    igraphdcopy_(n, &h__[(h_dim1 << 1) + 1], &c__1, &eig[1], &c__1);
    i__1 = *n - 1;
    igraphdcopy_(&i__1, &h__[h_dim1 + 2], &c__1, &workl[1], &c__1);
    igraphdstqrb_(n, &eig[1], &workl[1], &bounds[1], &workl[*n + 1], ierr);
    if (*ierr != 0) {
	goto L9000;
    }
    if (msglvl > 1) {
	igraphdvout_(&debug_1.logfil, n, &bounds[1], &debug_1.ndigit, "_seigt: las\
t row of the eigenvector matrix for H", (ftnlen)48);
    }

/*     %-----------------------------------------------% */
/*     | Finally determine the error bounds associated | */
/*     | with the n Ritz values of H.                  | */
/*     %-----------------------------------------------% */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	bounds[k] = *rnorm * (d__1 = bounds[k], abs(d__1));
/* L30: */
    }

    igraphsecond_(&t1);
    timing_1.tseigt += t1 - t0;

L9000:
    return 0;

/*     %---------------% */
/*     | End of dseigt | */
/*     %---------------% */

} /* igraphdseigt_ */

