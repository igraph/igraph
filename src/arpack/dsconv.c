/* dsconv.f -- translated by f2c (version 19991025).
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

static doublereal c_b3 = .66666666666666663;

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: dsconv */

/* \Description: */
/*  Convergence testing for the symmetric Arnoldi eigenvalue routine. */

/* \Usage: */
/*  call dsconv */
/*     ( N, RITZ, BOUNDS, TOL, NCONV ) */

/* \Arguments */
/*  N       Integer.  (INPUT) */
/*          Number of Ritz values to check for convergence. */

/*  RITZ    Double precision array of length N.  (INPUT) */
/*          The Ritz values to be checked for convergence. */

/*  BOUNDS  Double precision array of length N.  (INPUT) */
/*          Ritz estimates associated with the Ritz values in RITZ. */

/*  TOL     Double precision scalar.  (INPUT) */
/*          Desired relative accuracy for a Ritz value to be considered */
/*          "converged". */

/*  NCONV   Integer scalar.  (OUTPUT) */
/*          Number of "converged" Ritz values. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Routines called: */
/*     second  ARPACK utility routine for timing. */
/*     dlamch  LAPACK routine that determines machine constants. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2 */

/* \Remarks */
/*     1. Starting with version 2.4, this routine no longer uses the */
/*        Parlett strategy using the gap conditions. */

/* \EndLib */

/* ----------------------------------------------------------------------- */

/* Subroutine */ int igraphdsconv_(n, ritz, bounds, tol, nconv)
integer *n;
doublereal *ritz, *bounds, *tol;
integer *nconv;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double igraphpow_dd();

    /* Local variables */
    static doublereal eps23, temp;
    static integer i__;
    static real t0, t1;
    extern doublereal igraphdlamch_();
    extern /* Subroutine */ int igraphsecond_();


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


/*     %---------------% */
/*     | Local Scalars | */
/*     %---------------% */


/*     %-------------------% */
/*     | External routines | */
/*     %-------------------% */

/*     %---------------------% */
/*     | Intrinsic Functions | */
/*     %---------------------% */


/*     %-----------------------% */
/*     | Executable Statements | */
/*     %-----------------------% */

    /* Parameter adjustments */
    --bounds;
    --ritz;

    /* Function Body */
    igraphsecond_(&t0);

    eps23 = igraphdlamch_("Epsilon-Machine", (ftnlen)15);
    eps23 = igraphpow_dd(&eps23, &c_b3);

    *nconv = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        %-----------------------------------------------------% */
/*        | The i-th Ritz value is considered "converged"       | */
/*        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   | */
/*        %-----------------------------------------------------% */

/* Computing MAX */
	d__2 = eps23, d__3 = (d__1 = ritz[i__], abs(d__1));
	temp = max(d__2,d__3);
	if (bounds[i__] <= *tol * temp) {
	    ++(*nconv);
	}

/* L10: */
    }

    igraphsecond_(&t1);
    timing_1.tsconv += t1 - t0;

    return 0;

/*     %---------------% */
/*     | End of dsconv | */
/*     %---------------% */

} /* igraphdsconv_ */

