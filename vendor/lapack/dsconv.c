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

static doublereal c_b3 = .66666666666666663;

/* -----------------------------------------------------------------------   
   \BeginDoc   

   \Name: dsconv   

   \Description:   
    Convergence testing for the symmetric Arnoldi eigenvalue routine.   

   \Usage:   
    call dsconv   
       ( N, RITZ, BOUNDS, TOL, NCONV )   

   \Arguments   
    N       Integer.  (INPUT)   
            Number of Ritz values to check for convergence.   

    RITZ    Double precision array of length N.  (INPUT)   
            The Ritz values to be checked for convergence.   

    BOUNDS  Double precision array of length N.  (INPUT)   
            Ritz estimates associated with the Ritz values in RITZ.   

    TOL     Double precision scalar.  (INPUT)   
            Desired relative accuracy for a Ritz value to be considered   
            "converged".   

    NCONV   Integer scalar.  (OUTPUT)   
            Number of "converged" Ritz values.   

   \EndDoc   

   -----------------------------------------------------------------------   

   \BeginLib   

   \Routines called:   
       second  ARPACK utility routine for timing.   
       dlamch  LAPACK routine that determines machine constants.   

   \Author   
       Danny Sorensen               Phuong Vu   
       Richard Lehoucq              CRPC / Rice University   
       Dept. of Computational &     Houston, Texas   
       Applied Mathematics   
       Rice University   
       Houston, Texas   

   \SCCS Information: @(#)   
   FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2   

   \Remarks   
       1. Starting with version 2.4, this routine no longer uses the   
          Parlett strategy using the gap conditions.   

   \EndLib   

   -----------------------------------------------------------------------   

   Subroutine */ int igraphdsconv_(integer *n, doublereal *ritz, doublereal *bounds,
	 doublereal *tol, integer *nconv)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    integer i__;
    IGRAPH_F77_SAVE real t0, t1;
    doublereal eps23, temp;
    extern doublereal igraphdlamch_(char *);
    extern /* Subroutine */ int igraphsecond_(real *);
    real tsconv = 0;


/*     %----------------------------------------------------%   
       | Include files for debugging and timing information |   
       %----------------------------------------------------%   


       %------------------%   
       | Scalar Arguments |   
       %------------------%   


       %-----------------%   
       | Array Arguments |   
       %-----------------%   


       %---------------%   
       | Local Scalars |   
       %---------------%   


       %-------------------%   
       | External routines |   
       %-------------------%   

       %---------------------%   
       | Intrinsic Functions |   
       %---------------------%   


       %-----------------------%   
       | Executable Statements |   
       %-----------------------%   

       Parameter adjustments */
    --bounds;
    --ritz;

    /* Function Body */
    igraphsecond_(&t0);

    eps23 = igraphdlamch_("Epsilon-Machine");
    eps23 = pow_dd(&eps23, &c_b3);

    *nconv = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        %-----------------------------------------------------%   
          | The i-th Ritz value is considered "converged"       |   
          | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   |   
          %-----------------------------------------------------%   

   Computing MAX */
	d__2 = eps23, d__3 = (d__1 = ritz[i__], abs(d__1));
	temp = max(d__2,d__3);
	if (bounds[i__] <= *tol * temp) {
	    ++(*nconv);
	}

/* L10: */
    }

    igraphsecond_(&t1);
    tsconv += t1 - t0;

    return 0;

/*     %---------------%   
       | End of dsconv |   
       %---------------% */

} /* igraphdsconv_ */

