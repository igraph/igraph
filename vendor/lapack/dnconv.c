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

   \Name: dnconv   

   \Description:   
    Convergence testing for the nonsymmetric Arnoldi eigenvalue routine.   

   \Usage:   
    call dnconv   
       ( N, RITZR, RITZI, BOUNDS, TOL, NCONV )   

   \Arguments   
    N       Integer.  (INPUT)   
            Number of Ritz values to check for convergence.   

    RITZR,  Double precision arrays of length N.  (INPUT)   
    RITZI   Real and imaginary parts of the Ritz values to be checked   
            for convergence.   
    BOUNDS  Double precision array of length N.  (INPUT)   
            Ritz estimates for the Ritz values in RITZR and RITZI.   

    TOL     Double precision scalar.  (INPUT)   
            Desired backward error for a Ritz value to be considered   
            "converged".   

    NCONV   Integer scalar.  (OUTPUT)   
            Number of "converged" Ritz values.   

   \EndDoc   

   -----------------------------------------------------------------------   

   \BeginLib   

   \Local variables:   
       xxxxxx  real   

   \Routines called:   
       second  ARPACK utility routine for timing.   
       dlamch  LAPACK routine that determines machine constants.   
       dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.   

   \Author   
       Danny Sorensen               Phuong Vu   
       Richard Lehoucq              CRPC / Rice University   
       Dept. of Computational &     Houston, Texas   
       Applied Mathematics   
       Rice University   
       Houston, Texas   

   \Revision history:   
       xx/xx/92: Version ' 2.1'   

   \SCCS Information: @(#)   
   FILE: nconv.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2   

   \Remarks   
       1. xxxx   

   \EndLib   

   -----------------------------------------------------------------------   

   Subroutine */ int igraphdnconv_(integer *n, doublereal *ritzr, doublereal *ritzi,
	 doublereal *bounds, doublereal *tol, integer *nconv)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    integer i__;
    IGRAPH_F77_SAVE real t0, t1;
    doublereal eps23, temp;
    extern doublereal igraphdlapy2_(doublereal *, doublereal *), igraphdlamch_(char *);
    extern /* Subroutine */ int igraphsecond_(real *);
    real tnconv = 0.;


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


       %--------------------%   
       | External Functions |   
       %--------------------%   

       %-----------------------%   
       | Executable Statements |   
       %-----------------------%   

       %-------------------------------------------------------------%   
       | Convergence test: unlike in the symmetric code, I am not    |   
       | using things like refined error bounds and gap condition    |   
       | because I don't know the exact equivalent concept.          |   
       |                                                             |   
       | Instead the i-th Ritz value is considered "converged" when: |   
       |                                                             |   
       |     bounds(i) .le. ( TOL * | ritz | )                       |   
       |                                                             |   
       | for some appropriate choice of norm.                        |   
       %-------------------------------------------------------------%   

       Parameter adjustments */
    --bounds;
    --ritzi;
    --ritzr;

    /* Function Body */
    igraphsecond_(&t0);

/*     %---------------------------------%   
       | Get machine dependent constant. |   
       %---------------------------------% */

    eps23 = igraphdlamch_("Epsilon-Machine");
    eps23 = pow_dd(&eps23, &c_b3);

    *nconv = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__1 = eps23, d__2 = igraphdlapy2_(&ritzr[i__], &ritzi[i__]);
	temp = max(d__1,d__2);
	if (bounds[i__] <= *tol * temp) {
	    ++(*nconv);
	}
/* L20: */
    }

    igraphsecond_(&t1);
    tnconv += t1 - t0;

    return 0;

/*     %---------------%   
       | End of dnconv |   
       %---------------% */

} /* igraphdnconv_ */

