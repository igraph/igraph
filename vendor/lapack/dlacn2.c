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
static doublereal c_b11 = 1.;

/* > \brief \b DLACN2 estimates the 1-norm of a square matrix, using reverse communication for evaluating matr
ix-vector products.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLACN2 + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacn2.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacn2.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacn2.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DLACN2( N, V, X, ISGN, EST, KASE, ISAVE )   

         INTEGER            KASE, N   
         DOUBLE PRECISION   EST   
         INTEGER            ISGN( * ), ISAVE( 3 )   
         DOUBLE PRECISION   V( * ), X( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DLACN2 estimates the 1-norm of a square, real matrix A.   
   > Reverse communication is used for evaluating matrix-vector products.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >         The order of the matrix.  N >= 1.   
   > \endverbatim   
   >   
   > \param[out] V   
   > \verbatim   
   >          V is DOUBLE PRECISION array, dimension (N)   
   >         On the final return, V = A*W,  where  EST = norm(V)/norm(W)   
   >         (W is not returned).   
   > \endverbatim   
   >   
   > \param[in,out] X   
   > \verbatim   
   >          X is DOUBLE PRECISION array, dimension (N)   
   >         On an intermediate return, X should be overwritten by   
   >               A * X,   if KASE=1,   
   >               A**T * X,  if KASE=2,   
   >         and DLACN2 must be re-called with all the other parameters   
   >         unchanged.   
   > \endverbatim   
   >   
   > \param[out] ISGN   
   > \verbatim   
   >          ISGN is INTEGER array, dimension (N)   
   > \endverbatim   
   >   
   > \param[in,out] EST   
   > \verbatim   
   >          EST is DOUBLE PRECISION   
   >         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be   
   >         unchanged from the previous call to DLACN2.   
   >         On exit, EST is an estimate (a lower bound) for norm(A).   
   > \endverbatim   
   >   
   > \param[in,out] KASE   
   > \verbatim   
   >          KASE is INTEGER   
   >         On the initial call to DLACN2, KASE should be 0.   
   >         On an intermediate return, KASE will be 1 or 2, indicating   
   >         whether X should be overwritten by A * X  or A**T * X.   
   >         On the final return from DLACN2, KASE will again be 0.   
   > \endverbatim   
   >   
   > \param[in,out] ISAVE   
   > \verbatim   
   >          ISAVE is INTEGER array, dimension (3)   
   >         ISAVE is used to save variables between calls to DLACN2   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup doubleOTHERauxiliary   

   > \par Further Details:   
    =====================   
   >   
   > \verbatim   
   >   
   >  Originally named SONEST, dated March 16, 1988.   
   >   
   >  This is a thread safe version of DLACON, which uses the array ISAVE   
   >  in place of a SAVE statement, as follows:   
   >   
   >     DLACON     DLACN2   
   >      JUMP     ISAVE(1)   
   >      J        ISAVE(2)   
   >      ITER     ISAVE(3)   
   > \endverbatim   

   > \par Contributors:   
    ==================   
   >   
   >     Nick Higham, University of Manchester   

   > \par References:   
    ================   
   >   
   >  N.J. Higham, "FORTRAN codes for estimating the one-norm of   
   >  a real or complex matrix, with applications to condition estimation",   
   >  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.   
   >   
    =====================================================================   
   Subroutine */ int igraphdlacn2_(integer *n, doublereal *v, doublereal *x, 
	integer *isgn, doublereal *est, integer *kase, integer *isave)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);
    integer i_dnnt(doublereal *);

    /* Local variables */
    integer i__;
    doublereal temp;
    extern doublereal igraphdasum_(integer *, doublereal *, integer *);
    integer jlast;
    extern /* Subroutine */ int igraphdcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern integer igraphidamax_(integer *, doublereal *, integer *);
    doublereal altsgn, estold;


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       Parameter adjustments */
    --isave;
    --isgn;
    --x;
    --v;

    /* Function Body */
    if (*kase == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] = 1. / (doublereal) (*n);
/* L10: */
	}
	*kase = 1;
	isave[1] = 1;
	return 0;
    }

    switch (isave[1]) {
	case 1:  goto L20;
	case 2:  goto L40;
	case 3:  goto L70;
	case 4:  goto L110;
	case 5:  goto L140;
    }

/*     ................ ENTRY   (ISAVE( 1 ) = 1)   
       FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

L20:
    if (*n == 1) {
	v[1] = x[1];
	*est = abs(v[1]);
/*        ... QUIT */
	goto L150;
    }
    *est = igraphdasum_(n, &x[1], &c__1);

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = d_sign(&c_b11, &x[i__]);
	isgn[i__] = i_dnnt(&x[i__]);
/* L30: */
    }
    *kase = 2;
    isave[1] = 2;
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 2)   
       FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */

L40:
    isave[2] = igraphidamax_(n, &x[1], &c__1);
    isave[3] = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

L50:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = 0.;
/* L60: */
    }
    x[isave[2]] = 1.;
    *kase = 1;
    isave[1] = 3;
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 3)   
       X HAS BEEN OVERWRITTEN BY A*X. */

L70:
    igraphdcopy_(n, &x[1], &c__1, &v[1], &c__1);
    estold = *est;
    *est = igraphdasum_(n, &v[1], &c__1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = d_sign(&c_b11, &x[i__]);
	if (i_dnnt(&d__1) != isgn[i__]) {
	    goto L90;
	}
/* L80: */
    }
/*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
    goto L120;

L90:
/*     TEST FOR CYCLING. */
    if (*est <= estold) {
	goto L120;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = d_sign(&c_b11, &x[i__]);
	isgn[i__] = i_dnnt(&x[i__]);
/* L100: */
    }
    *kase = 2;
    isave[1] = 4;
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 4)   
       X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */

L110:
    jlast = isave[2];
    isave[2] = igraphidamax_(n, &x[1], &c__1);
    if (x[jlast] != (d__1 = x[isave[2]], abs(d__1)) && isave[3] < 5) {
	++isave[3];
	goto L50;
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

L120:
    altsgn = 1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) + 
		1.);
	altsgn = -altsgn;
/* L130: */
    }
    *kase = 1;
    isave[1] = 5;
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 5)   
       X HAS BEEN OVERWRITTEN BY A*X. */

L140:
    temp = igraphdasum_(n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
    if (temp > *est) {
	igraphdcopy_(n, &x[1], &c__1, &v[1], &c__1);
	*est = temp;
    }

L150:
    *kase = 0;
    return 0;

/*     End of DLACN2 */

} /* igraphdlacn2_ */

