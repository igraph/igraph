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

/* > \brief \b DLAGTF computes an LU factorization of a matrix T-λI, where T is a general tridiagonal matrix,
 and λ a scalar, using partial pivoting with row interchanges.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLAGTF + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlagtf.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlagtf.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlagtf.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DLAGTF( N, A, LAMBDA, B, C, TOL, D, IN, INFO )   

         INTEGER            INFO, N   
         DOUBLE PRECISION   LAMBDA, TOL   
         INTEGER            IN( * )   
         DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DLAGTF factorizes the matrix (T - lambda*I), where T is an n by n   
   > tridiagonal matrix and lambda is a scalar, as   
   >   
   >    T - lambda*I = PLU,   
   >   
   > where P is a permutation matrix, L is a unit lower tridiagonal matrix   
   > with at most one non-zero sub-diagonal elements per column and U is   
   > an upper triangular matrix with at most two non-zero super-diagonal   
   > elements per column.   
   >   
   > The factorization is obtained by Gaussian elimination with partial   
   > pivoting and implicit row scaling.   
   >   
   > The parameter LAMBDA is included in the routine so that DLAGTF may   
   > be used, in conjunction with DLAGTS, to obtain eigenvectors of T by   
   > inverse iteration.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix T.   
   > \endverbatim   
   >   
   > \param[in,out] A   
   > \verbatim   
   >          A is DOUBLE PRECISION array, dimension (N)   
   >          On entry, A must contain the diagonal elements of T.   
   >   
   >          On exit, A is overwritten by the n diagonal elements of the   
   >          upper triangular matrix U of the factorization of T.   
   > \endverbatim   
   >   
   > \param[in] LAMBDA   
   > \verbatim   
   >          LAMBDA is DOUBLE PRECISION   
   >          On entry, the scalar lambda.   
   > \endverbatim   
   >   
   > \param[in,out] B   
   > \verbatim   
   >          B is DOUBLE PRECISION array, dimension (N-1)   
   >          On entry, B must contain the (n-1) super-diagonal elements of   
   >          T.   
   >   
   >          On exit, B is overwritten by the (n-1) super-diagonal   
   >          elements of the matrix U of the factorization of T.   
   > \endverbatim   
   >   
   > \param[in,out] C   
   > \verbatim   
   >          C is DOUBLE PRECISION array, dimension (N-1)   
   >          On entry, C must contain the (n-1) sub-diagonal elements of   
   >          T.   
   >   
   >          On exit, C is overwritten by the (n-1) sub-diagonal elements   
   >          of the matrix L of the factorization of T.   
   > \endverbatim   
   >   
   > \param[in] TOL   
   > \verbatim   
   >          TOL is DOUBLE PRECISION   
   >          On entry, a relative tolerance used to indicate whether or   
   >          not the matrix (T - lambda*I) is nearly singular. TOL should   
   >          normally be chose as approximately the largest relative error   
   >          in the elements of T. For example, if the elements of T are   
   >          correct to about 4 significant figures, then TOL should be   
   >          set to about 5*10**(-4). If TOL is supplied as less than eps,   
   >          where eps is the relative machine precision, then the value   
   >          eps is used in place of TOL.   
   > \endverbatim   
   >   
   > \param[out] D   
   > \verbatim   
   >          D is DOUBLE PRECISION array, dimension (N-2)   
   >          On exit, D is overwritten by the (n-2) second super-diagonal   
   >          elements of the matrix U of the factorization of T.   
   > \endverbatim   
   >   
   > \param[out] IN   
   > \verbatim   
   >          IN is INTEGER array, dimension (N)   
   >          On exit, IN contains details of the permutation matrix P. If   
   >          an interchange occurred at the kth step of the elimination,   
   >          then IN(k) = 1, otherwise IN(k) = 0. The element IN(n)   
   >          returns the smallest positive integer j such that   
   >   
   >             abs( u(j,j) ).le. norm( (T - lambda*I)(j) )*TOL,   
   >   
   >          where norm( A(j) ) denotes the sum of the absolute values of   
   >          the jth row of the matrix A. If no such j exists then IN(n)   
   >          is returned as zero. If IN(n) is returned as positive, then a   
   >          diagonal element of U is small, indicating that   
   >          (T - lambda*I) is singular or nearly singular,   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
   >          = 0   : successful exit   
   >          .lt. 0: if INFO = -k, the kth argument had an illegal value   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup auxOTHERcomputational   

    =====================================================================   
   Subroutine */ int igraphdlagtf_(integer *n, doublereal *a, doublereal *lambda, 
	doublereal *b, doublereal *c__, doublereal *tol, doublereal *d__, 
	integer *in, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    integer k;
    doublereal tl, eps, piv1, piv2, temp, mult, scale1, scale2;
    extern doublereal igraphdlamch_(char *);
    extern /* Subroutine */ int igraphxerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


   =====================================================================   


       Parameter adjustments */
    --in;
    --d__;
    --c__;
    --b;
    --a;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	igraphxerbla_("DLAGTF", &i__1, (ftnlen)6);
	return 0;
    }

    if (*n == 0) {
	return 0;
    }

    a[1] -= *lambda;
    in[*n] = 0;
    if (*n == 1) {
	if (a[1] == 0.) {
	    in[1] = 1;
	}
	return 0;
    }

    eps = igraphdlamch_("Epsilon");

    tl = max(*tol,eps);
    scale1 = abs(a[1]) + abs(b[1]);
    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	a[k + 1] -= *lambda;
	scale2 = (d__1 = c__[k], abs(d__1)) + (d__2 = a[k + 1], abs(d__2));
	if (k < *n - 1) {
	    scale2 += (d__1 = b[k + 1], abs(d__1));
	}
	if (a[k] == 0.) {
	    piv1 = 0.;
	} else {
	    piv1 = (d__1 = a[k], abs(d__1)) / scale1;
	}
	if (c__[k] == 0.) {
	    in[k] = 0;
	    piv2 = 0.;
	    scale1 = scale2;
	    if (k < *n - 1) {
		d__[k] = 0.;
	    }
	} else {
	    piv2 = (d__1 = c__[k], abs(d__1)) / scale2;
	    if (piv2 <= piv1) {
		in[k] = 0;
		scale1 = scale2;
		c__[k] /= a[k];
		a[k + 1] -= c__[k] * b[k];
		if (k < *n - 1) {
		    d__[k] = 0.;
		}
	    } else {
		in[k] = 1;
		mult = a[k] / c__[k];
		a[k] = c__[k];
		temp = a[k + 1];
		a[k + 1] = b[k] - mult * temp;
		if (k < *n - 1) {
		    d__[k] = b[k + 1];
		    b[k + 1] = -mult * d__[k];
		}
		b[k] = temp;
		c__[k] = mult;
	    }
	}
	if (max(piv1,piv2) <= tl && in[*n] == 0) {
	    in[*n] = k;
	}
/* L10: */
    }
    if ((d__1 = a[*n], abs(d__1)) <= scale1 * tl && in[*n] == 0) {
	in[*n] = *n;
    }

    return 0;

/*     End of DLAGTF */

} /* igraphdlagtf_ */

