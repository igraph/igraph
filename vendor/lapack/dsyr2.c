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

/* > \brief \b DSYR2   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

    Definition:   
    ===========   

         SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)   

         DOUBLE PRECISION ALPHA   
         INTEGER INCX,INCY,LDA,N   
         CHARACTER UPLO   
         DOUBLE PRECISION A(LDA,*),X(*),Y(*)   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DSYR2  performs the symmetric rank 2 operation   
   >   
   >    A := alpha*x*y**T + alpha*y*x**T + A,   
   >   
   > where alpha is a scalar, x and y are n element vectors and A is an n   
   > by n symmetric matrix.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] UPLO   
   > \verbatim   
   >          UPLO is CHARACTER*1   
   >           On entry, UPLO specifies whether the upper or lower   
   >           triangular part of the array A is to be referenced as   
   >           follows:   
   >   
   >              UPLO = 'U' or 'u'   Only the upper triangular part of A   
   >                                  is to be referenced.   
   >   
   >              UPLO = 'L' or 'l'   Only the lower triangular part of A   
   >                                  is to be referenced.   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >           On entry, N specifies the order of the matrix A.   
   >           N must be at least zero.   
   > \endverbatim   
   >   
   > \param[in] ALPHA   
   > \verbatim   
   >          ALPHA is DOUBLE PRECISION.   
   >           On entry, ALPHA specifies the scalar alpha.   
   > \endverbatim   
   >   
   > \param[in] X   
   > \verbatim   
   >          X is DOUBLE PRECISION array, dimension at least   
   >           ( 1 + ( n - 1 )*abs( INCX ) ).   
   >           Before entry, the incremented array X must contain the n   
   >           element vector x.   
   > \endverbatim   
   >   
   > \param[in] INCX   
   > \verbatim   
   >          INCX is INTEGER   
   >           On entry, INCX specifies the increment for the elements of   
   >           X. INCX must not be zero.   
   > \endverbatim   
   >   
   > \param[in] Y   
   > \verbatim   
   >          Y is DOUBLE PRECISION array, dimension at least   
   >           ( 1 + ( n - 1 )*abs( INCY ) ).   
   >           Before entry, the incremented array Y must contain the n   
   >           element vector y.   
   > \endverbatim   
   >   
   > \param[in] INCY   
   > \verbatim   
   >          INCY is INTEGER   
   >           On entry, INCY specifies the increment for the elements of   
   >           Y. INCY must not be zero.   
   > \endverbatim   
   >   
   > \param[in,out] A   
   > \verbatim   
   >          A is DOUBLE PRECISION array, dimension ( LDA, N )   
   >           Before entry with  UPLO = 'U' or 'u', the leading n by n   
   >           upper triangular part of the array A must contain the upper   
   >           triangular part of the symmetric matrix and the strictly   
   >           lower triangular part of A is not referenced. On exit, the   
   >           upper triangular part of the array A is overwritten by the   
   >           upper triangular part of the updated matrix.   
   >           Before entry with UPLO = 'L' or 'l', the leading n by n   
   >           lower triangular part of the array A must contain the lower   
   >           triangular part of the symmetric matrix and the strictly   
   >           upper triangular part of A is not referenced. On exit, the   
   >           lower triangular part of the array A is overwritten by the   
   >           lower triangular part of the updated matrix.   
   > \endverbatim   
   >   
   > \param[in] LDA   
   > \verbatim   
   >          LDA is INTEGER   
   >           On entry, LDA specifies the first dimension of A as declared   
   >           in the calling (sub) program. LDA must be at least   
   >           max( 1, n ).   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \ingroup her2   

   > \par Further Details:   
    =====================   
   >   
   > \verbatim   
   >   
   >  Level 2 Blas routine.   
   >   
   >  -- Written on 22-October-1986.   
   >     Jack Dongarra, Argonne National Lab.   
   >     Jeremy Du Croz, Nag Central Office.   
   >     Sven Hammarling, Nag Central Office.   
   >     Richard Hanson, Sandia National Labs.   
   > \endverbatim   
   >   
    =====================================================================   
   Subroutine */ int igraphdsyr2_(char *uplo, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *y, integer *incy, 
	doublereal *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, ix, iy, jx, jy, kx, ky, info;
    doublereal temp1, temp2;
    extern logical igraphlsame_(char *, char *);
    extern /* Subroutine */ int igraphxerbla_(char *, integer *, ftnlen);


/*  -- Reference BLAS level2 routine --   
    -- Reference BLAS is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   


    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (! igraphlsame_(uplo, "U") && ! igraphlsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*n)) {
	info = 9;
    }
    if (info != 0) {
	igraphxerbla_("DSYR2 ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0.) {
	return 0;
    }

/*     Set up the start points in X and Y if the increments are not both   
       unity. */

    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A. */

    if (igraphlsame_(uplo, "U")) {

/*        Form  A  when A is stored in the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0. || y[j] != 0.) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
				temp1 + y[i__] * temp2;
/* L10: */
		    }
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0. || y[jy] != 0.) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = kx;
		    iy = ky;
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
				temp1 + y[iy] * temp2;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when A is stored in the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0. || y[j] != 0.) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
				temp1 + y[i__] * temp2;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0. || y[jy] != 0.) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = jx;
		    iy = jy;
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
				temp1 + y[iy] * temp2;
			ix += *incx;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    }

    return 0;

/*     End of DSYR2 */

} /* igraphdsyr2_ */

