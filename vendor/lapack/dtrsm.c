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

/* > \brief \b DTRSM   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

    Definition:   
    ===========   

         SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)   

         DOUBLE PRECISION ALPHA   
         INTEGER LDA,LDB,M,N   
         CHARACTER DIAG,SIDE,TRANSA,UPLO   
         DOUBLE PRECISION A(LDA,*),B(LDB,*)   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DTRSM  solves one of the matrix equations   
   >   
   >    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,   
   >   
   > where alpha is a scalar, X and B are m by n matrices, A is a unit, or   
   > non-unit,  upper or lower triangular matrix  and  op( A )  is one  of   
   >   
   >    op( A ) = A   or   op( A ) = A**T.   
   >   
   > The matrix X is overwritten on B.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] SIDE   
   > \verbatim   
   >          SIDE is CHARACTER*1   
   >           On entry, SIDE specifies whether op( A ) appears on the left   
   >           or right of X as follows:   
   >   
   >              SIDE = 'L' or 'l'   op( A )*X = alpha*B.   
   >   
   >              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.   
   > \endverbatim   
   >   
   > \param[in] UPLO   
   > \verbatim   
   >          UPLO is CHARACTER*1   
   >           On entry, UPLO specifies whether the matrix A is an upper or   
   >           lower triangular matrix as follows:   
   >   
   >              UPLO = 'U' or 'u'   A is an upper triangular matrix.   
   >   
   >              UPLO = 'L' or 'l'   A is a lower triangular matrix.   
   > \endverbatim   
   >   
   > \param[in] TRANSA   
   > \verbatim   
   >          TRANSA is CHARACTER*1   
   >           On entry, TRANSA specifies the form of op( A ) to be used in   
   >           the matrix multiplication as follows:   
   >   
   >              TRANSA = 'N' or 'n'   op( A ) = A.   
   >   
   >              TRANSA = 'T' or 't'   op( A ) = A**T.   
   >   
   >              TRANSA = 'C' or 'c'   op( A ) = A**T.   
   > \endverbatim   
   >   
   > \param[in] DIAG   
   > \verbatim   
   >          DIAG is CHARACTER*1   
   >           On entry, DIAG specifies whether or not A is unit triangular   
   >           as follows:   
   >   
   >              DIAG = 'U' or 'u'   A is assumed to be unit triangular.   
   >   
   >              DIAG = 'N' or 'n'   A is not assumed to be unit   
   >                                  triangular.   
   > \endverbatim   
   >   
   > \param[in] M   
   > \verbatim   
   >          M is INTEGER   
   >           On entry, M specifies the number of rows of B. M must be at   
   >           least zero.   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >           On entry, N specifies the number of columns of B.  N must be   
   >           at least zero.   
   > \endverbatim   
   >   
   > \param[in] ALPHA   
   > \verbatim   
   >          ALPHA is DOUBLE PRECISION.   
   >           On entry,  ALPHA specifies the scalar  alpha. When  alpha is   
   >           zero then  A is not referenced and  B need not be set before   
   >           entry.   
   > \endverbatim   
   >   
   > \param[in] A   
   > \verbatim   
   >          A is DOUBLE PRECISION array, dimension ( LDA, k ),   
   >           where k is m when SIDE = 'L' or 'l'   
   >             and k is n when SIDE = 'R' or 'r'.   
   >           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k   
   >           upper triangular part of the array  A must contain the upper   
   >           triangular matrix  and the strictly lower triangular part of   
   >           A is not referenced.   
   >           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k   
   >           lower triangular part of the array  A must contain the lower   
   >           triangular matrix  and the strictly upper triangular part of   
   >           A is not referenced.   
   >           Note that when  DIAG = 'U' or 'u',  the diagonal elements of   
   >           A  are not referenced either,  but are assumed to be  unity.   
   > \endverbatim   
   >   
   > \param[in] LDA   
   > \verbatim   
   >          LDA is INTEGER   
   >           On entry, LDA specifies the first dimension of A as declared   
   >           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then   
   >           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'   
   >           then LDA must be at least max( 1, n ).   
   > \endverbatim   
   >   
   > \param[in,out] B   
   > \verbatim   
   >          B is DOUBLE PRECISION array, dimension ( LDB, N )   
   >           Before entry,  the leading  m by n part of the array  B must   
   >           contain  the  right-hand  side  matrix  B,  and  on exit  is   
   >           overwritten by the solution matrix  X.   
   > \endverbatim   
   >   
   > \param[in] LDB   
   > \verbatim   
   >          LDB is INTEGER   
   >           On entry, LDB specifies the first dimension of B as declared   
   >           in  the  calling  (sub)  program.   LDB  must  be  at  least   
   >           max( 1, m ).   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \ingroup trsm   

   > \par Further Details:   
    =====================   
   >   
   > \verbatim   
   >   
   >  Level 3 Blas routine.   
   >   
   >   
   >  -- Written on 8-February-1989.   
   >     Jack Dongarra, Argonne National Laboratory.   
   >     Iain Duff, AERE Harwell.   
   >     Jeremy Du Croz, Numerical Algorithms Group Ltd.   
   >     Sven Hammarling, Numerical Algorithms Group Ltd.   
   > \endverbatim   
   >   
    =====================================================================   
   Subroutine */ int igraphdtrsm_(char *side, char *uplo, char *transa, char *diag, 
	integer *m, integer *n, doublereal *alpha, doublereal *a, integer *
	lda, doublereal *b, integer *ldb)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__, j, k, info;
    doublereal temp;
    logical lside;
    extern logical igraphlsame_(char *, char *);
    integer nrowa;
    logical upper;
    extern /* Subroutine */ int igraphxerbla_(char *, integer *, ftnlen);
    logical nounit;


/*  -- Reference BLAS level3 routine --   
    -- Reference BLAS is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   


    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    lside = igraphlsame_(side, "L");
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = igraphlsame_(diag, "N");
    upper = igraphlsame_(uplo, "U");

    info = 0;
    if (! lside && ! igraphlsame_(side, "R")) {
	info = 1;
    } else if (! upper && ! igraphlsame_(uplo, "L")) {
	info = 2;
    } else if (! igraphlsame_(transa, "N") && ! igraphlsame_(transa,
	     "T") && ! igraphlsame_(transa, "C")) {
	info = 3;
    } else if (! igraphlsame_(diag, "U") && ! igraphlsame_(diag, 
	    "N")) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < max(1,nrowa)) {
	info = 9;
    } else if (*ldb < max(1,*m)) {
	info = 11;
    }
    if (info != 0) {
	igraphxerbla_("DTRSM ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b[i__ + j * b_dim1] = 0.;
/* L10: */
	    }
/* L20: */
	}
	return 0;
    }

/*     Start the operations. */

    if (lside) {
	if (igraphlsame_(transa, "N")) {

/*           Form  B := alpha*inv( A )*B. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (*alpha != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L30: */
			}
		    }
		    for (k = *m; k >= 1; --k) {
			if (b[k + j * b_dim1] != 0.) {
			    if (nounit) {
				b[k + j * b_dim1] /= a[k + k * a_dim1];
			    }
			    i__2 = k - 1;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
					i__ + k * a_dim1];
/* L40: */
			    }
			}
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (*alpha != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L70: */
			}
		    }
		    i__2 = *m;
		    for (k = 1; k <= i__2; ++k) {
			if (b[k + j * b_dim1] != 0.) {
			    if (nounit) {
				b[k + j * b_dim1] /= a[k + k * a_dim1];
			    }
			    i__3 = *m;
			    for (i__ = k + 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
					i__ + k * a_dim1];
/* L80: */
			    }
			}
/* L90: */
		    }
/* L100: */
		}
	    }
	} else {

/*           Form  B := alpha*inv( A**T )*B. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = *alpha * b[i__ + j * b_dim1];
			i__3 = i__ - 1;
			for (k = 1; k <= i__3; ++k) {
			    temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
/* L110: */
			}
			if (nounit) {
			    temp /= a[i__ + i__ * a_dim1];
			}
			b[i__ + j * b_dim1] = temp;
/* L120: */
		    }
/* L130: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = *m; i__ >= 1; --i__) {
			temp = *alpha * b[i__ + j * b_dim1];
			i__2 = *m;
			for (k = i__ + 1; k <= i__2; ++k) {
			    temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
/* L140: */
			}
			if (nounit) {
			    temp /= a[i__ + i__ * a_dim1];
			}
			b[i__ + j * b_dim1] = temp;
/* L150: */
		    }
/* L160: */
		}
	    }
	}
    } else {
	if (igraphlsame_(transa, "N")) {

/*           Form  B := alpha*B*inv( A ). */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (*alpha != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L170: */
			}
		    }
		    i__2 = j - 1;
		    for (k = 1; k <= i__2; ++k) {
			if (a[k + j * a_dim1] != 0.) {
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
					i__ + k * b_dim1];
/* L180: */
			    }
			}
/* L190: */
		    }
		    if (nounit) {
			temp = 1. / a[j + j * a_dim1];
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
/* L200: */
			}
		    }
/* L210: */
		}
	    } else {
		for (j = *n; j >= 1; --j) {
		    if (*alpha != 1.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L220: */
			}
		    }
		    i__1 = *n;
		    for (k = j + 1; k <= i__1; ++k) {
			if (a[k + j * a_dim1] != 0.) {
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
					i__ + k * b_dim1];
/* L230: */
			    }
			}
/* L240: */
		    }
		    if (nounit) {
			temp = 1. / a[j + j * a_dim1];
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
/* L250: */
			}
		    }
/* L260: */
		}
	    }
	} else {

/*           Form  B := alpha*B*inv( A**T ). */

	    if (upper) {
		for (k = *n; k >= 1; --k) {
		    if (nounit) {
			temp = 1. / a[k + k * a_dim1];
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
/* L270: */
			}
		    }
		    i__1 = k - 1;
		    for (j = 1; j <= i__1; ++j) {
			if (a[j + k * a_dim1] != 0.) {
			    temp = a[j + k * a_dim1];
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] -= temp * b[i__ + k * 
					b_dim1];
/* L280: */
			    }
			}
/* L290: */
		    }
		    if (*alpha != 1.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
				    ;
/* L300: */
			}
		    }
/* L310: */
		}
	    } else {
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
		    if (nounit) {
			temp = 1. / a[k + k * a_dim1];
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
/* L320: */
			}
		    }
		    i__2 = *n;
		    for (j = k + 1; j <= i__2; ++j) {
			if (a[j + k * a_dim1] != 0.) {
			    temp = a[j + k * a_dim1];
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] -= temp * b[i__ + k * 
					b_dim1];
/* L330: */
			    }
			}
/* L340: */
		    }
		    if (*alpha != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
				    ;
/* L350: */
			}
		    }
/* L360: */
		}
	    }
	}
    }

    return 0;

/*     End of DTRSM */

} /* igraphdtrsm_ */

