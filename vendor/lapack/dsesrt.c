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

/* -----------------------------------------------------------------------   
   \BeginDoc   

   \Name: dsesrt   

   \Description:   
    Sort the array X in the order specified by WHICH and optionally   
    apply the permutation to the columns of the matrix A.   

   \Usage:   
    call dsesrt   
       ( WHICH, APPLY, N, X, NA, A, LDA)   

   \Arguments   
    WHICH   Character*2.  (Input)   
            'LM' -> X is sorted into increasing order of magnitude.   
            'SM' -> X is sorted into decreasing order of magnitude.   
            'LA' -> X is sorted into increasing order of algebraic.   
            'SA' -> X is sorted into decreasing order of algebraic.   

    APPLY   Logical.  (Input)   
            APPLY = .TRUE.  -> apply the sorted order to A.   
            APPLY = .FALSE. -> do not apply the sorted order to A.   

    N       Integer.  (INPUT)   
            Dimension of the array X.   

    X      Double precision array of length N.  (INPUT/OUTPUT)   
            The array to be sorted.   

    NA      Integer.  (INPUT)   
            Number of rows of the matrix A.   

    A      Double precision array of length NA by N.  (INPUT/OUTPUT)   

    LDA     Integer.  (INPUT)   
            Leading dimension of A.   

   \EndDoc   

   -----------------------------------------------------------------------   

   \BeginLib   

   \Routines   
       dswap  Level 1 BLAS that swaps the contents of two vectors.   

   \Authors   
       Danny Sorensen               Phuong Vu   
       Richard Lehoucq              CRPC / Rice University   
       Dept. of Computational &     Houston, Texas   
       Applied Mathematics   
       Rice University   
       Houston, Texas   

   \Revision history:   
       12/15/93: Version ' 2.1'.   
                 Adapted from the sort routine in LANSO and   
                 the ARPACK code dsortr   

   \SCCS Information: @(#)   
   FILE: sesrt.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2   

   \EndLib   

   -----------------------------------------------------------------------   

   Subroutine */ int igraphdsesrt_(char *which, logical *apply, integer *n, 
	doublereal *x, integer *na, doublereal *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    integer i__, j, igap;
    doublereal temp;
    extern /* Subroutine */ int igraphdswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);


/*     %------------------%   
       | Scalar Arguments |   
       %------------------%   


       %-----------------%   
       | Array Arguments |   
       %-----------------%   


       %---------------%   
       | Local Scalars |   
       %---------------%   


       %----------------------%   
       | External Subroutines |   
       %----------------------%   


       %-----------------------%   
       | Executable Statements |   
       %-----------------------%   

       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 0;
    a -= a_offset;

    /* Function Body */
    igap = *n / 2;

    if (s_cmp(which, "SA", (ftnlen)2, (ftnlen)2) == 0) {

/*        X is sorted into decreasing order of algebraic. */

L10:
	if (igap == 0) {
	    goto L9000;
	}
	i__1 = *n - 1;
	for (i__ = igap; i__ <= i__1; ++i__) {
	    j = i__ - igap;
L20:

	    if (j < 0) {
		goto L30;
	    }

	    if (x[j] < x[j + igap]) {
		temp = x[j];
		x[j] = x[j + igap];
		x[j + igap] = temp;
		if (*apply) {
		    igraphdswap_(na, &a[j * a_dim1 + 1], &c__1, &a[(j + igap) * 
			    a_dim1 + 1], &c__1);
		}
	    } else {
		goto L30;
	    }
	    j -= igap;
	    goto L20;
L30:
	    ;
	}
	igap /= 2;
	goto L10;

    } else if (s_cmp(which, "SM", (ftnlen)2, (ftnlen)2) == 0) {

/*        X is sorted into decreasing order of magnitude. */

L40:
	if (igap == 0) {
	    goto L9000;
	}
	i__1 = *n - 1;
	for (i__ = igap; i__ <= i__1; ++i__) {
	    j = i__ - igap;
L50:

	    if (j < 0) {
		goto L60;
	    }

	    if ((d__1 = x[j], abs(d__1)) < (d__2 = x[j + igap], abs(d__2))) {
		temp = x[j];
		x[j] = x[j + igap];
		x[j + igap] = temp;
		if (*apply) {
		    igraphdswap_(na, &a[j * a_dim1 + 1], &c__1, &a[(j + igap) * 
			    a_dim1 + 1], &c__1);
		}
	    } else {
		goto L60;
	    }
	    j -= igap;
	    goto L50;
L60:
	    ;
	}
	igap /= 2;
	goto L40;

    } else if (s_cmp(which, "LA", (ftnlen)2, (ftnlen)2) == 0) {

/*        X is sorted into increasing order of algebraic. */

L70:
	if (igap == 0) {
	    goto L9000;
	}
	i__1 = *n - 1;
	for (i__ = igap; i__ <= i__1; ++i__) {
	    j = i__ - igap;
L80:

	    if (j < 0) {
		goto L90;
	    }

	    if (x[j] > x[j + igap]) {
		temp = x[j];
		x[j] = x[j + igap];
		x[j + igap] = temp;
		if (*apply) {
		    igraphdswap_(na, &a[j * a_dim1 + 1], &c__1, &a[(j + igap) * 
			    a_dim1 + 1], &c__1);
		}
	    } else {
		goto L90;
	    }
	    j -= igap;
	    goto L80;
L90:
	    ;
	}
	igap /= 2;
	goto L70;

    } else if (s_cmp(which, "LM", (ftnlen)2, (ftnlen)2) == 0) {

/*        X is sorted into increasing order of magnitude. */

L100:
	if (igap == 0) {
	    goto L9000;
	}
	i__1 = *n - 1;
	for (i__ = igap; i__ <= i__1; ++i__) {
	    j = i__ - igap;
L110:

	    if (j < 0) {
		goto L120;
	    }

	    if ((d__1 = x[j], abs(d__1)) > (d__2 = x[j + igap], abs(d__2))) {
		temp = x[j];
		x[j] = x[j + igap];
		x[j + igap] = temp;
		if (*apply) {
		    igraphdswap_(na, &a[j * a_dim1 + 1], &c__1, &a[(j + igap) * 
			    a_dim1 + 1], &c__1);
		}
	    } else {
		goto L120;
	    }
	    j -= igap;
	    goto L110;
L120:
	    ;
	}
	igap /= 2;
	goto L100;
    }

L9000:
    return 0;

/*     %---------------%   
       | End of dsesrt |   
       %---------------% */

} /* igraphdsesrt_ */

