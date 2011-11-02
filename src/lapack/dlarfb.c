/*  -- translated by f2c (version 20100827).
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
static doublereal c_b14 = 1.;
static doublereal c_b25 = -1.;

/* Subroutine */ int igraphdlarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, doublereal *v, integer *
	ldv, doublereal *t, integer *ldt, doublereal *c__, integer *ldc, 
	doublereal *work, integer *ldwork)
{
    /* System generated locals */
    integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
	    work_offset, i__1, i__2;

    /* Local variables */
    integer i__, j;
    extern /* Subroutine */ int igraphdgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);
    extern logical igraphlsame_(char *, char *);
    integer lastc;
    extern /* Subroutine */ int igraphdcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), igraphdtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *);
    integer lastv;
    extern integer igraphiladlc_(integer *, integer *, doublereal *, integer *), 
	    igraphiladlr_(integer *, integer *, doublereal *, integer *);
    char transt[1];


/*  -- LAPACK auxiliary routine (version 3.3.1) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
    -- April 2011                                                      --   


    Purpose   
    =======   

    DLARFB applies a real block reflector H or its transpose H**T to a   
    real m by n matrix C, from either the left or the right.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply H or H**T from the Left   
            = 'R': apply H or H**T from the Right   

    TRANS   (input) CHARACTER*1   
            = 'N': apply H (No transpose)   
            = 'T': apply H**T (Transpose)   

    DIRECT  (input) CHARACTER*1   
            Indicates how H is formed from a product of elementary   
            reflectors   
            = 'F': H = H(1) H(2) . . . H(k) (Forward)   
            = 'B': H = H(k) . . . H(2) H(1) (Backward)   

    STOREV  (input) CHARACTER*1   
            Indicates how the vectors which define the elementary   
            reflectors are stored:   
            = 'C': Columnwise   
            = 'R': Rowwise   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    K       (input) INTEGER   
            The order of the matrix T (= the number of elementary   
            reflectors whose product defines the block reflector).   

    V       (input) DOUBLE PRECISION array, dimension   
                                  (LDV,K) if STOREV = 'C'   
                                  (LDV,M) if STOREV = 'R' and SIDE = 'L'   
                                  (LDV,N) if STOREV = 'R' and SIDE = 'R'   
            The matrix V. See Further Details.   

    LDV     (input) INTEGER   
            The leading dimension of the array V.   
            If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);   
            if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);   
            if STOREV = 'R', LDV >= K.   

    T       (input) DOUBLE PRECISION array, dimension (LDT,K)   
            The triangular k by k matrix T in the representation of the   
            block reflector.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= K.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= max(1,M).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)   

    LDWORK  (input) INTEGER   
            The leading dimension of the array WORK.   
            If SIDE = 'L', LDWORK >= max(1,N);   
            if SIDE = 'R', LDWORK >= max(1,M).   

    Further Details   
    ===============   

    The shape of the matrix V and the storage of the vectors which define   
    the H(i) is best illustrated by the following example with n = 5 and   
    k = 3. The elements equal to 1 are not stored; the corresponding   
    array elements are modified but restored on exit. The rest of the   
    array is not used.   

    DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':   

                 V = (  1       )                 V = (  1 v1 v1 v1 v1 )   
                     ( v1  1    )                     (     1 v2 v2 v2 )   
                     ( v1 v2  1 )                     (        1 v3 v3 )   
                     ( v1 v2 v3 )   
                     ( v1 v2 v3 )   

    DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':   

                 V = ( v1 v2 v3 )                 V = ( v1 v1  1       )   
                     ( v1 v2 v3 )                     ( v2 v2 v2  1    )   
                     (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )   
                     (     1 v3 )   
                     (        1 )   

    =====================================================================   


       Quick return if possible   

       Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    work_dim1 = *ldwork;
    work_offset = 1 + work_dim1;
    work -= work_offset;

    /* Function Body */
    if (*m <= 0 || *n <= 0) {
	return 0;
    }

    if (igraphlsame_(trans, "N")) {
	*(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transt = 'N';
    }

    if (igraphlsame_(storev, "C")) {

	if (igraphlsame_(direct, "F")) {

/*           Let  V =  ( V1 )    (first K rows)   
                       ( V2 )   
             where  V1  is unit lower triangular. */

	    if (igraphlsame_(side, "L")) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 )   
                                                      ( C2 )   

   Computing MAX */
		i__1 = *k, i__2 = igraphiladlr_(m, k, &v[v_offset], ldv);
		lastv = max(i__1,i__2);
		lastc = igraphiladlc_(&lastv, n, &c__[c_offset], ldc);

/*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)   

                W := C1**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    igraphdcopy_(&lastc, &c__[j + c_dim1], ldc, &work[j * work_dim1 
			    + 1], &c__1);
/* L10: */
		}

/*              W := W * V1 */

		igraphdtrmm_("Right", "Lower", "No transpose", "Unit", &lastc, k, &
			c_b14, &v[v_offset], ldv, &work[work_offset], ldwork);
		if (lastv > *k) {

/*                 W := W + C2**T *V2 */

		    i__1 = lastv - *k;
		    igraphdgemm_("Transpose", "No transpose", &lastc, k, &i__1, &
			    c_b14, &c__[*k + 1 + c_dim1], ldc, &v[*k + 1 + 
			    v_dim1], ldv, &c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T**T  or  W * T */

		igraphdtrmm_("Right", "Upper", transt, "Non-unit", &lastc, k, &
			c_b14, &t[t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - V * W**T */

		if (lastv > *k) {

/*                 C2 := C2 - V2 * W**T */

		    i__1 = lastv - *k;
		    igraphdgemm_("No transpose", "Transpose", &i__1, &lastc, k, &
			    c_b25, &v[*k + 1 + v_dim1], ldv, &work[
			    work_offset], ldwork, &c_b14, &c__[*k + 1 + 
			    c_dim1], ldc);
		}

/*              W := W * V1**T */

		igraphdtrmm_("Right", "Lower", "Transpose", "Unit", &lastc, k, &
			c_b14, &v[v_offset], ldv, &work[work_offset], ldwork);

/*              C1 := C1 - W**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = lastc;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
/* L20: */
		    }
/* L30: */
		}

	    } else if (igraphlsame_(side, "R")) {

/*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )   

   Computing MAX */
		i__1 = *k, i__2 = igraphiladlr_(n, k, &v[v_offset], ldv);
		lastv = max(i__1,i__2);
		lastc = igraphiladlr_(m, &lastv, &c__[c_offset], ldc);

/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)   

                W := C1 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    igraphdcopy_(&lastc, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
/* L40: */
		}

/*              W := W * V1 */

		igraphdtrmm_("Right", "Lower", "No transpose", "Unit", &lastc, k, &
			c_b14, &v[v_offset], ldv, &work[work_offset], ldwork);
		if (lastv > *k) {

/*                 W := W + C2 * V2 */

		    i__1 = lastv - *k;
		    igraphdgemm_("No transpose", "No transpose", &lastc, k, &i__1, &
			    c_b14, &c__[(*k + 1) * c_dim1 + 1], ldc, &v[*k + 
			    1 + v_dim1], ldv, &c_b14, &work[work_offset], 
			    ldwork);
		}

/*              W := W * T  or  W * T**T */

		igraphdtrmm_("Right", "Upper", trans, "Non-unit", &lastc, k, &c_b14,
			 &t[t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - W * V**T */

		if (lastv > *k) {

/*                 C2 := C2 - W * V2**T */

		    i__1 = lastv - *k;
		    igraphdgemm_("No transpose", "Transpose", &lastc, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v[*k + 1 + 
			    v_dim1], ldv, &c_b14, &c__[(*k + 1) * c_dim1 + 1],
			     ldc);
		}

/*              W := W * V1**T */

		igraphdtrmm_("Right", "Lower", "Transpose", "Unit", &lastc, k, &
			c_b14, &v[v_offset], ldv, &work[work_offset], ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = lastc;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
/* L50: */
		    }
/* L60: */
		}
	    }

	} else {

/*           Let  V =  ( V1 )   
                       ( V2 )    (last K rows)   
             where  V2  is unit upper triangular. */

	    if (igraphlsame_(side, "L")) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 )   
                                                      ( C2 )   

   Computing MAX */
		i__1 = *k, i__2 = igraphiladlr_(m, k, &v[v_offset], ldv);
		lastv = max(i__1,i__2);
		lastc = igraphiladlc_(&lastv, n, &c__[c_offset], ldc);

/*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)   

                W := C2**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    igraphdcopy_(&lastc, &c__[lastv - *k + j + c_dim1], ldc, &work[
			    j * work_dim1 + 1], &c__1);
/* L70: */
		}

/*              W := W * V2 */

		igraphdtrmm_("Right", "Upper", "No transpose", "Unit", &lastc, k, &
			c_b14, &v[lastv - *k + 1 + v_dim1], ldv, &work[
			work_offset], ldwork);
		if (lastv > *k) {

/*                 W := W + C1**T*V1 */

		    i__1 = lastv - *k;
		    igraphdgemm_("Transpose", "No transpose", &lastc, k, &i__1, &
			    c_b14, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T**T  or  W * T */

		igraphdtrmm_("Right", "Lower", transt, "Non-unit", &lastc, k, &
			c_b14, &t[t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - V * W**T */

		if (lastv > *k) {

/*                 C1 := C1 - V1 * W**T */

		    i__1 = lastv - *k;
		    igraphdgemm_("No transpose", "Transpose", &i__1, &lastc, k, &
			    c_b25, &v[v_offset], ldv, &work[work_offset], 
			    ldwork, &c_b14, &c__[c_offset], ldc);
		}

/*              W := W * V2**T */

		igraphdtrmm_("Right", "Upper", "Transpose", "Unit", &lastc, k, &
			c_b14, &v[lastv - *k + 1 + v_dim1], ldv, &work[
			work_offset], ldwork);

/*              C2 := C2 - W**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = lastc;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[lastv - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
/* L80: */
		    }
/* L90: */
		}

	    } else if (igraphlsame_(side, "R")) {

/*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )   

   Computing MAX */
		i__1 = *k, i__2 = igraphiladlr_(n, k, &v[v_offset], ldv);
		lastv = max(i__1,i__2);
		lastc = igraphiladlr_(m, &lastv, &c__[c_offset], ldc);

/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)   

                W := C2 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    igraphdcopy_(&lastc, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &
			    work[j * work_dim1 + 1], &c__1);
/* L100: */
		}

/*              W := W * V2 */

		igraphdtrmm_("Right", "Upper", "No transpose", "Unit", &lastc, k, &
			c_b14, &v[lastv - *k + 1 + v_dim1], ldv, &work[
			work_offset], ldwork);
		if (lastv > *k) {

/*                 W := W + C1 * V1 */

		    i__1 = lastv - *k;
		    igraphdgemm_("No transpose", "No transpose", &lastc, k, &i__1, &
			    c_b14, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T  or  W * T**T */

		igraphdtrmm_("Right", "Lower", trans, "Non-unit", &lastc, k, &c_b14,
			 &t[t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - W * V**T */

		if (lastv > *k) {

/*                 C1 := C1 - W * V1**T */

		    i__1 = lastv - *k;
		    igraphdgemm_("No transpose", "Transpose", &lastc, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v[v_offset], 
			    ldv, &c_b14, &c__[c_offset], ldc);
		}

/*              W := W * V2**T */

		igraphdtrmm_("Right", "Upper", "Transpose", "Unit", &lastc, k, &
			c_b14, &v[lastv - *k + 1 + v_dim1], ldv, &work[
			work_offset], ldwork);

/*              C2 := C2 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = lastc;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + (lastv - *k + j) * c_dim1] -= work[i__ + j *
				 work_dim1];
/* L110: */
		    }
/* L120: */
		}
	    }
	}

    } else if (igraphlsame_(storev, "R")) {

	if (igraphlsame_(direct, "F")) {

/*           Let  V =  ( V1  V2 )    (V1: first K columns)   
             where  V1  is unit upper triangular. */

	    if (igraphlsame_(side, "L")) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 )   
                                                      ( C2 )   

   Computing MAX */
		i__1 = *k, i__2 = igraphiladlc_(k, m, &v[v_offset], ldv);
		lastv = max(i__1,i__2);
		lastc = igraphiladlc_(&lastv, n, &c__[c_offset], ldc);

/*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)   

                W := C1**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    igraphdcopy_(&lastc, &c__[j + c_dim1], ldc, &work[j * work_dim1 
			    + 1], &c__1);
/* L130: */
		}

/*              W := W * V1**T */

		igraphdtrmm_("Right", "Upper", "Transpose", "Unit", &lastc, k, &
			c_b14, &v[v_offset], ldv, &work[work_offset], ldwork);
		if (lastv > *k) {

/*                 W := W + C2**T*V2**T */

		    i__1 = lastv - *k;
		    igraphdgemm_("Transpose", "Transpose", &lastc, k, &i__1, &c_b14,
			     &c__[*k + 1 + c_dim1], ldc, &v[(*k + 1) * v_dim1 
			    + 1], ldv, &c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T**T  or  W * T */

		igraphdtrmm_("Right", "Upper", transt, "Non-unit", &lastc, k, &
			c_b14, &t[t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - V**T * W**T */

		if (lastv > *k) {

/*                 C2 := C2 - V2**T * W**T */

		    i__1 = lastv - *k;
		    igraphdgemm_("Transpose", "Transpose", &i__1, &lastc, k, &c_b25,
			     &v[(*k + 1) * v_dim1 + 1], ldv, &work[
			    work_offset], ldwork, &c_b14, &c__[*k + 1 + 
			    c_dim1], ldc);
		}

/*              W := W * V1 */

		igraphdtrmm_("Right", "Upper", "No transpose", "Unit", &lastc, k, &
			c_b14, &v[v_offset], ldv, &work[work_offset], ldwork);

/*              C1 := C1 - W**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = lastc;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
/* L140: */
		    }
/* L150: */
		}

	    } else if (igraphlsame_(side, "R")) {

/*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )   

   Computing MAX */
		i__1 = *k, i__2 = igraphiladlc_(k, n, &v[v_offset], ldv);
		lastv = max(i__1,i__2);
		lastc = igraphiladlr_(m, &lastv, &c__[c_offset], ldc);

/*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)   

                W := C1 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    igraphdcopy_(&lastc, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
/* L160: */
		}

/*              W := W * V1**T */

		igraphdtrmm_("Right", "Upper", "Transpose", "Unit", &lastc, k, &
			c_b14, &v[v_offset], ldv, &work[work_offset], ldwork);
		if (lastv > *k) {

/*                 W := W + C2 * V2**T */

		    i__1 = lastv - *k;
		    igraphdgemm_("No transpose", "Transpose", &lastc, k, &i__1, &
			    c_b14, &c__[(*k + 1) * c_dim1 + 1], ldc, &v[(*k + 
			    1) * v_dim1 + 1], ldv, &c_b14, &work[work_offset],
			     ldwork);
		}

/*              W := W * T  or  W * T**T */

		igraphdtrmm_("Right", "Upper", trans, "Non-unit", &lastc, k, &c_b14,
			 &t[t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - W * V */

		if (lastv > *k) {

/*                 C2 := C2 - W * V2 */

		    i__1 = lastv - *k;
		    igraphdgemm_("No transpose", "No transpose", &lastc, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v[(*k + 1) * 
			    v_dim1 + 1], ldv, &c_b14, &c__[(*k + 1) * c_dim1 
			    + 1], ldc);
		}

/*              W := W * V1 */

		igraphdtrmm_("Right", "Upper", "No transpose", "Unit", &lastc, k, &
			c_b14, &v[v_offset], ldv, &work[work_offset], ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = lastc;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
/* L170: */
		    }
/* L180: */
		}

	    }

	} else {

/*           Let  V =  ( V1  V2 )    (V2: last K columns)   
             where  V2  is unit lower triangular. */

	    if (igraphlsame_(side, "L")) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 )   
                                                      ( C2 )   

   Computing MAX */
		i__1 = *k, i__2 = igraphiladlc_(k, m, &v[v_offset], ldv);
		lastv = max(i__1,i__2);
		lastc = igraphiladlc_(&lastv, n, &c__[c_offset], ldc);

/*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)   

                W := C2**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    igraphdcopy_(&lastc, &c__[lastv - *k + j + c_dim1], ldc, &work[
			    j * work_dim1 + 1], &c__1);
/* L190: */
		}

/*              W := W * V2**T */

		igraphdtrmm_("Right", "Lower", "Transpose", "Unit", &lastc, k, &
			c_b14, &v[(lastv - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork);
		if (lastv > *k) {

/*                 W := W + C1**T * V1**T */

		    i__1 = lastv - *k;
		    igraphdgemm_("Transpose", "Transpose", &lastc, k, &i__1, &c_b14,
			     &c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork);
		}

/*              W := W * T**T  or  W * T */

		igraphdtrmm_("Right", "Lower", transt, "Non-unit", &lastc, k, &
			c_b14, &t[t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - V**T * W**T */

		if (lastv > *k) {

/*                 C1 := C1 - V1**T * W**T */

		    i__1 = lastv - *k;
		    igraphdgemm_("Transpose", "Transpose", &i__1, &lastc, k, &c_b25,
			     &v[v_offset], ldv, &work[work_offset], ldwork, &
			    c_b14, &c__[c_offset], ldc);
		}

/*              W := W * V2 */

		igraphdtrmm_("Right", "Lower", "No transpose", "Unit", &lastc, k, &
			c_b14, &v[(lastv - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork);

/*              C2 := C2 - W**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = lastc;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[lastv - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
/* L200: */
		    }
/* L210: */
		}

	    } else if (igraphlsame_(side, "R")) {

/*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )   

   Computing MAX */
		i__1 = *k, i__2 = igraphiladlc_(k, n, &v[v_offset], ldv);
		lastv = max(i__1,i__2);
		lastc = igraphiladlr_(m, &lastv, &c__[c_offset], ldc);

/*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)   

                W := C2 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    igraphdcopy_(&lastc, &c__[(lastv - *k + j) * c_dim1 + 1], &c__1,
			     &work[j * work_dim1 + 1], &c__1);
/* L220: */
		}

/*              W := W * V2**T */

		igraphdtrmm_("Right", "Lower", "Transpose", "Unit", &lastc, k, &
			c_b14, &v[(lastv - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork);
		if (lastv > *k) {

/*                 W := W + C1 * V1**T */

		    i__1 = lastv - *k;
		    igraphdgemm_("No transpose", "Transpose", &lastc, k, &i__1, &
			    c_b14, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T  or  W * T**T */

		igraphdtrmm_("Right", "Lower", trans, "Non-unit", &lastc, k, &c_b14,
			 &t[t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - W * V */

		if (lastv > *k) {

/*                 C1 := C1 - W * V1 */

		    i__1 = lastv - *k;
		    igraphdgemm_("No transpose", "No transpose", &lastc, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v[v_offset], 
			    ldv, &c_b14, &c__[c_offset], ldc);
		}

/*              W := W * V2 */

		igraphdtrmm_("Right", "Lower", "No transpose", "Unit", &lastc, k, &
			c_b14, &v[(lastv - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = lastc;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + (lastv - *k + j) * c_dim1] -= work[i__ + j *
				 work_dim1];
/* L230: */
		    }
/* L240: */
		}

	    }

	}
    }

    return 0;

/*     End of DLARFB */

} /* igraphdlarfb_ */

