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

/* > \brief \b DLARFX applies an elementary reflector to a general rectangular matrix, with loop unrolling whe
n the reflector has order ≤ 10.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLARFX + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfx.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfx.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfx.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )   

         CHARACTER          SIDE   
         INTEGER            LDC, M, N   
         DOUBLE PRECISION   TAU   
         DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DLARFX applies a real elementary reflector H to a real m by n   
   > matrix C, from either the left or the right. H is represented in the   
   > form   
   >   
   >       H = I - tau * v * v**T   
   >   
   > where tau is a real scalar and v is a real vector.   
   >   
   > If tau = 0, then H is taken to be the unit matrix   
   >   
   > This version uses inline code if H has order < 11.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] SIDE   
   > \verbatim   
   >          SIDE is CHARACTER*1   
   >          = 'L': form  H * C   
   >          = 'R': form  C * H   
   > \endverbatim   
   >   
   > \param[in] M   
   > \verbatim   
   >          M is INTEGER   
   >          The number of rows of the matrix C.   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The number of columns of the matrix C.   
   > \endverbatim   
   >   
   > \param[in] V   
   > \verbatim   
   >          V is DOUBLE PRECISION array, dimension (M) if SIDE = 'L'   
   >                                     or (N) if SIDE = 'R'   
   >          The vector v in the representation of H.   
   > \endverbatim   
   >   
   > \param[in] TAU   
   > \verbatim   
   >          TAU is DOUBLE PRECISION   
   >          The value tau in the representation of H.   
   > \endverbatim   
   >   
   > \param[in,out] C   
   > \verbatim   
   >          C is DOUBLE PRECISION array, dimension (LDC,N)   
   >          On entry, the m by n matrix C.   
   >          On exit, C is overwritten by the matrix H * C if SIDE = 'L',   
   >          or C * H if SIDE = 'R'.   
   > \endverbatim   
   >   
   > \param[in] LDC   
   > \verbatim   
   >          LDC is INTEGER   
   >          The leading dimension of the array C. LDA >= (1,M).   
   > \endverbatim   
   >   
   > \param[out] WORK   
   > \verbatim   
   >          WORK is DOUBLE PRECISION array, dimension   
   >                      (N) if SIDE = 'L'   
   >                      or (M) if SIDE = 'R'   
   >          WORK is not referenced if H has order < 11.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date September 2012   

   > \ingroup doubleOTHERauxiliary   

    =====================================================================   
   Subroutine */ int igraphdlarfx_(char *side, integer *m, integer *n, doublereal *
	v, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1;

    /* Local variables */
    integer j;
    doublereal t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, v6, v7,
	     v8, v9, t10, v10, sum;
    extern /* Subroutine */ int igraphdlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *);
    extern logical igraphlsame_(char *, char *);


/*  -- LAPACK auxiliary routine (version 3.4.2) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       September 2012   


    =====================================================================   


       Parameter adjustments */
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    if (*tau == 0.) {
	return 0;
    }
    if (igraphlsame_(side, "L")) {

/*        Form  H * C, where H has order m. */

	switch (*m) {
	    case 1:  goto L10;
	    case 2:  goto L30;
	    case 3:  goto L50;
	    case 4:  goto L70;
	    case 5:  goto L90;
	    case 6:  goto L110;
	    case 7:  goto L130;
	    case 8:  goto L150;
	    case 9:  goto L170;
	    case 10:  goto L190;
	}

/*        Code for general M */

	igraphdlarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1]);
	goto L410;
L10:

/*        Special code for 1 x 1 Householder */

	t1 = 1. - *tau * v[1] * v[1];
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    c__[j * c_dim1 + 1] = t1 * c__[j * c_dim1 + 1];
/* L20: */
	}
	goto L410;
L30:

/*        Special code for 2 x 2 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2];
	    c__[j * c_dim1 + 1] -= sum * t1;
	    c__[j * c_dim1 + 2] -= sum * t2;
/* L40: */
	}
	goto L410;
L50:

/*        Special code for 3 x 3 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3];
	    c__[j * c_dim1 + 1] -= sum * t1;
	    c__[j * c_dim1 + 2] -= sum * t2;
	    c__[j * c_dim1 + 3] -= sum * t3;
/* L60: */
	}
	goto L410;
L70:

/*        Special code for 4 x 4 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4];
	    c__[j * c_dim1 + 1] -= sum * t1;
	    c__[j * c_dim1 + 2] -= sum * t2;
	    c__[j * c_dim1 + 3] -= sum * t3;
	    c__[j * c_dim1 + 4] -= sum * t4;
/* L80: */
	}
	goto L410;
L90:

/*        Special code for 5 x 5 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5];
	    c__[j * c_dim1 + 1] -= sum * t1;
	    c__[j * c_dim1 + 2] -= sum * t2;
	    c__[j * c_dim1 + 3] -= sum * t3;
	    c__[j * c_dim1 + 4] -= sum * t4;
	    c__[j * c_dim1 + 5] -= sum * t5;
/* L100: */
	}
	goto L410;
L110:

/*        Special code for 6 x 6 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6];
	    c__[j * c_dim1 + 1] -= sum * t1;
	    c__[j * c_dim1 + 2] -= sum * t2;
	    c__[j * c_dim1 + 3] -= sum * t3;
	    c__[j * c_dim1 + 4] -= sum * t4;
	    c__[j * c_dim1 + 5] -= sum * t5;
	    c__[j * c_dim1 + 6] -= sum * t6;
/* L120: */
	}
	goto L410;
L130:

/*        Special code for 7 x 7 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j * 
		    c_dim1 + 7];
	    c__[j * c_dim1 + 1] -= sum * t1;
	    c__[j * c_dim1 + 2] -= sum * t2;
	    c__[j * c_dim1 + 3] -= sum * t3;
	    c__[j * c_dim1 + 4] -= sum * t4;
	    c__[j * c_dim1 + 5] -= sum * t5;
	    c__[j * c_dim1 + 6] -= sum * t6;
	    c__[j * c_dim1 + 7] -= sum * t7;
/* L140: */
	}
	goto L410;
L150:

/*        Special code for 8 x 8 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	v8 = v[8];
	t8 = *tau * v8;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j * 
		    c_dim1 + 7] + v8 * c__[j * c_dim1 + 8];
	    c__[j * c_dim1 + 1] -= sum * t1;
	    c__[j * c_dim1 + 2] -= sum * t2;
	    c__[j * c_dim1 + 3] -= sum * t3;
	    c__[j * c_dim1 + 4] -= sum * t4;
	    c__[j * c_dim1 + 5] -= sum * t5;
	    c__[j * c_dim1 + 6] -= sum * t6;
	    c__[j * c_dim1 + 7] -= sum * t7;
	    c__[j * c_dim1 + 8] -= sum * t8;
/* L160: */
	}
	goto L410;
L170:

/*        Special code for 9 x 9 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	v8 = v[8];
	t8 = *tau * v8;
	v9 = v[9];
	t9 = *tau * v9;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j * 
		    c_dim1 + 7] + v8 * c__[j * c_dim1 + 8] + v9 * c__[j * 
		    c_dim1 + 9];
	    c__[j * c_dim1 + 1] -= sum * t1;
	    c__[j * c_dim1 + 2] -= sum * t2;
	    c__[j * c_dim1 + 3] -= sum * t3;
	    c__[j * c_dim1 + 4] -= sum * t4;
	    c__[j * c_dim1 + 5] -= sum * t5;
	    c__[j * c_dim1 + 6] -= sum * t6;
	    c__[j * c_dim1 + 7] -= sum * t7;
	    c__[j * c_dim1 + 8] -= sum * t8;
	    c__[j * c_dim1 + 9] -= sum * t9;
/* L180: */
	}
	goto L410;
L190:

/*        Special code for 10 x 10 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	v8 = v[8];
	t8 = *tau * v8;
	v9 = v[9];
	t9 = *tau * v9;
	v10 = v[10];
	t10 = *tau * v10;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j * 
		    c_dim1 + 7] + v8 * c__[j * c_dim1 + 8] + v9 * c__[j * 
		    c_dim1 + 9] + v10 * c__[j * c_dim1 + 10];
	    c__[j * c_dim1 + 1] -= sum * t1;
	    c__[j * c_dim1 + 2] -= sum * t2;
	    c__[j * c_dim1 + 3] -= sum * t3;
	    c__[j * c_dim1 + 4] -= sum * t4;
	    c__[j * c_dim1 + 5] -= sum * t5;
	    c__[j * c_dim1 + 6] -= sum * t6;
	    c__[j * c_dim1 + 7] -= sum * t7;
	    c__[j * c_dim1 + 8] -= sum * t8;
	    c__[j * c_dim1 + 9] -= sum * t9;
	    c__[j * c_dim1 + 10] -= sum * t10;
/* L200: */
	}
	goto L410;
    } else {

/*        Form  C * H, where H has order n. */

	switch (*n) {
	    case 1:  goto L210;
	    case 2:  goto L230;
	    case 3:  goto L250;
	    case 4:  goto L270;
	    case 5:  goto L290;
	    case 6:  goto L310;
	    case 7:  goto L330;
	    case 8:  goto L350;
	    case 9:  goto L370;
	    case 10:  goto L390;
	}

/*        Code for general N */

	igraphdlarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1]);
	goto L410;
L210:

/*        Special code for 1 x 1 Householder */

	t1 = 1. - *tau * v[1] * v[1];
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    c__[j + c_dim1] = t1 * c__[j + c_dim1];
/* L220: */
	}
	goto L410;
L230:

/*        Special code for 2 x 2 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)];
	    c__[j + c_dim1] -= sum * t1;
	    c__[j + (c_dim1 << 1)] -= sum * t2;
/* L240: */
	}
	goto L410;
L250:

/*        Special code for 3 x 3 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3];
	    c__[j + c_dim1] -= sum * t1;
	    c__[j + (c_dim1 << 1)] -= sum * t2;
	    c__[j + c_dim1 * 3] -= sum * t3;
/* L260: */
	}
	goto L410;
L270:

/*        Special code for 4 x 4 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)];
	    c__[j + c_dim1] -= sum * t1;
	    c__[j + (c_dim1 << 1)] -= sum * t2;
	    c__[j + c_dim1 * 3] -= sum * t3;
	    c__[j + (c_dim1 << 2)] -= sum * t4;
/* L280: */
	}
	goto L410;
L290:

/*        Special code for 5 x 5 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5];
	    c__[j + c_dim1] -= sum * t1;
	    c__[j + (c_dim1 << 1)] -= sum * t2;
	    c__[j + c_dim1 * 3] -= sum * t3;
	    c__[j + (c_dim1 << 2)] -= sum * t4;
	    c__[j + c_dim1 * 5] -= sum * t5;
/* L300: */
	}
	goto L410;
L310:

/*        Special code for 6 x 6 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6];
	    c__[j + c_dim1] -= sum * t1;
	    c__[j + (c_dim1 << 1)] -= sum * t2;
	    c__[j + c_dim1 * 3] -= sum * t3;
	    c__[j + (c_dim1 << 2)] -= sum * t4;
	    c__[j + c_dim1 * 5] -= sum * t5;
	    c__[j + c_dim1 * 6] -= sum * t6;
/* L320: */
	}
	goto L410;
L330:

/*        Special code for 7 x 7 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[
		    j + c_dim1 * 7];
	    c__[j + c_dim1] -= sum * t1;
	    c__[j + (c_dim1 << 1)] -= sum * t2;
	    c__[j + c_dim1 * 3] -= sum * t3;
	    c__[j + (c_dim1 << 2)] -= sum * t4;
	    c__[j + c_dim1 * 5] -= sum * t5;
	    c__[j + c_dim1 * 6] -= sum * t6;
	    c__[j + c_dim1 * 7] -= sum * t7;
/* L340: */
	}
	goto L410;
L350:

/*        Special code for 8 x 8 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	v8 = v[8];
	t8 = *tau * v8;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[
		    j + c_dim1 * 7] + v8 * c__[j + (c_dim1 << 3)];
	    c__[j + c_dim1] -= sum * t1;
	    c__[j + (c_dim1 << 1)] -= sum * t2;
	    c__[j + c_dim1 * 3] -= sum * t3;
	    c__[j + (c_dim1 << 2)] -= sum * t4;
	    c__[j + c_dim1 * 5] -= sum * t5;
	    c__[j + c_dim1 * 6] -= sum * t6;
	    c__[j + c_dim1 * 7] -= sum * t7;
	    c__[j + (c_dim1 << 3)] -= sum * t8;
/* L360: */
	}
	goto L410;
L370:

/*        Special code for 9 x 9 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	v8 = v[8];
	t8 = *tau * v8;
	v9 = v[9];
	t9 = *tau * v9;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[
		    j + c_dim1 * 7] + v8 * c__[j + (c_dim1 << 3)] + v9 * c__[
		    j + c_dim1 * 9];
	    c__[j + c_dim1] -= sum * t1;
	    c__[j + (c_dim1 << 1)] -= sum * t2;
	    c__[j + c_dim1 * 3] -= sum * t3;
	    c__[j + (c_dim1 << 2)] -= sum * t4;
	    c__[j + c_dim1 * 5] -= sum * t5;
	    c__[j + c_dim1 * 6] -= sum * t6;
	    c__[j + c_dim1 * 7] -= sum * t7;
	    c__[j + (c_dim1 << 3)] -= sum * t8;
	    c__[j + c_dim1 * 9] -= sum * t9;
/* L380: */
	}
	goto L410;
L390:

/*        Special code for 10 x 10 Householder */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	v8 = v[8];
	t8 = *tau * v8;
	v9 = v[9];
	t9 = *tau * v9;
	v10 = v[10];
	t10 = *tau * v10;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[
		    j + c_dim1 * 7] + v8 * c__[j + (c_dim1 << 3)] + v9 * c__[
		    j + c_dim1 * 9] + v10 * c__[j + c_dim1 * 10];
	    c__[j + c_dim1] -= sum * t1;
	    c__[j + (c_dim1 << 1)] -= sum * t2;
	    c__[j + c_dim1 * 3] -= sum * t3;
	    c__[j + (c_dim1 << 2)] -= sum * t4;
	    c__[j + c_dim1 * 5] -= sum * t5;
	    c__[j + c_dim1 * 6] -= sum * t6;
	    c__[j + c_dim1 * 7] -= sum * t7;
	    c__[j + (c_dim1 << 3)] -= sum * t8;
	    c__[j + c_dim1 * 9] -= sum * t9;
	    c__[j + c_dim1 * 10] -= sum * t10;
/* L400: */
	}
	goto L410;
    }
L410:
    return 0;

/*     End of DLARFX */

} /* igraphdlarfx_ */

