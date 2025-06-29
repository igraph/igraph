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

static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static integer c__1 = 1;

/* > \brief \b DLARF applies an elementary reflector to a general rectangular matrix.   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DLARF + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf.f
">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf.f
">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf.f
">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )   

         CHARACTER          SIDE   
         INTEGER            INCV, LDC, M, N   
         DOUBLE PRECISION   TAU   
         DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DLARF applies a real elementary reflector H to a real m by n matrix   
   > C, from either the left or the right. H is represented in the form   
   >   
   >       H = I - tau * v * v**T   
   >   
   > where tau is a real scalar and v is a real vector.   
   >   
   > If tau = 0, then H is taken to be the unit matrix.   
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
   >          V is DOUBLE PRECISION array, dimension   
   >                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'   
   >                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'   
   >          The vector v in the representation of H. V is not used if   
   >          TAU = 0.   
   > \endverbatim   
   >   
   > \param[in] INCV   
   > \verbatim   
   >          INCV is INTEGER   
   >          The increment between elements of v. INCV <> 0.   
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
   >          The leading dimension of the array C. LDC >= max(1,M).   
   > \endverbatim   
   >   
   > \param[out] WORK   
   > \verbatim   
   >          WORK is DOUBLE PRECISION array, dimension   
   >                         (N) if SIDE = 'L'   
   >                      or (M) if SIDE = 'R'   
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
   Subroutine */ int igraphdlarf_(char *side, integer *m, integer *n, doublereal *v,
	 integer *incv, doublereal *tau, doublereal *c__, integer *ldc, 
	doublereal *work)
{
    /* System generated locals */
    integer c_dim1, c_offset;
    doublereal d__1;

    /* Local variables */
    integer i__;
    logical applyleft;
    extern /* Subroutine */ int igraphdger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern logical igraphlsame_(char *, char *);
    extern /* Subroutine */ int igraphdgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    integer lastc, lastv;
    extern integer igraphiladlc_(integer *, integer *, doublereal *, integer *), 
	    igraphiladlr_(integer *, integer *, doublereal *, integer *);


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
    applyleft = igraphlsame_(side, "L");
    lastv = 0;
    lastc = 0;
    if (*tau != 0.) {
/*     Set up variables for scanning V.  LASTV begins pointing to the end   
       of V. */
	if (applyleft) {
	    lastv = *m;
	} else {
	    lastv = *n;
	}
	if (*incv > 0) {
	    i__ = (lastv - 1) * *incv + 1;
	} else {
	    i__ = 1;
	}
/*     Look for the last non-zero row in V. */
	while(lastv > 0 && v[i__] == 0.) {
	    --lastv;
	    i__ -= *incv;
	}
	if (applyleft) {
/*     Scan for the last non-zero column in C(1:lastv,:). */
	    lastc = igraphiladlc_(&lastv, n, &c__[c_offset], ldc);
	} else {
/*     Scan for the last non-zero row in C(:,1:lastv). */
	    lastc = igraphiladlr_(m, &lastv, &c__[c_offset], ldc);
	}
    }
/*     Note that lastc.eq.0 renders the BLAS operations null; no special   
       case is needed at this level. */
    if (applyleft) {

/*        Form  H * C */

	if (lastv > 0) {

/*           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1) */

	    igraphdgemv_("Transpose", &lastv, &lastc, &c_b4, &c__[c_offset], ldc, &
		    v[1], incv, &c_b5, &work[1], &c__1);

/*           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T */

	    d__1 = -(*tau);
	    igraphdger_(&lastv, &lastc, &d__1, &v[1], incv, &work[1], &c__1, &c__[
		    c_offset], ldc);
	}
    } else {

/*        Form  C * H */

	if (lastv > 0) {

/*           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1) */

	    igraphdgemv_("No transpose", &lastc, &lastv, &c_b4, &c__[c_offset], ldc,
		     &v[1], incv, &c_b5, &work[1], &c__1);

/*           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T */

	    d__1 = -(*tau);
	    igraphdger_(&lastc, &lastv, &d__1, &work[1], &c__1, &v[1], incv, &c__[
		    c_offset], ldc);
	}
    }
    return 0;

/*     End of DLARF */

} /* igraphdlarf_ */

