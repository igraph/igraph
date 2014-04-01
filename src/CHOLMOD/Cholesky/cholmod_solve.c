/* ========================================================================== */
/* === Cholesky/cholmod_solve =============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Cholesky Module.  Copyright (C) 2005-2013, Timothy A. Davis
 * The CHOLMOD/Cholesky Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Solve one of the following systems.  D is identity for an LL' factorization,
 * in which the D operation is skipped:
 *
 *      Ax=b        0: CHOLMOD_A     x = P' * (L' \ (D \ (L \ (P * b))))
 *      LDL'x=b     1: CHOLMOD_LDLt  x =      (L' \ (D \ (L \ (    b))))
 *      LDx=b       2: CHOLMOD_LD    x =      (     (D \ (L \ (    b))))
 *      DL'x=b      3: CHOLMOD_DLt   x =      (L' \ (D \ (    (    b))))
 *      Lx=b        4: CHOLMOD_L     x =      (     (    (L \ (    b))))
 *      L'x=b       5: CHOLMOD_Lt    x =      (L' \ (    (    (    b))))
 *      Dx=b        6: CHOLMOD_D     x =      (     (D \ (    (    b))))
 *      x=Pb        7: CHOLMOD_P     x =      (     (    (    (P * b))))
 *      x=P'b       8: CHOLMOD_Pt    x = P' * (     (    (    (    b))))
 *
 * The factorization can be simplicial LDL', simplicial LL', or supernodal LL'.
 * For an LL' factorization, D is the identity matrix.  Thus CHOLMOD_LD and
 * CHOLMOD_L solve the same system if an LL' factorization was performed,
 * for example.
 *
 * The supernodal solver uses BLAS routines dtrsv, dgemv, dtrsm, and dgemm,
 * or their complex counterparts ztrsv, zgemv, ztrsm, and zgemm.
 *
 * If both L and B are real, then X is returned real.  If either is complex
 * or zomplex, X is returned as either complex or zomplex, depending on the
 * Common->prefer_zomplex parameter.
 *
 * Supports any numeric xtype (pattern-only matrices not supported).
 *
 * This routine does not check to see if the diagonal of L or D is zero,
 * because sometimes a partial solve can be done with indefinite or singular
 * matrix.  If you wish to check in your own code, test L->minor.  If
 * L->minor == L->n, then the matrix has no zero diagonal entries.
 * If k = L->minor < L->n, then L(k,k) is zero for an LL' factorization, or
 * D(k,k) is zero for an LDL' factorization.
 *
 * This routine returns X as NULL only if it runs out of memory.  If L is
 * indefinite or singular, then X may contain Inf's or NaN's, but it will
 * exist on output.
 */

#ifndef NCHOLESKY

#include "cholmod_internal.h"
#include "cholmod_cholesky.h"

#ifndef NSUPERNODAL
#include "cholmod_supernodal.h"
#endif


/* ========================================================================== */
/* === TEMPLATE ============================================================= */
/* ========================================================================== */

#define REAL
#include "t_cholmod_solve.c"

#define COMPLEX
#include "t_cholmod_solve.c"

#define ZOMPLEX
#include "t_cholmod_solve.c"

/* ========================================================================== */
/* === Permutation macro ==================================================== */
/* ========================================================================== */

/* If Perm is NULL, it is interpretted as the identity permutation */

#define P(k) ((Perm == NULL) ? (k) : Perm [k])


/* ========================================================================== */
/* === perm ================================================================= */
/* ========================================================================== */

/* Y = B (P (1:nrow), k1 : min (k1+ncols,ncol)-1) where B is nrow-by-ncol.
 *
 * Creates a permuted copy of a contiguous set of columns of B.
 * Y is already allocated on input.  Y must be of sufficient size.  Let nk be
 * the number of columns accessed in B.  Y->xtype determines the complexity of
 * the result.
 *
 * If B is real and Y is complex (or zomplex), only the real part of B is
 * copied into Y.  The imaginary part of Y is set to zero.
 *
 * If B is complex (or zomplex) and Y is real, both the real and imaginary and
 * parts of B are returned in Y.  Y is returned as nrow-by-2*nk. The even
 * columns of Y contain the real part of B and the odd columns contain the
 * imaginary part of B.  Y->nzmax must be >= 2*nrow*nk.  Otherise, Y is
 * returned as nrow-by-nk with leading dimension nrow.  Y->nzmax must be >=
 * nrow*nk.
 *
 * The case where the input (B) is real and the output (Y) is zomplex is
 * not used.
 */

static void perm
(
    /* ---- input ---- */
    cholmod_dense *B,	/* input matrix B */
    Int *Perm,		/* optional input permutation (can be NULL) */
    Int k1,		/* first column of B to copy */
    Int ncols,		/* last column to copy is min(k1+ncols,B->ncol)-1 */
    /* ---- in/out --- */
    cholmod_dense *Y	/* output matrix Y, already allocated */
)
{
    double *Yx, *Yz, *Bx, *Bz ;
    Int k2, nk, p, k, j, nrow, ncol, d, dual, dj, j2 ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    ncol = B->ncol ;
    nrow = B->nrow ;
    k2 = MIN (k1+ncols, ncol) ;
    nk = MAX (k2 - k1, 0) ;
    dual = (Y->xtype == CHOLMOD_REAL && B->xtype != CHOLMOD_REAL) ? 2 : 1 ;
    d = B->d ;
    Bx = B->x ;
    Bz = B->z ;
    Yx = Y->x ;
    Yz = Y->z ;
    Y->nrow = nrow ;
    Y->ncol = dual*nk ;
    Y->d = nrow ;
    ASSERT (((Int) Y->nzmax) >= nrow*nk*dual) ;

    /* ---------------------------------------------------------------------- */
    /* Y = B (P (1:nrow), k1:k2-1) */
    /* ---------------------------------------------------------------------- */

    switch (Y->xtype)
    {

	case CHOLMOD_REAL:

	    switch (B->xtype)
	    {

		case CHOLMOD_REAL:
		    /* Y real, B real */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [k + j2] = Bx [p] ;		/* real */
			}
		    }
		    break ;

		case CHOLMOD_COMPLEX:
		    /* Y real, B complex. Y is nrow-by-2*nk */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * 2 * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [k + j2       ] = Bx [2*p  ] ;	/* real */
			    Yx [k + j2 + nrow] = Bx [2*p+1] ;	/* imag */
			}
		    }
		    break ;

		case CHOLMOD_ZOMPLEX:
		    /* Y real, B zomplex. Y is nrow-by-2*nk */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * 2 * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [k + j2       ] = Bx [p] ;	/* real */
			    Yx [k + j2 + nrow] = Bz [p] ;	/* imag */
			}
		    }
		    break ;

	    }
	    break ;

	case CHOLMOD_COMPLEX:

	    switch (B->xtype)
	    {

		case CHOLMOD_REAL:
		    /* Y complex, B real */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * 2 * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [2*k   + j2] = Bx [p] ;		/* real */
			    Yx [2*k+1 + j2] = 0 ;		/* imag */
			}
		    }
		    break ;

		case CHOLMOD_COMPLEX:
		    /* Y complex, B complex */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * 2 * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [2*k   + j2] = Bx [2*p  ] ;	/* real */
			    Yx [2*k+1 + j2] = Bx [2*p+1] ;	/* imag */
			}
		    }
		    break ;

		case CHOLMOD_ZOMPLEX:
		    /* Y complex, B zomplex */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * 2 * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [2*k   + j2] = Bx [p] ;		/* real */
			    Yx [2*k+1 + j2] = Bz [p] ;		/* imag */
			}
		    }
		    break ;

	    }
	    break ;

	case CHOLMOD_ZOMPLEX:

	    switch (B->xtype)
	    {

#if 0
		case CHOLMOD_REAL:
		    /* this case is not used */
		    break ;
#endif

		case CHOLMOD_COMPLEX:
		    /* Y zomplex, B complex */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [k + j2] = Bx [2*p  ] ;		/* real */
			    Yz [k + j2] = Bx [2*p+1] ;		/* imag */
			}
		    }
		    break ;

		case CHOLMOD_ZOMPLEX:
		    /* Y zomplex, B zomplex */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [k + j2] = Bx [p] ;		/* real */
			    Yz [k + j2] = Bz [p] ;		/* imag */
			}
		    }
		    break ;

	    }
	    break ;

    }
}


/* ========================================================================== */
/* === iperm ================================================================ */
/* ========================================================================== */

/* X (P (1:nrow), k1 : min (k1+ncols,ncol)-1) = Y where X is nrow-by-ncol.
 *
 * Copies and permutes Y into a contiguous set of columns of X.  X is already
 * allocated on input.  Y must be of sufficient size.  Let nk be the number
 * of columns accessed in X.  X->xtype determines the complexity of the result.
 *
 * If X is real and Y is complex (or zomplex), only the real part of B is
 * copied into X.  The imaginary part of Y is ignored.
 *
 * If X is complex (or zomplex) and Y is real, both the real and imaginary and
 * parts of Y are returned in X.  Y is nrow-by-2*nk. The even
 * columns of Y contain the real part of B and the odd columns contain the
 * imaginary part of B.  Y->nzmax must be >= 2*nrow*nk.  Otherise, Y is
 * nrow-by-nk with leading dimension nrow.  Y->nzmax must be >= nrow*nk.
 *
 * The case where the input (Y) is complex and the output (X) is real,
 * and the case where the input (Y) is zomplex and the output (X) is real,
 * are not used.
 */

static void iperm
(
    /* ---- input ---- */
    cholmod_dense *Y,	/* input matrix Y */
    Int *Perm,		/* optional input permutation (can be NULL) */
    Int k1,		/* first column of B to copy */
    Int ncols,		/* last column to copy is min(k1+ncols,B->ncol)-1 */
    /* ---- in/out --- */
    cholmod_dense *X	/* output matrix X, already allocated */
)
{
    double *Yx, *Yz, *Xx, *Xz ;
    Int k2, nk, p, k, j, nrow, ncol, d, dj, j2 ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    ncol = X->ncol ;
    nrow = X->nrow ;
    k2 = MIN (k1+ncols, ncol) ;
    nk = MAX (k2 - k1, 0) ;
    d = X->d ;
    Xx = X->x ;
    Xz = X->z ;
    Yx = Y->x ;
    Yz = Y->z ;
    ASSERT (((Int) Y->nzmax) >= nrow*nk*
	    ((X->xtype != CHOLMOD_REAL && Y->xtype == CHOLMOD_REAL) ? 2:1)) ;

    /* ---------------------------------------------------------------------- */
    /* X (P (1:nrow), k1:k2-1) = Y */
    /* ---------------------------------------------------------------------- */

    switch (Y->xtype)
    {

	case CHOLMOD_REAL:

	    switch (X->xtype)
	    {

		case CHOLMOD_REAL:
		    /* Y real, X real */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [p] = Yx [k + j2] ;		/* real */
			}
		    }
		    break ;

		case CHOLMOD_COMPLEX:
		    /* Y real, X complex. Y is nrow-by-2*nk */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * 2 * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [2*p  ] = Yx [k + j2       ] ;	/* real */
			    Xx [2*p+1] = Yx [k + j2 + nrow] ;	/* imag */
			}
		    }
		    break ;

		case CHOLMOD_ZOMPLEX:
		    /* Y real, X zomplex. Y is nrow-by-2*nk */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * 2 * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [p] = Yx [k + j2       ] ;	/* real */
			    Xz [p] = Yx [k + j2 + nrow] ;	/* imag */
			}
		    }
		    break ;

	    }
	    break ;

	case CHOLMOD_COMPLEX:

	    switch (X->xtype)
	    {

#if 0
		case CHOLMOD_REAL:
		    /* this case is not used */
		    break ;
#endif

		case CHOLMOD_COMPLEX:
		    /* Y complex, X complex */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * 2 * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [2*p  ] = Yx [2*k   + j2] ;	/* real */
			    Xx [2*p+1] = Yx [2*k+1 + j2] ;	/* imag */
			}
		    }
		    break ;

		case CHOLMOD_ZOMPLEX:
		    /* Y complex, X zomplex */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * 2 * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [p] = Yx [2*k   + j2] ;		/* real */
			    Xz [p] = Yx [2*k+1 + j2] ;		/* imag */
			}
		    }
		    break ;

	    }
	    break ;

	case CHOLMOD_ZOMPLEX:

	    switch (X->xtype)
	    {

#if 0
		case CHOLMOD_REAL:
		    /* this case is not used */
		    break ;
#endif

		case CHOLMOD_COMPLEX:
		    /* Y zomplex, X complex */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [2*p  ] = Yx [k + j2] ;		/* real */
			    Xx [2*p+1] = Yz [k + j2] ;		/* imag */
			}
		    }
		    break ;

		case CHOLMOD_ZOMPLEX:
		    /* Y zomplex, X zomplex */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = nrow * (j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [p] = Yx [k + j2] ;		/* real */
			    Xz [p] = Yz [k + j2] ;		/* imag */
			}
		    }
		    break ;

	    }
	    break ;

    }
}


/* ========================================================================== */
/* === ptrans =============================================================== */
/* ========================================================================== */

/* Y = B (P (1:nrow), k1 : min (k1+ncols,ncol)-1)' where B is nrow-by-ncol.
 *
 * Creates a permuted and transposed copy of a contiguous set of columns of B.
 * Y is already allocated on input.  Y must be of sufficient size.  Let nk be
 * the number of columns accessed in B.  Y->xtype determines the complexity of
 * the result.
 *
 * If B is real and Y is complex (or zomplex), only the real part of B is
 * copied into Y.  The imaginary part of Y is set to zero.
 *
 * If B is complex (or zomplex) and Y is real, both the real and imaginary and
 * parts of B are returned in Y.  Y is returned as 2*nk-by-nrow. The even
 * rows of Y contain the real part of B and the odd rows contain the
 * imaginary part of B.  Y->nzmax must be >= 2*nrow*nk.  Otherise, Y is
 * returned as nk-by-nrow with leading dimension nk.  Y->nzmax must be >=
 * nrow*nk.
 *
 * The array transpose is performed, not the complex conjugate transpose.
 */

static void ptrans
(
    /* ---- input ---- */
    cholmod_dense *B,	/* input matrix B */
    Int *Perm,		/* optional input permutation (can be NULL) */
    Int k1,		/* first column of B to copy */
    Int ncols,		/* last column to copy is min(k1+ncols,B->ncol)-1 */
    /* ---- in/out --- */
    cholmod_dense *Y	/* output matrix Y, already allocated */
)
{
    double *Yx, *Yz, *Bx, *Bz ;
    Int k2, nk, p, k, j, nrow, ncol, d, dual, dj, j2 ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    ncol = B->ncol ;
    nrow = B->nrow ;
    k2 = MIN (k1+ncols, ncol) ;
    nk = MAX (k2 - k1, 0) ;
    dual = (Y->xtype == CHOLMOD_REAL && B->xtype != CHOLMOD_REAL) ? 2 : 1 ;
    d = B->d ;
    Bx = B->x ;
    Bz = B->z ;
    Yx = Y->x ;
    Yz = Y->z ;
    Y->nrow = dual*nk ;
    Y->ncol = nrow ;
    Y->d = dual*nk ;
    ASSERT (((Int) Y->nzmax) >= nrow*nk*dual) ;

    /* ---------------------------------------------------------------------- */
    /* Y = B (P (1:nrow), k1:k2-1)' */
    /* ---------------------------------------------------------------------- */

    switch (Y->xtype)
    {

	case CHOLMOD_REAL:

	    switch (B->xtype)
	    {

		case CHOLMOD_REAL:
		    /* Y real, B real  */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = j-k1 ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [j2 + k*nk] = Bx [p] ;		/* real */
			}
		    }
		    break ;

		case CHOLMOD_COMPLEX:
		    /* Y real, B complex. Y is 2*nk-by-nrow */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = 2*(j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [j2   + k*2*nk] = Bx [2*p  ] ;	/* real */
			    Yx [j2+1 + k*2*nk] = Bx [2*p+1] ;	/* imag */
			}
		    }
		    break ;

		case CHOLMOD_ZOMPLEX:
		    /* Y real, B zomplex. Y is 2*nk-by-nrow */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = 2*(j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [j2   + k*2*nk] = Bx [p] ;	/* real */
			    Yx [j2+1 + k*2*nk] = Bz [p] ;	/* imag */
			}
		    }
		    break ;

	    }
	    break ;

	case CHOLMOD_COMPLEX:

	    switch (B->xtype)
	    {

		case CHOLMOD_REAL:
		    /* Y complex, B real  */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = 2*(j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [j2   + k*2*nk] = Bx [p] ;	/* real */
			    Yx [j2+1 + k*2*nk] = 0 ;		/* imag */
			}
		    }
		    break ;

		case CHOLMOD_COMPLEX:
		    /* Y complex, B complex  */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = 2*(j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [j2   + k*2*nk] = Bx [2*p  ] ;	/* real */
			    Yx [j2+1 + k*2*nk] = Bx [2*p+1] ;	/* imag */
			}
		    }
		    break ;

		case CHOLMOD_ZOMPLEX:
		    /* Y complex, B zomplex  */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = 2*(j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [j2   + k*2*nk] = Bx [p] ;	/* real */
			    Yx [j2+1 + k*2*nk] = Bz [p] ;	/* imag */
			}
		    }
		    break ;

	    }
	    break ;

	case CHOLMOD_ZOMPLEX:

	    switch (B->xtype)
	    {

		case CHOLMOD_REAL:
		    /* Y zomplex, B real  */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = j-k1 ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [j2 + k*nk] = Bx [p] ;		/* real */
			    Yz [j2 + k*nk] = 0 ;		/* imag */
			}
		    }
		    break ;

		case CHOLMOD_COMPLEX:
		    /* Y zomplex, B complex  */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = j-k1 ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [j2 + k*nk] = Bx [2*p  ] ;	/* real */
			    Yz [j2 + k*nk] = Bx [2*p+1] ;	/* imag */
			}
		    }
		    break ;

		case CHOLMOD_ZOMPLEX:
		    /* Y zomplex, B zomplex */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = j-k1 ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Yx [j2 + k*nk] = Bx [p] ;		/* real */
			    Yz [j2 + k*nk] = Bz [p] ;		/* imag */
			}
		    }
		    break ;

	    }
	    break ;

    }
}


/* ========================================================================== */
/* === iptrans ============================================================== */
/* ========================================================================== */

/* X (P (1:nrow), k1 : min (k1+ncols,ncol)-1) = Y' where X is nrow-by-ncol.
 *
 * Copies into a permuted and transposed contiguous set of columns of X.
 * X is already allocated on input.  Y must be of sufficient size.  Let nk be
 * the number of columns accessed in X.  X->xtype determines the complexity of
 * the result.
 *
 * If X is real and Y is complex (or zomplex), only the real part of Y is
 * copied into X.  The imaginary part of Y is ignored.
 *
 * If X is complex (or zomplex) and Y is real, both the real and imaginary and
 * parts of X are returned in Y.  Y is 2*nk-by-nrow. The even
 * rows of Y contain the real part of X and the odd rows contain the
 * imaginary part of X.  Y->nzmax must be >= 2*nrow*nk.  Otherise, Y is
 * nk-by-nrow with leading dimension nk.  Y->nzmax must be >= nrow*nk.
 *
 * The case where Y is complex or zomplex, and X is real, is not used.
 *
 * The array transpose is performed, not the complex conjugate transpose.
 */

static void iptrans
(
    /* ---- input ---- */
    cholmod_dense *Y,	/* input matrix Y */
    Int *Perm,		/* optional input permutation (can be NULL) */
    Int k1,		/* first column of X to copy into */
    Int ncols,		/* last column to copy is min(k1+ncols,X->ncol)-1 */
    /* ---- in/out --- */
    cholmod_dense *X	/* output matrix X, already allocated */
)
{
    double *Yx, *Yz, *Xx, *Xz ;
    Int k2, nk, p, k, j, nrow, ncol, d, dj, j2 ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    ncol = X->ncol ;
    nrow = X->nrow ;
    k2 = MIN (k1+ncols, ncol) ;
    nk = MAX (k2 - k1, 0) ;
    d = X->d ;
    Xx = X->x ;
    Xz = X->z ;
    Yx = Y->x ;
    Yz = Y->z ;
    ASSERT (((Int) Y->nzmax) >= nrow*nk*
	    ((X->xtype != CHOLMOD_REAL && Y->xtype == CHOLMOD_REAL) ? 2:1)) ;

    /* ---------------------------------------------------------------------- */
    /* X (P (1:nrow), k1:k2-1) = Y' */
    /* ---------------------------------------------------------------------- */

    switch (Y->xtype)
    {

	case CHOLMOD_REAL:

	    switch (X->xtype)
	    {

		case CHOLMOD_REAL:
		    /* Y real, X real  */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = j-k1 ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [p] = Yx [j2 + k*nk] ;		/* real */
			}
		    }
		    break ;

		case CHOLMOD_COMPLEX:
		    /* Y real, X complex. Y is 2*nk-by-nrow */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = 2*(j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [2*p  ] = Yx [j2   + k*2*nk] ;	/* real */
			    Xx [2*p+1] = Yx [j2+1 + k*2*nk] ;	/* imag */
			}
		    }
		    break ;

		case CHOLMOD_ZOMPLEX:
		    /* Y real, X zomplex. Y is 2*nk-by-nrow */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = 2*(j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [p] = Yx [j2   + k*2*nk] ;	/* real */
			    Xz [p] = Yx [j2+1 + k*2*nk] ;	/* imag */
			}
		    }
		    break ;

	    }
	    break ;

	case CHOLMOD_COMPLEX:

	    switch (X->xtype)
	    {

#if 0
		case CHOLMOD_REAL:
		    /* this case is not used */
		    break ;
#endif

		case CHOLMOD_COMPLEX:
		    /* Y complex, X complex  */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = 2*(j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [2*p  ] = Yx [j2   + k*2*nk] ;	/* real */
			    Xx [2*p+1] = Yx [j2+1 + k*2*nk] ;	/* imag */
			}
		    }
		    break ;

		case CHOLMOD_ZOMPLEX:
		    /* Y complex, X zomplex  */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = 2*(j-k1) ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [p] = Yx [j2   + k*2*nk] ;	/* real */
			    Xz [p] = Yx [j2+1 + k*2*nk] ;	/* imag */
			}
		    }
		    break ;

	    }
	    break ;

	case CHOLMOD_ZOMPLEX:

	    switch (X->xtype)
	    {

#if 0
		case CHOLMOD_REAL:
		    /* this case is not used */
		    break ;
#endif

		case CHOLMOD_COMPLEX:
		    /* Y zomplex, X complex  */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = j-k1 ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [2*p  ] = Yx [j2 + k*nk] ;	/* real */
			    Xx [2*p+1] = Yz [j2 + k*nk] ;	/* imag */
			}
		    }
		    break ;

		case CHOLMOD_ZOMPLEX:
		    /* Y zomplex, X zomplex */
		    for (j = k1 ; j < k2 ; j++)
		    {
			dj = d*j ;
			j2 = j-k1 ;
			for (k = 0 ; k < nrow ; k++)
			{
			    p = P(k) + dj ;
			    Xx [p] = Yx [j2 + k*nk] ;		/* real */
			    Xz [p] = Yz [j2 + k*nk] ;		/* imag */
			}
		    }
		    break ;

	    }
	    break ;

    }
}


/* ========================================================================== */
/* === cholmod_solve ======================================================== */
/* ========================================================================== */

/* Solve a linear system.
 *
 * The factorization can be simplicial LDL', simplicial LL', or supernodal LL'.
 * The Dx=b solve returns silently for the LL' factorizations (it is implicitly
 * identity).
 */

cholmod_dense *CHOLMOD(solve)
(
    /* ---- input ---- */
    int sys,		/* system to solve */
    cholmod_factor *L,	/* factorization to use */
    cholmod_dense *B,	/* right-hand-side */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_dense *Y = NULL, *X = NULL ;
    cholmod_dense *E = NULL ;
    int ok ;

    /* do the solve, allocating workspaces as needed  */
    ok = CHOLMOD (solve2) (sys, L, B, NULL, &X, NULL, &Y, &E, Common) ;

    /* free workspaces if allocated, and free result if an error occured */
    CHOLMOD(free_dense) (&Y, Common) ;
    CHOLMOD(free_dense) (&E, Common) ;
    if (!ok)
    {
        CHOLMOD(free_dense) (&X, Common) ;
    }
    return (X) ;
}


/* ========================================================================== */
/* === cholmod_solve2 ======================================================= */
/* ========================================================================== */

/* This function acts just like cholmod_solve, except that the solution X and
 * the internal workspace (Y and E) can be passed in preallocated.  If the
 * solution X or any required workspaces are not allocated on input, or if they
 * are the wrong size or type, then this function frees them and reallocates
 * them as the proper size and type.  Thus, if you have a sequence of solves to
 * do, you can let this function allocate X, Y, and E on the first call.
 * Subsequent calls to cholmod_solve2 can then reuse this space.  You must
 * then free the workspaces Y and E (and X if desired) when you are finished.
 * For example, the first call to cholmod_l_solve2, below, will solve the
 * requested system.  The next 2 calls (with different right-hand-sides but
 * the same value of "sys") will resuse the workspace and solution X from the
 * first call.  Finally, when all solves are done, you must free the workspaces
 * Y and E (otherwise you will have a memory leak), and you should also free X
 * when you are done with it.  Note that on input, X, Y, and E must be either
 * valid cholmod_dense matrices, or initialized to NULL.  You cannot pass in an
 * uninitialized X, Y, or E.
 *
 *      cholmod_dense *X = NULL, *Y = NULL, *E = NULL ;
 *      ...
 *      cholmod_l_solve2 (sys, L, B1, NULL, &X, NULL, &Y, &E, Common) ;
 *      cholmod_l_solve2 (sys, L, B2, NULL, &X, NULL, &Y, &E, Common) ;
 *      cholmod_l_solve2 (sys, L, B3, NULL, &X, NULL, &Y, &E, Common) ;
 *      cholmod_l_free_dense (&X, Common) ;
 *      cholmod_l_free_dense (&Y, Common) ;
 *      cholmod_l_free_dense (&E, Common) ;
 *
 * The equivalent when using cholmod_l_solve is:
 *
 *      cholmod_dense *X = NULL, *Y = NULL, *E = NULL ;
 *      ...
 *      X = cholmod_l_solve (sys, L, B1, Common) ;
 *      cholmod_l_free_dense (&X, Common) ;
 *      X = cholmod_l_solve (sys, L, B2, Common) ;
 *      cholmod_l_free_dense (&X, Common) ;
 *      X = cholmod_l_solve (sys, L, B3, Common) ;
 *      cholmod_l_free_dense (&X, Common) ;
 *
 * Both methods work fine, but in the 2nd method with cholmod_solve, the
 * internal workspaces (Y and E) are allocated and freed on each call.
 *
 * Bset is an optional sparse column (pattern only) that specifies a set
 * of row indices.  It is ignored if NULL, or if sys is CHOLMOD_P or
 * CHOLMOD_Pt.  If it is present and not ignored, B must be a dense column
 * vector, and only entries B(i) where i is in the pattern of Bset are
 * considered.  All others are treated as if they were zero (they are not
 * accessed).  L must be a simplicial factorization, not supernodal.  L is
 * converted from supernodal to simplicial if necessary.  The solution X is
 * defined only for entries in the output sparse pattern of Xset.
 * The xtype (real/complex/zomplex) of L and B must match.
 *
 * NOTE: If Bset is present and L is supernodal, it is converted to simplicial
 * on output.
 */

int CHOLMOD(solve2)         /* returns TRUE on success, FALSE on failure */
(
    /* ---- input ---- */
    int sys,		            /* system to solve */
    cholmod_factor *L,	            /* factorization to use */
    cholmod_dense *B,               /* right-hand-side */
    cholmod_sparse *Bset,
    /* ---- output --- */
    cholmod_dense **X_Handle,       /* solution, allocated if need be */
    cholmod_sparse **Xset_Handle,
    /* ---- workspace  */
    cholmod_dense **Y_Handle,       /* workspace, or NULL */
    cholmod_dense **E_Handle,       /* workspace, or NULL */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Yx, *Yz, *Bx, *Bz, *Xx, *Xz ;
    cholmod_dense *Y = NULL, *X = NULL ;
    cholmod_sparse *C, *Yset, C_header, Yset_header, *Xset ;
    Int *Perm = NULL, *IPerm = NULL ;
    Int n, nrhs, ncols, ctype, xtype, k1, nr, ytype, k, blen, p, i, d, nrow ;
    Int Cp [2], Ysetp [2], *Ci, *Yseti, ysetlen ;
    Int *Bsetp, *Bseti, *Bsetnz, *Xseti, *Xsetp, *Iwork ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_NULL (B, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (B, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    if (sys < CHOLMOD_A || sys > CHOLMOD_Pt)
    {
	ERROR (CHOLMOD_INVALID, "invalid system") ;
	return (FALSE) ;
    }
    DEBUG (CHOLMOD(dump_factor) (L, "L", Common)) ;
    DEBUG (CHOLMOD(dump_dense) (B, "B", Common)) ;
    nrhs = B->ncol ;
    n = (Int) L->n ;
    d = (Int) B->d ;
    nrow = (Int) B->nrow ;
    if (d < n || nrow != n)
    {
	ERROR (CHOLMOD_INVALID, "dimensions of L and B do not match") ;
	return (FALSE) ;
    }
    if (Bset)
    {
        if (nrhs != 1)
        {
            ERROR (CHOLMOD_INVALID, "Bset requires a single right-hand side") ;
            return (FALSE) ;
        }
        if (L->xtype != B->xtype)
        {
            ERROR (CHOLMOD_INVALID, "Bset requires xtype of L and B to match") ;
            return (FALSE) ;
        }
        DEBUG (CHOLMOD(dump_sparse) (Bset, "Bset", Common)) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    if ((sys == CHOLMOD_P || sys == CHOLMOD_Pt || sys == CHOLMOD_A)
	    && L->ordering != CHOLMOD_NATURAL)
    {
        /* otherwise, Perm is NULL, and the identity permutation is used */
	Perm = L->Perm ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate the result X (or resuse the space from a prior call) */
    /* ---------------------------------------------------------------------- */

    ctype = (Common->prefer_zomplex) ? CHOLMOD_ZOMPLEX : CHOLMOD_COMPLEX ;

    if (Bset)
    {
        xtype = L->xtype ;
    }
    else if (sys == CHOLMOD_P || sys == CHOLMOD_Pt)
    {
	/* x=Pb and x=P'b return X real if B is real; X is the preferred
	 * complex/zcomplex type if B is complex or zomplex */
	xtype = (B->xtype == CHOLMOD_REAL) ? CHOLMOD_REAL : ctype ;
    }
    else if (L->xtype == CHOLMOD_REAL && B->xtype == CHOLMOD_REAL)
    {
	/* X is real if both L and B are real */
	xtype = CHOLMOD_REAL ;
    }
    else
    {
	/* X is complex, use the preferred complex/zomplex type */
	xtype = ctype ;
    }

    /* ensure X has the right size and type */
    X = CHOLMOD(ensure_dense) (X_Handle, n, nrhs, n, xtype, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* solve using L, D, L', P, or some combination */
    /* ---------------------------------------------------------------------- */

    if (Bset)
    {

        /* ------------------------------------------------------------------ */
        /* solve for a subset of x, with a sparse b */
        /* ------------------------------------------------------------------ */

        Int save_realloc_state ;

#ifndef NSUPERNODAL
        /* convert a supernodal L to simplicial when using Bset */
        if (L->is_super)
        {
            /* Can only use Bset on a simplicial factorization.  The supernodal
             * factor L is converted to simplicial, leaving the xtype unchanged
             * (real, complex, or zomplex).  Since the supernodal factorization
             * is already LL', it is left in that form.   This conversion uses
             * the ll_super_to_simplicial_numeric function in
             * cholmod_change_factor.
             */
            CHOLMOD(change_factor) (
                CHOLMOD_REAL,   /* ignored, since L is already numeric */
                TRUE,           /* convert to LL' (no change to num. values) */
                FALSE,          /* convert to simplicial */
                FALSE,          /* do not pack the columns of L */
                FALSE,          /* (ignored) */
                L, Common) ;
            if (Common->status < CHOLMOD_OK)
            {
                /* out of memory, L is returned unchanged */
                return (FALSE) ;
            }
        }
#endif

        /* L, X, and B are all the same xtype */
        /* ensure Y is the the right size */
	Y = CHOLMOD(ensure_dense) (Y_Handle, 1, n, 1, L->xtype, Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
	    return (FALSE) ;
	}

        /* ------------------------------------------------------------------ */
        /* get the inverse permutation, constructing it if needed */
        /* ------------------------------------------------------------------ */

        DEBUG (CHOLMOD (dump_perm) (Perm,  n,n, "Perm",  Common)) ;

        if ((sys == CHOLMOD_A || sys == CHOLMOD_P) && Perm != NULL)
        {
            /* The inverse permutation IPerm is used for the c=Pb step,
               which is needed only for solving Ax=b or x=Pb.  No other
               steps should use IPerm */
            if (L->IPerm == NULL)
            {
                /* construct the inverse permutation.  This is done only once
                 * and then stored in L permanently.  */
                L->IPerm = CHOLMOD(malloc) (n, sizeof (Int), Common) ;
                if (Common->status < CHOLMOD_OK)
                {
                    /* out of memory */
                    return (FALSE) ;
                }
                IPerm = L->IPerm ;
                for (k = 0 ; k < n ; k++)
                {
                    IPerm [Perm [k]] = k ;
                }
            }
            /* x=A\b and x=Pb both need IPerm */
            IPerm = L->IPerm ;
        }

        if (sys == CHOLMOD_P)
        {
            /* x=Pb needs to turn off the subsequent x=P'b permutation */
            Perm = NULL ;
        }

        DEBUG (CHOLMOD (dump_perm) (Perm,  n,n, "Perm",  Common)) ;
        DEBUG (CHOLMOD (dump_perm) (IPerm, n,n, "IPerm", Common)) ;

        /* ------------------------------------------------------------------ */
        /* ensure Xset is the right size and type */
        /* ------------------------------------------------------------------ */

        /* Xset is n-by-1, nzmax >= n, pattern-only, packed, unsorted */
        Xset = *Xset_Handle ;
        if (Xset == NULL || (Int) Xset->nrow != n || (Int) Xset->ncol != 1 ||
            (Int) Xset->nzmax < n || Xset->itype != CHOLMOD_PATTERN)
        {
            /* this is done only once, for the 1st call to cholmod_solve */
            CHOLMOD(free_sparse) (Xset_Handle, Common) ;
            Xset = CHOLMOD(allocate_sparse) (n, 1, n, FALSE, TRUE, 0,
                CHOLMOD_PATTERN, Common) ;
            *Xset_Handle = Xset ;
        }
        Xset->sorted = FALSE ;
        Xset->stype = 0 ;
        if (Common->status < CHOLMOD_OK)
        {
            /* out of memory */
            return (FALSE) ;
        }

        /* -------------------------------------------------------------- */
        /* ensure Flag of size n, and 3*n Int workspace is available */
        /* -------------------------------------------------------------- */

        /* does no work if prior calls already allocated enough space */
        CHOLMOD(allocate_work) (n, 3*n, 0, Common) ;
        if (Common->status < CHOLMOD_OK)
        {
            /* out of memory */
            return (FALSE) ;
        }

        /* [ use Iwork (n:3n-1) for Ci and Yseti */
        Iwork = Common->Iwork ;
        /* Iwork (0:n-1) is not used because it is used by check_perm,
           print_perm, check_sparse, and print_sparse */
        Ci = Iwork + n ;
        Yseti = Ci + n ;

        /* reallocating workspace would break Ci and Yseti */
        save_realloc_state = Common->no_workspace_reallocate ;
        Common->no_workspace_reallocate = TRUE ;

        /* -------------------------------------------------------------- */
        /* C = permuted Bset, to correspond to the permutation of L */
        /* -------------------------------------------------------------- */

        /* C = IPerm (Bset) */
        DEBUG (CHOLMOD(dump_sparse) (Bset, "Bset", Common)) ;

        Bsetp = Bset->p ;
        Bseti = Bset->i ;
        Bsetnz = Bset->nz ;
        blen = (Bset->packed) ? Bsetp [1] : Bsetnz [0] ;

        /* C = spones (P*B) or C = spones (B) if IPerm is NULL */
        C = &C_header ;
        C->nrow = n ;
        C->ncol = 1 ;
        C->nzmax = n ;
        C->packed = TRUE ;
        C->stype = 0 ;
        C->itype = ITYPE ;
        C->xtype = CHOLMOD_PATTERN ;
        C->dtype = CHOLMOD_DOUBLE ;
        C->nz = NULL ;
        C->p = Cp ;
        C->i = Ci ;
        C->x = NULL ;
        C->z = NULL ;
        C->sorted = FALSE ;
        Cp [0] = 0 ;
        Cp [1] = blen ;
        for (p = 0 ; p < blen ; p++)
        {
            Int iold = Bseti [p] ;
            Ci [p] = IPerm ? IPerm [iold] : iold ;
        }
        DEBUG (CHOLMOD (dump_sparse) (C, "C", Common)) ;

        /* create a sparse column Yset from Iwork (n:2n-1) */
        Yset = &Yset_header ;
        Yset->nrow = n ;
        Yset->ncol = 1 ;
        Yset->nzmax = n ;
        Yset->packed = TRUE ;
        Yset->stype = 0 ;
        Yset->itype = ITYPE ;
        Yset->xtype = CHOLMOD_PATTERN ;
        Yset->dtype = CHOLMOD_DOUBLE ;
        Yset->nz = NULL ;
        Yset->p = Ysetp ;
        Yset->i = Yseti ;
        Yset->x = NULL ;
        Yset->z = NULL ;
        Yset->sorted = FALSE ;
        Ysetp [0] = 0 ;
        Ysetp [1] = 0 ;
        DEBUG (CHOLMOD (dump_sparse) (Yset, "Yset empty", Common)) ;

        /* -------------------------------------------------------------- */
        /* Yset = nonzero pattern of L\C, or just C itself */
        /* -------------------------------------------------------------- */

        /* this takes O(ysetlen) time  */
        if (sys == CHOLMOD_P || sys == CHOLMOD_Pt || sys == CHOLMOD_D)
        {
            Ysetp [1] = blen ;
            for (p = 0 ; p < blen ; p++)
            {
                Yseti [p] = Ci [p] ;
            }
        }
        else
        {
            if (!CHOLMOD(lsolve_pattern) (C, L, Yset, Common))
            {
                Common->no_workspace_reallocate = save_realloc_state ;
                return (FALSE) ;
            }
        }
        DEBUG (CHOLMOD (dump_sparse) (Yset, "Yset", Common)) ;

        /* -------------------------------------------------------------- */
        /* clear the parts of Y that we will use in the solve */
        /* -------------------------------------------------------------- */

        Yx = Y->x ;
        Yz = Y->z ;
        ysetlen = Ysetp [1] ;

        switch (L->xtype)
        {

            case CHOLMOD_REAL:
                for (p = 0 ; p < ysetlen ; p++)
                {
                    i = Yseti [p] ;
                    Yx [i] = 0 ;
                }
                break ;

            case CHOLMOD_COMPLEX:
                for (p = 0 ; p < ysetlen ; p++)
                {
                    i = Yseti [p] ;
                    Yx [2*i  ] = 0 ;
                    Yx [2*i+1] = 0 ;
                }
                break ;

            case CHOLMOD_ZOMPLEX:
                for (p = 0 ; p < ysetlen ; p++)
                {
                    i = Yseti [p] ;
                    Yx [i] = 0 ;
                    Yz [i] = 0 ;
                }
                break ;
        }

        DEBUG (CHOLMOD (dump_dense) (Y, "Y (Yset) = 0", Common)) ;

        /* -------------------------------------------------------------- */
        /* scatter and permute B into Y */
        /* -------------------------------------------------------------- */

        /* Y (C) = B (Bset) */
        Bx = B->x ;
        Bz = B->z ;

        switch (L->xtype)
        {

            case CHOLMOD_REAL:
                for (p = 0 ; p < blen ; p++)
                {
                    Int iold = Bseti [p] ;
                    Int inew = Ci [p] ;
                    Yx [inew] = Bx [iold] ;
                }
                break ;

            case CHOLMOD_COMPLEX:
                for (p = 0 ; p < blen ; p++)
                {
                    Int iold = Bseti [p] ;
                    Int inew = Ci [p] ;
                    Yx [2*inew  ] = Bx [2*iold  ] ;
                    Yx [2*inew+1] = Bx [2*iold+1] ;
                }
                break ;

            case CHOLMOD_ZOMPLEX:
                for (p = 0 ; p < blen ; p++)
                {
                    Int iold = Bseti [p] ;
                    Int inew = Ci [p] ;
                    Yx [inew] = Bx [iold] ;
                    Yz [inew] = Bz [iold] ;
                }
                break ;
        }

        DEBUG (CHOLMOD (dump_dense) (Y, "Y (C) = B (Bset)", Common)) ;

        /* -------------------------------------------------------------- */
        /* solve Y = (L' \ (L \ Y'))', or other system, with template */
        /* -------------------------------------------------------------- */

        /* the solve only iterates over columns in Yseti [0...ysetlen-1] */

        if (! (sys == CHOLMOD_P || sys == CHOLMOD_Pt))
        {
            switch (L->xtype)
            {
                case CHOLMOD_REAL:
                    r_simplicial_solver (sys, L, Y, Yseti, ysetlen) ;
                    break ;

                case CHOLMOD_COMPLEX:
                    c_simplicial_solver (sys, L, Y, Yseti, ysetlen) ;
                    break ;

                case CHOLMOD_ZOMPLEX:
                    z_simplicial_solver (sys, L, Y, Yseti, ysetlen) ;
                    break ;
            }
        }

        DEBUG (CHOLMOD (dump_dense) (Y, "Y after solve", Common)) ;

        /* -------------------------------------------------------------- */
        /* X = P'*Y, but only for rows in Yset, and create Xset */
        /* -------------------------------------------------------------- */

        /* X (Perm (Yset)) = Y (Yset) */
        Xx = X->x ;
        Xz = X->z ;
        Xseti = Xset->i ;
        Xsetp = Xset->p ;

        switch (L->xtype)
        {

            case CHOLMOD_REAL:
                for (p = 0 ; p < ysetlen ; p++)
                {
                    Int inew = Yseti [p] ;
                    Int iold = Perm ? Perm [inew] : inew ;
                    Xx [iold] = Yx [inew] ;
                    Xseti [p] = iold ;
                }
                break ;

            case CHOLMOD_COMPLEX:
                for (p = 0 ; p < ysetlen ; p++)
                {
                    Int inew = Yseti [p] ;
                    Int iold = Perm ? Perm [inew] : inew ;
                    Xx [2*iold  ] = Yx [2*inew] ;
                    Xx [2*iold+1] = Yx [2*inew+1] ;
                    Xseti [p] = iold ;
                }
                break ;

            case CHOLMOD_ZOMPLEX:
                for (p = 0 ; p < ysetlen ; p++)
                {
                    Int inew = Yseti [p] ;
                    Int iold = Perm ? Perm [inew] : inew ;
                    Xx [iold] = Yx [inew] ;
                    Xz [iold] = Yz [inew] ;
                    Xseti [p] = iold ;
                }
                break ;
        }

        Xsetp [0] = 0 ;
        Xsetp [1] = ysetlen ;

        DEBUG (CHOLMOD(dump_sparse) (Xset, "Xset", Common)) ;
        DEBUG (CHOLMOD(dump_dense) (X, "X", Common)) ;
        Common->no_workspace_reallocate = save_realloc_state ;
        /* done using Iwork (n:3n-1) for Ci and Yseti ] */

    }
    else if (sys == CHOLMOD_P)
    {

	/* ------------------------------------------------------------------ */
	/* x = P*b */
	/* ------------------------------------------------------------------ */

	perm (B, Perm, 0, nrhs, X) ;

    }
    else if (sys == CHOLMOD_Pt)
    {

	/* ------------------------------------------------------------------ */
	/* x = P'*b */
	/* ------------------------------------------------------------------ */

	iperm (B, Perm, 0, nrhs, X) ;

    }
    else if (L->is_super)
    {

	/* ------------------------------------------------------------------ */
	/* solve using a supernodal LL' factorization */
	/* ------------------------------------------------------------------ */

#ifndef NSUPERNODAL
	/* allocate workspace */
	cholmod_dense *E ;
	Int dual ;
        Common->blas_ok = TRUE ;
	dual = (L->xtype == CHOLMOD_REAL && B->xtype != CHOLMOD_REAL) ? 2 : 1 ;
	Y = CHOLMOD(ensure_dense) (Y_Handle, n, dual*nrhs, n, L->xtype, Common);
	E = CHOLMOD(ensure_dense) (E_Handle, dual*nrhs, L->maxesize, dual*nrhs,
		L->xtype, Common) ;

	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
            return (FALSE) ;
	}

	perm (B, Perm, 0, nrhs, Y) ;			    /* Y = P*B */

	if (sys == CHOLMOD_A || sys == CHOLMOD_LDLt)
	{
	    CHOLMOD(super_lsolve) (L, Y, E, Common) ;	    /* Y = L\Y */
            CHOLMOD(super_ltsolve) (L, Y, E, Common) ;	    /* Y = L'\Y*/
	}
	else if (sys == CHOLMOD_L || sys == CHOLMOD_LD)
	{
	    CHOLMOD(super_lsolve) (L, Y, E, Common) ;	    /* Y = L\Y */
	}
	else if (sys == CHOLMOD_Lt || sys == CHOLMOD_DLt)
	{
	    CHOLMOD(super_ltsolve) (L, Y, E, Common) ;      /* Y = L'\Y*/
	}

	iperm (Y, Perm, 0, nrhs, X) ;			    /* X = P'*Y */

	if (CHECK_BLAS_INT && !Common->blas_ok)
	{
	    /* Integer overflow in the BLAS.  This is probably impossible,
	     * since the BLAS were used to create the supernodal factorization.
	     * It might be possible for the calls to the BLAS to differ between
	     * factorization and forward/backsolves, however.  This statement
	     * is untested; it does not appear in the compiled code if
             * CHECK_BLAS_INT is true (when the same integer is used in
             * CHOLMOD and the BLAS. */
	    return (FALSE) ;
	}

#else
	/* CHOLMOD Supernodal module not installed */
	ERROR (CHOLMOD_NOT_INSTALLED,"Supernodal module not installed") ;
#endif

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* solve using a simplicial LL' or LDL' factorization */
	/* ------------------------------------------------------------------ */

        if (L->xtype == CHOLMOD_REAL && B->xtype == CHOLMOD_REAL)
	{
	    /* L, B, and Y are all real */
	    /* solve with up to 4 columns of B at a time */
            ncols = 4 ;
            nr = MAX (4, nrhs) ;
	    ytype = CHOLMOD_REAL ;
	}
	else if (L->xtype == CHOLMOD_REAL)
	{
            /* L is real and B is complex or zomplex */
	    /* solve with one column of B (real/imag), at a time */
	    ncols = 1 ;
	    nr = 2 ;
	    ytype = CHOLMOD_REAL ;
	}
	else
	{
	    /* L is complex or zomplex, B is real/complex/zomplex, Y has the
	     * same complexity as L.  Solve with one column of B at a time. */
	    ncols = 1 ;
	    nr = 1 ;
	    ytype = L->xtype ;
	}

	Y = CHOLMOD(ensure_dense) (Y_Handle, nr, n, nr, ytype, Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
	    return (FALSE) ;
	}

        for (k1 = 0 ; k1 < nrhs ; k1 += ncols)
        {

            /* -------------------------------------------------------------- */
            /* Y = B (P, k1:k1+ncols-1)' = (P * B (:,...))' */
            /* -------------------------------------------------------------- */

            ptrans (B, Perm, k1, ncols, Y) ;

            /* -------------------------------------------------------------- */
            /* solve Y = (L' \ (L \ Y'))', or other system, with template */
            /* -------------------------------------------------------------- */

            switch (L->xtype)
            {
                case CHOLMOD_REAL:
                    r_simplicial_solver (sys, L, Y, NULL, 0) ;
                    break ;

                case CHOLMOD_COMPLEX:
                    c_simplicial_solver (sys, L, Y, NULL, 0) ;
                    break ;

                case CHOLMOD_ZOMPLEX:
                    z_simplicial_solver (sys, L, Y, NULL, 0) ;
                    break ;
            }

            /* -------------------------------------------------------------- */
            /* X (P, k1:k2+ncols-1) = Y' */
            /* -------------------------------------------------------------- */

            iptrans (Y, Perm, k1, ncols, X) ;
        }
    }

    /*
    printf ("bye from solve2\n") ;
    */
    DEBUG (CHOLMOD(dump_dense) (X, "X result", Common)) ;
    return (TRUE) ;
}
#endif
