/* ========================================================================== */
/* === MatrixOps/cholmod_scale ============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/MatrixOps Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/MatrixOps Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* scale a matrix:  A = diag(s)*A, A*diag(s), s*A, or diag(s)*A*diag(s)
 *
 * A can be of any type (packed/unpacked, upper/lower/unsymmetric).
 * The symmetry of A is ignored; all entries in the matrix are modified.
 *
 * If A is m-by-n unsymmetric but scaled symmtrically, the result is
 * A = diag (s (1:m)) * A * diag (s (1:n)).
 *
 * Note: diag(s) should be interpretted as spdiags(s,0,n,n) where n=length(s).
 *
 * Row or column scaling of a symmetric matrix still results in a symmetric
 * matrix, since entries are still ignored by other routines.
 * For example, when row-scaling a symmetric matrix where just the upper
 * triangular part is stored (and lower triangular entries ignored)
 * A = diag(s)*triu(A) is performed, where the result A is also
 * symmetric-upper.  This has the effect of modifying the implicit lower
 * triangular part.  In MATLAB notation:
 *
 *	U = diag(s)*triu(A) ;
 *	L = tril (U',-1)
 *	A = L + U ;
 *
 * The scale parameter determines the kind of scaling to perform:
 *
 *	 CHOLMOD_SCALAR: s[0]*A
 *	 CHOLMOD_ROW:    diag(s)*A
 *	 CHOLMOD_COL:    A*diag(s)
 *	 CHOLMOD_SYM:    diag(s)*A*diag(s)
 *
 * The size of S depends on the scale parameter:
 *
 *	 CHOLMOD_SCALAR: size 1
 *	 CHOLMOD_ROW:    size nrow-by-1 or 1-by-nrow
 *	 CHOLMOD_COL:    size ncol-by-1 or 1-by-ncol
 *	 CHOLMOD_SYM:    size max(nrow,ncol)-by-1, or 1-by-max(nrow,ncol)
 *
 * workspace: none
 *
 * Only real matrices are supported.
 */

#ifndef NMATRIXOPS

#include "cholmod_internal.h"
#include "cholmod_matrixops.h"


/* ========================================================================== */
/* === cholmod_scale ======================================================== */
/* ========================================================================== */

int CHOLMOD(scale)
(
    /* ---- input ---- */
    cholmod_dense *S,	/* scale factors (scalar or vector) */
    int scale,		/* type of scaling to compute */
    /* ---- in/out --- */
    cholmod_sparse *A,	/* matrix to scale */
    /* --------------- */
    cholmod_common *Common
)
{
    double t ;
    double *Ax, *s ;
    Int *Ap, *Anz, *Ai ;
    Int packed, j, ncol, nrow, p, pend, sncol, snrow, nn, ok ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (S, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_REAL, CHOLMOD_REAL, FALSE) ;
    RETURN_IF_XTYPE_INVALID (S, CHOLMOD_REAL, CHOLMOD_REAL, FALSE) ;
    ncol = A->ncol ;
    nrow = A->nrow ;
    sncol = S->ncol ;
    snrow = S->nrow ;
    if (scale == CHOLMOD_SCALAR)
    {
	ok = (snrow == 1 && sncol == 1) ;
    }
    else if (scale == CHOLMOD_ROW)
    {
	ok = (snrow == nrow && sncol == 1) || (snrow == 1 && sncol == nrow) ;
    }
    else if (scale == CHOLMOD_COL)
    {
	ok = (snrow == ncol && sncol == 1) || (snrow == 1 && sncol == ncol) ;
    }
    else if (scale == CHOLMOD_SYM)
    {
	nn = MAX (nrow, ncol) ;
	ok = (snrow == nn && sncol == 1) || (snrow == 1 && sncol == nn) ;
    }
    else
    {
	/* scale invalid */
	ERROR (CHOLMOD_INVALID, "invalid scaling option") ;
	return (FALSE) ;
    }
    if (!ok)
    {
	/* S is wrong size */
	ERROR (CHOLMOD_INVALID, "invalid scale factors") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    packed = A->packed ;
    s = S->x ;

    /* ---------------------------------------------------------------------- */
    /* scale the matrix */
    /* ---------------------------------------------------------------------- */

    if (scale == CHOLMOD_ROW)
    {

	/* ------------------------------------------------------------------ */
	/* A = diag(s)*A, row scaling */
	/* ------------------------------------------------------------------ */

	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		Ax [p] *= s [Ai [p]] ;
	    }
	}

    }
    else if (scale == CHOLMOD_COL)
    {

	/* ------------------------------------------------------------------ */
	/* A = A*diag(s), column scaling */
	/* ------------------------------------------------------------------ */

	for (j = 0 ; j < ncol ; j++)
	{
	    t = s [j] ;
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		Ax [p] *= t ;
	    }
	}

    }
    else if (scale == CHOLMOD_SYM)
    {

	/* ------------------------------------------------------------------ */
	/* A = diag(s)*A*diag(s), symmetric scaling */
	/* ------------------------------------------------------------------ */

	for (j = 0 ; j < ncol ; j++)
	{
	    t = s [j] ;
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		Ax [p] *= t * s [Ai [p]] ;
	    }
	}

    }
    else if (scale == CHOLMOD_SCALAR)
    {

	/* ------------------------------------------------------------------ */
	/* A = s[0] * A, scalar scaling */
	/* ------------------------------------------------------------------ */

	t = s [0] ;
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		Ax [p] *= t ;
	    }
	}
    }

    ASSERT (CHOLMOD(dump_sparse) (A, "A scaled", Common) >= 0) ;
    return (TRUE) ;
}
#endif
