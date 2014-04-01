/* ========================================================================== */
/* === MatrixOps/cholmod_vertcat ============================================ */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/MatrixOps Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/MatrixOps Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* Vertical concatenation, C = [A ; B] in MATLAB notation.
 *
 * A and B can be up/lo/unsym; C is unsymmetric and packed.
 * A and B must have the same number of columns.
 * C is sorted if both A and B are sorted.
 *
 * workspace: Iwork (max (A->nrow, A->ncol, B->nrow, B->ncol)).
 *	allocates temporary copies of A and B if they are symmetric.
 *
 * Only pattern and real matrices are supported.  Complex and zomplex matrices
 * are supported only if "values" is FALSE.
 */

#ifndef NMATRIXOPS

#include "cholmod_internal.h"
#include "cholmod_matrixops.h"


/* ========================================================================== */
/* === cholmod_vertcat ====================================================== */
/* ========================================================================== */

cholmod_sparse *CHOLMOD(vertcat)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* left matrix to concatenate */
    cholmod_sparse *B,	/* right matrix to concatenate */
    int values,		/* if TRUE compute the numerical values of C */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Ax, *Bx, *Cx ;
    Int *Ap, *Ai, *Anz, *Bp, *Bi, *Bnz, *Cp, *Ci ;
    cholmod_sparse *C, *A2, *B2 ;
    Int apacked, bpacked, anrow, bnrow, ncol, nrow, anz, bnz, nz, j, p, pend,
	pdest ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    RETURN_IF_NULL (B, NULL) ;
    values = values &&
	(A->xtype != CHOLMOD_PATTERN) && (B->xtype != CHOLMOD_PATTERN) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN,
	    values ? CHOLMOD_REAL : CHOLMOD_ZOMPLEX, NULL) ;
    RETURN_IF_XTYPE_INVALID (B, CHOLMOD_PATTERN,
	    values ? CHOLMOD_REAL : CHOLMOD_ZOMPLEX, NULL) ;
    if (A->ncol != B->ncol)
    {
	/* A and B must have the same number of columns */
	ERROR (CHOLMOD_INVALID, "A and B must have same # of columns") ;
	return (NULL) ;
    }
    /* A and B must have the same numerical type if values is TRUE (both must
     * be CHOLMOD_REAL, this is implicitly checked above) */
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    anrow = A->nrow ;
    bnrow = B->nrow ;
    ncol = A->ncol ;
    CHOLMOD(allocate_work) (0, MAX3 (anrow, bnrow, ncol), 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    /* convert A to unsymmetric, if necessary */
    A2 = NULL ;
    if (A->stype != 0)
    {
	/* workspace: Iwork (max (A->nrow,A->ncol)) */
	A2 = CHOLMOD(copy) (A, 0, values, Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
	    return (NULL) ;
	}
	A = A2 ;
    }

    /* convert B to unsymmetric, if necessary */
    B2 = NULL ;
    if (B->stype != 0)
    {
	/* workspace: Iwork (max (B->nrow,B->ncol)) */
	B2 = CHOLMOD(copy) (B, 0, values, Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
	    CHOLMOD(free_sparse) (&A2, Common) ;
	    return (NULL) ;
	}
	B = B2 ;
    }

    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    apacked = A->packed ;

    Bp  = B->p ;
    Bnz = B->nz ;
    Bi  = B->i ;
    Bx  = B->x ;
    bpacked = B->packed ;

    /* ---------------------------------------------------------------------- */
    /* allocate C */
    /* ---------------------------------------------------------------------- */

    anz = CHOLMOD(nnz) (A, Common) ;
    bnz = CHOLMOD(nnz) (B, Common) ;
    nrow = anrow + bnrow ;
    nz = anz + bnz ;

    C = CHOLMOD(allocate_sparse) (nrow, ncol, nz, A->sorted && B->sorted, TRUE,
	    0, values ? A->xtype : CHOLMOD_PATTERN, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	CHOLMOD(free_sparse) (&A2, Common) ;
	CHOLMOD(free_sparse) (&B2, Common) ;
	return (NULL) ;
    }
    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* ---------------------------------------------------------------------- */
    /* C = [A ; B] */
    /* ---------------------------------------------------------------------- */

    pdest = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	/* attach A(:,j) as the first part of C(:,j) */
	p = Ap [j] ;
	pend = (apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
	Cp [j] = pdest ;
	for ( ; p < pend ; p++)
	{
	    Ci [pdest] = Ai [p] ;
	    if (values)
	    {
		Cx [pdest] = Ax [p] ;
	    }
	    pdest++ ;
	}

	/* attach B(:,j) as the second part of C(:,j) */
	p = Bp [j] ;
	pend = (bpacked) ? (Bp [j+1]) : (p + Bnz [j]) ;
	for ( ; p < pend ; p++)
	{
	    Ci [pdest] = Bi [p] + anrow ;
	    if (values)
	    {
		Cx [pdest] = Bx [p] ;
	    }
	    pdest++ ;
	}
    }
    Cp [ncol] = pdest ;
    ASSERT (pdest == nz) ;

    /* ---------------------------------------------------------------------- */
    /* free the unsymmetric copies of A and B, and return C */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(free_sparse) (&A2, Common) ;
    CHOLMOD(free_sparse) (&B2, Common) ;
    return (C) ;
}
#endif
