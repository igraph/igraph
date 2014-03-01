/* ========================================================================== */
/* === Core/cholmod_add ===================================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2006,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* C = alpha*A + beta*B, or spones(A+B).  Result is packed, with sorted or
 * unsorted columns.  This routine is much faster and takes less memory if C
 * is allowed to have unsorted columns.
 *
 * If A and B are both symmetric (in upper form) then C is the same.  Likewise,
 * if A and B are both symmetric (in lower form) then C is the same.
 * Otherwise, C is unsymmetric.  A and B must have the same dimension.
 *
 * workspace: Flag (nrow), W (nrow) if values, Iwork (max (nrow,ncol)).
 *	allocates temporary copies for A and B if they are symmetric.
 *	allocates temporary copy of C if it is to be returned sorted.
 *
 * A and B can have an xtype of pattern or real.  Complex or zomplex cases
 * are supported only if the "values" input parameter is FALSE.
 */

#include "cholmod_internal.h"
#include "cholmod_core.h"

cholmod_sparse *CHOLMOD(add)
(
    /* ---- input ---- */
    cholmod_sparse *A,	    /* matrix to add */
    cholmod_sparse *B,	    /* matrix to add */
    double alpha [2],	    /* scale factor for A */
    double beta [2],	    /* scale factor for B */
    int values,		    /* if TRUE compute the numerical values of C */
    int sorted,		    /* if TRUE, sort columns of C */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Ax, *Bx, *Cx, *W ;
    Int apacked, up, lo, nrow, ncol, bpacked, nzmax, pa, paend, pb, pbend, i,
	j, p, mark, nz ;
    Int *Ap, *Ai, *Anz, *Bp, *Bi, *Bnz, *Flag, *Cp, *Ci ;
    cholmod_sparse *A2, *B2, *C ;

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
    if (A->nrow != B->nrow || A->ncol != B->ncol)
    {
	/* A and B must have the same dimensions */
	ERROR (CHOLMOD_INVALID, "A and B dimesions do not match") ;
	return (NULL) ;
    }
    /* A and B must have the same numerical type if values is TRUE (both must
     * be CHOLMOD_REAL, this is implicitly checked above) */

    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    ncol = A->ncol ;
    CHOLMOD(allocate_work) (nrow, MAX (nrow,ncol), values ? nrow : 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    if (nrow <= 1)
    {
	/* C will be implicitly sorted, so no need to sort it here */
	sorted = FALSE ;
    }

    /* convert A or B to unsymmetric, if necessary */
    A2 = NULL ;
    B2 = NULL ;

    if (A->stype != B->stype)
    {
	if (A->stype)
	{
	    /* workspace: Iwork (max (nrow,ncol)) */
	    A2 = CHOLMOD(copy) (A, 0, values, Common) ;
	    if (Common->status < CHOLMOD_OK)
	    {
		return (NULL) ;	    /* out of memory */
	    }
	    A = A2 ;
	}
	if (B->stype)
	{
	    /* workspace: Iwork (max (nrow,ncol)) */
	    B2 = CHOLMOD(copy) (B, 0, values, Common) ;
	    if (Common->status < CHOLMOD_OK)
	    {
		CHOLMOD(free_sparse) (&A2, Common) ;
		return (NULL) ;	    /* out of memory */
	    }
	    B = B2 ;
	}
    }

    /* get the A matrix */
    ASSERT (A->stype == B->stype) ;
    up = (A->stype > 0) ;
    lo = (A->stype < 0) ;

    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    apacked = A->packed ;

    /* get the B matrix */
    Bp  = B->p ;
    Bnz = B->nz ;
    Bi  = B->i ;
    Bx  = B->x ;
    bpacked = B->packed ;

    /* get workspace */
    W = Common->Xwork ;	    /* size nrow, used if values is TRUE */
    Flag = Common->Flag ;   /* size nrow, Flag [0..nrow-1] < mark on input */

    /* ---------------------------------------------------------------------- */
    /* allocate the result C */
    /* ---------------------------------------------------------------------- */

    /* If integer overflow occurs, nzmax < 0 and the allocate fails properly
     * (likewise in most other matrix manipulation routines). */

    nzmax = CHOLMOD(nnz) (A, Common) + CHOLMOD(nnz) (B, Common) ;

    C = CHOLMOD(allocate_sparse) (nrow, ncol, nzmax, FALSE, TRUE,
	    SIGN (A->stype), values ? A->xtype : CHOLMOD_PATTERN, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	CHOLMOD(free_sparse) (&A2, Common) ;
	CHOLMOD(free_sparse) (&B2, Common) ;
	return (NULL) ;	    /* out of memory */
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* ---------------------------------------------------------------------- */
    /* compute C = alpha*A + beta*B */
    /* ---------------------------------------------------------------------- */

    nz = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	Cp [j] = nz ;

	/* clear the Flag array */
	/* mark = CHOLMOD(clear_flag) (Common) ; */
	CHOLMOD_CLEAR_FLAG (Common) ;
	mark = Common->mark ;

	/* scatter B into W */
	pb = Bp [j] ;
	pbend = (bpacked) ? (Bp [j+1]) : (pb + Bnz [j]) ;
	for (p = pb ; p < pbend ; p++)
	{
	    i = Bi [p] ;
	    if ((up && i > j) || (lo && i < j))
	    {
		continue ;
	    }
	    Flag [i] = mark ;
	    if (values)
	    {
		W [i] = beta [0] * Bx [p] ;
	    }
	}

	/* add A and gather from W into C(:,j) */
	pa = Ap [j] ;
	paend = (apacked) ? (Ap [j+1]) : (pa + Anz [j]) ;
	for (p = pa ; p < paend ; p++)
	{
	    i = Ai [p] ;
	    if ((up && i > j) || (lo && i < j))
	    {
		continue ;
	    }
	    Flag [i] = EMPTY ;
	    Ci [nz] = i ;
	    if (values)
	    {
		Cx [nz] = W [i] + alpha [0] * Ax [p] ;
		W [i] = 0 ;
	    }
	    nz++ ;
	}

	/* gather remaining entries into C(:,j), using pattern of B */
	for (p = pb ; p < pbend ; p++)
	{
	    i = Bi [p] ;
	    if ((up && i > j) || (lo && i < j))
	    {
		continue ;
	    }
	    if (Flag [i] == mark)
	    {
		Ci [nz] = i ;
		if (values)
		{
		    Cx [nz] = W [i] ;
		    W [i] = 0 ;
		}
		nz++ ;
	    }
	}
    }

    Cp [ncol] = nz ;

    /* ---------------------------------------------------------------------- */
    /* reduce C in size and free temporary matrices */
    /* ---------------------------------------------------------------------- */

    ASSERT (MAX (1,nz) <= C->nzmax) ;
    CHOLMOD(reallocate_sparse) (nz, C, Common) ;
    ASSERT (Common->status >= CHOLMOD_OK) ;

    /* clear the Flag array */
    mark = CHOLMOD(clear_flag) (Common) ;

    CHOLMOD(free_sparse) (&A2, Common) ;
    CHOLMOD(free_sparse) (&B2, Common) ;

    /* ---------------------------------------------------------------------- */
    /* sort C, if requested */
    /* ---------------------------------------------------------------------- */

    if (sorted)
    {
	/* workspace: Iwork (max (nrow,ncol)) */
	if (!CHOLMOD(sort) (C, Common))
	{
	    CHOLMOD(free_sparse) (&C, Common) ;
	    if (Common->status < CHOLMOD_OK)
	    {
		return (NULL) ;		/* out of memory */
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    ASSERT (CHOLMOD(dump_sparse) (C, "add", Common) >= 0) ;
    return (C) ;
}
