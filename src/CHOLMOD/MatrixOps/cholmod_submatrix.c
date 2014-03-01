/* ========================================================================== */
/* === MatrixOps/cholmod_submatrix ========================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/MatrixOps Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/MatrixOps Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* C = A (rset,cset), where C becomes length(rset)-by-length(cset) in dimension.
 * rset and cset can have duplicate entries.  A and C must be unsymmetric.   C
 * is packed.  If the sorted flag is TRUE on input, or rset is sorted and A is
 * sorted, then C is sorted; otherwise C is unsorted.
 *
 * A NULL rset or cset means "[ ]" in MATLAB notation.
 * If the length of rset or cset is negative, it denotes ":" in MATLAB notation.
 *
 * For permuting a matrix, this routine is an alternative to cholmod_ptranspose
 * (which permutes and transposes a matrix and can work on symmetric matrices).
 *
 * The time taken by this routine is O(A->nrow) if the Common workspace needs
 * to be initialized, plus O(C->nrow + C->ncol + nnz (A (:,cset))).  Thus, if C
 * is small and the workspace is not initialized, the time can be dominated by
 * the call to cholmod_allocate_work.  However, once the workspace is
 * allocated, subsequent calls take less time.
 *
 * workspace:  Iwork (max (A->nrow + length (rset), length (cset))).
 *	allocates temporary copy of C if it is to be returned sorted.
 *
 * Future work:  A common case occurs where A has sorted columns, and rset is in
 * the form lo:hi in MATLAB notation.  This routine could exploit that case
 * to run even faster when the matrix is sorted, particularly when lo is small.
 *
 * Only pattern and real matrices are supported.  Complex and zomplex matrices
 * are supported only when "values" is FALSE.
 */

#ifndef NMATRIXOPS

#include "cholmod_internal.h"
#include "cholmod_matrixops.h"

/* ========================================================================== */
/* === check_subset ========================================================= */
/* ========================================================================== */

/* Check the rset or cset, and return TRUE if valid, FALSE if invalid */

static int check_subset (Int *set, Int len, Int n)
{
    Int k ;
    if (set == NULL)
    {
	return (TRUE) ;
    }
    for (k = 0 ; k < len ; k++)
    {
	if (set [k] < 0 || set [k] >= n)
	{
	    return (FALSE) ;
	}
    }
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_submatrix ==================================================== */
/* ========================================================================== */

cholmod_sparse *CHOLMOD(submatrix)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to subreference */
    Int *rset,		/* set of row indices, duplicates OK */
    SuiteSparse_long rsize,	/* size of rset, or -1 for ":" */
    Int *cset,		/* set of column indices, duplicates OK */
    SuiteSparse_long csize,	/* size of cset, or -1 for ":" */
    int values,		/* if TRUE compute the numerical values of C */
    int sorted,		/* if TRUE then return C with sorted columns */
    /* --------------- */
    cholmod_common *Common
)
{
    double aij = 0 ;
    double *Ax, *Cx ;
    Int *Ap, *Ai, *Anz, *Ci, *Cp, *Head, *Rlen, *Rnext, *Iwork ;
    cholmod_sparse *C ;
    Int packed, ancol, anrow, cnrow, cncol, nnz, i, j, csorted, ilast, p,
	pend, pdest, ci, cj, head, nr, nc ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    values = (values && (A->xtype != CHOLMOD_PATTERN)) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN,
	    values ? CHOLMOD_REAL : CHOLMOD_ZOMPLEX, NULL) ;
    if (A->stype != 0)
    {
	/* A must be unsymmetric */
	ERROR (CHOLMOD_INVALID, "symmetric upper or lower case not supported") ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    ancol = A->ncol ;
    anrow = A->nrow ;
    nr = rsize ;
    nc = csize ;
    if (rset == NULL)
    {
	/* nr = 0 denotes rset = [ ], nr < 0 denotes rset = 0:anrow-1 */
	nr = (nr < 0) ? (-1) : 0 ;
    }
    if (cset == NULL)
    {
	/* nr = 0 denotes cset = [ ], nr < 0 denotes cset = 0:ancol-1 */
	nc = (nc < 0) ? (-1) : 0 ;
    }
    cnrow = (nr < 0) ? anrow : nr ;  /* negative rset means rset = 0:anrow-1 */
    cncol = (nc < 0) ? ancol : nc ;  /* negative cset means cset = 0:ancol-1 */

    if (nr < 0 && nc < 0)
    {

	/* ------------------------------------------------------------------ */
	/* C = A (:,:), use cholmod_copy instead */
	/* ------------------------------------------------------------------ */

	/* workspace: Iwork (max (C->nrow,C->ncol)) */
	PRINT1 (("submatrix C = A (:,:)\n")) ;
	C = CHOLMOD(copy) (A, 0, values, Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
	    return (NULL) ;
	}
	return (C) ;
    }
    PRINT1 (("submatrix nr "ID" nc "ID" Cnrow "ID" Cncol "ID""
	    "  Anrow "ID" Ancol "ID"\n", nr, nc, cnrow, cncol, anrow, ancol)) ;

    /* s = MAX3 (anrow+MAX(0,nr), cncol, cnrow) ; */
    s = CHOLMOD(add_size_t) (anrow, MAX (0,nr), &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (NULL) ;
    }
    s = MAX3 (s, ((size_t) cncol), ((size_t) cnrow)) ;

    CHOLMOD(allocate_work) (anrow, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	return (NULL) ;
    }

    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    packed = A->packed ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Head  = Common->Head ;	    /* size anrow */
    Iwork = Common->Iwork ;
    Rlen  = Iwork ;		    /* size anrow (i/i/l) */
    Rnext = Iwork + anrow ;	    /* size nr (i/i/l), not used if nr < 0 */

    /* ---------------------------------------------------------------------- */
    /* construct inverse of rset and compute nnz (C) */
    /* ---------------------------------------------------------------------- */

    PRINT1 (("nr "ID" nc "ID"\n", nr, nc)) ;
    PRINT1 (("anrow "ID" ancol "ID"\n", anrow, ancol)) ;
    PRINT1 (("cnrow "ID" cncol "ID"\n", cnrow, cncol)) ;
    DEBUG (for (i = 0 ; i < nr ; i++) PRINT2 (("rset ["ID"] = "ID"\n",
		    i, rset [i])));
    DEBUG (for (i = 0 ; i < nc ; i++) PRINT2 (("cset ["ID"] = "ID"\n",
		    i, cset [i])));

    /* C is sorted if A and rset are sorted, or if C has one row or less */
    csorted = A->sorted || (cnrow <= 1) ;

    if (!check_subset (rset, nr, anrow))
    {
	ERROR (CHOLMOD_INVALID, "invalid rset") ;
	return (NULL) ;
    }

    if (!check_subset (cset, nc, ancol))
    {
	ERROR (CHOLMOD_INVALID, "invalid cset") ;
	return (NULL) ;
    }

    nnz = 0 ;
    if (nr < 0)
    {
	/* C = A (:,cset) where cset = [ ] or cset is not empty */
	ASSERT (IMPLIES (cncol > 0, cset != NULL)) ;
	for (cj = 0 ; cj < cncol ; cj++)
	{
	    /* construct column cj of C, which is column j of A */
	    j = cset [cj] ;
	    nnz += (packed) ? (Ap [j+1] - Ap [j]) : MAX (0, Anz [j]) ;
	}
    }
    else
    {
	/* C = A (rset,cset), where rset is not empty but cset might be empty */
	/* create link lists in reverse order to preserve natural order */
	ilast = anrow ;
	for (ci = nr-1 ; ci >= 0 ; ci--)
	{
	    /* row i of A becomes row ci of C; add ci to ith link list */
	    i = rset [ci] ;
	    head = Head [i] ;
	    Rlen [i] = (head == EMPTY) ? 1 : (Rlen [i] + 1) ;
	    Rnext [ci] = head ;
	    Head [i] = ci ;
	    if (i > ilast)
	    {
		/* row indices in columns of C will not be sorted */
		csorted = FALSE ;
	    }
	    ilast = i ;
	}

#ifndef NDEBUG
	for (i = 0 ; i < anrow ; i++)
	{
	    Int k = 0 ;
	    Int rlen = (Head [i] != EMPTY) ? Rlen [i] : -1 ;
	    PRINT1 (("Row "ID" Rlen "ID": ", i, rlen)) ;
	    for (ci = Head [i] ; ci != EMPTY ; ci = Rnext [ci])
	    {
		k++ ;
		PRINT2 ((""ID" ", ci)) ;
	    }
	    PRINT1 (("\n")) ;
	    ASSERT (IMPLIES (Head [i] != EMPTY, k == Rlen [i])) ;
	}
#endif

	/* count nonzeros in C */
	for (cj = 0 ; cj < cncol ; cj++)
	{
	    /* count rows in column cj of C, which is column j of A */
	    j = (nc < 0) ? cj : (cset [cj]) ;
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		/* row i of A becomes multiple rows (ci) of C */
		i = Ai [p] ;
		ASSERT (i >= 0 && i < anrow) ;
		if (Head [i] != EMPTY)
		{
		    nnz += Rlen [i] ;
		}
	    }
	}
    }
    PRINT1 (("nnz (C) "ID"\n", nnz)) ;

    /* rset and cset are now valid */
    DEBUG (CHOLMOD(dump_subset) (rset, rsize, anrow, "rset", Common)) ;
    DEBUG (CHOLMOD(dump_subset) (cset, csize, ancol, "cset", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* allocate C */
    /* ---------------------------------------------------------------------- */

    C = CHOLMOD(allocate_sparse) (cnrow, cncol, nnz, csorted, TRUE, 0,
	    values ? A->xtype : CHOLMOD_PATTERN, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	for (i = 0 ; i < anrow ; i++)
	{
	    Head [i] = EMPTY ;
	}
	ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;
	return (NULL) ;
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* ---------------------------------------------------------------------- */
    /* C = A (rset,cset) */
    /* ---------------------------------------------------------------------- */

    pdest = 0 ;
    if (nnz == 0)
    {
	/* C has no nonzeros */
	for (cj = 0 ; cj <= cncol ; cj++)
	{
	    Cp [cj] = 0 ;
	}
    }
    else if (nr < 0)
    {
	/* C = A (:,cset), where cset is not empty */
	for (cj = 0 ; cj < cncol ; cj++)
	{
	    /* construct column cj of C, which is column j of A */
	    PRINT1 (("construct cj = j = "ID"\n", cj)) ;
	    j = cset [cj] ;
	    Cp [cj] = pdest ;
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		Ci [pdest] = Ai [p] ;
		if (values)
		{
		    Cx [pdest] = Ax [p] ;
		}
		pdest++ ;
		ASSERT (pdest <= nnz) ;
	    }
	}
    }
    else
    {
	/* C = A (rset,cset), where rset is not empty but cset might be empty */
	for (cj = 0 ; cj < cncol ; cj++)
	{
	    /* construct column cj of C, which is column j of A */
	    PRINT1 (("construct cj = "ID"\n", cj)) ;
	    j = (nc < 0) ? cj : (cset [cj]) ;
	    PRINT1 (("cj = "ID"\n", j)) ;
	    Cp [cj] = pdest ;
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		/* row (Ai [p]) of A becomes multiple rows (ci) of C */
		PRINT2 (("i: "ID" becomes: ", Ai [p])) ;
		if (values)
		{
		    aij = Ax [p] ;
		}
		for (ci = Head [Ai [p]] ; ci != EMPTY ; ci = Rnext [ci])
		{
		    PRINT3 ((""ID" ", ci)) ;
		    Ci [pdest] = ci ;
		    if (values)
		    {
			Cx [pdest] = aij ;
		    }
		    pdest++ ;
		    ASSERT (pdest <= nnz) ;
		}
		PRINT2 (("\n")) ;
	    }
	}
    }
    Cp [cncol] = pdest ;
    ASSERT (nnz == pdest) ;

    /* ---------------------------------------------------------------------- */
    /* clear workspace */
    /* ---------------------------------------------------------------------- */

    for (ci = 0 ; ci < nr ; ci++)
    {
	Head [rset [ci]] = EMPTY ;
    }

    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* sort C, if requested */
    /* ---------------------------------------------------------------------- */

    ASSERT (CHOLMOD(dump_sparse) (C , "C before sort", Common) >= 0) ;

    if (sorted && !csorted)
    {
	/* workspace: Iwork (max (C->nrow,C->ncol)) */
	if (!CHOLMOD(sort) (C, Common))
	{
	    /* out of memory */
	    CHOLMOD(free_sparse) (&C, Common) ;
	    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;
	    return (NULL) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    ASSERT (CHOLMOD(dump_sparse) (C , "Final C", Common) >= 0) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;
    return (C) ;
}
#endif
