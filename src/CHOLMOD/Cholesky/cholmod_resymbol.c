/* ========================================================================== */
/* === Cholesky/cholmod_resymbol ============================================ */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Cholesky Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/Cholesky Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Recompute the symbolic pattern of L.  Entries not in the symbolic pattern
 * are dropped.  L->Perm can be used (or not) to permute the input matrix A.
 *
 * These routines are used after a supernodal factorization is converted into
 * a simplicial one, to remove zero entries that were added due to relaxed
 * supernode amalgamation.  They can also be used after a series of downdates
 * to remove entries that would no longer be present if the matrix were
 * factorized from scratch.  A downdate (cholmod_updown) does not remove any
 * entries from L.
 *
 * workspace: Flag (nrow), Head (nrow+1),
 *	if symmetric:   Iwork (2*nrow)
 *	if unsymmetric: Iwork (2*nrow+ncol).
 *	Allocates up to 2 copies of its input matrix A (pattern only).
 */

#ifndef NCHOLESKY

#include "cholmod_internal.h"
#include "cholmod_cholesky.h"

/* ========================================================================== */
/* === cholmod_resymbol ===================================================== */
/* ========================================================================== */

/* Remove entries from L that are not in the factorization of P*A*P', P*A*A'*P',
 * or P*F*F'*P' (depending on A->stype and whether fset is NULL or not).
 *
 * cholmod_resymbol is the same as cholmod_resymbol_noperm, except that it
 * first permutes A according to L->Perm.  A can be upper/lower/unsymmetric,
 * in contrast to cholmod_resymbol_noperm (which can be lower or unsym). */

int CHOLMOD(resymbol)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int pack,		/* if TRUE, pack the columns of L */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factorization, entries pruned on output */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_sparse *H, *F, *G ;
    Int stype, nrow, ncol ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    Common->status = CHOLMOD_OK ;
    if (L->is_super)
    {
	/* cannot operate on a supernodal factorization */
	ERROR (CHOLMOD_INVALID, "cannot operate on supernodal L") ;
	return (FALSE) ;
    }
    if (L->n != A->nrow)
    {
	/* dimensions must agree */
	ERROR (CHOLMOD_INVALID, "A and L dimensions do not match") ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    stype = A->stype ;
    nrow = A->nrow ;
    ncol = A->ncol ;

    /* s = 2*nrow + (stype ? 0 : ncol) */
    s = CHOLMOD(mult_size_t) (nrow, 2, &ok) ;
    s = CHOLMOD(add_size_t) (s, (stype ? 0 : ncol), &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    CHOLMOD(allocate_work) (nrow, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* permute the input matrix if necessary */
    /* ---------------------------------------------------------------------- */

    H = NULL ;
    G = NULL ;

    if (stype > 0)
    {
	if (L->ordering == CHOLMOD_NATURAL)
	{
	    /* F = triu(A)' */
	    /* workspace: Iwork (nrow) */
	    G = CHOLMOD(ptranspose) (A, 0, NULL, NULL, 0, Common) ;
	}
	else
	{
	    /* F = triu(A(p,p))' */
	    /* workspace: Iwork (2*nrow) */
	    G = CHOLMOD(ptranspose) (A, 0, L->Perm, NULL, 0, Common) ;
	}
	F = G ;
    }
    else if (stype < 0)
    {
	if (L->ordering == CHOLMOD_NATURAL)
	{
	    F = A ;
	}
	else
	{
	    /* G = triu(A(p,p))' */
	    /* workspace: Iwork (2*nrow) */
	    G = CHOLMOD(ptranspose) (A, 0, L->Perm, NULL, 0, Common) ;
	    /* H = G' */
	    /* workspace: Iwork (nrow) */
	    H = CHOLMOD(ptranspose) (G, 0, NULL, NULL, 0, Common) ;
	    F = H ;
	}
    }
    else
    {
	if (L->ordering == CHOLMOD_NATURAL)
	{
	    F = A ;
	}
	else
	{
	    /* G = A(p,f)' */
	    /* workspace: Iwork (nrow if no fset; MAX (nrow,ncol) if fset)*/
	    G = CHOLMOD(ptranspose) (A, 0, L->Perm, fset, fsize, Common) ;
	    /* H = G' */
	    /* workspace: Iwork (ncol) */
	    H = CHOLMOD(ptranspose) (G, 0, NULL, NULL, 0, Common) ;
	    F = H ;
	}
    }

    /* No need to check for failure here.  cholmod_resymbol_noperm will return
     * FALSE if F is NULL. */

    /* ---------------------------------------------------------------------- */
    /* resymbol */
    /* ---------------------------------------------------------------------- */

    ok = CHOLMOD(resymbol_noperm) (F, fset, fsize, pack, L, Common) ;

    /* ---------------------------------------------------------------------- */
    /* free the temporary matrices, if they exist */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(free_sparse) (&H, Common) ;
    CHOLMOD(free_sparse) (&G, Common) ;
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_resymbol_noperm ============================================== */
/* ========================================================================== */

/* Redo symbolic LDL' or LL' factorization of I + F*F' or I+A, where F=A(:,f).
 *
 * L already exists, but is a superset of the true dynamic pattern (simple
 * column downdates and row deletions haven't pruned anything).  Just redo the
 * symbolic factorization and drop entries that are no longer there.  The
 * diagonal is not modified.  The number of nonzeros in column j of L
 * (L->nz[j]) can decrease.  The column pointers (L->p[j]) remain unchanged if
 * pack is FALSE or if L is not monotonic.  Otherwise, the columns of L are
 * packed in place.
 *
 * For the symmetric case, the columns of the lower triangular part of A
 * are accessed by column.  NOTE that this the transpose of the general case.
 *
 * For the unsymmetric case, F=A(:,f) is accessed by column.
 *
 * A need not be sorted, and can be packed or unpacked.  If L->Perm is not
 * identity, then A must already be permuted according to the permutation used
 * to factorize L.  The advantage of using this routine is that it does not
 * need to create permuted copies of A first.
 *
 * This routine can be called if L is only partially factored via cholmod_rowfac
 * since all it does is prune.  If an entry is in F*F' or A, but not in L, it
 * isn't added to L.
 *
 * L must be simplicial LDL' or LL'; it cannot be supernodal or symbolic.
 *
 * The set f is held in fset and fsize.
 *	fset = NULL means ":" in MATLAB. fset is ignored.
 *	fset != NULL means f = fset [0..fset-1].
 *	fset != NULL and fsize = 0 means f is the empty set.
 *	There can be no duplicates in fset.
 *	Common->status is set to CHOLMOD_INVALID if fset is invalid.
 *
 * workspace: Flag (nrow), Head (nrow+1),
 *	if symmetric:   Iwork (2*nrow)
 *	if unsymmetric: Iwork (2*nrow+ncol).
 *	Unlike cholmod_resymbol, this routine does not allocate any temporary
 *	copies of its input matrix.
 */

int CHOLMOD(resymbol_noperm)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int pack,		/* if TRUE, pack the columns of L */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factorization, entries pruned on output */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Lx, *Lz ;
    Int i, j, k, row, parent, p, pend, pdest, ncol, apacked, sorted, nrow, nf,
	use_fset, mark, jj, stype, xtype ;
    Int *Ap, *Ai, *Anz, *Li, *Lp, *Lnz, *Flag, *Head, *Link, *Anext, *Iwork ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    ncol = A->ncol ;
    nrow = A->nrow ;
    stype = A->stype ;
    ASSERT (IMPLIES (stype != 0, nrow == ncol)) ;
    if (stype > 0)
    {
	/* symmetric, with upper triangular part, not supported */
	ERROR (CHOLMOD_INVALID, "symmetric upper not supported ") ;
	return (FALSE) ;
    }
    if (L->is_super)
    {
	/* cannot operate on a supernodal or symbolic factorization */
	ERROR (CHOLMOD_INVALID, "cannot operate on supernodal L") ;
	return (FALSE) ;
    }
    if (L->n != A->nrow)
    {
	/* dimensions must agree */
	ERROR (CHOLMOD_INVALID, "A and L dimensions do not match") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* s = 2*nrow + (stype ? 0 : ncol) */
    s = CHOLMOD(mult_size_t) (nrow, 2, &ok) ;
    if (stype != 0)
    {
	s = CHOLMOD(add_size_t) (s, ncol, &ok) ;
    }
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    CHOLMOD(allocate_work) (nrow, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;	/* out of memory */
    }
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ai = A->i ;
    Ap = A->p ;
    Anz = A->nz ;
    apacked = A->packed ;
    sorted = A->sorted ;

    Li = L->i ;
    Lx = L->x ;
    Lz = L->z ;
    Lp = L->p ;
    Lnz = L->nz ;
    xtype = L->xtype ;

    /* If L is monotonic on input, then it can be packed or
     * unpacked on output, depending on the pack input parameter. */

    /* cannot pack a non-monotonic matrix */
    if (!(L->is_monotonic))
    {
	pack = FALSE ;
    }

    ASSERT (L->nzmax >= (size_t) (Lp [L->n])) ;

    pdest = 0 ;

    PRINT1 (("\n\n===================== Resymbol pack %d Apacked %d\n",
	pack, A->packed)) ;
    ASSERT (CHOLMOD(dump_sparse) (A, "ReSymbol A:", Common) >= 0) ;
    DEBUG (CHOLMOD(dump_factor) (L, "ReSymbol initial L (i, x):", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Flag  = Common->Flag ;	/* size nrow */
    Head  = Common->Head ;	/* size nrow+1 */
    Iwork = Common->Iwork ;
    Link  = Iwork ;		/* size nrow (i/i/l) [ */
    Lnz   = Iwork + nrow ;	/* size nrow (i/i/l), if L not packed */
    Anext = Iwork + 2*((size_t) nrow) ;	/* size ncol (i/i/l), unsym. only */
    for (j = 0 ; j < nrow ; j++)
    {
	Link [j] = EMPTY ;
    }

    /* use Lnz in L itself */
    Lnz = L->nz ;
    ASSERT (Lnz != NULL) ;

    /* ---------------------------------------------------------------------- */
    /* for the unsymmetric case, queue each column of A (:,f) */
    /* ---------------------------------------------------------------------- */

    /* place each column of the basis set on the link list corresponding to */
    /* the smallest row index in that column */

    if (stype == 0)
    {
	use_fset = (fset != NULL) ;
	if (use_fset)
	{
	    nf = fsize ;
	    /* This is the only O(ncol) loop in cholmod_resymbol.
	     * It is required only to check the fset. */
	    for (j = 0 ; j < ncol ; j++)
	    {
		Anext [j] = -2 ;
	    }
	    for (jj = 0 ; jj < nf ; jj++)
	    {
		j = fset [jj] ;
		if (j < 0 || j > ncol || Anext [j] != -2)
		{
		    /* out-of-range or duplicate entry in fset */
		    ERROR (CHOLMOD_INVALID, "fset invalid") ;
		    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;
		    return (FALSE) ;
		}
		/* flag column j as having been seen */
		Anext [j] = EMPTY ;
	    }
	    /* the fset is now valid */
	    ASSERT (CHOLMOD(dump_perm) (fset, nf, ncol, "fset", Common)) ;
	}
	else
	{
	    nf = ncol ;
	}
	for (jj = 0 ; jj < nf ; jj++)
	{
	    j = (use_fset) ? (fset [jj]) : jj ;
	    /* column j is the fset; find the smallest row (if any) */
	    p = Ap [j] ;
	    pend = (apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
	    if (pend > p)
	    {
		k = Ai [p] ;
		if (!sorted)
		{
		    for ( ; p < pend ; p++)
		    {
			k = MIN (k, Ai [p]) ;
		    }
		}
		/* place column j on link list k */
		ASSERT (k >= 0 && k < nrow) ;
		Anext [j] = Head [k] ;
		Head [k] = j ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* recompute symbolic LDL' factorization */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < nrow ; k++)
    {

#ifndef NDEBUG
	PRINT1 (("\n\n================== Initial column k = "ID"\n", k)) ;
	for (p = Lp [k] ; p < Lp [k] + Lnz [k] ; p++)
	{
	    PRINT1 ((" row: "ID"  value: ", Li [p])) ;
	    PRINT1 (("\n")) ;
	}
	PRINT1 (("Recomputing LDL, column k = "ID"\n", k)) ;
#endif

	/* ------------------------------------------------------------------ */
	/* compute column k of I+F*F' or I+A */
	/* ------------------------------------------------------------------ */

	/* flag the diagonal entry */
	/* mark = CHOLMOD(clear_flag) (Common) ; */
	CHOLMOD_CLEAR_FLAG (Common) ;
	mark = Common->mark ;

	Flag [k] = mark ;
	PRINT1 (("	row: "ID" (diagonal)\n", k)) ;

	if (stype != 0)
	{
	    /* merge column k of A into Flag (lower triangular part only) */
	    p = Ap [k] ;
	    pend = (apacked) ? (Ap [k+1]) : (p + Anz [k]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i > k)
		{
		    Flag [i] = mark ;
		}
	    }
	}
	else
	{
	    /* for each column j whos first row index is in row k */
	    for (j = Head [k] ; j != EMPTY ; j = Anext [j])
	    {
		/* merge column j of A into Flag */
		PRINT1 (("	---- A column "ID"\n", j)) ;
		p = Ap [j] ;
		pend = (apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
		PRINT1 (("  length "ID"  adding\n", pend-p)) ;
		for ( ; p < pend ; p++)
		{
#ifndef NDEBUG
		    ASSERT (Ai [p] >= k && Ai [p] < nrow) ;
		    if (Flag [Ai [p]] < mark) PRINT1 ((" row "ID"\n", Ai [p])) ;
#endif
		    Flag [Ai [p]] = mark ;
		}
	    }
	    /* clear the kth link list */
	    Head [k] = EMPTY ;
	}

	/* ------------------------------------------------------------------ */
	/* compute pruned pattern of kth column of L = union of children */
	/* ------------------------------------------------------------------ */

	/* for each column j of L whose parent is k */
	for (j = Link [k] ; j != EMPTY ; j = Link [j])
	{
	    /* merge column j of L into Flag */
	    PRINT1 (("	---- L column "ID"\n", k)) ;
	    ASSERT (j < k) ;
	    ASSERT (Lnz [j] > 0) ;
	    p = Lp [j] ;
	    pend = p + Lnz [j] ;
	    ASSERT (Li [p] == j && Li [p+1] == k) ;
	    p++ ;	    /* skip past the diagonal entry */
	    for ( ; p < pend ; p++)
	    {
		/* add to pattern */
		ASSERT (Li [p] >= k && Li [p] < nrow) ;
		Flag [Li [p]] = mark ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* prune the kth column of L */
	/* ------------------------------------------------------------------ */

	PRINT1 (("Final column of L:\n")) ;
	p = Lp [k] ;
	pend = p + Lnz [k] ;

	if (pack)
	{
	    /* shift column k upwards */
	    Lp [k] = pdest ;
	}
	else
	{
	    /* leave column k in place, just reduce Lnz [k] */
	    pdest = p ;
	}

	for ( ; p < pend ; p++)
	{
	    ASSERT (pdest < pend) ;
	    ASSERT (pdest <= p) ;
	    row = Li [p] ;
	    ASSERT (row >= k && row < nrow) ;
	    if (Flag [row] == mark)
	    {
		/* keep this entry */
		Li [pdest] = row ;
		if (xtype == CHOLMOD_REAL)
		{
		    Lx [pdest] = Lx [p] ;
		}
		else if (xtype == CHOLMOD_COMPLEX)
		{
		    Lx [2*pdest  ] = Lx [2*p  ] ;
		    Lx [2*pdest+1] = Lx [2*p+1] ;
		}
		else if (xtype == CHOLMOD_ZOMPLEX)
		{
		    Lx [pdest] = Lx [p] ;
		    Lz [pdest] = Lz [p] ;
		}
		pdest++ ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* prepare this column for its parent */
	/* ------------------------------------------------------------------ */

	Lnz [k] = pdest - Lp [k] ;

	PRINT1 ((" L("ID") length "ID"\n", k, Lnz [k])) ;
	ASSERT (Lnz [k] > 0) ;

	/* parent is the first entry in the column after the diagonal */
	parent = (Lnz [k] > 1) ? (Li [Lp [k] + 1]) : EMPTY ;

	PRINT1 (("parent ("ID") = "ID"\n", k, parent)) ;
	ASSERT ((parent > k && parent < nrow) || (parent == EMPTY)) ;

	if (parent != EMPTY)
	{
	    Link [k] = Link [parent] ;
	    Link [parent] = k ;
	}
    }

    /* done using Iwork for Link, Lnz (if needed), and Anext ] */

    /* ---------------------------------------------------------------------- */
    /* convert L to packed, if requested */
    /* ---------------------------------------------------------------------- */

    if (pack)
    {
	/* finalize Lp */
	Lp [nrow] = pdest ;
	/* Shrink L to be just large enough.  It cannot fail. */
	/* workspace: none */
	ASSERT ((size_t) (Lp [nrow]) <= L->nzmax) ;
	CHOLMOD(reallocate_factor) (Lp [nrow], L, Common) ;
	ASSERT (Common->status >= CHOLMOD_OK) ;
    }

    /* ---------------------------------------------------------------------- */
    /* clear workspace */
    /* ---------------------------------------------------------------------- */

    /* CHOLMOD(clear_flag) (Common) ; */
    CHOLMOD_CLEAR_FLAG (Common) ;

    DEBUG (CHOLMOD(dump_factor) (L, "ReSymbol final L (i, x):", Common)) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;
    return (TRUE) ;
}
#endif
