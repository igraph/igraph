/* ========================================================================== */
/* === Core/cholmod_factor ================================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2013,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Core utility routines for the cholmod_factor object:
 *
 * The data structure for an LL' or LDL' factorization is too complex to
 * describe in one sentence.  This object can hold the symbolic analysis alone,
 * or in combination with a "simplicial" (similar to a sparse matrix) or
 * "supernodal" form of the numerical factorization.  Only the routine to free
 * a factor is primary, since a factor object is created by the factorization
 * routine (cholmod_factorize).  It must be freed with cholmod_free_factor.
 *
 * Primary routine:
 * ----------------
 * cholmod_free_factor		free a factor
 *
 * Secondary routines:
 * -------------------
 * cholmod_allocate_factor	allocate a symbolic factor (LL' or LDL')
 * cholmod_reallocate_factor	change the # entries in a factor 
 * cholmod_change_factor	change the type of factor (e.g., LDL' to LL')
 * cholmod_pack_factor		pack the columns of a factor
 * cholmod_reallocate_column	resize a single column of a factor
 * cholmod_factor_to_sparse	create a sparse matrix copy of a factor
 * cholmod_copy_factor		create a copy of a factor
 *
 * Note that there is no cholmod_sparse_to_factor routine to create a factor
 * as a copy of a sparse matrix.  It could be done, after a fashion, but a
 * lower triangular sparse matrix would not necessarily have a chordal graph,
 * which would break the many CHOLMOD routines that rely on this property.
 *
 * The cholmod_factor_to_sparse routine is provided so that matrix operations
 * in the MatrixOps module may be applied to L.  Those operations operate on
 * cholmod_sparse objects, and they are not guaranteed to maintain the chordal
 * property of L.  Such a modified L cannot be safely converted back to a
 * cholmod_factor object.
 */

#include "cholmod_internal.h"
#include "cholmod_core.h"


/* ========================================================================== */
/* === cholmod_allocate_factor ============================================== */
/* ========================================================================== */

/* Allocate a simplicial symbolic factor, with L->Perm and L->ColCount allocated
 * and initialized to "empty" values (Perm [k] = k, and ColCount[k] = 1).
 * The integer and numerical parts of L are not allocated.  L->xtype is returned
 * as CHOLMOD_PATTERN and L->is_super are returned as FALSE.  L->is_ll is also
 * returned FALSE, but this may be modified when the matrix is factorized.
 *
 * This is sufficient (but far from ideal) for input to cholmod_factorize,
 * since the simplicial LL' or LDL' factorization (cholmod_rowfac) can
 * reallocate the columns of L as needed.  The primary purpose of this routine
 * is to allocate space for a symbolic factorization, for the "expert" user to
 * do his or her own symbolic analysis.  The typical user should use
 * cholmod_analyze instead of this routine.
 *
 * workspace: none
 */

cholmod_factor *CHOLMOD(allocate_factor)
(
    /* ---- input ---- */
    size_t n,		/* L is n-by-n */
    /* --------------- */
    cholmod_common *Common
)
{
    Int j ;
    Int *Perm, *ColCount ;
    cholmod_factor *L ;
    int ok = TRUE ;

    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;

    /* ensure the dimension does not cause integer overflow */
    (void) CHOLMOD(add_size_t) (n, 2, &ok) ;
    if (!ok || n > Int_max)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (NULL) ;
    }

    L = CHOLMOD(malloc) (sizeof (cholmod_factor), 1, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }
    L->n = n ;
    L->is_ll = FALSE ;
    L->is_super = FALSE ;
    L->is_monotonic = TRUE ;
    L->itype = ITYPE ;
    L->xtype = CHOLMOD_PATTERN ;
    L->dtype = DTYPE ;

    /* allocate the purely symbolic part of L */
    L->ordering = CHOLMOD_NATURAL ;
    L->Perm = CHOLMOD(malloc) (n, sizeof (Int), Common) ;
    L->IPerm = NULL ;       /* only created by cholmod_solve2 when needed */
    L->ColCount = CHOLMOD(malloc) (n, sizeof (Int), Common) ;

    /* simplicial part of L is empty */
    L->nzmax = 0 ;
    L->p = NULL ;
    L->i = NULL ;
    L->x = NULL ;
    L->z = NULL ;
    L->nz = NULL ;
    L->next = NULL ;
    L->prev = NULL ;

    /* supernodal part of L is also empty */
    L->nsuper = 0 ;
    L->ssize = 0 ;
    L->xsize = 0 ;
    L->maxesize = 0 ;
    L->maxcsize = 0 ;
    L->super = NULL ;
    L->pi = NULL ;
    L->px = NULL ;
    L->s = NULL ;

    /* L has not been factorized */
    L->minor = n ;

    if (Common->status < CHOLMOD_OK)
    {
	CHOLMOD(free_factor) (&L, Common) ;
	return (NULL) ;		/* out of memory */
    }

    /* initialize Perm and ColCount */
    Perm = L->Perm ;
    for (j = 0 ; j < ((Int) n) ; j++)
    {
	Perm [j] = j ;
    }
    ColCount = L->ColCount ;
    for (j = 0 ; j < ((Int) n) ; j++)
    {
	ColCount [j] = 1 ;
    }

    return (L) ;
}


/* ========================================================================== */
/* === cholmod_free_factor ================================================== */
/* ========================================================================== */

/* Free a factor object.
 *
 * workspace: none
 */

int CHOLMOD(free_factor)
(
    /* ---- in/out --- */
    cholmod_factor **LHandle,	/* factor to free, NULL on output */
    /* --------------- */
    cholmod_common *Common
)
{
    Int n, lnz, xs, ss, s ;
    cholmod_factor *L ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    if (LHandle == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }
    L = *LHandle ;
    if (L == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }

    n = L->n ;
    lnz = L->nzmax ;
    s = L->nsuper + 1 ;
    xs = (L->is_super) ? ((Int) (L->xsize)) : (lnz) ;
    ss = L->ssize ;

    /* symbolic part of L */
    CHOLMOD(free) (n,   sizeof (Int), L->Perm,     Common) ;
    CHOLMOD(free) (n,   sizeof (Int), L->IPerm,    Common) ;
    CHOLMOD(free) (n,   sizeof (Int), L->ColCount, Common) ;

    /* simplicial form of L */
    CHOLMOD(free) (n+1, sizeof (Int), L->p,        Common) ;
    CHOLMOD(free) (lnz, sizeof (Int), L->i,        Common) ;
    CHOLMOD(free) (n,   sizeof (Int), L->nz,       Common) ;
    CHOLMOD(free) (n+2, sizeof (Int), L->next,     Common) ;
    CHOLMOD(free) (n+2, sizeof (Int), L->prev,     Common) ;

    /* supernodal form of L */
    CHOLMOD(free) (s,   sizeof (Int), L->pi,       Common) ;
    CHOLMOD(free) (s,   sizeof (Int), L->px,       Common) ;
    CHOLMOD(free) (s,   sizeof (Int), L->super,    Common) ;
    CHOLMOD(free) (ss,  sizeof (Int), L->s,        Common) ;

    /* numerical values for both simplicial and supernodal L */
    if (L->xtype == CHOLMOD_REAL)
    {
	CHOLMOD(free) (xs, sizeof (double), L->x, Common) ;
    }
    else if (L->xtype == CHOLMOD_COMPLEX)
    {
	CHOLMOD(free) (xs, 2*sizeof (double), L->x, Common) ;
    }
    else if (L->xtype == CHOLMOD_ZOMPLEX)
    {
	CHOLMOD(free) (xs, sizeof (double), L->x, Common) ;
	CHOLMOD(free) (xs, sizeof (double), L->z, Common) ;
    }

    *LHandle = CHOLMOD(free) (1, sizeof (cholmod_factor), (*LHandle), Common) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_reallocate_factor ============================================ */
/* ========================================================================== */

/* Change the size of L->i and L->x, or allocate them if their current size
 * is zero.  L must be simplicial.
 *
 * workspace: none
 */

int CHOLMOD(reallocate_factor)
(
    /* ---- input ---- */
    size_t nznew,	/* new # of entries in L */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    /* --------------- */
    cholmod_common *Common
)
{
    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    PRINT1 (("realloc factor: xtype %d\n", L->xtype)) ;
    if (L->is_super)
    {
	/* L must be simplicial, and not symbolic */
	ERROR (CHOLMOD_INVALID, "L invalid") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;
    PRINT1 (("realloc factor %g to %g\n", (double) L->nzmax, (double) nznew)) ;

    /* ---------------------------------------------------------------------- */
    /* resize (or allocate) the L->i and L->x components of the factor */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(realloc_multiple) (nznew, 1, L->xtype, &(L->i), NULL,
	    &(L->x), &(L->z), &(L->nzmax), Common) ;
    return (Common->status == CHOLMOD_OK) ;
}


/* ========================================================================== */
/* === cholmod_reallocate_column =========================================== */
/* ========================================================================== */

/* Column j needs more space, reallocate it at the end of L->i and L->x.
 * If the reallocation fails, the factor is converted to a simplicial
 * symbolic factor (no pattern, just L->Perm and L->ColCount).
 *
 * workspace: none
 */

int CHOLMOD(reallocate_column)
(
    /* ---- input ---- */
    size_t j,		/* the column to reallocate */
    size_t need,	/* required size of column j */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    /* --------------- */
    cholmod_common *Common
)
{
    double xneed ;
    double *Lx, *Lz ;
    Int *Lp, *Lprev, *Lnext, *Li, *Lnz ;
    Int n, pold, pnew, len, k, tail ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    if (L->is_super)
    {
	ERROR (CHOLMOD_INVALID, "L must be simplicial") ;
	return (FALSE) ;
    }
    n = L->n ;
    if (j >= L->n || need == 0)
    {
	ERROR (CHOLMOD_INVALID, "j invalid") ;
	return (FALSE) ;	    /* j out of range */
    }
    Common->status = CHOLMOD_OK ;

    DEBUG (CHOLMOD(dump_factor) (L, "start colrealloc", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* increase the size of L if needed */
    /* ---------------------------------------------------------------------- */

    /* head = n+1 ; */
    tail = n ;
    Lp = L->p ;
    Lnz = L->nz ;
    Lprev = L->prev ;
    Lnext = L->next ;

    ASSERT (Lnz != NULL) ;
    ASSERT (Lnext != NULL && Lprev != NULL) ;
    PRINT1 (("col %g need %g\n", (double) j, (double) need)) ;

    /* column j cannot have more than n-j entries if all entries are present */
    need = MIN (need, n-j) ;

    /* compute need in double to avoid integer overflow */
    if (Common->grow1 >= 1.0)
    {
	xneed = (double) need ;
	xneed = Common->grow1 * xneed + Common->grow2 ;
	xneed = MIN (xneed, n-j) ;
	need = (Int) xneed ;
    }
    PRINT1 (("really new need %g current %g\n", (double) need,
	    (double) (Lp [Lnext [j]] - Lp [j]))) ;
    ASSERT (need >= 1 && need <= n-j) ;

    if (Lp [Lnext [j]] - Lp [j] >= (Int) need)
    {
	/* no need to reallocate the column, it's already big enough */
	PRINT1 (("colrealloc: quick return %g %g\n",
	    (double) (Lp [Lnext [j]] - Lp [j]), (double) need)) ;
	return (TRUE) ;

    }

    if (Lp [tail] + need > L->nzmax)
    {
	/* use double to avoid integer overflow */
	xneed = (double) need ;
	if (Common->grow0 < 1.2)	    /* fl. pt. compare, false if NaN */
	{
	    /* if grow0 is less than 1.2 or NaN, don't use it */
	    xneed = 1.2 * (((double) L->nzmax) + xneed + 1) ;
	}
	else
	{
	    xneed = Common->grow0 * (((double) L->nzmax) + xneed + 1) ;
	}
	if (xneed > Size_max ||
		!CHOLMOD(reallocate_factor) ((Int) xneed, L, Common))
	{
	    /* out of memory, convert to simplicial symbolic */
	    CHOLMOD(change_factor) (CHOLMOD_PATTERN, L->is_ll, FALSE, TRUE,
		    TRUE, L, Common) ;
	    ERROR (CHOLMOD_OUT_OF_MEMORY, "out of memory; L now symbolic") ;
	    return (FALSE) ;	    /* out of memory */
	}
	PRINT1 (("\n=== GROW L from %g to %g\n",
		    (double) L->nzmax, (double) xneed)) ;
	/* pack all columns so that each column has at most grow2 free space */
	CHOLMOD(pack_factor) (L, Common) ;
	ASSERT (Common->status == CHOLMOD_OK) ;
	Common->nrealloc_factor++ ;
    }

    /* ---------------------------------------------------------------------- */
    /* reallocate the column */
    /* ---------------------------------------------------------------------- */

    Common->nrealloc_col++ ;

    Li = L->i ;
    Lx = L->x ;
    Lz = L->z ;

    /* remove j from its current position in the list */
    Lnext [Lprev [j]] = Lnext [j] ;
    Lprev [Lnext [j]] = Lprev [j] ;

    /* place j at the end of the list */
    Lnext [Lprev [tail]] = j ;
    Lprev [j] = Lprev [tail] ;
    Lnext [j] = n ;
    Lprev [tail] = j ;

    /* L is no longer monotonic; columns are out-of-order */
    L->is_monotonic = FALSE ;

    /* allocate space for column j */
    pold = Lp [j] ;
    pnew = Lp [tail] ;
    Lp [j] = pnew  ;
    Lp [tail] += need ;

    /* copy column j to the new space */
    len = Lnz [j] ;
    for (k = 0 ; k < len ; k++)
    {
	Li [pnew + k] = Li [pold + k] ;
    }

    if (L->xtype == CHOLMOD_REAL)
    {
	for (k = 0 ; k < len ; k++)
	{
	    Lx [pnew + k] = Lx [pold + k] ;
	}
    }
    else if (L->xtype == CHOLMOD_COMPLEX)
    {
	for (k = 0 ; k < len ; k++)
	{
	    Lx [2*(pnew + k)  ] = Lx [2*(pold + k)  ] ;
	    Lx [2*(pnew + k)+1] = Lx [2*(pold + k)+1] ;
	}
    }
    else if (L->xtype == CHOLMOD_ZOMPLEX)
    {
	for (k = 0 ; k < len ; k++)
	{
	    Lx [pnew + k] = Lx [pold + k] ;
	    Lz [pnew + k] = Lz [pold + k] ;
	}
    }

    DEBUG (CHOLMOD(dump_factor) (L, "colrealloc done", Common)) ;

    /* successful reallocation of column j of L */
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_pack_factor ================================================== */
/* ========================================================================== */

/* Pack the columns of a simplicial LDL' or LL' factor.  This can be followed
 * by a call to cholmod_reallocate_factor to reduce the size of L to the exact
 * size required by the factor, if desired.  Alternatively, you can leave the
 * size of L->i and L->x the same, to allow space for future updates/rowadds.
 *
 * Each column is reduced in size so that it has at most Common->grow2 free
 * space at the end of the column.
 *
 * Does nothing and returns silently if given any other type of factor.
 *
 * Does NOT force the columns of L to be monotonic.  It thus differs from
 * cholmod_change_factor (xtype, -, FALSE, TRUE, TRUE, L, Common), which
 * packs the columns and ensures that they appear in monotonic order.
 */

int CHOLMOD(pack_factor)
(
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Lx, *Lz ;
    Int *Lp, *Li, *Lnz, *Lnext ;
    Int pnew, j, k, pold, len, n, head, tail, grow2 ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    Common->status = CHOLMOD_OK ;
    DEBUG (CHOLMOD(dump_factor) (L, "start pack", Common)) ;
    PRINT1 (("PACK factor %d\n", L->is_super)) ;

    if (L->xtype == CHOLMOD_PATTERN || L->is_super)
    {
	/* nothing to do unless L is simplicial numeric */
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* pack */
    /* ---------------------------------------------------------------------- */

    grow2 = Common->grow2 ;
    PRINT1 (("\nPACK grow2 "ID"\n", grow2)) ;

    pnew = 0 ;
    n = L->n ;
    Lp = L->p ;
    Li = L->i ;
    Lx = L->x ;
    Lz = L->z ;
    Lnz = L->nz ;
    Lnext = L->next ;

    head = n+1 ;
    tail = n ;

    for (j = Lnext [head] ; j != tail ; j = Lnext [j])
    {
	/* pack column j */
	pold = Lp [j] ;
	len = Lnz [j] ;
	ASSERT (len > 0) ;
	PRINT2 (("col "ID" pnew "ID" pold "ID"\n", j, pnew, pold)) ;
	if (pnew < pold)
	{
	    PRINT2 (("    pack this column\n")) ;

	    for (k = 0 ; k < len ; k++)
	    {
		Li [pnew + k] = Li [pold + k] ;
	    }

	    if (L->xtype == CHOLMOD_REAL)
	    {
		for (k = 0 ; k < len ; k++)
		{
		    Lx [pnew + k] = Lx [pold + k] ;
		}
	    }
	    else if (L->xtype == CHOLMOD_COMPLEX)
	    {
		for (k = 0 ; k < len ; k++)
		{
		    Lx [2*(pnew + k)  ] = Lx [2*(pold + k)  ] ;
		    Lx [2*(pnew + k)+1] = Lx [2*(pold + k)+1] ;
		}
	    }
	    else if (L->xtype == CHOLMOD_ZOMPLEX)
	    {
		for (k = 0 ; k < len ; k++)
		{
		    Lx [pnew + k] = Lx [pold + k] ;
		    Lz [pnew + k] = Lz [pold + k] ;
		}
	    }

	    Lp [j] = pnew ;
	}
	len = MIN (len + grow2, n - j) ;
	pnew = MIN (Lp [j] + len, Lp [Lnext [j]]) ;
    }
    PRINT2 (("final pnew = "ID"\n", pnew)) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_factor_to_sparse ============================================= */
/* ========================================================================== */

/* Constructs a column-oriented sparse matrix containing the pattern and values
 * of a simplicial or supernodal numerical factor, and then converts the factor
 * into a simplicial symbolic factor.  If L is already packed, monotonic,
 * and simplicial (which is the case when cholmod_factorize uses the simplicial
 * Cholesky factorization algorithm) then this routine requires only O(1)
 * memory and takes O(1) time.
 *
 * Only operates on numeric factors (real, complex, or zomplex).  Does not
 * change the numeric L->xtype (the resulting sparse matrix has the same xtype
 * as L).  If this routine fails, L is left unmodified.
 */

cholmod_sparse *CHOLMOD(factor_to_sparse)
(
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to copy, converted to symbolic on output */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_sparse *Lsparse ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (L, NULL) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, NULL) ;
    Common->status = CHOLMOD_OK ;
    DEBUG (CHOLMOD(dump_factor) (L, "start convert to matrix", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* convert to packed, monotonic, simplicial, numeric */
    /* ---------------------------------------------------------------------- */

    /* leave as LL or LDL' */
    if (!CHOLMOD(change_factor) (L->xtype, L->is_ll, FALSE, TRUE, TRUE, L,
		Common))
    {
	ERROR (CHOLMOD_INVALID, "cannot convert L") ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* create Lsparse */
    /* ---------------------------------------------------------------------- */

    /* allocate the header for Lsparse, the sparse matrix version of L */
    Lsparse = CHOLMOD(malloc) (sizeof (cholmod_sparse), 1, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;		/* out of memory */
    }

    /* transfer the contents from L to Lsparse */
    Lsparse->nrow = L->n ;
    Lsparse->ncol = L->n ;
    Lsparse->p = L->p ;
    Lsparse->i = L->i ;
    Lsparse->x = L->x ;
    Lsparse->z = L->z ;
    Lsparse->nz = NULL ;
    Lsparse->stype = 0 ;
    Lsparse->itype = L->itype ;
    Lsparse->xtype = L->xtype ;
    Lsparse->dtype = L->dtype ;
    Lsparse->sorted = TRUE ;
    Lsparse->packed = TRUE ;
    Lsparse->nzmax = L->nzmax ;
    ASSERT (CHOLMOD(dump_sparse) (Lsparse, "Lsparse", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* convert L to symbolic, but do not free contents transfered to Lsparse */
    /* ---------------------------------------------------------------------- */

    L->p = NULL ;
    L->i = NULL ;
    L->x = NULL ;
    L->z = NULL ;
    L->xtype = CHOLMOD_PATTERN ;
    CHOLMOD(change_factor) (CHOLMOD_PATTERN, FALSE, FALSE, TRUE, TRUE, L,
	    Common) ;

    return (Lsparse) ;
}


/* ========================================================================== */
/* === cholmod_copy_factor ================================================== */
/* ========================================================================== */

/* Create an exact copy of a factor, with one exception:
 *
 * Entries in unused space are not copied (they might not be initialized,
 *	and copying them would cause program checkers such as purify and
 *	valgrind to complain).
 *
 * Note that a supernodal L cannot be zomplex.
 */

cholmod_factor *CHOLMOD(copy_factor)
(
    /* ---- input ---- */
    cholmod_factor *L,	/* factor to copy */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_factor *L2 ;
    double *Lx, *L2x, *Lz, *L2z ;
    Int *Perm, *ColCount, *Lp, *Li, *Lnz, *Lnext, *Lprev, *Lsuper, *Lpi, *Lpx,
	*Ls, *Perm2, *ColCount2, *L2p, *L2i, *L2nz, *L2next, *L2prev, *L2super,
	*L2pi, *L2px, *L2s ;
    Int n, j, p, pend, s, xsize, ssize, nsuper ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (L, NULL) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, NULL) ;
    Common->status = CHOLMOD_OK ;
    DEBUG (CHOLMOD(dump_factor) (L, "start copy", Common)) ;

    n = L->n ;

    /* ---------------------------------------------------------------------- */
    /* allocate a simplicial symbolic factor  */
    /* ---------------------------------------------------------------------- */

    /* allocates L2->Perm and L2->ColCount */
    L2 = CHOLMOD(allocate_factor) (n, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }
    ASSERT (L2->xtype == CHOLMOD_PATTERN && !(L2->is_super)) ;

    Perm = L->Perm ;
    ColCount = L->ColCount ;
    Perm2 = L2->Perm ;
    ColCount2 = L2->ColCount ;
    L2->ordering = L->ordering ;

    for (j = 0 ; j < n ; j++)
    {
	Perm2 [j] = Perm [j] ;
    }
    for (j = 0 ; j < n ; j++)
    {
	ColCount2 [j] = ColCount [j] ;
    }
    L2->is_ll = L->is_ll ;

    /* ---------------------------------------------------------------------- */
    /* copy the rest of the factor */
    /* ---------------------------------------------------------------------- */

    if (L->xtype != CHOLMOD_PATTERN && !(L->super))
    {

	/* ------------------------------------------------------------------ */
	/* allocate a simplicial numeric factor */
	/* ------------------------------------------------------------------ */

	/* allocate L2->p, L2->nz, L2->prev, L2->next, L2->i, and L2->x.
	 * packed = -1 so that cholmod_change_factor allocates space of
	 * size L2->nzmax */
	L2->nzmax = L->nzmax ;
	if (!CHOLMOD(change_factor) (L->xtype, L->is_ll, FALSE, -1, TRUE,
		    L2, Common))
	{
	    CHOLMOD(free_factor) (&L2, Common) ;
	    return (NULL) ;	/* out of memory */
	}
	ASSERT (MAX (1, L->nzmax) == L2->nzmax) ;

	/* ------------------------------------------------------------------ */
	/* copy the contents of a simplicial numeric factor */
	/* ------------------------------------------------------------------ */

	Lp = L->p ;
	Li = L->i ;
	Lx = L->x ;
	Lz = L->z ;
	Lnz = L->nz ;
	Lnext = L->next ;
	Lprev = L->prev ;

	L2p = L2->p ;
	L2i = L2->i ;
	L2x = L2->x ;
	L2z = L2->z ;
	L2nz = L2->nz ;
	L2next = L2->next ;
	L2prev = L2->prev ;
	L2->xtype = L->xtype ;
	L2->dtype = L->dtype ;

	for (j = 0 ; j <= n ; j++)
	{
	    L2p [j] = Lp [j] ;
	}

	for (j = 0 ; j < n+2 ; j++)
	{
	    L2prev [j] = Lprev [j] ;
	}

	for (j = 0 ; j < n+2 ; j++)
	{
	    L2next [j] = Lnext [j] ;
	}

	for (j = 0 ; j < n ; j++)
	{
	    L2nz [j] = Lnz [j] ;
	}

	for (j = 0 ; j < n ; j++)
	{
	    p = Lp [j] ;
	    pend = p + Lnz [j] ;
	    for ( ; p < pend ; p++)
	    {
		L2i [p] = Li [p] ;
	    }
	    p = Lp [j] ;

	    if (L->xtype == CHOLMOD_REAL)
	    {
		for ( ; p < pend ; p++)
		{
		    L2x [p] = Lx [p] ;
		}
	    }
	    else if (L->xtype == CHOLMOD_COMPLEX)
	    {
		for ( ; p < pend ; p++)
		{
		    L2x [2*p  ] = Lx [2*p  ] ;
		    L2x [2*p+1] = Lx [2*p+1] ;
		}
	    }
	    else if (L->xtype == CHOLMOD_ZOMPLEX)
	    {
		for ( ; p < pend ; p++)
		{
		    L2x [p] = Lx [p] ;
		    L2z [p] = Lz [p] ;
		}
	    }

	}

    }
    else if (L->is_super)
    {

	/* ------------------------------------------------------------------ */
	/* copy a supernodal factor */
	/* ------------------------------------------------------------------ */

	xsize = L->xsize ;
	ssize = L->ssize ;
	nsuper = L->nsuper ;

	L2->xsize = xsize ;
	L2->ssize = ssize ;
	L2->nsuper = nsuper ;

	/* allocate L2->super, L2->pi, L2->px, and L2->s.  Allocate L2->x if
	 * L is numeric */
	if (!CHOLMOD(change_factor) (L->xtype, TRUE, TRUE, TRUE, TRUE, L2,
		    Common))
	{
	    CHOLMOD(free_factor) (&L2, Common) ;
	    return (NULL) ;	/* out of memory */
	}

	ASSERT (L2->s != NULL) ;

	/* ------------------------------------------------------------------ */
	/* copy the contents of a supernodal factor */
	/* ------------------------------------------------------------------ */

	Lsuper = L->super ;
	Lpi = L->pi ;
	Lpx = L->px ;
	Ls = L->s ;
	Lx = L->x ;

	L2super = L2->super ;
	L2pi = L2->pi ;
	L2px = L2->px ;
	L2s = L2->s ;
	L2x = L2->x ;

	L2->maxcsize = L->maxcsize ;
	L2->maxesize = L->maxesize ;

	for (s = 0 ; s <= nsuper ; s++)
	{
	    L2super [s] = Lsuper [s] ;
	}
	for (s = 0 ; s <= nsuper ; s++)
	{
	    L2pi [s] = Lpi [s] ;
	}
	for (s = 0 ; s <= nsuper ; s++)
	{
	    L2px [s] = Lpx [s] ;
	}

	L2s [0] = 0 ;
	for (p = 0 ; p < ssize ; p++)
	{
	    L2s [p] = Ls [p] ;
	}

	if (L->xtype == CHOLMOD_REAL)
	{
	    for (p = 0 ; p < xsize ; p++)
	    {
		L2x [p] = Lx [p] ;
	    }
	}
	else if (L->xtype == CHOLMOD_COMPLEX)
	{
	    for (p = 0 ; p < 2*xsize ; p++)
	    {
		L2x [p] = Lx [p] ;
	    }
	}
    }

    L2->minor = L->minor ;
    L2->is_monotonic = L->is_monotonic ;

    DEBUG (CHOLMOD(dump_factor) (L2, "L2 got copied", Common)) ;
    ASSERT (L2->xtype == L->xtype && L2->is_super == L->is_super) ;
    return (L2) ;
}
