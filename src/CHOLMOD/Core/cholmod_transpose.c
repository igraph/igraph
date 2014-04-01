/* ========================================================================== */
/* === Core/cholmod_transpose =============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2006,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Core utility routines for the cholmod_sparse object to
 * compute the transpose or permuted transpose of a matrix:
 *
 * Primary routines:
 * -----------------
 * cholmod_transpose		transpose sparse matrix
 * cholmod_ptranspose		transpose and permute sparse matrix
 * cholmod_sort			sort row indices in each column of sparse matrix
 *
 * Secondary routines:
 * -------------------
 * cholmod_transpose_unsym	transpose unsymmetric sparse matrix
 * cholmod_transpose_sym	transpose symmetric sparse matrix
 *
 * All xtypes (pattern, real, complex, and zomplex) are supported.
 *
 * ---------------------------------------
 * Unsymmetric case: A->stype is zero.
 * ---------------------------------------
 *
 * Computes F = A', F = A (:,f)' or F = A (p,f)', except that the indexing by
 * f does not work the same as the MATLAB notation (see below).  A->stype
 * is zero, which denotes that both the upper and lower triangular parts of
 * A are present (and used).  A may in fact be symmetric in pattern and/or
 * value; A->stype just denotes which part of A are stored.  A may be
 * rectangular.
 *
 * p is a permutation of 0:m-1, and f is a subset of 0:n-1, where A is m-by-n.
 * There can be no duplicate entries in p or f.
 *
 * The set f is held in fset and fsize.
 *	fset = NULL means ":" in MATLAB. fsize is ignored.
 *	fset != NULL means f = fset [0..fsize-1].
 *	fset != NULL and fsize = 0 means f is the empty set.
 *
 * Columns not in the set f are considered to be zero.  That is,
 * if A is 5-by-10 then F = A (:,[3 4])' is not 2-by-5, but 10-by-5, and rows
 * 3 and 4 of F are equal to columns 3 and 4 of A (the other rows of F are
 * zero).  More precisely, in MATLAB notation:
 *
 *	[m n] = size (A) ;
 *	F = A ;
 *	notf = ones (1,n) ;
 *	notf (f) = 0 ;
 *	F (:, find (notf)) = 0
 *	F = F'
 *
 * If you want the MATLAB equivalent F=A(p,f) operation, use cholmod_submatrix
 * instead (which does not compute the transpose).
 *
 * F->nzmax must be large enough to hold the matrix F.  It is not modified.
 * If F->nz is present then F->nz [j] = # of entries in column j of F.
 *
 * A can be sorted or unsorted, with packed or unpacked columns.
 *
 * If f is present and not sorted in ascending order, then F is unsorted
 * (that is, it may contain columns whose row indices do not appear in
 * ascending order).  Otherwise, F is sorted (the row indices in each
 * column of F appear in strictly ascending order).
 *
 * F is returned in packed or unpacked form, depending on F->packed on input.
 * If F->packed is false, then F is returned in unpacked form (F->nz must be
 * present).  Each row i of F is large enough to hold all the entries in row i
 * of A, even if f is provided.  That is, F->i and
 * F->x [F->p [i] .. F->p [i] + F->nz [i] - 1] contain all entries in A (i,f),
 * but F->p [i+1] - F->p [i] is equal to the number of nonzeros in A (i,:),
 * not just A (i,f).
 *
 * The cholmod_transpose_unsym routine is the only operation in CHOLMOD that
 * can produce an unpacked matrix.
 *
 * ---------------------------------------
 * Symmetric case: A->stype is nonzero.
 * ---------------------------------------
 *
 * Computes F = A' or F = A(p,p)', the transpose or permuted transpose, where
 * A->stype is nonzero.
 *
 * If A->stype > 0, then A is a symmetric matrix where just the upper part
 * of the matrix is stored.  Entries in the lower triangular part may be
 * present, but are ignored.  A must be square.  If F=A', then F is returned
 * sorted; otherwise F is unsorted for the F=A(p,p)' case.
 *
 * There can be no duplicate entries in p.
 * The fset and fsize parameters are not used.
 *
 * Three kinds of transposes are available, depending on the "values" parameter:
 * 0: do not transpose the numerical values; create a CHOLMOD_PATTERN matrix
 * 1: array transpose
 * 2: complex conjugate transpose (same as 2 if input is real or pattern)
 *
 * -----------------------------------------------------------------------------
 *
 * For cholmod_transpose_unsym and cholmod_transpose_sym, the output matrix
 * F must already be pre-allocated by the caller, with the correct dimensions.
 * If F is not valid or has the wrong dimensions, it is not modified.
 * Otherwise, if F is too small, the transpose is not computed; the contents
 * of F->p contain the column pointers of the resulting matrix, where
 * F->p [F->ncol] > F->nzmax.  In this case, the remaining contents of F are
 * not modified.  F can still be properly free'd with cholmod_free_sparse.
 */

#include "cholmod_internal.h"
#include "cholmod_core.h"


/* ========================================================================== */
/* === TEMPLATE ============================================================= */
/* ========================================================================== */

#define PATTERN
#include "t_cholmod_transpose.c"
#define REAL
#include "t_cholmod_transpose.c"
#define COMPLEX
#include "t_cholmod_transpose.c"
#define COMPLEX
#define NCONJUGATE
#include "t_cholmod_transpose.c"
#define ZOMPLEX
#include "t_cholmod_transpose.c"
#define ZOMPLEX
#define NCONJUGATE
#include "t_cholmod_transpose.c"


/* ========================================================================== */
/* === cholmod_transpose_unsym ============================================== */
/* ========================================================================== */

/* Compute F = A', A (:,f)', or A (p,f)', where A is unsymmetric and F is
 * already allocated.  See cholmod_transpose for a simpler routine.
 *
 * workspace:
 * Iwork (MAX (nrow,ncol)) if fset is present
 * Iwork (nrow) if fset is NULL
 *
 * The xtype of A and F must match, unless values is zero or F->xtype is
 * CHOLMOD_PATTERN (in which case only the pattern of A is transpose into F).
 */

int CHOLMOD(transpose_unsym)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to transpose */
    int values,		/* 2: complex conj. transpose, 1: array transpose,
			   0: do not transpose the numerical values */
    Int *Perm,		/* size nrow, if present (can be NULL) */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* ---- output --- */
    cholmod_sparse *F,	/* F = A', A(:,f)', or A(p,f)' */
    /* --------------- */
    cholmod_common *Common
)
{
    Int *Fp, *Fnz, *Ap, *Ai, *Anz, *Wi ;
    Int nrow, ncol, permute, use_fset, Apacked, Fpacked, p, pend,
	i, j, k, Fsorted, nf, jj, jlast ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (F, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (F, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    if (A->nrow != F->ncol || A->ncol != F->nrow)
    {
	ERROR (CHOLMOD_INVALID, "F has the wrong dimensions") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nf = fsize ;
    use_fset = (fset != NULL) ;
    nrow = A->nrow ;
    ncol = A->ncol ;

    Ap = A->p ;		/* size A->ncol+1, column pointers of A */
    Ai = A->i ;		/* size nz = Ap [A->ncol], row indices of A */
    Anz = A->nz ;
    Apacked = A->packed ;
    ASSERT (IMPLIES (!Apacked, Anz != NULL)) ;

    permute = (Perm != NULL) ;

    Fp = F->p ;		/* size A->nrow+1, row pointers of F */
    Fnz = F->nz ;
    Fpacked = F->packed ;
    ASSERT (IMPLIES (!Fpacked, Fnz != NULL)) ;

    nf = (use_fset) ? nf : ncol ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* s = nrow + ((fset != NULL) ? ncol : 0) */
    s = CHOLMOD(add_size_t) (nrow, ((fset != NULL) ? ncol : 0), &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    CHOLMOD(allocate_work) (0, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;	/* out of memory */
    }

    Wi = Common->Iwork ;	/* size nrow (i/l/l) */

    /* ---------------------------------------------------------------------- */
    /* check Perm and fset */
    /* ---------------------------------------------------------------------- */

    if (permute)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    Wi [i] = 1 ;
	}
	for (k = 0 ; k < nrow ; k++)
	{
	    i = Perm [k] ;
	    if (i < 0 || i > nrow || Wi [i] == 0)
	    {
		ERROR (CHOLMOD_INVALID, "invalid permutation") ;
		return (FALSE) ;
	    }
	    Wi [i] = 0 ;
	}
    }

    if (use_fset)
    {
	for (j = 0 ; j < ncol ; j++)
	{
	    Wi [j] = 1 ;
	}
	for (k = 0 ; k < nf ; k++)
	{
	    j = fset [k] ;
	    if (j < 0 || j > ncol || Wi [j] == 0)
	    {
		ERROR (CHOLMOD_INVALID, "invalid fset") ;
		return (FALSE) ;
	    }
	    Wi [j] = 0 ;
	}
    }

    /* Perm and fset are now valid */
    ASSERT (CHOLMOD(dump_perm) (Perm, nrow, nrow, "Perm", Common)) ;
    ASSERT (CHOLMOD(dump_perm) (fset, nf, ncol, "fset", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* count the entries in each row of A or A(:,f) */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < nrow ; i++)
    {
	Wi [i] = 0 ;
    }

    jlast = EMPTY ;
    Fsorted = TRUE ;

    if (use_fset)
    {
	/* count entries in each row of A(:,f) */
	for (jj = 0 ; jj < nf ; jj++)
	{
	    j = fset [jj] ;
	    if (j <= jlast)
	    {
		Fsorted = FALSE ;
	    }
	    p = Ap [j] ;
	    pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		Wi [Ai [p]]++ ;
	    }
	    jlast = j ;
	}

	/* save the nz counts if F is unpacked, and recount all of A */
	if (!Fpacked)
	{
	    if (permute)
	    {
		for (i = 0 ; i < nrow ; i++)
		{
		    Fnz [i] = Wi [Perm [i]] ;
		}
	    }
	    else
	    {
		for (i = 0 ; i < nrow ; i++)
		{
		    Fnz [i] = Wi [i] ;
		}
	    }
	    for (i = 0 ; i < nrow ; i++)
	    {
		Wi [i] = 0 ;
	    }

	    /* count entries in each row of A */
	    for (j = 0 ; j < ncol ; j++)
	    {
		p = Ap [j] ;
		pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    Wi [Ai [p]]++ ;
		}
	    }
	}

    }
    else
    {

	/* count entries in each row of A */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		Wi [Ai [p]]++ ;
	    }
	}

	/* save the nz counts if F is unpacked */
	if (!Fpacked)
	{
	    if (permute)
	    {
		for (i = 0 ; i < nrow ; i++)
		{
		    Fnz [i] = Wi [Perm [i]] ;
		}
	    }
	    else
	    {
		for (i = 0 ; i < nrow ; i++)
		{
		    Fnz [i] = Wi [i] ;
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* compute the row pointers */
    /* ---------------------------------------------------------------------- */

    p = 0 ;
    if (permute)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    Fp [i] = p ;
	    p += Wi [Perm [i]] ;
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    Wi [Perm [i]] = Fp [i] ;
	}
    }
    else
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    Fp [i] = p ;
	    p += Wi [i] ;
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    Wi [i] = Fp [i] ;
	}
    }
    Fp [nrow] = p ;

    if (p > (Int) (F->nzmax))
    {
	ERROR (CHOLMOD_INVALID, "F is too small") ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* transpose matrix, using template routine */
    /* ---------------------------------------------------------------------- */

    ok = FALSE ;
    if (values == 0 || F->xtype == CHOLMOD_PATTERN)
    {
	ok = p_cholmod_transpose_unsym (A, Perm, fset, nf, F, Common) ;
    }
    else if (F->xtype == CHOLMOD_REAL)
    {
	ok = r_cholmod_transpose_unsym (A, Perm, fset, nf, F, Common) ;
    }
    else if (F->xtype == CHOLMOD_COMPLEX)
    {
	if (values == 1)
	{
	    /* array transpose */
	    ok = ct_cholmod_transpose_unsym (A, Perm, fset, nf, F, Common) ;
	}
	else
	{
	    /* complex conjugate transpose */
	    ok = c_cholmod_transpose_unsym (A, Perm, fset, nf, F, Common) ;
	}
    }
    else if (F->xtype == CHOLMOD_ZOMPLEX)
    {
	if (values == 1)
	{
	    /* array transpose */
	    ok = zt_cholmod_transpose_unsym (A, Perm, fset, nf, F, Common) ;
	}
	else
	{
	    /* complex conjugate transpose */
	    ok = z_cholmod_transpose_unsym (A, Perm, fset, nf, F, Common) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* finalize result F */
    /* ---------------------------------------------------------------------- */

    if (ok)
    {
	F->sorted = Fsorted ;
    }
    ASSERT (CHOLMOD(dump_sparse) (F, "output F unsym", Common) >= 0) ;
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_transpose_sym ================================================ */
/* ========================================================================== */

/* Compute F = A' or A (p,p)', where A is symmetric and F is already allocated.
 * See cholmod_transpose for a simpler routine.
 *
 * workspace:  Iwork (nrow) if Perm NULL, Iwork (2*nrow) if Perm non-NULL.
 */

int CHOLMOD(transpose_sym)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to transpose */
    int values,		/* 2: complex conj. transpose, 1: array transpose,
			   0: do not transpose the numerical values */
    Int *Perm,		/* size nrow, if present (can be NULL) */
    /* ---- output --- */
    cholmod_sparse *F,	/* F = A' or A(p,p)' */
    /* --------------- */
    cholmod_common *Common
)
{
    Int *Ap, *Anz, *Ai, *Fp, *Wi, *Pinv, *Iwork ;
    Int p, pend, packed, upper, permute, jold, n, i, j, k, iold ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (F, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (F, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    if (A->nrow != A->ncol || A->stype == 0)
    {
	/* this routine handles square symmetric matrices only */
	ERROR (CHOLMOD_INVALID, "matrix must be symmetric") ;
	return (FALSE) ;
    }
    if (A->nrow != F->ncol || A->ncol != F->nrow)
    {
	ERROR (CHOLMOD_INVALID, "F has the wrong dimensions") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    permute = (Perm != NULL) ;
    n = A->nrow ;
    Ap = A->p ;		/* size A->ncol+1, column pointers of A */
    Ai = A->i ;		/* size nz = Ap [A->ncol], row indices of A */
    Anz = A->nz ;
    packed = A->packed ;
    ASSERT (IMPLIES (!packed, Anz != NULL)) ;
    upper = (A->stype > 0) ;

    Fp = F->p ;		/* size A->nrow+1, row pointers of F */

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* s = (Perm != NULL) ? 2*n : n */
    s = CHOLMOD(add_size_t) (n, ((Perm != NULL) ? n : 0), &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    CHOLMOD(allocate_work) (0, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;	/* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Wi   = Iwork ;	    /* size n (i/l/l) */
    Pinv = Iwork + n ;	    /* size n (i/i/l) , unused if Perm NULL */

    /* ---------------------------------------------------------------------- */
    /* check Perm and construct inverse permutation */
    /* ---------------------------------------------------------------------- */

    if (permute)
    {
	for (i = 0 ; i < n ; i++)
	{
	    Pinv [i] = EMPTY ;
	}
	for (k = 0 ; k < n ; k++)
	{
	    i = Perm [k] ;
	    if (i < 0 || i > n || Pinv [i] != EMPTY)
	    {
		ERROR (CHOLMOD_INVALID, "invalid permutation") ;
		return (FALSE) ;
	    }
	    Pinv [i] = k ;
	}
    }

    /* Perm is now valid */
    ASSERT (CHOLMOD(dump_perm) (Perm, n, n, "Perm", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* count the entries in each row of F */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < n ; i++)
    {
	Wi [i] = 0 ;
    }

    if (packed)
    {
	if (permute)
	{
	    if (upper)
	    {
		/* packed, permuted, upper */
		for (j = 0 ; j < n ; j++)
		{
		    jold = Perm [j] ;
		    pend = Ap [jold+1] ;
		    for (p = Ap [jold] ; p < pend ; p++)
		    {
			iold = Ai [p] ;
			if (iold <= jold)
			{
			    i = Pinv [iold] ;
			    Wi [MIN (i, j)]++ ;
			}
		    }
		}
	    }
	    else
	    {
		/* packed, permuted, lower */
		for (j = 0 ; j < n ; j++)
		{
		    jold = Perm [j] ;
		    pend = Ap [jold+1] ;
		    for (p = Ap [jold] ; p < pend ; p++)
		    {
			iold = Ai [p] ;
			if (iold >= jold)
			{
			    i = Pinv [iold] ;
			    Wi [MAX (i, j)]++ ;
			}
		    }
		}
	    }
	}
	else
	{
	    if (upper)
	    {
		/* packed, unpermuted, upper */
		for (j = 0 ; j < n ; j++)
		{
		    pend = Ap [j+1] ;
		    for (p = Ap [j] ; p < pend ; p++)
		    {
			i = Ai [p] ;
			if (i <= j)
			{
			    Wi [i]++ ;
			}
		    }
		}
	    }
	    else
	    {
		/* packed, unpermuted, lower */
		for (j = 0 ; j < n ; j++)
		{
		    pend = Ap [j+1] ;
		    for (p = Ap [j] ; p < pend ; p++)
		    {
			i = Ai [p] ;
			if (i >= j)
			{
			    Wi [i]++ ;
			}
		    }
		}
	    }
	}
    }
    else
    {
	if (permute)
	{
	    if (upper)
	    {
		/* unpacked, permuted, upper */
		for (j = 0 ; j < n ; j++)
		{
		    jold = Perm [j] ;
		    p = Ap [jold] ;
		    pend = p + Anz [jold] ;
		    for ( ; p < pend ; p++)
		    {
			iold = Ai [p] ;
			if (iold <= jold)
			{
			    i = Pinv [iold] ;
			    Wi [MIN (i, j)]++ ;
			}
		    }
		}
	    }
	    else
	    {
		/* unpacked, permuted, lower */
		for (j = 0 ; j < n ; j++)
		{
		    jold = Perm [j] ;
		    p = Ap [jold] ;
		    pend = p + Anz [jold] ;
		    for ( ; p < pend ; p++)
		    {
			iold = Ai [p] ;
			if (iold >= jold)
			{
			    i = Pinv [iold] ;
			    Wi [MAX (i, j)]++ ;
			}
		    }
		}
	    }
	}
	else
	{
	    if (upper)
	    {
		/* unpacked, unpermuted, upper */
		for (j = 0 ; j < n ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			i = Ai [p] ;
			if (i <= j)
			{
			    Wi [i]++ ;
			}
		    }
		}
	    }
	    else
	    {
		/* unpacked, unpermuted, lower */
		for (j = 0 ; j < n ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			i = Ai [p] ;
			if (i >= j)
			{
			    Wi [i]++ ;
			}
		    }
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* compute the row pointers */
    /* ---------------------------------------------------------------------- */

    p = 0 ;
    for (i = 0 ; i < n ; i++)
    {
	Fp [i] = p ;
	p += Wi [i] ;
    }
    Fp [n] = p ;
    for (i = 0 ; i < n ; i++)
    {
	Wi [i] = Fp [i] ;
    }

    if (p > (Int) (F->nzmax))
    {
	ERROR (CHOLMOD_INVALID, "F is too small") ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* transpose matrix, using template routine */
    /* ---------------------------------------------------------------------- */

    ok = FALSE ;
    if (values == 0 || F->xtype == CHOLMOD_PATTERN)
    {
	PRINT2 (("\n:::: p_transpose_sym Perm %p\n", Perm)) ;
	ok = p_cholmod_transpose_sym (A, Perm, F, Common) ;
    }
    else if (F->xtype == CHOLMOD_REAL)
    {
	PRINT2 (("\n:::: r_transpose_sym Perm %p\n", Perm)) ;
	ok = r_cholmod_transpose_sym (A, Perm, F, Common) ;
    }
    else if (F->xtype == CHOLMOD_COMPLEX)
    {
	if (values == 1)
	{
	    /* array transpose */
	    PRINT2 (("\n:::: ct_transpose_sym Perm %p\n", Perm)) ;
	    ok = ct_cholmod_transpose_sym (A, Perm, F, Common) ;
	}
	else
	{
	    /* complex conjugate transpose */
	    PRINT2 (("\n:::: c_transpose_sym Perm %p\n", Perm)) ;
	    ok = c_cholmod_transpose_sym (A, Perm, F, Common) ;
	}
    }
    else if (F->xtype == CHOLMOD_ZOMPLEX)
    {
	if (values == 1)
	{
	    /* array transpose */
	    PRINT2 (("\n:::: zt_transpose_sym Perm %p\n", Perm)) ;
	    ok = zt_cholmod_transpose_sym (A, Perm, F, Common) ;
	}
	else
	{
	    /* complex conjugate transpose */
	    PRINT2 (("\n:::: z_transpose_sym Perm %p\n", Perm)) ;
	    ok = z_cholmod_transpose_sym (A, Perm, F, Common) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* finalize result F */
    /* ---------------------------------------------------------------------- */

    /* F is sorted if there is no permutation vector */
    if (ok)
    {
	F->sorted = !permute ;
	F->packed = TRUE ;
	F->stype = - SIGN (A->stype) ;	/* flip the stype */
	ASSERT (CHOLMOD(dump_sparse) (F, "output F sym", Common) >= 0) ;
    }
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_transpose ==================================================== */
/* ========================================================================== */

/* Returns A'.  See also cholmod_ptranspose below. */

cholmod_sparse *CHOLMOD(transpose)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to transpose */
    int values,		/* 2: complex conj. transpose, 1: array transpose,
			   0: do not transpose the numerical values
			   (returns its result as CHOLMOD_PATTERN) */
    /* --------------- */
    cholmod_common *Common
)
{
    return (CHOLMOD(ptranspose) (A, values, NULL, NULL, 0, Common)) ;
}


/* ========================================================================== */
/* === cholmod_ptranspose =================================================== */
/* ========================================================================== */

/* Return A' or A(p,p)' if A is symmetric.  Return A', A(:,f)', or A(p,f)' if
 * A is unsymmetric.
 *
 * workspace:
 * Iwork (MAX (nrow,ncol)) if unsymmetric and fset is non-NULL
 * Iwork (nrow) if unsymmetric and fset is NULL
 * Iwork (2*nrow) if symmetric and Perm is non-NULL.
 * Iwork (nrow) if symmetric and Perm is NULL.
 *
 * A simple worst-case upper bound on the workspace is nrow+ncol.
 */

cholmod_sparse *CHOLMOD(ptranspose)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to transpose */
    int values,		/* 2: complex conj. transpose, 1: array transpose,
			   0: do not transpose the numerical values */
    Int *Perm,		/* if non-NULL, F = A(p,f) or A(p,p) */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* --------------- */
    cholmod_common *Common
)
{
    Int *Ap, *Anz ;
    cholmod_sparse *F ;
    Int nrow, ncol, use_fset, j, jj, fnz, packed, stype, nf, xtype ;
    size_t ineed ;
    int ok = TRUE ;

    nf = fsize ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, NULL) ;
    stype = A->stype ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    ncol = A->ncol ;

    if (stype != 0)
    {
	use_fset = FALSE ;
	if (Perm != NULL)
	{
	    ineed = CHOLMOD(mult_size_t) (A->nrow, 2, &ok) ;
	}
	else
	{
	    ineed = A->nrow ;
	}
    }
    else
    {
	use_fset = (fset != NULL) ;
	if (use_fset)
	{
	    ineed = MAX (A->nrow, A->ncol) ;
	}
	else
	{
	    ineed = A->nrow ;
	}
    }

    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (NULL) ;
    }

    CHOLMOD(allocate_work) (0, ineed, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Anz = A->nz ;
    packed = A->packed ;
    ASSERT (IMPLIES (!packed, Anz != NULL)) ;
    xtype = values ? A->xtype : CHOLMOD_PATTERN ;

    /* ---------------------------------------------------------------------- */
    /* allocate F */
    /* ---------------------------------------------------------------------- */

    /* determine # of nonzeros in F */
    if (stype != 0)
    {
	/* F=A' or F=A(p,p)', fset is ignored */
	fnz = CHOLMOD(nnz) (A, Common) ;
    }
    else
    {
	nf = (use_fset) ? nf : ncol ;
	if (use_fset)
	{
	    fnz = 0 ;
	    /* F=A(:,f)' or F=A(p,f)' */
	    for (jj = 0 ; jj < nf ; jj++)
	    {
		/* The fset is not yet checked; it will be thoroughly checked
		 * in cholmod_transpose_unsym.  For now, just make sure we don't
		 * access Ap and Anz out of bounds. */
		j = fset [jj] ;
		if (j >= 0 && j < ncol)
		{
		    fnz += packed ? (Ap [j+1] - Ap [j]) : MAX (0, Anz [j]) ;
		}
	    }
	}
	else
	{
	    /* F=A' or F=A(p,:)' */
	    fnz = CHOLMOD(nnz) (A, Common) ;
	}
    }

    /* F is ncol-by-nrow, fnz nonzeros, sorted unless f is present and unsorted,
     * packed, of opposite stype as A, and with/without numerical values */
    F = CHOLMOD(allocate_sparse) (ncol, nrow, fnz, TRUE, TRUE, -SIGN(stype),
	    xtype, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* transpose and optionally permute the matrix A */
    /* ---------------------------------------------------------------------- */

    if (stype != 0)
    {
	/* F = A (p,p)', using upper or lower triangular part of A only */
	ok = CHOLMOD(transpose_sym) (A, values, Perm, F, Common) ;
    }
    else
    {
	/* F = A (p,f)' */
	ok = CHOLMOD(transpose_unsym) (A, values, Perm, fset, nf, F, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* return the matrix F, or NULL if an error occured */
    /* ---------------------------------------------------------------------- */

    if (!ok)
    {
	CHOLMOD(free_sparse) (&F, Common) ;
    }
    return (F) ;
}


/* ========================================================================== */
/* === cholmod_sort ========================================================= */
/* ========================================================================== */

/* Sort the columns of A, in place.  Returns A in packed form, even if it
 * starts as unpacked.  Removes entries in the ignored part of a symmetric
 * matrix.
 *
 * workspace: Iwork (max (nrow,ncol)).  Allocates additional workspace for a
 * temporary copy of A'.
 */

int CHOLMOD(sort)
(
    /* ---- in/out --- */
    cholmod_sparse *A,	/* matrix to sort */
    /* --------------- */
    cholmod_common *Common
)
{
    Int *Ap ;
    cholmod_sparse *F ;
    Int anz, ncol, nrow, stype ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    Common->status = CHOLMOD_OK ;
    nrow = A->nrow ;
    if (nrow <= 1)
    {
	/* a 1-by-n sparse matrix must be sorted */
	A->sorted = TRUE ;
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    ncol = A->ncol ;
    CHOLMOD(allocate_work) (0, MAX (nrow, ncol), 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;	/* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    anz = CHOLMOD(nnz) (A, Common) ;
    stype = A->stype ;

    /* ---------------------------------------------------------------------- */
    /* sort the columns of the matrix */
    /* ---------------------------------------------------------------------- */

    /* allocate workspace for transpose: ncol-by-nrow, same # of nonzeros as A,
     * sorted, packed, same stype as A, and of the same numeric type as A. */
    F = CHOLMOD(allocate_sparse) (ncol, nrow, anz, TRUE, TRUE, stype,
	    A->xtype, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;	/* out of memory */
    }

    if (stype != 0)
    {
	/* F = A', upper or lower triangular part only */
	CHOLMOD(transpose_sym) (A, 1, NULL, F, Common) ;
	A->packed = TRUE ;
	/* A = F' */
	CHOLMOD(transpose_sym) (F, 1, NULL, A, Common) ;
    }
    else
    {
	/* F = A' */
	CHOLMOD(transpose_unsym) (A, 1, NULL, NULL, 0, F, Common) ;
	A->packed = TRUE ;
	/* A = F' */
	CHOLMOD(transpose_unsym) (F, 1, NULL, NULL, 0, A, Common) ;
    }

    ASSERT (A->sorted && A->packed) ;
    ASSERT (CHOLMOD(dump_sparse) (A, "Asorted", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* reduce A in size, if needed.  This must succeed. */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    anz = Ap [ncol] ;
    ASSERT ((size_t) anz <= A->nzmax) ;
    CHOLMOD(reallocate_sparse) (anz, A, Common) ;
    ASSERT (Common->status >= CHOLMOD_OK) ;

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(free_sparse) (&F, Common) ;
    return (TRUE) ;
}
