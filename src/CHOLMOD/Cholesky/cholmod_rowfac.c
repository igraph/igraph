/* ========================================================================== */
/* === Cholesky/cholmod_rowfac ============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Cholesky Module.  Copyright (C) 2005-2013, Timothy A. Davis
 * The CHOLMOD/Cholesky Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Full or incremental numerical LDL' or LL' factorization (simplicial, not
 * supernodal) cholmod_factorize is the "easy" wrapper for this code, but it
 * does not provide access to incremental factorization.
 *
 * cholmod_rowfac computes the full or incremental LDL' or LL' factorization of
 * A+beta*I (where A is symmetric) or A*F+beta*I (where A and F are unsymmetric
 * and only the upper triangular part of A*F+beta*I is used).  It computes
 * L (and D, for LDL') one row at a time.  beta is real.
 *
 * A is nrow-by-ncol or nrow-by-nrow.  In "packed" form it is a conventional
 * column-oriented sparse matrix.  Row indices of column j are in
 * Ai [Ap [j] ... Ap [j+1]-1] and values in the same locations of Ax.
 * will be faster if A has sorted columns.  In "unpacked" form the column
 * of A ends at Ap [j] + Anz [j] - 1 instead of Ap [j+1] - 1.
 *
 * Row indices in each column of A can be sorted or unsorted, but the routine
 * routine works fastest if A is sorted, or if only triu(A) is provided
 * for the symmetric case.
 *
 * The unit-diagonal nrow-by-nrow output matrix L is returned in "unpacked"
 * column form, with row indices of column j in Li [Lp [j] ...
 * Lp [j] + Lnz [j] - 1] and values in the same location in Lx.  The row
 * indices in each column of L are in sorted order.  The unit diagonal of L
 * is not stored.
 *
 * L can be a simplicial symbolic or numeric (L->is_super must be FALSE).
 * A symbolic factor is converted immediately into a numeric factor containing
 * the identity matrix.
 *
 * For a full factorization, kstart = 0 and kend = nrow.  The existing nonzero
 * entries (numerical values in L->x and L->z for the zomplex case, and indices
 * in L->i), if any, are overwritten.
 *
 * To compute an incremental factorization, select kstart and kend as the range
 * of rows of L you wish to compute.  A correct factorization will be computed
 * only if all descendants of all nodes k = kstart to kend-1 in the etree have
 * been factorized by a prior call to this routine, and if rows kstart to kend-1
 * have not been factorized.  This condition is NOT checked on input.
 *
 * ---------------
 * Symmetric case:
 * ---------------
 *
 *	The factorization (in MATLAB notation) is:
 *
 *	S = beta*I + A
 *	S = triu (S) + triu (S,1)'
 *	L*D*L' = S, or L*L' = S
 *
 *	A is a conventional sparse matrix in compressed column form.  Only the
 *	diagonal and upper triangular part of A is accessed; the lower
 *	triangular part is ignored and assumed to be equal to the upper
 *	triangular part.  For an incremental factorization, only columns kstart
 *	to kend-1 of A are accessed.  F is not used.
 *
 * ---------------
 * Unsymmetric case:
 * ---------------
 *
 *	The factorization (in MATLAB notation) is:
 *
 *	S = beta*I + A*F
 *	S = triu (S) + triu (S,1)'
 *	L*D*L' = S, or L*L' = S
 *
 *	The typical case is F=A'.  Alternatively, if F=A(:,f)', then this
 *	routine factorizes S = beta*I + A(:,f)*A(:,f)'.
 *
 *	All of A and F are accessed, but only the upper triangular part of A*F
 *	is used.  F must be of size A->ncol by A->nrow.  F is used for the
 *	unsymmetric case only.  F can be packed or unpacked and it need not be
 *	sorted.
 *
 *	For a complete factorization of beta*I + A*A',
 *	this routine performs a number of flops exactly equal to:
 *
 *	sum (for each column j of A) of (Anz (j)^2 + Anz (j)), to form S
 *	+
 *	sum (for each column j of L) of (Lnz (j)^2 + 3*Lnz (j)), to factorize S
 *
 *	where Anz (j) is the number of nonzeros in column j of A, and Lnz (j)
 *	is the number of nonzero in column j of L below the diagonal.
 *
 *
 * workspace: Flag (nrow), W (nrow if real, 2*nrow if complex/zomplex),
 * Iwork (nrow)
 *
 * Supports any xtype, except a pattern-only input matrix A cannot be
 * factorized.
 */

#ifndef NCHOLESKY

#include "cholmod_internal.h"
#include "cholmod_cholesky.h"

/* ========================================================================== */
/* === subtree ============================================================== */
/* ========================================================================== */

/* Compute the nonzero pattern of the sparse triangular solve Lx=b, where L in
 * this case is L(0:k-1,0:k-1), and b is a column of A.  This is done by
 * traversing the kth row-subtree of the elimination tree of L, starting from
 * each nonzero entry in b.  The pattern is returned postordered, and is valid
 * for a subsequent numerical triangular solve of Lx=b.  The elimination tree
 * can be provided in a Parent array, or extracted from the pattern of L itself.
 *
 * The pattern of x = inv(L)*b is returned in Stack [top...].
 * Also scatters b, or a multiple of b, into the work vector W.
 *
 * The SCATTER macro is defines how the numerical values of A or A*A' are to be
 * scattered.
 *
 * PARENT(i) is a macro the defines how the etree is accessed.  It is either:
 *	#define PARENT(i) Parent [i]
 *	#define PARENT(i) (Lnz [i] > 1) ? (Li [Lp [i] + 1]) : EMPTY
 */

#define SUBTREE \
    for ( ; p < pend ; p++) \
    { \
	i = Ai [p] ; \
	if (i <= k) \
	{ \
	    /* scatter the column of A, or A*A' into Wx and Wz */ \
	    SCATTER ; \
	    /* start at node i and traverse up the subtree, stop at node k */ \
	    for (len = 0 ; i < k && i != EMPTY && Flag [i] < mark ; i = parent) \
	    { \
		/* L(k,i) is nonzero, and seen for the first time */ \
		Stack [len++] = i ;	    /* place i on the stack */ \
		Flag [i] = mark ;	    /* mark i as visited */ \
		parent = PARENT (i) ;   /* traverse up the etree to the parent */ \
	    } \
	    /* move the path down to the bottom of the stack */ \
	    while (len > 0) \
	    { \
		Stack [--top] = Stack [--len] ; \
	    } \
	} \
	else if (sorted) \
	{ \
	    break ; \
	} \
    }


/* ========================================================================== */
/* === TEMPLATE ============================================================= */
/* ========================================================================== */

#define REAL
#include "t_cholmod_rowfac.c"
#define COMPLEX
#include "t_cholmod_rowfac.c"
#define ZOMPLEX
#include "t_cholmod_rowfac.c"

#define MASK
#define REAL
#include "t_cholmod_rowfac.c"
#define COMPLEX
#include "t_cholmod_rowfac.c"
#define ZOMPLEX
#include "t_cholmod_rowfac.c"
#undef MASK


/* ========================================================================== */
/* === cholmod_row_subtree ================================================== */
/* ========================================================================== */

/* Compute the nonzero pattern of the solution to the lower triangular system
 * L(0:k-1,0:k-1) * x = A (0:k-1,k) if A is symmetric, or
 * L(0:k-1,0:k-1) * x = A (0:k-1,:) * A (:,k)' if A is unsymmetric.
 * This gives the nonzero pattern of row k of L (excluding the diagonal).
 * The pattern is returned postordered.
 *
 * The symmetric case requires A to be in symmetric-upper form.
 *
 * The result is returned in R, a pre-allocated sparse matrix of size nrow-by-1,
 * with R->nzmax >= nrow.  R is assumed to be packed (Rnz [0] is not updated);
 * the number of entries in R is given by Rp [0].
 *
 * FUTURE WORK:  a very minor change to this routine could allow it to compute
 * the nonzero pattern of x for any system Lx=b.  The SUBTREE macro would need
 * to change, to eliminate its dependence on k.
 *
 * workspace: Flag (nrow)
 */

int CHOLMOD(row_subtree)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    cholmod_sparse *F,	/* used for A*A' case only. F=A' or A(:,f)' */
    size_t krow,	/* row k of L */
    Int *Parent,	/* elimination tree */
    /* ---- output --- */
    cholmod_sparse *R,	/* pattern of L(k,:), 1-by-n with R->nzmax >= n */
    /* --------------- */
    cholmod_common *Common
)
{
    Int *Rp, *Stack, *Flag, *Ap, *Ai, *Anz, *Fp, *Fi, *Fnz ;
    Int p, pend, parent, t, stype, nrow, k, pf, pfend, Fpacked, packed,
	sorted, top, len, i, mark ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (R, FALSE) ;
    RETURN_IF_NULL (Parent, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (R, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    stype = A->stype ;
    if (stype == 0)
    {
	RETURN_IF_NULL (F, FALSE) ;
	RETURN_IF_XTYPE_INVALID (F, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    }
    if (krow >= A->nrow)
    {
	ERROR (CHOLMOD_INVALID, "subtree: k invalid") ;
	return (FALSE) ;
    }
    if (R->ncol != 1 || A->nrow != R->nrow || A->nrow > R->nzmax)
    {
	ERROR (CHOLMOD_INVALID, "subtree: R invalid") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    CHOLMOD(allocate_work) (nrow, 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    if (stype > 0)
    {
	/* symmetric upper case: F is not needed.  It may be NULL */
	Fp = NULL ;
	Fi = NULL ;
	Fnz = NULL ;
	Fpacked = TRUE ;
    }
    else if (stype == 0)
    {
	/* unsymmetric case: F is required. */
	Fp = F->p ;
	Fi = F->i ;
	Fnz = F->nz ;
	Fpacked = F->packed ;
    }
    else
    {
	/* symmetric lower triangular form not supported */
	ERROR (CHOLMOD_INVALID, "symmetric lower not supported") ;
	return (FALSE) ;
    }

    Ap = A->p ;
    Ai = A->i ;
    Anz = A->nz ;
    packed = A->packed ;
    sorted = A->sorted ;

    k = krow ;
    Stack = R->i ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Flag = Common->Flag ;	/* size nrow, Flag [i] < mark must hold */
    /* mark = CHOLMOD(clear_flag) (Common) ; */
    CHOLMOD_CLEAR_FLAG (Common) ;
    mark = Common->mark ;

    /* ---------------------------------------------------------------------- */
    /* compute the pattern of L(k,:) */
    /* ---------------------------------------------------------------------- */

    top = nrow ;		/* Stack is empty */
    Flag [k] = mark ;		/* do not include diagonal entry in Stack */

#define SCATTER			/* do not scatter numerical values */
#define PARENT(i) Parent [i]	/* use Parent for etree */

    if (stype != 0)
    {
	/* scatter kth col of triu (A), get pattern L(k,:) */
	p = Ap [k] ;
	pend = (packed) ? (Ap [k+1]) : (p + Anz [k]) ;
	SUBTREE ;
    }
    else
    {
	/* scatter kth col of triu (beta*I+AA'), get pattern L(k,:) */
	pf = Fp [k] ;
	pfend = (Fpacked) ? (Fp [k+1]) : (pf + Fnz [k]) ;
	for ( ; pf < pfend ; pf++)
	{
	    /* get nonzero entry F (t,k) */
	    t = Fi [pf] ;
	    p = Ap [t] ;
	    pend = (packed) ? (Ap [t+1]) : (p + Anz [t]) ;
	    SUBTREE ;
	}
    }

#undef SCATTER
#undef PARENT

    /* shift the stack upwards, to the first part of R */
    len = nrow - top ;
    for (i = 0 ; i < len ; i++)
    {
	Stack [i] = Stack [top + i] ;
    }

    Rp = R->p ;
    Rp [0] = 0 ;
    Rp [1] = len ;
    R->sorted = FALSE ;

    CHOLMOD(clear_flag) (Common) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_lsolve_pattern =============================================== */
/* ========================================================================== */

/* Compute the nonzero pattern of Y=L\B.  L must be simplicial, and B
 * must be a single sparse column vector with B->stype = 0.  The values of
 * B are not used; it just specifies a nonzero pattern.  The pattern of
 * Y is not sorted, but is in topological order instead (suitable for a
 * sparse forward/backsolve).
 */

int CHOLMOD(lsolve_pattern)
(
    /* ---- input ---- */
    cholmod_sparse *B,	/* sparse right-hand-side (a single sparse column) */
    cholmod_factor *L,	/* the factor L from which parent(i) is derived */
    /* ---- output --- */
    cholmod_sparse *Yset,   /* pattern of Y=L\B, n-by-1 with Y->nzmax >= n */
    /* --------------- */
    cholmod_common *Common
)
{
    size_t krow ;
    RETURN_IF_NULL (B, FALSE) ;
    krow = B->nrow ;
    return (CHOLMOD(row_lsubtree) (B, NULL, 0, krow, L, Yset, Common)) ;
}


/* ========================================================================== */
/* === cholmod_row_lsubtree ================================================= */
/* ========================================================================== */

/* Identical to cholmod_row_subtree, except that the elimination tree is
 * obtained from L itself, as the first off-diagonal entry in each column.
 * L must be simplicial, not supernodal.
 *
 * If krow = A->nrow, then A must be a single sparse column vector, (A->stype
 * must be zero), and the nonzero pattern of x=L\b is computed, where b=A(:,0)
 * is the single sparse right-hand-side.  The inputs Fi and fnz are ignored.
 * See CHOLMOD(lsolve_pattern) above for a simpler interface for this case.
 */

int CHOLMOD(row_lsubtree)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    Int *Fi, size_t fnz,    /* nonzero pattern of kth row of A', not required
			     * for the symmetric case.  Need not be sorted. */
    size_t krow,	/* row k of L */
    cholmod_factor *L,	/* the factor L from which parent(i) is derived */
    /* ---- output --- */
    cholmod_sparse *R,	/* pattern of L(k,:), n-by-1 with R->nzmax >= n */
    /* --------------- */
    cholmod_common *Common
)
{
    Int *Rp, *Stack, *Flag, *Ap, *Ai, *Anz, *Lp, *Li, *Lnz ;
    Int p, pend, parent, t, stype, nrow, k, pf, packed, sorted, top, len, i,
	mark, ka ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (R, FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (R, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;

    nrow = A->nrow ;
    stype = A->stype ;
    if (stype < 0)
    {
	/* symmetric lower triangular form not supported */
	ERROR (CHOLMOD_INVALID, "symmetric lower not supported") ;
	return (FALSE) ;
    }

    if (krow > nrow)
    {
        ERROR (CHOLMOD_INVALID, "lsubtree: krow invalid") ;
        return (FALSE) ;
    }
    else if (krow == nrow)
    {
        /* find pattern of x=L\b where b=A(:,0) */
        k = nrow ;      /* compute all of the result; don't stop in SUBTREE */
        ka = 0 ;        /* use column A(:,0) */
        if (stype != 0 || A->ncol != 1)
        {
            /* A must be unsymmetric (it's a single sparse column vector) */
            ERROR (CHOLMOD_INVALID, "lsubtree: A invalid") ;
            return (FALSE) ;
        }
    }
    else
    {
        /* find pattern of L(k,:) using A(:,k) and Fi if A unsymmetric */
        k = krow ;      /* which row of L to compute */
        ka = k ;        /* which column of A to use */
        if (stype == 0)
        {
            RETURN_IF_NULL (Fi, FALSE) ;
        }
    }

    if (R->ncol != 1 || nrow != R->nrow || nrow > R->nzmax || ka >= A->ncol)
    {
	ERROR (CHOLMOD_INVALID, "lsubtree: R invalid") ;
	return (FALSE) ;
    }
    if (L->is_super)
    {
	ERROR (CHOLMOD_INVALID, "lsubtree: L invalid (cannot be supernodal)") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(allocate_work) (nrow, 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    Anz = A->nz ;
    packed = A->packed ;
    sorted = A->sorted ;

    Stack = R->i ;

    Lp = L->p ;
    Li = L->i ;
    Lnz = L->nz ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Flag = Common->Flag ;	/* size nrow, Flag [i] < mark must hold */
    mark = CHOLMOD(clear_flag) (Common) ;

    /* ---------------------------------------------------------------------- */
    /* compute the pattern of L(k,:) */
    /* ---------------------------------------------------------------------- */

    top = nrow ;		/* Stack is empty */
    if (k < nrow)
    {
        Flag [k] = mark ;       /* do not include diagonal entry in Stack */
    }

#define SCATTER			/* do not scatter numerical values */
#define PARENT(i) (Lnz [i] > 1) ? (Li [Lp [i] + 1]) : EMPTY

    if (krow == nrow || stype != 0)
    {
	/* scatter kth col of triu (A), get pattern L(k,:) */
	p = Ap [ka] ;
	pend = (packed) ? (Ap [ka+1]) : (p + Anz [ka]) ;
	SUBTREE ;
    }
    else
    {
	/* scatter kth col of triu (beta*I+AA'), get pattern L(k,:) */
	for (pf = 0 ; pf < (Int) fnz ; pf++)
	{
	    /* get nonzero entry F (t,k) */
	    t = Fi [pf] ;
	    p = Ap [t] ;
	    pend = (packed) ? (Ap [t+1]) : (p + Anz [t]) ;
	    SUBTREE ;
	}
    }

#undef SCATTER
#undef PARENT

    /* shift the stack upwards, to the first part of R */
    len = nrow - top ;
    for (i = 0 ; i < len ; i++)
    {
	Stack [i] = Stack [top + i] ;
    }

    Rp = R->p ;
    Rp [0] = 0 ;
    Rp [1] = len ;
    R->sorted = FALSE ;

    CHOLMOD(clear_flag) (Common) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_rowfac ======================================================= */
/* ========================================================================== */

/* This is the incremental factorization for general purpose usage. */

int CHOLMOD(rowfac)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to factorize */
    cholmod_sparse *F,	/* used for A*A' case only. F=A' or A(:,f)' */
    double beta [2],	/* factorize beta*I+A or beta*I+AA' */
    size_t kstart,	/* first row to factorize */
    size_t kend,	/* last row to factorize is kend-1 */
    /* ---- in/out --- */
    cholmod_factor *L,
    /* --------------- */
    cholmod_common *Common
)
{
    return (CHOLMOD(rowfac_mask) (A, F, beta, kstart, kend, NULL, NULL, L,
	Common)) ;
}


/* ========================================================================== */
/* === cholmod_rowfac_mask ================================================== */
/* ========================================================================== */

/* This is meant for use in LPDASA only. */

int CHOLMOD(rowfac_mask)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to factorize */
    cholmod_sparse *F,	/* used for A*A' case only. F=A' or A(:,f)' */
    double beta [2],	/* factorize beta*I+A or beta*I+AA' */
    size_t kstart,	/* first row to factorize */
    size_t kend,	/* last row to factorize is kend-1 */
    Int *mask,		/* size A->nrow. if mask[i] >= 0 row i is set to zero */
    Int *RLinkUp,	/* size A->nrow. link list of rows to compute */
    /* ---- in/out --- */
    cholmod_factor *L,
    /* --------------- */
    cholmod_common *Common
)
{
    Int n ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    if (L->xtype != CHOLMOD_PATTERN && A->xtype != L->xtype)
    {
	ERROR (CHOLMOD_INVALID, "xtype of A and L do not match") ;
	return (FALSE) ;
    }
    if (L->is_super)
    {
	ERROR (CHOLMOD_INVALID, "can only do simplicial factorization");
	return (FALSE) ;
    }
    if (A->stype == 0)
    {
	RETURN_IF_NULL (F, FALSE) ;
	if (A->xtype != F->xtype)
	{
	    ERROR (CHOLMOD_INVALID, "xtype of A and F do not match") ;
	    return (FALSE) ;
	}
    }
    if (A->stype < 0)
    {
	/* symmetric lower triangular form not supported */
	ERROR (CHOLMOD_INVALID, "symmetric lower not supported") ;
	return (FALSE) ;
    }
    if (kend > L->n)
    {
	ERROR (CHOLMOD_INVALID, "kend invalid") ;
	return (FALSE) ;
    }
    if (A->nrow != L->n)
    {
	ERROR (CHOLMOD_INVALID, "dimensions of A and L do not match") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;
    Common->rowfacfl = 0 ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* Xwork is of size n for the real case, 2*n for complex/zomplex */
    n = L->n  ;

    /* s = ((A->xtype != CHOLMOD_REAL) ? 2:1)*n */
    s = CHOLMOD(mult_size_t) (n, ((A->xtype != CHOLMOD_REAL) ? 2:1), &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    CHOLMOD(allocate_work) (n, n, s, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, A->nrow, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* factorize the matrix, using template routine */
    /* ---------------------------------------------------------------------- */

    if (RLinkUp == NULL)
    {

	switch (A->xtype)
	{
	    case CHOLMOD_REAL:
		ok = r_cholmod_rowfac (A, F, beta, kstart, kend, L, Common) ;
		break ;

	    case CHOLMOD_COMPLEX:
		ok = c_cholmod_rowfac (A, F, beta, kstart, kend, L, Common) ;
		break ;

	    case CHOLMOD_ZOMPLEX:
		ok = z_cholmod_rowfac (A, F, beta, kstart, kend, L, Common) ;
		break ;
	}

    }
    else
    {

	switch (A->xtype)
	{
	    case CHOLMOD_REAL:
		ok = r_cholmod_rowfac_mask (A, F, beta, kstart, kend,
		    mask, RLinkUp, L, Common) ;
		break ;

	    case CHOLMOD_COMPLEX:
		ok = c_cholmod_rowfac_mask (A, F, beta, kstart, kend,
		    mask, RLinkUp, L, Common) ;
		break ;

	    case CHOLMOD_ZOMPLEX:
		ok = z_cholmod_rowfac_mask (A, F, beta, kstart, kend,
		    mask, RLinkUp, L, Common) ;
		break ;
	}
    }

    return (ok) ;
}
#endif
