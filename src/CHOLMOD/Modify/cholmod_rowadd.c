/* ========================================================================== */
/* === Modify/cholmod_rowadd ================================================ */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Modify Module.
 * Copyright (C) 2005-2006, Timothy A. Davis and William W. Hager.
 * The CHOLMOD/Modify Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* Adds a row and column to an LDL' factorization, and optionally updates the
 * solution to Lx=b.
 *
 * workspace: Flag (nrow), Head (nrow+1), W (2*nrow), Iwork (2*nrow)
 *
 * Only real matrices are supported.  A symbolic L is converted into a
 * numeric identity matrix before the row is added.
 */

#ifndef NMODIFY

#include "cholmod_internal.h"
#include "cholmod_modify.h"


/* ========================================================================== */
/* === cholmod_rowadd ======================================================= */
/* ========================================================================== */

/* cholmod_rowadd adds a row to the LDL' factorization.  It computes the kth
 * row and kth column of L, and then updates the submatrix L (k+1:n,k+1:n)
 * accordingly.  The kth row and column of L should originally be equal to the
 * kth row and column of the identity matrix (they are treated as such, if they
 * are not).  The kth row/column of L is computed as the factorization of the
 * kth row/column of the matrix to factorize, which is provided as a single
 * n-by-1 sparse matrix R.  The sparse vector R need not be sorted.
 */

int CHOLMOD(rowadd)
(
    /* ---- input ---- */
    size_t k,		/* row/column index to add */
    cholmod_sparse *R,	/* row/column of matrix to factorize (n-by-1) */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    /* --------------- */
    cholmod_common *Common
)
{
    double bk [2] ;
    bk [0] = 0. ;
    bk [1] = 0. ;
    return (CHOLMOD(rowadd_mark) (k, R, bk, NULL, L, NULL, NULL, Common)) ;
}


/* ========================================================================== */
/* === cholmod_rowadd_solve ================================================= */
/* ========================================================================== */

/* Does the same as cholmod_rowadd, and also updates the solution to Lx=b
 * See cholmod_updown for a description of how Lx=b is updated.  There is on
 * additional parameter:  bk specifies the new kth entry of b.
 */

int CHOLMOD(rowadd_solve)
(
    /* ---- input ---- */
    size_t k,		/* row/column index to add */
    cholmod_sparse *R,	/* row/column of matrix to factorize (n-by-1) */
    double bk [2],	/* kth entry of the right-hand-side b */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    cholmod_dense *X,	/* solution to Lx=b (size n-by-1) */
    cholmod_dense *DeltaB,  /* change in b, zero on output */
    /* --------------- */
    cholmod_common *Common
)
{
    return (CHOLMOD(rowadd_mark) (k, R, bk, NULL, L, X, DeltaB, Common)) ;
}


/* ========================================================================== */
/* === icomp ================================================================ */
/* ========================================================================== */

/* for sorting by qsort */
static int icomp (Int *i, Int *j)
{
    if (*i < *j)
    {
	return (-1) ;
    }
    else
    {
	return (1) ;
    }
}


/* ========================================================================== */
/* === cholmod_rowadd_mark ================================================== */
/* ========================================================================== */

/* Does the same as cholmod_rowadd_solve, except only part of L is used in
 * the update/downdate of the solution to Lx=b.  This routine is an "expert"
 * routine.  It is meant for use in LPDASA only.  */

int CHOLMOD(rowadd_mark)
(
    /* ---- input ---- */
    size_t kadd,	/* row/column index to add */
    cholmod_sparse *R,	/* row/column of matrix to factorize (n-by-1) */
    double bk [2],	/* kth entry of the right hand side, b */
    Int *colmark,	/* Int array of size 1.  See cholmod_updown.c */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    cholmod_dense *X,	/* solution to Lx=b (size n-by-1) */
    cholmod_dense *DeltaB,  /* change in b, zero on output */
    /* --------------- */
    cholmod_common *Common
)
{
    double dk, yj, l_kj, lx, l_ij, sqrt_dk, dj, xk, rnz, fl ;
    double *Lx, *W, *Cx, *Rx, *Xx, *Nx ;
    Int *Li, *Lp, *Lnz, *Flag, *Stack, *Ci, *Rj, *Rp, *Lnext, *Iwork, *Rnz ;
    cholmod_sparse *C, Cmatrix ;
    Int i, j, p, pend, top, len, kk, li, lnz, mark, k, n, parent, Cp [2],
	do_solve, do_update ;
    size_t s ;
    int ok = TRUE ;
    DEBUG (Int lastrow) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_NULL (R, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_PATTERN, CHOLMOD_REAL, FALSE) ;
    RETURN_IF_XTYPE_INVALID (R, CHOLMOD_REAL, CHOLMOD_REAL, FALSE) ;
    n = L->n ;
    k = kadd ;
    if (kadd >= L->n || k < 0)
    {
	ERROR (CHOLMOD_INVALID, "k invalid") ;
	return (FALSE) ;
    }
    if (R->ncol != 1 || R->nrow != L->n)
    {
	ERROR (CHOLMOD_INVALID, "R invalid") ;
	return (FALSE) ;
    }
    Rj = R->i ;
    Rx = R->x ;
    Rp = R->p ;
    Rnz = R->nz ;
    rnz = (R->packed) ? (Rp [1]) : (Rnz [0]) ;
    do_solve = (X != NULL) && (DeltaB != NULL) ;
    if (do_solve)
    {
	RETURN_IF_XTYPE_INVALID (X, CHOLMOD_REAL, CHOLMOD_REAL, FALSE) ;
	RETURN_IF_XTYPE_INVALID (DeltaB, CHOLMOD_REAL, CHOLMOD_REAL, FALSE) ;
	Xx = X->x ;
	Nx = DeltaB->x ;
	if (X->nrow != L->n || X->ncol != 1 || DeltaB->nrow != L->n ||
		DeltaB->ncol != 1 || Xx == NULL || Nx == NULL)
	{
	    ERROR (CHOLMOD_INVALID, "X and/or DeltaB invalid") ;
	    return (FALSE) ;
	}
    }
    else
    {
	Xx = NULL ;
	Nx = NULL ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* s = 2*n */
    s = CHOLMOD(mult_size_t) (n, 2, &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    CHOLMOD(allocate_work) (n, s, s, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, s, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* convert to simplicial numeric LDL' factor, if not already */
    /* ---------------------------------------------------------------------- */

    if (L->xtype == CHOLMOD_PATTERN || L->is_super || L->is_ll) 
    {
	/* can only update/downdate a simplicial LDL' factorization */
	CHOLMOD(change_factor) (CHOLMOD_REAL, FALSE, FALSE, FALSE, FALSE, L,
		Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory, L is returned unchanged */
	    return (FALSE) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    /* inputs, not modified on output: */
    Lp = L->p ;		/* size n+1.  input, not modified on output */

    /* outputs, contents defined on input for incremental case only: */
    Lnz = L->nz ;	/* size n */
    Li = L->i ;		/* size L->nzmax.  Can change in size. */
    Lx = L->x ;		/* size L->nzmax.  Can change in size. */
    Lnext = L->next ;	/* size n+2 */

    ASSERT (L->nz != NULL) ;

    PRINT1 (("rowadd:\n")) ;
    fl = 0 ;

#if 0
#ifndef NDEBUG
    /* column k of L should be zero, except for the diagonal.  This test is
     * overly cautious. */
    for (p = Lp [k] + 1 ; p < Lp [k] + Lnz [k] ; p++) ASSERT (Lx [p] == 0) ;
#endif
#endif

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Flag = Common->Flag ;   /* size n */
    W = Common->Xwork ;     /* size n */
    Cx = W + n ;	    /* size n (use 2nd column of Xwork for C) */
    Iwork = Common->Iwork ;
    Stack = Iwork ;	    /* size n (i/i/l), also in cholmod_updown */
    Ci = Iwork + n ;	    /* size n (i/i/l) */
    /* NOTE: cholmod_updown uses Iwork [0..n-1] (i/i/l) as Stack as well */

    mark = Common->mark ;

    /* copy Rj/Rx into W/Ci */
    for (p = 0 ; p < rnz ; p++)
    {
	i = Rj [p] ;
	ASSERT (i >= 0 && i < n) ;
	W [i] = Rx [p] ;
	Ci [p] = i ;
    }

    /* At this point, W [Ci [0..rnz-1]] holds the sparse vector to add */
    /* The nonzero pattern of column W is held in Ci (it may be unsorted). */

    /* ---------------------------------------------------------------------- */
    /* symbolic factorization to get pattern of kth row of L */
    /* ---------------------------------------------------------------------- */

    DEBUG (for (p = 0 ; p < rnz ; p++)
	    PRINT1 (("C ("ID",%g)\n", Ci [p], W [Ci [p]]))) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* flag the diagonal */
    Flag [k] = mark ;

    /* find the union of all the paths */
    top = n ;
    lnz = 0 ;	/* # of nonzeros in column k of L, excluding diagonal */
    for (p = 0 ; p < rnz ; p++)
    {
	i = Ci [p] ;

	if (i < k)
	{

	    /* walk from i = entry in Ci to root (and stop if i marked)*/
	    PRINT2 (("\nwalk from i = "ID" towards k = "ID"\n", i, k)) ;
	    len = 0 ;

	    /* walk up tree, but stop if we go below the diagonal */
	    while (i < k && i != EMPTY && Flag [i] < mark)
	    {
		PRINT2 (("   Add "ID" to path\n", i)) ;
		ASSERT (i >= 0 && i < k) ;
		Stack [len++] = i ;	/* place i on the stack */
		Flag [i] = mark ;		/* mark i as visited */
		/* parent is the first entry in the column after the diagonal */
		ASSERT (Lnz [i] > 0) ;
		parent = (Lnz [i] > 1) ? (Li [Lp [i] + 1]) : EMPTY ;
		PRINT2 (("                      parent: "ID"\n", parent)) ;
		i = parent ;	/* go up the tree */
	    }
	    ASSERT (len <= top) ;

	    /* move the path down to the bottom of the stack */
	    /* this shifts Stack [0..len-1] down to [ ... oldtop-1] */
	    while (len > 0)
	    {
		Stack [--top] = Stack [--len] ;
	    }
	}
	else if (i > k)
	{
	    /* prune the diagonal and upper triangular entries from Ci */
	    Ci [lnz++] = i ;
	    Flag [i] = mark ;
	}
    }

#ifndef NDEBUG
    PRINT1 (("length of S after prune: "ID"\n", lnz)) ;
    for (p = 0 ; p < lnz ; p++)
    {
	PRINT1 (("After prune Ci ["ID"] = "ID"\n", p, Ci [p])) ;
	ASSERT (Ci [p] > k) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* ensure each column of L has enough space to grow */
    /* ---------------------------------------------------------------------- */

    for (kk = top ; kk < n ; kk++)
    {
	/* could skip this if we knew column j already included row k */
	j = Stack [kk] ;
	if (Lp [j] + Lnz [j] >= Lp [Lnext [j]])
	{
	    PRINT1 (("Col "ID" realloc, old Lnz "ID"\n", j, Lnz [j])) ;
	    if (!CHOLMOD(reallocate_column) (j, Lnz [j] + 1, L, Common))
	    {
		/* out of memory, L is now simplicial symbolic */
		/* CHOLMOD(clear_flag) (Common) ; */
		CHOLMOD_CLEAR_FLAG (Common) ;
		for (i = 0 ; i < n ; i++)
		{
		    W [i] = 0 ;
		}
		return (FALSE) ;
	    }
	    /* L->i and L->x may have moved */
	    Li = L->i ;
	    Lx = L->x ;
	}
	ASSERT (Lp [j] + Lnz [j] < Lp [Lnext [j]]
	    || (Lp [Lnext [j]] - Lp [j] == n-j)) ;
    }

    /* ---------------------------------------------------------------------- */
    /* compute kth row of L and store in column form */
    /* ---------------------------------------------------------------------- */

    /* solve L (1:k-1, 1:k-1) * y (1:k-1) = b (1:k-1) */
    /* where b (1:k) is in W and Ci */

    /* L (k, 1:k-1) = y (1:k-1) ./ D (1:k-1) */
    /* D (k) = B (k,k) - L (k, 1:k-1) * y (1:k-1) */

    PRINT2 (("\nForward solve: "ID" to "ID"\n", top, n)) ;
    ASSERT (Lnz [k] >= 1 && Li [Lp [k]] == k) ;
    DEBUG (for (i = top ; i < n ; i++) PRINT2 ((" Path: "ID"\n", Stack [i]))) ;

    dk = W [k] ;
    W [k] = 0.0 ;

    /* if do_solve: compute x (k) = b (k) - L (k, 1:k-1) * x (1:k-1) */
    xk = bk [0] ;
    PRINT2 (("B [k] = %g\n", xk)) ;

    for (kk = top ; kk < n ; kk++)
    {
	j = Stack [kk] ;
	i = j ;
	PRINT2 (("Forward solve col j = "ID":\n", j)) ;
	ASSERT (j >= 0 && j < k) ;

	/* forward solve using L (j+1:k-1,j) */
	yj = W [j] ;
	W [j] = 0.0 ;
	p = Lp [j] ;
	pend = p + Lnz [j] ;
	ASSERT (Lnz [j] > 0) ;
	dj = Lx [p++] ;
	for ( ; p < pend ; p++)
	{
	    i = Li [p] ;
	    PRINT2 (("    row "ID"\n", i)) ;
	    ASSERT (i > j) ;
	    ASSERT (i < n) ;
	    /* stop at row k */
	    if (i >= k)
	    {
		break ;
	    }
	    W [i] -= Lx [p] * yj ;
	}

	/* each iteration of the above for loop did 2 flops, and 3 flops
	 * are done below.  so: fl += 2 * (Lp [j] - p - 1) + 3 becomes: */
	fl += 2 * (Lp [j] - p) + 1 ;

	/* scale L (k,1:k-1) and compute dot product for D (k,k) */
	l_kj = yj / dj ;
	dk -= l_kj * yj ;

	/* compute dot product for X(k) */
	if (do_solve)
	{
	    xk -= l_kj * Xx [j] ;
	}

	/* store l_kj in the jth column of L */
	/* and shift the rest of the column down */

	li = k ;
	lx = l_kj ;

	if (i == k)
	{
	    /* no need to modify the nonzero pattern of L, since it already
	     * contains row index k. */
	    ASSERT (Li [p] == k) ;
	    Lx [p] = l_kj ;

	    for (p++ ; p < pend ; p++)
	    {
		i    = Li [p] ;
		l_ij = Lx [p] ;
		ASSERT (i > k && i < n) ;
		PRINT2 (("   apply to row "ID" of column k of L\n", i)) ;

		/* add to the pattern of the kth column of L */
		if (Flag [i] < mark)
		{
		    PRINT2 (("   add Ci["ID"] = "ID"\n", lnz, i)) ;
		    ASSERT (i > k) ;
		    Ci [lnz++] = i ;
		    Flag [i] = mark ;
		}

		/* apply the update to the kth column of L */
		/* yj is equal to l_kj * d_j */

		W [i] -= l_ij * yj ;
	    }

	}
	else
	{

	    PRINT2 (("Shift col j = "ID", apply saxpy to col k of L\n", j)) ;
	    for ( ; p < pend ; p++)
	    {
		/* swap (Li [p],Lx [p]) with (li,lx) */
		i    = Li [p] ;
		l_ij = Lx [p] ;
		Li [p] = li ;
		Lx [p] = lx ;
		li = i ;
		lx = l_ij ;
		ASSERT (i > k && i < n) ;
		PRINT2 (("   apply to row "ID" of column k of L\n", i)) ;

		/* add to the pattern of the kth column of L */
		if (Flag [i] < mark)
		{
		    PRINT2 (("   add Ci["ID"] = "ID"\n", lnz, i)) ;
		    ASSERT (i > k) ;
		    Ci [lnz++] = i ;
		    Flag [i] = mark ;
		}

		/* apply the update to the kth column of L */
		/* yj is equal to l_kj * d_j */

		W [i] -= l_ij * yj ;
	    }

	    /* store the last value in the jth column of L */
	    Li [p] = li ;
	    Lx [p] = lx ;
	    Lnz [j]++ ;

	}
    }

    /* ---------------------------------------------------------------------- */
    /* merge C with the pattern of the existing column of L */
    /* ---------------------------------------------------------------------- */

    /* This column should be zero, but it may contain explicit zero entries.
     * These entries should be kept, not dropped. */
    p = Lp [k] ;
    pend = p + Lnz [k] ;
    for (p++ ; p < pend ; p++)
    {
	i = Li [p] ;
	/* add to the pattern of the kth column of L */
	if (Flag [i] < mark)
	{
	    PRINT2 (("   add Ci["ID"] = "ID" from existing col k\n", lnz, i)) ;
	    ASSERT (i > k) ;
	    Ci [lnz++] = i ;
	    Flag [i] = mark ;
	}
    }

    /* ---------------------------------------------------------------------- */

    if (do_solve)
    {
	Xx [k] = xk ;
	PRINT2 (("Xx [k] = %g\n", Xx [k])) ;
    }

    /* ---------------------------------------------------------------------- */
    /* ensure abs (dk) >= dbound, if dbound is given */
    /* ---------------------------------------------------------------------- */

    dk = (IS_GT_ZERO (Common->dbound)) ? (CHOLMOD(dbound) (dk, Common)) : dk ;

    PRINT2 (("D [k = "ID"] = %g\n", k, dk)) ;

    /* ---------------------------------------------------------------------- */
    /* store the kth column of L */
    /* ---------------------------------------------------------------------- */

    /* ensure the new column of L has enough space */
    if (Lp [k] + lnz + 1 > Lp [Lnext [k]])
    {
	PRINT1 (("New Col "ID" realloc, old Lnz "ID"\n", k, Lnz [k])) ;
	if (!CHOLMOD(reallocate_column) (k, lnz + 1, L, Common))
	{
	    /* out of memory, L is now simplicial symbolic */
	    CHOLMOD(clear_flag) (Common) ;
	    for (i = 0 ; i < n ; i++)
	    {
		W [i] = 0 ;
	    }
	    return (FALSE) ;
	}
	/* L->i and L->x may have moved */
	Li = L->i ;
	Lx = L->x ;
    }
    ASSERT (Lp [k] + lnz + 1 <= Lp [Lnext [k]]) ;

#ifndef NDEBUG
    PRINT2 (("\nPrior to sort: lnz "ID" (excluding diagonal)\n", lnz)) ;
    for (kk = 0 ; kk < lnz ; kk++)
    {
	i = Ci [kk] ;
	PRINT2 (("L ["ID"] kept: "ID" %e\n", kk, i, W [i] / dk)) ;
    }
#endif

    /* sort Ci */
    qsort (Ci, lnz, sizeof (Int), (int (*) (const void *, const void *)) icomp);

    /* store the kth column of L */
    DEBUG (lastrow = k) ;
    p = Lp [k] ;
    Lx [p++] = dk ;
    Lnz [k] = lnz + 1 ;
    fl += lnz ;
    for (kk = 0 ; kk < lnz ; kk++, p++)
    {
	i = Ci [kk] ;
	PRINT2 (("L ["ID"] after sort: "ID", %e\n", kk, i, W [i] / dk)) ;
	ASSERT (i > lastrow) ;
	Li [p] = i ;
	Lx [p] = W [i] / dk ;
	W [i] = 0.0 ;
	DEBUG (lastrow = i) ;
    }

    /* compute DeltaB for updown (in DeltaB) */
    if (do_solve)
    {
	p = Lp [k] ;
	pend = p + Lnz [k] ;
	for (p++ ; p < pend ; p++)
	{
	    ASSERT (Li [p] > k) ;
	    Nx [Li [p]] -= Lx [p] * xk ;
	}
    }

    /* clear the flag for the update */
    mark = CHOLMOD(clear_flag) (Common) ;

    /* workspaces are now cleared */
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 2*n, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* update/downdate */
    /* ---------------------------------------------------------------------- */

    /* update or downdate L (k+1:n, k+1:n) with the vector
     * C = L (:,k) * sqrt (abs (D [k])).
     * Do a numeric update if D[k] < 0, numeric downdate otherwise.
     */

    ok = TRUE ;
    Common->modfl = 0 ;

    PRINT1 (("rowadd update lnz = "ID"\n", lnz)) ;
    if (lnz > 0)
    {
	do_update = IS_LT_ZERO (dk) ;
	if (do_update)
	{
	    dk = -dk ;
	}
	sqrt_dk = sqrt (dk) ;
	p = Lp [k] + 1 ;
	for (kk = 0 ; kk < lnz ; kk++, p++)
	{
	    Cx [kk] = Lx [p] * sqrt_dk ;
	}
	fl += lnz + 1 ;

	/* create a n-by-1 sparse matrix to hold the single column */
	C = &Cmatrix ;
	C->nrow = n ;
	C->ncol = 1 ;
	C->nzmax = lnz ;
	C->sorted = TRUE ;
	C->packed = TRUE ;
	C->p = Cp ;
	C->i = Ci ;
	C->x = Cx ;
	C->nz = NULL ;
	C->itype = L->itype ;
	C->xtype = L->xtype ;
	C->dtype = L->dtype ;
	C->z = NULL ;
	C->stype = 0 ;

	Cp [0] = 0 ;
	Cp [1] = lnz ;

	/* numeric downdate if dk > 0, and optional Lx=b change */
	/* workspace: Flag (nrow), Head (nrow+1), W (nrow), Iwork (2*nrow) */
	ok = CHOLMOD(updown_mark) (do_update ? (1) : (0), C, colmark,
		L, X, DeltaB, Common) ;

	/* clear workspace */
	for (kk = 0 ; kk < lnz ; kk++)
	{
	    Cx [kk] = 0 ;
	}
    }

    Common->modfl += fl ;

    DEBUG (CHOLMOD(dump_factor) (L, "LDL factorization, L:", Common)) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 2*n, Common)) ;
    return (ok) ;
}
#endif
