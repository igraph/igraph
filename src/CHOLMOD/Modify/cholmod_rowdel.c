/* ========================================================================== */
/* === Modify/cholmod_rowdel ================================================ */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Modify Module.
 * Copyright (C) 2005-2006, Timothy A. Davis and William W. Hager.
 * The CHOLMOD/Modify Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* Deletes a row and column from an LDL' factorization.  The row and column k
 * is set to the kth row and column of the identity matrix.  Optionally
 * downdates the solution to Lx=b.
 *
 * workspace: Flag (nrow), Head (nrow+1), W (nrow*2), Iwork (2*nrow)
 *
 * Only real matrices are supported (exception: since only the pattern of R
 * is used, it can have any valid xtype).
 */

#ifndef NMODIFY

#include "cholmod_internal.h"
#include "cholmod_modify.h"


/* ========================================================================== */
/* === cholmod_rowdel ======================================================= */
/* ========================================================================== */

/* Sets the kth row and column of L to be the kth row and column of the identity
 * matrix, and updates L(k+1:n,k+1:n) accordingly.   To reduce the running time,
 * the caller can optionally provide the nonzero pattern (or an upper bound) of
 * kth row of L, as the sparse n-by-1 vector R.  Provide R as NULL if you want
 * CHOLMOD to determine this itself, which is easier for the caller, but takes
 * a little more time.
 */

int CHOLMOD(rowdel)
(
    /* ---- input ---- */
    size_t k,		/* row/column index to delete */
    cholmod_sparse *R,	/* NULL, or the nonzero pattern of kth row of L */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    /* --------------- */
    cholmod_common *Common
)
{
    double yk [2] ;
    yk [0] = 0. ;
    yk [1] = 0. ;
    return (CHOLMOD(rowdel_mark) (k, R, yk, NULL, L, NULL, NULL, Common)) ;
}


/* ========================================================================== */
/* === cholmod_rowdel_solve ================================================= */
/* ========================================================================== */

/* Does the same as cholmod_rowdel, but also downdates the solution to Lx=b.
 * When row/column k of A is "deleted" from the system A*y=b, this can induce
 * a change to x, in addition to changes arising when L and b are modified.
 * If this is the case, the kth entry of y is required as input (yk) */

int CHOLMOD(rowdel_solve)
(
    /* ---- input ---- */
    size_t k,		/* row/column index to delete */
    cholmod_sparse *R,	/* NULL, or the nonzero pattern of kth row of L */
    double yk [2],	/* kth entry in the solution to A*y=b */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    cholmod_dense *X,	/* solution to Lx=b (size n-by-1) */
    cholmod_dense *DeltaB,  /* change in b, zero on output */
    /* --------------- */
    cholmod_common *Common
)
{
    return (CHOLMOD(rowdel_mark) (k, R, yk, NULL, L, X, DeltaB, Common)) ;
}


/* ========================================================================== */
/* === cholmod_rowdel_mark ================================================== */
/* ========================================================================== */

/* Does the same as cholmod_rowdel_solve, except only part of L is used in
 * the update/downdate of the solution to Lx=b.  This routine is an "expert"
 * routine.  It is meant for use in LPDASA only.
 *
 * if R == NULL then columns 0:k-1 of L are searched for row k.  Otherwise, it
 * searches columns in the set defined by the pattern of the first column of R.
 * This is meant to be the pattern of row k of L (a superset of that pattern is
 * OK too).  R must be a permutation of a subset of 0:k-1.
 */

int CHOLMOD(rowdel_mark)
(
    /* ---- input ---- */
    size_t kdel,	/* row/column index to delete */
    cholmod_sparse *R,	/* NULL, or the nonzero pattern of kth row of L */
    double yk [2],	/* kth entry in the solution to A*y=b */
    Int *colmark,	/* Int array of size 1.  See cholmod_updown.c */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    cholmod_dense *X,	/* solution to Lx=b (size n-by-1) */
    cholmod_dense *DeltaB,  /* change in b, zero on output */
    /* --------------- */
    cholmod_common *Common
)
{
    double dk, sqrt_dk, xk, dj, fl ;
    double *Lx, *Cx, *W, *Xx, *Nx ;
    Int *Li, *Lp, *Lnz, *Ci, *Rj, *Rp, *Iwork ;
    cholmod_sparse *C, Cmatrix ;
    Int j, p, pend, kk, lnz, n, Cp [2], do_solve, do_update, left, k,
	right, middle, i, klast, given_row, rnz ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_PATTERN, CHOLMOD_REAL, FALSE) ;
    n = L->n ;
    k = kdel ;
    if (kdel >= L->n || k < 0)
    {
	ERROR (CHOLMOD_INVALID, "k invalid") ;
	return (FALSE) ;
    }
    if (R == NULL)
    {
	Rj = NULL ;
	rnz = EMPTY ;
    }
    else
    {
	RETURN_IF_XTYPE_INVALID (R, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
	if (R->ncol != 1 || R->nrow != L->n)
	{
	    ERROR (CHOLMOD_INVALID, "R invalid") ;
	    return (FALSE) ;
	}
	Rj = R->i ;
	Rp = R->p ;
	rnz = Rp [1] ;
    }
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
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 2*n, Common)) ;

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
    Lp = L->p ;		/* size n+1 */

    /* outputs, contents defined on input for incremental case only: */
    Lnz = L->nz ;	/* size n */
    Li = L->i ;		/* size L->nzmax.  Can change in size. */
    Lx = L->x ;		/* size L->nzmax.  Can change in size. */

    ASSERT (L->nz != NULL) ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    W = Common->Xwork ; 	/* size n, used only in cholmod_updown */
    Cx = W + n ;		/* use 2nd column of Xwork for C (size n) */
    Iwork = Common->Iwork ;
    Ci = Iwork + n ;		/* size n (i/i/l) */
    /* NOTE: cholmod_updown uses Iwork [0..n-1] (i/i/l) as Stack */

    /* ---------------------------------------------------------------------- */
    /* prune row k from all columns of L */
    /* ---------------------------------------------------------------------- */

    given_row = (rnz >= 0) ;
    klast = given_row ? rnz : k ;
    PRINT2 (("given_row "ID"\n", given_row)) ;

    for (kk = 0 ; kk < klast ; kk++)
    {
	/* either search j = 0:k-1 or j = Rj [0:rnz-1] */
	j = given_row ? (Rj [kk]) : (kk) ;

	if (j < 0 || j >= k)
	{
	    ERROR (CHOLMOD_INVALID, "R invalid") ;
	    return (FALSE) ;
	}

	PRINT2 (("Prune col j = "ID":\n", j)) ;

	lnz = Lnz [j] ;
	dj = Lx [Lp [j]] ;
	ASSERT (Lnz [j] > 0 && Li [Lp [j]] == j) ;

	if (lnz > 1)
	{
	    left = Lp [j] ;
	    pend = left + lnz ;
	    right = pend - 1 ;

	    i = Li [right] ;

	    if (i < k)
	    {
		/* row k is not in column j */
		continue ;
	    }
	    else if (i == k)
	    {
		/* k is the last row index in this column (quick delete) */
		if (do_solve)
		{
		    Xx [j] -= yk [0] * dj * Lx [right] ;
		}
		Lx [right] = 0 ;
	    }
	    else
	    {
		/* binary search for row k in column j */
		PRINT2 (("\nBinary search: lnz "ID" k = "ID"\n", lnz, k)) ;
		while (left < right)
		{
		    middle = (left + right) / 2 ;
		    PRINT2 (("left "ID" right "ID" middle "ID": ["ID" "ID""
			""ID"]\n", left, right, middle,
			Li [left], Li [middle], Li [right])) ;
		    if (k > Li [middle])
		    {
			left = middle + 1 ;
		    }
		    else
		    {
			right = middle ;
		    }
		}
		ASSERT (left >= Lp [j] && left < pend) ;

#ifndef NDEBUG
		/* brute force, linear-time search */
		{
		    Int p3 = Lp [j] ;
		    i = EMPTY ;
		    PRINT2 (("Brute force:\n")) ;
		    for ( ; p3 < pend ; p3++)
		    {
			i = Li [p3] ;
			PRINT2 (("p "ID" ["ID"]\n", p3, i)) ;
			if (i >= k)
			{
			    break ;
			}
		    }
		    if (i == k)
		    {
			ASSERT (k == Li [p3]) ;
			ASSERT (p3 == left) ;
		    }
		}
#endif

		if (k == Li [left])
		{
		    if (do_solve)
		    {
			Xx [j] -= yk [0] * dj * Lx [left] ;
		    }
		    /* found row k in column j.  Prune it from the column.*/
		    Lx [left] = 0 ;
		}
	    }
	}
    }

#ifndef NDEBUG
    /* ensure that row k has been deleted from the matrix L */
    for (j = 0 ; j < k ; j++)
    {
	Int lasti ;
	lasti = EMPTY ;
	p = Lp [j] ;
	pend = p + Lnz [j] ;
	/* look for row k in column j */
	PRINT1 (("Pruned column "ID"\n", j)) ;
	for ( ; p < pend ; p++)
	{
	    i = Li [p] ;
	    PRINT2 ((" "ID"", i)) ;
	    PRINT2 ((" %g\n", Lx [p])) ;
	    ASSERT (IMPLIES (i == k, Lx [p] == 0)) ;
	    ASSERT (i > lasti) ;
	    lasti = i ;
	}
	PRINT1 (("\n")) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* set diagonal and clear column k of L */
    /* ---------------------------------------------------------------------- */

    lnz = Lnz [k] - 1 ;
    ASSERT (Lnz [k] > 0) ;

    /* ---------------------------------------------------------------------- */
    /* update/downdate */
    /* ---------------------------------------------------------------------- */

    /* update or downdate L (k+1:n, k+1:n) with the vector
     * C = L (:,k) * sqrt (abs (D [k]))
     * Do a numeric update if D[k] > 0, numeric downdate otherwise.
     */

    PRINT1 (("rowdel downdate lnz = "ID"\n", lnz)) ;

    /* store the new unit diagonal */
    p = Lp [k] ;
    pend = p + lnz + 1 ;
    dk = Lx [p] ;
    Lx [p++] = 1 ;
    PRINT2 (("D [k = "ID"] = %g\n", k, dk)) ;
    ok = TRUE ;
    fl = 0 ;

    if (lnz > 0)
    {
	/* compute DeltaB for updown (in DeltaB) */
	if (do_solve)
	{
	    xk = Xx [k] - yk [0] * dk ;
	    for ( ; p < pend ; p++)
	    {
		Nx [Li [p]] += Lx [p] * xk ;
	    }
	}

	do_update = IS_GT_ZERO (dk) ;
	if (!do_update)
	{
	    dk = -dk ;
	}
	sqrt_dk = sqrt (dk) ;
	p = Lp [k] + 1 ;
	for (kk = 0 ; kk < lnz ; kk++, p++)
	{
	    Ci [kk] = Li [p] ;
	    Cx [kk] = Lx [p] * sqrt_dk ;
	    Lx [p] = 0 ;		/* clear column k */
	}
	fl = lnz + 1 ;

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

	/* numeric update if dk > 0, and with Lx=b change */
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

    if (do_solve)
    {
	/* kth equation becomes identity, so X(k) is now Y(k) */
	Xx [k] = yk [0] ;
    }

    DEBUG (CHOLMOD(dump_factor) (L, "LDL factorization, L:", Common)) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 2*n, Common)) ;
    return (ok) ;
}
#endif
