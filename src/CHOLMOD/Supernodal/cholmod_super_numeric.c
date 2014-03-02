/* ========================================================================== */
/* === Supernodal/cholmod_super_numeric ===================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Supernodal Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/Supernodal Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* Computes the Cholesky factorization of A+beta*I or A*F+beta*I.  Only the
 * the lower triangular part of A+beta*I or A*F+beta*I is accessed.  The
 * matrices A and F must already be permuted according to the fill-reduction
 * permutation L->Perm.  cholmod_factorize is an "easy" wrapper for this code
 * which applies that permutation.  beta is real.
 *
 * Symmetric case: A is a symmetric (lower) matrix.  F is not accessed.
 * With a fill-reducing permutation, A(p,p) should be passed instead, where is
 * p is L->Perm.
 *
 * Unsymmetric case: A is unsymmetric, and F must be present.  Normally, F=A'.
 * With a fill-reducing permutation, A(p,f) and A(p,f)' should be passed as A
 * and F, respectively, where f is a list of the subset of the columns of A.
 *
 * The input factorization L must be supernodal (L->is_super is TRUE).  It can
 * either be symbolic or numeric.  In the first case, L has been analyzed by
 * cholmod_analyze or cholmod_super_symbolic, but the matrix has not yet been
 * numerically factorized.  The numerical values are allocated here and the
 * factorization is computed.  In the second case, a prior matrix has been
 * analyzed and numerically factorized, and a new matrix is being factorized.
 * The numerical values of L are replaced with the new numerical factorization.
 *
 * L->is_ll is ignored, and set to TRUE.  This routine always computes an LL'
 * factorization.  Supernodal LDL' factorization is not (yet) supported.
 * FUTURE WORK: perform a supernodal LDL' factorization if L->is_ll is FALSE.
 *
 * Uses BLAS routines dsyrk, dgemm, dtrsm, and the LAPACK routine dpotrf.
 * The supernodal solver uses BLAS routines dtrsv, dgemv, dtrsm, and dgemm.
 *
 * If the matrix is not positive definite the routine returns TRUE, but sets
 * Common->status to CHOLMOD_NOT_POSDEF and L->minor is set to the column at
 * which the failure occurred.  The supernode containing the non-positive
 * diagonal entry is set to zero (this includes columns to the left of L->minor
 * in the same supernode), as are all subsequent supernodes.
 *
 * workspace: Flag (nrow), Head (nrow+1), Iwork (2*nrow + 4*nsuper).
 *	Allocates temporary space of size L->maxcsize * sizeof(double)
 *	(twice that for the complex/zomplex case).
 *
 * If L is supernodal symbolic on input, it is converted to a supernodal numeric
 * factor on output, with an xtype of real if A is real, or complex if A is
 * complex or zomplex.  If L is supernodal numeric on input, its xtype must
 * match A (except that L can be complex and A zomplex).  The xtype of A and F
 * must match.
 */

#ifndef NSUPERNODAL

#include "cholmod_internal.h"
#include "cholmod_supernodal.h"
#include "igraph_blas_internal.h"
#include "igraph_lapack_internal.h"

/* ========================================================================== */
/* === TEMPLATE codes for GPU and regular numeric factorization ============= */
/* ========================================================================== */

#ifdef GPU_BLAS
#define REAL
#include "t_cholmod_gpu.c"
#define COMPLEX
#include "t_cholmod_gpu.c"
#define ZOMPLEX
#include "t_cholmod_gpu.c"
#endif

#define REAL
#include "t_cholmod_super_numeric.c"
/* #define COMPLEX */
/* #include "t_cholmod_super_numeric.c" */
/* #define ZOMPLEX */
/* #include "t_cholmod_super_numeric.c" */

/* ========================================================================== */
/* === cholmod_super_numeric ================================================ */
/* ========================================================================== */

/* Returns TRUE if successful, or if the matrix is not positive definite.
 * Returns FALSE if out of memory, inputs are invalid, or other fatal error
 * occurs.
 */

int CHOLMOD(super_numeric)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to factorize */
    cholmod_sparse *F,	/* F = A' or A(:,f)' */
    double beta [2],	/* beta*I is added to diagonal of matrix to factorize */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factorization */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_dense *C ;
    Int *Super, *Map, *SuperMap ;
    size_t maxcsize ;
    Int nsuper, n, i, k, s, stype, nrow ;
    int ok = TRUE, symbolic ;
    size_t t, w ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_PATTERN, CHOLMOD_COMPLEX, FALSE) ;
    stype = A->stype ;
    if (stype < 0)
    {
	if (A->nrow != A->ncol || A->nrow != L->n)
	{
	    ERROR (CHOLMOD_INVALID, "invalid dimensions") ;
	    return (FALSE) ;
	}
    }
    else if (stype == 0)
    {
	if (A->nrow != L->n)
	{
	    ERROR (CHOLMOD_INVALID, "invalid dimensions") ;
	    return (FALSE) ;
	}
	RETURN_IF_NULL (F, FALSE) ;
	RETURN_IF_XTYPE_INVALID (F, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
	if (A->nrow != F->ncol || A->ncol != F->nrow || F->stype != 0)
	{
	    ERROR (CHOLMOD_INVALID, "F invalid") ;
	    return (FALSE) ;
	}
	if (A->xtype != F->xtype)
	{
	    ERROR (CHOLMOD_INVALID, "A and F must have same xtype") ;
	    return (FALSE) ;
	}
    }
    else
    {
	/* symmetric upper case not suppored */
	ERROR (CHOLMOD_INVALID, "symmetric upper case not supported") ;
	return (FALSE) ;
    }
    if (!(L->is_super))
    {
	ERROR (CHOLMOD_INVALID, "L not supernodal") ;
	return (FALSE) ;
    }
    if (L->xtype != CHOLMOD_PATTERN)
    {
	if (! ((A->xtype == CHOLMOD_REAL    && L->xtype == CHOLMOD_REAL)
	    || (A->xtype == CHOLMOD_COMPLEX && L->xtype == CHOLMOD_COMPLEX)
	    || (A->xtype == CHOLMOD_ZOMPLEX && L->xtype == CHOLMOD_COMPLEX)))
	{
	    ERROR (CHOLMOD_INVALID, "complex type mismatch") ;
	    return (FALSE) ;
	}
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace in Common */
    /* ---------------------------------------------------------------------- */

    nsuper = L->nsuper ;
    maxcsize = L->maxcsize ;
    nrow = A->nrow ;
    n = nrow ;

    PRINT1 (("nsuper "ID" maxcsize %g\n", nsuper, (double) maxcsize)) ;
    ASSERT (nsuper >= 0 && maxcsize > 0) ;

    /* w = 2*n + 4*nsuper */
    w = CHOLMOD(mult_size_t) (n, 2, &ok) ;
    t = CHOLMOD(mult_size_t) (nsuper, 4, &ok) ;
    w = CHOLMOD(add_size_t) (w, t, &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    CHOLMOD(allocate_work) (n, w, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get the current factor L and allocate numerical part, if needed */
    /* ---------------------------------------------------------------------- */

    Super = L->super ;
    symbolic = (L->xtype == CHOLMOD_PATTERN) ;
    if (symbolic)
    {
	/* convert to supernodal numeric by allocating L->x */
	CHOLMOD(change_factor) (
		(A->xtype == CHOLMOD_REAL) ? CHOLMOD_REAL : CHOLMOD_COMPLEX,
		TRUE, TRUE, TRUE, TRUE, L, Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* the factor L remains in symbolic supernodal form */
	    return (FALSE) ;
	}
    }
    ASSERT (L->dtype == DTYPE) ;
    ASSERT (L->xtype == CHOLMOD_REAL || L->xtype == CHOLMOD_COMPLEX) ;

    /* supernodal LDL' is not supported */
    L->is_ll = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* get more workspace */
    /* ---------------------------------------------------------------------- */

    C = CHOLMOD(allocate_dense) (maxcsize, 1, maxcsize, L->xtype, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	int status = Common->status ;
	if (symbolic)
	{
	    /* Change L back to symbolic, since the numeric values are not
	     * initialized.  This cannot fail. */
	    CHOLMOD(change_factor) (CHOLMOD_PATTERN, TRUE, TRUE, TRUE, TRUE,
		    L, Common) ;
	}
	/* the factor L is now back to the form it had on input */
	Common->status = status ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    SuperMap = Common->Iwork ;		/* size n (i/i/l) */
    Map = Common->Flag ;    /* size n, use Flag as workspace for Map array */
    for (i = 0 ; i < n ; i++)
    {
	Map [i] = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* find the mapping of nodes to relaxed supernodes */
    /* ---------------------------------------------------------------------- */

    /* SuperMap [k] = s if column k is contained in supernode s */
    for (s = 0 ; s < nsuper ; s++)
    {
	PRINT1 (("Super ["ID"] "ID" ncols "ID"\n",
		    s, Super[s], Super[s+1]-Super[s]));
	for (k = Super [s] ; k < Super [s+1] ; k++)
	{
	    SuperMap [k] = s ;
	    PRINT2 (("relaxed SuperMap ["ID"] = "ID"\n", k, SuperMap [k])) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* supernodal numerical factorization, using template routine */
    /* ---------------------------------------------------------------------- */

    switch (A->xtype)
    {
	case CHOLMOD_REAL:
	    ok = r_cholmod_super_numeric (A, F, beta, L, C, Common) ;
	    break ;

	/* case CHOLMOD_COMPLEX: */
	/*     ok = c_cholmod_super_numeric (A, F, beta, L, C, Common) ; */
	/*     break ; */

	/* case CHOLMOD_ZOMPLEX: */
	/*     /\* This operates on complex L, not zomplex *\/ */
	/*     ok = z_cholmod_super_numeric (A, F, beta, L, C, Common) ; */
	/*     break ; */
    }

    /* ---------------------------------------------------------------------- */
    /* clear Common workspace, free temp workspace C, and return */
    /* ---------------------------------------------------------------------- */

    /* Flag array was used as workspace, clear it */
    Common->mark = EMPTY ;
    /* CHOLMOD(clear_flag) (Common) ; */
    CHOLMOD_CLEAR_FLAG (Common) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;
    CHOLMOD(free_dense) (&C, Common) ;
    return (ok) ;
}
#endif
