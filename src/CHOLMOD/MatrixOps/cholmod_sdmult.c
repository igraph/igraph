/* ========================================================================== */
/* === MatrixOps/cholmod_sdmult ============================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/MatrixOps Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/MatrixOps Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* Sparse matrix times dense matrix:
 * Y = alpha*(A*X) + beta*Y or Y = alpha*(A'*X) + beta*Y,
 * where A is sparse and X and Y are dense.
 *
 * when using A,  X has A->ncol columns and Y has A->nrow rows
 * when using A', X has A->nrow columns and Y has A->ncol rows
 *
 * workspace: none in Common.  Temporary workspace of size 4*(X->nrow) is used
 * if A is stored in symmetric form and X has four columns or more.  If the
 * workspace is not available, a slower method is used instead that requires
 * no workspace.
 *
 * transpose = 0: use A
 * otherwise, use A'  (complex conjugate transpose)
 *
 * transpose is ignored if the matrix is symmetric or Hermitian.
 * (the array transpose A.' is not supported).
 *
 * Supports real, complex, and zomplex matrices, but the xtypes of A, X, and Y
 * must all match.
 */

#ifndef NMATRIXOPS

#include "cholmod_internal.h"
#include "cholmod_matrixops.h"


/* ========================================================================== */
/* === TEMPLATE ============================================================= */
/* ========================================================================== */

#define REAL
#include "t_cholmod_sdmult.c"
#define COMPLEX
#include "t_cholmod_sdmult.c"
#define ZOMPLEX
#include "t_cholmod_sdmult.c"

/* ========================================================================== */
/* === cholmod_sdmult ======================================================= */
/* ========================================================================== */

int CHOLMOD(sdmult)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* sparse matrix to multiply */
    int transpose,	/* use A if 0, otherwise use A' */
    double alpha [2],   /* scale factor for A */
    double beta [2],    /* scale factor for Y */
    cholmod_dense *X,	/* dense matrix to multiply */
    /* ---- in/out --- */
    cholmod_dense *Y,	/* resulting dense matrix */
    /* --------------- */
    cholmod_common *Common
)
{
    double *w ;
    size_t nx, ny ;
    Int e ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (X, FALSE) ;
    RETURN_IF_NULL (Y, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (X, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (Y, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    ny = transpose ? A->ncol : A->nrow ;	/* required length of Y */
    nx = transpose ? A->nrow : A->ncol ;	/* required length of X */
    if (X->nrow != nx || X->ncol != Y->ncol || Y->nrow != ny)
    {
	/* X and/or Y have the wrong dimension */
	ERROR (CHOLMOD_INVALID, "X and/or Y have wrong dimensions") ;
	return (FALSE) ;
    }
    if (A->xtype != X->xtype || A->xtype != Y->xtype)
    {
	ERROR (CHOLMOD_INVALID, "A, X, and Y must have same xtype") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace, if required */
    /* ---------------------------------------------------------------------- */

    w = NULL ;
    e = (A->xtype == CHOLMOD_REAL ? 1:2) ;
    if (A->stype && X->ncol >= 4)
    {
	w = CHOLMOD(malloc) (nx, 4*e*sizeof (double), Common) ;
    }
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* Y = alpha*op(A)*X + beta*Y via template routine */
    /* ---------------------------------------------------------------------- */

    ASSERT (CHOLMOD(dump_sparse) (A, "A", Common) >= 0) ;
    DEBUG (CHOLMOD(dump_dense) (X, "X", Common)) ;
    DEBUG (if (IS_NONZERO (beta [0])
	   || (IS_NONZERO (beta [1]) && A->xtype != CHOLMOD_REAL))
	    CHOLMOD(dump_dense) (Y, "Y", Common)) ;

    switch (A->xtype)
    {

	case CHOLMOD_REAL:
	    r_cholmod_sdmult (A, transpose, alpha, beta, X, Y, w) ;
	    break ;

	case CHOLMOD_COMPLEX:
	    c_cholmod_sdmult (A, transpose, alpha, beta, X, Y, w) ;
	    break ;

	case CHOLMOD_ZOMPLEX:
	    z_cholmod_sdmult (A, transpose, alpha, beta, X, Y, w) ;
	    break ;
    }

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(free) (4*nx, e*sizeof (double), w, Common) ;
    DEBUG (CHOLMOD(dump_dense) (Y, "Y", Common)) ;
    return (TRUE) ;
}
#endif
