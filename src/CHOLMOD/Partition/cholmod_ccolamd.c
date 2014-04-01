/* ========================================================================== */
/* === Partition/cholmod_ccolamd ============================================ */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Partition Module. 
 * Copyright (C) 2005-2013, Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Partition Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* CHOLMOD interface to the CCOLAMD ordering routine.  Finds a permutation
 * p such that the Cholesky factorization of PAA'P' is sparser than AA'.
 * The column etree is found and postordered, and the ccolamd ordering is then
 * combined with its postordering.  A must be unsymmetric.
 *
 * workspace: Iwork (MAX (nrow,ncol))
 *	Allocates a copy of its input matrix, which is
 *	then used as CCOLAMD's workspace.
 *
 * Supports any xtype (pattern, real, complex, or zomplex).
 */

#ifndef NCAMD

#include "cholmod_internal.h"
#include "ccolamd.h"
#include "cholmod_camd.h"

#if (CCOLAMD_VERSION < CCOLAMD_VERSION_CODE (2,5))
#error "CCOLAMD v2.0 or later is required"
#endif

/* ========================================================================== */
/* === ccolamd_interface ==================================================== */
/* ========================================================================== */

/* Order with ccolamd */

static int ccolamd_interface
(
    cholmod_sparse *A,
    size_t alen,
    Int *Perm,
    Int *Cmember,
    Int *fset,
    Int fsize,
    cholmod_sparse *C,
    cholmod_common *Common
)
{
    double knobs [CCOLAMD_KNOBS] ;
    Int *Cp = NULL ;
    Int ok, k, nrow, ncol, stats [CCOLAMD_STATS] ;

    nrow = A->nrow ;
    ncol = A->ncol ;

    /* ---------------------------------------------------------------------- */
    /* copy (and transpose) the input matrix A into the ccolamd workspace */
    /* ---------------------------------------------------------------------- */

    /* C = A (:,f)', which also packs A if needed. */
    /* workspace: Iwork (nrow if no fset; MAX (nrow,ncol) if fset non-NULL) */
    ok = CHOLMOD(transpose_unsym) (A, 0, NULL, fset, fsize, C, Common) ;

    /* ---------------------------------------------------------------------- */
    /* order the matrix (destroys the contents of C->i and C->p) */
    /* ---------------------------------------------------------------------- */

    /* get parameters */
#ifdef LONG
    ccolamd_l_set_defaults (knobs) ;
#else
    ccolamd_set_defaults (knobs) ;
#endif

    if (Common->current < 0 || Common->current >= CHOLMOD_MAXMETHODS)
    {
	/* this is the CHOLMOD default, not the CCOLAMD default */
	knobs [CCOLAMD_DENSE_ROW] = -1 ;
    }
    else
    {
	/* get the knobs from the Common parameters */
	knobs [CCOLAMD_DENSE_COL] =Common->method[Common->current].prune_dense ;
	knobs [CCOLAMD_DENSE_ROW] =Common->method[Common->current].prune_dense2;
	knobs [CCOLAMD_AGGRESSIVE]=Common->method[Common->current].aggressive ;
	knobs [CCOLAMD_LU]        =Common->method[Common->current].order_for_lu;
    }

    if (ok)
    {

#ifdef LONG
	ccolamd_l (ncol, nrow, alen, C->i, C->p, knobs, stats, Cmember) ;
#else
	ccolamd (ncol, nrow, alen, C->i, C->p, knobs, stats, Cmember) ;
#endif

	ok = stats [CCOLAMD_STATUS] ;

	ok = (ok == CCOLAMD_OK || ok == CCOLAMD_OK_BUT_JUMBLED) ;
	/* permutation returned in C->p, if the ordering succeeded */
	Cp = C->p ;
	for (k = 0 ; k < nrow ; k++)
	{
	    Perm [k] = Cp [k] ;
	}
    }

    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_ccolamd ====================================================== */
/* ========================================================================== */

/* Order AA' or A(:,f)*A(:,f)' using CCOLAMD. */

int CHOLMOD(ccolamd)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    Int *Cmember,	/* size A->nrow.  Cmember [i] = c if row i is in the
			 * constraint set c.  c must be >= 0.  The # of
			 * constraint sets is max (Cmember) + 1.  If Cmember is
			 * NULL, then it is interpretted as Cmember [i] = 0 for
			 * all i */
    /* ---- output --- */
    Int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_sparse *C ;
    Int ok, nrow, ncol ;
    size_t alen ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Perm, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    if (A->stype != 0)
    {
	ERROR (CHOLMOD_INVALID, "matrix must be unsymmetric") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    ncol = A->ncol ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

#ifdef LONG
    alen = ccolamd_l_recommended (A->nzmax, ncol, nrow) ;
#else
    alen = ccolamd_recommended (A->nzmax, ncol, nrow) ;
#endif

    if (alen == 0)
    {
	ERROR (CHOLMOD_TOO_LARGE, "matrix invalid or too large") ;
	return (FALSE) ;
    }

    CHOLMOD(allocate_work) (0, MAX (nrow,ncol), 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }

    C = CHOLMOD(allocate_sparse) (ncol, nrow, alen, TRUE, TRUE, 0,
	    CHOLMOD_PATTERN, Common) ;

    /* ---------------------------------------------------------------------- */
    /* order with ccolamd */
    /* ---------------------------------------------------------------------- */

    ok = ccolamd_interface (A, alen, Perm, Cmember, fset, fsize, C, Common) ;

    /* ---------------------------------------------------------------------- */
    /* free the workspace and return result */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(free_sparse) (&C, Common) ;
    return (ok) ;
}
#endif
