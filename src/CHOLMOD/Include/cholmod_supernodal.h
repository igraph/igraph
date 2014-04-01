/* ========================================================================== */
/* === Include/cholmod_supernodal.h ========================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod_supernodal.h.
 * Copyright (C) 2005-2006, Timothy A. Davis
 * CHOLMOD/Include/cholmod_supernodal.h is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* CHOLMOD Supernodal module.
 *
 * Supernodal analysis, factorization, and solve.  The simplest way to use
 * these routines is via the Cholesky module.  It does not provide any
 * fill-reducing orderings, but does accept the orderings computed by the
 * Cholesky module.  It does not require the Cholesky module itself, however.
 *
 * Primary routines:
 * -----------------
 * cholmod_super_symbolic	supernodal symbolic analysis
 * cholmod_super_numeric	supernodal numeric factorization
 * cholmod_super_lsolve		supernodal Lx=b solve
 * cholmod_super_ltsolve	supernodal L'x=b solve
 *
 * Prototypes for the BLAS and LAPACK routines that CHOLMOD uses are listed
 * below, including how they are used in CHOLMOD.
 *
 * BLAS routines:
 * --------------
 * dtrsv	solve Lx=b or L'x=b, L non-unit diagonal, x and b stride-1
 * dtrsm	solve LX=B or L'X=b, L non-unit diagonal
 * dgemv	y=y-A*x or y=y-A'*x (x and y stride-1)
 * dgemm	C=A*B', C=C-A*B, or C=C-A'*B
 * dsyrk	C=tril(A*A')
 *
 * LAPACK routines:
 * ----------------
 * dpotrf	LAPACK: A=chol(tril(A))
 *
 * Requires the Core module, and two external packages: LAPACK and the BLAS.
 * Optionally used by the Cholesky module.
 */

#ifndef CHOLMOD_SUPERNODAL_H
#define CHOLMOD_SUPERNODAL_H

#include "cholmod_core.h"

/* -------------------------------------------------------------------------- */
/* cholmod_super_symbolic */
/* -------------------------------------------------------------------------- */

/* Analyzes A, AA', or A(:,f)*A(:,f)' in preparation for a supernodal numeric
 * factorization.  The user need not call this directly; cholmod_analyze is
 * a "simple" wrapper for this routine.
 */

int cholmod_super_symbolic
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    cholmod_sparse *F,	/* F = A' or A(:,f)' */
    int *Parent,	/* elimination tree */
    /* ---- in/out --- */
    cholmod_factor *L,	/* simplicial symbolic on input,
			 * supernodal symbolic on output */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_super_symbolic (cholmod_sparse *, cholmod_sparse *,
    SuiteSparse_long *, cholmod_factor *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_super_symbolic2 */
/* -------------------------------------------------------------------------- */

/* Analyze for supernodal Cholesky or multifrontal QR.  CHOLMOD itself always
 * analyzes for supernodal Cholesky, of course.  This "for_cholesky = TRUE"
 * option is used by SuiteSparseQR only.   Added for V1.7 */

int cholmod_super_symbolic2
(
    /* ---- input ---- */
    int for_cholesky,   /* Cholesky if TRUE, QR if FALSE */
    cholmod_sparse *A,	/* matrix to analyze */
    cholmod_sparse *F,	/* F = A' or A(:,f)' */
    int *Parent,	/* elimination tree */
    /* ---- in/out --- */
    cholmod_factor *L,	/* simplicial symbolic on input,
			 * supernodal symbolic on output */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_super_symbolic2 (int, cholmod_sparse *, cholmod_sparse *,
    SuiteSparse_long *, cholmod_factor *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_super_numeric */
/* -------------------------------------------------------------------------- */

/* Computes the numeric LL' factorization of A, AA', or A(:,f)*A(:,f)' using
 * a BLAS-based supernodal method.  The user need not call this directly;
 * cholmod_factorize is a "simple" wrapper for this routine.
 */

int cholmod_super_numeric
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to factorize */
    cholmod_sparse *F,	/* F = A' or A(:,f)' */
    double beta [2],	/* beta*I is added to diagonal of matrix to factorize */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factorization */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_super_numeric (cholmod_sparse *, cholmod_sparse *, double *,
    cholmod_factor *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_super_lsolve */
/* -------------------------------------------------------------------------- */

/* Solve Lx=b where L is from a supernodal numeric factorization.  The user
 * need not call this routine directly.  cholmod_solve is a "simple" wrapper
 * for this routine. */

int cholmod_super_lsolve
(
    /* ---- input ---- */
    cholmod_factor *L,	/* factor to use for the forward solve */
    /* ---- output ---- */
    cholmod_dense *X,	/* b on input, solution to Lx=b on output */
    /* ---- workspace   */
    cholmod_dense *E,	/* workspace of size nrhs*(L->maxesize) */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_super_lsolve (cholmod_factor *, cholmod_dense *, cholmod_dense *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_super_ltsolve */
/* -------------------------------------------------------------------------- */

/* Solve L'x=b where L is from a supernodal numeric factorization.  The user
 * need not call this routine directly.  cholmod_solve is a "simple" wrapper
 * for this routine. */

int cholmod_super_ltsolve
(
    /* ---- input ---- */
    cholmod_factor *L,	/* factor to use for the backsolve */
    /* ---- output ---- */
    cholmod_dense *X,	/* b on input, solution to L'x=b on output */
    /* ---- workspace   */
    cholmod_dense *E,	/* workspace of size nrhs*(L->maxesize) */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_super_ltsolve (cholmod_factor *, cholmod_dense *, cholmod_dense *,
    cholmod_common *) ;

#endif
