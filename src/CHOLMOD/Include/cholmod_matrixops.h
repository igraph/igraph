/* ========================================================================== */
/* === Include/cholmod_matrixops.h ========================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod_matrixops.h.
 * Copyright (C) 2005-2006, Timothy A. Davis
 * CHOLMOD/Include/cholmod_matrixops.h is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* CHOLMOD MatrixOps module.
 *
 * Basic operations on sparse and dense matrices.
 *
 * cholmod_drop		    A = entries in A with abs. value >= tol
 * cholmod_norm_dense	    s = norm (X), 1-norm, inf-norm, or 2-norm
 * cholmod_norm_sparse	    s = norm (A), 1-norm or inf-norm
 * cholmod_horzcat	    C = [A,B]
 * cholmod_scale	    A = diag(s)*A, A*diag(s), s*A or diag(s)*A*diag(s)
 * cholmod_sdmult	    Y = alpha*(A*X) + beta*Y or alpha*(A'*X) + beta*Y
 * cholmod_ssmult	    C = A*B
 * cholmod_submatrix	    C = A (i,j), where i and j are arbitrary vectors
 * cholmod_vertcat	    C = [A ; B]
 *
 * A, B, C: sparse matrices (cholmod_sparse)
 * X, Y: dense matrices (cholmod_dense)
 * s: scalar or vector
 *
 * Requires the Core module.  Not required by any other CHOLMOD module.
 */

#ifndef CHOLMOD_MATRIXOPS_H
#define CHOLMOD_MATRIXOPS_H

#include "cholmod_core.h"

/* -------------------------------------------------------------------------- */
/* cholmod_drop:  drop entries with small absolute value */
/* -------------------------------------------------------------------------- */

int cholmod_drop
(
    /* ---- input ---- */
    double tol,		/* keep entries with absolute value > tol */
    /* ---- in/out --- */
    cholmod_sparse *A,	/* matrix to drop entries from */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_drop (double, cholmod_sparse *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_norm_dense:  s = norm (X), 1-norm, inf-norm, or 2-norm */
/* -------------------------------------------------------------------------- */

double cholmod_norm_dense
(
    /* ---- input ---- */
    cholmod_dense *X,	/* matrix to compute the norm of */
    int norm,		/* type of norm: 0: inf. norm, 1: 1-norm, 2: 2-norm */
    /* --------------- */
    cholmod_common *Common
) ;

double cholmod_l_norm_dense (cholmod_dense *, int, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_norm_sparse:  s = norm (A), 1-norm or inf-norm */
/* -------------------------------------------------------------------------- */

double cholmod_norm_sparse
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to compute the norm of */
    int norm,		/* type of norm: 0: inf. norm, 1: 1-norm */
    /* --------------- */
    cholmod_common *Common
) ;

double cholmod_l_norm_sparse (cholmod_sparse *, int, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_horzcat:  C = [A,B] */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_horzcat
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* left matrix to concatenate */
    cholmod_sparse *B,	/* right matrix to concatenate */
    int values,		/* if TRUE compute the numerical values of C */
    /* --------------- */
    cholmod_common *Common
) ;

cholmod_sparse *cholmod_l_horzcat (cholmod_sparse *, cholmod_sparse *, int,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_scale:  A = diag(s)*A, A*diag(s), s*A or diag(s)*A*diag(s) */
/* -------------------------------------------------------------------------- */

/* scaling modes, selected by the scale input parameter: */
#define CHOLMOD_SCALAR 0	/* A = s*A */
#define CHOLMOD_ROW 1		/* A = diag(s)*A */
#define CHOLMOD_COL 2		/* A = A*diag(s) */
#define CHOLMOD_SYM 3		/* A = diag(s)*A*diag(s) */

int cholmod_scale
(
    /* ---- input ---- */
    cholmod_dense *S,	/* scale factors (scalar or vector) */
    int scale,		/* type of scaling to compute */
    /* ---- in/out --- */
    cholmod_sparse *A,	/* matrix to scale */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_scale (cholmod_dense *, int, cholmod_sparse *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_sdmult:  Y = alpha*(A*X) + beta*Y or alpha*(A'*X) + beta*Y */
/* -------------------------------------------------------------------------- */

/* Sparse matrix times dense matrix */

int cholmod_sdmult
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* sparse matrix to multiply */
    int transpose,	/* use A if 0, or A' otherwise */
    double alpha [2],   /* scale factor for A */
    double beta [2],    /* scale factor for Y */
    cholmod_dense *X,	/* dense matrix to multiply */
    /* ---- in/out --- */
    cholmod_dense *Y,	/* resulting dense matrix */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_sdmult (cholmod_sparse *, int, double *, double *,
    cholmod_dense *, cholmod_dense *Y, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_ssmult:  C = A*B */
/* -------------------------------------------------------------------------- */

/* Sparse matrix times sparse matrix */

cholmod_sparse *cholmod_ssmult
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* left matrix to multiply */
    cholmod_sparse *B,	/* right matrix to multiply */
    int stype,		/* requested stype of C */
    int values,		/* TRUE: do numerical values, FALSE: pattern only */
    int sorted,		/* if TRUE then return C with sorted columns */
    /* --------------- */
    cholmod_common *Common
) ;

cholmod_sparse *cholmod_l_ssmult (cholmod_sparse *, cholmod_sparse *, int, int,
    int, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_submatrix:  C = A (r,c), where i and j are arbitrary vectors */
/* -------------------------------------------------------------------------- */

/* rsize < 0 denotes ":" in MATLAB notation, or more precisely 0:(A->nrow)-1.
 * In this case, r can be NULL.  An rsize of zero, or r = NULL and rsize >= 0,
 * denotes "[ ]" in MATLAB notation (the empty set).
 * Similar rules hold for csize.
 */

cholmod_sparse *cholmod_submatrix
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to subreference */
    int *rset,		/* set of row indices, duplicates OK */
    SuiteSparse_long rsize,	/* size of r; rsize < 0 denotes ":" */
    int *cset,		/* set of column indices, duplicates OK */
    SuiteSparse_long csize,	/* size of c; csize < 0 denotes ":" */
    int values,		/* if TRUE compute the numerical values of C */
    int sorted,		/* if TRUE then return C with sorted columns */
    /* --------------- */
    cholmod_common *Common
) ;

cholmod_sparse *cholmod_l_submatrix (cholmod_sparse *, SuiteSparse_long *,
    SuiteSparse_long, SuiteSparse_long *, SuiteSparse_long, int, int,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_vertcat:  C = [A ; B] */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_vertcat
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* left matrix to concatenate */
    cholmod_sparse *B,	/* right matrix to concatenate */
    int values,		/* if TRUE compute the numerical values of C */
    /* --------------- */
    cholmod_common *Common
) ;

cholmod_sparse *cholmod_l_vertcat (cholmod_sparse *, cholmod_sparse *, int,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_symmetry: determine if a sparse matrix is symmetric */
/* -------------------------------------------------------------------------- */

int cholmod_symmetry
(
    /* ---- input ---- */
    cholmod_sparse *A,
    int option,
    /* ---- output ---- */
    int *xmatched,
    int *pmatched,
    int *nzoffdiag,
    int *nzdiag,
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_symmetry (cholmod_sparse *, int, SuiteSparse_long *,
    SuiteSparse_long *, SuiteSparse_long *, SuiteSparse_long *,
    cholmod_common *) ;

#endif
