/* ========================================================================== */
/* === Include/cholmod_modify.h ============================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod_modify.h.
 * Copyright (C) 2005-2006, Timothy A. Davis and William W. Hager
 * CHOLMOD/Include/cholmod_modify.h is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* CHOLMOD Modify module.
 *
 * Sparse Cholesky modification routines: update / downdate / rowadd / rowdel.
 * Can also modify a corresponding solution to Lx=b when L is modified.  This
 * module is most useful when applied on a Cholesky factorization computed by
 * the Cholesky module, but it does not actually require the Cholesky module.
 * The Core module can create an identity Cholesky factorization (LDL' where
 * L=D=I) that can then by modified by these routines.
 *
 * Primary routines:
 * -----------------
 *
 * cholmod_updown	    multiple rank update/downdate
 * cholmod_rowadd	    add a row to an LDL' factorization
 * cholmod_rowdel	    delete a row from an LDL' factorization
 *
 * Secondary routines:
 * -------------------
 *
 * cholmod_updown_solve	    update/downdate, and modify solution to Lx=b
 * cholmod_updown_mark	    update/downdate, and modify solution to partial Lx=b
 * cholmod_updown_mask	    update/downdate for LPDASA
 * cholmod_rowadd_solve	    add a row, and update solution to Lx=b
 * cholmod_rowadd_mark	    add a row, and update solution to partial Lx=b
 * cholmod_rowdel_solve	    delete a row, and downdate Lx=b
 * cholmod_rowdel_mark	    delete a row, and downdate solution to partial Lx=b
 *
 * Requires the Core module.  Not required by any other CHOLMOD module.
 */

#ifndef CHOLMOD_MODIFY_H
#define CHOLMOD_MODIFY_H

#include "cholmod_core.h"

/* -------------------------------------------------------------------------- */
/* cholmod_updown:  multiple rank update/downdate */
/* -------------------------------------------------------------------------- */

/* Compute the new LDL' factorization of LDL'+CC' (an update) or LDL'-CC'
 * (a downdate).  The factor object L need not be an LDL' factorization; it
 * is converted to one if it isn't. */

int cholmod_updown 
(
    /* ---- input ---- */
    int update,		/* TRUE for update, FALSE for downdate */
    cholmod_sparse *C,	/* the incoming sparse update */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_updown (int, cholmod_sparse *, cholmod_factor *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_updown_solve:  update/downdate, and modify solution to Lx=b */
/* -------------------------------------------------------------------------- */

/* Does the same as cholmod_updown, except that it also updates/downdates the
 * solution to Lx=b+DeltaB.  x and b must be n-by-1 dense matrices.  b is not
 * need as input to this routine, but a sparse change to b is (DeltaB).  Only
 * entries in DeltaB corresponding to columns modified in L are accessed; the
 * rest must be zero. */

int cholmod_updown_solve
(
    /* ---- input ---- */
    int update,		/* TRUE for update, FALSE for downdate */
    cholmod_sparse *C,	/* the incoming sparse update */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    cholmod_dense *X,	/* solution to Lx=b (size n-by-1) */
    cholmod_dense *DeltaB,  /* change in b, zero on output */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_updown_solve (int, cholmod_sparse *, cholmod_factor *,
    cholmod_dense *, cholmod_dense *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_updown_mark:  update/downdate, and modify solution to partial Lx=b */
/* -------------------------------------------------------------------------- */

/* Does the same as cholmod_updown_solve, except only part of L is used in
 * the update/downdate of the solution to Lx=b.  This routine is an "expert"
 * routine.  It is meant for use in LPDASA only.  See cholmod_updown.c for
 * a description of colmark. */

int cholmod_updown_mark
(
    /* ---- input ---- */
    int update,		/* TRUE for update, FALSE for downdate */
    cholmod_sparse *C,	/* the incoming sparse update */
    int *colmark,	/* int array of size n.  See cholmod_updown.c */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    cholmod_dense *X,	/* solution to Lx=b (size n-by-1) */
    cholmod_dense *DeltaB,  /* change in b, zero on output */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_updown_mark (int, cholmod_sparse *, SuiteSparse_long *,
    cholmod_factor *, cholmod_dense *, cholmod_dense *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_updown_mask:  update/downdate, for LPDASA */
/* -------------------------------------------------------------------------- */

/* Does the same as cholmod_updown_mark, except has an additional "mask"
 * argument.  This routine is an "expert" routine.  It is meant for use in
 * LPDASA only.  See cholmod_updown.c for a description of mask. */

int cholmod_updown_mask
(
    /* ---- input ---- */
    int update,		/* TRUE for update, FALSE for downdate */
    cholmod_sparse *C,	/* the incoming sparse update */
    int *colmark,	/* int array of size n.  See cholmod_updown.c */
    int *mask,		/* size n */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    cholmod_dense *X,	/* solution to Lx=b (size n-by-1) */
    cholmod_dense *DeltaB,  /* change in b, zero on output */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_updown_mask (int, cholmod_sparse *, SuiteSparse_long *,
    SuiteSparse_long *, cholmod_factor *, cholmod_dense *, cholmod_dense *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_rowadd:  add a row to an LDL' factorization (a rank-2 update) */
/* -------------------------------------------------------------------------- */

/* cholmod_rowadd adds a row to the LDL' factorization.  It computes the kth
 * row and kth column of L, and then updates the submatrix L (k+1:n,k+1:n)
 * accordingly.  The kth row and column of L must originally be equal to the
 * kth row and column of the identity matrix.  The kth row/column of L is
 * computed as the factorization of the kth row/column of the matrix to
 * factorize, which is provided as a single n-by-1 sparse matrix R. */

int cholmod_rowadd 
(
    /* ---- input ---- */
    size_t k,		/* row/column index to add */
    cholmod_sparse *R,	/* row/column of matrix to factorize (n-by-1) */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_rowadd (size_t, cholmod_sparse *, cholmod_factor *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_rowadd_solve:  add a row, and update solution to Lx=b */
/* -------------------------------------------------------------------------- */

/* Does the same as cholmod_rowadd, and also updates the solution to Lx=b
 * See cholmod_updown for a description of how Lx=b is updated.  There is on
 * additional parameter:  bk specifies the new kth entry of b. */

int cholmod_rowadd_solve
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
) ;

int cholmod_l_rowadd_solve (size_t, cholmod_sparse *, double *,
    cholmod_factor *, cholmod_dense *, cholmod_dense *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_rowadd_mark:  add a row, and update solution to partial Lx=b */
/* -------------------------------------------------------------------------- */

/* Does the same as cholmod_rowadd_solve, except only part of L is used in
 * the update/downdate of the solution to Lx=b.  This routine is an "expert"
 * routine.  It is meant for use in LPDASA only.  */

int cholmod_rowadd_mark
(
    /* ---- input ---- */
    size_t k,		/* row/column index to add */
    cholmod_sparse *R,	/* row/column of matrix to factorize (n-by-1) */
    double bk [2],	/* kth entry of the right hand side, b */
    int *colmark,	/* int array of size n.  See cholmod_updown.c */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    cholmod_dense *X,	/* solution to Lx=b (size n-by-1) */
    cholmod_dense *DeltaB,  /* change in b, zero on output */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_rowadd_mark (size_t, cholmod_sparse *, double *,
    SuiteSparse_long *, cholmod_factor *, cholmod_dense *, cholmod_dense *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_rowdel:  delete a row from an LDL' factorization (a rank-2 update) */
/* -------------------------------------------------------------------------- */

/* Sets the kth row and column of L to be the kth row and column of the identity
 * matrix, and updates L(k+1:n,k+1:n) accordingly.   To reduce the running time,
 * the caller can optionally provide the nonzero pattern (or an upper bound) of
 * kth row of L, as the sparse n-by-1 vector R.  Provide R as NULL if you want
 * CHOLMOD to determine this itself, which is easier for the caller, but takes
 * a little more time.
 */

int cholmod_rowdel 
(
    /* ---- input ---- */
    size_t k,		/* row/column index to delete */
    cholmod_sparse *R,	/* NULL, or the nonzero pattern of kth row of L */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_rowdel (size_t, cholmod_sparse *, cholmod_factor *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_rowdel_solve:  delete a row, and downdate Lx=b */
/* -------------------------------------------------------------------------- */

/* Does the same as cholmod_rowdel, but also downdates the solution to Lx=b.
 * When row/column k of A is "deleted" from the system A*y=b, this can induce
 * a change to x, in addition to changes arising when L and b are modified.
 * If this is the case, the kth entry of y is required as input (yk) */

int cholmod_rowdel_solve
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
) ;

int cholmod_l_rowdel_solve (size_t, cholmod_sparse *, double *,
    cholmod_factor *, cholmod_dense *, cholmod_dense *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_rowdel_mark:  delete a row, and downdate solution to partial Lx=b */
/* -------------------------------------------------------------------------- */

/* Does the same as cholmod_rowdel_solve, except only part of L is used in
 * the update/downdate of the solution to Lx=b.  This routine is an "expert"
 * routine.  It is meant for use in LPDASA only.  */

int cholmod_rowdel_mark
(
    /* ---- input ---- */
    size_t k,		/* row/column index to delete */
    cholmod_sparse *R,	/* NULL, or the nonzero pattern of kth row of L */
    double yk [2],	/* kth entry in the solution to A*y=b */
    int *colmark,	/* int array of size n.  See cholmod_updown.c */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    cholmod_dense *X,	/* solution to Lx=b (size n-by-1) */
    cholmod_dense *DeltaB,  /* change in b, zero on output */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_rowdel_mark (size_t, cholmod_sparse *, double *,
    SuiteSparse_long *, cholmod_factor *, cholmod_dense *, cholmod_dense *,
    cholmod_common *) ;

#endif
