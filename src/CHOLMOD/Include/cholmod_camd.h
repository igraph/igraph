/* ========================================================================== */
/* === Include/cholmod_camd.h =============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod_camd.h.
 * Copyright (C) 2005-2013, Univ. of Florida.  Author: Timothy A. Davis
 * CHOLMOD/Include/cholmod_partition.h is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* CHOLMOD Partition module, interface to CAMD, CCOLAMD, and CSYMAMD
 *
 * An interface to CCOLAMD and CSYMAMD, constrained minimum degree ordering
 * methods which order a matrix following constraints determined via nested
 * dissection.
 *
 * These functions do not require METIS.  They are installed unless NCAMD
 * is defined:
 * cholmod_ccolamd		interface to CCOLAMD ordering
 * cholmod_csymamd		interface to CSYMAMD ordering
 * cholmod_camd			interface to CAMD ordering
 *
 * Requires the Core and Cholesky modules, and two packages: CAMD,
 * and CCOLAMD.  Used by functions in the Partition Module.
 */

#ifndef CHOLMOD_CAMD_H
#define CHOLMOD_CAMD_H

#include "cholmod_core.h"

/* -------------------------------------------------------------------------- */
/* cholmod_ccolamd */
/* -------------------------------------------------------------------------- */

/* Order AA' or A(:,f)*A(:,f)' using CCOLAMD. */

int cholmod_ccolamd
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int *Cmember,	/* size A->nrow.  Cmember [i] = c if row i is in the
			 * constraint set c.  c must be >= 0.  The # of
			 * constraint sets is max (Cmember) + 1.  If Cmember is
			 * NULL, then it is interpretted as Cmember [i] = 0 for
			 * all i */
    /* ---- output --- */
    int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_ccolamd (cholmod_sparse *, SuiteSparse_long *, size_t,
    SuiteSparse_long *, SuiteSparse_long *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_csymamd */
/* -------------------------------------------------------------------------- */

/* Order A using CSYMAMD. */

int cholmod_csymamd
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    /* ---- output --- */
    int *Cmember,	/* size nrow.  see cholmod_ccolamd above */
    int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_csymamd (cholmod_sparse *, SuiteSparse_long *,
    SuiteSparse_long *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_camd */
/* -------------------------------------------------------------------------- */

/* Order A using CAMD. */

int cholmod_camd
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* ---- output --- */
    int *Cmember,	/* size nrow.  see cholmod_ccolamd above */
    int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_camd (cholmod_sparse *, SuiteSparse_long *, size_t,
    SuiteSparse_long *, SuiteSparse_long *, cholmod_common *) ;

#endif
