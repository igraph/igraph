/* ========================================================================== */
/* === Include/cholmod_check.h ============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod_check.h.  Copyright (C) 2005-2006, Timothy A. Davis
 * CHOLMOD/Include/cholmod_check.h is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* CHOLMOD Check module.
 *
 * Routines that check and print the 5 basic data types in CHOLMOD, and 3 kinds
 * of integer vectors (subset, perm, and parent), and read in matrices from a
 * file:
 *
 * cholmod_check_common	    check/print the Common object
 * cholmod_print_common
 *
 * cholmod_check_sparse	    check/print a sparse matrix in column-oriented form
 * cholmod_print_sparse
 *
 * cholmod_check_dense	    check/print a dense matrix
 * cholmod_print_dense
 *
 * cholmod_check_factor	    check/print a Cholesky factorization
 * cholmod_print_factor
 *
 * cholmod_check_triplet    check/print a sparse matrix in triplet form
 * cholmod_print_triplet
 *
 * cholmod_check_subset	    check/print a subset (integer vector in given range)
 * cholmod_print_subset
 *
 * cholmod_check_perm	    check/print a permutation (an integer vector)
 * cholmod_print_perm
 *
 * cholmod_check_parent	    check/print an elimination tree (an integer vector)
 * cholmod_print_parent
 *
 * cholmod_read_triplet	    read a matrix in triplet form (any Matrix Market
 *			    "coordinate" format, or a generic triplet format).
 *
 * cholmod_read_sparse	    read a matrix in sparse form (same file format as
 *			    cholmod_read_triplet).
 *
 * cholmod_read_dense	    read a dense matrix (any Matrix Market "array"
 *			    format, or a generic dense format).
 *
 * cholmod_write_sparse	    write a sparse matrix to a Matrix Market file.
 *
 * cholmod_write_dense	    write a dense matrix to a Matrix Market file.
 *
 * cholmod_print_common and cholmod_check_common are the only two routines that
 * you may call after calling cholmod_finish.
 *
 * Requires the Core module.  Not required by any CHOLMOD module, except when
 * debugging is enabled (in which case all modules require the Check module).
 *
 * See cholmod_read.c for a description of the file formats supported by the
 * cholmod_read_* routines.
 */

#ifndef CHOLMOD_CHECK_H
#define CHOLMOD_CHECK_H

#include "cholmod_core.h"
#include <stdio.h>

/* -------------------------------------------------------------------------- */
/* cholmod_check_common:  check the Common object */
/* -------------------------------------------------------------------------- */

int cholmod_check_common
(
    cholmod_common *Common
) ;

int cholmod_l_check_common (cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_common:  print the Common object */
/* -------------------------------------------------------------------------- */

int cholmod_print_common
(
    /* ---- input ---- */
    const char *name,	/* printed name of Common object */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_common (const char *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_gpu_stats:  print the GPU / CPU statistics */
/* -------------------------------------------------------------------------- */

int cholmod_gpu_stats   (cholmod_common *) ;
int cholmod_l_gpu_stats (cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_sparse:  check a sparse matrix */
/* -------------------------------------------------------------------------- */

int cholmod_check_sparse
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* sparse matrix to check */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_sparse (cholmod_sparse *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_sparse */
/* -------------------------------------------------------------------------- */

int cholmod_print_sparse
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* sparse matrix to print */
    const char *name,	/* printed name of sparse matrix */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_sparse (cholmod_sparse *, const char *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_dense:  check a dense matrix */
/* -------------------------------------------------------------------------- */

int cholmod_check_dense
(
    /* ---- input ---- */
    cholmod_dense *X,	/* dense matrix to check */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_dense (cholmod_dense *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_dense:  print a dense matrix */
/* -------------------------------------------------------------------------- */

int cholmod_print_dense
(
    /* ---- input ---- */
    cholmod_dense *X,	/* dense matrix to print */
    const char *name,	/* printed name of dense matrix */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_dense (cholmod_dense *, const char *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_factor:  check a factor */
/* -------------------------------------------------------------------------- */

int cholmod_check_factor
(
    /* ---- input ---- */
    cholmod_factor *L,	/* factor to check */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_factor (cholmod_factor *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_factor:  print a factor */
/* -------------------------------------------------------------------------- */

int cholmod_print_factor
(
    /* ---- input ---- */
    cholmod_factor *L,	/* factor to print */
    const char *name,	/* printed name of factor */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_factor (cholmod_factor *, const char *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_triplet:  check a sparse matrix in triplet form */
/* -------------------------------------------------------------------------- */

int cholmod_check_triplet
(
    /* ---- input ---- */
    cholmod_triplet *T,	/* triplet matrix to check */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_triplet (cholmod_triplet *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_triplet:  print a triplet matrix */
/* -------------------------------------------------------------------------- */

int cholmod_print_triplet
(
    /* ---- input ---- */
    cholmod_triplet *T,	/* triplet matrix to print */
    const char *name,	/* printed name of triplet matrix */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_triplet (cholmod_triplet *, const char *, cholmod_common *);

/* -------------------------------------------------------------------------- */
/* cholmod_check_subset:  check a subset */
/* -------------------------------------------------------------------------- */

int cholmod_check_subset
(
    /* ---- input ---- */
    int *Set,		/* Set [0:len-1] is a subset of 0:n-1.  Duplicates OK */
    SuiteSparse_long len,   /* size of Set (an integer array) */
    size_t n,		/* 0:n-1 is valid range */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_subset (SuiteSparse_long *, SuiteSparse_long, size_t,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_subset:  print a subset */
/* -------------------------------------------------------------------------- */

int cholmod_print_subset
(
    /* ---- input ---- */
    int *Set,		/* Set [0:len-1] is a subset of 0:n-1.  Duplicates OK */
    SuiteSparse_long len,   /* size of Set (an integer array) */
    size_t n,		/* 0:n-1 is valid range */
    const char *name,	/* printed name of Set */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_subset (SuiteSparse_long *, SuiteSparse_long, size_t,
    const char *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_perm:  check a permutation */
/* -------------------------------------------------------------------------- */

int cholmod_check_perm
(
    /* ---- input ---- */
    int *Perm,		/* Perm [0:len-1] is a permutation of subset of 0:n-1 */
    size_t len,		/* size of Perm (an integer array) */
    size_t n,		/* 0:n-1 is valid range */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_perm (SuiteSparse_long *, size_t, size_t, cholmod_common *);

/* -------------------------------------------------------------------------- */
/* cholmod_print_perm:  print a permutation vector */
/* -------------------------------------------------------------------------- */

int cholmod_print_perm
(
    /* ---- input ---- */
    int *Perm,		/* Perm [0:len-1] is a permutation of subset of 0:n-1 */
    size_t len,		/* size of Perm (an integer array) */
    size_t n,		/* 0:n-1 is valid range */
    const char *name,	/* printed name of Perm */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_perm (SuiteSparse_long *, size_t, size_t, const char *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_parent:  check an elimination tree */
/* -------------------------------------------------------------------------- */

int cholmod_check_parent
(
    /* ---- input ---- */
    int *Parent,	/* Parent [0:n-1] is an elimination tree */
    size_t n,		/* size of Parent */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_parent (SuiteSparse_long *, size_t, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_parent */
/* -------------------------------------------------------------------------- */

int cholmod_print_parent
(
    /* ---- input ---- */
    int *Parent,	/* Parent [0:n-1] is an elimination tree */
    size_t n,		/* size of Parent */
    const char *name,	/* printed name of Parent */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_parent (SuiteSparse_long *, size_t, const char *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_read_sparse: read a sparse matrix from a file */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_read_sparse
(
    /* ---- input ---- */
    FILE *f,		/* file to read from, must already be open */
    /* --------------- */
    cholmod_common *Common
) ;

cholmod_sparse *cholmod_l_read_sparse (FILE *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_read_triplet: read a triplet matrix from a file */
/* -------------------------------------------------------------------------- */

cholmod_triplet *cholmod_read_triplet
(
    /* ---- input ---- */
    FILE *f,		/* file to read from, must already be open */
    /* --------------- */
    cholmod_common *Common
) ;

cholmod_triplet *cholmod_l_read_triplet (FILE *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_read_dense: read a dense matrix from a file */
/* -------------------------------------------------------------------------- */

cholmod_dense *cholmod_read_dense
(
    /* ---- input ---- */
    FILE *f,		/* file to read from, must already be open */
    /* --------------- */
    cholmod_common *Common
) ;

cholmod_dense *cholmod_l_read_dense (FILE *, cholmod_common *) ; 

/* -------------------------------------------------------------------------- */
/* cholmod_read_matrix: read a sparse or dense matrix from a file */
/* -------------------------------------------------------------------------- */

void *cholmod_read_matrix
(
    /* ---- input ---- */
    FILE *f,		/* file to read from, must already be open */
    int prefer,		/* If 0, a sparse matrix is always return as a
			 *	cholmod_triplet form.  It can have any stype
			 *	(symmetric-lower, unsymmetric, or
			 *	symmetric-upper).
			 * If 1, a sparse matrix is returned as an unsymmetric
			 *	cholmod_sparse form (A->stype == 0), with both
			 *	upper and lower triangular parts present.
			 *	This is what the MATLAB mread mexFunction does,
			 *	since MATLAB does not have an stype.
			 * If 2, a sparse matrix is returned with an stype of 0
			 *	or 1 (unsymmetric, or symmetric with upper part
			 *	stored).
			 * This argument has no effect for dense matrices.
			 */
    /* ---- output---- */
    int *mtype,		/* CHOLMOD_TRIPLET, CHOLMOD_SPARSE or CHOLMOD_DENSE */
    /* --------------- */
    cholmod_common *Common
) ;

void *cholmod_l_read_matrix (FILE *, int, int *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_write_sparse: write a sparse matrix to a file */
/* -------------------------------------------------------------------------- */

int cholmod_write_sparse
(
    /* ---- input ---- */
    FILE *f,		    /* file to write to, must already be open */
    cholmod_sparse *A,	    /* matrix to print */
    cholmod_sparse *Z,	    /* optional matrix with pattern of explicit zeros */
    const char *comments,    /* optional filename of comments to include */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_write_sparse (FILE *, cholmod_sparse *, cholmod_sparse *,
    const char *c, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_write_dense: write a dense matrix to a file */
/* -------------------------------------------------------------------------- */

int cholmod_write_dense
(
    /* ---- input ---- */
    FILE *f,		    /* file to write to, must already be open */
    cholmod_dense *X,	    /* matrix to print */
    const char *comments,    /* optional filename of comments to include */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_write_dense (FILE *, cholmod_dense *, const char *,
    cholmod_common *) ;
#endif
