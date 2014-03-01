/* ========================================================================== */
/* === Include/cholmod_cholesky.h =========================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod_cholesky.h. Copyright (C) 2005-2013, Timothy A. Davis
 * CHOLMOD/Include/cholmod_cholesky.h is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* CHOLMOD Cholesky module.
 *
 * Sparse Cholesky routines: analysis, factorization, and solve.
 *
 * The primary routines are all that a user requires to order, analyze, and
 * factorize a sparse symmetric positive definite matrix A (or A*A'), and
 * to solve Ax=b (or A*A'x=b).  The primary routines rely on the secondary
 * routines, the CHOLMOD Core module, and the AMD and COLAMD packages.  They
 * make optional use of the CHOLMOD Supernodal and Partition modules, the
 * METIS package, and the CCOLAMD package.
 *
 * Primary routines:
 * -----------------
 *
 * cholmod_analyze		order and analyze (simplicial or supernodal)
 * cholmod_factorize		simplicial or supernodal Cholesky factorization
 * cholmod_solve		solve a linear system (simplicial or supernodal)
 * cholmod_solve2		like cholmod_solve, but reuse workspace
 * cholmod_spsolve		solve a linear system (sparse x and b)
 *
 * Secondary routines:
 * ------------------
 *
 * cholmod_analyze_p		analyze, with user-provided permutation or f set
 * cholmod_factorize_p		factorize, with user-provided permutation or f
 * cholmod_analyze_ordering	analyze a fill-reducing ordering
 * cholmod_etree		find the elimination tree
 * cholmod_rowcolcounts		compute the row/column counts of L
 * cholmod_amd			order using AMD
 * cholmod_colamd		order using COLAMD
 * cholmod_rowfac		incremental simplicial factorization
 * cholmod_rowfac_mask		rowfac, specific to LPDASA
 * cholmod_row_subtree		find the nonzero pattern of a row of L
 * cholmod_resymbol		recompute the symbolic pattern of L
 * cholmod_resymbol_noperm	recompute the symbolic pattern of L, no L->Perm
 * cholmod_postorder		postorder a tree
 *
 * Requires the Core module, and two packages: AMD and COLAMD.
 * Optionally uses the Supernodal and Partition modules.
 * Required by the Partition module.
 */

#ifndef CHOLMOD_CHOLESKY_H
#define CHOLMOD_CHOLESKY_H

#include "cholmod_config.h"
#include "cholmod_core.h"

#ifndef NPARTITION
#include "cholmod_partition.h"
#endif

#ifndef NSUPERNODAL
#include "cholmod_supernodal.h"
#endif

/* -------------------------------------------------------------------------- */
/* cholmod_analyze:  order and analyze (simplicial or supernodal) */
/* -------------------------------------------------------------------------- */

/* Orders and analyzes A, AA', PAP', or PAA'P' and returns a symbolic factor
 * that can later be passed to cholmod_factorize. */

cholmod_factor *cholmod_analyze 
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order and analyze */
    /* --------------- */
    cholmod_common *Common
) ;

cholmod_factor *cholmod_l_analyze (cholmod_sparse *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_analyze_p:  analyze, with user-provided permutation or f set */
/* -------------------------------------------------------------------------- */

/* Orders and analyzes A, AA', PAP', PAA'P', FF', or PFF'P and returns a
 * symbolic factor that can later be passed to cholmod_factorize, where
 * F = A(:,fset) if fset is not NULL and A->stype is zero.
 * UserPerm is tried if non-NULL.  */

cholmod_factor *cholmod_analyze_p
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order and analyze */
    int *UserPerm,	/* user-provided permutation, size A->nrow */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* --------------- */
    cholmod_common *Common
) ;

cholmod_factor *cholmod_l_analyze_p (cholmod_sparse *, SuiteSparse_long *,
    SuiteSparse_long *, size_t, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_analyze_p2:  analyze for sparse Cholesky or sparse QR */
/* -------------------------------------------------------------------------- */

cholmod_factor *cholmod_analyze_p2
(
    /* ---- input ---- */
    int for_cholesky,   /* if TRUE, then analyze for Cholesky; else for QR */
    cholmod_sparse *A,	/* matrix to order and analyze */
    int *UserPerm,	/* user-provided permutation, size A->nrow */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* --------------- */
    cholmod_common *Common
) ;

cholmod_factor *cholmod_l_analyze_p2 (int, cholmod_sparse *, SuiteSparse_long *,
    SuiteSparse_long *, size_t, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_factorize:  simplicial or supernodal Cholesky factorization */
/* -------------------------------------------------------------------------- */

/* Factorizes PAP' (or PAA'P' if A->stype is 0), using a factor obtained
 * from cholmod_analyze.  The analysis can be re-used simply by calling this
 * routine a second time with another matrix.  A must have the same nonzero
 * pattern as that passed to cholmod_analyze. */

int cholmod_factorize 
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to factorize */
    /* ---- in/out --- */
    cholmod_factor *L,	/* resulting factorization */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_factorize (cholmod_sparse *, cholmod_factor *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_factorize_p:  factorize, with user-provided permutation or fset */
/* -------------------------------------------------------------------------- */

/* Same as cholmod_factorize, but with more options. */

int cholmod_factorize_p
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to factorize */
    double beta [2],	/* factorize beta*I+A or beta*I+A'*A */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* ---- in/out --- */
    cholmod_factor *L,	/* resulting factorization */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_factorize_p (cholmod_sparse *, double *, SuiteSparse_long *,
    size_t, cholmod_factor *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_solve:  solve a linear system (simplicial or supernodal) */
/* -------------------------------------------------------------------------- */

/* Solves one of many linear systems with a dense right-hand-side, using the
 * factorization from cholmod_factorize (or as modified by any other CHOLMOD
 * routine).  D is identity for LL' factorizations. */

#define CHOLMOD_A    0		/* solve Ax=b */
#define CHOLMOD_LDLt 1		/* solve LDL'x=b */
#define CHOLMOD_LD   2		/* solve LDx=b */
#define CHOLMOD_DLt  3		/* solve DL'x=b */
#define CHOLMOD_L    4		/* solve Lx=b */
#define CHOLMOD_Lt   5		/* solve L'x=b */
#define CHOLMOD_D    6		/* solve Dx=b */
#define CHOLMOD_P    7		/* permute x=Px */
#define CHOLMOD_Pt   8		/* permute x=P'x */

cholmod_dense *cholmod_solve	/* returns the solution X */
(
    /* ---- input ---- */
    int sys,		/* system to solve */
    cholmod_factor *L,	/* factorization to use */
    cholmod_dense *B,	/* right-hand-side */
    /* --------------- */
    cholmod_common *Common
) ;

cholmod_dense *cholmod_l_solve (int, cholmod_factor *, cholmod_dense *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_solve2:  like cholmod_solve, but with reusable workspace */
/* -------------------------------------------------------------------------- */

int cholmod_solve2     /* returns TRUE on success, FALSE on failure */
(
    /* ---- input ---- */
    int sys,		            /* system to solve */
    cholmod_factor *L,	            /* factorization to use */
    cholmod_dense *B,               /* right-hand-side */
    cholmod_sparse *Bset,
    /* ---- output --- */
    cholmod_dense **X_Handle,       /* solution, allocated if need be */
    cholmod_sparse **Xset_Handle,
    /* ---- workspace  */
    cholmod_dense **Y_Handle,       /* workspace, or NULL */
    cholmod_dense **E_Handle,       /* workspace, or NULL */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_solve2 (int, cholmod_factor *, cholmod_dense *, cholmod_sparse *,
    cholmod_dense **, cholmod_sparse **, cholmod_dense **, cholmod_dense **,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_spsolve:  solve a linear system with a sparse right-hand-side */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_spsolve
(
    /* ---- input ---- */
    int sys,		/* system to solve */
    cholmod_factor *L,	/* factorization to use */
    cholmod_sparse *B,	/* right-hand-side */
    /* --------------- */
    cholmod_common *Common
) ;

cholmod_sparse *cholmod_l_spsolve (int, cholmod_factor *, cholmod_sparse *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_etree: find the elimination tree of A or A'*A */
/* -------------------------------------------------------------------------- */

int cholmod_etree
(
    /* ---- input ---- */
    cholmod_sparse *A,
    /* ---- output --- */
    int *Parent,	/* size ncol.  Parent [j] = p if p is the parent of j */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_etree (cholmod_sparse *, SuiteSparse_long *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_rowcolcounts: compute the row/column counts of L */
/* -------------------------------------------------------------------------- */

int cholmod_rowcolcounts
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int *Parent,	/* size nrow.  Parent [i] = p if p is the parent of i */
    int *Post,		/* size nrow.  Post [k] = i if i is the kth node in
			 * the postordered etree. */
    /* ---- output --- */
    int *RowCount,	/* size nrow. RowCount [i] = # entries in the ith row of
			 * L, including the diagonal. */
    int *ColCount,	/* size nrow. ColCount [i] = # entries in the ith
			 * column of L, including the diagonal. */
    int *First,		/* size nrow.  First [i] = k is the least postordering
			 * of any descendant of i. */
    int *Level,		/* size nrow.  Level [i] is the length of the path from
			 * i to the root, with Level [root] = 0. */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_rowcolcounts (cholmod_sparse *, SuiteSparse_long *, size_t,
    SuiteSparse_long *, SuiteSparse_long *, SuiteSparse_long *,
    SuiteSparse_long *, SuiteSparse_long *, SuiteSparse_long *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_analyze_ordering:  analyze a fill-reducing ordering */
/* -------------------------------------------------------------------------- */

int cholmod_analyze_ordering
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    int ordering,	/* ordering method used */
    int *Perm,		/* size n, fill-reducing permutation to analyze */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* ---- output --- */
    int *Parent,	/* size n, elimination tree */
    int *Post,		/* size n, postordering of elimination tree */
    int *ColCount,	/* size n, nnz in each column of L */
    /* ---- workspace  */
    int *First,		/* size nworkspace for cholmod_postorder */
    int *Level,		/* size n workspace for cholmod_postorder */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_analyze_ordering (cholmod_sparse *, int, SuiteSparse_long *,
    SuiteSparse_long *, size_t, SuiteSparse_long *, SuiteSparse_long *,
    SuiteSparse_long *, SuiteSparse_long *, SuiteSparse_long *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_amd:  order using AMD */
/* -------------------------------------------------------------------------- */

/* Finds a permutation P to reduce fill-in in the factorization of P*A*P'
 * or P*A*A'P' */

int cholmod_amd
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* ---- output --- */
    int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_amd (cholmod_sparse *, SuiteSparse_long *, size_t,
    SuiteSparse_long *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_colamd:  order using COLAMD */
/* -------------------------------------------------------------------------- */

/* Finds a permutation P to reduce fill-in in the factorization of P*A*A'*P'.
 * Orders F*F' where F = A (:,fset) if fset is not NULL */

int cholmod_colamd
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int postorder,	/* if TRUE, follow with a coletree postorder */
    /* ---- output --- */
    int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_colamd (cholmod_sparse *, SuiteSparse_long *, size_t, int,
    SuiteSparse_long *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_rowfac:  incremental simplicial factorization */
/* -------------------------------------------------------------------------- */

/* Partial or complete simplicial factorization.  Rows and columns kstart:kend-1
 * of L and D must be initially equal to rows/columns kstart:kend-1 of the
 * identity matrix.   Row k can only be factorized if all descendants of node
 * k in the elimination tree have been factorized. */

int cholmod_rowfac
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to factorize */
    cholmod_sparse *F,	/* used for A*A' case only. F=A' or A(:,fset)' */
    double beta [2],	/* factorize beta*I+A or beta*I+A'*A */
    size_t kstart,	/* first row to factorize */
    size_t kend,	/* last row to factorize is kend-1 */
    /* ---- in/out --- */
    cholmod_factor *L,
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_rowfac (cholmod_sparse *, cholmod_sparse *, double *, size_t,
    size_t, cholmod_factor *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_rowfac_mask:  incremental simplicial factorization */
/* -------------------------------------------------------------------------- */

/* cholmod_rowfac_mask is a version of cholmod_rowfac that is specific to
 * LPDASA.  It is unlikely to be needed by any other application. */

int cholmod_rowfac_mask
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to factorize */
    cholmod_sparse *F,	/* used for A*A' case only. F=A' or A(:,fset)' */
    double beta [2],	/* factorize beta*I+A or beta*I+A'*A */
    size_t kstart,	/* first row to factorize */
    size_t kend,	/* last row to factorize is kend-1 */
    int *mask,		/* if mask[i] >= 0, then set row i to zero */
    int *RLinkUp,	/* link list of rows to compute */
    /* ---- in/out --- */
    cholmod_factor *L,
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_rowfac_mask (cholmod_sparse *, cholmod_sparse *, double *, size_t,
    size_t, SuiteSparse_long *, SuiteSparse_long *, cholmod_factor *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_row_subtree:  find the nonzero pattern of a row of L */
/* -------------------------------------------------------------------------- */

/* Find the nonzero pattern of x for the system Lx=b where L = (0:k-1,0:k-1)
 * and b = kth column of A or A*A' (rows 0 to k-1 only) */

int cholmod_row_subtree
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    cholmod_sparse *F,	/* used for A*A' case only. F=A' or A(:,fset)' */
    size_t k,		/* row k of L */
    int *Parent,	/* elimination tree */
    /* ---- output --- */
    cholmod_sparse *R,	/* pattern of L(k,:), n-by-1 with R->nzmax >= n */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_row_subtree (cholmod_sparse *, cholmod_sparse *, size_t,
    SuiteSparse_long *, cholmod_sparse *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_lsolve_pattern: find the nonzero pattern of x=L\b */
/* -------------------------------------------------------------------------- */

int cholmod_lsolve_pattern
(
    /* ---- input ---- */
    cholmod_sparse *B,	/* sparse right-hand-side (a single sparse column) */
    cholmod_factor *L,	/* the factor L from which parent(i) is derived */
    /* ---- output --- */
    cholmod_sparse *X,	/* pattern of X=L\B, n-by-1 with X->nzmax >= n */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_lsolve_pattern (cholmod_sparse *, cholmod_factor *,
    cholmod_sparse *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_row_lsubtree:  find the nonzero pattern of a row of L */
/* -------------------------------------------------------------------------- */

/* Identical to cholmod_row_subtree, except that it finds the elimination tree
 * from L itself. */

int cholmod_row_lsubtree
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    int *Fi, size_t fnz,    /* nonzero pattern of kth row of A', not required
			     * for the symmetric case.  Need not be sorted. */
    size_t k,		/* row k of L */
    cholmod_factor *L,	/* the factor L from which parent(i) is derived */
    /* ---- output --- */
    cholmod_sparse *R,	/* pattern of L(k,:), n-by-1 with R->nzmax >= n */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_row_lsubtree (cholmod_sparse *, SuiteSparse_long *, size_t,
    size_t, cholmod_factor *, cholmod_sparse *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_resymbol:  recompute the symbolic pattern of L */
/* -------------------------------------------------------------------------- */

/* Remove entries from L that are not in the factorization of P*A*P', P*A*A'*P',
 * or P*F*F'*P' (depending on A->stype and whether fset is NULL or not).
 *
 * cholmod_resymbol is the same as cholmod_resymbol_noperm, except that it
 * first permutes A according to L->Perm.  A can be upper/lower/unsymmetric,
 * in contrast to cholmod_resymbol_noperm (which can be lower or unsym). */

int cholmod_resymbol 
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int pack,		/* if TRUE, pack the columns of L */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factorization, entries pruned on output */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_resymbol (cholmod_sparse *, SuiteSparse_long *, size_t, int,
    cholmod_factor *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_resymbol_noperm:  recompute the symbolic pattern of L, no L->Perm */
/* -------------------------------------------------------------------------- */

/* Remove entries from L that are not in the factorization of A, A*A',
 * or F*F' (depending on A->stype and whether fset is NULL or not). */

int cholmod_resymbol_noperm
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int pack,		/* if TRUE, pack the columns of L */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factorization, entries pruned on output */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_resymbol_noperm (cholmod_sparse *, SuiteSparse_long *, size_t, int,
    cholmod_factor *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_rcond:  compute rough estimate of reciprocal of condition number */
/* -------------------------------------------------------------------------- */

double cholmod_rcond	    /* return min(diag(L)) / max(diag(L)) */
(
    /* ---- input ---- */
    cholmod_factor *L,
    /* --------------- */
    cholmod_common *Common
) ;

double cholmod_l_rcond (cholmod_factor *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_postorder: Compute the postorder of a tree */
/* -------------------------------------------------------------------------- */

SuiteSparse_long cholmod_postorder	/* return # of nodes postordered */
(
    /* ---- input ---- */
    int *Parent,	/* size n. Parent [j] = p if p is the parent of j */
    size_t n,
    int *Weight_p,	/* size n, optional. Weight [j] is weight of node j */
    /* ---- output --- */
    int *Post,		/* size n. Post [k] = j is kth in postordered tree */
    /* --------------- */
    cholmod_common *Common
) ;

SuiteSparse_long cholmod_l_postorder (SuiteSparse_long *, size_t,
    SuiteSparse_long *, SuiteSparse_long *, cholmod_common *) ;

#endif
