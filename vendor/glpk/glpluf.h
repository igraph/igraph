/* glpluf.h (LU-factorization) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
*  2009, 2010 Andrew Makhorin, Department for Applied Informatics,
*  Moscow Aviation Institute, Moscow, Russia. All rights reserved.
*  E-mail: <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#ifndef GLPLUF_H
#define GLPLUF_H

/***********************************************************************
*  The structure LUF defines LU-factorization of a square matrix A and
*  is the following quartet:
*
*     [A] = (F, V, P, Q),                                            (1)
*
*  where F and V are such matrices that
*
*     A = F * V,                                                     (2)
*
*  and P and Q are such permutation matrices that the matrix
*
*     L = P * F * inv(P)                                             (3)
*
*  is lower triangular with unity diagonal, and the matrix
*
*     U = P * V * Q                                                  (4)
*
*  is upper triangular. All the matrices have the order n.
*
*  Matrices F and V are stored in row- and column-wise sparse format
*  as row and column linked lists of non-zero elements. Unity elements
*  on the main diagonal of matrix F are not stored. Pivot elements of
*  matrix V (which correspond to diagonal elements of matrix U) are
*  stored separately in an ordinary array.
*
*  Permutation matrices P and Q are stored in ordinary arrays in both
*  row- and column-like formats.
*
*  Matrices L and U are completely defined by matrices F, V, P, and Q
*  and therefore not stored explicitly.
*
*  The factorization (1)-(4) is a version of LU-factorization. Indeed,
*  from (3) and (4) it follows that:
*
*     F = inv(P) * L * P,
*
*     U = inv(P) * U * inv(Q),
*
*  and substitution into (2) leads to:
*
*     A = F * V = inv(P) * L * U * inv(Q).
*
*  For more details see the program documentation. */

typedef struct LUF LUF;

struct LUF
{     /* LU-factorization of a square matrix */
      int n_max;
      /* maximal value of n (increased automatically, if necessary) */
      int n;
      /* the order of matrices A, F, V, P, Q */
      int valid;
      /* the factorization is valid only if this flag is set */
      /*--------------------------------------------------------------*/
      /* matrix F in row-wise format */
      int *fr_ptr; /* int fr_ptr[1+n_max]; */
      /* fr_ptr[i], i = 1,...,n, is a pointer to the first element of
         i-th row in SVA */
      int *fr_len; /* int fr_len[1+n_max]; */
      /* fr_len[i], i = 1,...,n, is the number of elements in i-th row
         (except unity diagonal element) */
      /*--------------------------------------------------------------*/
      /* matrix F in column-wise format */
      int *fc_ptr; /* int fc_ptr[1+n_max]; */
      /* fc_ptr[j], j = 1,...,n, is a pointer to the first element of
         j-th column in SVA */
      int *fc_len; /* int fc_len[1+n_max]; */
      /* fc_len[j], j = 1,...,n, is the number of elements in j-th
         column (except unity diagonal element) */
      /*--------------------------------------------------------------*/
      /* matrix V in row-wise format */
      int *vr_ptr; /* int vr_ptr[1+n_max]; */
      /* vr_ptr[i], i = 1,...,n, is a pointer to the first element of
         i-th row in SVA */
      int *vr_len; /* int vr_len[1+n_max]; */
      /* vr_len[i], i = 1,...,n, is the number of elements in i-th row
         (except pivot element) */
      int *vr_cap; /* int vr_cap[1+n_max]; */
      /* vr_cap[i], i = 1,...,n, is the capacity of i-th row, i.e.
         maximal number of elements which can be stored in the row
         without relocating it, vr_cap[i] >= vr_len[i] */
      double *vr_piv; /* double vr_piv[1+n_max]; */
      /* vr_piv[p], p = 1,...,n, is the pivot element v[p,q] which
         corresponds to a diagonal element of matrix U = P*V*Q */
      /*--------------------------------------------------------------*/
      /* matrix V in column-wise format */
      int *vc_ptr; /* int vc_ptr[1+n_max]; */
      /* vc_ptr[j], j = 1,...,n, is a pointer to the first element of
         j-th column in SVA */
      int *vc_len; /* int vc_len[1+n_max]; */
      /* vc_len[j], j = 1,...,n, is the number of elements in j-th
         column (except pivot element) */
      int *vc_cap; /* int vc_cap[1+n_max]; */
      /* vc_cap[j], j = 1,...,n, is the capacity of j-th column, i.e.
         maximal number of elements which can be stored in the column
         without relocating it, vc_cap[j] >= vc_len[j] */
      /*--------------------------------------------------------------*/
      /* matrix P */
      int *pp_row; /* int pp_row[1+n_max]; */
      /* pp_row[i] = j means that P[i,j] = 1 */
      int *pp_col; /* int pp_col[1+n_max]; */
      /* pp_col[j] = i means that P[i,j] = 1 */
      /* if i-th row or column of matrix F is i'-th row or column of
         matrix L, or if i-th row of matrix V is i'-th row of matrix U,
         then pp_row[i'] = i and pp_col[i] = i' */
      /*--------------------------------------------------------------*/
      /* matrix Q */
      int *qq_row; /* int qq_row[1+n_max]; */
      /* qq_row[i] = j means that Q[i,j] = 1 */
      int *qq_col; /* int qq_col[1+n_max]; */
      /* qq_col[j] = i means that Q[i,j] = 1 */
      /* if j-th column of matrix V is j'-th column of matrix U, then
         qq_row[j] = j' and qq_col[j'] = j */
      /*--------------------------------------------------------------*/
      /* the Sparse Vector Area (SVA) is a set of locations used to
         store sparse vectors representing rows and columns of matrices
         F and V; each location is a doublet (ind, val), where ind is
         an index, and val is a numerical value of a sparse vector
         element; in the whole each sparse vector is a set of adjacent
         locations defined by a pointer to the first element and the
         number of elements; these pointer and number are stored in the
         corresponding matrix data structure (see above); the left part
         of SVA is used to store rows and columns of matrix V, and its
         right part is used to store rows and columns of matrix F; the
         middle part of SVA contains free (unused) locations */
      int sv_size;
      /* the size of SVA, in locations; all locations are numbered by
         integers 1, ..., n, and location 0 is not used; if necessary,
         the SVA size is automatically increased */
      int sv_beg, sv_end;
      /* SVA partitioning pointers:
         locations from 1 to sv_beg-1 belong to the left part
         locations from sv_beg to sv_end-1 belong to the middle part
         locations from sv_end to sv_size belong to the right part
         the size of the middle part is (sv_end - sv_beg) */
      int *sv_ind; /* sv_ind[1+sv_size]; */
      /* sv_ind[k], 1 <= k <= sv_size, is the index field of k-th
         location */
      double *sv_val; /* sv_val[1+sv_size]; */
      /* sv_val[k], 1 <= k <= sv_size, is the value field of k-th
         location */
      /*--------------------------------------------------------------*/
      /* in order to efficiently defragment the left part of SVA there
         is a doubly linked list of rows and columns of matrix V, where
         rows are numbered by 1, ..., n, while columns are numbered by
         n+1, ..., n+n, that allows uniquely identifying each row and
         column of V by only one integer; in this list rows and columns
         are ordered by ascending their pointers vr_ptr and vc_ptr */
      int sv_head;
      /* the number of leftmost row/column */
      int sv_tail;
      /* the number of rightmost row/column */
      int *sv_prev; /* int sv_prev[1+n_max+n_max]; */
      /* sv_prev[k], k = 1,...,n+n, is the number of a row/column which
         precedes k-th row/column */
      int *sv_next; /* int sv_next[1+n_max+n_max]; */
      /* sv_next[k], k = 1,...,n+n, is the number of a row/column which
         succedes k-th row/column */
      /*--------------------------------------------------------------*/
      /* working segment (used only during factorization) */
      double *vr_max; /* int vr_max[1+n_max]; */
      /* vr_max[i], 1 <= i <= n, is used only if i-th row of matrix V
         is active (i.e. belongs to the active submatrix), and is the
         largest magnitude of elements in i-th row; if vr_max[i] < 0,
         the largest magnitude is not known yet and should be computed
         by the pivoting routine */
      /*--------------------------------------------------------------*/
      /* in order to efficiently implement Markowitz strategy and Duff
         search technique there are two families {R[0], R[1], ..., R[n]}
         and {C[0], C[1], ..., C[n]}; member R[k] is the set of active
         rows of matrix V, which have k non-zeros, and member C[k] is
         the set of active columns of V, which have k non-zeros in the
         active submatrix (i.e. in the active rows); each set R[k] and
         C[k] is implemented as a separate doubly linked list */
      int *rs_head; /* int rs_head[1+n_max]; */
      /* rs_head[k], 0 <= k <= n, is the number of first active row,
         which has k non-zeros */
      int *rs_prev; /* int rs_prev[1+n_max]; */
      /* rs_prev[i], 1 <= i <= n, is the number of previous row, which
         has the same number of non-zeros as i-th row */
      int *rs_next; /* int rs_next[1+n_max]; */
      /* rs_next[i], 1 <= i <= n, is the number of next row, which has
         the same number of non-zeros as i-th row */
      int *cs_head; /* int cs_head[1+n_max]; */
      /* cs_head[k], 0 <= k <= n, is the number of first active column,
         which has k non-zeros (in the active rows) */
      int *cs_prev; /* int cs_prev[1+n_max]; */
      /* cs_prev[j], 1 <= j <= n, is the number of previous column,
         which has the same number of non-zeros (in the active rows) as
         j-th column */
      int *cs_next; /* int cs_next[1+n_max]; */
      /* cs_next[j], 1 <= j <= n, is the number of next column, which
         has the same number of non-zeros (in the active rows) as j-th
         column */
      /* (end of working segment) */
      /*--------------------------------------------------------------*/
      /* working arrays */
      int *flag; /* int flag[1+n_max]; */
      /* integer working array */
      double *work; /* double work[1+n_max]; */
      /* floating-point working array */
      /*--------------------------------------------------------------*/
      /* control parameters */
      int new_sva;
      /* new required size of the sparse vector area, in locations; set
         automatically by the factorizing routine */
      double piv_tol;
      /* threshold pivoting tolerance, 0 < piv_tol < 1; element v[i,j]
         of the active submatrix fits to be pivot if it satisfies to the
         stability criterion |v[i,j]| >= piv_tol * max |v[i,*]|, i.e. if
         it is not very small in the magnitude among other elements in
         the same row; decreasing this parameter gives better sparsity
         at the expense of numerical accuracy and vice versa */
      int piv_lim;
      /* maximal allowable number of pivot candidates to be considered;
         if piv_lim pivot candidates have been considered, the pivoting
         routine terminates the search with the best candidate found */
      int suhl;
      /* if this flag is set, the pivoting routine applies a heuristic
         proposed by Uwe Suhl: if a column of the active submatrix has
         no eligible pivot candidates (i.e. all its elements do not
         satisfy to the stability criterion), the routine excludes it
         from futher consideration until it becomes column singleton;
         in many cases this allows reducing the time needed for pivot
         searching */
      double eps_tol;
      /* epsilon tolerance; each element of the active submatrix, whose
         magnitude is less than eps_tol, is replaced by exact zero */
      double max_gro;
      /* maximal allowable growth of elements of matrix V during all
         the factorization process; if on some eliminaion step the ratio
         big_v / max_a (see below) becomes greater than max_gro, matrix
         A is considered as ill-conditioned (assuming that the pivoting
         tolerance piv_tol has an appropriate value) */
      /*--------------------------------------------------------------*/
      /* some statistics */
      int nnz_a;
      /* the number of non-zeros in matrix A */
      int nnz_f;
      /* the number of non-zeros in matrix F (except diagonal elements,
         which are not stored) */
      int nnz_v;
      /* the number of non-zeros in matrix V (except its pivot elements,
         which are stored in a separate array) */
      double max_a;
      /* the largest magnitude of elements of matrix A */
      double big_v;
      /* the largest magnitude of elements of matrix V appeared in the
         active submatrix during all the factorization process */
      int rank;
      /* estimated rank of matrix A */
};

/* return codes: */
#define LUF_ESING    1  /* singular matrix */
#define LUF_ECOND    2  /* ill-conditioned matrix */

#define luf_create_it _glp_luf_create_it
LUF *luf_create_it(void);
/* create LU-factorization */

#define luf_defrag_sva _glp_luf_defrag_sva
void luf_defrag_sva(LUF *luf);
/* defragment the sparse vector area */

#define luf_enlarge_row _glp_luf_enlarge_row
int luf_enlarge_row(LUF *luf, int i, int cap);
/* enlarge row capacity */

#define luf_enlarge_col _glp_luf_enlarge_col
int luf_enlarge_col(LUF *luf, int j, int cap);
/* enlarge column capacity */

#define luf_factorize _glp_luf_factorize
int luf_factorize(LUF *luf, int n, int (*col)(void *info, int j,
      int ind[], double val[]), void *info);
/* compute LU-factorization */

#define luf_f_solve _glp_luf_f_solve
void luf_f_solve(LUF *luf, int tr, double x[]);
/* solve system F*x = b or F'*x = b */

#define luf_v_solve _glp_luf_v_solve
void luf_v_solve(LUF *luf, int tr, double x[]);
/* solve system V*x = b or V'*x = b */

#define luf_a_solve _glp_luf_a_solve
void luf_a_solve(LUF *luf, int tr, double x[]);
/* solve system A*x = b or A'*x = b */

#define luf_delete_it _glp_luf_delete_it
void luf_delete_it(LUF *luf);
/* delete LU-factorization */

#endif

/* eof */
