/* sgf.h (sparse Gaussian factorizer) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2012-2013 Free Software Foundation, Inc.
*  Written by Andrew Makhorin <mao@gnu.org>.
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

#ifndef SGF_H
#define SGF_H

#include "luf.h"

typedef struct SGF SGF;

struct SGF
{     /* sparse Gaussian factorizer workspace */
      LUF *luf;
      /* LU-factorization being computed */
      /*--------------------------------------------------------------*/
      /* to efficiently choose pivot elements according to Markowitz
       * strategy, the search technique proposed by Iain Duff is used;
       * it is based on using two families of sets {R[0], ..., R[n]}
       * and {C[0], ..., C[n]}, where R[k] and C[k], 0 <= k <= n, are,
       * respectively, sets of rows and columns of the active submatrix
       * of matrix V having k non-zeros (i.e. whose length is k); each
       * set R[k] and C[k] is implemented as a doubly linked list */
      int *rs_head; /* int rs_head[1+n]; */
      /* rs_head[k], 0 <= k <= n, is the number of first row, which
       * has k non-zeros in the active submatrix */
      int *rs_prev; /* int rs_prev[1+n]; */
      /* rs_prev[0] is not used;
       * rs_prev[i], 1 <= i <= n, is the number of previous row, which
       * has the same number of non-zeros as i-th row;
       * rs_prev[i] < 0 means that i-th row is inactive */
      int *rs_next; /* int rs_next[1+n]; */
      /* rs_next[0] is not used;
       * rs_next[i], 1 <= i <= n, is the number of next row, which has
       * the same number of non-zeros as i-th row;
       * rs_next[i] < 0 means that i-th row is inactive */
      int *cs_head; /* int cs_head[1+n]; */
      /* cs_head[k], 0 <= k <= n, is the number of first column, which
       * has k non-zeros in the active submatrix */
      int *cs_prev; /* int cs_prev[1+n]; */
      /* cs_prev[0] is not used;
       * cs_prev[j], 1 <= j <= n, is the number of previous column,
       * which has the same number of non-zeros as j-th column;
       * cs_prev[j] < 0 means that j-th column is inactive */
      int *cs_next; /* int cs_next[1+n]; */
      /* cs_next[0] is not used;
       * cs_next[j], 1 <= j <= n, is the number of next column, which
       * has the same number of non-zeros as j-th column;
       * cs_next[j] < 0 means that j-th column is inactive */
      /* NOTE: cs_prev[j] = cs_next[j] = j means that j-th column was
       *       temporarily removed from corresponding set C[k] by the
       *       pivoting routine according to Uwe Suhl's heuristic */
      /*--------------------------------------------------------------*/
      /* working arrays */
      double *vr_max; /* int vr_max[1+n]; */
      /* vr_max[0] is not used;
       * vr_max[i], 1 <= i <= n, is used only if i-th row of matrix V
       * is active (i.e. belongs to the active submatrix), and is the
       * largest magnitude of elements in that row; if vr_max[i] < 0,
       * the largest magnitude is unknown yet */
      char *flag; /* char flag[1+n]; */
      /* boolean working array */
      double *work; /* double work[1+n]; */
      /* floating-point working array */
      /*--------------------------------------------------------------*/
      /* control parameters */
      int updat;
      /* if this flag is set, the matrix V is assumed to be updatable;
       * in this case factorized (non-active) part of V is stored in
       * the left part of SVA rather than in its right part */
      double piv_tol;
      /* threshold pivoting tolerance, 0 < piv_tol < 1; element v[i,j]
       * of the active submatrix fits to be pivot if it satisfies to
       * the stability criterion |v[i,j]| >= piv_tol * max |v[i,*]|,
       * i.e. if it is not very small in the magnitude among other
       * elements in the same row; decreasing this parameter gives
       * better sparsity at the expense of numerical accuracy and vice
       * versa */
      int piv_lim;
      /* maximal allowable number of pivot candidates to be considered;
       * if piv_lim pivot candidates have been considered, the pivoting
       * routine terminates the search with the best candidate found */
      int suhl;
      /* if this flag is set, the pivoting routine applies a heuristic
       * proposed by Uwe Suhl: if a column of the active submatrix has
       * no eligible pivot candidates (i.e. all its elements do not
       * satisfy to the stability criterion), the routine excludes it
       * from futher consideration until it becomes column singleton;
       * in many cases this allows reducing the time needed to choose
       * the pivot */
      double eps_tol;
      /* epsilon tolerance; each element of the active submatrix, whose
       * magnitude is less than eps_tol, is replaced by exact zero */
#if 0 /* FIXME */
      double den_lim;
      /* density limit; if the density of the active submatrix reaches
       * this limit, the factorization routine switches from sparse to
       * dense mode */
#endif
};

#define sgf_activate_row(i) \
      do \
      {  int len = vr_len[i]; \
         rs_prev[i] = 0; \
         rs_next[i] = rs_head[len]; \
         if (rs_next[i] != 0) \
            rs_prev[rs_next[i]] = i; \
         rs_head[len] = i; \
      } while (0)
/* include i-th row of matrix V in active set R[len] */

#define sgf_deactivate_row(i) \
      do \
      {  if (rs_prev[i] == 0) \
            rs_head[vr_len[i]] = rs_next[i]; \
         else \
            rs_next[rs_prev[i]] = rs_next[i]; \
         if (rs_next[i] == 0) \
            ; \
         else \
            rs_prev[rs_next[i]] = rs_prev[i]; \
         rs_prev[i] = rs_next[i] = -1; \
      } while (0)
/* remove i-th row of matrix V from active set R[len] */

#define sgf_activate_col(j) \
      do \
      {  int len = vc_len[j]; \
         cs_prev[j] = 0; \
         cs_next[j] = cs_head[len]; \
         if (cs_next[j] != 0) \
            cs_prev[cs_next[j]] = j; \
         cs_head[len] = j; \
      } while (0)
/* include j-th column of matrix V in active set C[len] */

#define sgf_deactivate_col(j) \
      do \
      {  if (cs_prev[j] == 0) \
            cs_head[vc_len[j]] = cs_next[j]; \
         else \
            cs_next[cs_prev[j]] = cs_next[j]; \
         if (cs_next[j] == 0) \
            ; \
         else \
            cs_prev[cs_next[j]] = cs_prev[j]; \
         cs_prev[j] = cs_next[j] = -1; \
      } while (0)
/* remove j-th column of matrix V from active set C[len] */

#define sgf_reduce_nuc _glp_sgf_reduce_nuc
int sgf_reduce_nuc(LUF *luf, int *k1, int *k2, int cnt[/*1+n*/],
      int list[/*1+n*/]);
/* initial reordering to minimize nucleus size */

#define sgf_singl_phase _glp_sgf_singl_phase
int sgf_singl_phase(LUF *luf, int k1, int k2, int updat,
      int ind[/*1+n*/], double val[/*1+n*/]);
/* compute LU-factorization (singleton phase) */

#define sgf_choose_pivot _glp_sgf_choose_pivot
int sgf_choose_pivot(SGF *sgf, int *p, int *q);
/* choose pivot element v[p,q] */

#define sgf_eliminate _glp_sgf_eliminate
int sgf_eliminate(SGF *sgf, int p, int q);
/* perform gaussian elimination */

#define sgf_dense_lu _glp_sgf_dense_lu
int sgf_dense_lu(int n, double a[], int r[], int c[], double eps);
/* compute dense LU-factorization with full pivoting */

#define sgf_dense_phase _glp_sgf_dense_phase
int sgf_dense_phase(LUF *luf, int k, int updat);
/* compute LU-factorization (dense phase) */

#define sgf_factorize _glp_sgf_factorize
int sgf_factorize(SGF *sgf, int singl);
/* compute LU-factorization (main routine) */

#endif

/* eof */
