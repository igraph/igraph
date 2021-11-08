/* luf.h (sparse LU-factorization) */

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

#ifndef LUF_H
#define LUF_H

#include "sva.h"

/***********************************************************************
*  The structure LUF describes sparse LU-factorization.
*
*  The LU-factorization has the following format:
*
*     A = F * V = P * L * U * Q,                                     (1)
*
*     F = P * L * P',                                                (2)
*
*     V = P * U * Q,                                                 (3)
*
*  where A is a given (unsymmetric) square matrix, F and V are matrix
*  factors actually computed, L is a lower triangular matrix with unity
*  diagonal, U is an upper triangular matrix, P and Q are permutation
*  matrices, P' is a matrix transposed to P. All the matrices have the
*  same order n.
*
*  Matrices F and V are stored in both row- and column-wise sparse
*  formats in the associated sparse vector area (SVA). Unity diagonal
*  elements of matrix F are not stored. Pivot elements of matrix V
*  (which correspond to diagonal elements of matrix U) are stored in
*  a separate ordinary array.
*
*  Permutation matrices P and Q are stored in ordinary arrays in both
*  row- and column-like formats.
*
*  Matrices L and U are completely defined by matrices F, V, P, and Q,
*  and therefore not stored explicitly. */

typedef struct LUF LUF;

struct LUF
{     /* sparse LU-factorization */
      int n;
      /* order of matrices A, F, V, P, Q */
      SVA *sva;
      /* associated sparse vector area (SVA) used to store rows and
       * columns of matrices F and V; note that different objects may
       * share the same SVA */
      /*--------------------------------------------------------------*/
      /* matrix F in row-wise format */
      /* during the factorization process this object is not used */
      int fr_ref;
      /* reference number of sparse vector in SVA, which is the first
       * row of matrix F */
#if 0 + 0
      int *fr_ptr = &sva->ptr[fr_ref-1];
      /* fr_ptr[0] is not used;
       * fr_ptr[i], 1 <= i <= n, is pointer to i-th row in SVA */
      int *fr_len = &sva->len[fr_ref-1];
      /* fr_len[0] is not used;
       * fr_len[i], 1 <= i <= n, is length of i-th row */
#endif
      /*--------------------------------------------------------------*/
      /* matrix F in column-wise format */
      /* during the factorization process this object is constructed
       * by columns */
      int fc_ref;
      /* reference number of sparse vector in SVA, which is the first
       * column of matrix F */
#if 0 + 0
      int *fc_ptr = &sva->ptr[fc_ref-1];
      /* fc_ptr[0] is not used;
       * fc_ptr[j], 1 <= j <= n, is pointer to j-th column in SVA */
      int *fc_len = &sva->len[fc_ref-1];
      /* fc_len[0] is not used;
       * fc_len[j], 1 <= j <= n, is length of j-th column */
#endif
      /*--------------------------------------------------------------*/
      /* matrix V in row-wise format */
      int vr_ref;
      /* reference number of sparse vector in SVA, which is the first
       * row of matrix V */
#if 0 + 0
      int *vr_ptr = &sva->ptr[vr_ref-1];
      /* vr_ptr[0] is not used;
       * vr_ptr[i], 1 <= i <= n, is pointer to i-th row in SVA */
      int *vr_len = &sva->len[vr_ref-1];
      /* vr_len[0] is not used;
       * vr_len[i], 1 <= i <= n, is length of i-th row */
      int *vr_cap = &sva->cap[vr_ref-1];
      /* vr_cap[0] is not used;
       * vr_cap[i], 1 <= i <= n, is capacity of i-th row */
#endif
      double *vr_piv; /* double vr_piv[1+n]; */
      /* vr_piv[0] is not used;
       * vr_piv[i], 1 <= i <= n, is pivot element of i-th row */
      /*--------------------------------------------------------------*/
      /* matrix V in column-wise format */
      /* during the factorization process this object contains only the
       * patterns (row indices) of columns of the active submatrix */
      int vc_ref;
      /* reference number of sparse vector in SVA, which is the first
       * column of matrix V */
#if 0 + 0
      int *vc_ptr = &sva->ptr[vc_ref-1];
      /* vc_ptr[0] is not used;
       * vc_ptr[j], 1 <= j <= n, is pointer to j-th column in SVA */
      int *vc_len = &sva->len[vc_ref-1];
      /* vc_len[0] is not used;
       * vc_len[j], 1 <= j <= n, is length of j-th column */
      int *vc_cap = &sva->cap[vc_ref-1];
      /* vc_cap[0] is not used;
       * vc_cap[j], 1 <= j <= n, is capacity of j-th column */
#endif
      /*--------------------------------------------------------------*/
      /* matrix P */
      int *pp_ind; /* int pp_ind[1+n]; */
      /* pp_ind[i] = j means that P[i,j] = 1 */
      int *pp_inv; /* int pp_inv[1+n]; */
      /* pp_inv[j] = i means that P[i,j] = 1 */
      /* if i-th row or column of matrix F is i'-th row or column of
       * matrix L, or if i-th row of matrix V is i'-th row of matrix U,
       * then pp_ind[i] = i' and pp_inv[i'] = i */
      /*--------------------------------------------------------------*/
      /* matrix Q */
      int *qq_ind; /* int qq_ind[1+n]; */
      /* qq_ind[i] = j means that Q[i,j] = 1 */
      int *qq_inv; /* int qq_inv[1+n]; */
      /* qq_inv[j] = i means that Q[i,j] = 1 */
      /* if j-th column of matrix V is j'-th column of matrix U, then
       * qq_ind[j'] = j and qq_inv[j] = j' */
};

#define luf_swap_u_rows(i1, i2) \
      do \
      {  int j1, j2; \
         j1 = pp_inv[i1], j2 = pp_inv[i2]; \
         pp_ind[j1] = i2, pp_inv[i2] = j1; \
         pp_ind[j2] = i1, pp_inv[i1] = j2; \
      } while (0)
/* swap rows i1 and i2 of matrix U = P'* V * Q' */

#define luf_swap_u_cols(j1, j2) \
      do \
      {  int i1, i2; \
         i1 = qq_ind[j1], i2 = qq_ind[j2]; \
         qq_ind[j1] = i2, qq_inv[i2] = j1; \
         qq_ind[j2] = i1, qq_inv[i1] = j2; \
      } while (0)
/* swap columns j1 and j2 of matrix U = P'* V * Q' */

#define luf_store_v_cols _glp_luf_store_v_cols
int luf_store_v_cols(LUF *luf, int (*col)(void *info, int j, int ind[],
      double val[]), void *info, int ind[], double val[]);
/* store matrix V = A in column-wise format */

#define luf_check_all _glp_luf_check_all
void luf_check_all(LUF *luf, int k);
/* check LU-factorization before k-th elimination step */

#define luf_build_v_rows _glp_luf_build_v_rows
void luf_build_v_rows(LUF *luf, int len[/*1+n*/]);
/* build matrix V in row-wise format */

#define luf_build_f_rows _glp_luf_build_f_rows
void luf_build_f_rows(LUF *luf, int len[/*1+n*/]);
/* build matrix F in row-wise format */

#define luf_build_v_cols _glp_luf_build_v_cols
void luf_build_v_cols(LUF *luf, int updat, int len[/*1+n*/]);
/* build matrix V in column-wise format */

#define luf_check_f_rc _glp_luf_check_f_rc
void luf_check_f_rc(LUF *luf);
/* check rows and columns of matrix F */

#define luf_check_v_rc _glp_luf_check_v_rc
void luf_check_v_rc(LUF *luf);
/* check rows and columns of matrix V */

#define luf_f_solve _glp_luf_f_solve
void luf_f_solve(LUF *luf, double x[/*1+n*/]);
/* solve system F * x = b */

#define luf_ft_solve _glp_luf_ft_solve
void luf_ft_solve(LUF *luf, double x[/*1+n*/]);
/* solve system F' * x = b */

#define luf_v_solve _glp_luf_v_solve
void luf_v_solve(LUF *luf, double b[/*1+n*/], double x[/*1+n*/]);
/* solve system V * x = b */

#define luf_vt_solve _glp_luf_vt_solve
void luf_vt_solve(LUF *luf, double b[/*1+n*/], double x[/*1+n*/]);
/* solve system V' * x = b */

#define luf_vt_solve1 _glp_luf_vt_solve1
void luf_vt_solve1(LUF *luf, double e[/*1+n*/], double y[/*1+n*/]);
/* solve system V' * y = e' to cause growth in y */

#define luf_estimate_norm _glp_luf_estimate_norm
double luf_estimate_norm(LUF *luf, double w1[/*1+n*/], double
      w2[/*1+n*/]);
/* estimate 1-norm of inv(A) */

#endif

/* eof */
