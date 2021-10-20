/* btf.h (sparse block triangular LU-factorization) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2013-2014 Free Software Foundation, Inc.
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

#ifndef BTF_H
#define BTF_H

#include "sva.h"

/***********************************************************************
*  The structure BTF describes BT-factorization, which is sparse block
*  triangular LU-factorization.
*
*  The BT-factorization has the following format:
*
*     A = P * A~ * Q,                                                (1)
*
*  where A is a given (unsymmetric) square matrix, A~ is an upper block
*  triangular matrix (see below), P and Q are permutation matrices. All
*  the matrices have the same order n.
*
*  The matrix A~, which is a permuted version of the original matrix A,
*  has the following structure:
*
*     A~[1,1]  A~[1,2]  ...  A~[1,num-1]      A~[1,num]
*
*              A~[2,2]  ...  A~[2,num-1]      A~[2,num]
*
*                       . . . . . . . . .                            (2)
*
*                        A~[num-1,num-1]  A~[num-1,num]
*
*                                           A~[num,num]
*
*  where A~[i,j] is a submatrix called a "block," num is the number of
*  blocks. Each diagonal block A~[k,k] is a non-singular square matrix,
*  and each subdiagonal block A~[i,j], i > j, is a zero submatrix, thus
*  A~ is an upper block triangular matrix.
*
*  Permutation matrices P and Q are stored in ordinary arrays in both
*  row- and column-like formats.
*
*  The original matrix A is stored in both row- and column-wise sparse
*  formats in the associated sparse vector area (SVA). Should note that
*  elements of all diagonal blocks A~[k,k] in matrix A are set to zero
*  (i.e. removed), so only elements of non-diagonal blocks are stored.
*
*  Each diagonal block A~[k,k], 1 <= k <= num, is stored in the form of
*  LU-factorization (see the module LUF). */

typedef struct BTF BTF;

struct BTF
{     /* sparse block triangular LU-factorization */
      int n;
      /* order of matrices A, A~, P, Q */
      SVA *sva;
      /* associated sparse vector area used to store rows and columns
       * of matrix A as well as sparse vectors for LU-factorizations of
       * all diagonal blocks A~[k,k] */
      /*--------------------------------------------------------------*/
      /* matrix P */
      int *pp_ind; /* int pp_ind[1+n]; */
      /* pp_ind[i] = j means that P[i,j] = 1 */
      int *pp_inv; /* int pp_inv[1+n]; */
      /* pp_inv[j] = i means that P[i,j] = 1 */
      /* if i-th row of matrix A is i'-th row of matrix A~, then
       * pp_ind[i] = i' and pp_inv[i'] = i */
      /*--------------------------------------------------------------*/
      /* matrix Q */
      int *qq_ind; /* int qq_ind[1+n]; */
      /* qq_ind[i] = j means that Q[i,j] = 1 */
      int *qq_inv; /* int qq_inv[1+n]; */
      /* qq_inv[j] = i means that Q[i,j] = 1 */
      /* if j-th column of matrix A is j'-th column of matrix A~, then
       * qq_ind[j'] = j and qq_inv[j] = j' */
      /*--------------------------------------------------------------*/
      /* block triangular structure of matrix A~ */
      int num;
      /* number of diagonal blocks, 1 <= num <= n */
      int *beg; /* int beg[1+num+1]; */
      /* beg[0] is not used;
       * beg[k], 1 <= k <= num, is index of first row/column of k-th
       * block of matrix A~;
       * beg[num+1] is always n+1;
       * note that order (size) of k-th diagonal block can be computed
       * as beg[k+1] - beg[k] */
      /*--------------------------------------------------------------*/
      /* original matrix A in row-wise format */
      /* NOTE: elements of all diagonal blocks A~[k,k] are removed */
      int ar_ref;
      /* reference number of sparse vector in SVA, which is the first
       * row of matrix A */
#if 0 + 0
      int *ar_ptr = &sva->ptr[ar_ref-1];
      /* ar_ptr[0] is not used;
       * ar_ptr[i], 1 <= i <= n, is pointer to i-th row in SVA */
      int *ar_len = &sva->ptr[ar_ref-1];
      /* ar_len[0] is not used;
       * ar_len[i], 1 <= i <= n, is length of i-th row */
#endif
      /*--------------------------------------------------------------*/
      /* original matrix A in column-wise format */
      /* NOTE: elements of all diagonal blocks A~[k,k] are removed */
      int ac_ref;
      /* reference number of sparse vector in SVA, which is the first
       * column of matrix A */
#if 0 + 0
      int *ac_ptr = &sva->ptr[ac_ref-1];
      /* ac_ptr[0] is not used;
       * ac_ptr[j], 1 <= j <= n, is pointer to j-th column in SVA */
      int *ac_len = &sva->ptr[ac_ref-1];
      /* ac_len[0] is not used;
       * ac_len[j], 1 <= j <= n, is length of j-th column */
#endif
      /*--------------------------------------------------------------*/
      /* LU-factorizations of diagonal blocks A~[k,k] */
      /* to decrease overhead expenses similar arrays for all LUFs are
       * packed into a single array; for example, elements fr_ptr[1],
       * ..., fr_ptr[n1], where n1 = beg[2] - beg[1], are related to
       * LUF for first diagonal block A~[1,1], elements fr_ptr[n1+1],
       * ..., fr_ptr[n1+n2], where n2 = beg[3] - beg[2], are related to
       * LUF for second diagonal block A~[2,2], etc.; in other words,
       * elements related to LUF for k-th diagonal block A~[k,k] have
       * indices beg[k], beg[k]+1, ..., beg[k+1]-1 */
      /* for details about LUF see description of the LUF module */
      int fr_ref;
      /* reference number of sparse vector in SVA, which is the first
         row of matrix F for first diagonal block A~[1,1] */
      int fc_ref;
      /* reference number of sparse vector in SVA, which is the first
         column of matrix F for first diagonal block A~[1,1] */
      int vr_ref;
      /* reference number of sparse vector in SVA, which is the first
         row of matrix V for first diagonal block A~[1,1] */
      double *vr_piv; /* double vr_piv[1+n]; */
      /* vr_piv[0] is not used;
         vr_piv[1,...,n] are pivot elements for all diagonal blocks */
      int vc_ref;
      /* reference number of sparse vector in SVA, which is the first
         column of matrix V for first diagonal block A~[1,1] */
      int *p1_ind; /* int p1_ind[1+n]; */
      int *p1_inv; /* int p1_inv[1+n]; */
      int *q1_ind; /* int q1_ind[1+n]; */
      int *q1_inv; /* int q1_inv[1+n]; */
      /* permutation matrices P and Q for all diagonal blocks */
};

#define btf_store_a_cols _glp_btf_store_a_cols
int btf_store_a_cols(BTF *btf, int (*col)(void *info, int j, int ind[],
      double val[]), void *info, int ind[], double val[]);
/* store pattern of matrix A in column-wise format */

#define btf_make_blocks _glp_btf_make_blocks
int btf_make_blocks(BTF *btf);
/* permutations to block triangular form */

#define btf_check_blocks _glp_btf_check_blocks
void btf_check_blocks(BTF *btf);
/* check structure of matrix A~ */

#define btf_build_a_rows _glp_btf_build_a_rows
void btf_build_a_rows(BTF *btf, int len[/*1+n*/]);
/* build matrix A in row-wise format */

#define btf_a_solve _glp_btf_a_solve
void btf_a_solve(BTF *btf, double b[/*1+n*/], double x[/*1+n*/],
      double w1[/*1+n*/], double w2[/*1+n*/]);
/* solve system A * x = b */

#define btf_at_solve _glp_btf_at_solve
void btf_at_solve(BTF *btf, double b[/*1+n*/], double x[/*1+n*/],
      double w1[/*1+n*/], double w2[/*1+n*/]);
/* solve system A'* x = b */

#define btf_at_solve1 _glp_btf_at_solve1
void btf_at_solve1(BTF *btf, double e[/*1+n*/], double y[/*1+n*/],
      double w1[/*1+n*/], double w2[/*1+n*/]);
/* solve system A'* y = e' to cause growth in y */

#define btf_estimate_norm _glp_btf_estimate_norm
double btf_estimate_norm(BTF *btf, double w1[/*1+n*/], double
      w2[/*1+n*/], double w3[/*1+n*/], double w4[/*1+n*/]);
/* estimate 1-norm of inv(A) */

#endif

/* eof */
