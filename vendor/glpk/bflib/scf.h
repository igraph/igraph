/* scf.h (sparse updatable Schur-complement-based factorization) */

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

#ifndef SCF_H
#define SCF_H

#include "btf.h"
#include "ifu.h"
#include "luf.h"

/***********************************************************************
*  The structure SCF describes sparse updatable factorization based on
*  Schur complement.
*
*  The SCF-factorization has the following format:
*
*     ( A   A1~ )     ( A0  A1 )       ( R0    ) ( S0  S )
*     (         ) = P (        ) Q = P (       ) (       ) Q,        (1)
*     ( A2~ A3~ )     ( A2  A3 )       ( R   I ) (     C )
*
*  where:
*
*  A is current (unsymmetric) square matrix (not stored);
*
*  A1~, A2~, A3~ are some additional matrices (not stored);
*
*  A0 is initial (unsymmetric) square matrix (not stored);
*
*  A1, A2, A3 are some additional matrices (not stored);
*
*  R0 and S0 are matrices that define factorization of the initial
*  matrix A0 = R0 * S0 (stored in an invertable form);
*
*  R is a matrix defined from R * S0 = A2, so R = A2 * inv(S0) (stored
*  in row-wise sparse format);
*
*  S is a matrix defined from R0 * S = A1, so S = inv(R0) * A1 (stored
*  in column-wise sparse format);
*
*  C is Schur complement (to matrix A0) defined from R * S + C = A3,
*  so C = A3 - R * S = A3 - A2 * inv(A0) * A1 (stored in an invertable
*  form).
*
*  P, Q are permutation matrices (stored in both row- and column-like
*  formats). */

typedef struct SCF SCF;

struct SCF
{     /* Schur-complement-based factorization */
      int n;
      /* order of current matrix A */
      /*--------------------------------------------------------------*/
      /* initial matrix A0 = R0 * S0 of order n0 in invertable form */
      int n0;
      /* order of matrix A0 */
      int type;
      /* type of factorization used:
       * 1 - LU-factorization (R0 = F0, S0 = V0)
       * 2 - BT-factorization (R0 = I, S0 = A0) */
      union
      {  LUF *luf; /* type = 1 */
         BTF *btf; /* type = 2 */
      }  a0;
      /* factorization of matrix A0 */
      /*--------------------------------------------------------------*/
      /* augmented matrix (A0, A1; A2, A3) of order n0+nn */
      int nn_max;
      /* maximal number of additional rows and columns in the augmented
       * matrix (this limits the number of updates) */
      int nn;
      /* current number of additional rows and columns in the augmented
       * matrix, 0 <= nn <= nn_max */
      SVA *sva;
      /* associated sparse vector area (SVA) used to store rows of
       * matrix R and columns of matrix S */
      /*--------------------------------------------------------------*/
      /* nn*n0-matrix R in row-wise format */
      int rr_ref;
      /* reference number of sparse vector in SVA, which is the first
       * row of matrix R */
#if 0 + 0
      int *rr_ptr = &sva->ptr[rr_ref-1];
      /* rr_ptr[0] is not used;
       * rr_ptr[i], 1 <= i <= nn, is pointer to i-th row in SVA;
       * rr_ptr[nn+1,...,nn_max] are reserved locations */
      int *rr_len = &sva->len[rr_ref-1];
      /* rr_len[0] is not used;
       * rr_len[i], 1 <= i <= nn, is length of i-th row;
       * rr_len[nn+1,...,nn_max] are reserved locations */
#endif
      /*--------------------------------------------------------------*/
      /* n0*nn-matrix S in column-wise format */
      int ss_ref;
      /* reference number of sparse vector in SVA, which is the first
       * column of matrix S */
#if 0 + 0
      int *ss_ptr = &sva->ptr[ss_ref-1];
      /* ss_ptr[0] is not used;
       * ss_ptr[j], 1 <= j <= nn, is pointer to j-th column in SVA;
       * ss_ptr[nn+1,...,nn_max] are reserved locations */
      int *ss_len = &sva->len[ss_ref-1];
      /* ss_len[0] is not used;
       * ss_len[j], 1 <= j <= nn, is length of j-th column;
       * ss_len[nn+1,...,nn_max] are reserved locations */
#endif
      /*--------------------------------------------------------------*/
      /* Schur complement C of order nn in invertable form */
      IFU ifu;
      /* IFU-factorization of matrix C */
      /*--------------------------------------------------------------*/
      /* permutation matrix P of order n0+nn */
      int *pp_ind; /* int pp_ind[1+n0+nn_max]; */
      /* pp_ind[i] = j means that P[i,j] = 1 */
      int *pp_inv; /* int pp_inv[1+n0+nn_max]; */
      /* pp_inv[j] = i means that P[i,j] = 1 */
      /*--------------------------------------------------------------*/
      /* permutation matrix Q of order n0+nn */
      int *qq_ind; /* int qq_ind[1+n0+nn_max]; */
      /* qq_ind[i] = j means that Q[i,j] = 1 */
      int *qq_inv; /* int qq_inv[1+n0+nn_max]; */
      /* qq_inv[j] = i means that Q[i,j] = 1 */
};

#define scf_swap_q_cols(j1, j2) \
      do \
      {  int i1, i2; \
         i1 = qq_inv[j1], i2 = qq_inv[j2]; \
         qq_ind[i1] = j2, qq_inv[j2] = i1; \
         qq_ind[i2] = j1, qq_inv[j1] = i2; \
      }  while (0)
/* swap columns j1 and j2 of permutation matrix Q */

#define scf_r0_solve _glp_scf_r0_solve
void scf_r0_solve(SCF *scf, int tr, double x[/*1+n0*/]);
/* solve system R0 * x = b or R0'* x = b */

#define scf_s0_solve _glp_scf_s0_solve
void scf_s0_solve(SCF *scf, int tr, double x[/*1+n0*/],
      double w1[/*1+n0*/], double w2[/*1+n0*/], double w3[/*1+n0*/]);
/* solve system S0 * x = b or S0'* x = b */

#define scf_r_prod _glp_scf_r_prod
void scf_r_prod(SCF *scf, double y[/*1+nn*/], double a, const double
      x[/*1+n0*/]);
/* compute product y := y + alpha * R * x */

#define scf_rt_prod _glp_scf_rt_prod
void scf_rt_prod(SCF *scf, double y[/*1+n0*/], double a, const double
      x[/*1+nn*/]);
/* compute product y := y + alpha * R'* x */

#define scf_s_prod _glp_scf_s_prod
void scf_s_prod(SCF *scf, double y[/*1+n0*/], double a, const double
      x[/*1+nn*/]);
/* compute product y := y + alpha * S * x */

#define scf_st_prod _glp_scf_st_prod
void scf_st_prod(SCF *scf, double y[/*1+nn*/], double a, const double
      x[/*1+n0*/]);
/* compute product y := y + alpha * S'* x */

#define scf_a_solve _glp_scf_a_solve
void scf_a_solve(SCF *scf, double x[/*1+n*/],
      double w[/*1+n0+nn*/], double work1[/*1+max(n0,nn)*/],
      double work2[/*1+n*/], double work3[/*1+n*/]);
/* solve system A * x = b */

#define scf_at_solve _glp_scf_at_solve
void scf_at_solve(SCF *scf, double x[/*1+n*/],
      double w[/*1+n0+nn*/], double work1[/*1+max(n0,nn)*/],
      double work2[/*1+n*/], double work3[/*1+n*/]);
/* solve system A'* x = b */

#define scf_add_r_row _glp_scf_add_r_row
void scf_add_r_row(SCF *scf, const double w[/*1+n0*/]);
/* add new row to matrix R */

#define scf_add_s_col _glp_scf_add_s_col
void scf_add_s_col(SCF *scf, const double v[/*1+n0*/]);
/* add new column to matrix S */

#define scf_update_aug _glp_scf_update_aug
int scf_update_aug(SCF *scf, double b[/*1+n0*/], double d[/*1+n0*/],
      double f[/*1+nn*/], double g[/*1+nn*/], double h, int upd,
      double w1[/*1+n0*/], double w2[/*1+n0*/], double w3[/*1+n0*/]);
/* update factorization of augmented matrix */

#endif

/* eof */
