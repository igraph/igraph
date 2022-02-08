/* spxlp.h */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2015 Free Software Foundation, Inc.
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

#ifndef SPXLP_H
#define SPXLP_H

#include "bfd.h"

/***********************************************************************
*  The structure SPXLP describes LP problem and its current basis.
*
*  It is assumed that LP problem has the following formulation (this is
*  so called "working format"):
*
*     z = c'* x + c0 -> min                                          (1)
*
*     A * x = b                                                      (2)
*
*     l <= x <= u                                                    (3)
*
*  where:
*
*  x = (x[k]) is a n-vector of variables;
*
*  z is an objective function;
*
*  c = (c[k]) is a n-vector of objective coefficients;
*
*  c0 is a constant term of the objective function;
*
*  A = (a[i,k]) is a mxn-matrix of constraint coefficients;
*
*  b = (b[i]) is a m-vector of right-hand sides;
*
*  l = (l[k]) is a n-vector of lower bounds of variables;
*
*  u = (u[k]) is a n-vector of upper bounds of variables.
*
*  If variable x[k] has no lower (upper) bound, it is formally assumed
*  that l[k] = -inf (u[k] = +inf). Variable having no bounds is called
*  free (unbounded) variable. If l[k] = u[k], variable x[k] is assumed
*  to be fixed.
*
*  It is also assumed that matrix A has full row rank: rank(A) = m,
*  i.e. all its rows are linearly independent, so m <= n.
*
*  The (current) basis is defined by an appropriate permutation matrix
*  P of order n such that:
*
*             ( xB )
*     P * x = (    ),                                                (4)
*             ( xN )
*
*  where xB = (xB[i]) is a m-vector of basic variables, xN = (xN[j]) is
*  a (n-m)-vector of non-basic variables. If a non-basic variable xN[j]
*  has both lower and upper bounds, there is used an additional flag to
*  indicate which bound is active.
*
*  From (2) and (4) it follows that:
*
*     A * P'* P * x = b   <=>   B * xB + N * xN = b,                 (5)
*
*  where P' is a matrix transposed to P, and
*
*     A * P' = (B | N).                                              (6)
*
*  Here B is the basis matrix, which is a square non-singular matrix
*  of order m composed from columns of matrix A that correspond to
*  basic variables xB, and N is a mx(n-m) matrix composed from columns
*  of matrix A that correspond to non-basic variables xN. */

typedef struct SPXLP SPXLP;

struct SPXLP
{     /* LP problem data and its (current) basis */
      int m;
      /* number of equality constraints, m > 0 */
      int n;
      /* number of variables, n >= m */
      int nnz;
      /* number of non-zeros in constraint matrix A */
      /*--------------------------------------------------------------*/
      /* mxn-matrix A of constraint coefficients in sparse column-wise
       * format */
      int *A_ptr; /* int A_ptr[1+n+1]; */
      /* A_ptr[0] is not used;
       * A_ptr[k], 1 <= k <= n, is starting position of k-th column in
       * arrays A_ind and A_val; note that A_ptr[1] is always 1;
       * A_ptr[n+1] indicates the position after the last element in
       * arrays A_ind and A_val, i.e. A_ptr[n+1] = nnz+1, where nnz is
       * the number of non-zero elements in matrix A;
       * the length of k-th column (the number of non-zero elements in
       * that column) can be calculated as A_ptr[k+1] - A_ptr[k] */
      int *A_ind; /* int A_ind[1+nnz]; */
      /* row indices */
      double *A_val; /* double A_val[1+nnz]; */
      /* non-zero element values (constraint coefficients) */
      /*--------------------------------------------------------------*/
      /* principal vectors of LP formulation */
      double *b; /* double b[1+m]; */
      /* b[0] is not used;
       * b[i], 1 <= i <= m, is the right-hand side of i-th equality
       * constraint */
      double *c; /* double c[1+n]; */
      /* c[0] is the constant term of the objective function;
       * c[k], 1 <= k <= n, is the objective function coefficient at
       * variable x[k] */
      double *l; /* double l[1+n]; */
      /* l[0] is not used;
       * l[k], 1 <= k <= n, is the lower bound of variable x[k];
       * if x[k] has no lower bound, l[k] = -DBL_MAX */
      double *u; /* double u[1+n]; */
      /* u[0] is not used;
       * u[k], 1 <= k <= n, is the upper bound of variable u[k];
       * if x[k] has no upper bound, u[k] = +DBL_MAX;
       * note that l[k] = u[k] means that x[k] is fixed variable */
      /*--------------------------------------------------------------*/
      /* LP basis */
      int *head; /* int head[1+n]; */
      /* basis header, which is permutation matrix P (4):
       * head[0] is not used;
       * head[i] = k means that xB[i] = x[k], 1 <= i <= m;
       * head[m+j] = k, means that xN[j] = x[k], 1 <= j <= n-m */
      char *flag; /* char flag[1+n-m]; */
      /* flags of non-basic variables:
       * flag[0] is not used;
       * flag[j], 1 <= j <= n-m, indicates that non-basic variable
       * xN[j] is non-fixed and has its upper bound active */
      /*--------------------------------------------------------------*/
      /* basis matrix B of order m stored in factorized form */
      int valid;
      /* factorization validity flag */
      BFD *bfd;
      /* driver to factorization of the basis matrix */
};

#define spx_factorize _glp_spx_factorize
int spx_factorize(SPXLP *lp);
/* compute factorization of current basis matrix */

#define spx_eval_beta _glp_spx_eval_beta
void spx_eval_beta(SPXLP *lp, double beta[/*1+m*/]);
/* compute values of basic variables */

#define spx_eval_obj _glp_spx_eval_obj
double spx_eval_obj(SPXLP *lp, const double beta[/*1+m*/]);
/* compute value of objective function */

#define spx_eval_pi _glp_spx_eval_pi
void spx_eval_pi(SPXLP *lp, double pi[/*1+m*/]);
/* compute simplex multipliers */

#define spx_eval_dj _glp_spx_eval_dj
double spx_eval_dj(SPXLP *lp, const double pi[/*1+m*/], int j);
/* compute reduced cost of j-th non-basic variable */

#define spx_eval_tcol _glp_spx_eval_tcol
void spx_eval_tcol(SPXLP *lp, int j, double tcol[/*1+m*/]);
/* compute j-th column of simplex table */

#define spx_eval_rho _glp_spx_eval_rho
void spx_eval_rho(SPXLP *lp, int i, double rho[/*1+m*/]);
/* compute i-th row of basis matrix inverse */

#if 1 /* 31/III-2016 */
#define spx_eval_rho_s _glp_spx_eval_rho_s
void spx_eval_rho_s(SPXLP *lp, int i, FVS *rho);
/* sparse version of spx_eval_rho */
#endif

#define spx_eval_tij _glp_spx_eval_tij
double spx_eval_tij(SPXLP *lp, const double rho[/*1+m*/], int j);
/* compute element T[i,j] of simplex table */

#define spx_eval_trow _glp_spx_eval_trow
void spx_eval_trow(SPXLP *lp, const double rho[/*1+m*/], double
      trow[/*1+n-m*/]);
/* compute i-th row of simplex table */

#define spx_update_beta _glp_spx_update_beta
void spx_update_beta(SPXLP *lp, double beta[/*1+m*/], int p,
      int p_flag, int q, const double tcol[/*1+m*/]);
/* update values of basic variables */

#if 1 /* 30/III-2016 */
#define spx_update_beta_s _glp_spx_update_beta_s
void spx_update_beta_s(SPXLP *lp, double beta[/*1+m*/], int p,
      int p_flag, int q, const FVS *tcol);
/* sparse version of spx_update_beta */
#endif

#define spx_update_d _glp_spx_update_d
double spx_update_d(SPXLP *lp, double d[/*1+n-m*/], int p, int q,
      const double trow[/*1+n-m*/], const double tcol[/*1+m*/]);
/* update reduced costs of non-basic variables */

#if 1 /* 30/III-2016 */
#define spx_update_d_s _glp_spx_update_d_s
double spx_update_d_s(SPXLP *lp, double d[/*1+n-m*/], int p, int q,
      const FVS *trow, const FVS *tcol);
/* sparse version of spx_update_d */
#endif

#define spx_change_basis _glp_spx_change_basis
void spx_change_basis(SPXLP *lp, int p, int p_flag, int q);
/* change current basis to adjacent one */

#define spx_update_invb _glp_spx_update_invb
int spx_update_invb(SPXLP *lp, int i, int k);
/* update factorization of basis matrix */

#endif

/* eof */
