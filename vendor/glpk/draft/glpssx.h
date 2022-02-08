/* glpssx.h (simplex method, rational arithmetic) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2003-2013 Free Software Foundation, Inc.
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

#ifndef GLPSSX_H
#define GLPSSX_H

#include "bfx.h"
#include "env.h"
#if 1 /* 25/XI-2017 */
#include "glpk.h"
#endif

typedef struct SSX SSX;

struct SSX
{     /* simplex solver workspace */
/*----------------------------------------------------------------------
// LP PROBLEM DATA
//
// It is assumed that LP problem has the following statement:
//
//    minimize (or maximize)
//
//       z = c[1]*x[1] + ... + c[m+n]*x[m+n] + c[0]                  (1)
//
//    subject to equality constraints
//
//       x[1] - a[1,1]*x[m+1] - ... - a[1,n]*x[m+n] = 0
//
//          .  .  .  .  .  .  .                                      (2)
//
//       x[m] - a[m,1]*x[m+1] + ... - a[m,n]*x[m+n] = 0
//
//    and bounds of variables
//
//         l[1] <= x[1]   <= u[1]
//
//          .  .  .  .  .  .  .                                      (3)
//
//       l[m+n] <= x[m+n] <= u[m+n]
//
// where:
// x[1], ..., x[m]      - auxiliary variables;
// x[m+1], ..., x[m+n]  - structural variables;
// z                    - objective function;
// c[1], ..., c[m+n]    - coefficients of the objective function;
// c[0]                 - constant term of the objective function;
// a[1,1], ..., a[m,n]  - constraint coefficients;
// l[1], ..., l[m+n]    - lower bounds of variables;
// u[1], ..., u[m+n]    - upper bounds of variables.
//
// Bounds of variables can be finite as well as inifinite. Besides,
// lower and upper bounds can be equal to each other. So the following
// five types of variables are possible:
//
//    Bounds of variable      Type of variable
//    -------------------------------------------------
//    -inf <  x[k] <  +inf    Free (unbounded) variable
//    l[k] <= x[k] <  +inf    Variable with lower bound
//    -inf <  x[k] <= u[k]    Variable with upper bound
//    l[k] <= x[k] <= u[k]    Double-bounded variable
//    l[k] =  x[k] =  u[k]    Fixed variable
//
// Using vector-matrix notations the LP problem (1)-(3) can be written
// as follows:
//
//    minimize (or maximize)
//
//       z = c * x + c[0]                                            (4)
//
//    subject to equality constraints
//
//       xR - A * xS = 0                                             (5)
//
//    and bounds of variables
//
//       l <= x <= u                                                 (6)
//
// where:
// xR                   - vector of auxiliary variables;
// xS                   - vector of structural variables;
// x = (xR, xS)         - vector of all variables;
// z                    - objective function;
// c                    - vector of objective coefficients;
// c[0]                 - constant term of the objective function;
// A                    - matrix of constraint coefficients (has m rows
//                        and n columns);
// l                    - vector of lower bounds of variables;
// u                    - vector of upper bounds of variables.
//
// The simplex method makes no difference between auxiliary and
// structural variables, so it is convenient to think the system of
// equality constraints (5) written in a homogeneous form:
//
//    (I | -A) * x = 0,                                              (7)
//
// where (I | -A) is an augmented (m+n)xm constraint matrix, I is mxm
// unity matrix whose columns correspond to auxiliary variables, and A
// is the original mxn constraint matrix whose columns correspond to
// structural variables. Note that only the matrix A is stored.
----------------------------------------------------------------------*/
      int m;
      /* number of rows (auxiliary variables), m > 0 */
      int n;
      /* number of columns (structural variables), n > 0 */
      int *type; /* int type[1+m+n]; */
      /* type[0] is not used;
         type[k], 1 <= k <= m+n, is the type of variable x[k]: */
#define SSX_FR          0     /* free (unbounded) variable */
#define SSX_LO          1     /* variable with lower bound */
#define SSX_UP          2     /* variable with upper bound */
#define SSX_DB          3     /* double-bounded variable */
#define SSX_FX          4     /* fixed variable */
      mpq_t *lb; /* mpq_t lb[1+m+n]; alias: l */
      /* lb[0] is not used;
         lb[k], 1 <= k <= m+n, is an lower bound of variable x[k];
         if x[k] has no lower bound, lb[k] is zero */
      mpq_t *ub; /* mpq_t ub[1+m+n]; alias: u */
      /* ub[0] is not used;
         ub[k], 1 <= k <= m+n, is an upper bound of variable x[k];
         if x[k] has no upper bound, ub[k] is zero;
         if x[k] is of fixed type, ub[k] is equal to lb[k] */
      int dir;
      /* optimization direction (sense of the objective function): */
#define SSX_MIN         0     /* minimization */
#define SSX_MAX         1     /* maximization */
      mpq_t *coef; /* mpq_t coef[1+m+n]; alias: c */
      /* coef[0] is a constant term of the objective function;
         coef[k], 1 <= k <= m+n, is a coefficient of the objective
         function at variable x[k];
         note that auxiliary variables also may have non-zero objective
         coefficients */
      int *A_ptr; /* int A_ptr[1+n+1]; */
      int *A_ind; /* int A_ind[A_ptr[n+1]]; */
      mpq_t *A_val; /* mpq_t A_val[A_ptr[n+1]]; */
      /* constraint matrix A (see (5)) in storage-by-columns format */
/*----------------------------------------------------------------------
// LP BASIS AND CURRENT BASIC SOLUTION
//
// The LP basis is defined by the following partition of the augmented
// constraint matrix (7):
//
//    (B | N) = (I | -A) * Q,                                        (8)
//
// where B is a mxm non-singular basis matrix whose columns correspond
// to basic variables xB, N is a mxn matrix whose columns correspond to
// non-basic variables xN, and Q is a permutation (m+n)x(m+n) matrix.
//
// From (7) and (8) it follows that
//
//    (I | -A) * x = (I | -A) * Q * Q' * x = (B | N) * (xB, xN),
//
// therefore
//
//    (xB, xN) = Q' * x,                                             (9)
//
// where x is the vector of all variables in the original order, xB is
// a vector of basic variables, xN is a vector of non-basic variables,
// Q' = inv(Q) is a matrix transposed to Q.
//
// Current values of non-basic variables xN[j], j = 1, ..., n, are not
// stored; they are defined implicitly by their statuses as follows:
//
//    0,             if xN[j] is free variable
//    lN[j],         if xN[j] is on its lower bound                 (10)
//    uN[j],         if xN[j] is on its upper bound
//    lN[j] = uN[j], if xN[j] is fixed variable
//
// where lN[j] and uN[j] are lower and upper bounds of xN[j].
//
// Current values of basic variables xB[i], i = 1, ..., m, are computed
// as follows:
//
//    beta = - inv(B) * N * xN,                                     (11)
//
// where current values of xN are defined by (10).
//
// Current values of simplex multipliers pi[i], i = 1, ..., m (which
// are values of Lagrange multipliers for equality constraints (7) also
// called shadow prices) are computed as follows:
//
//    pi = inv(B') * cB,                                            (12)
//
// where B' is a matrix transposed to B, cB is a vector of objective
// coefficients at basic variables xB.
//
// Current values of reduced costs d[j], j = 1, ..., n, (which are
// values of Langrange multipliers for active inequality constraints
// corresponding to non-basic variables) are computed as follows:
//
//    d = cN - N' * pi,                                             (13)
//
// where N' is a matrix transposed to N, cN is a vector of objective
// coefficients at non-basic variables xN.
----------------------------------------------------------------------*/
      int *stat; /* int stat[1+m+n]; */
      /* stat[0] is not used;
         stat[k], 1 <= k <= m+n, is the status of variable x[k]: */
#define SSX_BS          0     /* basic variable */
#define SSX_NL          1     /* non-basic variable on lower bound */
#define SSX_NU          2     /* non-basic variable on upper bound */
#define SSX_NF          3     /* non-basic free variable */
#define SSX_NS          4     /* non-basic fixed variable */
      int *Q_row; /* int Q_row[1+m+n]; */
      /* matrix Q in row-like format;
         Q_row[0] is not used;
         Q_row[i] = j means that q[i,j] = 1 */
      int *Q_col; /* int Q_col[1+m+n]; */
      /* matrix Q in column-like format;
         Q_col[0] is not used;
         Q_col[j] = i means that q[i,j] = 1 */
      /* if k-th column of the matrix (I | A) is k'-th column of the
         matrix (B | N), then Q_row[k] = k' and Q_col[k'] = k;
         if x[k] is xB[i], then Q_row[k] = i and Q_col[i] = k;
         if x[k] is xN[j], then Q_row[k] = m+j and Q_col[m+j] = k */
      BFX *binv;
      /* invertable form of the basis matrix B */
      mpq_t *bbar; /* mpq_t bbar[1+m]; alias: beta */
      /* bbar[0] is a value of the objective function;
         bbar[i], 1 <= i <= m, is a value of basic variable xB[i] */
      mpq_t *pi; /* mpq_t pi[1+m]; */
      /* pi[0] is not used;
         pi[i], 1 <= i <= m, is a simplex multiplier corresponding to
         i-th row (equality constraint) */
      mpq_t *cbar; /* mpq_t cbar[1+n]; alias: d */
      /* cbar[0] is not used;
         cbar[j], 1 <= j <= n, is a reduced cost of non-basic variable
         xN[j] */
/*----------------------------------------------------------------------
// SIMPLEX TABLE
//
// Due to (8) and (9) the system of equality constraints (7) for the
// current basis can be written as follows:
//
//    xB = A~ * xN,                                                 (14)
//
// where
//
//    A~ = - inv(B) * N                                             (15)
//
// is a mxn matrix called the simplex table.
//
// The revised simplex method uses only two components of A~, namely,
// pivot column corresponding to non-basic variable xN[q] chosen to
// enter the basis, and pivot row corresponding to basic variable xB[p]
// chosen to leave the basis.
//
// Pivot column alfa_q is q-th column of A~, so
//
//    alfa_q = A~ * e[q] = - inv(B) * N * e[q] = - inv(B) * N[q],   (16)
//
// where N[q] is q-th column of the matrix N.
//
// Pivot row alfa_p is p-th row of A~ or, equivalently, p-th column of
// A~', a matrix transposed to A~, so
//
//    alfa_p = A~' * e[p] = - N' * inv(B') * e[p] = - N' * rho_p,   (17)
//
// where (*)' means transposition, and
//
//    rho_p = inv(B') * e[p],                                       (18)
//
// is p-th column of inv(B') or, that is the same, p-th row of inv(B).
----------------------------------------------------------------------*/
      int p;
      /* number of basic variable xB[p], 1 <= p <= m, chosen to leave
         the basis */
      mpq_t *rho; /* mpq_t rho[1+m]; */
      /* p-th row of the inverse inv(B); see (18) */
      mpq_t *ap; /* mpq_t ap[1+n]; */
      /* p-th row of the simplex table; see (17) */
      int q;
      /* number of non-basic variable xN[q], 1 <= q <= n, chosen to
         enter the basis */
      mpq_t *aq; /* mpq_t aq[1+m]; */
      /* q-th column of the simplex table; see (16) */
/*--------------------------------------------------------------------*/
      int q_dir;
      /* direction in which non-basic variable xN[q] should change on
         moving to the adjacent vertex of the polyhedron:
         +1 means that xN[q] increases
         -1 means that xN[q] decreases */
      int p_stat;
      /* non-basic status which should be assigned to basic variable
         xB[p] when it has left the basis and become xN[q] */
      mpq_t delta;
      /* actual change of xN[q] in the adjacent basis (it has the same
         sign as q_dir) */
/*--------------------------------------------------------------------*/
#if 1 /* 25/XI-2017 */
      int msg_lev;
      /* verbosity level:
         GLP_MSG_OFF no output
         GLP_MSG_ERR report errors and warnings
         GLP_MSG_ON  normal output
         GLP_MSG_ALL highest verbosity */
#endif
      int it_lim;
      /* simplex iterations limit; if this value is positive, it is
         decreased by one each time when one simplex iteration has been
         performed, and reaching zero value signals the solver to stop
         the search; negative value means no iterations limit */
      int it_cnt;
      /* simplex iterations count; this count is increased by one each
         time when one simplex iteration has been performed */
      double tm_lim;
      /* searching time limit, in seconds; if this value is positive,
         it is decreased each time when one simplex iteration has been
         performed by the amount of time spent for the iteration, and
         reaching zero value signals the solver to stop the search;
         negative value means no time limit */
      double out_frq;
      /* output frequency, in seconds; this parameter specifies how
         frequently the solver sends information about the progress of
         the search to the standard output */
#if 0 /* 10/VI-2013 */
      glp_long tm_beg;
#else
      double tm_beg;
#endif
      /* starting time of the search, in seconds; the total time of the
         search is the difference between xtime() and tm_beg */
#if 0 /* 10/VI-2013 */
      glp_long tm_lag;
#else
      double tm_lag;
#endif
      /* the most recent time, in seconds, at which the progress of the
         the search was displayed */
};

#define ssx_create            _glp_ssx_create
#define ssx_factorize         _glp_ssx_factorize
#define ssx_get_xNj           _glp_ssx_get_xNj
#define ssx_eval_bbar         _glp_ssx_eval_bbar
#define ssx_eval_pi           _glp_ssx_eval_pi
#define ssx_eval_dj           _glp_ssx_eval_dj
#define ssx_eval_cbar         _glp_ssx_eval_cbar
#define ssx_eval_rho          _glp_ssx_eval_rho
#define ssx_eval_row          _glp_ssx_eval_row
#define ssx_eval_col          _glp_ssx_eval_col
#define ssx_chuzc             _glp_ssx_chuzc
#define ssx_chuzr             _glp_ssx_chuzr
#define ssx_update_bbar       _glp_ssx_update_bbar
#define ssx_update_pi         _glp_ssx_update_pi
#define ssx_update_cbar       _glp_ssx_update_cbar
#define ssx_change_basis      _glp_ssx_change_basis
#define ssx_delete            _glp_ssx_delete

#define ssx_phase_I           _glp_ssx_phase_I
#define ssx_phase_II          _glp_ssx_phase_II
#define ssx_driver            _glp_ssx_driver

SSX *ssx_create(int m, int n, int nnz);
/* create simplex solver workspace */

int ssx_factorize(SSX *ssx);
/* factorize the current basis matrix */

void ssx_get_xNj(SSX *ssx, int j, mpq_t x);
/* determine value of non-basic variable */

void ssx_eval_bbar(SSX *ssx);
/* compute values of basic variables */

void ssx_eval_pi(SSX *ssx);
/* compute values of simplex multipliers */

void ssx_eval_dj(SSX *ssx, int j, mpq_t dj);
/* compute reduced cost of non-basic variable */

void ssx_eval_cbar(SSX *ssx);
/* compute reduced costs of all non-basic variables */

void ssx_eval_rho(SSX *ssx);
/* compute p-th row of the inverse */

void ssx_eval_row(SSX *ssx);
/* compute pivot row of the simplex table */

void ssx_eval_col(SSX *ssx);
/* compute pivot column of the simplex table */

void ssx_chuzc(SSX *ssx);
/* choose pivot column */

void ssx_chuzr(SSX *ssx);
/* choose pivot row */

void ssx_update_bbar(SSX *ssx);
/* update values of basic variables */

void ssx_update_pi(SSX *ssx);
/* update simplex multipliers */

void ssx_update_cbar(SSX *ssx);
/* update reduced costs of non-basic variables */

void ssx_change_basis(SSX *ssx);
/* change current basis to adjacent one */

void ssx_delete(SSX *ssx);
/* delete simplex solver workspace */

int ssx_phase_I(SSX *ssx);
/* find primal feasible solution */

int ssx_phase_II(SSX *ssx);
/* find optimal solution */

int ssx_driver(SSX *ssx);
/* base driver to exact simplex method */

#endif

/* eof */
