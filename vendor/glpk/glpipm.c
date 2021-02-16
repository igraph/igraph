/* glpipm.c */

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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wsometimes-uninitialized"
#endif

#include "glpipm.h"
#include "glpmat.h"

#define ITER_MAX 100
/* maximal number of iterations */

struct csa
{     /* common storage area */
      /*--------------------------------------------------------------*/
      /* LP data */
      int m;
      /* number of rows (equality constraints) */
      int n;
      /* number of columns (structural variables) */
      int *A_ptr; /* int A_ptr[1+m+1]; */
      int *A_ind; /* int A_ind[A_ptr[m+1]]; */
      double *A_val; /* double A_val[A_ptr[m+1]]; */
      /* mxn-matrix A in storage-by-rows format */
      double *b; /* double b[1+m]; */
      /* m-vector b of right-hand sides */
      double *c; /* double c[1+n]; */
      /* n-vector c of objective coefficients; c[0] is constant term of
         the objective function */
      /*--------------------------------------------------------------*/
      /* LP solution */
      double *x; /* double x[1+n]; */
      double *y; /* double y[1+m]; */
      double *z; /* double z[1+n]; */
      /* current point in primal-dual space; the best point on exit */
      /*--------------------------------------------------------------*/
      /* control parameters */
      const glp_iptcp *parm;
      /*--------------------------------------------------------------*/
      /* working arrays and variables */
      double *D; /* double D[1+n]; */
      /* diagonal nxn-matrix D = X*inv(Z), where X = diag(x[j]) and
         Z = diag(z[j]) */
      int *P; /* int P[1+m+m]; */
      /* permutation mxm-matrix P used to minimize fill-in in Cholesky
         factorization */
      int *S_ptr; /* int S_ptr[1+m+1]; */
      int *S_ind; /* int S_ind[S_ptr[m+1]]; */
      double *S_val; /* double S_val[S_ptr[m+1]]; */
      double *S_diag; /* double S_diag[1+m]; */
      /* symmetric mxm-matrix S = P*A*D*A'*P' whose upper triangular
         part without diagonal elements is stored in S_ptr, S_ind, and
         S_val in storage-by-rows format, diagonal elements are stored
         in S_diag */
      int *U_ptr; /* int U_ptr[1+m+1]; */
      int *U_ind; /* int U_ind[U_ptr[m+1]]; */
      double *U_val; /* double U_val[U_ptr[m+1]]; */
      double *U_diag; /* double U_diag[1+m]; */
      /* upper triangular mxm-matrix U defining Cholesky factorization
         S = U'*U; its non-diagonal elements are stored in U_ptr, U_ind,
         U_val in storage-by-rows format, diagonal elements are stored
         in U_diag */
      int iter;
      /* iteration number (0, 1, 2, ...); iter = 0 corresponds to the
         initial point */
      double obj;
      /* current value of the objective function */
      double rpi;
      /* relative primal infeasibility rpi = ||A*x-b||/(1+||b||) */
      double rdi;
      /* relative dual infeasibility rdi = ||A'*y+z-c||/(1+||c||) */
      double gap;
      /* primal-dual gap = |c'*x-b'*y|/(1+|c'*x|) which is a relative
         difference between primal and dual objective functions */
      double phi;
      /* merit function phi = ||A*x-b||/max(1,||b||) +
                            + ||A'*y+z-c||/max(1,||c||) +
                            + |c'*x-b'*y|/max(1,||b||,||c||) */
      double mu;
      /* duality measure mu = x'*z/n (used as barrier parameter) */
      double rmu;
      /* rmu = max(||A*x-b||,||A'*y+z-c||)/mu */
      double rmu0;
      /* the initial value of rmu on iteration 0 */
      double *phi_min; /* double phi_min[1+ITER_MAX]; */
      /* phi_min[k] = min(phi[k]), where phi[k] is the value of phi on
         k-th iteration, 0 <= k <= iter */
      int best_iter;
      /* iteration number, on which the value of phi reached its best
         (minimal) value */
      double *best_x; /* double best_x[1+n]; */
      double *best_y; /* double best_y[1+m]; */
      double *best_z; /* double best_z[1+n]; */
      /* best point (in the sense of the merit function phi) which has
         been reached on iteration iter_best */
      double best_obj;
      /* objective value at the best point */
      double *dx_aff; /* double dx_aff[1+n]; */
      double *dy_aff; /* double dy_aff[1+m]; */
      double *dz_aff; /* double dz_aff[1+n]; */
      /* affine scaling direction */
      double alfa_aff_p, alfa_aff_d;
      /* maximal primal and dual stepsizes in affine scaling direction,
         on which x and z are still non-negative */
      double mu_aff;
      /* duality measure mu_aff = x_aff'*z_aff/n in the boundary point
         x_aff' = x+alfa_aff_p*dx_aff, z_aff' = z+alfa_aff_d*dz_aff */
      double sigma;
      /* Mehrotra's heuristic parameter (0 <= sigma <= 1) */
      double *dx_cc; /* double dx_cc[1+n]; */
      double *dy_cc; /* double dy_cc[1+m]; */
      double *dz_cc; /* double dz_cc[1+n]; */
      /* centering corrector direction */
      double *dx; /* double dx[1+n]; */
      double *dy; /* double dy[1+m]; */
      double *dz; /* double dz[1+n]; */
      /* final combined direction dx = dx_aff+dx_cc, dy = dy_aff+dy_cc,
         dz = dz_aff+dz_cc */
      double alfa_max_p;
      double alfa_max_d;
      /* maximal primal and dual stepsizes in combined direction, on
         which x and z are still non-negative */
};

/***********************************************************************
*  initialize - allocate and initialize common storage area
*
*  This routine allocates and initializes the common storage area (CSA)
*  used by interior-point method routines. */

static void initialize(struct csa *csa)
{     int m = csa->m;
      int n = csa->n;
      int i;
      if (csa->parm->msg_lev >= GLP_MSG_ALL)
         xprintf("Matrix A has %d non-zeros\n", csa->A_ptr[m+1]-1);
      csa->D = xcalloc(1+n, sizeof(double));
      /* P := I */
      csa->P = xcalloc(1+m+m, sizeof(int));
      for (i = 1; i <= m; i++) csa->P[i] = csa->P[m+i] = i;
      /* S := A*A', symbolically */
      csa->S_ptr = xcalloc(1+m+1, sizeof(int));
      csa->S_ind = adat_symbolic(m, n, csa->P, csa->A_ptr, csa->A_ind,
         csa->S_ptr);
      if (csa->parm->msg_lev >= GLP_MSG_ALL)
         xprintf("Matrix S = A*A' has %d non-zeros (upper triangle)\n",
            csa->S_ptr[m+1]-1 + m);
      /* determine P using specified ordering algorithm */
      if (csa->parm->ord_alg == GLP_ORD_NONE)
      {  if (csa->parm->msg_lev >= GLP_MSG_ALL)
            xprintf("Original ordering is being used\n");
         for (i = 1; i <= m; i++)
            csa->P[i] = csa->P[m+i] = i;
      }
      else if (csa->parm->ord_alg == GLP_ORD_QMD)
      {  if (csa->parm->msg_lev >= GLP_MSG_ALL)
            xprintf("Minimum degree ordering (QMD)...\n");
         min_degree(m, csa->S_ptr, csa->S_ind, csa->P);
      }
      else if (csa->parm->ord_alg == GLP_ORD_AMD)
      {  if (csa->parm->msg_lev >= GLP_MSG_ALL)
            xprintf("Approximate minimum degree ordering (AMD)...\n");
         amd_order1(m, csa->S_ptr, csa->S_ind, csa->P);
      }
      else if (csa->parm->ord_alg == GLP_ORD_SYMAMD)
      {  if (csa->parm->msg_lev >= GLP_MSG_ALL)
            xprintf("Approximate minimum degree ordering (SYMAMD)...\n")
               ;
         symamd_ord(m, csa->S_ptr, csa->S_ind, csa->P);
      }
      else
         xassert(csa != csa);
      /* S := P*A*A'*P', symbolically */
      xfree(csa->S_ind);
      csa->S_ind = adat_symbolic(m, n, csa->P, csa->A_ptr, csa->A_ind,
         csa->S_ptr);
      csa->S_val = xcalloc(csa->S_ptr[m+1], sizeof(double));
      csa->S_diag = xcalloc(1+m, sizeof(double));
      /* compute Cholesky factorization S = U'*U, symbolically */
      if (csa->parm->msg_lev >= GLP_MSG_ALL)
         xprintf("Computing Cholesky factorization S = L*L'...\n");
      csa->U_ptr = xcalloc(1+m+1, sizeof(int));
      csa->U_ind = chol_symbolic(m, csa->S_ptr, csa->S_ind, csa->U_ptr);
      if (csa->parm->msg_lev >= GLP_MSG_ALL)
         xprintf("Matrix L has %d non-zeros\n", csa->U_ptr[m+1]-1 + m);
      csa->U_val = xcalloc(csa->U_ptr[m+1], sizeof(double));
      csa->U_diag = xcalloc(1+m, sizeof(double));
      csa->iter = 0;
      csa->obj = 0.0;
      csa->rpi = 0.0;
      csa->rdi = 0.0;
      csa->gap = 0.0;
      csa->phi = 0.0;
      csa->mu = 0.0;
      csa->rmu = 0.0;
      csa->rmu0 = 0.0;
      csa->phi_min = xcalloc(1+ITER_MAX, sizeof(double));
      csa->best_iter = 0;
      csa->best_x = xcalloc(1+n, sizeof(double));
      csa->best_y = xcalloc(1+m, sizeof(double));
      csa->best_z = xcalloc(1+n, sizeof(double));
      csa->best_obj = 0.0;
      csa->dx_aff = xcalloc(1+n, sizeof(double));
      csa->dy_aff = xcalloc(1+m, sizeof(double));
      csa->dz_aff = xcalloc(1+n, sizeof(double));
      csa->alfa_aff_p = 0.0;
      csa->alfa_aff_d = 0.0;
      csa->mu_aff = 0.0;
      csa->sigma = 0.0;
      csa->dx_cc = xcalloc(1+n, sizeof(double));
      csa->dy_cc = xcalloc(1+m, sizeof(double));
      csa->dz_cc = xcalloc(1+n, sizeof(double));
      csa->dx = csa->dx_aff;
      csa->dy = csa->dy_aff;
      csa->dz = csa->dz_aff;
      csa->alfa_max_p = 0.0;
      csa->alfa_max_d = 0.0;
      return;
}

/***********************************************************************
*  A_by_vec - compute y = A*x
*
*  This routine computes matrix-vector product y = A*x, where A is the
*  constraint matrix. */

static void A_by_vec(struct csa *csa, double x[], double y[])
{     /* compute y = A*x */
      int m = csa->m;
      int *A_ptr = csa->A_ptr;
      int *A_ind = csa->A_ind;
      double *A_val = csa->A_val;
      int i, t, beg, end;
      double temp;
      for (i = 1; i <= m; i++)
      {  temp = 0.0;
         beg = A_ptr[i], end = A_ptr[i+1];
         for (t = beg; t < end; t++) temp += A_val[t] * x[A_ind[t]];
         y[i] = temp;
      }
      return;
}

/***********************************************************************
*  AT_by_vec - compute y = A'*x
*
*  This routine computes matrix-vector product y = A'*x, where A' is a
*  matrix transposed to the constraint matrix A. */

static void AT_by_vec(struct csa *csa, double x[], double y[])
{     /* compute y = A'*x, where A' is transposed to A */
      int m = csa->m;
      int n = csa->n;
      int *A_ptr = csa->A_ptr;
      int *A_ind = csa->A_ind;
      double *A_val = csa->A_val;
      int i, j, t, beg, end;
      double temp;
      for (j = 1; j <= n; j++) y[j] = 0.0;
      for (i = 1; i <= m; i++)
      {  temp = x[i];
         if (temp == 0.0) continue;
         beg = A_ptr[i], end = A_ptr[i+1];
         for (t = beg; t < end; t++) y[A_ind[t]] += A_val[t] * temp;
      }
      return;
}

/***********************************************************************
*  decomp_NE - numeric factorization of matrix S = P*A*D*A'*P'
*
*  This routine implements numeric phase of Cholesky factorization of
*  the matrix S = P*A*D*A'*P', which is a permuted matrix of the normal
*  equation system. Matrix D is assumed to be already computed. */

static void decomp_NE(struct csa *csa)
{     adat_numeric(csa->m, csa->n, csa->P, csa->A_ptr, csa->A_ind,
         csa->A_val, csa->D, csa->S_ptr, csa->S_ind, csa->S_val,
         csa->S_diag);
      chol_numeric(csa->m, csa->S_ptr, csa->S_ind, csa->S_val,
         csa->S_diag, csa->U_ptr, csa->U_ind, csa->U_val, csa->U_diag);
      return;
}

/***********************************************************************
*  solve_NE - solve normal equation system
*
*  This routine solves the normal equation system:
*
*     A*D*A'*y = h.
*
*  It is assumed that the matrix A*D*A' has been previously factorized
*  by the routine decomp_NE.
*
*  On entry the array y contains the vector of right-hand sides h. On
*  exit this array contains the computed vector of unknowns y.
*
*  Once the vector y has been computed the routine checks for numeric
*  stability. If the residual vector:
*
*     r = A*D*A'*y - h
*
*  is relatively small, the routine returns zero, otherwise non-zero is
*  returned. */

static int solve_NE(struct csa *csa, double y[])
{     int m = csa->m;
      int n = csa->n;
      int *P = csa->P;
      int i, j, ret = 0;
      double *h, *r, *w;
      /* save vector of right-hand sides h */
      h = xcalloc(1+m, sizeof(double));
      for (i = 1; i <= m; i++) h[i] = y[i];
      /* solve normal equation system (A*D*A')*y = h */
      /* since S = P*A*D*A'*P' = U'*U, then A*D*A' = P'*U'*U*P, so we
         have inv(A*D*A') = P'*inv(U)*inv(U')*P */
      /* w := P*h */
      w = xcalloc(1+m, sizeof(double));
      for (i = 1; i <= m; i++) w[i] = y[P[i]];
      /* w := inv(U')*w */
      ut_solve(m, csa->U_ptr, csa->U_ind, csa->U_val, csa->U_diag, w);
      /* w := inv(U)*w */
      u_solve(m, csa->U_ptr, csa->U_ind, csa->U_val, csa->U_diag, w);
      /* y := P'*w */
      for (i = 1; i <= m; i++) y[i] = w[P[m+i]];
      xfree(w);
      /* compute residual vector r = A*D*A'*y - h */
      r = xcalloc(1+m, sizeof(double));
      /* w := A'*y */
      w = xcalloc(1+n, sizeof(double));
      AT_by_vec(csa, y, w);
      /* w := D*w */
      for (j = 1; j <= n; j++) w[j] *= csa->D[j];
      /* r := A*w */
      A_by_vec(csa, w, r);
      xfree(w);
      /* r := r - h */
      for (i = 1; i <= m; i++) r[i] -= h[i];
      /* check for numeric stability */
      for (i = 1; i <= m; i++)
      {  if (fabs(r[i]) / (1.0 + fabs(h[i])) > 1e-4)
         {  ret = 1;
            break;
         }
      }
      xfree(h);
      xfree(r);
      return ret;
}

/***********************************************************************
*  solve_NS - solve Newtonian system
*
*  This routine solves the Newtonian system:
*
*     A*dx               = p
*
*           A'*dy +   dz = q
*
*     Z*dx        + X*dz = r
*
*  where X = diag(x[j]), Z = diag(z[j]), by reducing it to the normal
*  equation system:
*
*     (A*inv(Z)*X*A')*dy = A*inv(Z)*(X*q-r)+p
*
*  (it is assumed that the matrix A*inv(Z)*X*A' has been factorized by
*  the routine decomp_NE).
*
*  Once vector dy has been computed the routine computes vectors dx and
*  dz as follows:
*
*     dx = inv(Z)*(X*(A'*dy-q)+r)
*
*     dz = inv(X)*(r-Z*dx)
*
*  The routine solve_NS returns the same code which was reported by the
*  routine solve_NE (see above). */

static int solve_NS(struct csa *csa, double p[], double q[], double r[],
      double dx[], double dy[], double dz[])
{     int m = csa->m;
      int n = csa->n;
      double *x = csa->x;
      double *z = csa->z;
      int i, j, ret;
      double *w = dx;
      /* compute the vector of right-hand sides A*inv(Z)*(X*q-r)+p for
         the normal equation system */
      for (j = 1; j <= n; j++)
         w[j] = (x[j] * q[j] - r[j]) / z[j];
      A_by_vec(csa, w, dy);
      for (i = 1; i <= m; i++) dy[i] += p[i];
      /* solve the normal equation system to compute vector dy */
      ret = solve_NE(csa, dy);
      /* compute vectors dx and dz */
      AT_by_vec(csa, dy, dx);
      for (j = 1; j <= n; j++)
      {  dx[j] = (x[j] * (dx[j] - q[j]) + r[j]) / z[j];
         dz[j] = (r[j] - z[j] * dx[j]) / x[j];
      }
      return ret;
}

/***********************************************************************
*  initial_point - choose initial point using Mehrotra's heuristic
*
*  This routine chooses a starting point using a heuristic proposed in
*  the paper:
*
*  S. Mehrotra. On the implementation of a primal-dual interior point
*  method. SIAM J. on Optim., 2(4), pp. 575-601, 1992.
*
*  The starting point x in the primal space is chosen as a solution of
*  the following least squares problem:
*
*     minimize    ||x||
*
*     subject to  A*x = b
*
*  which can be computed explicitly as follows:
*
*     x = A'*inv(A*A')*b
*
*  Similarly, the starting point (y, z) in the dual space is chosen as
*  a solution of the following least squares problem:
*
*     minimize    ||z||
*
*     subject to  A'*y + z = c
*
*  which can be computed explicitly as follows:
*
*     y = inv(A*A')*A*c
*
*     z = c - A'*y
*
*  However, some components of the vectors x and z may be non-positive
*  or close to zero, so the routine uses a Mehrotra's heuristic to find
*  a more appropriate starting point. */

static void initial_point(struct csa *csa)
{     int m = csa->m;
      int n = csa->n;
      double *b = csa->b;
      double *c = csa->c;
      double *x = csa->x;
      double *y = csa->y;
      double *z = csa->z;
      double *D = csa->D;
      int i, j;
      double dp, dd, ex, ez, xz;
      /* factorize A*A' */
      for (j = 1; j <= n; j++) D[j] = 1.0;
      decomp_NE(csa);
      /* x~ = A'*inv(A*A')*b */
      for (i = 1; i <= m; i++) y[i] = b[i];
      solve_NE(csa, y);
      AT_by_vec(csa, y, x);
      /* y~ = inv(A*A')*A*c */
      A_by_vec(csa, c, y);
      solve_NE(csa, y);
      /* z~ = c - A'*y~ */
      AT_by_vec(csa, y,z);
      for (j = 1; j <= n; j++) z[j] = c[j] - z[j];
      /* use Mehrotra's heuristic in order to choose more appropriate
         starting point with positive components of vectors x and z */
      dp = dd = 0.0;
      for (j = 1; j <= n; j++)
      {  if (dp < -1.5 * x[j]) dp = -1.5 * x[j];
         if (dd < -1.5 * z[j]) dd = -1.5 * z[j];
      }
      /* note that b = 0 involves x = 0, and c = 0 involves y = 0 and
         z = 0, so we need to be careful */
      if (dp == 0.0) dp = 1.5;
      if (dd == 0.0) dd = 1.5;
      ex = ez = xz = 0.0;
      for (j = 1; j <= n; j++)
      {  ex += (x[j] + dp);
         ez += (z[j] + dd);
         xz += (x[j] + dp) * (z[j] + dd);
      }
      dp += 0.5 * (xz / ez);
      dd += 0.5 * (xz / ex);
      for (j = 1; j <= n; j++)
      {  x[j] += dp;
         z[j] += dd;
         xassert(x[j] > 0.0 && z[j] > 0.0);
      }
      return;
}

/***********************************************************************
*  basic_info - perform basic computations at the current point
*
*  This routine computes the following quantities at the current point:
*
*  1) value of the objective function:
*
*     F = c'*x + c[0]
*
*  2) relative primal infeasibility:
*
*     rpi = ||A*x-b|| / (1+||b||)
*
*  3) relative dual infeasibility:
*
*     rdi = ||A'*y+z-c|| / (1+||c||)
*
*  4) primal-dual gap (relative difference between the primal and the
*     dual objective function values):
*
*     gap = |c'*x-b'*y| / (1+|c'*x|)
*
*  5) merit function:
*
*     phi = ||A*x-b|| / max(1,||b||) + ||A'*y+z-c|| / max(1,||c||) +
*
*         + |c'*x-b'*y| / max(1,||b||,||c||)
*
*  6) duality measure:
*
*     mu = x'*z / n
*
*  7) the ratio of infeasibility to mu:
*
*     rmu = max(||A*x-b||,||A'*y+z-c||) / mu
*
*  where ||*|| denotes euclidian norm, *' denotes transposition. */

static void basic_info(struct csa *csa)
{     int m = csa->m;
      int n = csa->n;
      double *b = csa->b;
      double *c = csa->c;
      double *x = csa->x;
      double *y = csa->y;
      double *z = csa->z;
      int i, j;
      double norm1, bnorm, norm2, cnorm, cx, by, *work, temp;
      /* compute value of the objective function */
      temp = c[0];
      for (j = 1; j <= n; j++) temp += c[j] * x[j];
      csa->obj = temp;
      /* norm1 = ||A*x-b|| */
      work = xcalloc(1+m, sizeof(double));
      A_by_vec(csa, x, work);
      norm1 = 0.0;
      for (i = 1; i <= m; i++)
         norm1 += (work[i] - b[i]) * (work[i] - b[i]);
      norm1 = sqrt(norm1);
      xfree(work);
      /* bnorm = ||b|| */
      bnorm = 0.0;
      for (i = 1; i <= m; i++) bnorm += b[i] * b[i];
      bnorm = sqrt(bnorm);
      /* compute relative primal infeasibility */
      csa->rpi = norm1 / (1.0 + bnorm);
      /* norm2 = ||A'*y+z-c|| */
      work = xcalloc(1+n, sizeof(double));
      AT_by_vec(csa, y, work);
      norm2 = 0.0;
      for (j = 1; j <= n; j++)
         norm2 += (work[j] + z[j] - c[j]) * (work[j] + z[j] - c[j]);
      norm2 = sqrt(norm2);
      xfree(work);
      /* cnorm = ||c|| */
      cnorm = 0.0;
      for (j = 1; j <= n; j++) cnorm += c[j] * c[j];
      cnorm = sqrt(cnorm);
      /* compute relative dual infeasibility */
      csa->rdi = norm2 / (1.0 + cnorm);
      /* by = b'*y */
      by = 0.0;
      for (i = 1; i <= m; i++) by += b[i] * y[i];
      /* cx = c'*x */
      cx = 0.0;
      for (j = 1; j <= n; j++) cx += c[j] * x[j];
      /* compute primal-dual gap */
      csa->gap = fabs(cx - by) / (1.0 + fabs(cx));
      /* compute merit function */
      csa->phi = 0.0;
      csa->phi += norm1 / (bnorm > 1.0 ? bnorm : 1.0);
      csa->phi += norm2 / (cnorm > 1.0 ? cnorm : 1.0);
      temp = 1.0;
      if (temp < bnorm) temp = bnorm;
      if (temp < cnorm) temp = cnorm;
      csa->phi += fabs(cx - by) / temp;
      /* compute duality measure */
      temp = 0.0;
      for (j = 1; j <= n; j++) temp += x[j] * z[j];
      csa->mu = temp / (double)n;
      /* compute the ratio of infeasibility to mu */
      csa->rmu = (norm1 > norm2 ? norm1 : norm2) / csa->mu;
      return;
}

/***********************************************************************
*  make_step - compute next point using Mehrotra's technique
*
*  This routine computes the next point using the predictor-corrector
*  technique proposed in the paper:
*
*  S. Mehrotra. On the implementation of a primal-dual interior point
*  method. SIAM J. on Optim., 2(4), pp. 575-601, 1992.
*
*  At first, the routine computes so called affine scaling (predictor)
*  direction (dx_aff,dy_aff,dz_aff) which is a solution of the system:
*
*     A*dx_aff                       = b - A*x
*
*               A'*dy_aff +   dz_aff = c - A'*y - z
*
*     Z*dx_aff            + X*dz_aff = - X*Z*e
*
*  where (x,y,z) is the current point, X = diag(x[j]), Z = diag(z[j]),
*  e = (1,...,1)'.
*
*  Then, the routine computes the centering parameter sigma, using the
*  following Mehrotra's heuristic:
*
*     alfa_aff_p = inf{0 <= alfa <= 1 | x+alfa*dx_aff >= 0}
*
*     alfa_aff_d = inf{0 <= alfa <= 1 | z+alfa*dz_aff >= 0}
*
*     mu_aff = (x+alfa_aff_p*dx_aff)'*(z+alfa_aff_d*dz_aff)/n
*
*     sigma = (mu_aff/mu)^3
*
*  where alfa_aff_p is the maximal stepsize along the affine scaling
*  direction in the primal space, alfa_aff_d is the maximal stepsize
*  along the same direction in the dual space.
*
*  After determining sigma the routine computes so called centering
*  (corrector) direction (dx_cc,dy_cc,dz_cc) which is the solution of
*  the system:
*
*     A*dx_cc                     = 0
*
*              A'*dy_cc +   dz_cc = 0
*
*     Z*dx_cc           + X*dz_cc = sigma*mu*e - X*Z*e
*
*  Finally, the routine computes the combined direction
*
*     (dx,dy,dz) = (dx_aff,dy_aff,dz_aff) + (dx_cc,dy_cc,dz_cc)
*
*  and determines maximal primal and dual stepsizes along the combined
*  direction:
*
*     alfa_max_p = inf{0 <= alfa <= 1 | x+alfa*dx >= 0}
*
*     alfa_max_d = inf{0 <= alfa <= 1 | z+alfa*dz >= 0}
*
*  In order to prevent the next point to be too close to the boundary
*  of the positive ortant, the routine decreases maximal stepsizes:
*
*     alfa_p = gamma_p * alfa_max_p
*
*     alfa_d = gamma_d * alfa_max_d
*
*  where gamma_p and gamma_d are scaling factors, and computes the next
*  point:
*
*     x_new = x + alfa_p * dx
*
*     y_new = y + alfa_d * dy
*
*     z_new = z + alfa_d * dz
*
*  which becomes the current point on the next iteration. */

static int make_step(struct csa *csa)
{     int m = csa->m;
      int n = csa->n;
      double *b = csa->b;
      double *c = csa->c;
      double *x = csa->x;
      double *y = csa->y;
      double *z = csa->z;
      double *dx_aff = csa->dx_aff;
      double *dy_aff = csa->dy_aff;
      double *dz_aff = csa->dz_aff;
      double *dx_cc = csa->dx_cc;
      double *dy_cc = csa->dy_cc;
      double *dz_cc = csa->dz_cc;
      double *dx = csa->dx;
      double *dy = csa->dy;
      double *dz = csa->dz;
      int i, j, ret = 0;
      double temp, gamma_p, gamma_d, *p, *q, *r;
      /* allocate working arrays */
      p = xcalloc(1+m, sizeof(double));
      q = xcalloc(1+n, sizeof(double));
      r = xcalloc(1+n, sizeof(double));
      /* p = b - A*x */
      A_by_vec(csa, x, p);
      for (i = 1; i <= m; i++) p[i] = b[i] - p[i];
      /* q = c - A'*y - z */
      AT_by_vec(csa, y,q);
      for (j = 1; j <= n; j++) q[j] = c[j] - q[j] - z[j];
      /* r = - X * Z * e */
      for (j = 1; j <= n; j++) r[j] = - x[j] * z[j];
      /* solve the first Newtonian system */
      if (solve_NS(csa, p, q, r, dx_aff, dy_aff, dz_aff))
      {  ret = 1;
         goto done;
      }
      /* alfa_aff_p = inf{0 <= alfa <= 1 | x + alfa*dx_aff >= 0} */
      /* alfa_aff_d = inf{0 <= alfa <= 1 | z + alfa*dz_aff >= 0} */
      csa->alfa_aff_p = csa->alfa_aff_d = 1.0;
      for (j = 1; j <= n; j++)
      {  if (dx_aff[j] < 0.0)
         {  temp = - x[j] / dx_aff[j];
            if (csa->alfa_aff_p > temp) csa->alfa_aff_p = temp;
         }
         if (dz_aff[j] < 0.0)
         {  temp = - z[j] / dz_aff[j];
            if (csa->alfa_aff_d > temp) csa->alfa_aff_d = temp;
         }
      }
      /* mu_aff = (x+alfa_aff_p*dx_aff)' * (z+alfa_aff_d*dz_aff) / n */
      temp = 0.0;
      for (j = 1; j <= n; j++)
         temp += (x[j] + csa->alfa_aff_p * dx_aff[j]) *
                 (z[j] + csa->alfa_aff_d * dz_aff[j]);
      csa->mu_aff = temp / (double)n;
      /* sigma = (mu_aff/mu)^3 */
      temp = csa->mu_aff / csa->mu;
      csa->sigma = temp * temp * temp;
      /* p = 0 */
      for (i = 1; i <= m; i++) p[i] = 0.0;
      /* q = 0 */
      for (j = 1; j <= n; j++) q[j] = 0.0;
      /* r = sigma * mu * e - X * Z * e */
      for (j = 1; j <= n; j++)
         r[j] = csa->sigma * csa->mu - dx_aff[j] * dz_aff[j];
      /* solve the second Newtonian system with the same coefficients
         but with altered right-hand sides */
      if (solve_NS(csa, p, q, r, dx_cc, dy_cc, dz_cc))
      {  ret = 1;
         goto done;
      }
      /* (dx,dy,dz) = (dx_aff,dy_aff,dz_aff) + (dx_cc,dy_cc,dz_cc) */
      for (j = 1; j <= n; j++) dx[j] = dx_aff[j] + dx_cc[j];
      for (i = 1; i <= m; i++) dy[i] = dy_aff[i] + dy_cc[i];
      for (j = 1; j <= n; j++) dz[j] = dz_aff[j] + dz_cc[j];
      /* alfa_max_p = inf{0 <= alfa <= 1 | x + alfa*dx >= 0} */
      /* alfa_max_d = inf{0 <= alfa <= 1 | z + alfa*dz >= 0} */
      csa->alfa_max_p = csa->alfa_max_d = 1.0;
      for (j = 1; j <= n; j++)
      {  if (dx[j] < 0.0)
         {  temp = - x[j] / dx[j];
            if (csa->alfa_max_p > temp) csa->alfa_max_p = temp;
         }
         if (dz[j] < 0.0)
         {  temp = - z[j] / dz[j];
            if (csa->alfa_max_d > temp) csa->alfa_max_d = temp;
         }
      }
      /* determine scale factors (not implemented yet) */
      gamma_p = 0.90;
      gamma_d = 0.90;
      /* compute the next point */
      for (j = 1; j <= n; j++)
      {  x[j] += gamma_p * csa->alfa_max_p * dx[j];
         xassert(x[j] > 0.0);
      }
      for (i = 1; i <= m; i++)
         y[i] += gamma_d * csa->alfa_max_d * dy[i];
      for (j = 1; j <= n; j++)
      {  z[j] += gamma_d * csa->alfa_max_d * dz[j];
         xassert(z[j] > 0.0);
      }
done: /* free working arrays */
      xfree(p);
      xfree(q);
      xfree(r);
      return ret;
}

/***********************************************************************
*  terminate - deallocate common storage area
*
*  This routine frees all memory allocated to the common storage area
*  used by interior-point method routines. */

static void terminate(struct csa *csa)
{     xfree(csa->D);
      xfree(csa->P);
      xfree(csa->S_ptr);
      xfree(csa->S_ind);
      xfree(csa->S_val);
      xfree(csa->S_diag);
      xfree(csa->U_ptr);
      xfree(csa->U_ind);
      xfree(csa->U_val);
      xfree(csa->U_diag);
      xfree(csa->phi_min);
      xfree(csa->best_x);
      xfree(csa->best_y);
      xfree(csa->best_z);
      xfree(csa->dx_aff);
      xfree(csa->dy_aff);
      xfree(csa->dz_aff);
      xfree(csa->dx_cc);
      xfree(csa->dy_cc);
      xfree(csa->dz_cc);
      return;
}

/***********************************************************************
*  ipm_main - main interior-point method routine
*
*  This is a main routine of the primal-dual interior-point method.
*
*  The routine ipm_main returns one of the following codes:
*
*  0 - optimal solution found;
*  1 - problem has no feasible (primal or dual) solution;
*  2 - no convergence;
*  3 - iteration limit exceeded;
*  4 - numeric instability on solving Newtonian system.
*
*  In case of non-zero return code the routine returns the best point,
*  which has been reached during optimization. */

static int ipm_main(struct csa *csa)
{     int m = csa->m;
      int n = csa->n;
      int i, j, status;
      double temp;
      /* choose initial point using Mehrotra's heuristic */
      if (csa->parm->msg_lev >= GLP_MSG_ALL)
         xprintf("Guessing initial point...\n");
      initial_point(csa);
      /* main loop starts here */
      if (csa->parm->msg_lev >= GLP_MSG_ALL)
         xprintf("Optimization begins...\n");
      for (;;)
      {  /* perform basic computations at the current point */
         basic_info(csa);
         /* save initial value of rmu */
         if (csa->iter == 0) csa->rmu0 = csa->rmu;
         /* accumulate values of min(phi[k]) and save the best point */
         xassert(csa->iter <= ITER_MAX);
         if (csa->iter == 0 || csa->phi_min[csa->iter-1] > csa->phi)
         {  csa->phi_min[csa->iter] = csa->phi;
            csa->best_iter = csa->iter;
            for (j = 1; j <= n; j++) csa->best_x[j] = csa->x[j];
            for (i = 1; i <= m; i++) csa->best_y[i] = csa->y[i];
            for (j = 1; j <= n; j++) csa->best_z[j] = csa->z[j];
            csa->best_obj = csa->obj;
         }
         else
            csa->phi_min[csa->iter] = csa->phi_min[csa->iter-1];
         /* display information at the current point */
         if (csa->parm->msg_lev >= GLP_MSG_ON)
            xprintf("%3d: obj = %17.9e; rpi = %8.1e; rdi = %8.1e; gap ="
               " %8.1e\n", csa->iter, csa->obj, csa->rpi, csa->rdi,
               csa->gap);
         /* check if the current point is optimal */
         if (csa->rpi < 1e-8 && csa->rdi < 1e-8 && csa->gap < 1e-8)
         {  if (csa->parm->msg_lev >= GLP_MSG_ALL)
               xprintf("OPTIMAL SOLUTION FOUND\n");
            status = 0;
            break;
         }
         /* check if the problem has no feasible solution */
         temp = 1e5 * csa->phi_min[csa->iter];
         if (temp < 1e-8) temp = 1e-8;
         if (csa->phi >= temp)
         {  if (csa->parm->msg_lev >= GLP_MSG_ALL)
               xprintf("PROBLEM HAS NO FEASIBLE PRIMAL/DUAL SOLUTION\n")
                  ;
            status = 1;
            break;
         }
         /* check for very slow convergence or divergence */
         if (((csa->rpi >= 1e-8 || csa->rdi >= 1e-8) && csa->rmu /
               csa->rmu0 >= 1e6) ||
               (csa->iter >= 30 && csa->phi_min[csa->iter] >= 0.5 *
               csa->phi_min[csa->iter - 30]))
         {  if (csa->parm->msg_lev >= GLP_MSG_ALL)
               xprintf("NO CONVERGENCE; SEARCH TERMINATED\n");
            status = 2;
            break;
         }
         /* check for maximal number of iterations */
         if (csa->iter == ITER_MAX)
         {  if (csa->parm->msg_lev >= GLP_MSG_ALL)
               xprintf("ITERATION LIMIT EXCEEDED; SEARCH TERMINATED\n");
            status = 3;
            break;
         }
         /* start the next iteration */
         csa->iter++;
         /* factorize normal equation system */
         for (j = 1; j <= n; j++) csa->D[j] = csa->x[j] / csa->z[j];
         decomp_NE(csa);
         /* compute the next point using Mehrotra's predictor-corrector
            technique */
         if (make_step(csa))
         {  if (csa->parm->msg_lev >= GLP_MSG_ALL)
               xprintf("NUMERIC INSTABILITY; SEARCH TERMINATED\n");
            status = 4;
            break;
         }
      }
      /* restore the best point */
      if (status != 0)
      {  for (j = 1; j <= n; j++) csa->x[j] = csa->best_x[j];
         for (i = 1; i <= m; i++) csa->y[i] = csa->best_y[i];
         for (j = 1; j <= n; j++) csa->z[j] = csa->best_z[j];
         if (csa->parm->msg_lev >= GLP_MSG_ALL)
            xprintf("Best point %17.9e was reached on iteration %d\n",
               csa->best_obj, csa->best_iter);
      }
      /* return to the calling program */
      return status;
}

/***********************************************************************
*  NAME
*
*  ipm_solve - core LP solver based on the interior-point method
*
*  SYNOPSIS
*
*  #include "glpipm.h"
*  int ipm_solve(glp_prob *P, const glp_iptcp *parm);
*
*  DESCRIPTION
*
*  The routine ipm_solve is a core LP solver based on the primal-dual
*  interior-point method.
*
*  The routine assumes the following standard formulation of LP problem
*  to be solved:
*
*     minimize
*
*        F = c[0] + c[1]*x[1] + c[2]*x[2] + ... + c[n]*x[n]
*
*     subject to linear constraints
*
*        a[1,1]*x[1] + a[1,2]*x[2] + ... + a[1,n]*x[n] = b[1]
*
*        a[2,1]*x[1] + a[2,2]*x[2] + ... + a[2,n]*x[n] = b[2]
*
*              . . . . . .
*
*        a[m,1]*x[1] + a[m,2]*x[2] + ... + a[m,n]*x[n] = b[m]
*
*     and non-negative variables
*
*        x[1] >= 0, x[2] >= 0, ..., x[n] >= 0
*
*  where:
*  F                    is the objective function;
*  x[1], ..., x[n]      are (structural) variables;
*  c[0]                 is a constant term of the objective function;
*  c[1], ..., c[n]      are objective coefficients;
*  a[1,1], ..., a[m,n]  are constraint coefficients;
*  b[1], ..., b[n]      are right-hand sides.
*
*  The solution is three vectors x, y, and z, which are stored by the
*  routine in the arrays x, y, and z, respectively. These vectors
*  correspond to the best primal-dual point found during optimization.
*  They are approximate solution of the following system (which is the
*  Karush-Kuhn-Tucker optimality conditions):
*
*     A*x      = b      (primal feasibility condition)
*
*     A'*y + z = c      (dual feasibility condition)
*
*     x'*z     = 0      (primal-dual complementarity condition)
*
*     x >= 0, z >= 0    (non-negativity condition)
*
*  where:
*  x[1], ..., x[n]      are primal (structural) variables;
*  y[1], ..., y[m]      are dual variables (Lagrange multipliers) for
*                       equality constraints;
*  z[1], ..., z[n]      are dual variables (Lagrange multipliers) for
*                       non-negativity constraints.
*
*  RETURNS
*
*  0  LP has been successfully solved.
*
*  GLP_ENOCVG
*     No convergence.
*
*  GLP_EITLIM
*     Iteration limit exceeded.
*
*  GLP_EINSTAB
*     Numeric instability on solving Newtonian system.
*
*  In case of non-zero return code the routine returns the best point,
*  which has been reached during optimization. */

int ipm_solve(glp_prob *P, const glp_iptcp *parm)
{     struct csa _dsa, *csa = &_dsa;
      int m = P->m;
      int n = P->n;
      int nnz = P->nnz;
      GLPROW *row;
      GLPCOL *col;
      GLPAIJ *aij;
      int i, j, loc, ret, *A_ind, *A_ptr;
      double dir, *A_val, *b, *c, *x, *y, *z;
      xassert(m > 0);
      xassert(n > 0);
      /* allocate working arrays */
      A_ptr = xcalloc(1+m+1, sizeof(int));
      A_ind = xcalloc(1+nnz, sizeof(int));
      A_val = xcalloc(1+nnz, sizeof(double));
      b = xcalloc(1+m, sizeof(double));
      c = xcalloc(1+n, sizeof(double));
      x = xcalloc(1+n, sizeof(double));
      y = xcalloc(1+m, sizeof(double));
      z = xcalloc(1+n, sizeof(double));
      /* prepare rows and constraint coefficients */
      loc = 1;
      for (i = 1; i <= m; i++)
      {  row = P->row[i];
         xassert(row->type == GLP_FX);
         b[i] = row->lb * row->rii;
         A_ptr[i] = loc;
         for (aij = row->ptr; aij != NULL; aij = aij->r_next)
         {  A_ind[loc] = aij->col->j;
            A_val[loc] = row->rii * aij->val * aij->col->sjj;
            loc++;
         }
      }
      A_ptr[m+1] = loc;
      xassert(loc-1 == nnz);
      /* prepare columns and objective coefficients */
      if (P->dir == GLP_MIN)
         dir = +1.0;
      else if (P->dir == GLP_MAX)
         dir = -1.0;
      else
         xassert(P != P);
      c[0] = dir * P->c0;
      for (j = 1; j <= n; j++)
      {  col = P->col[j];
         xassert(col->type == GLP_LO && col->lb == 0.0);
         c[j] = dir * col->coef * col->sjj;
      }
      /* allocate and initialize the common storage area */
      csa->m = m;
      csa->n = n;
      csa->A_ptr = A_ptr;
      csa->A_ind = A_ind;
      csa->A_val = A_val;
      csa->b = b;
      csa->c = c;
      csa->x = x;
      csa->y = y;
      csa->z = z;
      csa->parm = parm;
      initialize(csa);
      /* solve LP with the interior-point method */
      ret = ipm_main(csa);
      /* deallocate the common storage area */
      terminate(csa);
      /* determine solution status */
      if (ret == 0)
      {  /* optimal solution found */
         P->ipt_stat = GLP_OPT;
         ret = 0;
      }
      else if (ret == 1)
      {  /* problem has no feasible (primal or dual) solution */
         P->ipt_stat = GLP_NOFEAS;
         ret = 0;
      }
      else if (ret == 2)
      {  /* no convergence */
         P->ipt_stat = GLP_INFEAS;
         ret = GLP_ENOCVG;
      }
      else if (ret == 3)
      {  /* iteration limit exceeded */
         P->ipt_stat = GLP_INFEAS;
         ret = GLP_EITLIM;
      }
      else if (ret == 4)
      {  /* numeric instability on solving Newtonian system */
         P->ipt_stat = GLP_INFEAS;
         ret = GLP_EINSTAB;
      }
      else
         xassert(ret != ret);
      /* store row solution components */
      for (i = 1; i <= m; i++)
      {  row = P->row[i];
         row->pval = row->lb;
         row->dval = dir * y[i] * row->rii;
      }
      /* store column solution components */
      P->ipt_obj = P->c0;
      for (j = 1; j <= n; j++)
      {  col = P->col[j];
         col->pval = x[j] * col->sjj;
         col->dval = dir * z[j] / col->sjj;
         P->ipt_obj += col->coef * col->pval;
      }
      /* free working arrays */
      xfree(A_ptr);
      xfree(A_ind);
      xfree(A_val);
      xfree(b);
      xfree(c);
      xfree(x);
      xfree(y);
      xfree(z);
      return ret;
}

/* eof */
