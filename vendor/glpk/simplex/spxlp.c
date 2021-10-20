/* spxlp.c */

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

#include "env.h"
#include "spxlp.h"

/***********************************************************************
*  spx_factorize - compute factorization of current basis matrix
*
*  This routine computes factorization of the current basis matrix B.
*
*  If the factorization has been successfully computed, the routine
*  validates it and returns zero. Otherwise, the routine invalidates
*  the factorization and returns the code provided by the factorization
*  driver (bfd_factorize). */

static int jth_col(void *info, int j, int ind[], double val[])
{     /* provide column B[j] */
      SPXLP *lp = info;
      int m = lp->m;
      int *A_ptr = lp->A_ptr;
      int *head = lp->head;
      int k, ptr, len;
      xassert(1 <= j && j <= m);
      k = head[j]; /* x[k] = xB[j] */
      ptr = A_ptr[k];
      len = A_ptr[k+1] - ptr;
      memcpy(&ind[1], &lp->A_ind[ptr], len * sizeof(int));
      memcpy(&val[1], &lp->A_val[ptr], len * sizeof(double));
      return len;
}

int spx_factorize(SPXLP *lp)
{     int ret;
      ret = bfd_factorize(lp->bfd, lp->m, jth_col, lp);
      lp->valid = (ret == 0);
      return ret;
}

/***********************************************************************
*  spx_eval_beta - compute current values of basic variables
*
*  This routine computes vector beta = (beta[i]) of current values of
*  basic variables xB = (xB[i]). (Factorization of the current basis
*  matrix should be valid.)
*
*  First the routine computes a modified vector of right-hand sides:
*
*                         n-m
*     y = b - N * f = b - sum N[j] * f[j],
*                         j=1
*
*  where b = (b[i]) is the original vector of right-hand sides, N is
*  a matrix composed from columns of the original constraint matrix A,
*  which (columns) correspond to non-basic variables, f = (f[j]) is the
*  vector of active bounds of non-basic variables xN = (xN[j]),
*  N[j] = A[k] is a column of matrix A corresponding to non-basic
*  variable xN[j] = x[k], f[j] is current active bound lN[j] = l[k] or
*  uN[j] = u[k] of non-basic variable xN[j] = x[k]. The matrix-vector
*  product N * f is computed as a linear combination of columns of N,
*  so if f[j] = 0, column N[j] can be skipped.
*
*  Then the routine performs FTRAN to compute the vector beta:
*
*     beta = inv(B) * y.
*
*  On exit the routine stores components of the vector beta to array
*  locations beta[1], ..., beta[m]. */

void spx_eval_beta(SPXLP *lp, double beta[/*1+m*/])
{     int m = lp->m;
      int n = lp->n;
      int *A_ptr = lp->A_ptr;
      int *A_ind = lp->A_ind;
      double *A_val = lp->A_val;
      double *b = lp->b;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      int j, k, ptr, end;
      double fj, *y;
      /* compute y = b - N * xN */
      /* y := b */
      y = beta;
      memcpy(&y[1], &b[1], m * sizeof(double));
      /* y := y - N * f */
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         /* f[j] := active bound of xN[j] */
         fj = flag[j] ? u[k] : l[k];
         if (fj == 0.0 || fj == -DBL_MAX)
         {  /* either xN[j] has zero active bound or it is unbounded;
             * in the latter case its value is assumed to be zero */
            continue;
         }
         /* y := y - N[j] * f[j] */
         ptr = A_ptr[k];
         end = A_ptr[k+1];
         for (; ptr < end; ptr++)
            y[A_ind[ptr]] -= A_val[ptr] * fj;
      }
      /* compute beta = inv(B) * y */
      xassert(lp->valid);
      bfd_ftran(lp->bfd, beta);
      return;
}

/***********************************************************************
*  spx_eval_obj - compute current value of objective function
*
*  This routine computes the value of the objective function in the
*  current basic solution:
*
*     z = cB'* beta + cN'* f + c[0] =
*
*          m                    n-m
*       = sum cB[i] * beta[i] + sum cN[j] * f[j] + c[0],
*         i=1                   j=1
*
*  where cB = (cB[i]) is the vector of objective coefficients at basic
*  variables, beta = (beta[i]) is the vector of current values of basic
*  variables, cN = (cN[j]) is the vector of objective coefficients at
*  non-basic variables, f = (f[j]) is the vector of current active
*  bounds of non-basic variables, c[0] is the constant term of the
*  objective function.
*
*  It as assumed that components of the vector beta are stored in the
*  array locations beta[1], ..., beta[m]. */

double spx_eval_obj(SPXLP *lp, const double beta[/*1+m*/])
{     int m = lp->m;
      int n = lp->n;
      double *c = lp->c;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      int i, j, k;
      double fj, z;
      /* compute z = cB'* beta + cN'* f + c0 */
      /* z := c0 */
      z = c[0];
      /* z := z + cB'* beta */
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
         z += c[k] * beta[i];
      }
      /* z := z + cN'* f */
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         /* f[j] := active bound of xN[j] */
         fj = flag[j] ? u[k] : l[k];
         if (fj == 0.0 || fj == -DBL_MAX)
         {  /* either xN[j] has zero active bound or it is unbounded;
             * in the latter case its value is assumed to be zero */
            continue;
         }
         z += c[k] * fj;
      }
      return z;
}

/***********************************************************************
*  spx_eval_pi - compute simplex multipliers in current basis
*
*  This routine computes vector pi = (pi[i]) of simplex multipliers in
*  the current basis. (Factorization of the current basis matrix should
*  be valid.)
*
*  The vector pi is computed by performing BTRAN:
*
*     pi = inv(B') * cB,
*
*  where cB = (cB[i]) is the vector of objective coefficients at basic
*  variables xB = (xB[i]).
*
*  On exit components of vector pi are stored in the array locations
*  pi[1], ..., pi[m]. */

void spx_eval_pi(SPXLP *lp, double pi[/*1+m*/])
{     int m = lp->m;
      double *c = lp->c;
      int *head = lp->head;
      int i;
      double *cB;
      /* construct cB */
      cB = pi;
      for (i = 1; i <= m; i++)
         cB[i] = c[head[i]];
      /* compute pi = inv(B) * cB */
      bfd_btran(lp->bfd, pi);
      return;
}

/***********************************************************************
*  spx_eval_dj - compute reduced cost of j-th non-basic variable
*
*  This routine computes reduced cost d[j] of non-basic variable
*  xN[j] = x[k], 1 <= j <= n-m, in the current basic solution:
*
*     d[j] = c[k] - A'[k] * pi,
*
*  where c[k] is the objective coefficient at x[k], A[k] is k-th column
*  of the constraint matrix, pi is the vector of simplex multipliers in
*  the current basis.
*
*  It as assumed that components of the vector pi are stored in the
*  array locations pi[1], ..., pi[m]. */

double spx_eval_dj(SPXLP *lp, const double pi[/*1+m*/], int j)
{     int m = lp->m;
      int n = lp->n;
      int *A_ptr = lp->A_ptr;
      int *A_ind = lp->A_ind;
      double *A_val = lp->A_val;
      int k, ptr, end;
      double dj;
      xassert(1 <= j && j <= n-m);
      k = lp->head[m+j]; /* x[k] = xN[j] */
      /* dj := c[k] */
      dj = lp->c[k];
      /* dj := dj - A'[k] * pi */
      ptr = A_ptr[k];
      end = A_ptr[k+1];
      for (; ptr < end; ptr++)
         dj -= A_val[ptr] * pi[A_ind[ptr]];
      return dj;
}

/***********************************************************************
*  spx_eval_tcol - compute j-th column of simplex table
*
*  This routine computes j-th column of the current simplex table
*  T = (T[i,j]) = - inv(B) * N, 1 <= j <= n-m. (Factorization of the
*  current basis matrix should be valid.)
*
*  The simplex table column is computed by performing FTRAN:
*
*     tcol = - inv(B) * N[j],
*
*  where B is the current basis matrix, N[j] = A[k] is a column of the
*  constraint matrix corresponding to non-basic variable xN[j] = x[k].
*
*  On exit components of the simplex table column are stored in the
*  array locations tcol[1], ... tcol[m]. */

void spx_eval_tcol(SPXLP *lp, int j, double tcol[/*1+m*/])
{     int m = lp->m;
      int n = lp->n;
      int *A_ptr = lp->A_ptr;
      int *A_ind = lp->A_ind;
      double *A_val = lp->A_val;
      int *head = lp->head;
      int i, k, ptr, end;
      xassert(1 <= j && j <= n-m);
      k = head[m+j]; /* x[k] = xN[j] */
      /* compute tcol = - inv(B) * N[j] */
      for (i = 1; i <= m; i++)
         tcol[i] = 0.0;
      ptr = A_ptr[k];
      end = A_ptr[k+1];
      for (; ptr < end; ptr++)
         tcol[A_ind[ptr]] = -A_val[ptr];
      bfd_ftran(lp->bfd, tcol);
      return;
}

/***********************************************************************
*  spx_eval_rho - compute i-th row of basis matrix inverse
*
*  This routine computes i-th row of the matrix inv(B), where B is
*  the current basis matrix, 1 <= i <= m. (Factorization of the current
*  basis matrix should be valid.)
*
*  The inverse row is computed by performing BTRAN:
*
*     rho = inv(B') * e[i],
*
*  where e[i] is i-th column of unity matrix.
*
*  On exit components of the row are stored in the array locations
*  row[1], ..., row[m]. */

void spx_eval_rho(SPXLP *lp, int i, double rho[/*1+m*/])
{     int m = lp->m;
      int j;
      xassert(1 <= i && i <= m);
      /* compute rho = inv(B') * e[i] */
      for (j = 1; j <= m; j++)
         rho[j] = 0.0;
      rho[i] = 1.0;
      bfd_btran(lp->bfd, rho);
      return;
}

#if 1 /* 31/III-2016 */
void spx_eval_rho_s(SPXLP *lp, int i, FVS *rho)
{     /* sparse version of spx_eval_rho */
      int m = lp->m;
      xassert(1 <= i && i <= m);
      /* compute rho = inv(B') * e[i] */
      xassert(rho->n == m);
      fvs_clear_vec(rho);
      rho->nnz = 1;
      rho->ind[1] = i;
      rho->vec[i] = 1.0;
      bfd_btran_s(lp->bfd, rho);
      return;
}
#endif

/***********************************************************************
*  spx_eval_tij - compute element T[i,j] of simplex table
*
*  This routine computes element T[i,j] of the current simplex table
*  T = - inv(B) * N, 1 <= i <= m, 1 <= j <= n-m, with the following
*  formula:
*
*     T[i,j] = - N'[j] * rho,                                        (1)
*
*  where N[j] = A[k] is a column of the constraint matrix corresponding
*  to non-basic variable xN[j] = x[k], rho is i-th row of the inverse
*  matrix inv(B).
*
*  It as assumed that components of the inverse row rho = (rho[j]) are
*  stored in the array locations rho[1], ..., rho[m]. */

double spx_eval_tij(SPXLP *lp, const double rho[/*1+m*/], int j)
{     int m = lp->m;
      int n = lp->n;
      int *A_ptr = lp->A_ptr;
      int *A_ind = lp->A_ind;
      double *A_val = lp->A_val;
      int k, ptr, end;
      double tij;
      xassert(1 <= j && j <= n-m);
      k = lp->head[m+j]; /* x[k] = xN[j] */
      /* compute t[i,j] = - N'[j] * pi */
      tij = 0.0;
      ptr = A_ptr[k];
      end = A_ptr[k+1];
      for (; ptr < end; ptr++)
         tij -= A_val[ptr] * rho[A_ind[ptr]];
      return tij;
}

/***********************************************************************
*  spx_eval_trow - compute i-th row of simplex table
*
*  This routine computes i-th row of the current simplex table
*  T = (T[i,j]) = - inv(B) * N, 1 <= i <= m.
*
*  Elements of the row T[i] = (T[i,j]), j = 1, ..., n-m, are computed
*  directly with the routine spx_eval_tij.
*
*  The vector rho = (rho[j]), which is i-th row of the basis inverse
*  inv(B), should be previously computed with the routine spx_eval_rho.
*  It is assumed that elements of this vector are stored in the array
*  locations rho[1], ..., rho[m].
*
*  On exit components of the simplex table row are stored in the array
*  locations trow[1], ... trow[n-m].
*
*  NOTE: For testing/debugging only. */

void spx_eval_trow(SPXLP *lp, const double rho[/*1+m*/], double
      trow[/*1+n-m*/])
{     int m = lp->m;
      int n = lp->n;
      int j;
      for (j = 1; j <= n-m; j++)
         trow[j] = spx_eval_tij(lp, rho, j);
      return;
}

/***********************************************************************
*  spx_update_beta - update values of basic variables
*
*  This routine updates the vector beta = (beta[i]) of values of basic
*  variables xB = (xB[i]) for the adjacent basis.
*
*  On entry to the routine components of the vector beta in the current
*  basis should be placed in array locations beta[1], ..., beta[m].
*
*  The parameter 1 <= p <= m specifies basic variable xB[p] which
*  becomes non-basic variable xN[q] in the adjacent basis. The special
*  case p < 0 means that non-basic variable xN[q] goes from its current
*  active bound to opposite one in the adjacent basis.
*
*  If the flag p_flag is set, the active bound of xB[p] in the adjacent
*  basis is set to its upper bound. (In this case xB[p] should have its
*  upper bound and should not be fixed.)
*
*  The parameter 1 <= q <= n-m specifies non-basic variable xN[q] which
*  becomes basic variable xB[p] in the adjacent basis (if 1 <= p <= m),
*  or goes to its opposite bound (if p < 0). (In the latter case xN[q]
*  should have both lower and upper bounds and should not be fixed.)
*
*  It is assumed that the array tcol contains elements of q-th (pivot)
*  column T[q] of the simple table in locations tcol[1], ..., tcol[m].
*  (This column should be computed for the current basis.)
*
*  First, the routine determines the increment of basic variable xB[p]
*  in the adjacent basis (but only if 1 <= p <= m):
*
*                   (       - beta[p], if -inf < xB[p] < +inf
*                   (
*     delta xB[p] = { lB[p] - beta[p], if p_flag = 0
*                   (
*                   ( uB[p] - beta[p], if p_flag = 1
*
*  where beta[p] is the value of xB[p] in the current basis, lB[p] and
*  uB[p] are its lower and upper bounds. Then, the routine determines
*  the increment of non-basic variable xN[q] in the adjacent basis:
*
*                   ( delta xB[p] / T[p,q], if 1 <= p <= m
*                   (
*     delta xN[q] = { uN[q] - lN[q],        if p < 0 and f[q] = lN[q]
*                   (
*                   ( lN[q] - uN[q],        if p < 0 and f[q] = uN[q]
*
*  where T[p,q] is the pivot element of the simplex table, f[q] is the
*  active bound of xN[q] in the current basis.
*
*  If 1 <= p <= m, in the adjacent basis xN[q] becomes xB[p], so:
*
*     new beta[p] = f[q] + delta xN[q].
*
*  Values of other basic variables xB[i] for 1 <= i <= m, i != p, are
*  updated as follows:
*
*     new beta[i] = beta[i] + T[i,q] * delta xN[q].
*
*  On exit the routine stores updated components of the vector beta to
*  the same locations, where the input vector beta was stored. */

void spx_update_beta(SPXLP *lp, double beta[/*1+m*/], int p,
      int p_flag, int q, const double tcol[/*1+m*/])
{     int m = lp->m;
      int n = lp->n;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      int i, k;
      double delta_p, delta_q;
      if (p < 0)
      {  /* special case: xN[q] goes to its opposite bound */
         xassert(1 <= q && q <= n-m);
         /* xN[q] should be double-bounded variable */
         k = head[m+q]; /* x[k] = xN[q] */
         xassert(l[k] != -DBL_MAX && u[k] != +DBL_MAX && l[k] != u[k]);
         /* determine delta xN[q] */
         if (flag[q])
         {  /* xN[q] goes from its upper bound to its lower bound */
            delta_q = l[k] - u[k];
         }
         else
         {  /* xN[q] goes from its lower bound to its upper bound */
            delta_q = u[k] - l[k];
         }
      }
      else
      {  /* xB[p] leaves the basis, xN[q] enters the basis */
         xassert(1 <= p && p <= m);
         xassert(1 <= q && q <= n-m);
         /* determine delta xB[p] */
         k = head[p]; /* x[k] = xB[p] */
         if (p_flag)
         {  /* xB[p] goes to its upper bound */
            xassert(l[k] != u[k] && u[k] != +DBL_MAX);
            delta_p = u[k] - beta[p];
         }
         else if (l[k] == -DBL_MAX)
         {  /* unbounded xB[p] becomes non-basic (unusual case) */
            xassert(u[k] == +DBL_MAX);
            delta_p = 0.0 - beta[p];
         }
         else
         {  /* xB[p] goes to its lower bound or becomes fixed */
            delta_p = l[k] - beta[p];
         }
         /* determine delta xN[q] */
         delta_q = delta_p / tcol[p];
         /* compute new beta[p], which is the value of xN[q] in the
          * adjacent basis */
         k = head[m+q]; /* x[k] = xN[q] */
         if (flag[q])
         {  /* xN[q] has its upper bound active */
            xassert(l[k] != u[k] && u[k] != +DBL_MAX);
            beta[p] = u[k] + delta_q;
         }
         else if (l[k] == -DBL_MAX)
         {  /* xN[q] is non-basic unbounded variable */
            xassert(u[k] == +DBL_MAX);
            beta[p] = 0.0 + delta_q;
         }
         else
         {  /* xN[q] has its lower bound active or is fixed (latter
             * case is unusual) */
            beta[p] = l[k] + delta_q;
         }
      }
      /* compute new beta[i] for all i != p */
      for (i = 1; i <= m; i++)
      {  if (i != p)
            beta[i] += tcol[i] * delta_q;
      }
      return;
}

#if 1 /* 30/III-2016 */
void spx_update_beta_s(SPXLP *lp, double beta[/*1+m*/], int p,
      int p_flag, int q, const FVS *tcol)
{     /* sparse version of spx_update_beta */
      int m = lp->m;
      int n = lp->n;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      int nnz = tcol->nnz;
      int *ind = tcol->ind;
      double *vec = tcol->vec;
      int i, k;
      double delta_p, delta_q;
      xassert(tcol->n == m);
      if (p < 0)
      {  /* special case: xN[q] goes to its opposite bound */
#if 0 /* 11/VI-2017 */
         /* FIXME: not tested yet */
         xassert(0);
#endif
         xassert(1 <= q && q <= n-m);
         /* xN[q] should be double-bounded variable */
         k = head[m+q]; /* x[k] = xN[q] */
         xassert(l[k] != -DBL_MAX && u[k] != +DBL_MAX && l[k] != u[k]);
         /* determine delta xN[q] */
         if (flag[q])
         {  /* xN[q] goes from its upper bound to its lower bound */
            delta_q = l[k] - u[k];
         }
         else
         {  /* xN[q] goes from its lower bound to its upper bound */
            delta_q = u[k] - l[k];
         }
      }
      else
      {  /* xB[p] leaves the basis, xN[q] enters the basis */
         xassert(1 <= p && p <= m);
         xassert(1 <= q && q <= n-m);
         /* determine delta xB[p] */
         k = head[p]; /* x[k] = xB[p] */
         if (p_flag)
         {  /* xB[p] goes to its upper bound */
            xassert(l[k] != u[k] && u[k] != +DBL_MAX);
            delta_p = u[k] - beta[p];
         }
         else if (l[k] == -DBL_MAX)
         {  /* unbounded xB[p] becomes non-basic (unusual case) */
            xassert(u[k] == +DBL_MAX);
            delta_p = 0.0 - beta[p];
         }
         else
         {  /* xB[p] goes to its lower bound or becomes fixed */
            delta_p = l[k] - beta[p];
         }
         /* determine delta xN[q] */
         delta_q = delta_p / vec[p];
         /* compute new beta[p], which is the value of xN[q] in the
          * adjacent basis */
         k = head[m+q]; /* x[k] = xN[q] */
         if (flag[q])
         {  /* xN[q] has its upper bound active */
            xassert(l[k] != u[k] && u[k] != +DBL_MAX);
            beta[p] = u[k] + delta_q;
         }
         else if (l[k] == -DBL_MAX)
         {  /* xN[q] is non-basic unbounded variable */
            xassert(u[k] == +DBL_MAX);
            beta[p] = 0.0 + delta_q;
         }
         else
         {  /* xN[q] has its lower bound active or is fixed (latter
             * case is unusual) */
            beta[p] = l[k] + delta_q;
         }
      }
      /* compute new beta[i] for all i != p */
      for (k = 1; k <= nnz; k++)
      {  i = ind[k];
         if (i != p)
            beta[i] += vec[i] * delta_q;
      }
      return;
}
#endif

/***********************************************************************
*  spx_update_d - update reduced costs of non-basic variables
*
*  This routine updates the vector d = (d[j]) of reduced costs of
*  non-basic variables xN = (xN[j]) for the adjacent basis.
*
*  On entry to the routine components of the vector d in the current
*  basis should be placed in locations d[1], ..., d[n-m].
*
*  The parameter 1 <= p <= m specifies basic variable xB[p] which
*  becomes non-basic variable xN[q] in the adjacent basis.
*
*  The parameter 1 <= q <= n-m specified non-basic variable xN[q] which
*  becomes basic variable xB[p] in the adjacent basis.
*
*  It is assumed that the array trow contains elements of p-th (pivot)
*  row T'[p] of the simplex table in locations trow[1], ..., trow[n-m].
*  It is also assumed that the array tcol contains elements of q-th
*  (pivot) column T[q] of the simple table in locations tcol[1], ...,
*  tcol[m]. (These row and column should be computed for the current
*  basis.)
*
*  First, the routine computes more accurate reduced cost d[q] in the
*  current basis using q-th column of the simplex table:
*
*                    n-m
*     d[q] = cN[q] + sum t[i,q] * cB[i],
*                    i=1
*
*  where cN[q] and cB[i] are objective coefficients at variables xN[q]
*  and xB[i], resp. The routine also computes the relative error:
*
*     e = |d[q] - d'[q]| / (1 + |d[q]|),
*
*  where d'[q] is the reduced cost of xN[q] on entry to the routine,
*  and returns e on exit. (If e happens to be large enough, the calling
*  program may compute the reduced costs directly, since other reduced
*  costs also may be inaccurate.)
*
*  In the adjacent basis xB[p] becomes xN[q], so:
*
*     new d[q] = d[q] / T[p,q],
*
*  where T[p,q] is the pivot element of the simplex table (it is taken
*  from column T[q] as more accurate). Reduced costs of other non-basic
*  variables xN[j] for 1 <= j <= n-m, j != q, are updated as follows:
*
*     new d[j] = d[j] + T[p,j] * new d[q].
*
*  On exit the routine stores updated components of the vector d to the
*  same locations, where the input vector d was stored. */

double spx_update_d(SPXLP *lp, double d[/*1+n-m*/], int p, int q,
      const double trow[/*1+n-m*/], const double tcol[/*1+m*/])
{     int m = lp->m;
      int n = lp->n;
      double *c = lp->c;
      int *head = lp->head;
      int i, j, k;
      double dq, e;
      xassert(1 <= p && p <= m);
      xassert(1 <= q && q <= n);
      /* compute d[q] in current basis more accurately */
      k = head[m+q]; /* x[k] = xN[q] */
      dq = c[k];
      for (i = 1; i <= m; i++)
         dq += tcol[i] * c[head[i]];
      /* compute relative error in d[q] */
      e = fabs(dq - d[q]) / (1.0 + fabs(dq));
      /* compute new d[q], which is the reduced cost of xB[p] in the
       * adjacent basis */
      d[q] = (dq /= tcol[p]);
      /* compute new d[j] for all j != q */
      for (j = 1; j <= n-m; j++)
      {  if (j != q)
            d[j] -= trow[j] * dq;
      }
      return e;
}

#if 1 /* 30/III-2016 */
double spx_update_d_s(SPXLP *lp, double d[/*1+n-m*/], int p, int q,
      const FVS *trow, const FVS *tcol)
{     /* sparse version of spx_update_d */
      int m = lp->m;
      int n = lp->n;
      double *c = lp->c;
      int *head = lp->head;
      int trow_nnz = trow->nnz;
      int *trow_ind = trow->ind;
      double *trow_vec = trow->vec;
      int tcol_nnz = tcol->nnz;
      int *tcol_ind = tcol->ind;
      double *tcol_vec = tcol->vec;
      int i, j, k;
      double dq, e;
      xassert(1 <= p && p <= m);
      xassert(1 <= q && q <= n);
      xassert(trow->n == n-m);
      xassert(tcol->n == m);
      /* compute d[q] in current basis more accurately */
      k = head[m+q]; /* x[k] = xN[q] */
      dq = c[k];
      for (k = 1; k <= tcol_nnz; k++)
      {  i = tcol_ind[k];
         dq += tcol_vec[i] * c[head[i]];
      }
      /* compute relative error in d[q] */
      e = fabs(dq - d[q]) / (1.0 + fabs(dq));
      /* compute new d[q], which is the reduced cost of xB[p] in the
       * adjacent basis */
      d[q] = (dq /= tcol_vec[p]);
      /* compute new d[j] for all j != q */
      for (k = 1; k <= trow_nnz; k++)
      {  j = trow_ind[k];
         if (j != q)
            d[j] -= trow_vec[j] * dq;
      }
      return e;
}
#endif

/***********************************************************************
*  spx_change_basis - change current basis to adjacent one
*
*  This routine changes the current basis to the adjacent one making
*  necessary changes in lp->head and lp->flag members.
*
*  The parameters p, p_flag, and q have the same meaning as for the
*  routine spx_update_beta. */

void spx_change_basis(SPXLP *lp, int p, int p_flag, int q)
{     int m = lp->m;
      int n = lp->n;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      int k;
      if (p < 0)
      {  /* special case: xN[q] goes to its opposite bound */
         xassert(1 <= q && q <= n-m);
         /* xN[q] should be double-bounded variable */
         k = head[m+q]; /* x[k] = xN[q] */
         xassert(l[k] != -DBL_MAX && u[k] != +DBL_MAX && l[k] != u[k]);
         /* change active bound flag */
         flag[q] = 1 - flag[q];
      }
      else
      {  /* xB[p] leaves the basis, xN[q] enters the basis */
         xassert(1 <= p && p <= m);
         xassert(p_flag == 0 || p_flag == 1);
         xassert(1 <= q && q <= n-m);
         k = head[p]; /* xB[p] = x[k] */
         if (p_flag)
         {  /* xB[p] goes to its upper bound */
            xassert(l[k] != u[k] && u[k] != +DBL_MAX);
         }
         /* swap xB[p] and xN[q] in the basis */
         head[p] = head[m+q], head[m+q] = k;
         /* and set active bound flag for new xN[q] */
         lp->flag[q] = p_flag;
      }
      return;
}

/***********************************************************************
*  spx_update_invb - update factorization of basis matrix
*
*  This routine updates factorization of the basis matrix B when i-th
*  column of B is replaced by k-th column of the constraint matrix A.
*
*  The parameter 1 <= i <= m specifies the number of column of matrix B
*  to be replaced by a new column.
*
*  The parameter 1 <= k <= n specifies the number of column of matrix A
*  to be used for replacement.
*
*  If the factorization has been successfully updated, the routine
*  validates it and returns zero. Otherwise, the routine invalidates
*  the factorization and returns the code provided by the factorization
*  driver (bfd_update). */

int spx_update_invb(SPXLP *lp, int i, int k)
{     int m = lp->m;
      int n = lp->n;
      int *A_ptr = lp->A_ptr;
      int *A_ind = lp->A_ind;
      double *A_val = lp->A_val;
      int ptr, len, ret;
      xassert(1 <= i && i <= m);
      xassert(1 <= k && k <= n);
      ptr = A_ptr[k];
      len = A_ptr[k+1] - ptr;
      ret = bfd_update(lp->bfd, i, len, &A_ind[ptr-1], &A_val[ptr-1]);
      lp->valid = (ret == 0);
      return ret;
}

/* eof */
