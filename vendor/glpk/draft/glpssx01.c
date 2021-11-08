/* glpssx01.c (simplex method, rational arithmetic) */

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

#include "env.h"
#include "glpssx.h"
#define xfault xerror

/*----------------------------------------------------------------------
// ssx_create - create simplex solver workspace.
//
// This routine creates the workspace used by simplex solver routines,
// and returns a pointer to it.
//
// Parameters m, n, and nnz specify, respectively, the number of rows,
// columns, and non-zero constraint coefficients.
//
// This routine only allocates the memory for the workspace components,
// so the workspace needs to be saturated by data. */

SSX *ssx_create(int m, int n, int nnz)
{     SSX *ssx;
      int i, j, k;
      if (m < 1)
         xfault("ssx_create: m = %d; invalid number of rows\n", m);
      if (n < 1)
         xfault("ssx_create: n = %d; invalid number of columns\n", n);
      if (nnz < 0)
         xfault("ssx_create: nnz = %d; invalid number of non-zero const"
            "raint coefficients\n", nnz);
      ssx = xmalloc(sizeof(SSX));
      ssx->m = m;
      ssx->n = n;
      ssx->type = xcalloc(1+m+n, sizeof(int));
      ssx->lb = xcalloc(1+m+n, sizeof(mpq_t));
      for (k = 1; k <= m+n; k++) mpq_init(ssx->lb[k]);
      ssx->ub = xcalloc(1+m+n, sizeof(mpq_t));
      for (k = 1; k <= m+n; k++) mpq_init(ssx->ub[k]);
      ssx->coef = xcalloc(1+m+n, sizeof(mpq_t));
      for (k = 0; k <= m+n; k++) mpq_init(ssx->coef[k]);
      ssx->A_ptr = xcalloc(1+n+1, sizeof(int));
      ssx->A_ptr[n+1] = nnz+1;
      ssx->A_ind = xcalloc(1+nnz, sizeof(int));
      ssx->A_val = xcalloc(1+nnz, sizeof(mpq_t));
      for (k = 1; k <= nnz; k++) mpq_init(ssx->A_val[k]);
      ssx->stat = xcalloc(1+m+n, sizeof(int));
      ssx->Q_row = xcalloc(1+m+n, sizeof(int));
      ssx->Q_col = xcalloc(1+m+n, sizeof(int));
      ssx->binv = bfx_create_binv();
      ssx->bbar = xcalloc(1+m, sizeof(mpq_t));
      for (i = 0; i <= m; i++) mpq_init(ssx->bbar[i]);
      ssx->pi = xcalloc(1+m, sizeof(mpq_t));
      for (i = 1; i <= m; i++) mpq_init(ssx->pi[i]);
      ssx->cbar = xcalloc(1+n, sizeof(mpq_t));
      for (j = 1; j <= n; j++) mpq_init(ssx->cbar[j]);
      ssx->rho = xcalloc(1+m, sizeof(mpq_t));
      for (i = 1; i <= m; i++) mpq_init(ssx->rho[i]);
      ssx->ap = xcalloc(1+n, sizeof(mpq_t));
      for (j = 1; j <= n; j++) mpq_init(ssx->ap[j]);
      ssx->aq = xcalloc(1+m, sizeof(mpq_t));
      for (i = 1; i <= m; i++) mpq_init(ssx->aq[i]);
      mpq_init(ssx->delta);
      return ssx;
}

/*----------------------------------------------------------------------
// ssx_factorize - factorize the current basis matrix.
//
// This routine computes factorization of the current basis matrix B
// and returns the singularity flag. If the matrix B is non-singular,
// the flag is zero, otherwise non-zero. */

static int basis_col(void *info, int j, int ind[], mpq_t val[])
{     /* this auxiliary routine provides row indices and numeric values
         of non-zero elements in j-th column of the matrix B */
      SSX *ssx = info;
      int m = ssx->m;
      int n = ssx->n;
      int *A_ptr = ssx->A_ptr;
      int *A_ind = ssx->A_ind;
      mpq_t *A_val = ssx->A_val;
      int *Q_col = ssx->Q_col;
      int k, len, ptr;
      xassert(1 <= j && j <= m);
      k = Q_col[j]; /* x[k] = xB[j] */
      xassert(1 <= k && k <= m+n);
      /* j-th column of the matrix B is k-th column of the augmented
         constraint matrix (I | -A) */
      if (k <= m)
      {  /* it is a column of the unity matrix I */
         len = 1, ind[1] = k, mpq_set_si(val[1], 1, 1);
      }
      else
      {  /* it is a column of the original constraint matrix -A */
         len = 0;
         for (ptr = A_ptr[k-m]; ptr < A_ptr[k-m+1]; ptr++)
         {  len++;
            ind[len] = A_ind[ptr];
            mpq_neg(val[len], A_val[ptr]);
         }
      }
      return len;
}

int ssx_factorize(SSX *ssx)
{     int ret;
      ret = bfx_factorize(ssx->binv, ssx->m, basis_col, ssx);
      return ret;
}

/*----------------------------------------------------------------------
// ssx_get_xNj - determine value of non-basic variable.
//
// This routine determines the value of non-basic variable xN[j] in the
// current basic solution defined as follows:
//
//    0,             if xN[j] is free variable
//    lN[j],         if xN[j] is on its lower bound
//    uN[j],         if xN[j] is on its upper bound
//    lN[j] = uN[j], if xN[j] is fixed variable
//
// where lN[j] and uN[j] are lower and upper bounds of xN[j]. */

void ssx_get_xNj(SSX *ssx, int j, mpq_t x)
{     int m = ssx->m;
      int n = ssx->n;
      mpq_t *lb = ssx->lb;
      mpq_t *ub = ssx->ub;
      int *stat = ssx->stat;
      int *Q_col = ssx->Q_col;
      int k;
      xassert(1 <= j && j <= n);
      k = Q_col[m+j]; /* x[k] = xN[j] */
      xassert(1 <= k && k <= m+n);
      switch (stat[k])
      {  case SSX_NL:
            /* xN[j] is on its lower bound */
            mpq_set(x, lb[k]); break;
         case SSX_NU:
            /* xN[j] is on its upper bound */
            mpq_set(x, ub[k]); break;
         case SSX_NF:
            /* xN[j] is free variable */
            mpq_set_si(x, 0, 1); break;
         case SSX_NS:
            /* xN[j] is fixed variable */
            mpq_set(x, lb[k]); break;
         default:
            xassert(stat != stat);
      }
      return;
}

/*----------------------------------------------------------------------
// ssx_eval_bbar - compute values of basic variables.
//
// This routine computes values of basic variables xB in the current
// basic solution as follows:
//
//    beta = - inv(B) * N * xN,
//
// where B is the basis matrix, N is the matrix of non-basic columns,
// xN is a vector of current values of non-basic variables. */

void ssx_eval_bbar(SSX *ssx)
{     int m = ssx->m;
      int n = ssx->n;
      mpq_t *coef = ssx->coef;
      int *A_ptr = ssx->A_ptr;
      int *A_ind = ssx->A_ind;
      mpq_t *A_val = ssx->A_val;
      int *Q_col = ssx->Q_col;
      mpq_t *bbar = ssx->bbar;
      int i, j, k, ptr;
      mpq_t x, temp;
      mpq_init(x);
      mpq_init(temp);
      /* bbar := 0 */
      for (i = 1; i <= m; i++)
         mpq_set_si(bbar[i], 0, 1);
      /* bbar := - N * xN = - N[1] * xN[1] - ... - N[n] * xN[n] */
      for (j = 1; j <= n; j++)
      {  ssx_get_xNj(ssx, j, x);
         if (mpq_sgn(x) == 0) continue;
         k = Q_col[m+j]; /* x[k] = xN[j] */
         if (k <= m)
         {  /* N[j] is a column of the unity matrix I */
            mpq_sub(bbar[k], bbar[k], x);
         }
         else
         {  /* N[j] is a column of the original constraint matrix -A */
            for (ptr = A_ptr[k-m]; ptr < A_ptr[k-m+1]; ptr++)
            {  mpq_mul(temp, A_val[ptr], x);
               mpq_add(bbar[A_ind[ptr]], bbar[A_ind[ptr]], temp);
            }
         }
      }
      /* bbar := inv(B) * bbar */
      bfx_ftran(ssx->binv, bbar, 0);
#if 1
      /* compute value of the objective function */
      /* bbar[0] := c[0] */
      mpq_set(bbar[0], coef[0]);
      /* bbar[0] := bbar[0] + sum{i in B} cB[i] * xB[i] */
      for (i = 1; i <= m; i++)
      {  k = Q_col[i]; /* x[k] = xB[i] */
         if (mpq_sgn(coef[k]) == 0) continue;
         mpq_mul(temp, coef[k], bbar[i]);
         mpq_add(bbar[0], bbar[0], temp);
      }
      /* bbar[0] := bbar[0] + sum{j in N} cN[j] * xN[j] */
      for (j = 1; j <= n; j++)
      {  k = Q_col[m+j]; /* x[k] = xN[j] */
         if (mpq_sgn(coef[k]) == 0) continue;
         ssx_get_xNj(ssx, j, x);
         mpq_mul(temp, coef[k], x);
         mpq_add(bbar[0], bbar[0], temp);
      }
#endif
      mpq_clear(x);
      mpq_clear(temp);
      return;
}

/*----------------------------------------------------------------------
// ssx_eval_pi - compute values of simplex multipliers.
//
// This routine computes values of simplex multipliers (shadow prices)
// pi in the current basic solution as follows:
//
//    pi = inv(B') * cB,
//
// where B' is a matrix transposed to the basis matrix B, cB is a vector
// of objective coefficients at basic variables xB. */

void ssx_eval_pi(SSX *ssx)
{     int m = ssx->m;
      mpq_t *coef = ssx->coef;
      int *Q_col = ssx->Q_col;
      mpq_t *pi = ssx->pi;
      int i;
      /* pi := cB */
      for (i = 1; i <= m; i++) mpq_set(pi[i], coef[Q_col[i]]);
      /* pi := inv(B') * cB */
      bfx_btran(ssx->binv, pi);
      return;
}

/*----------------------------------------------------------------------
// ssx_eval_dj - compute reduced cost of non-basic variable.
//
// This routine computes reduced cost d[j] of non-basic variable xN[j]
// in the current basic solution as follows:
//
//    d[j] = cN[j] - N[j] * pi,
//
// where cN[j] is an objective coefficient at xN[j], N[j] is a column
// of the augmented constraint matrix (I | -A) corresponding to xN[j],
// pi is the vector of simplex multipliers (shadow prices). */

void ssx_eval_dj(SSX *ssx, int j, mpq_t dj)
{     int m = ssx->m;
      int n = ssx->n;
      mpq_t *coef = ssx->coef;
      int *A_ptr = ssx->A_ptr;
      int *A_ind = ssx->A_ind;
      mpq_t *A_val = ssx->A_val;
      int *Q_col = ssx->Q_col;
      mpq_t *pi = ssx->pi;
      int k, ptr, end;
      mpq_t temp;
      mpq_init(temp);
      xassert(1 <= j && j <= n);
      k = Q_col[m+j]; /* x[k] = xN[j] */
      xassert(1 <= k && k <= m+n);
      /* j-th column of the matrix N is k-th column of the augmented
         constraint matrix (I | -A) */
      if (k <= m)
      {  /* it is a column of the unity matrix I */
         mpq_sub(dj, coef[k], pi[k]);
      }
      else
      {  /* it is a column of the original constraint matrix -A */
         mpq_set(dj, coef[k]);
         for (ptr = A_ptr[k-m], end = A_ptr[k-m+1]; ptr < end; ptr++)
         {  mpq_mul(temp, A_val[ptr], pi[A_ind[ptr]]);
            mpq_add(dj, dj, temp);
         }
      }
      mpq_clear(temp);
      return;
}

/*----------------------------------------------------------------------
// ssx_eval_cbar - compute reduced costs of all non-basic variables.
//
// This routine computes the vector of reduced costs pi in the current
// basic solution for all non-basic variables, including fixed ones. */

void ssx_eval_cbar(SSX *ssx)
{     int n = ssx->n;
      mpq_t *cbar = ssx->cbar;
      int j;
      for (j = 1; j <= n; j++)
         ssx_eval_dj(ssx, j, cbar[j]);
      return;
}

/*----------------------------------------------------------------------
// ssx_eval_rho - compute p-th row of the inverse.
//
// This routine computes p-th row of the matrix inv(B), where B is the
// current basis matrix.
//
// p-th row of the inverse is computed using the following formula:
//
//    rho = inv(B') * e[p],
//
// where B' is a matrix transposed to B, e[p] is a unity vector, which
// contains one in p-th position. */

void ssx_eval_rho(SSX *ssx)
{     int m = ssx->m;
      int p = ssx->p;
      mpq_t *rho = ssx->rho;
      int i;
      xassert(1 <= p && p <= m);
      /* rho := 0 */
      for (i = 1; i <= m; i++) mpq_set_si(rho[i], 0, 1);
      /* rho := e[p] */
      mpq_set_si(rho[p], 1, 1);
      /* rho := inv(B') * rho */
      bfx_btran(ssx->binv, rho);
      return;
}

/*----------------------------------------------------------------------
// ssx_eval_row - compute pivot row of the simplex table.
//
// This routine computes p-th (pivot) row of the current simplex table
// A~ = - inv(B) * N using the following formula:
//
//    A~[p] = - N' * inv(B') * e[p] = - N' * rho[p],
//
// where N' is a matrix transposed to the matrix N, rho[p] is p-th row
// of the inverse inv(B). */

void ssx_eval_row(SSX *ssx)
{     int m = ssx->m;
      int n = ssx->n;
      int *A_ptr = ssx->A_ptr;
      int *A_ind = ssx->A_ind;
      mpq_t *A_val = ssx->A_val;
      int *Q_col = ssx->Q_col;
      mpq_t *rho = ssx->rho;
      mpq_t *ap = ssx->ap;
      int j, k, ptr;
      mpq_t temp;
      mpq_init(temp);
      for (j = 1; j <= n; j++)
      {  /* ap[j] := - N'[j] * rho (inner product) */
         k = Q_col[m+j]; /* x[k] = xN[j] */
         if (k <= m)
            mpq_neg(ap[j], rho[k]);
         else
         {  mpq_set_si(ap[j], 0, 1);
            for (ptr = A_ptr[k-m]; ptr < A_ptr[k-m+1]; ptr++)
            {  mpq_mul(temp, A_val[ptr], rho[A_ind[ptr]]);
               mpq_add(ap[j], ap[j], temp);
            }
         }
      }
      mpq_clear(temp);
      return;
}

/*----------------------------------------------------------------------
// ssx_eval_col - compute pivot column of the simplex table.
//
// This routine computes q-th (pivot) column of the current simplex
// table A~ = - inv(B) * N using the following formula:
//
//    A~[q] = - inv(B) * N[q],
//
// where N[q] is q-th column of the matrix N corresponding to chosen
// non-basic variable xN[q]. */

void ssx_eval_col(SSX *ssx)
{     int m = ssx->m;
      int n = ssx->n;
      int *A_ptr = ssx->A_ptr;
      int *A_ind = ssx->A_ind;
      mpq_t *A_val = ssx->A_val;
      int *Q_col = ssx->Q_col;
      int q = ssx->q;
      mpq_t *aq = ssx->aq;
      int i, k, ptr;
      xassert(1 <= q && q <= n);
      /* aq := 0 */
      for (i = 1; i <= m; i++) mpq_set_si(aq[i], 0, 1);
      /* aq := N[q] */
      k = Q_col[m+q]; /* x[k] = xN[q] */
      if (k <= m)
      {  /* N[q] is a column of the unity matrix I */
         mpq_set_si(aq[k], 1, 1);
      }
      else
      {  /* N[q] is a column of the original constraint matrix -A */
         for (ptr = A_ptr[k-m]; ptr < A_ptr[k-m+1]; ptr++)
            mpq_neg(aq[A_ind[ptr]], A_val[ptr]);
      }
      /* aq := inv(B) * aq */
      bfx_ftran(ssx->binv, aq, 1);
      /* aq := - aq */
      for (i = 1; i <= m; i++) mpq_neg(aq[i], aq[i]);
      return;
}

/*----------------------------------------------------------------------
// ssx_chuzc - choose pivot column.
//
// This routine chooses non-basic variable xN[q] whose reduced cost
// indicates possible improving of the objective function to enter it
// in the basis.
//
// Currently the standard (textbook) pricing is used, i.e. that
// non-basic variable is preferred which has greatest reduced cost (in
// magnitude).
//
// If xN[q] has been chosen, the routine stores its number q and also
// sets the flag q_dir that indicates direction in which xN[q] has to
// change (+1 means increasing, -1 means decreasing).
//
// If the choice cannot be made, because the current basic solution is
// dual feasible, the routine sets the number q to 0. */

void ssx_chuzc(SSX *ssx)
{     int m = ssx->m;
      int n = ssx->n;
      int dir = (ssx->dir == SSX_MIN ? +1 : -1);
      int *Q_col = ssx->Q_col;
      int *stat = ssx->stat;
      mpq_t *cbar = ssx->cbar;
      int j, k, s, q, q_dir;
      double best, temp;
      /* nothing is chosen so far */
      q = 0, q_dir = 0, best = 0.0;
      /* look through the list of non-basic variables */
      for (j = 1; j <= n; j++)
      {  k = Q_col[m+j]; /* x[k] = xN[j] */
         s = dir * mpq_sgn(cbar[j]);
         if ((stat[k] == SSX_NF || stat[k] == SSX_NL) && s < 0 ||
             (stat[k] == SSX_NF || stat[k] == SSX_NU) && s > 0)
         {  /* reduced cost of xN[j] indicates possible improving of
               the objective function */
            temp = fabs(mpq_get_d(cbar[j]));
            xassert(temp != 0.0);
            if (q == 0 || best < temp)
               q = j, q_dir = - s, best = temp;
         }
      }
      ssx->q = q, ssx->q_dir = q_dir;
      return;
}

/*----------------------------------------------------------------------
// ssx_chuzr - choose pivot row.
//
// This routine looks through elements of q-th column of the simplex
// table and chooses basic variable xB[p] which should leave the basis.
//
// The choice is based on the standard (textbook) ratio test.
//
// If xB[p] has been chosen, the routine stores its number p and also
// sets its non-basic status p_stat which should be assigned to xB[p]
// when it has left the basis and become xN[q].
//
// Special case p < 0 means that xN[q] is double-bounded variable and
// it reaches its opposite bound before any basic variable does that,
// so the current basis remains unchanged.
//
// If the choice cannot be made, because xN[q] can infinitely change in
// the feasible direction, the routine sets the number p to 0. */

void ssx_chuzr(SSX *ssx)
{     int m = ssx->m;
      int n = ssx->n;
      int *type = ssx->type;
      mpq_t *lb = ssx->lb;
      mpq_t *ub = ssx->ub;
      int *Q_col = ssx->Q_col;
      mpq_t *bbar = ssx->bbar;
      int q = ssx->q;
      mpq_t *aq = ssx->aq;
      int q_dir = ssx->q_dir;
      int i, k, s, t, p, p_stat;
      mpq_t teta, temp;
      mpq_init(teta);
      mpq_init(temp);
      xassert(1 <= q && q <= n);
      xassert(q_dir == +1 || q_dir == -1);
      /* nothing is chosen so far */
      p = 0, p_stat = 0;
      /* look through the list of basic variables */
      for (i = 1; i <= m; i++)
      {  s = q_dir * mpq_sgn(aq[i]);
         if (s < 0)
         {  /* xB[i] decreases */
            k = Q_col[i]; /* x[k] = xB[i] */
            t = type[k];
            if (t == SSX_LO || t == SSX_DB || t == SSX_FX)
            {  /* xB[i] has finite lower bound */
               mpq_sub(temp, bbar[i], lb[k]);
               mpq_div(temp, temp, aq[i]);
               mpq_abs(temp, temp);
               if (p == 0 || mpq_cmp(teta, temp) > 0)
               {  p = i;
                  p_stat = (t == SSX_FX ? SSX_NS : SSX_NL);
                  mpq_set(teta, temp);
               }
            }
         }
         else if (s > 0)
         {  /* xB[i] increases */
            k = Q_col[i]; /* x[k] = xB[i] */
            t = type[k];
            if (t == SSX_UP || t == SSX_DB || t == SSX_FX)
            {  /* xB[i] has finite upper bound */
               mpq_sub(temp, bbar[i], ub[k]);
               mpq_div(temp, temp, aq[i]);
               mpq_abs(temp, temp);
               if (p == 0 || mpq_cmp(teta, temp) > 0)
               {  p = i;
                  p_stat = (t == SSX_FX ? SSX_NS : SSX_NU);
                  mpq_set(teta, temp);
               }
            }
         }
         /* if something has been chosen and the ratio test indicates
            exact degeneracy, the search can be finished */
         if (p != 0 && mpq_sgn(teta) == 0) break;
      }
      /* if xN[q] is double-bounded, check if it can reach its opposite
         bound before any basic variable */
      k = Q_col[m+q]; /* x[k] = xN[q] */
      if (type[k] == SSX_DB)
      {  mpq_sub(temp, ub[k], lb[k]);
         if (p == 0 || mpq_cmp(teta, temp) > 0)
         {  p = -1;
            p_stat = -1;
            mpq_set(teta, temp);
         }
      }
      ssx->p = p;
      ssx->p_stat = p_stat;
      /* if xB[p] has been chosen, determine its actual change in the
         adjacent basis (it has the same sign as q_dir) */
      if (p != 0)
      {  xassert(mpq_sgn(teta) >= 0);
         if (q_dir > 0)
            mpq_set(ssx->delta, teta);
         else
            mpq_neg(ssx->delta, teta);
      }
      mpq_clear(teta);
      mpq_clear(temp);
      return;
}

/*----------------------------------------------------------------------
// ssx_update_bbar - update values of basic variables.
//
// This routine recomputes the current values of basic variables for
// the adjacent basis.
//
// The simplex table for the current basis is the following:
//
//    xB[i] = sum{j in 1..n} alfa[i,j] * xN[q],  i = 1,...,m
//
// therefore
//
//    delta xB[i] = alfa[i,q] * delta xN[q],  i = 1,...,m
//
// where delta xN[q] = xN.new[q] - xN[q] is the change of xN[q] in the
// adjacent basis, and delta xB[i] = xB.new[i] - xB[i] is the change of
// xB[i]. This gives formulae for recomputing values of xB[i]:
//
//    xB.new[p] = xN[q] + delta xN[q]
//
// (because xN[q] becomes xB[p] in the adjacent basis), and
//
//    xB.new[i] = xB[i] + alfa[i,q] * delta xN[q],  i != p
//
// for other basic variables. */

void ssx_update_bbar(SSX *ssx)
{     int m = ssx->m;
      int n = ssx->n;
      mpq_t *bbar = ssx->bbar;
      mpq_t *cbar = ssx->cbar;
      int p = ssx->p;
      int q = ssx->q;
      mpq_t *aq = ssx->aq;
      int i;
      mpq_t temp;
      mpq_init(temp);
      xassert(1 <= q && q <= n);
      if (p < 0)
      {  /* xN[q] is double-bounded and goes to its opposite bound */
         /* nop */;
      }
      else
      {  /* xN[q] becomes xB[p] in the adjacent basis */
         /* xB.new[p] = xN[q] + delta xN[q] */
         xassert(1 <= p && p <= m);
         ssx_get_xNj(ssx, q, temp);
         mpq_add(bbar[p], temp, ssx->delta);
      }
      /* update values of other basic variables depending on xN[q] */
      for (i = 1; i <= m; i++)
      {  if (i == p) continue;
         /* xB.new[i] = xB[i] + alfa[i,q] * delta xN[q] */
         if (mpq_sgn(aq[i]) == 0) continue;
         mpq_mul(temp, aq[i], ssx->delta);
         mpq_add(bbar[i], bbar[i], temp);
      }
#if 1
      /* update value of the objective function */
      /* z.new = z + d[q] * delta xN[q] */
      mpq_mul(temp, cbar[q], ssx->delta);
      mpq_add(bbar[0], bbar[0], temp);
#endif
      mpq_clear(temp);
      return;
}

/*----------------------------------------------------------------------
-- ssx_update_pi - update simplex multipliers.
--
-- This routine recomputes the vector of simplex multipliers for the
-- adjacent basis. */

void ssx_update_pi(SSX *ssx)
{     int m = ssx->m;
      int n = ssx->n;
      mpq_t *pi = ssx->pi;
      mpq_t *cbar = ssx->cbar;
      int p = ssx->p;
      int q = ssx->q;
      mpq_t *aq = ssx->aq;
      mpq_t *rho = ssx->rho;
      int i;
      mpq_t new_dq, temp;
      mpq_init(new_dq);
      mpq_init(temp);
      xassert(1 <= p && p <= m);
      xassert(1 <= q && q <= n);
      /* compute d[q] in the adjacent basis */
      mpq_div(new_dq, cbar[q], aq[p]);
      /* update the vector of simplex multipliers */
      for (i = 1; i <= m; i++)
      {  if (mpq_sgn(rho[i]) == 0) continue;
         mpq_mul(temp, new_dq, rho[i]);
         mpq_sub(pi[i], pi[i], temp);
      }
      mpq_clear(new_dq);
      mpq_clear(temp);
      return;
}

/*----------------------------------------------------------------------
// ssx_update_cbar - update reduced costs of non-basic variables.
//
// This routine recomputes the vector of reduced costs of non-basic
// variables for the adjacent basis. */

void ssx_update_cbar(SSX *ssx)
{     int m = ssx->m;
      int n = ssx->n;
      mpq_t *cbar = ssx->cbar;
      int p = ssx->p;
      int q = ssx->q;
      mpq_t *ap = ssx->ap;
      int j;
      mpq_t temp;
      mpq_init(temp);
      xassert(1 <= p && p <= m);
      xassert(1 <= q && q <= n);
      /* compute d[q] in the adjacent basis */
      /* d.new[q] = d[q] / alfa[p,q] */
      mpq_div(cbar[q], cbar[q], ap[q]);
      /* update reduced costs of other non-basic variables */
      for (j = 1; j <= n; j++)
      {  if (j == q) continue;
         /* d.new[j] = d[j] - (alfa[p,j] / alfa[p,q]) * d[q] */
         if (mpq_sgn(ap[j]) == 0) continue;
         mpq_mul(temp, ap[j], cbar[q]);
         mpq_sub(cbar[j], cbar[j], temp);
      }
      mpq_clear(temp);
      return;
}

/*----------------------------------------------------------------------
// ssx_change_basis - change current basis to adjacent one.
//
// This routine changes the current basis to the adjacent one swapping
// basic variable xB[p] and non-basic variable xN[q]. */

void ssx_change_basis(SSX *ssx)
{     int m = ssx->m;
      int n = ssx->n;
      int *type = ssx->type;
      int *stat = ssx->stat;
      int *Q_row = ssx->Q_row;
      int *Q_col = ssx->Q_col;
      int p = ssx->p;
      int q = ssx->q;
      int p_stat = ssx->p_stat;
      int k, kp, kq;
      if (p < 0)
      {  /* special case: xN[q] goes to its opposite bound */
         xassert(1 <= q && q <= n);
         k = Q_col[m+q]; /* x[k] = xN[q] */
         xassert(type[k] == SSX_DB);
         switch (stat[k])
         {  case SSX_NL:
               stat[k] = SSX_NU;
               break;
            case SSX_NU:
               stat[k] = SSX_NL;
               break;
            default:
               xassert(stat != stat);
         }
      }
      else
      {  /* xB[p] leaves the basis, xN[q] enters the basis */
         xassert(1 <= p && p <= m);
         xassert(1 <= q && q <= n);
         kp = Q_col[p];   /* x[kp] = xB[p] */
         kq = Q_col[m+q]; /* x[kq] = xN[q] */
         /* check non-basic status of xB[p] which becomes xN[q] */
         switch (type[kp])
         {  case SSX_FR:
               xassert(p_stat == SSX_NF);
               break;
            case SSX_LO:
               xassert(p_stat == SSX_NL);
               break;
            case SSX_UP:
               xassert(p_stat == SSX_NU);
               break;
            case SSX_DB:
               xassert(p_stat == SSX_NL || p_stat == SSX_NU);
               break;
            case SSX_FX:
               xassert(p_stat == SSX_NS);
               break;
            default:
               xassert(type != type);
         }
         /* swap xB[p] and xN[q] */
         stat[kp] = (char)p_stat, stat[kq] = SSX_BS;
         Q_row[kp] = m+q, Q_row[kq] = p;
         Q_col[p] = kq, Q_col[m+q] = kp;
         /* update factorization of the basis matrix */
         if (bfx_update(ssx->binv, p))
         {  if (ssx_factorize(ssx))
               xassert(("Internal error: basis matrix is singular", 0));
         }
      }
      return;
}

/*----------------------------------------------------------------------
// ssx_delete - delete simplex solver workspace.
//
// This routine deletes the simplex solver workspace freeing all the
// memory allocated to this object. */

void ssx_delete(SSX *ssx)
{     int m = ssx->m;
      int n = ssx->n;
      int nnz = ssx->A_ptr[n+1]-1;
      int i, j, k;
      xfree(ssx->type);
      for (k = 1; k <= m+n; k++) mpq_clear(ssx->lb[k]);
      xfree(ssx->lb);
      for (k = 1; k <= m+n; k++) mpq_clear(ssx->ub[k]);
      xfree(ssx->ub);
      for (k = 0; k <= m+n; k++) mpq_clear(ssx->coef[k]);
      xfree(ssx->coef);
      xfree(ssx->A_ptr);
      xfree(ssx->A_ind);
      for (k = 1; k <= nnz; k++) mpq_clear(ssx->A_val[k]);
      xfree(ssx->A_val);
      xfree(ssx->stat);
      xfree(ssx->Q_row);
      xfree(ssx->Q_col);
      bfx_delete_binv(ssx->binv);
      for (i = 0; i <= m; i++) mpq_clear(ssx->bbar[i]);
      xfree(ssx->bbar);
      for (i = 1; i <= m; i++) mpq_clear(ssx->pi[i]);
      xfree(ssx->pi);
      for (j = 1; j <= n; j++) mpq_clear(ssx->cbar[j]);
      xfree(ssx->cbar);
      for (i = 1; i <= m; i++) mpq_clear(ssx->rho[i]);
      xfree(ssx->rho);
      for (j = 1; j <= n; j++) mpq_clear(ssx->ap[j]);
      xfree(ssx->ap);
      for (i = 1; i <= m; i++) mpq_clear(ssx->aq[i]);
      xfree(ssx->aq);
      mpq_clear(ssx->delta);
      xfree(ssx);
      return;
}

/* eof */
