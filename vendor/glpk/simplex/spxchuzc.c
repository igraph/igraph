/* spxchuzc.c */

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
#include "spxchuzc.h"

/***********************************************************************
*  spx_chuzc_sel - select eligible non-basic variables
*
*  This routine selects eligible non-basic variables xN[j], whose
*  reduced costs d[j] have "wrong" sign, i.e. changing such xN[j] in
*  feasible direction improves (decreases) the objective function.
*
*  Reduced costs of non-basic variables should be placed in the array
*  locations d[1], ..., d[n-m].
*
*  Non-basic variable xN[j] is considered eligible if:
*
*     d[j] <= -eps[j] and xN[j] can increase
*
*     d[j] >= +eps[j] and xN[j] can decrease
*
*  for
*
*     eps[j] = tol + tol1 * |cN[j]|,
*
*  where cN[j] is the objective coefficient at xN[j], tol and tol1 are
*  specified tolerances.
*
*  On exit the routine stores indices j of eligible non-basic variables
*  xN[j] to the array locations list[1], ..., list[num] and returns the
*  number of such variables 0 <= num <= n-m. (If the parameter list is
*  specified as NULL, no indices are stored.) */

int spx_chuzc_sel(SPXLP *lp, const double d[/*1+n-m*/], double tol,
      double tol1, int list[/*1+n-m*/])
{     int m = lp->m;
      int n = lp->n;
      double *c = lp->c;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      int j, k, num;
      double ck, eps;
      num = 0;
      /* walk thru list of non-basic variables */
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         if (l[k] == u[k])
         {  /* xN[j] is fixed variable; skip it */
            continue;
         }
         /* determine absolute tolerance eps[j] */
         ck = c[k];
         eps = tol + tol1 * (ck >= 0.0 ? +ck : -ck);
         /* check if xN[j] is eligible */
         if (d[j] <= -eps)
         {  /* xN[j] should be able to increase */
            if (flag[j])
            {  /* but its upper bound is active */
               continue;
            }
         }
         else if (d[j] >= +eps)
         {  /* xN[j] should be able to decrease */
            if (!flag[j] && l[k] != -DBL_MAX)
            {  /* but its lower bound is active */
               continue;
            }
         }
         else /* -eps < d[j] < +eps */
         {  /* xN[j] does not affect the objective function within the
             * specified tolerance */
            continue;
         }
         /* xN[j] is eligible non-basic variable */
         num++;
         if (list != NULL)
            list[num] = j;
      }
      return num;
}

/***********************************************************************
*  spx_chuzc_std - choose non-basic variable (Dantzig's rule)
*
*  This routine chooses most eligible non-basic variable xN[q]
*  according to Dantzig's ("standard") rule:
*
*     d[q] =   max |d[j]|,
*            j in J
*
*  where J <= {1, ..., n-m} is the set of indices of eligible non-basic
*  variables, d[j] is the reduced cost of non-basic variable xN[j] in
*  the current basis.
*
*  Reduced costs of non-basic variables should be placed in the array
*  locations d[1], ..., d[n-m].
*
*  Indices of eligible non-basic variables j in J should be placed in
*  the array locations list[1], ..., list[num], where num = |J| > 0 is
*  the total number of such variables.
*
*  On exit the routine returns q, the index of the non-basic variable
*  xN[q] chosen. */

int spx_chuzc_std(SPXLP *lp, const double d[/*1+n-m*/], int num,
      const int list[])
{     int m = lp->m;
      int n = lp->n;
      int j, q, t;
      double abs_dj, abs_dq;
      xassert(0 < num && num <= n-m);
      q = 0, abs_dq = -1.0;
      for (t = 1; t <= num; t++)
      {  j = list[t];
         abs_dj = (d[j] >= 0.0 ? +d[j] : -d[j]);
         if (abs_dq < abs_dj)
            q = j, abs_dq = abs_dj;
      }
      xassert(q != 0);
      return q;
}

/***********************************************************************
*  spx_alloc_se - allocate pricing data block
*
*  This routine allocates the memory for arrays used in the pricing
*  data block. */

void spx_alloc_se(SPXLP *lp, SPXSE *se)
{     int m = lp->m;
      int n = lp->n;
      se->valid = 0;
      se->refsp = talloc(1+n, char);
      se->gamma = talloc(1+n-m, double);
      se->work = talloc(1+m, double);
      return;
}

/***********************************************************************
*  spx_reset_refsp - reset reference space
*
*  This routine resets (re-initializes) the reference space composing
*  it from variables which are non-basic in the current basis, and sets
*  all weights gamma[j] to 1. */

void spx_reset_refsp(SPXLP *lp, SPXSE *se)
{     int m = lp->m;
      int n = lp->n;
      int *head = lp->head;
      char *refsp = se->refsp;
      double *gamma = se->gamma;
      int j, k;
      se->valid = 1;
      memset(&refsp[1], 0, n * sizeof(char));
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         refsp[k] = 1;
         gamma[j] = 1.0;
      }
      return;
}

/***********************************************************************
*  spx_eval_gamma_j - compute projected steepest edge weight directly
*
*  This routine computes projected steepest edge weight gamma[j],
*  1 <= j <= n-m, for the current basis directly with the formula:
*
*                            m
*     gamma[j] = delta[j] + sum eta[i] * T[i,j]**2,
*                           i=1
*
*  where T[i,j] is element of the current simplex table, and
*
*                ( 1, if xB[i] is in the reference space
*     eta[i]   = {
*                ( 0, otherwise
*
*                ( 1, if xN[j] is in the reference space
*     delta[j] = {
*                ( 0, otherwise
*
*  NOTE: For testing/debugging only. */

double spx_eval_gamma_j(SPXLP *lp, SPXSE *se, int j)
{     int m = lp->m;
      int n = lp->n;
      int *head = lp->head;
      char *refsp = se->refsp;
      double *tcol = se->work;
      int i, k;
      double gamma_j;
      xassert(se->valid);
      xassert(1 <= j && j <= n-m);
      k = head[m+j]; /* x[k] = xN[j] */
      gamma_j = (refsp[k] ? 1.0 : 0.0);
      spx_eval_tcol(lp, j, tcol);
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
         if (refsp[k])
            gamma_j += tcol[i] * tcol[i];
      }
      return gamma_j;
}

/***********************************************************************
*  spx_chuzc_pse - choose non-basic variable (projected steepest edge)
*
*  This routine chooses most eligible non-basic variable xN[q]
*  according to the projected steepest edge method:
*
*      d[q]**2           d[j]**2
*     -------- =   max  -------- ,
*     gamma[q]   j in J gamma[j]
*
*  where J <= {1, ..., n-m} is the set of indices of eligible non-basic
*  variable, d[j] is the reduced cost of non-basic variable xN[j] in
*  the current basis, gamma[j] is the projected steepest edge weight.
*
*  Reduced costs of non-basic variables should be placed in the array
*  locations d[1], ..., d[n-m].
*
*  Indices of eligible non-basic variables j in J should be placed in
*  the array locations list[1], ..., list[num], where num = |J| > 0 is
*  the total number of such variables.
*
*  On exit the routine returns q, the index of the non-basic variable
*  xN[q] chosen. */

int spx_chuzc_pse(SPXLP *lp, SPXSE *se, const double d[/*1+n-m*/],
      int num, const int list[])
{     int m = lp->m;
      int n = lp->n;
      double *gamma = se->gamma;
      int j, q, t;
      double best, temp;
      xassert(se->valid);
      xassert(0 < num && num <= n-m);
      q = 0, best = -1.0;
      for (t = 1; t <= num; t++)
      {  j = list[t];
         /* FIXME */
         if (gamma[j] < DBL_EPSILON)
            temp = 0.0;
         else
            temp = (d[j] * d[j]) / gamma[j];
         if (best < temp)
            q = j, best = temp;
      }
      xassert(q != 0);
      return q;
}

/***********************************************************************
*  spx_update_gamma - update projected steepest edge weights exactly
*
*  This routine updates the vector gamma = (gamma[j]) of projected
*  steepest edge weights exactly, for the adjacent basis.
*
*  On entry to the routine the content of the se object should be valid
*  and should correspond to the current basis.
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
*  For details about the formulae used see the program documentation.
*
*  The routine also computes the relative error:
*
*     e = |gamma[q] - gamma'[q]| / (1 + |gamma[q]|),
*
*  where gamma'[q] is the weight for xN[q] on entry to the routine,
*  and returns e on exit. (If e happens to be large enough, the calling
*  program may reset the reference space, since other weights also may
*  be inaccurate.) */

double spx_update_gamma(SPXLP *lp, SPXSE *se, int p, int q,
      const double trow[/*1+n-m*/], const double tcol[/*1+m*/])
{     int m = lp->m;
      int n = lp->n;
      int *head = lp->head;
      char *refsp = se->refsp;
      double *gamma = se->gamma;
      double *u = se->work;
      int i, j, k, ptr, end;
      double gamma_q, delta_q, e, r, s, t1, t2;
      xassert(se->valid);
      xassert(1 <= p && p <= m);
      xassert(1 <= q && q <= n-m);
      /* compute gamma[q] in current basis more accurately; also
       * compute auxiliary vector u */
      k = head[m+q]; /* x[k] = xN[q] */
      gamma_q = delta_q = (refsp[k] ? 1.0 : 0.0);
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
         if (refsp[k])
         {  gamma_q += tcol[i] * tcol[i];
            u[i] = tcol[i];
         }
         else
            u[i] = 0.0;
      }
      bfd_btran(lp->bfd, u);
      /* compute relative error in gamma[q] */
      e = fabs(gamma_q - gamma[q]) / (1.0 + gamma_q);
      /* compute new gamma[q] */
      gamma[q] = gamma_q / (tcol[p] * tcol[p]);
      /* compute new gamma[j] for all j != q */
      for (j = 1; j <= n-m; j++)
      {  if (j == q)
            continue;
         if (-1e-9 < trow[j] && trow[j] < +1e-9)
         {  /* T[p,j] is close to zero; gamma[j] is not changed */
            continue;
         }
         /* compute r[j] = T[p,j] / T[p,q] */
         r = trow[j] / tcol[p];
         /* compute inner product s[j] = N'[j] * u, where N[j] = A[k]
          * is constraint matrix column corresponding to xN[j] */
         s = 0.0;
         k = head[m+j]; /* x[k] = xN[j] */
         ptr = lp->A_ptr[k];
         end = lp->A_ptr[k+1];
         for (; ptr < end; ptr++)
            s += lp->A_val[ptr] * u[lp->A_ind[ptr]];
         /* compute new gamma[j] */
         t1 = gamma[j] + r * (r * gamma_q + s + s);
         t2 = (refsp[k] ? 1.0 : 0.0) + delta_q * r * r;
         gamma[j] = (t1 >= t2 ? t1 : t2);
      }
      return e;
}

/***********************************************************************
*  spx_free_se - deallocate pricing data block
*
*  This routine deallocates the memory used for arrays in the pricing
*  data block. */

void spx_free_se(SPXLP *lp, SPXSE *se)
{     xassert(lp == lp);
      tfree(se->refsp);
      tfree(se->gamma);
      tfree(se->work);
      return;
}

/* eof */
