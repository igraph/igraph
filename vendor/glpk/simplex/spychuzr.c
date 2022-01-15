/* spychuzr.c */

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
#include "spychuzr.h"

/***********************************************************************
*  spy_chuzr_sel - select eligible basic variables
*
*  This routine selects eligible basic variables xB[i], whose value
*  beta[i] violates corresponding lower lB[i] or upper uB[i] bound.
*  Positive bound violation rp[i] = lb[i] - beta[i] > 0 is the reduced
*  cost of non-basic dual variable lambda^+B[i] >= 0, so increasing it
*  increases the dual objective. Similarly, negative bound violation
*  rn[i] = ub[i] - beta[i] < 0 is the reduced cost of non-basic dual
*  variable lambda^-B[i] <= 0, so decreasing it also increases the dual
*  objective.
*
*  Current values of basic variables should be placed in the array
*  locations beta[1], ..., beta[m].
*
*  Basic variable xB[i] is considered eligible, if:
*
*     beta[i] <= lB[i] - eps1[i], or
*
*     beta[i] >= uB[i] + eps2[i],
*
*  for
*
*     eps1[i] = tol + tol1 * |lB[i]|,
*
*     eps2[i] = tol + tol2 * |uB[i]|,
*
*  where lB[i] and uB[i] are, resp., lower and upper bounds of xB[i],
*  tol and tol1 are specified tolerances.
*
*  On exit the routine stores indices i of eligible basic variables
*  xB[i] to the array locations list[1], ..., list[num] and returns the
*  number of such variables 0 <= num <= m. (If the parameter list is
*  specified as NULL, no indices are stored.) */

int spy_chuzr_sel(SPXLP *lp, const double beta[/*1+m*/], double tol,
      double tol1, int list[/*1+m*/])
{     int m = lp->m;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      int i, k, num;
      double lk, uk, eps;
      num = 0;
      /* walk thru list of basic variables */
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
         lk = l[k], uk = u[k];
         /* check if xB[i] is eligible */
         if (beta[i] < lk)
         {  /* determine absolute tolerance eps1[i] */
            eps = tol + tol1 * (lk >= 0.0 ? +lk : -lk);
            if (beta[i] < lk - eps)
            {  /* lower bound is violated */
               num++;
               if (list != NULL)
                  list[num] = i;
            }
         }
         else if (beta[i] > uk)
         {  /* determine absolute tolerance eps2[i] */
            eps = tol + tol1 * (uk >= 0.0 ? +uk : -uk);
            if (beta[i] > uk + eps)
            {  /* upper bound is violated */
               num++;
               if (list != NULL)
                  list[num] = i;
            }
         }
      }
      return num;
}

/***********************************************************************
*  spy_chuzr_std - choose basic variable (dual Dantzig's rule)
*
*  This routine chooses most eligible basic variable xB[p] according
*  to dual Dantzig's ("standard") rule:
*
*     r[p] =   max  |r[i]|,
*            i in I
*
*            ( lB[i] - beta[i], if beta[i] < lB[i]
*            (
*     r[i] = { 0,               if lB[i] <= beta[i] <= uB[i]
*            (
*            ( uB[i] - beta[i], if beta[i] > uB[i]
*
*  where I <= {1, ..., m} is the set of indices of eligible basic
*  variables, beta[i] is current value of xB[i], lB[i] and uB[i] are,
*  resp., lower and upper bounds of xB[i], r[i] is bound violation.
*
*  Current values of basic variables should be placed in the array
*  locations beta[1], ..., beta[m].
*
*  Indices of eligible basic variables i in I should be placed in the
*  array locations list[1], ..., list[num], where num = |J| > 0 is the
*  total number of such variables.
*
*  On exit the routine returns p, the index of the basic variable xB[p]
*  chosen. */

int spy_chuzr_std(SPXLP *lp, const double beta[/*1+m*/], int num,
      const int list[])
{     int m = lp->m;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      int i, k, p, t;
      double abs_ri, abs_rp;
      xassert(0 < num && num <= m);
      p = 0, abs_rp = -1.0;
      for (t = 1; t <= num; t++)
      {  i = list[t];
         k = head[i]; /* x[k] = xB[i] */
         if (beta[i] < l[k])
            abs_ri = l[k] - beta[i];
         else if (beta[i] > u[k])
            abs_ri = beta[i] - u[k];
         else
            xassert(t != t);
         if (abs_rp < abs_ri)
            p = i, abs_rp = abs_ri;
      }
      xassert(p != 0);
      return p;
}

/***********************************************************************
*  spy_alloc_se - allocate dual pricing data block
*
*  This routine allocates the memory for arrays used in the dual
*  pricing data block. */

void spy_alloc_se(SPXLP *lp, SPYSE *se)
{     int m = lp->m;
      int n = lp->n;
#if 1 /* 30/III-2016 */
      int i;
#endif
      se->valid = 0;
      se->refsp = talloc(1+n, char);
      se->gamma = talloc(1+m, double);
      se->work = talloc(1+m, double);
#if 1 /* 30/III-2016 */
      se->u.n = m;
      se->u.nnz = 0;
      se->u.ind = talloc(1+m, int);
      se->u.vec = talloc(1+m, double);
      for (i = 1; i <= m; i++)
         se->u.vec[i] = 0.0;
#endif
      return;
}

/***********************************************************************
*  spy_reset_refsp - reset dual reference space
*
*  This routine resets (re-initializes) the dual reference space
*  composing it from dual variables which are non-basic (corresponding
*  to basic primal variables) in the current basis, and sets all
*  weights gamma[i] to 1. */

void spy_reset_refsp(SPXLP *lp, SPYSE *se)
{     int m = lp->m;
      int n = lp->n;
      int *head = lp->head;
      char *refsp = se->refsp;
      double *gamma = se->gamma;
      int i, k;
      se->valid = 1;
      memset(&refsp[1], 0, n * sizeof(char));
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
         refsp[k] = 1;
         gamma[i] = 1.0;
      }
      return;
}

/***********************************************************************
*  spy_eval_gamma_i - compute dual proj. steepest edge weight directly
*
*  This routine computes dual projected steepest edge weight gamma[i],
*  1 <= i <= m, for the current basis directly with the formula:
*
*                           n-m
*     gamma[i] = delta[i] + sum eta[j] * T[i,j]**2,
*                           j=1
*
*  where T[i,j] is element of the current simplex table, and
*
*                ( 1, if lambdaN[j] is in the reference space
*     eta[j]   = {
*                ( 0, otherwise
*
*                ( 1, if lambdaB[i] is in the reference space
*     delta[i] = {
*                ( 0, otherwise
*
*  Dual basic variable lambdaN[j] corresponds to primal non-basic
*  variable xN[j], and dual non-basic variable lambdaB[j] corresponds
*  to primal basic variable xB[i].
*
*  NOTE: For testing/debugging only. */

double spy_eval_gamma_i(SPXLP *lp, SPYSE *se, int i)
{     int m = lp->m;
      int n = lp->n;
      int *head = lp->head;
      char *refsp = se->refsp;
      double *rho = se->work;
      int j, k;
      double gamma_i, t_ij;
      xassert(se->valid);
      xassert(1 <= i && i <= m);
      k = head[i]; /* x[k] = xB[i] */
      gamma_i = (refsp[k] ? 1.0 : 0.0);
      spx_eval_rho(lp, i, rho);
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         if (refsp[k])
         {  t_ij = spx_eval_tij(lp, rho, j);
            gamma_i += t_ij * t_ij;
         }
      }
      return gamma_i;
}

/***********************************************************************
*  spy_chuzr_pse - choose basic variable (dual projected steepest edge)
*
*  This routine chooses most eligible basic variable xB[p] according
*  to the dual projected steepest edge method:
*
*      r[p]**2           r[i]**2
*     -------- =   max  -------- ,
*     gamma[p]   i in I gamma[i]
*
*            ( lB[i] - beta[i], if beta[i] < lB[i]
*            (
*     r[i] = { 0,               if lB[i] <= beta[i] <= uB[i]
*            (
*            ( uB[i] - beta[i], if beta[i] > uB[i]
*
*  where I <= {1, ..., m} is the set of indices of eligible basic
*  variables, beta[i] is current value of xB[i], lB[i] and uB[i] are,
*  resp., lower and upper bounds of xB[i], r[i] is bound violation.
*
*  Current values of basic variables should be placed in the array
*  locations beta[1], ..., beta[m].
*
*  Indices of eligible basic variables i in I should be placed in the
*  array locations list[1], ..., list[num], where num = |J| > 0 is the
*  total number of such variables.
*
*  On exit the routine returns p, the index of the basic variable xB[p]
*  chosen. */

int spy_chuzr_pse(SPXLP *lp, SPYSE *se, const double beta[/*1+m*/],
      int num, const int list[])
{     int m = lp->m;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      double *gamma = se->gamma;
      int i, k, p, t;
      double best, ri, temp;
      xassert(0 < num && num <= m);
      p = 0, best = -1.0;
      for (t = 1; t <= num; t++)
      {  i = list[t];
         k = head[i]; /* x[k] = xB[i] */
         if (beta[i] < l[k])
            ri = l[k] - beta[i];
         else if (beta[i] > u[k])
            ri = u[k] - beta[i];
         else
            xassert(t != t);
         /* FIXME */
         if (gamma[i] < DBL_EPSILON)
            temp = 0.0;
         else
            temp = (ri * ri) / gamma[i];
         if (best < temp)
            p = i, best = temp;
      }
      xassert(p != 0);
      return p;
}

/***********************************************************************
*  spy_update_gamma - update dual proj. steepest edge weights exactly
*
*  This routine updates the vector gamma = (gamma[i]) of dual projected
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
*     e = |gamma[p] - gamma'[p]| / (1 + |gamma[p]|),
*
*  where gamma'[p] is the weight for lambdaB[p] (which is dual
*  non-basic variable corresponding to xB[p]) on entry to the routine,
*  and returns e on exit. (If e happens to be large enough, the calling
*  program may reset the reference space, since other weights also may
*  be inaccurate.) */

double spy_update_gamma(SPXLP *lp, SPYSE *se, int p, int q,
      const double trow[/*1+n-m*/], const double tcol[/*1+m*/])
{     int m = lp->m;
      int n = lp->n;
      int *head = lp->head;
      char *refsp = se->refsp;
      double *gamma = se->gamma;
      double *u = se->work;
      int i, j, k, ptr, end;
      double gamma_p, delta_p, e, r, t1, t2;
      xassert(se->valid);
      xassert(1 <= p && p <= m);
      xassert(1 <= q && q <= n-m);
      /* compute gamma[p] in current basis more accurately; also
       * compute auxiliary vector u */
      k = head[p]; /* x[k] = xB[p] */
      gamma_p = delta_p = (refsp[k] ? 1.0 : 0.0);
      for (i = 1; i <= m; i++)
         u[i] = 0.0;
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         if (refsp[k] && trow[j] != 0.0)
         {  gamma_p += trow[j] * trow[j];
            /* u := u + T[p,j] * N[j], where N[j] = A[k] is constraint
             * matrix column corresponding to xN[j] */
            ptr = lp->A_ptr[k];
            end = lp->A_ptr[k+1];
            for (; ptr < end; ptr++)
               u[lp->A_ind[ptr]] += trow[j] * lp->A_val[ptr];
         }
      }
      bfd_ftran(lp->bfd, u);
      /* compute relative error in gamma[p] */
      e = fabs(gamma_p - gamma[p]) / (1.0 + gamma_p);
      /* compute new gamma[p] */
      gamma[p] = gamma_p / (tcol[p] * tcol[p]);
      /* compute new gamma[i] for all i != p */
      for (i = 1; i <= m; i++)
      {  if (i == p)
            continue;
         /* compute r[i] = T[i,q] / T[p,q] */
         r = tcol[i] / tcol[p];
         /* compute new gamma[i] */
         t1 = gamma[i] + r * (r * gamma_p + u[i] + u[i]);
         k = head[i]; /* x[k] = xB[i] */
         t2 = (refsp[k] ? 1.0 : 0.0) + delta_p * r * r;
         gamma[i] = (t1 >= t2 ? t1 : t2);
      }
      return e;
}

#if 1 /* 30/III-2016 */
double spy_update_gamma_s(SPXLP *lp, SPYSE *se, int p, int q,
      const FVS *trow, const FVS *tcol)
{     /* sparse version of spy_update_gamma */
      int m = lp->m;
      int n = lp->n;
      int *head = lp->head;
      char *refsp = se->refsp;
      double *gamma = se->gamma;
      double *u = se->work;
      int trow_nnz = trow->nnz;
      int *trow_ind = trow->ind;
      double *trow_vec = trow->vec;
      int tcol_nnz = tcol->nnz;
      int *tcol_ind = tcol->ind;
      double *tcol_vec = tcol->vec;
      int i, j, k, t, ptr, end;
      double gamma_p, delta_p, e, r, t1, t2;
      xassert(se->valid);
      xassert(1 <= p && p <= m);
      xassert(1 <= q && q <= n-m);
      /* compute gamma[p] in current basis more accurately; also
       * compute auxiliary vector u */
      k = head[p]; /* x[k] = xB[p] */
      gamma_p = delta_p = (refsp[k] ? 1.0 : 0.0);
      for (i = 1; i <= m; i++)
         u[i] = 0.0;
      for (t = 1; t <= trow_nnz; t++)
      {  j = trow_ind[t];
         k = head[m+j]; /* x[k] = xN[j] */
         if (refsp[k])
         {  gamma_p += trow_vec[j] * trow_vec[j];
            /* u := u + T[p,j] * N[j], where N[j] = A[k] is constraint
             * matrix column corresponding to xN[j] */
            ptr = lp->A_ptr[k];
            end = lp->A_ptr[k+1];
            for (; ptr < end; ptr++)
               u[lp->A_ind[ptr]] += trow_vec[j] * lp->A_val[ptr];
         }
      }
      bfd_ftran(lp->bfd, u);
      /* compute relative error in gamma[p] */
      e = fabs(gamma_p - gamma[p]) / (1.0 + gamma_p);
      /* compute new gamma[p] */
      gamma[p] = gamma_p / (tcol_vec[p] * tcol_vec[p]);
      /* compute new gamma[i] for all i != p */
      for (t = 1; t <= tcol_nnz; t++)
      {  i = tcol_ind[t];
         if (i == p)
            continue;
         /* compute r[i] = T[i,q] / T[p,q] */
         r = tcol_vec[i] / tcol_vec[p];
         /* compute new gamma[i] */
         t1 = gamma[i] + r * (r * gamma_p + u[i] + u[i]);
         k = head[i]; /* x[k] = xB[i] */
         t2 = (refsp[k] ? 1.0 : 0.0) + delta_p * r * r;
         gamma[i] = (t1 >= t2 ? t1 : t2);
      }
      return e;
}
#endif

/***********************************************************************
*  spy_free_se - deallocate dual pricing data block
*
*  This routine deallocates the memory used for arrays in the dual
*  pricing data block. */

void spy_free_se(SPXLP *lp, SPYSE *se)
{     xassert(lp == lp);
      tfree(se->refsp);
      tfree(se->gamma);
      tfree(se->work);
#if 1 /* 30/III-2016 */
      tfree(se->u.ind);
      tfree(se->u.vec);
#endif
      return;
}

/* eof */
