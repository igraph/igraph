/* glpapi10.c (solution checking routines) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000-2013 Free Software Foundation, Inc.
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
#include "prob.h"

void glp_check_kkt(glp_prob *P, int sol, int cond, double *_ae_max,
      int *_ae_ind, double *_re_max, int *_re_ind)
{     /* check feasibility and optimality conditions */
      int m = P->m;
      int n = P->n;
      GLPROW *row;
      GLPCOL *col;
      GLPAIJ *aij;
      int i, j, ae_ind, re_ind;
      double e, sp, sn, t, ae_max, re_max;
      if (!(sol == GLP_SOL || sol == GLP_IPT || sol == GLP_MIP))
         xerror("glp_check_kkt: sol = %d; invalid solution indicator\n",
            sol);
      if (!(cond == GLP_KKT_PE || cond == GLP_KKT_PB ||
            cond == GLP_KKT_DE || cond == GLP_KKT_DB ||
            cond == GLP_KKT_CS))
         xerror("glp_check_kkt: cond = %d; invalid condition indicator "
            "\n", cond);
      ae_max = re_max = 0.0;
      ae_ind = re_ind = 0;
      if (cond == GLP_KKT_PE)
      {  /* xR - A * xS = 0 */
         for (i = 1; i <= m; i++)
         {  row = P->row[i];
            sp = sn = 0.0;
            /* t := xR[i] */
            if (sol == GLP_SOL)
               t = row->prim;
            else if (sol == GLP_IPT)
               t = row->pval;
            else if (sol == GLP_MIP)
               t = row->mipx;
            else
               xassert(sol != sol);
            if (t >= 0.0) sp += t; else sn -= t;
            for (aij = row->ptr; aij != NULL; aij = aij->r_next)
            {  col = aij->col;
               /* t := - a[i,j] * xS[j] */
               if (sol == GLP_SOL)
                  t = - aij->val * col->prim;
               else if (sol == GLP_IPT)
                  t = - aij->val * col->pval;
               else if (sol == GLP_MIP)
                  t = - aij->val * col->mipx;
               else
                  xassert(sol != sol);
               if (t >= 0.0) sp += t; else sn -= t;
            }
            /* absolute error */
            e = fabs(sp - sn);
            if (ae_max < e)
               ae_max = e, ae_ind = i;
            /* relative error */
            e /= (1.0 + sp + sn);
            if (re_max < e)
               re_max = e, re_ind = i;
         }
      }
      else if (cond == GLP_KKT_PB)
      {  /* lR <= xR <= uR */
         for (i = 1; i <= m; i++)
         {  row = P->row[i];
            /* t := xR[i] */
            if (sol == GLP_SOL)
               t = row->prim;
            else if (sol == GLP_IPT)
               t = row->pval;
            else if (sol == GLP_MIP)
               t = row->mipx;
            else
               xassert(sol != sol);
            /* check lower bound */
            if (row->type == GLP_LO || row->type == GLP_DB ||
                row->type == GLP_FX)
            {  if (t < row->lb)
               {  /* absolute error */
                  e = row->lb - t;
                  if (ae_max < e)
                     ae_max = e, ae_ind = i;
                  /* relative error */
                  e /= (1.0 + fabs(row->lb));
                  if (re_max < e)
                     re_max = e, re_ind = i;
               }
            }
            /* check upper bound */
            if (row->type == GLP_UP || row->type == GLP_DB ||
                row->type == GLP_FX)
            {  if (t > row->ub)
               {  /* absolute error */
                  e = t - row->ub;
                  if (ae_max < e)
                     ae_max = e, ae_ind = i;
                  /* relative error */
                  e /= (1.0 + fabs(row->ub));
                  if (re_max < e)
                     re_max = e, re_ind = i;
               }
            }
         }
         /* lS <= xS <= uS */
         for (j = 1; j <= n; j++)
         {  col = P->col[j];
            /* t := xS[j] */
            if (sol == GLP_SOL)
               t = col->prim;
            else if (sol == GLP_IPT)
               t = col->pval;
            else if (sol == GLP_MIP)
               t = col->mipx;
            else
               xassert(sol != sol);
            /* check lower bound */
            if (col->type == GLP_LO || col->type == GLP_DB ||
                col->type == GLP_FX)
            {  if (t < col->lb)
               {  /* absolute error */
                  e = col->lb - t;
                  if (ae_max < e)
                     ae_max = e, ae_ind = m+j;
                  /* relative error */
                  e /= (1.0 + fabs(col->lb));
                  if (re_max < e)
                     re_max = e, re_ind = m+j;
               }
            }
            /* check upper bound */
            if (col->type == GLP_UP || col->type == GLP_DB ||
                col->type == GLP_FX)
            {  if (t > col->ub)
               {  /* absolute error */
                  e = t - col->ub;
                  if (ae_max < e)
                     ae_max = e, ae_ind = m+j;
                  /* relative error */
                  e /= (1.0 + fabs(col->ub));
                  if (re_max < e)
                     re_max = e, re_ind = m+j;
               }
            }
         }
      }
      else if (cond == GLP_KKT_DE)
      {  /* A' * (lambdaR - cR) + (lambdaS - cS) = 0 */
         for (j = 1; j <= n; j++)
         {  col = P->col[j];
            sp = sn = 0.0;
            /* t := lambdaS[j] - cS[j] */
            if (sol == GLP_SOL)
               t = col->dual - col->coef;
            else if (sol == GLP_IPT)
               t = col->dval - col->coef;
            else
               xassert(sol != sol);
            if (t >= 0.0) sp += t; else sn -= t;
            for (aij = col->ptr; aij != NULL; aij = aij->c_next)
            {  row = aij->row;
               /* t := a[i,j] * (lambdaR[i] - cR[i]) */
               if (sol == GLP_SOL)
                  t = aij->val * row->dual;
               else if (sol == GLP_IPT)
                  t = aij->val * row->dval;
               else
                  xassert(sol != sol);
               if (t >= 0.0) sp += t; else sn -= t;
            }
            /* absolute error */
            e = fabs(sp - sn);
            if (ae_max < e)
               ae_max = e, ae_ind = m+j;
            /* relative error */
            e /= (1.0 + sp + sn);
            if (re_max < e)
               re_max = e, re_ind = m+j;
         }
      }
      else if (cond == GLP_KKT_DB)
      {  /* check lambdaR */
         for (i = 1; i <= m; i++)
         {  row = P->row[i];
            /* t := lambdaR[i] */
            if (sol == GLP_SOL)
               t = row->dual;
            else if (sol == GLP_IPT)
               t = row->dval;
            else
               xassert(sol != sol);
            /* correct sign */
            if (P->dir == GLP_MIN)
               t = + t;
            else if (P->dir == GLP_MAX)
               t = - t;
            else
               xassert(P != P);
            /* check for positivity */
#if 1 /* 08/III-2013 */
            /* the former check was correct */
            /* the bug reported by David Price is related to violation
               of complementarity slackness, not to this condition */
            if (row->type == GLP_FR || row->type == GLP_LO)
#else
            if (row->stat == GLP_NF || row->stat == GLP_NL)
#endif
            {  if (t < 0.0)
               {  e = - t;
                  if (ae_max < e)
                     ae_max = re_max = e, ae_ind = re_ind = i;
               }
            }
            /* check for negativity */
#if 1 /* 08/III-2013 */
            /* see comment above */
            if (row->type == GLP_FR || row->type == GLP_UP)
#else
            if (row->stat == GLP_NF || row->stat == GLP_NU)
#endif
            {  if (t > 0.0)
               {  e = + t;
                  if (ae_max < e)
                     ae_max = re_max = e, ae_ind = re_ind = i;
               }
            }
         }
         /* check lambdaS */
         for (j = 1; j <= n; j++)
         {  col = P->col[j];
            /* t := lambdaS[j] */
            if (sol == GLP_SOL)
               t = col->dual;
            else if (sol == GLP_IPT)
               t = col->dval;
            else
               xassert(sol != sol);
            /* correct sign */
            if (P->dir == GLP_MIN)
               t = + t;
            else if (P->dir == GLP_MAX)
               t = - t;
            else
               xassert(P != P);
            /* check for positivity */
#if 1 /* 08/III-2013 */
            /* see comment above */
            if (col->type == GLP_FR || col->type == GLP_LO)
#else
            if (col->stat == GLP_NF || col->stat == GLP_NL)
#endif
            {  if (t < 0.0)
               {  e = - t;
                  if (ae_max < e)
                     ae_max = re_max = e, ae_ind = re_ind = m+j;
               }
            }
            /* check for negativity */
#if 1 /* 08/III-2013 */
            /* see comment above */
            if (col->type == GLP_FR || col->type == GLP_UP)
#else
            if (col->stat == GLP_NF || col->stat == GLP_NU)
#endif
            {  if (t > 0.0)
               {  e = + t;
                  if (ae_max < e)
                     ae_max = re_max = e, ae_ind = re_ind = m+j;
               }
            }
         }
      }
      else
         xassert(cond != cond);
      if (_ae_max != NULL) *_ae_max = ae_max;
      if (_ae_ind != NULL) *_ae_ind = ae_ind;
      if (_re_max != NULL) *_re_max = re_max;
      if (_re_ind != NULL) *_re_ind = re_ind;
      return;
}

/* eof */
