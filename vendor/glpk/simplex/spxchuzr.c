/* spxchuzr.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2015-2018 Free Software Foundation, Inc.
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
#include "spxchuzr.h"

/***********************************************************************
*  spx_chuzr_std - choose basic variable (textbook ratio test)
*
*  This routine implements an improved textbook ratio test to choose
*  basic variable xB[p].
*
*  The parameter phase specifies the search phase:
*
*  1 - searching for feasible basic solution. In this case the routine
*      uses artificial bounds of basic variables that correspond to
*      breakpoints of the penalty function:
*
*                 ( lB[i], if cB[i] = 0
*                 (
*        lB'[i] = { uB[i], if cB[i] > 0
*                 (
*                 (  -inf, if cB[i] < 0
*
*                 ( uB[i], if cB[i] = 0
*                 (
*        uB'[i] = {  +inf, if cB[i] > 0
*                 (
*                 ( lB[i], if cB[i] < 0
*
*      where lB[i] and uB[i] are original bounds of variable xB[i],
*      cB[i] is the penalty (objective) coefficient of that variable.
*
*  2 - searching for optimal basic solution. In this case the routine
*      uses original bounds of basic variables.
*
*  Current values of basic variables should be placed in the array
*  locations beta[1], ..., beta[m].
*
*  The parameter 1 <= q <= n-m specifies the index of non-basic
*  variable xN[q] chosen.
*
*  The parameter s specifies the direction in which xN[q] changes:
*  s = +1.0 means xN[q] increases, and s = -1.0 means xN[q] decreases.
*  (Thus, the corresponding ray parameter is theta = s (xN[q] - f[q]),
*  where f[q] is the active bound of xN[q] in the current basis.)
*
*  Elements of q-th simplex table column T[q] = (t[i,q]) corresponding
*  to non-basic variable xN[q] should be placed in the array locations
*  tcol[1], ..., tcol[m].
*
*  The parameter tol_piv specifies a tolerance for elements of the
*  simplex table column T[q]. If |t[i,q]| < tol_piv, basic variable
*  xB[i] is skipped, i.e. it is assumed that it does not depend on the
*  ray parameter theta.
*
*  The parameters tol and tol1 specify tolerances used to increase the
*  choice freedom by simulating an artificial degeneracy as follows.
*  If beta[i] <= lB[i] + delta[i], where delta[i] = tol + tol1 |lB[i]|,
*  it is assumed that beta[i] is exactly the same as lB[i]. Similarly,
*  if beta[i] >= uB[i] - delta[i], where delta[i] = tol + tol1 |uB[i]|,
*  it is assumed that beta[i] is exactly the same as uB[i].
*
*  The routine determines the index 1 <= p <= m of basic variable xB[p]
*  that reaches its (lower or upper) bound first on increasing the ray
*  parameter theta, stores the bound flag (0 - lower bound or fixed
*  value, 1 - upper bound) to the location pointed to by the pointer
*  p_flag, and returns the index p. If non-basic variable xN[q] is
*  double-bounded and reaches its opposite bound first, the routine
*  returns (-1). And if the ray parameter may increase unlimitedly, the
*  routine returns zero.
*
*  Should note that the bound flag stored to the location pointed to by
*  p_flag corresponds to the original (not artficial) bound of variable
*  xB[p] and defines the active bound flag lp->flag[q] to be set in the
*  adjacent basis for that basic variable. */

int spx_chuzr_std(SPXLP *lp, int phase, const double beta[/*1+m*/],
      int q, double s, const double tcol[/*1+m*/], int *p_flag,
      double tol_piv, double tol, double tol1)
{     int m = lp->m;
      int n = lp->n;
      double *c = lp->c;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      int i, i_flag, k, p;
      double alfa, biga, delta, lk, uk, teta, teta_min;
      xassert(phase == 1 || phase == 2);
      xassert(1 <= q && q <= n-m);
      xassert(s == +1.0 || s == -1.0);
      /* determine initial teta_min */
      k = head[m+q]; /* x[k] = xN[q] */
      if (l[k] == -DBL_MAX || u[k] == +DBL_MAX)
      {  /* xN[q] has no opposite bound */
         p = 0, *p_flag = 0, teta_min = DBL_MAX, biga = 0.0;
      }
      else
      {  /* xN[q] have both lower and upper bounds */
         p = -1, *p_flag = 0, teta_min = fabs(l[k] - u[k]), biga = 1.0;
      }
      /* walk thru the list of basic variables */
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
         /* determine alfa such that delta xB[i] = alfa * teta */
         alfa = s * tcol[i];
         if (alfa <= -tol_piv)
         {  /* xB[i] decreases */
            /* determine actual lower bound of xB[i] */
            if (phase == 1 && c[k] < 0.0)
            {  /* xB[i] has no actual lower bound */
               continue;
            }
            else if (phase == 1 && c[k] > 0.0)
            {  /* actual lower bound of xB[i] is its upper bound */
               lk = u[k];
               xassert(lk != +DBL_MAX);
               i_flag = 1;
            }
            else
            {  /* actual lower bound of xB[i] is its original bound */
               lk = l[k];
               if (lk == -DBL_MAX)
                  continue;
               i_flag = 0;
            }
            /* determine teta on which xB[i] reaches its lower bound */
            delta = tol + tol1 * (lk >= 0.0 ? +lk : -lk);
            if (beta[i] <= lk + delta)
               teta = 0.0;
            else
               teta = (lk - beta[i]) / alfa;
         }
         else if (alfa >= +tol_piv)
         {  /* xB[i] increases */
            /* determine actual upper bound of xB[i] */
            if (phase == 1 && c[k] < 0.0)
            {  /* actual upper bound of xB[i] is its lower bound */
               uk = l[k];
               xassert(uk != -DBL_MAX);
               i_flag = 0;
            }
            else if (phase == 1 && c[k] > 0.0)
            {  /* xB[i] has no actual upper bound */
               continue;
            }
            else
            {  /* actual upper bound of xB[i] is its original bound */
               uk = u[k];
               if (uk == +DBL_MAX)
                  continue;
               i_flag = 1;
            }
            /* determine teta on which xB[i] reaches its upper bound */
            delta = tol + tol1 * (uk >= 0.0 ? +uk : -uk);
            if (beta[i] >= uk - delta)
               teta = 0.0;
            else
               teta = (uk - beta[i]) / alfa;
         }
         else
         {  /* xB[i] does not depend on teta */
            continue;
         }
         /* choose basic variable xB[p] for which teta is minimal */
         xassert(teta >= 0.0);
         alfa = (alfa >= 0.0 ? +alfa : -alfa);
         if (teta_min > teta || (teta_min == teta && biga < alfa))
            p = i, *p_flag = i_flag, teta_min = teta, biga = alfa;
      }
      /* if xB[p] is fixed variable, adjust its bound flag */
      if (p > 0)
      {  k = head[p];
         if (l[k] == u[k])
            *p_flag = 0;
      }
      return p;
}

/***********************************************************************
*  spx_chuzr_harris - choose basic variable (Harris' ratio test)
*
*  This routine implements Harris' ratio test to choose basic variable
*  xB[p].
*
*  All the parameters, except tol and tol1, as well as the returned
*  value have the same meaning as for the routine spx_chuzr_std (see
*  above).
*
*  The parameters tol and tol1 specify tolerances on bound violations
*  for basic variables. For the lower bound of basic variable xB[i] the
*  tolerance is delta[i] = tol + tol1 |lB[i]|, and for the upper bound
*  the tolerance is delta[i] = tol + tol1 |uB[i]|. */

int spx_chuzr_harris(SPXLP *lp, int phase, const double beta[/*1+m*/],
      int q, double s, const double tcol[/*1+m*/], int *p_flag,
      double tol_piv, double tol, double tol1)
{     int m = lp->m;
      int n = lp->n;
      double *c = lp->c;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      int i, i_flag, k, p;
      double alfa, biga, delta, lk, uk, teta, teta_min;
      xassert(phase == 1 || phase == 2);
      xassert(1 <= q && q <= n-m);
      xassert(s == +1.0 || s == -1.0);
      /*--------------------------------------------------------------*/
      /* first pass: determine teta_min for relaxed bounds            */
      /*--------------------------------------------------------------*/
      teta_min = DBL_MAX;
      /* walk thru the list of basic variables */
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
         /* determine alfa such that delta xB[i] = alfa * teta */
         alfa = s * tcol[i];
         if (alfa <= -tol_piv)
         {  /* xB[i] decreases */
            /* determine actual lower bound of xB[i] */
            if (phase == 1 && c[k] < 0.0)
            {  /* xB[i] has no actual lower bound */
               continue;
            }
            else if (phase == 1 && c[k] > 0.0)
            {  /* actual lower bound of xB[i] is its upper bound */
               lk = u[k];
               xassert(lk != +DBL_MAX);
            }
            else
            {  /* actual lower bound of xB[i] is its original bound */
               lk = l[k];
               if (lk == -DBL_MAX)
                  continue;
            }
            /* determine teta on which xB[i] reaches its relaxed lower
             * bound */
            delta = tol + tol1 * (lk >= 0.0 ? +lk : -lk);
            if (beta[i] < lk)
               teta = - delta / alfa;
            else
               teta = ((lk - delta) - beta[i]) / alfa;
         }
         else if (alfa >= +tol_piv)
         {  /* xB[i] increases */
            /* determine actual upper bound of xB[i] */
            if (phase == 1 && c[k] < 0.0)
            {  /* actual upper bound of xB[i] is its lower bound */
               uk = l[k];
               xassert(uk != -DBL_MAX);
            }
            else if (phase == 1 && c[k] > 0.0)
            {  /* xB[i] has no actual upper bound */
               continue;
            }
            else
            {  /* actual upper bound of xB[i] is its original bound */
               uk = u[k];
               if (uk == +DBL_MAX)
                  continue;
            }
            /* determine teta on which xB[i] reaches its relaxed upper
             * bound */
            delta = tol + tol1 * (uk >= 0.0 ? +uk : -uk);
            if (beta[i] > uk)
               teta = + delta / alfa;
            else
               teta = ((uk + delta) - beta[i]) / alfa;
         }
         else
         {  /* xB[i] does not depend on teta */
            continue;
         }
         xassert(teta >= 0.0);
         if (teta_min > teta)
            teta_min = teta;
      }
      /*--------------------------------------------------------------*/
      /* second pass: choose basic variable xB[p]                     */
      /*--------------------------------------------------------------*/
      k = head[m+q]; /* x[k] = xN[q] */
      if (l[k] != -DBL_MAX && u[k] != +DBL_MAX)
      {  /* xN[q] has both lower and upper bounds */
         if (fabs(l[k] - u[k]) <= teta_min)
         {  /* and reaches its opposite bound */
            p = -1, *p_flag = 0;
            goto done;
         }
      }
      if (teta_min == DBL_MAX)
      {  /* teta may increase unlimitedly */
         p = 0, *p_flag = 0;
         goto done;
      }
      /* nothing is chosen so far */
      p = 0, *p_flag = 0, biga = 0.0;
      /* walk thru the list of basic variables */
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
         /* determine alfa such that delta xB[i] = alfa * teta */
         alfa = s * tcol[i];
         if (alfa <= -tol_piv)
         {  /* xB[i] decreases */
            /* determine actual lower bound of xB[i] */
            if (phase == 1 && c[k] < 0.0)
            {  /* xB[i] has no actual lower bound */
               continue;
            }
            else if (phase == 1 && c[k] > 0.0)
            {  /* actual lower bound of xB[i] is its upper bound */
               lk = u[k];
               xassert(lk != +DBL_MAX);
               i_flag = 1;
            }
            else
            {  /* actual lower bound of xB[i] is its original bound */
               lk = l[k];
               if (lk == -DBL_MAX)
                  continue;
               i_flag = 0;
            }
            /* determine teta on which xB[i] reaches its lower bound */
            teta = (lk - beta[i]) / alfa;
         }
         else if (alfa >= +tol_piv)
         {  /* xB[i] increases */
            /* determine actual upper bound of xB[i] */
            if (phase == 1 && c[k] < 0.0)
            {  /* actual upper bound of xB[i] is its lower bound */
               uk = l[k];
               xassert(uk != -DBL_MAX);
               i_flag = 0;
            }
            else if (phase == 1 && c[k] > 0.0)
            {  /* xB[i] has no actual upper bound */
               continue;
            }
            else
            {  /* actual upper bound of xB[i] is its original bound */
               uk = u[k];
               if (uk == +DBL_MAX)
                  continue;
               i_flag = 1;
            }
            /* determine teta on which xB[i] reaches its upper bound */
            teta = (uk - beta[i]) / alfa;
         }
         else
         {  /* xB[i] does not depend on teta */
            continue;
         }
         /* choose basic variable for which teta is not greater than
          * teta_min determined for relaxed bounds and which has best
          * (largest in magnitude) pivot */
         alfa = (alfa >= 0.0 ? +alfa : -alfa);
         if (teta <= teta_min && biga < alfa)
            p = i, *p_flag = i_flag, biga = alfa;
      }
      /* something must be chosen */
      xassert(1 <= p && p <= m);
      /* if xB[p] is fixed variable, adjust its bound flag */
      k = head[p];
      if (l[k] == u[k])
         *p_flag = 0;
done: return p;
}

#if 1 /* 22/VI-2017 */
/***********************************************************************
*  spx_ls_eval_bp - determine penalty function break points
*
*  This routine determines break points of the penalty function (which
*  is the sum of primal infeasibilities).
*
*  The parameters lp, beta, q, dq, tcol, and tol_piv have the same
*  meaning as for the routine spx_chuzr_std (see above).
*
*  The routine stores the break-points determined to the array elements
*  bp[1], ..., bp[nbp] in *arbitrary* order, where 0 <= nbp <= 2*m+1 is
*  the number of break-points returned by the routine on exit. */

int spx_ls_eval_bp(SPXLP *lp, const double beta[/*1+m*/],
      int q, double dq, const double tcol[/*1+m*/], double tol_piv,
      SPXBP bp[/*1+2*m+1*/])
{     int m = lp->m;
      int n = lp->n;
      double *c = lp->c;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      int i, k, nbp;
      double s, alfa;
      xassert(1 <= q && q <= n-m);
      xassert(dq != 0.0);
      s = (dq < 0.0 ? +1.0 : -1.0);
      nbp = 0;
      /* if chosen non-basic variable xN[q] is double-bounded, include
       * it in the list, because it can cross its opposite bound */
      k = head[m+q]; /* x[k] = xN[q] */
      if (l[k] != -DBL_MAX && u[k] != +DBL_MAX)
      {  nbp++;
         bp[nbp].i = 0;
         xassert(l[k] < u[k]); /* xN[q] cannot be fixed */
         bp[nbp].teta = u[k] - l[k];
         bp[nbp].dc = s;
      }
      /* build the list of all basic variables xB[i] that can cross
       * their bound(s) for the ray parameter 0 <= teta < teta_max */
      for (i = 1; i <= m; i++)
      {  k = head[i]; /* x[k] = xB[i] */
         xassert(l[k] <= u[k]);
         /* determine alfa such that (delta xB[i]) = alfa * teta */
         alfa = s * tcol[i];
         if (alfa >= +tol_piv)
         {  /* xB[i] increases on increasing teta */
            if (l[k] == u[k])
            {  /* xB[i] is fixed at lB[i] = uB[i] */
               if (c[k] <= 0.0)
               {  /* increasing xB[i] can cross its fixed value lB[i],
                   * because currently xB[i] <= lB[i] */
                  nbp++;
                  bp[nbp].i = +i;
                  bp[nbp].teta = (l[k] - beta[i]) / alfa;
                  /* if xB[i] > lB[i] then cB[i] = +1 */
                  bp[nbp].dc = +1.0 - c[k];
               }
            }
            else
            {  if (l[k] != -DBL_MAX && c[k] < 0.0)
               {  /* increasing xB[i] can cross its lower bound lB[i],
                   * because currently xB[i] < lB[i] */
                  nbp++;
                  bp[nbp].i = +i;
                  bp[nbp].teta = (l[k] - beta[i]) / alfa;
                  bp[nbp].dc = +1.0;
               }
               if (u[k] != +DBL_MAX && c[k] <= 0.0)
               {  /* increasing xB[i] can cross its upper bound uB[i],
                   * because currently xB[i] does not violate it */
                  nbp++;
                  bp[nbp].i = -i;
                  bp[nbp].teta = (u[k] - beta[i]) / alfa;
                  bp[nbp].dc = +1.0;
               }
            }
         }
         else if (alfa <= -tol_piv)
         {  /* xB[i] decreases on increasing teta */
            if (l[k] == u[k])
            {  /* xB[i] is fixed at lB[i] = uB[i] */
               if (c[k] >= 0.0)
               {  /* decreasing xB[i] can cross its fixed value lB[i],
                   * because currently xB[i] >= lB[i] */
                  nbp++;
                  bp[nbp].i = +i;
                  bp[nbp].teta = (l[k] - beta[i]) / alfa;
                  /* if xB[i] < lB[i] then cB[i] = -1 */
                  bp[nbp].dc = -1.0 - c[k];
               }
            }
            else
            {  if (l[k] != -DBL_MAX && c[k] >= 0.0)
               {  /* decreasing xB[i] can cross its lower bound lB[i],
                   * because currently xB[i] does not violate it */
                  nbp++;
                  bp[nbp].i = +i;
                  bp[nbp].teta = (l[k] - beta[i]) / alfa;
                  bp[nbp].dc = -1.0;
               }
               if (u[k] != +DBL_MAX && c[k] > 0.0)
               {  /* decreasing xB[i] can cross its upper bound uB[i],
                   * because currently xB[i] > uB[i] */
                  nbp++;
                  bp[nbp].i = -i;
                  bp[nbp].teta = (u[k] - beta[i]) / alfa;
                  bp[nbp].dc = -1.0;
               }
            }
         }
         else
         {  /* xB[i] does not depend on teta within a tolerance */
            continue;
         }
         /* teta < 0 may happen only due to round-off errors when the
          * current value of xB[i] is *close* to its (lower or upper)
          * bound; in this case we replace teta by exact zero */
         if (bp[nbp].teta < 0.0)
            bp[nbp].teta = 0.0;
      }
      xassert(nbp <= 2*m+1);
      return nbp;
}
#endif

#if 1 /* 22/VI-2017 */
/***********************************************************************
*  spx_ls_select_bp - select and process penalty function break points
*
*  This routine selects a next portion of the penalty function break
*  points and processes them.
*
*  On entry to the routine it is assumed that break points bp[1], ...,
*  bp[num] are already processed, and slope is the penalty function
*  slope to the right of the last processed break point bp[num].
*  (Initially, when num = 0, slope should be specified as -fabs(d[q]),
*  where d[q] is the reduced cost of chosen non-basic variable xN[q].)
*
*  The routine selects break points among bp[num+1], ..., bp[nbp], for
*  which teta <= teta_lim, and moves these break points to the array
*  elements bp[num+1], ..., bp[num1], where num <= num1 <= 2*m+1 is the
*  new number of processed break points returned by the routine on
*  exit. Then the routine sorts the break points by ascending teta and
*  computes the change of the penalty function relative to its value at
*  teta = 0.
*
*  On exit the routine also replaces the parameter slope with a new
*  value that corresponds to the new last break-point bp[num1]. */

static int CDECL fcmp(const void *v1, const void *v2)
{     const SPXBP *p1 = v1, *p2 = v2;
      if (p1->teta < p2->teta)
         return -1;
      else if (p1->teta > p2->teta)
         return +1;
      else
         return 0;
}

int spx_ls_select_bp(SPXLP *lp, const double tcol[/*1+m*/],
      int nbp, SPXBP bp[/*1+m+m+1*/], int num, double *slope, double
      teta_lim)
{     int m = lp->m;
      int i, t, num1;
      double teta, dz;
      xassert(0 <= num && num <= nbp && nbp <= m+m+1);
      /* select a new portion of break points */
      num1 = num;
      for (t = num+1; t <= nbp; t++)
      {  if (bp[t].teta <= teta_lim)
         {  /* move break point to the beginning of the new portion */
            num1++;
            i = bp[num1].i, teta = bp[num1].teta, dz = bp[num1].dc;
            bp[num1].i = bp[t].i, bp[num1].teta = bp[t].teta,
            bp[num1].dc = bp[t].dc;
            bp[t].i = i, bp[t].teta = teta, bp[t].dc = dz;
         }
      }
      /* sort new break points bp[num+1], ..., bp[num1] by ascending
       * the ray parameter teta */
      if (num1 - num > 1)
         qsort(&bp[num+1], num1 - num, sizeof(SPXBP), fcmp);
      /* calculate the penalty function change at the new break points
       * selected */
      for (t = num+1; t <= num1; t++)
      {  /* calculate the penalty function change relative to its value
          * at break point bp[t-1] */
         dz = (*slope) * (bp[t].teta - (t == 1 ? 0.0 : bp[t-1].teta));
         /* calculate the penalty function change relative to its value
          * at teta = 0 */
         bp[t].dz = (t == 1 ? 0.0 : bp[t-1].dz) + dz;
         /* calculate a new slope of the penalty function to the right
          * of the current break point bp[t] */
         i = (bp[t].i >= 0 ? bp[t].i : -bp[t].i);
         xassert(0 <= i && i <= m);
         if (i == 0)
            *slope += fabs(1.0 * bp[t].dc);
         else
            *slope += fabs(tcol[i] * bp[t].dc);
      }
      return num1;
}
#endif

/* eof */
