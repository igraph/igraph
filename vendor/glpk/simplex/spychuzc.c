/* spychuzc.c */

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
#include "spychuzc.h"

/***********************************************************************
*  spy_chuzc_std - choose non-basic variable (dual textbook ratio test)
*
*  This routine implements an improved dual textbook ratio test to
*  choose non-basic variable xN[q].
*
*  Current reduced costs of non-basic variables should be placed in the
*  array locations d[1], ..., d[n-m]. Note that d[j] is a value of dual
*  basic variable lambdaN[j] in the current basis.
*
#if 0 (* 14/III-2016 *)
*  The parameter s specifies the sign of bound violation for basic
*  variable xB[p] chosen: s = +1.0 means that xB[p] violates its lower
*  bound, so dual non-basic variable lambdaB[p] = lambda^+B[p]
*  increases, and s = -1.0 means that xB[p] violates its upper bound,
*  so dual non-basic variable lambdaB[p] = lambda^-B[p] decreases.
*  (Thus, the dual ray parameter theta = s * lambdaB[p] >= 0.)
#else
*  The parameter r specifies the bound violation for basic variable
*  xB[p] chosen:
*
*  r = lB[p] - beta[p] > 0 means that xB[p] violates its lower bound,
*  so dual non-basic variable lambdaB[p] = lambda^+B[p] increases; and
*
*  r = uB[p] - beta[p] < 0 means that xB[p] violates its upper bound,
*  so dual non-basic variable lambdaB[p] = lambda^-B[p] decreases.
*
*  (Note that r is the dual reduced cost of lambdaB[p].)
#endif
*
*  Elements of p-th simplex table row t[p] = (t[p,j]) corresponding
*  to basic variable xB[p] should be placed in the array locations
*  trow[1], ..., trow[n-m].
*
*  The parameter tol_piv specifies a tolerance for elements of the
*  simplex table row t[p]. If |t[p,j]| < tol_piv, dual basic variable
*  lambdaN[j] is skipped, i.e. it is assumed that it does not depend on
*  the dual ray parameter theta.
*
*  The parameters tol and tol1 specify tolerances used to increase the
*  choice freedom by simulating an artificial degeneracy as follows.
*  If lambdaN[j] = lambda^+N[j] >= 0 and d[j] <= +delta[j], or if
*  lambdaN[j] = lambda^-N[j] <= 0 and d[j] >= -delta[j], where
*  delta[j] = tol + tol1 * |cN[j]|, cN[j] is objective coefficient at
*  xN[j], then it is assumed that reduced cost d[j] is equal to zero.
*
*  The routine determines the index 1 <= q <= n-m of non-basic variable
*  xN[q], for which corresponding dual basic variable lambda^+N[j] or
*  lambda^-N[j] reaches its zero bound first on increasing the dual ray
*  parameter theta, and returns p on exit. And if theta may increase
*  unlimitedly, the routine returns zero. */

int spy_chuzc_std(SPXLP *lp, const double d[/*1+n-m*/],
#if 0 /* 14/III-2016 */
      double s, const double trow[/*1+n-m*/], double tol_piv,
#else
      double r, const double trow[/*1+n-m*/], double tol_piv,
#endif
      double tol, double tol1)
{     int m = lp->m;
      int n = lp->n;
      double *c = lp->c;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      int j, k, q;
      double alfa, biga, delta, teta, teta_min;
#if 0 /* 14/III-2016 */
      xassert(s == +1.0 || s == -1.0);
#else
      double s;
      xassert(r != 0.0);
      s = (r > 0.0 ? +1.0 : -1.0);
#endif
      /* nothing is chosen so far */
      q = 0, teta_min = DBL_MAX, biga = 0.0;
      /* walk thru the list of non-basic variables */
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         /* if xN[j] is fixed variable, skip it */
         if (l[k] == u[k])
            continue;
         alfa = s * trow[j];
         if (alfa >= +tol_piv && !flag[j])
         {  /* xN[j] is either free or has its lower bound active, so
             * lambdaN[j] = d[j] >= 0 decreases down to zero */
            delta = tol + tol1 * (c[k] >= 0.0 ? +c[k] : -c[k]);
            /* determine theta on which lambdaN[j] reaches zero */
            teta = (d[j] < +delta ? 0.0 : d[j] / alfa);
         }
         else if (alfa <= -tol_piv && (l[k] == -DBL_MAX || flag[j]))
         {  /* xN[j] is either free or has its upper bound active, so
             * lambdaN[j] = d[j] <= 0 increases up to zero */
            delta = tol + tol1 * (c[k] >= 0.0 ? +c[k] : -c[k]);
            /* determine theta on which lambdaN[j] reaches zero */
            teta = (d[j] > -delta ? 0.0 : d[j] / alfa);
         }
         else
         {  /* lambdaN[j] cannot reach zero on increasing theta */
            continue;
         }
         /* choose non-basic variable xN[q] by corresponding dual basic
          * variable lambdaN[q] for which theta is minimal */
         xassert(teta >= 0.0);
         alfa = (alfa >= 0.0 ? +alfa : -alfa);
         if (teta_min > teta || (teta_min == teta && biga < alfa))
            q = j, teta_min = teta, biga = alfa;
      }
      return q;
}

/***********************************************************************
*  spy_chuzc_harris - choose non-basic var. (dual Harris' ratio test)
*
*  This routine implements dual Harris' ratio test to choose non-basic
*  variable xN[q].
*
*  All the parameters, except tol and tol1, as well as the returned
*  value have the same meaning as for the routine spx_chuzr_std (see
*  above).
*
*  The parameters tol and tol1 specify tolerances on zero bound
*  violations for reduced costs of non-basic variables. For reduced
*  cost d[j] the tolerance is delta[j] = tol + tol1 |cN[j]|, where
*  cN[j] is objective coefficient at non-basic variable xN[j]. */

int spy_chuzc_harris(SPXLP *lp, const double d[/*1+n-m*/],
#if 0 /* 14/III-2016 */
      double s, const double trow[/*1+n-m*/], double tol_piv,
#else
      double r, const double trow[/*1+n-m*/], double tol_piv,
#endif
      double tol, double tol1)
{     int m = lp->m;
      int n = lp->n;
      double *c = lp->c;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      int j, k, q;
      double alfa, biga, delta, teta, teta_min;
#if 0 /* 14/III-2016 */
      xassert(s == +1.0 || s == -1.0);
#else
      double s;
      xassert(r != 0.0);
      s = (r > 0.0 ? +1.0 : -1.0);
#endif
      /*--------------------------------------------------------------*/
      /* first pass: determine teta_min for relaxed bounds            */
      /*--------------------------------------------------------------*/
      teta_min = DBL_MAX;
      /* walk thru the list of non-basic variables */
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         /* if xN[j] is fixed variable, skip it */
         if (l[k] == u[k])
            continue;
         alfa = s * trow[j];
         if (alfa >= +tol_piv && !flag[j])
         {  /* xN[j] is either free or has its lower bound active, so
             * lambdaN[j] = d[j] >= 0 decreases down to zero */
            delta = tol + tol1 * (c[k] >= 0.0 ? +c[k] : -c[k]);
            /* determine theta on which lambdaN[j] reaches -delta */
            teta = ((d[j] < 0.0 ? 0.0 : d[j]) + delta) / alfa;
         }
         else if (alfa <= -tol_piv && (l[k] == -DBL_MAX || flag[j]))
         {  /* xN[j] is either free or has its upper bound active, so
             * lambdaN[j] = d[j] <= 0 increases up to zero */
            delta = tol + tol1 * (c[k] >= 0.0 ? +c[k] : -c[k]);
            /* determine theta on which lambdaN[j] reaches +delta */
            teta = ((d[j] > 0.0 ? 0.0 : d[j]) - delta) / alfa;
         }
         else
         {  /* lambdaN[j] cannot reach zero on increasing theta */
            continue;
         }
         xassert(teta >= 0.0);
         if (teta_min > teta)
            teta_min = teta;
      }
      /*--------------------------------------------------------------*/
      /* second pass: choose non-basic variable xN[q]                 */
      /*--------------------------------------------------------------*/
      if (teta_min == DBL_MAX)
      {  /* theta may increase unlimitedly */
         q = 0;
         goto done;
      }
      /* nothing is chosen so far */
      q = 0, biga = 0.0;
      /* walk thru the list of non-basic variables */
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         /* if xN[j] is fixed variable, skip it */
         if (l[k] == u[k])
            continue;
         alfa = s * trow[j];
         if (alfa >= +tol_piv && !flag[j])
         {  /* xN[j] is either free or has its lower bound active, so
             * lambdaN[j] = d[j] >= 0 decreases down to zero */
            /* determine theta on which lambdaN[j] reaches zero */
            teta = d[j] / alfa;
         }
         else if (alfa <= -tol_piv && (l[k] == -DBL_MAX || flag[j]))
         {  /* xN[j] is either free or has its upper bound active, so
             * lambdaN[j] = d[j] <= 0 increases up to zero */
            /* determine theta on which lambdaN[j] reaches zero */
            teta = d[j] / alfa;
         }
         else
         {  /* lambdaN[j] cannot reach zero on increasing theta */
            continue;
         }
         /* choose non-basic variable for which theta is not greater
          * than theta_min determined for relaxed bounds and which has
          * best (largest in magnitude) pivot */
         alfa = (alfa >= 0.0 ? +alfa : -alfa);
         if (teta <= teta_min && biga < alfa)
            q = j, biga = alfa;
      }
      /* something must be chosen */
      xassert(1 <= q && q <= n-m);
done: return q;
}

#if 0 /* 23/III-2016 */
/***********************************************************************
*  spy_eval_bp - determine dual objective function break-points
*
*  This routine determines the dual objective function break-points.
*
*  The parameters lp, d, r, trow, and tol_piv have the same meaning as
*  for the routine spx_chuzc_std (see above).
*
*  On exit the routine stores the break-points determined to the array
*  elements bp[1], ..., bp[num], where 0 <= num <= n-m is the number of
*  break-points returned by the routine.
*
*  The break-points stored in the array bp are ordered by ascending
*  the ray parameter teta >= 0. The break-points numbered 1, ..., num-1
*  always correspond to non-basic non-fixed variables xN[j] of primal
*  LP having both lower and upper bounds while the last break-point
*  numbered num may correspond to a non-basic variable having only one
*  lower or upper bound, if such variable prevents further increasing
*  of the ray parameter teta. Besides, the routine includes in the
*  array bp only the break-points that correspond to positive increment
*  of the dual objective. */

static int CDECL fcmp(const void *v1, const void *v2)
{     const SPYBP *p1 = v1, *p2 = v2;
      if (p1->teta < p2->teta)
         return -1;
      else if (p1->teta > p2->teta)
         return +1;
      else
         return 0;
}

int spy_eval_bp(SPXLP *lp, const double d[/*1+n-m*/],
      double r, const double trow[/*1+n-m*/], double tol_piv,
      SPYBP bp[/*1+n-m*/])
{     int m = lp->m;
      int n = lp->n;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      int j, j_max, k, t, nnn, num;
      double s, alfa, teta, teta_max, dz, v;
      xassert(r != 0.0);
      s = (r > 0.0 ? +1.0 : -1.0);
      /* build the list of all dual basic variables lambdaN[j] that
       * can reach zero on increasing the ray parameter teta >= 0 */
      num = 0;
      /* walk thru the list of non-basic variables */
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         /* if xN[j] is fixed variable, skip it */
         if (l[k] == u[k])
            continue;
         alfa = s * trow[j];
         if (alfa >= +tol_piv && !flag[j])
         {  /* xN[j] is either free or has its lower bound active, so
             * lambdaN[j] = d[j] >= 0 decreases down to zero */
            /* determine teta[j] on which lambdaN[j] reaches zero */
            teta = (d[j] < 0.0 ? 0.0 : d[j] / alfa);
         }
         else if (alfa <= -tol_piv && (l[k] == -DBL_MAX || flag[j]))
         {  /* xN[j] is either free or has its upper bound active, so
             * lambdaN[j] = d[j] <= 0 increases up to zero */
            /* determine teta[j] on which lambdaN[j] reaches zero */
            teta = (d[j] > 0.0 ? 0.0 : d[j] / alfa);
         }
         else
         {  /* lambdaN[j] cannot reach zero on increasing teta */
            continue;
         }
         /* add lambdaN[j] to the list */
         num++;
         bp[num].j = j;
         bp[num].teta = teta;
      }
      if (num == 0)
      {  /* dual unboundedness */
         goto done;
      }
      /* determine "blocking" dual basic variable lambdaN[j_max] that
       * prevents increasing teta more than teta_max */
      j_max = 0, teta_max = DBL_MAX;
      for (t = 1; t <= num; t++)
      {  j = bp[t].j;
         k = head[m+j]; /* x[k] = xN[j] */
         if (l[k] == -DBL_MAX || u[k] == +DBL_MAX)
         {  /* lambdaN[j] cannot intersect zero */
            if (j_max == 0
               || teta_max > bp[t].teta
               || (teta_max == bp[t].teta
                  && fabs(trow[j_max]) < fabs(trow[j])))
               j_max = j, teta_max = bp[t].teta;
         }
      }
      /* keep in the list only dual basic variables lambdaN[j] that
       * correspond to primal double-bounded variables xN[j] and whose
       * teta[j] is not greater than teta_max */
      nnn = 0;
      for (t = 1; t <= num; t++)
      {  j = bp[t].j;
         k = head[m+j]; /* x[k] = xN[j] */
         if (l[k] != -DBL_MAX && u[k] != +DBL_MAX
            && bp[t].teta <= teta_max)
         {  nnn++;
            bp[nnn].j = j;
            bp[nnn].teta = bp[t].teta;
         }
      }
      num = nnn;
      /* sort break-points by ascending teta[j] */
      qsort(&bp[1], num, sizeof(SPYBP), fcmp);
      /* add lambdaN[j_max] to the end of the list */
      if (j_max != 0)
      {  xassert(num < n-m);
         num++;
         bp[num].j = j_max;
         bp[num].teta = teta_max;
      }
      /* compute increments of the dual objective at all break-points
       * (relative to its value at teta = 0) */
      dz = 0.0;      /* dual objective increment */
      v = fabs(r);   /* dual objective slope d zeta / d teta */
      for (t = 1; t <= num; t++)
      {  /* compute increment at current break-point */
         dz += v * (bp[t].teta - (t == 1 ? 0.0 : bp[t-1].teta));
         if (dz < 0.001)
         {  /* break-point with non-positive increment reached */
            num = t - 1;
            break;
         }
         bp[t].dz = dz;
         /* compute next slope on the right to current break-point */
         if (t < num)
         {  j = bp[t].j;
            k = head[m+j]; /* x[k] = xN[j] */
            xassert(-DBL_MAX < l[k] && l[k] < u[k] && u[k] < +DBL_MAX);
            v -= fabs(trow[j]) * (u[k] - l[k]);
         }
      }
done: return num;
}
#endif

/***********************************************************************
*  spy_ls_eval_bp - determine dual objective function break-points
*
*  This routine determines the dual objective function break-points.
*
*  The parameters lp, d, r, trow, and tol_piv have the same meaning as
*  for the routine spx_chuzc_std (see above).
*
*  The routine stores the break-points determined to the array elements
*  bp[1], ..., bp[nbp] in *arbitrary* order, where 0 <= nbp <= n-m is
*  the number of break-points returned by the routine on exit. */

int spy_ls_eval_bp(SPXLP *lp, const double d[/*1+n-m*/],
      double r, const double trow[/*1+n-m*/], double tol_piv,
      SPYBP bp[/*1+n-m*/])
{     int m = lp->m;
      int n = lp->n;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      char *flag = lp->flag;
      int j, k, t, nnn, nbp;
      double s, alfa, teta, teta_max;
      xassert(r != 0.0);
      s = (r > 0.0 ? +1.0 : -1.0);
      /* build the list of all dual basic variables lambdaN[j] that
       * can reach zero on increasing the ray parameter teta >= 0 */
      nnn = 0, teta_max = DBL_MAX;
      /* walk thru the list of non-basic variables */
      for (j = 1; j <= n-m; j++)
      {  k = head[m+j]; /* x[k] = xN[j] */
         /* if xN[j] is fixed variable, skip it */
         if (l[k] == u[k])
            continue;
         alfa = s * trow[j];
         if (alfa >= +tol_piv && !flag[j])
         {  /* xN[j] is either free or has its lower bound active, so
             * lambdaN[j] = d[j] >= 0 decreases down to zero */
            /* determine teta[j] on which lambdaN[j] reaches zero */
            teta = (d[j] < 0.0 ? 0.0 : d[j] / alfa);
            /* if xN[j] has no upper bound, lambdaN[j] cannot become
             * negative and thereby blocks further increasing teta */
            if (u[k] == +DBL_MAX && teta_max > teta)
               teta_max = teta;
         }
         else if (alfa <= -tol_piv && (l[k] == -DBL_MAX || flag[j]))
         {  /* xN[j] is either free or has its upper bound active, so
             * lambdaN[j] = d[j] <= 0 increases up to zero */
            /* determine teta[j] on which lambdaN[j] reaches zero */
            teta = (d[j] > 0.0 ? 0.0 : d[j] / alfa);
            /* if xN[j] has no lower bound, lambdaN[j] cannot become
             * positive and thereby blocks further increasing teta */
            if (l[k] == -DBL_MAX && teta_max > teta)
               teta_max = teta;
         }
         else
         {  /* lambdaN[j] cannot reach zero on increasing teta */
            continue;
         }
         /* add lambdaN[j] to the list */
         nnn++;
         bp[nnn].j = j;
         bp[nnn].teta = teta;
      }
      /* remove from the list all dual basic variables lambdaN[j], for
       * which teta[j] > teta_max */
      nbp = 0;
      for (t = 1; t <= nnn; t++)
      {  if (bp[t].teta <= teta_max + 1e-6)
         {  nbp++;
            bp[nbp].j = bp[t].j;
            bp[nbp].teta = bp[t].teta;
         }
      }
      return nbp;
}

/***********************************************************************
*  spy_ls_select_bp - select and process dual objective break-points
*
*  This routine selects a next portion of the dual objective function
*  break-points and processes them.
*
*  On entry to the routine it is assumed that break-points bp[1], ...,
*  bp[num] are already processed, and slope is the dual objective slope
*  to the right of the last processed break-point bp[num]. (Initially,
*  when num = 0, slope should be specified as fabs(r), where r has the
*  same meaning as above.)
*
*  The routine selects break-points among bp[num+1], ..., bp[nbp], for
*  which teta <= teta_lim, and moves these break-points to the array
*  elements bp[num+1], ..., bp[num1], where num <= num1 <= n-m is the
*  new number of processed break-points returned by the routine on
*  exit. Then the routine sorts these break-points by ascending teta
*  and computes the change of the dual objective function relative to
*  its value at teta = 0.
*
*  On exit the routine also replaces the parameter slope with a new
*  value that corresponds to the new last break-point bp[num1]. */

static int CDECL fcmp(const void *v1, const void *v2)
{     const SPYBP *p1 = v1, *p2 = v2;
      if (p1->teta < p2->teta)
         return -1;
      else if (p1->teta > p2->teta)
         return +1;
      else
         return 0;
}

int spy_ls_select_bp(SPXLP *lp, const double trow[/*1+n-m*/],
      int nbp, SPYBP bp[/*1+n-m*/], int num, double *slope, double
      teta_lim)
{     int m = lp->m;
      int n = lp->n;
      double *l = lp->l;
      double *u = lp->u;
      int *head = lp->head;
      int j, k, t, num1;
      double teta, dz;
      xassert(0 <= num && num <= nbp && nbp <= n-m);
      /* select a new portion of break-points */
      num1 = num;
      for (t = num+1; t <= nbp; t++)
      {  if (bp[t].teta <= teta_lim)
         {  /* move break-point to the beginning of the new portion */
            num1++;
            j = bp[num1].j, teta = bp[num1].teta;
            bp[num1].j = bp[t].j, bp[num1].teta = bp[t].teta;
            bp[t].j = j, bp[t].teta = teta;
         }
      }
      /* sort new break-points bp[num+1], ..., bp[num1] by ascending
       * the ray parameter teta */
      if (num1 - num > 1)
         qsort(&bp[num+1], num1 - num, sizeof(SPYBP), fcmp);
      /* calculate the dual objective change at the new break-points */
      for (t = num+1; t <= num1; t++)
      {  /* calculate the dual objective change relative to its value
          * at break-point bp[t-1] */
         if (*slope == -DBL_MAX)
            dz = -DBL_MAX;
         else
            dz = (*slope) *
               (bp[t].teta - (t == 1 ? 0.0 : bp[t-1].teta));
         /* calculate the dual objective change relative to its value
          * at teta = 0 */
         if (dz == -DBL_MAX)
            bp[t].dz = -DBL_MAX;
         else
            bp[t].dz = (t == 1 ? 0.0 : bp[t-1].dz) + dz;
         /* calculate a new slope of the dual objective to the right of
          * the current break-point bp[t] */
         if (*slope != -DBL_MAX)
         {  j = bp[t].j;
            k = head[m+j]; /* x[k] = xN[j] */
            if (l[k] == -DBL_MAX || u[k] == +DBL_MAX)
               *slope = -DBL_MAX; /* blocking break-point reached */
            else
            {  xassert(l[k] < u[k]);
               *slope -= fabs(trow[j]) * (u[k] - l[k]);
            }
         }
      }
      return num1;
}

/* eof */
