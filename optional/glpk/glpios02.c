/* glpios02.c (preprocess current subproblem) */

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

#include "glpios.h"

/***********************************************************************
*  prepare_row_info - prepare row info to determine implied bounds
*
*  Given a row (linear form)
*
*      n
*     sum a[j] * x[j]                                                (1)
*     j=1
*
*  and bounds of columns (variables)
*
*     l[j] <= x[j] <= u[j]                                           (2)
*
*  this routine computes f_min, j_min, f_max, j_max needed to determine
*  implied bounds.
*
*  ALGORITHM
*
*  Let J+ = {j : a[j] > 0} and J- = {j : a[j] < 0}.
*
*  Parameters f_min and j_min are computed as follows:
*
*  1) if there is no x[k] such that k in J+ and l[k] = -inf or k in J-
*     and u[k] = +inf, then
*
*     f_min :=   sum   a[j] * l[j] +   sum   a[j] * u[j]
*              j in J+               j in J-
*                                                                    (3)
*     j_min := 0
*
*  2) if there is exactly one x[k] such that k in J+ and l[k] = -inf
*     or k in J- and u[k] = +inf, then
*
*     f_min :=   sum       a[j] * l[j] +   sum       a[j] * u[j]
*              j in J+\{k}               j in J-\{k}
*                                                                    (4)
*     j_min := k
*
*  3) if there are two or more x[k] such that k in J+ and l[k] = -inf
*     or k in J- and u[k] = +inf, then
*
*     f_min := -inf
*                                                                    (5)
*     j_min := 0
*
*  Parameters f_max and j_max are computed in a similar way as follows:
*
*  1) if there is no x[k] such that k in J+ and u[k] = +inf or k in J-
*     and l[k] = -inf, then
*
*     f_max :=   sum   a[j] * u[j] +   sum   a[j] * l[j]
*              j in J+               j in J-
*                                                                    (6)
*     j_max := 0
*
*  2) if there is exactly one x[k] such that k in J+ and u[k] = +inf
*     or k in J- and l[k] = -inf, then
*
*     f_max :=   sum       a[j] * u[j] +   sum       a[j] * l[j]
*              j in J+\{k}               j in J-\{k}
*                                                                    (7)
*     j_max := k
*
*  3) if there are two or more x[k] such that k in J+ and u[k] = +inf
*     or k in J- and l[k] = -inf, then
*
*     f_max := +inf
*                                                                    (8)
*     j_max := 0                                                      */

struct f_info
{     int j_min, j_max;
      double f_min, f_max;
};

static void prepare_row_info(int n, const double a[], const double l[],
      const double u[], struct f_info *f)
{     int j, j_min, j_max;
      double f_min, f_max;
      xassert(n >= 0);
      /* determine f_min and j_min */
      f_min = 0.0, j_min = 0;
      for (j = 1; j <= n; j++)
      {  if (a[j] > 0.0)
         {  if (l[j] == -DBL_MAX)
            {  if (j_min == 0)
                  j_min = j;
               else
               {  f_min = -DBL_MAX, j_min = 0;
                  break;
               }
            }
            else
               f_min += a[j] * l[j];
         }
         else if (a[j] < 0.0)
         {  if (u[j] == +DBL_MAX)
            {  if (j_min == 0)
                  j_min = j;
               else
               {  f_min = -DBL_MAX, j_min = 0;
                  break;
               }
            }
            else
               f_min += a[j] * u[j];
         }
         else
            xassert(a != a);
      }
      f->f_min = f_min, f->j_min = j_min;
      /* determine f_max and j_max */
      f_max = 0.0, j_max = 0;
      for (j = 1; j <= n; j++)
      {  if (a[j] > 0.0)
         {  if (u[j] == +DBL_MAX)
            {  if (j_max == 0)
                  j_max = j;
               else
               {  f_max = +DBL_MAX, j_max = 0;
                  break;
               }
            }
            else
               f_max += a[j] * u[j];
         }
         else if (a[j] < 0.0)
         {  if (l[j] == -DBL_MAX)
            {  if (j_max == 0)
                  j_max = j;
               else
               {  f_max = +DBL_MAX, j_max = 0;
                  break;
               }
            }
            else
               f_max += a[j] * l[j];
         }
         else
            xassert(a != a);
      }
      f->f_max = f_max, f->j_max = j_max;
      return;
}

/***********************************************************************
*  row_implied_bounds - determine row implied bounds
*
*  Given a row (linear form)
*
*      n
*     sum a[j] * x[j]
*     j=1
*
*  and bounds of columns (variables)
*
*     l[j] <= x[j] <= u[j]
*
*  this routine determines implied bounds of the row.
*
*  ALGORITHM
*
*  Let J+ = {j : a[j] > 0} and J- = {j : a[j] < 0}.
*
*  The implied lower bound of the row is computed as follows:
*
*     L' :=   sum   a[j] * l[j] +   sum   a[j] * u[j]                (9)
*           j in J+               j in J-
*
*  and as it follows from (3), (4), and (5):
*
*     L' := if j_min = 0 then f_min else -inf                       (10)
*
*  The implied upper bound of the row is computed as follows:
*
*     U' :=   sum   a[j] * u[j] +   sum   a[j] * l[j]               (11)
*           j in J+               j in J-
*
*  and as it follows from (6), (7), and (8):
*
*     U' := if j_max = 0 then f_max else +inf                       (12)
*
*  The implied bounds are stored in locations LL and UU. */

static void row_implied_bounds(const struct f_info *f, double *LL,
      double *UU)
{     *LL = (f->j_min == 0 ? f->f_min : -DBL_MAX);
      *UU = (f->j_max == 0 ? f->f_max : +DBL_MAX);
      return;
}

/***********************************************************************
*  col_implied_bounds - determine column implied bounds
*
*  Given a row (constraint)
*
*           n
*     L <= sum a[j] * x[j] <= U                                     (13)
*          j=1
*
*  and bounds of columns (variables)
*
*     l[j] <= x[j] <= u[j]
*
*  this routine determines implied bounds of variable x[k].
*
*  It is assumed that if L != -inf, the lower bound of the row can be
*  active, and if U != +inf, the upper bound of the row can be active.
*
*  ALGORITHM
*
*  From (13) it follows that
*
*     L <= sum a[j] * x[j] + a[k] * x[k] <= U
*          j!=k
*  or
*
*     L - sum a[j] * x[j] <= a[k] * x[k] <= U - sum a[j] * x[j]
*         j!=k                                  j!=k
*
*  Thus, if the row lower bound L can be active, implied lower bound of
*  term a[k] * x[k] can be determined as follows:
*
*     ilb(a[k] * x[k]) = min(L - sum a[j] * x[j]) =
*                                j!=k
*                                                                   (14)
*                      = L - max sum a[j] * x[j]
*                            j!=k
*
*  where, as it follows from (6), (7), and (8)
*
*                           / f_max - a[k] * u[k], j_max = 0, a[k] > 0
*                           |
*                           | f_max - a[k] * l[k], j_max = 0, a[k] < 0
*     max sum a[j] * x[j] = {
*         j!=k              | f_max,               j_max = k
*                           |
*                           \ +inf,                j_max != 0
*
*  and if the upper bound U can be active, implied upper bound of term
*  a[k] * x[k] can be determined as follows:
*
*     iub(a[k] * x[k]) = max(U - sum a[j] * x[j]) =
*                                j!=k
*                                                                   (15)
*                      = U - min sum a[j] * x[j]
*                            j!=k
*
*  where, as it follows from (3), (4), and (5)
*
*                           / f_min - a[k] * l[k], j_min = 0, a[k] > 0
*                           |
*                           | f_min - a[k] * u[k], j_min = 0, a[k] < 0
*     min sum a[j] * x[j] = {
*         j!=k              | f_min,               j_min = k
*                           |
*                           \ -inf,                j_min != 0
*
*  Since
*
*     ilb(a[k] * x[k]) <= a[k] * x[k] <= iub(a[k] * x[k])
*
*  implied lower and upper bounds of x[k] are determined as follows:
*
*     l'[k] := if a[k] > 0 then ilb / a[k] else ulb / a[k]          (16)
*
*     u'[k] := if a[k] > 0 then ulb / a[k] else ilb / a[k]          (17)
*
*  The implied bounds are stored in locations ll and uu. */

static void col_implied_bounds(const struct f_info *f, int n,
      const double a[], double L, double U, const double l[],
      const double u[], int k, double *ll, double *uu)
{     double ilb, iub;
      xassert(n >= 0);
      xassert(1 <= k && k <= n);
      /* determine implied lower bound of term a[k] * x[k] (14) */
      if (L == -DBL_MAX || f->f_max == +DBL_MAX)
         ilb = -DBL_MAX;
      else if (f->j_max == 0)
      {  if (a[k] > 0.0)
         {  xassert(u[k] != +DBL_MAX);
            ilb = L - (f->f_max - a[k] * u[k]);
         }
         else if (a[k] < 0.0)
         {  xassert(l[k] != -DBL_MAX);
            ilb = L - (f->f_max - a[k] * l[k]);
         }
         else
            xassert(a != a);
      }
      else if (f->j_max == k)
         ilb = L - f->f_max;
      else
         ilb = -DBL_MAX;
      /* determine implied upper bound of term a[k] * x[k] (15) */
      if (U == +DBL_MAX || f->f_min == -DBL_MAX)
         iub = +DBL_MAX;
      else if (f->j_min == 0)
      {  if (a[k] > 0.0)
         {  xassert(l[k] != -DBL_MAX);
            iub = U - (f->f_min - a[k] * l[k]);
         }
         else if (a[k] < 0.0)
         {  xassert(u[k] != +DBL_MAX);
            iub = U - (f->f_min - a[k] * u[k]);
         }
         else
            xassert(a != a);
      }
      else if (f->j_min == k)
         iub = U - f->f_min;
      else
         iub = +DBL_MAX;
      /* determine implied bounds of x[k] (16) and (17) */
#if 1
      /* do not use a[k] if it has small magnitude to prevent wrong
         implied bounds; for example, 1e-15 * x1 >= x2 + x3, where
         x1 >= -10, x2, x3 >= 0, would lead to wrong conclusion that
         x1 >= 0 */
      if (fabs(a[k]) < 1e-6)
         *ll = -DBL_MAX, *uu = +DBL_MAX; else
#endif
      if (a[k] > 0.0)
      {  *ll = (ilb == -DBL_MAX ? -DBL_MAX : ilb / a[k]);
         *uu = (iub == +DBL_MAX ? +DBL_MAX : iub / a[k]);
      }
      else if (a[k] < 0.0)
      {  *ll = (iub == +DBL_MAX ? -DBL_MAX : iub / a[k]);
         *uu = (ilb == -DBL_MAX ? +DBL_MAX : ilb / a[k]);
      }
      else
         xassert(a != a);
      return;
}

/***********************************************************************
*  check_row_bounds - check and relax original row bounds
*
*  Given a row (constraint)
*
*           n
*     L <= sum a[j] * x[j] <= U
*          j=1
*
*  and bounds of columns (variables)
*
*     l[j] <= x[j] <= u[j]
*
*  this routine checks the original row bounds L and U for feasibility
*  and redundancy. If the original lower bound L or/and upper bound U
*  cannot be active due to bounds of variables, the routine remove them
*  replacing by -inf or/and +inf, respectively.
*
*  If no primal infeasibility is detected, the routine returns zero,
*  otherwise non-zero. */

static int check_row_bounds(const struct f_info *f, double *L_,
      double *U_)
{     int ret = 0;
      double L = *L_, U = *U_, LL, UU;
      /* determine implied bounds of the row */
      row_implied_bounds(f, &LL, &UU);
      /* check if the original lower bound is infeasible */
      if (L != -DBL_MAX)
      {  double eps = 1e-3 * (1.0 + fabs(L));
         if (UU < L - eps)
         {  ret = 1;
            goto done;
         }
      }
      /* check if the original upper bound is infeasible */
      if (U != +DBL_MAX)
      {  double eps = 1e-3 * (1.0 + fabs(U));
         if (LL > U + eps)
         {  ret = 1;
            goto done;
         }
      }
      /* check if the original lower bound is redundant */
      if (L != -DBL_MAX)
      {  double eps = 1e-12 * (1.0 + fabs(L));
         if (LL > L - eps)
         {  /* it cannot be active, so remove it */
            *L_ = -DBL_MAX;
         }
      }
      /* check if the original upper bound is redundant */
      if (U != +DBL_MAX)
      {  double eps = 1e-12 * (1.0 + fabs(U));
         if (UU < U + eps)
         {  /* it cannot be active, so remove it */
            *U_ = +DBL_MAX;
         }
      }
done: return ret;
}

/***********************************************************************
*  check_col_bounds - check and tighten original column bounds
*
*  Given a row (constraint)
*
*           n
*     L <= sum a[j] * x[j] <= U
*          j=1
*
*  and bounds of columns (variables)
*
*     l[j] <= x[j] <= u[j]
*
*  for column (variable) x[j] this routine checks the original column
*  bounds l[j] and u[j] for feasibility and redundancy. If the original
*  lower bound l[j] or/and upper bound u[j] cannot be active due to
*  bounds of the constraint and other variables, the routine tighten
*  them replacing by corresponding implied bounds, if possible.
*
*  NOTE: It is assumed that if L != -inf, the row lower bound can be
*        active, and if U != +inf, the row upper bound can be active.
*
*  The flag means that variable x[j] is required to be integer.
*
*  New actual bounds for x[j] are stored in locations lj and uj.
*
*  If no primal infeasibility is detected, the routine returns zero,
*  otherwise non-zero. */

static int check_col_bounds(const struct f_info *f, int n,
      const double a[], double L, double U, const double l[],
      const double u[], int flag, int j, double *_lj, double *_uj)
{     int ret = 0;
      double lj, uj, ll, uu;
      xassert(n >= 0);
      xassert(1 <= j && j <= n);
      lj = l[j], uj = u[j];
      /* determine implied bounds of the column */
      col_implied_bounds(f, n, a, L, U, l, u, j, &ll, &uu);
      /* if x[j] is integral, round its implied bounds */
      if (flag)
      {  if (ll != -DBL_MAX)
            ll = (ll - floor(ll) < 1e-3 ? floor(ll) : ceil(ll));
         if (uu != +DBL_MAX)
            uu = (ceil(uu) - uu < 1e-3 ? ceil(uu) : floor(uu));
      }
      /* check if the original lower bound is infeasible */
      if (lj != -DBL_MAX)
      {  double eps = 1e-3 * (1.0 + fabs(lj));
         if (uu < lj - eps)
         {  ret = 1;
            goto done;
         }
      }
      /* check if the original upper bound is infeasible */
      if (uj != +DBL_MAX)
      {  double eps = 1e-3 * (1.0 + fabs(uj));
         if (ll > uj + eps)
         {  ret = 1;
            goto done;
         }
      }
      /* check if the original lower bound is redundant */
      if (ll != -DBL_MAX)
      {  double eps = 1e-3 * (1.0 + fabs(ll));
         if (lj < ll - eps)
         {  /* it cannot be active, so tighten it */
            lj = ll;
         }
      }
      /* check if the original upper bound is redundant */
      if (uu != +DBL_MAX)
      {  double eps = 1e-3 * (1.0 + fabs(uu));
         if (uj > uu + eps)
         {  /* it cannot be active, so tighten it */
            uj = uu;
         }
      }
      /* due to round-off errors it may happen that lj > uj (although
         lj < uj + eps, since no primal infeasibility is detected), so
         adjuct the new actual bounds to provide lj <= uj */
      if (!(lj == -DBL_MAX || uj == +DBL_MAX))
      {  double t1 = fabs(lj), t2 = fabs(uj);
         double eps = 1e-10 * (1.0 + (t1 <= t2 ? t1 : t2));
         if (lj > uj - eps)
         {  if (lj == l[j])
               uj = lj;
            else if (uj == u[j])
               lj = uj;
            else if (t1 <= t2)
               uj = lj;
            else
               lj = uj;
         }
      }
      *_lj = lj, *_uj = uj;
done: return ret;
}

/***********************************************************************
*  check_efficiency - check if change in column bounds is efficient
*
*  Given the original bounds of a column l and u and its new actual
*  bounds l' and u' (possibly tighten by the routine check_col_bounds)
*  this routine checks if the change in the column bounds is efficient
*  enough. If so, the routine returns non-zero, otherwise zero.
*
*  The flag means that the variable is required to be integer. */

static int check_efficiency(int flag, double l, double u, double ll,
      double uu)
{     int eff = 0;
      /* check efficiency for lower bound */
      if (l < ll)
      {  if (flag || l == -DBL_MAX)
            eff++;
         else
         {  double r;
            if (u == +DBL_MAX)
               r = 1.0 + fabs(l);
            else
               r = 1.0 + (u - l);
            if (ll - l >= 0.25 * r)
               eff++;
         }
      }
      /* check efficiency for upper bound */
      if (u > uu)
      {  if (flag || u == +DBL_MAX)
            eff++;
         else
         {  double r;
            if (l == -DBL_MAX)
               r = 1.0 + fabs(u);
            else
               r = 1.0 + (u - l);
            if (u - uu >= 0.25 * r)
               eff++;
         }
      }
      return eff;
}

/***********************************************************************
*  basic_preprocessing - perform basic preprocessing
*
*  This routine performs basic preprocessing of the specified MIP that
*  includes relaxing some row bounds and tightening some column bounds.
*
*  On entry the arrays L and U contains original row bounds, and the
*  arrays l and u contains original column bounds:
*
*  L[0] is the lower bound of the objective row;
*  L[i], i = 1,...,m, is the lower bound of i-th row;
*  U[0] is the upper bound of the objective row;
*  U[i], i = 1,...,m, is the upper bound of i-th row;
*  l[0] is not used;
*  l[j], j = 1,...,n, is the lower bound of j-th column;
*  u[0] is not used;
*  u[j], j = 1,...,n, is the upper bound of j-th column.
*
*  On exit the arrays L, U, l, and u contain new actual bounds of rows
*  and column in the same locations.
*
*  The parameters nrs and num specify an initial list of rows to be
*  processed:
*
*  nrs is the number of rows in the initial list, 0 <= nrs <= m+1;
*  num[0] is not used;
*  num[1,...,nrs] are row numbers (0 means the objective row).
*
*  The parameter max_pass specifies the maximal number of times that
*  each row can be processed, max_pass > 0.
*
*  If no primal infeasibility is detected, the routine returns zero,
*  otherwise non-zero. */

static int basic_preprocessing(glp_prob *mip, double L[], double U[],
      double l[], double u[], int nrs, const int num[], int max_pass)
{     int m = mip->m;
      int n = mip->n;
      struct f_info f;
      int i, j, k, len, size, ret = 0;
      int *ind, *list, *mark, *pass;
      double *val, *lb, *ub;
      xassert(0 <= nrs && nrs <= m+1);
      xassert(max_pass > 0);
      /* allocate working arrays */
      ind = xcalloc(1+n, sizeof(int));
      list = xcalloc(1+m+1, sizeof(int));
      mark = xcalloc(1+m+1, sizeof(int));
      memset(&mark[0], 0, (m+1) * sizeof(int));
      pass = xcalloc(1+m+1, sizeof(int));
      memset(&pass[0], 0, (m+1) * sizeof(int));
      val = xcalloc(1+n, sizeof(double));
      lb = xcalloc(1+n, sizeof(double));
      ub = xcalloc(1+n, sizeof(double));
      /* initialize the list of rows to be processed */
      size = 0;
      for (k = 1; k <= nrs; k++)
      {  i = num[k];
         xassert(0 <= i && i <= m);
         /* duplicate row numbers are not allowed */
         xassert(!mark[i]);
         list[++size] = i, mark[i] = 1;
      }
      xassert(size == nrs);
      /* process rows in the list until it becomes empty */
      while (size > 0)
      {  /* get a next row from the list */
         i = list[size--], mark[i] = 0;
         /* increase the row processing count */
         pass[i]++;
         /* if the row is free, skip it */
         if (L[i] == -DBL_MAX && U[i] == +DBL_MAX) continue;
         /* obtain coefficients of the row */
         len = 0;
         if (i == 0)
         {  for (j = 1; j <= n; j++)
            {  GLPCOL *col = mip->col[j];
               if (col->coef != 0.0)
                  len++, ind[len] = j, val[len] = col->coef;
            }
         }
         else
         {  GLPROW *row = mip->row[i];
            GLPAIJ *aij;
            for (aij = row->ptr; aij != NULL; aij = aij->r_next)
               len++, ind[len] = aij->col->j, val[len] = aij->val;
         }
         /* determine lower and upper bounds of columns corresponding
            to non-zero row coefficients */
         for (k = 1; k <= len; k++)
            j = ind[k], lb[k] = l[j], ub[k] = u[j];
         /* prepare the row info to determine implied bounds */
         prepare_row_info(len, val, lb, ub, &f);
         /* check and relax bounds of the row */
         if (check_row_bounds(&f, &L[i], &U[i]))
         {  /* the feasible region is empty */
            ret = 1;
            goto done;
         }
         /* if the row became free, drop it */
         if (L[i] == -DBL_MAX && U[i] == +DBL_MAX) continue;
         /* process columns having non-zero coefficients in the row */
         for (k = 1; k <= len; k++)
         {  GLPCOL *col;
            int flag, eff;
            double ll, uu;
            /* take a next column in the row */
            j = ind[k], col = mip->col[j];
            flag = col->kind != GLP_CV;
            /* check and tighten bounds of the column */
            if (check_col_bounds(&f, len, val, L[i], U[i], lb, ub,
                flag, k, &ll, &uu))
            {  /* the feasible region is empty */
               ret = 1;
               goto done;
            }
            /* check if change in the column bounds is efficient */
            eff = check_efficiency(flag, l[j], u[j], ll, uu);
            /* set new actual bounds of the column */
            l[j] = ll, u[j] = uu;
            /* if the change is efficient, add all rows affected by the
               corresponding column, to the list */
            if (eff > 0)
            {  GLPAIJ *aij;
               for (aij = col->ptr; aij != NULL; aij = aij->c_next)
               {  int ii = aij->row->i;
                  /* if the row was processed maximal number of times,
                     skip it */
                  if (pass[ii] >= max_pass) continue;
                  /* if the row is free, skip it */
                  if (L[ii] == -DBL_MAX && U[ii] == +DBL_MAX) continue;
                  /* put the row into the list */
                  if (mark[ii] == 0)
                  {  xassert(size <= m);
                     list[++size] = ii, mark[ii] = 1;
                  }
               }
            }
         }
      }
done: /* free working arrays */
      xfree(ind);
      xfree(list);
      xfree(mark);
      xfree(pass);
      xfree(val);
      xfree(lb);
      xfree(ub);
      return ret;
}

/***********************************************************************
*  NAME
*
*  ios_preprocess_node - preprocess current subproblem
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  int ios_preprocess_node(glp_tree *tree, int max_pass);
*
*  DESCRIPTION
*
*  The routine ios_preprocess_node performs basic preprocessing of the
*  current subproblem.
*
*  RETURNS
*
*  If no primal infeasibility is detected, the routine returns zero,
*  otherwise non-zero. */

int ios_preprocess_node(glp_tree *tree, int max_pass)
{     glp_prob *mip = tree->mip;
      int m = mip->m;
      int n = mip->n;
      int i, j, nrs, *num, ret = 0;
      double *L, *U, *l, *u;
      /* the current subproblem must exist */
      xassert(tree->curr != NULL);
      /* determine original row bounds */
      L = xcalloc(1+m, sizeof(double));
      U = xcalloc(1+m, sizeof(double));
      switch (mip->mip_stat)
      {  case GLP_UNDEF:
            L[0] = -DBL_MAX, U[0] = +DBL_MAX;
            break;
         case GLP_FEAS:
            switch (mip->dir)
            {  case GLP_MIN:
                  L[0] = -DBL_MAX, U[0] = mip->mip_obj - mip->c0;
                  break;
               case GLP_MAX:
                  L[0] = mip->mip_obj - mip->c0, U[0] = +DBL_MAX;
                  break;
               default:
                  xassert(mip != mip);
            }
            break;
         default:
            xassert(mip != mip);
      }
      for (i = 1; i <= m; i++)
      {  L[i] = glp_get_row_lb(mip, i);
         U[i] = glp_get_row_ub(mip, i);
      }
      /* determine original column bounds */
      l = xcalloc(1+n, sizeof(double));
      u = xcalloc(1+n, sizeof(double));
      for (j = 1; j <= n; j++)
      {  l[j] = glp_get_col_lb(mip, j);
         u[j] = glp_get_col_ub(mip, j);
      }
      /* build the initial list of rows to be analyzed */
      nrs = m + 1;
      num = xcalloc(1+nrs, sizeof(int));
      for (i = 1; i <= nrs; i++) num[i] = i - 1;
      /* perform basic preprocessing */
      if (basic_preprocessing(mip , L, U, l, u, nrs, num, max_pass))
      {  ret = 1;
         goto done;
      }
      /* set new actual (relaxed) row bounds */
      for (i = 1; i <= m; i++)
      {  /* consider only non-active rows to keep dual feasibility */
         if (glp_get_row_stat(mip, i) == GLP_BS)
         {  if (L[i] == -DBL_MAX && U[i] == +DBL_MAX)
               glp_set_row_bnds(mip, i, GLP_FR, 0.0, 0.0);
            else if (U[i] == +DBL_MAX)
               glp_set_row_bnds(mip, i, GLP_LO, L[i], 0.0);
            else if (L[i] == -DBL_MAX)
               glp_set_row_bnds(mip, i, GLP_UP, 0.0, U[i]);
         }
      }
      /* set new actual (tightened) column bounds */
      for (j = 1; j <= n; j++)
      {  int type;
         if (l[j] == -DBL_MAX && u[j] == +DBL_MAX)
            type = GLP_FR;
         else if (u[j] == +DBL_MAX)
            type = GLP_LO;
         else if (l[j] == -DBL_MAX)
            type = GLP_UP;
         else if (l[j] != u[j])
            type = GLP_DB;
         else
            type = GLP_FX;
         glp_set_col_bnds(mip, j, type, l[j], u[j]);
      }
done: /* free working arrays and return */
      xfree(L);
      xfree(U);
      xfree(l);
      xfree(u);
      xfree(num);
      return ret;
}

/* eof */
