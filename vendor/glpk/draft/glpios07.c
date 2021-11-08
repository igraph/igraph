/* glpios07.c (mixed cover cut generator) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2005-2018 Free Software Foundation, Inc.
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
#include "ios.h"

/*----------------------------------------------------------------------
-- COVER INEQUALITIES
--
-- Consider the set of feasible solutions to 0-1 knapsack problem:
--
--    sum a[j]*x[j] <= b,                                            (1)
--  j in J
--
--    x[j] is binary,                                                (2)
--
-- where, wlog, we assume that a[j] > 0 (since 0-1 variables can be
-- complemented) and a[j] <= b (since a[j] > b implies x[j] = 0).
--
-- A set C within J is called a cover if
--
--    sum a[j] > b.                                                  (3)
--  j in C
--
-- For any cover C the inequality
--
--    sum x[j] <= |C| - 1                                            (4)
--  j in C
--
-- is called a cover inequality and is valid for (1)-(2).
--
-- MIXED COVER INEQUALITIES
--
-- Consider the set of feasible solutions to mixed knapsack problem:
--
--    sum a[j]*x[j] + y <= b,                                        (5)
--  j in J
--
--    x[j] is binary,                                                (6)
--
--    0 <= y <= u is continuous,                                     (7)
--
-- where again we assume that a[j] > 0.
--
-- Let C within J be some set. From (1)-(4) it follows that
--
--    sum a[j] > b - y                                               (8)
--  j in C
--
-- implies
--
--    sum x[j] <= |C| - 1.                                           (9)
--  j in C
--
-- Thus, we need to modify the inequality (9) in such a way that it be
-- a constraint only if the condition (8) is satisfied.
--
-- Consider the following inequality:
--
--    sum x[j] <= |C| - t.                                          (10)
--  j in C
--
-- If 0 < t <= 1, then (10) is equivalent to (9), because all x[j] are
-- binary variables. On the other hand, if t <= 0, (10) being satisfied
-- for any values of x[j] is not a constraint.
--
-- Let
--
--    t' = sum a[j] + y - b.                                        (11)
--       j in C
--
-- It is understood that the condition t' > 0 is equivalent to (8).
-- Besides, from (6)-(7) it follows that t' has an implied upper bound:
--
--    t'max = sum a[j] + u - b.                                     (12)
--          j in C
--
-- This allows to express the parameter t having desired properties:
--
--    t = t' / t'max.                                               (13)
--
-- In fact, t <= 1 by definition, and t > 0 being equivalent to t' > 0
-- is equivalent to (8).
--
-- Thus, the inequality (10), where t is given by formula (13) is valid
-- for (5)-(7).
--
-- Note that if u = 0, then y = 0, so t = 1, and the conditions (8) and
-- (10) is transformed to the conditions (3) and (4).
--
-- GENERATING MIXED COVER CUTS
--
-- To generate a mixed cover cut in the form (10) we need to find such
-- set C which satisfies to the inequality (8) and for which, in turn,
-- the inequality (10) is violated in the current point.
--
-- Substituting t from (13) to (10) gives:
--
--                        1
--    sum x[j] <= |C| - -----  (sum a[j] + y - b),                  (14)
--  j in C              t'max j in C
--
-- and finally we have the cut inequality in the standard form:
--
--    sum x[j] + alfa * y <= beta,                                  (15)
--  j in C
--
-- where:
--
--    alfa = 1 / t'max,                                             (16)
--
--    beta = |C| - alfa *  (sum a[j] - b).                          (17)
--                        j in C                                      */

#if 1
#define MAXTRY 1000
#else
#define MAXTRY 10000
#endif

static int cover2(int n, double a[], double b, double u, double x[],
      double y, int cov[], double *_alfa, double *_beta)
{     /* try to generate mixed cover cut using two-element cover */
      int i, j, try = 0, ret = 0;
      double eps, alfa, beta, temp, rmax = 0.001;
      eps = 0.001 * (1.0 + fabs(b));
      for (i = 0+1; i <= n; i++)
      for (j = i+1; j <= n; j++)
      {  /* C = {i, j} */
         try++;
         if (try > MAXTRY) goto done;
         /* check if condition (8) is satisfied */
         if (a[i] + a[j] + y > b + eps)
         {  /* compute parameters for inequality (15) */
            temp = a[i] + a[j] - b;
            alfa = 1.0 / (temp + u);
            beta = 2.0 - alfa * temp;
            /* compute violation of inequality (15) */
            temp = x[i] + x[j] + alfa * y - beta;
            /* choose C providing maximum violation */
            if (rmax < temp)
            {  rmax = temp;
               cov[1] = i;
               cov[2] = j;
               *_alfa = alfa;
               *_beta = beta;
               ret = 1;
            }
         }
      }
done: return ret;
}

static int cover3(int n, double a[], double b, double u, double x[],
      double y, int cov[], double *_alfa, double *_beta)
{     /* try to generate mixed cover cut using three-element cover */
      int i, j, k, try = 0, ret = 0;
      double eps, alfa, beta, temp, rmax = 0.001;
      eps = 0.001 * (1.0 + fabs(b));
      for (i = 0+1; i <= n; i++)
      for (j = i+1; j <= n; j++)
      for (k = j+1; k <= n; k++)
      {  /* C = {i, j, k} */
         try++;
         if (try > MAXTRY) goto done;
         /* check if condition (8) is satisfied */
         if (a[i] + a[j] + a[k] + y > b + eps)
         {  /* compute parameters for inequality (15) */
            temp = a[i] + a[j] + a[k] - b;
            alfa = 1.0 / (temp + u);
            beta = 3.0 - alfa * temp;
            /* compute violation of inequality (15) */
            temp = x[i] + x[j] + x[k] + alfa * y - beta;
            /* choose C providing maximum violation */
            if (rmax < temp)
            {  rmax = temp;
               cov[1] = i;
               cov[2] = j;
               cov[3] = k;
               *_alfa = alfa;
               *_beta = beta;
               ret = 1;
            }
         }
      }
done: return ret;
}

static int cover4(int n, double a[], double b, double u, double x[],
      double y, int cov[], double *_alfa, double *_beta)
{     /* try to generate mixed cover cut using four-element cover */
      int i, j, k, l, try = 0, ret = 0;
      double eps, alfa, beta, temp, rmax = 0.001;
      eps = 0.001 * (1.0 + fabs(b));
      for (i = 0+1; i <= n; i++)
      for (j = i+1; j <= n; j++)
      for (k = j+1; k <= n; k++)
      for (l = k+1; l <= n; l++)
      {  /* C = {i, j, k, l} */
         try++;
         if (try > MAXTRY) goto done;
         /* check if condition (8) is satisfied */
         if (a[i] + a[j] + a[k] + a[l] + y > b + eps)
         {  /* compute parameters for inequality (15) */
            temp = a[i] + a[j] + a[k] + a[l] - b;
            alfa = 1.0 / (temp + u);
            beta = 4.0 - alfa * temp;
            /* compute violation of inequality (15) */
            temp = x[i] + x[j] + x[k] + x[l] + alfa * y - beta;
            /* choose C providing maximum violation */
            if (rmax < temp)
            {  rmax = temp;
               cov[1] = i;
               cov[2] = j;
               cov[3] = k;
               cov[4] = l;
               *_alfa = alfa;
               *_beta = beta;
               ret = 1;
            }
         }
      }
done: return ret;
}

static int cover(int n, double a[], double b, double u, double x[],
      double y, int cov[], double *alfa, double *beta)
{     /* try to generate mixed cover cut;
         input (see (5)):
         n        is the number of binary variables;
         a[1:n]   are coefficients at binary variables;
         b        is the right-hand side;
         u        is upper bound of continuous variable;
         x[1:n]   are values of binary variables at current point;
         y        is value of continuous variable at current point;
         output (see (15), (16), (17)):
         cov[1:r] are indices of binary variables included in cover C,
                  where r is the set cardinality returned on exit;
         alfa     coefficient at continuous variable;
         beta     is the right-hand side; */
      int j;
      /* perform some sanity checks */
      xassert(n >= 2);
      for (j = 1; j <= n; j++) xassert(a[j] > 0.0);
#if 1 /* ??? */
      xassert(b > -1e-5);
#else
      xassert(b > 0.0);
#endif
      xassert(u >= 0.0);
      for (j = 1; j <= n; j++) xassert(0.0 <= x[j] && x[j] <= 1.0);
      xassert(0.0 <= y && y <= u);
      /* try to generate mixed cover cut */
      if (cover2(n, a, b, u, x, y, cov, alfa, beta)) return 2;
      if (cover3(n, a, b, u, x, y, cov, alfa, beta)) return 3;
      if (cover4(n, a, b, u, x, y, cov, alfa, beta)) return 4;
      return 0;
}

/*----------------------------------------------------------------------
-- lpx_cover_cut - generate mixed cover cut.
--
-- SYNOPSIS
--
-- int lpx_cover_cut(LPX *lp, int len, int ind[], double val[],
--    double work[]);
--
-- DESCRIPTION
--
-- The routine lpx_cover_cut generates a mixed cover cut for a given
-- row of the MIP problem.
--
-- The given row of the MIP problem should be explicitly specified in
-- the form:
--
--    sum{j in J} a[j]*x[j] <= b.                                    (1)
--
-- On entry indices (ordinal numbers) of structural variables, which
-- have non-zero constraint coefficients, should be placed in locations
-- ind[1], ..., ind[len], and corresponding constraint coefficients
-- should be placed in locations val[1], ..., val[len]. The right-hand
-- side b should be stored in location val[0].
--
-- The working array work should have at least nb locations, where nb
-- is the number of binary variables in (1).
--
-- The routine generates a mixed cover cut in the same form as (1) and
-- stores the cut coefficients and right-hand side in the same way as
-- just described above.
--
-- RETURNS
--
-- If the cutting plane has been successfully generated, the routine
-- returns 1 <= len' <= n, which is the number of non-zero coefficients
-- in the inequality constraint. Otherwise, the routine returns zero. */

static int lpx_cover_cut(glp_prob *lp, int len, int ind[],
      double val[], double work[])
{     int cov[1+4], j, k, nb, newlen, r;
      double f_min, f_max, alfa, beta, u, *x = work, y;
      /* substitute and remove fixed variables */
      newlen = 0;
      for (k = 1; k <= len; k++)
      {  j = ind[k];
         if (glp_get_col_type(lp, j) == GLP_FX)
            val[0] -= val[k] * glp_get_col_lb(lp, j);
         else
         {  newlen++;
            ind[newlen] = ind[k];
            val[newlen] = val[k];
         }
      }
      len = newlen;
      /* move binary variables to the beginning of the list so that
         elements 1, 2, ..., nb correspond to binary variables, and
         elements nb+1, nb+2, ..., len correspond to rest variables */
      nb = 0;
      for (k = 1; k <= len; k++)
      {  j = ind[k];
         if (glp_get_col_kind(lp, j) == GLP_BV)
         {  /* binary variable */
            int ind_k;
            double val_k;
            nb++;
            ind_k = ind[nb], val_k = val[nb];
            ind[nb] = ind[k], val[nb] = val[k];
            ind[k] = ind_k, val[k] = val_k;
         }
      }
      /* now the specified row has the form:
         sum a[j]*x[j] + sum a[j]*y[j] <= b,
         where x[j] are binary variables, y[j] are rest variables */
      /* at least two binary variables are needed */
      if (nb < 2) return 0;
      /* compute implied lower and upper bounds for sum a[j]*y[j] */
      f_min = f_max = 0.0;
      for (k = nb+1; k <= len; k++)
      {  j = ind[k];
         /* both bounds must be finite */
         if (glp_get_col_type(lp, j) != GLP_DB) return 0;
         if (val[k] > 0.0)
         {  f_min += val[k] * glp_get_col_lb(lp, j);
            f_max += val[k] * glp_get_col_ub(lp, j);
         }
         else
         {  f_min += val[k] * glp_get_col_ub(lp, j);
            f_max += val[k] * glp_get_col_lb(lp, j);
         }
      }
      /* sum a[j]*x[j] + sum a[j]*y[j] <= b ===>
         sum a[j]*x[j] + (sum a[j]*y[j] - f_min) <= b - f_min ===>
         sum a[j]*x[j] + y <= b - f_min,
         where y = sum a[j]*y[j] - f_min;
         note that 0 <= y <= u, u = f_max - f_min */
      /* determine upper bound of y */
      u = f_max - f_min;
      /* determine value of y at the current point */
      y = 0.0;
      for (k = nb+1; k <= len; k++)
      {  j = ind[k];
         y += val[k] * glp_get_col_prim(lp, j);
      }
      y -= f_min;
      if (y < 0.0) y = 0.0;
      if (y > u) y = u;
      /* modify the right-hand side b */
      val[0] -= f_min;
      /* now the transformed row has the form:
         sum a[j]*x[j] + y <= b, where 0 <= y <= u */
      /* determine values of x[j] at the current point */
      for (k = 1; k <= nb; k++)
      {  j = ind[k];
         x[k] = glp_get_col_prim(lp, j);
         if (x[k] < 0.0) x[k] = 0.0;
         if (x[k] > 1.0) x[k] = 1.0;
      }
      /* if a[j] < 0, replace x[j] by its complement 1 - x'[j] */
      for (k = 1; k <= nb; k++)
      {  if (val[k] < 0.0)
         {  ind[k] = - ind[k];
            val[k] = - val[k];
            val[0] += val[k];
            x[k] = 1.0 - x[k];
         }
      }
      /* try to generate a mixed cover cut for the transformed row */
      r = cover(nb, val, val[0], u, x, y, cov, &alfa, &beta);
      if (r == 0) return 0;
      xassert(2 <= r && r <= 4);
      /* now the cut is in the form:
         sum{j in C} x[j] + alfa * y <= beta */
      /* store the right-hand side beta */
      ind[0] = 0, val[0] = beta;
      /* restore the original ordinal numbers of x[j] */
      for (j = 1; j <= r; j++) cov[j] = ind[cov[j]];
      /* store cut coefficients at binary variables complementing back
         the variables having negative row coefficients */
      xassert(r <= nb);
      for (k = 1; k <= r; k++)
      {  if (cov[k] > 0)
         {  ind[k] = +cov[k];
            val[k] = +1.0;
         }
         else
         {  ind[k] = -cov[k];
            val[k] = -1.0;
            val[0] -= 1.0;
         }
      }
      /* substitute y = sum a[j]*y[j] - f_min */
      for (k = nb+1; k <= len; k++)
      {  r++;
         ind[r] = ind[k];
         val[r] = alfa * val[k];
      }
      val[0] += alfa * f_min;
      xassert(r <= len);
      len = r;
      return len;
}

/*----------------------------------------------------------------------
-- lpx_eval_row - compute explictily specified row.
--
-- SYNOPSIS
--
-- double lpx_eval_row(LPX *lp, int len, int ind[], double val[]);
--
-- DESCRIPTION
--
-- The routine lpx_eval_row computes the primal value of an explicitly
-- specified row using current values of structural variables.
--
-- The explicitly specified row may be thought as a linear form:
--
--    y = a[1]*x[m+1] + a[2]*x[m+2] + ... + a[n]*x[m+n],
--
-- where y is an auxiliary variable for this row, a[j] are coefficients
-- of the linear form, x[m+j] are structural variables.
--
-- On entry column indices and numerical values of non-zero elements of
-- the row should be stored in locations ind[1], ..., ind[len] and
-- val[1], ..., val[len], where len is the number of non-zero elements.
-- The array ind and val are not changed on exit.
--
-- RETURNS
--
-- The routine returns a computed value of y, the auxiliary variable of
-- the specified row. */

static double lpx_eval_row(glp_prob *lp, int len, int ind[],
      double val[])
{     int n = glp_get_num_cols(lp);
      int j, k;
      double sum = 0.0;
      if (len < 0)
         xerror("lpx_eval_row: len = %d; invalid row length\n", len);
      for (k = 1; k <= len; k++)
      {  j = ind[k];
         if (!(1 <= j && j <= n))
            xerror("lpx_eval_row: j = %d; column number out of range\n",
               j);
         sum += val[k] * glp_get_col_prim(lp, j);
      }
      return sum;
}

/***********************************************************************
*  NAME
*
*  ios_cov_gen - generate mixed cover cuts
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_cov_gen(glp_tree *tree);
*
*  DESCRIPTION
*
*  The routine ios_cov_gen generates mixed cover cuts for the current
*  point and adds them to the cut pool. */

void ios_cov_gen(glp_tree *tree)
{     glp_prob *prob = tree->mip;
      int m = glp_get_num_rows(prob);
      int n = glp_get_num_cols(prob);
      int i, k, type, kase, len, *ind;
      double r, *val, *work;
      xassert(glp_get_status(prob) == GLP_OPT);
      /* allocate working arrays */
      ind = xcalloc(1+n, sizeof(int));
      val = xcalloc(1+n, sizeof(double));
      work = xcalloc(1+n, sizeof(double));
      /* look through all rows */
      for (i = 1; i <= m; i++)
      for (kase = 1; kase <= 2; kase++)
      {  type = glp_get_row_type(prob, i);
         if (kase == 1)
         {  /* consider rows of '<=' type */
            if (!(type == GLP_UP || type == GLP_DB)) continue;
            len = glp_get_mat_row(prob, i, ind, val);
            val[0] = glp_get_row_ub(prob, i);
         }
         else
         {  /* consider rows of '>=' type */
            if (!(type == GLP_LO || type == GLP_DB)) continue;
            len = glp_get_mat_row(prob, i, ind, val);
            for (k = 1; k <= len; k++) val[k] = - val[k];
            val[0] = - glp_get_row_lb(prob, i);
         }
         /* generate mixed cover cut:
            sum{j in J} a[j] * x[j] <= b */
         len = lpx_cover_cut(prob, len, ind, val, work);
         if (len == 0) continue;
         /* at the current point the cut inequality is violated, i.e.
            sum{j in J} a[j] * x[j] - b > 0 */
         r = lpx_eval_row(prob, len, ind, val) - val[0];
         if (r < 1e-3) continue;
         /* add the cut to the cut pool */
         glp_ios_add_row(tree, NULL, GLP_RF_COV, 0, len, ind, val,
            GLP_UP, val[0]);
      }
      /* free working arrays */
      xfree(ind);
      xfree(val);
      xfree(work);
      return;
}

/* eof */
