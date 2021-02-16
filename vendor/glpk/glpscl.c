/* glpscl.c */

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

#include "glpapi.h"

/***********************************************************************
*  min_row_aij - determine minimal |a[i,j]| in i-th row
*
*  This routine returns minimal magnitude of (non-zero) constraint
*  coefficients in i-th row of the constraint matrix.
*
*  If the parameter scaled is zero, the original constraint matrix A is
*  assumed. Otherwise, the scaled constraint matrix R*A*S is assumed.
*
*  If i-th row of the matrix is empty, the routine returns 1. */

static double min_row_aij(glp_prob *lp, int i, int scaled)
{     GLPAIJ *aij;
      double min_aij, temp;
      xassert(1 <= i && i <= lp->m);
      min_aij = 1.0;
      for (aij = lp->row[i]->ptr; aij != NULL; aij = aij->r_next)
      {  temp = fabs(aij->val);
         if (scaled) temp *= (aij->row->rii * aij->col->sjj);
         if (aij->r_prev == NULL || min_aij > temp)
            min_aij = temp;
      }
      return min_aij;
}

/***********************************************************************
*  max_row_aij - determine maximal |a[i,j]| in i-th row
*
*  This routine returns maximal magnitude of (non-zero) constraint
*  coefficients in i-th row of the constraint matrix.
*
*  If the parameter scaled is zero, the original constraint matrix A is
*  assumed. Otherwise, the scaled constraint matrix R*A*S is assumed.
*
*  If i-th row of the matrix is empty, the routine returns 1. */

static double max_row_aij(glp_prob *lp, int i, int scaled)
{     GLPAIJ *aij;
      double max_aij, temp;
      xassert(1 <= i && i <= lp->m);
      max_aij = 1.0;
      for (aij = lp->row[i]->ptr; aij != NULL; aij = aij->r_next)
      {  temp = fabs(aij->val);
         if (scaled) temp *= (aij->row->rii * aij->col->sjj);
         if (aij->r_prev == NULL || max_aij < temp)
            max_aij = temp;
      }
      return max_aij;
}

/***********************************************************************
*  min_col_aij - determine minimal |a[i,j]| in j-th column
*
*  This routine returns minimal magnitude of (non-zero) constraint
*  coefficients in j-th column of the constraint matrix.
*
*  If the parameter scaled is zero, the original constraint matrix A is
*  assumed. Otherwise, the scaled constraint matrix R*A*S is assumed.
*
*  If j-th column of the matrix is empty, the routine returns 1. */

static double min_col_aij(glp_prob *lp, int j, int scaled)
{     GLPAIJ *aij;
      double min_aij, temp;
      xassert(1 <= j && j <= lp->n);
      min_aij = 1.0;
      for (aij = lp->col[j]->ptr; aij != NULL; aij = aij->c_next)
      {  temp = fabs(aij->val);
         if (scaled) temp *= (aij->row->rii * aij->col->sjj);
         if (aij->c_prev == NULL || min_aij > temp)
            min_aij = temp;
      }
      return min_aij;
}

/***********************************************************************
*  max_col_aij - determine maximal |a[i,j]| in j-th column
*
*  This routine returns maximal magnitude of (non-zero) constraint
*  coefficients in j-th column of the constraint matrix.
*
*  If the parameter scaled is zero, the original constraint matrix A is
*  assumed. Otherwise, the scaled constraint matrix R*A*S is assumed.
*
*  If j-th column of the matrix is empty, the routine returns 1. */

static double max_col_aij(glp_prob *lp, int j, int scaled)
{     GLPAIJ *aij;
      double max_aij, temp;
      xassert(1 <= j && j <= lp->n);
      max_aij = 1.0;
      for (aij = lp->col[j]->ptr; aij != NULL; aij = aij->c_next)
      {  temp = fabs(aij->val);
         if (scaled) temp *= (aij->row->rii * aij->col->sjj);
         if (aij->c_prev == NULL || max_aij < temp)
            max_aij = temp;
      }
      return max_aij;
}

/***********************************************************************
*  min_mat_aij - determine minimal |a[i,j]| in constraint matrix
*
*  This routine returns minimal magnitude of (non-zero) constraint
*  coefficients in the constraint matrix.
*
*  If the parameter scaled is zero, the original constraint matrix A is
*  assumed. Otherwise, the scaled constraint matrix R*A*S is assumed.
*
*  If the matrix is empty, the routine returns 1. */

static double min_mat_aij(glp_prob *lp, int scaled)
{     int i;
      double min_aij, temp;
      min_aij = 1.0;
      for (i = 1; i <= lp->m; i++)
      {  temp = min_row_aij(lp, i, scaled);
         if (i == 1 || min_aij > temp)
            min_aij = temp;
      }
      return min_aij;
}

/***********************************************************************
*  max_mat_aij - determine maximal |a[i,j]| in constraint matrix
*
*  This routine returns maximal magnitude of (non-zero) constraint
*  coefficients in the constraint matrix.
*
*  If the parameter scaled is zero, the original constraint matrix A is
*  assumed. Otherwise, the scaled constraint matrix R*A*S is assumed.
*
*  If the matrix is empty, the routine returns 1. */

static double max_mat_aij(glp_prob *lp, int scaled)
{     int i;
      double max_aij, temp;
      max_aij = 1.0;
      for (i = 1; i <= lp->m; i++)
      {  temp = max_row_aij(lp, i, scaled);
         if (i == 1 || max_aij < temp)
            max_aij = temp;
      }
      return max_aij;
}

/***********************************************************************
*  eq_scaling - perform equilibration scaling
*
*  This routine performs equilibration scaling of rows and columns of
*  the constraint matrix.
*
*  If the parameter flag is zero, the routine scales rows at first and
*  then columns. Otherwise, the routine scales columns and then rows.
*
*  Rows are scaled as follows:
*
*                         n
*     a'[i,j] = a[i,j] / max |a[i,j]|,  i = 1,...,m.
*                        j=1
*
*  This makes the infinity (maximum) norm of each row of the matrix
*  equal to 1.
*
*  Columns are scaled as follows:
*
*                         n
*     a'[i,j] = a[i,j] / max |a[i,j]|,  j = 1,...,n.
*                        i=1
*
*  This makes the infinity (maximum) norm of each column of the matrix
*  equal to 1. */

static void eq_scaling(glp_prob *lp, int flag)
{     int i, j, pass;
      double temp;
      xassert(flag == 0 || flag == 1);
      for (pass = 0; pass <= 1; pass++)
      {  if (pass == flag)
         {  /* scale rows */
            for (i = 1; i <= lp->m; i++)
            {  temp = max_row_aij(lp, i, 1);
               glp_set_rii(lp, i, glp_get_rii(lp, i) / temp);
            }
         }
         else
         {  /* scale columns */
            for (j = 1; j <= lp->n; j++)
            {  temp = max_col_aij(lp, j, 1);
               glp_set_sjj(lp, j, glp_get_sjj(lp, j) / temp);
            }
         }
      }
      return;
}

/***********************************************************************
*  gm_scaling - perform geometric mean scaling
*
*  This routine performs geometric mean scaling of rows and columns of
*  the constraint matrix.
*
*  If the parameter flag is zero, the routine scales rows at first and
*  then columns. Otherwise, the routine scales columns and then rows.
*
*  Rows are scaled as follows:
*
*     a'[i,j] = a[i,j] / sqrt(alfa[i] * beta[i]),  i = 1,...,m,
*
*  where:
*                n                        n
*     alfa[i] = min |a[i,j]|,  beta[i] = max |a[i,j]|.
*               j=1                      j=1
*
*  This allows decreasing the ratio beta[i] / alfa[i] for each row of
*  the matrix.
*
*  Columns are scaled as follows:
*
*     a'[i,j] = a[i,j] / sqrt(alfa[j] * beta[j]),  j = 1,...,n,
*
*  where:
*                m                        m
*     alfa[j] = min |a[i,j]|,  beta[j] = max |a[i,j]|.
*               i=1                      i=1
*
*  This allows decreasing the ratio beta[j] / alfa[j] for each column
*  of the matrix. */

static void gm_scaling(glp_prob *lp, int flag)
{     int i, j, pass;
      double temp;
      xassert(flag == 0 || flag == 1);
      for (pass = 0; pass <= 1; pass++)
      {  if (pass == flag)
         {  /* scale rows */
            for (i = 1; i <= lp->m; i++)
            {  temp = min_row_aij(lp, i, 1) * max_row_aij(lp, i, 1);
               glp_set_rii(lp, i, glp_get_rii(lp, i) / sqrt(temp));
            }
         }
         else
         {  /* scale columns */
            for (j = 1; j <= lp->n; j++)
            {  temp = min_col_aij(lp, j, 1) * max_col_aij(lp, j, 1);
               glp_set_sjj(lp, j, glp_get_sjj(lp, j) / sqrt(temp));
            }
         }
      }
      return;
}

/***********************************************************************
*  max_row_ratio - determine worst scaling "quality" for rows
*
*  This routine returns the worst scaling "quality" for rows of the
*  currently scaled constraint matrix:
*
*              m
*     ratio = max ratio[i],
*             i=1
*  where:
*                 n              n
*     ratio[i] = max |a[i,j]| / min |a[i,j]|,  1 <= i <= m,
*                j=1            j=1
*
*  is the scaling "quality" of i-th row. */

static double max_row_ratio(glp_prob *lp)
{     int i;
      double ratio, temp;
      ratio = 1.0;
      for (i = 1; i <= lp->m; i++)
      {  temp = max_row_aij(lp, i, 1) / min_row_aij(lp, i, 1);
         if (i == 1 || ratio < temp) ratio = temp;
      }
      return ratio;
}

/***********************************************************************
*  max_col_ratio - determine worst scaling "quality" for columns
*
*  This routine returns the worst scaling "quality" for columns of the
*  currently scaled constraint matrix:
*
*              n
*     ratio = max ratio[j],
*             j=1
*  where:
*                 m              m
*     ratio[j] = max |a[i,j]| / min |a[i,j]|,  1 <= j <= n,
*                i=1            i=1
*
*  is the scaling "quality" of j-th column. */

static double max_col_ratio(glp_prob *lp)
{     int j;
      double ratio, temp;
      ratio = 1.0;
      for (j = 1; j <= lp->n; j++)
      {  temp = max_col_aij(lp, j, 1) / min_col_aij(lp, j, 1);
         if (j == 1 || ratio < temp) ratio = temp;
      }
      return ratio;
}

/***********************************************************************
*  gm_iterate - perform iterative geometric mean scaling
*
*  This routine performs iterative geometric mean scaling of rows and
*  columns of the constraint matrix.
*
*  The parameter it_max specifies the maximal number of iterations.
*  Recommended value of it_max is 15.
*
*  The parameter tau specifies a minimal improvement of the scaling
*  "quality" on each iteration, 0 < tau < 1. It means than the scaling
*  process continues while the following condition is satisfied:
*
*     ratio[k] <= tau * ratio[k-1],
*
*  where ratio = max |a[i,j]| / min |a[i,j]| is the scaling "quality"
*  to be minimized, k is the iteration number. Recommended value of tau
*  is 0.90. */

static void gm_iterate(glp_prob *lp, int it_max, double tau)
{     int k, flag;
      double ratio = 0.0, r_old;
      /* if the scaling "quality" for rows is better than for columns,
         the rows are scaled first; otherwise, the columns are scaled
         first */
      flag = (max_row_ratio(lp) > max_col_ratio(lp));
      for (k = 1; k <= it_max; k++)
      {  /* save the scaling "quality" from previous iteration */
         r_old = ratio;
         /* determine the current scaling "quality" */
         ratio = max_mat_aij(lp, 1) / min_mat_aij(lp, 1);
#if 0
         xprintf("k = %d; ratio = %g\n", k, ratio);
#endif
         /* if improvement is not enough, terminate scaling */
         if (k > 1 && ratio > tau * r_old) break;
         /* otherwise, perform another iteration */
         gm_scaling(lp, flag);
      }
      return;
}

/***********************************************************************
*  NAME
*
*  scale_prob - scale problem data
*
*  SYNOPSIS
*
*  #include "glpscl.h"
*  void scale_prob(glp_prob *lp, int flags);
*
*  DESCRIPTION
*
*  The routine scale_prob performs automatic scaling of problem data
*  for the specified problem object. */

static void scale_prob(glp_prob *lp, int flags)
{     static const char *fmt =
         "%s: min|aij| = %10.3e  max|aij| = %10.3e  ratio = %10.3e\n";
      double min_aij, max_aij, ratio;
      xprintf("Scaling...\n");
      /* cancel the current scaling effect */
      glp_unscale_prob(lp);
      /* report original scaling "quality" */
      min_aij = min_mat_aij(lp, 1);
      max_aij = max_mat_aij(lp, 1);
      ratio = max_aij / min_aij;
      xprintf(fmt, " A", min_aij, max_aij, ratio);
      /* check if the problem is well scaled */
      if (min_aij >= 0.10 && max_aij <= 10.0)
      {  xprintf("Problem data seem to be well scaled\n");
         /* skip scaling, if required */
         if (flags & GLP_SF_SKIP) goto done;
      }
      /* perform iterative geometric mean scaling, if required */
      if (flags & GLP_SF_GM)
      {  gm_iterate(lp, 15, 0.90);
         min_aij = min_mat_aij(lp, 1);
         max_aij = max_mat_aij(lp, 1);
         ratio = max_aij / min_aij;
         xprintf(fmt, "GM", min_aij, max_aij, ratio);
      }
      /* perform equilibration scaling, if required */
      if (flags & GLP_SF_EQ)
      {  eq_scaling(lp, max_row_ratio(lp) > max_col_ratio(lp));
         min_aij = min_mat_aij(lp, 1);
         max_aij = max_mat_aij(lp, 1);
         ratio = max_aij / min_aij;
         xprintf(fmt, "EQ", min_aij, max_aij, ratio);
      }
      /* round scale factors to nearest power of two, if required */
      if (flags & GLP_SF_2N)
      {  int i, j;
         for (i = 1; i <= lp->m; i++)
            glp_set_rii(lp, i, round2n(glp_get_rii(lp, i)));
         for (j = 1; j <= lp->n; j++)
            glp_set_sjj(lp, j, round2n(glp_get_sjj(lp, j)));
         min_aij = min_mat_aij(lp, 1);
         max_aij = max_mat_aij(lp, 1);
         ratio = max_aij / min_aij;
         xprintf(fmt, "2N", min_aij, max_aij, ratio);
      }
done: return;
}

/***********************************************************************
*  NAME
*
*  glp_scale_prob - scale problem data
*
*  SYNOPSIS
*
*  void glp_scale_prob(glp_prob *lp, int flags);
*
*  DESCRIPTION
*
*  The routine glp_scale_prob performs automatic scaling of problem
*  data for the specified problem object.
*
*  The parameter flags specifies scaling options used by the routine.
*  Options can be combined with the bitwise OR operator and may be the
*  following:
*
*  GLP_SF_GM      perform geometric mean scaling;
*  GLP_SF_EQ      perform equilibration scaling;
*  GLP_SF_2N      round scale factors to nearest power of two;
*  GLP_SF_SKIP    skip scaling, if the problem is well scaled.
*
*  The parameter flags may be specified as GLP_SF_AUTO, in which case
*  the routine chooses scaling options automatically. */

void glp_scale_prob(glp_prob *lp, int flags)
{     if (flags & ~(GLP_SF_GM | GLP_SF_EQ | GLP_SF_2N | GLP_SF_SKIP |
                    GLP_SF_AUTO))
         xerror("glp_scale_prob: flags = 0x%02X; invalid scaling option"
            "s\n", flags);
      if (flags & GLP_SF_AUTO)
         flags = (GLP_SF_GM | GLP_SF_EQ | GLP_SF_SKIP);
      scale_prob(lp, flags);
      return;
}

/* eof */
