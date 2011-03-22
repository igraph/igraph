/* glpini02.c */

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

struct var
{     /* structural variable */
      int j;
      /* ordinal number */
      double q;
      /* penalty value */
};

static int fcmp(const void *ptr1, const void *ptr2)
{     /* this routine is passed to the qsort() function */
      struct var *col1 = (void *)ptr1, *col2 = (void *)ptr2;
      if (col1->q < col2->q) return -1;
      if (col1->q > col2->q) return +1;
      return 0;
}

static int get_column(glp_prob *lp, int j, int ind[], double val[])
{     /* Bixby's algorithm assumes that the constraint matrix is scaled
         such that the maximum absolute value in every non-zero row and
         column is 1 */
      int k, len;
      double big;
      len = glp_get_mat_col(lp, j, ind, val);
      big = 0.0;
      for (k = 1; k <= len; k++)
         if (big < fabs(val[k])) big = fabs(val[k]);
      if (big == 0.0) big = 1.0;
      for (k = 1; k <= len; k++) val[k] /= big;
      return len;
}

static void cpx_basis(glp_prob *lp)
{     /* main routine */
      struct var *C, *C2, *C3, *C4;
      int m, n, i, j, jk, k, l, ll, t, n2, n3, n4, type, len, *I, *r,
         *ind;
      double alpha, gamma, cmax, temp, *v, *val;
      xprintf("Constructing initial basis...\n");
      /* determine the number of rows and columns */
      m = glp_get_num_rows(lp);
      n = glp_get_num_cols(lp);
      /* allocate working arrays */
      C = xcalloc(1+n, sizeof(struct var));
      I = xcalloc(1+m, sizeof(int));
      r = xcalloc(1+m, sizeof(int));
      v = xcalloc(1+m, sizeof(double));
      ind = xcalloc(1+m, sizeof(int));
      val = xcalloc(1+m, sizeof(double));
      /* make all auxiliary variables non-basic */
      for (i = 1; i <= m; i++)
      {  if (glp_get_row_type(lp, i) != GLP_DB)
            glp_set_row_stat(lp, i, GLP_NS);
         else if (fabs(glp_get_row_lb(lp, i)) <=
                  fabs(glp_get_row_ub(lp, i)))
            glp_set_row_stat(lp, i, GLP_NL);
         else
            glp_set_row_stat(lp, i, GLP_NU);
      }
      /* make all structural variables non-basic */
      for (j = 1; j <= n; j++)
      {  if (glp_get_col_type(lp, j) != GLP_DB)
            glp_set_col_stat(lp, j, GLP_NS);
         else if (fabs(glp_get_col_lb(lp, j)) <=
                  fabs(glp_get_col_ub(lp, j)))
            glp_set_col_stat(lp, j, GLP_NL);
         else
            glp_set_col_stat(lp, j, GLP_NU);
      }
      /* C2 is a set of free structural variables */
      n2 = 0, C2 = C + 0;
      for (j = 1; j <= n; j++)
      {  type = glp_get_col_type(lp, j);
         if (type == GLP_FR)
         {  n2++;
            C2[n2].j = j;
            C2[n2].q = 0.0;
         }
      }
      /* C3 is a set of structural variables having excatly one (lower
         or upper) bound */
      n3 = 0, C3 = C2 + n2;
      for (j = 1; j <= n; j++)
      {  type = glp_get_col_type(lp, j);
         if (type == GLP_LO)
         {  n3++;
            C3[n3].j = j;
            C3[n3].q = + glp_get_col_lb(lp, j);
         }
         else if (type == GLP_UP)
         {  n3++;
            C3[n3].j = j;
            C3[n3].q = - glp_get_col_ub(lp, j);
         }
      }
      /* C4 is a set of structural variables having both (lower and
         upper) bounds */
      n4 = 0, C4 = C3 + n3;
      for (j = 1; j <= n; j++)
      {  type = glp_get_col_type(lp, j);
         if (type == GLP_DB)
         {  n4++;
            C4[n4].j = j;
            C4[n4].q = glp_get_col_lb(lp, j) - glp_get_col_ub(lp, j);
         }
      }
      /* compute gamma = max{|c[j]|: 1 <= j <= n} */
      gamma = 0.0;
      for (j = 1; j <= n; j++)
      {  temp = fabs(glp_get_obj_coef(lp, j));
         if (gamma < temp) gamma = temp;
      }
      /* compute cmax */
      cmax = (gamma == 0.0 ? 1.0 : 1000.0 * gamma);
      /* compute final penalty for all structural variables within sets
         C2, C3, and C4 */
      switch (glp_get_obj_dir(lp))
      {  case GLP_MIN: temp = +1.0; break;
         case GLP_MAX: temp = -1.0; break;
         default: xassert(lp != lp);
      }
      for (k = 1; k <= n2+n3+n4; k++)
      {  j = C[k].j;
         C[k].q += (temp * glp_get_obj_coef(lp, j)) / cmax;
      }
      /* sort structural variables within C2, C3, and C4 in ascending
         order of penalty value */
      qsort(C2+1, n2, sizeof(struct var), fcmp);
      for (k = 1; k < n2; k++) xassert(C2[k].q <= C2[k+1].q);
      qsort(C3+1, n3, sizeof(struct var), fcmp);
      for (k = 1; k < n3; k++) xassert(C3[k].q <= C3[k+1].q);
      qsort(C4+1, n4, sizeof(struct var), fcmp);
      for (k = 1; k < n4; k++) xassert(C4[k].q <= C4[k+1].q);
      /*** STEP 1 ***/
      for (i = 1; i <= m; i++)
      {  type = glp_get_row_type(lp, i);
         if (type != GLP_FX)
         {  /* row i is either free or inequality constraint */
            glp_set_row_stat(lp, i, GLP_BS);
            I[i] = 1;
            r[i] = 1;
         }
         else
         {  /* row i is equality constraint */
            I[i] = 0;
            r[i] = 0;
         }
         v[i] = +DBL_MAX;
      }
      /*** STEP 2 ***/
      for (k = 1; k <= n2+n3+n4; k++)
      {  jk = C[k].j;
         len = get_column(lp, jk, ind, val);
         /* let alpha = max{|A[l,jk]|: r[l] = 0} and let l' be such
            that alpha = |A[l',jk]| */
         alpha = 0.0, ll = 0;
         for (t = 1; t <= len; t++)
         {  l = ind[t];
            if (r[l] == 0 && alpha < fabs(val[t]))
               alpha = fabs(val[t]), ll = l;
         }
         if (alpha >= 0.99)
         {  /* B := B union {jk} */
            glp_set_col_stat(lp, jk, GLP_BS);
            I[ll] = 1;
            v[ll] = alpha;
            /* r[l] := r[l] + 1 for all l such that |A[l,jk]| != 0 */
            for (t = 1; t <= len; t++)
            {  l = ind[t];
               if (val[t] != 0.0) r[l]++;
            }
            /* continue to the next k */
            continue;
         }
         /* if |A[l,jk]| > 0.01 * v[l] for some l, continue to the
            next k */
         for (t = 1; t <= len; t++)
         {  l = ind[t];
            if (fabs(val[t]) > 0.01 * v[l]) break;
         }
         if (t <= len) continue;
         /* otherwise, let alpha = max{|A[l,jk]|: I[l] = 0} and let l'
            be such that alpha = |A[l',jk]| */
         alpha = 0.0, ll = 0;
         for (t = 1; t <= len; t++)
         {  l = ind[t];
            if (I[l] == 0 && alpha < fabs(val[t]))
               alpha = fabs(val[t]), ll = l;
         }
         /* if alpha = 0, continue to the next k */
         if (alpha == 0.0) continue;
         /* B := B union {jk} */
         glp_set_col_stat(lp, jk, GLP_BS);
         I[ll] = 1;
         v[ll] = alpha;
         /* r[l] := r[l] + 1 for all l such that |A[l,jk]| != 0 */
         for (t = 1; t <= len; t++)
         {  l = ind[t];
            if (val[t] != 0.0) r[l]++;
         }
      }
      /*** STEP 3 ***/
      /* add an artificial variable (auxiliary variable for equality
         constraint) to cover each remaining uncovered row */
      for (i = 1; i <= m; i++)
         if (I[i] == 0) glp_set_row_stat(lp, i, GLP_BS);
      /* free working arrays */
      xfree(C);
      xfree(I);
      xfree(r);
      xfree(v);
      xfree(ind);
      xfree(val);
      return;
}

/***********************************************************************
*  NAME
*
*  glp_cpx_basis - construct Bixby's initial LP basis
*
*  SYNOPSIS
*
*  void glp_cpx_basis(glp_prob *lp);
*
*  DESCRIPTION
*
*  The routine glp_cpx_basis constructs an advanced initial basis for
*  the specified problem object.
*
*  The routine is based on Bixby's algorithm described in the paper:
*
*  Robert E. Bixby. Implementing the Simplex Method: The Initial Basis.
*  ORSA Journal on Computing, Vol. 4, No. 3, 1992, pp. 267-84. */

void glp_cpx_basis(glp_prob *lp)
{     if (lp->m == 0 || lp->n == 0)
         glp_std_basis(lp);
      else
         cpx_basis(lp);
      return;
}

/* eof */
