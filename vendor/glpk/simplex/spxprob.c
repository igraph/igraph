/* spxprob.c */

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
#include "spxprob.h"

/***********************************************************************
*  spx_init_lp - initialize working LP object
*
*  This routine determines the number of equality constraints m, the
*  number of variables n, and the number of non-zero elements nnz in
*  the constraint matrix for the working LP, which corresponds to the
*  original LP, and stores these dimensions to the working LP object.
*  (The working LP object should be allocated by the calling routine.)
*
*  If the flag excl is set, the routine assumes that non-basic fixed
*  variables will be excluded from the working LP. */

void spx_init_lp(SPXLP *lp, glp_prob *P, int excl)
{     int i, j, m, n, nnz;
      m = P->m;
      xassert(m > 0);
      n = 0;
      nnz = P->nnz;
      xassert(P->valid);
      /* scan rows of original LP */
      for (i = 1; i <= m; i++)
      {  GLPROW *row = P->row[i];
         if (excl && row->stat == GLP_NS)
         {  /* skip non-basic fixed auxiliary variable */
            /* nop */
         }
         else
         {  /* include auxiliary variable in working LP */
            n++;
            nnz++; /* unity column */
         }
      }
      /* scan columns of original LP */
      for (j = 1; j <= P->n; j++)
      {  GLPCOL *col = P->col[j];
         if (excl && col->stat == GLP_NS)
         {  /* skip non-basic fixed structural variable */
            GLPAIJ *aij;
            for (aij = col->ptr; aij != NULL; aij = aij->c_next)
               nnz--;
         }
         else
         {  /* include structural variable in working LP */
            n++;
         }
      }
      /* initialize working LP data block */
      memset(lp, 0, sizeof(SPXLP));
      lp->m = m;
      xassert(n > 0);
      lp->n = n;
      lp->nnz = nnz;
      return;
}

/***********************************************************************
*  spx_alloc_lp - allocate working LP arrays
*
*  This routine allocates the memory for all arrays in the working LP
*  object. */

void spx_alloc_lp(SPXLP *lp)
{     int m = lp->m;
      int n = lp->n;
      int nnz = lp->nnz;
      lp->A_ptr = talloc(1+n+1, int);
      lp->A_ind = talloc(1+nnz, int);
      lp->A_val = talloc(1+nnz, double);
      lp->b = talloc(1+m, double);
      lp->c = talloc(1+n, double);
      lp->l = talloc(1+n, double);
      lp->u = talloc(1+n, double);
      lp->head = talloc(1+n, int);
      lp->flag = talloc(1+n-m, char);
      return;
}

/***********************************************************************
*  spx_build_lp - convert original LP to working LP
*
*  This routine converts components (except the current basis) of the
*  original LP to components of the working LP and perform scaling of
*  these components. Also, if the original LP is maximization, the
*  routine changes the signs of the objective coefficients and constant
*  term to opposite ones.
*
*  If the flag excl is set, original non-basic fixed variables are
*  *not* included in the working LP. Otherwise, all (auxiliary and
*  structural) original variables are included in the working LP. Note
*  that this flag should have the same value as it has in a call to the
*  routine spx_init_lp.
*
*  If the flag shift is set, the routine shift bounds of variables
*  included in the working LP to make at least one bound to be zero.
*  If a variable has both lower and upper bounds, the bound having
*  smaller magnitude is shifted to zero.
*
*  On exit the routine stores information about correspondence between
*  numbers of variables in the original and working LPs to the array
*  map, which should have 1+P->m+P->n locations (location [0] is not
*  used), where P->m is the numbers of rows and P->n is the number of
*  columns in the original LP:
*
*  map[i] = +k, 1 <= i <= P->m, means that i-th auxiliary variable of
*  the original LP corresponds to variable x[k] of the working LP;
*
*  map[i] = -k, 1 <= i <= P->m, means that i-th auxiliary variable of
*  the original LP corresponds to variable x[k] of the working LP, and
*  the upper bound of that variable was shifted to zero;
*
*  map[i] = 0, 1 <= i <= P->m, means that i-th auxiliary variable of
*  the original LP was excluded from the working LP;
*
*  map[P->m+j], 1 <= j <= P->n, has the same sense as above, however,
*  for j-th structural variable of the original LP. */

void spx_build_lp(SPXLP *lp, glp_prob *P, int excl, int shift,
      int map[/*1+P->m+P->n*/])
{     int m = lp->m;
      int n = lp->n;
      int nnz = lp->nnz;
      int *A_ptr = lp->A_ptr;
      int *A_ind = lp->A_ind;
      double *A_val = lp->A_val;
      double *b = lp->b;
      double *c = lp->c;
      double *l = lp->l;
      double *u = lp->u;
      int i, j, k, kk, ptr, end;
      double dir, delta;
      /* working LP is always minimization */
      switch (P->dir)
      {  case GLP_MIN:
            dir = +1.0;
            break;
         case GLP_MAX:
            dir = -1.0;
            break;
         default:
            xassert(P != P);
      }
      /* initialize constant term of the objective */
      c[0] = dir * P->c0;
      k = 0; /* number of variable in working LP */
      ptr = 1; /* current available position in A_ind/A_val */
      /* process rows of original LP */
      xassert(P->m == m);
      for (i = 1; i <= m; i++)
      {  GLPROW *row = P->row[i];
         if (excl && row->stat == GLP_NS)
         {  /* i-th auxiliary variable is non-basic and fixed */
            /* substitute its scaled value in working LP */
            xassert(row->type == GLP_FX);
            map[i] = 0;
            b[i] = - row->lb * row->rii;
         }
         else
         {  /* include i-th auxiliary variable in working LP */
            map[i] = ++k;
            /* setup k-th column of working constraint matrix which is
             * i-th column of unity matrix */
            A_ptr[k] = ptr;
            A_ind[ptr] = i;
            A_val[ptr] = 1.0;
            ptr++;
            /* initialize right-hand side of i-th equality constraint
             * and setup zero objective coefficient at variable x[k] */
            b[i] = c[k] = 0.0;
            /* setup scaled bounds of variable x[k] */
            switch (row->type)
            {  case GLP_FR:
                  l[k] = -DBL_MAX, u[k] = +DBL_MAX;
                  break;
               case GLP_LO:
                  l[k] = row->lb * row->rii, u[k] = +DBL_MAX;
                  break;
               case GLP_UP:
                  l[k] = -DBL_MAX, u[k] = row->ub * row->rii;
                  break;
               case GLP_DB:
                  l[k] = row->lb * row->rii, u[k] = row->ub * row->rii;
                  xassert(l[k] != u[k]);
                  break;
               case GLP_FX:
                  l[k] = u[k] = row->lb * row->rii;
                  break;
               default:
                  xassert(row != row);
            }
         }
      }
      /* process columns of original LP */
      for (j = 1; j <= P->n; j++)
      {  GLPCOL *col = P->col[j];
         GLPAIJ *aij;
         if (excl && col->stat == GLP_NS)
         {  /* j-th structural variable is non-basic and fixed */
            /* substitute its scaled value in working LP */
            xassert(col->type == GLP_FX);
            map[m+j] = 0;
            if (col->lb != 0.0)
            {  /* (note that sjj scale factor is cancelled) */
               for (aij = col->ptr; aij != NULL; aij = aij->c_next)
                  b[aij->row->i] +=
                     (aij->row->rii * aij->val) * col->lb;
               c[0] += (dir * col->coef) * col->lb;
            }
         }
         else
         {  /* include j-th structural variable in working LP */
            map[m+j] = ++k;
            /* setup k-th column of working constraint matrix which is
             * scaled j-th column of original constraint matrix (-A) */
            A_ptr[k] = ptr;
            for (aij = col->ptr; aij != NULL; aij = aij->c_next)
            {  A_ind[ptr] = aij->row->i;
               A_val[ptr] = - aij->row->rii * aij->val * col->sjj;
               ptr++;
            }
            /* setup scaled objective coefficient at variable x[k] */
            c[k] = dir * col->coef * col->sjj;
            /* setup scaled bounds of variable x[k] */
            switch (col->type)
            {  case GLP_FR:
                  l[k] = -DBL_MAX, u[k] = +DBL_MAX;
                  break;
               case GLP_LO:
                  l[k] = col->lb / col->sjj, u[k] = +DBL_MAX;
                  break;
               case GLP_UP:
                  l[k] = -DBL_MAX, u[k] = col->ub / col->sjj;
                  break;
               case GLP_DB:
                  l[k] = col->lb / col->sjj, u[k] = col->ub / col->sjj;
                  xassert(l[k] != u[k]);
                  break;
               case GLP_FX:
                  l[k] = u[k] = col->lb / col->sjj;
                  break;
               default:
                  xassert(col != col);
            }
         }
      }
      xassert(k == n);
      xassert(ptr == nnz+1);
      A_ptr[n+1] = ptr;
      /* shift bounds of all variables of working LP (optionally) */
      if (shift)
      {  for (kk = 1; kk <= m+P->n; kk++)
         {  k = map[kk];
            if (k == 0)
            {  /* corresponding original variable was excluded */
               continue;
            }
            /* shift bounds of variable x[k] */
            if (l[k] == -DBL_MAX && u[k] == +DBL_MAX)
            {  /* x[k] is unbounded variable */
               delta = 0.0;
            }
            else if (l[k] != -DBL_MAX && u[k] == +DBL_MAX)
            {  /* shift lower bound to zero */
               delta = l[k];
               l[k] = 0.0;
            }
            else if (l[k] == -DBL_MAX && u[k] != +DBL_MAX)
            {  /* shift upper bound to zero */
               map[kk] = -k;
               delta = u[k];
               u[k] = 0.0;
            }
            else if (l[k] != u[k])
            {  /* x[k] is double bounded variable */
               if (fabs(l[k]) <= fabs(u[k]))
               {  /* shift lower bound to zero */
                  delta = l[k];
                  l[k] = 0.0, u[k] -= delta;
               }
               else
               {  /* shift upper bound to zero */
                  map[kk] = -k;
                  delta = u[k];
                  l[k] -= delta, u[k] = 0.0;
               }
               xassert(l[k] != u[k]);
            }
            else
            {  /* shift fixed value to zero */
               delta = l[k];
               l[k] = u[k] = 0.0;
            }
            /* substitute x[k] = x'[k] + delta into all constraints
             * and the objective function of working LP */
            if (delta != 0.0)
            {  ptr = A_ptr[k];
               end = A_ptr[k+1];
               for (; ptr < end; ptr++)
                  b[A_ind[ptr]] -= A_val[ptr] * delta;
               c[0] += c[k] * delta;
            }
         }
      }
      return;
}

/***********************************************************************
*  spx_build_basis - convert original LP basis to working LP basis
*
*  This routine converts the current basis of the original LP to
*  corresponding initial basis of the working LP, and moves the basis
*  factorization driver from the original LP object to the working LP
*  object.
*
*  The array map should contain information provided by the routine
*  spx_build_lp. */

void spx_build_basis(SPXLP *lp, glp_prob *P, const int map[])
{     int m = lp->m;
      int n = lp->n;
      int *head = lp->head;
      char *flag = lp->flag;
      int i, j, k, ii, jj;
      /* original basis factorization should be valid that guarantees
       * the basis is correct */
      xassert(P->m == m);
      xassert(P->valid);
      /* initialize basis header for working LP */
      memset(&head[1], 0, m * sizeof(int));
      jj = 0;
      /* scan rows of original LP */
      xassert(P->m == m);
      for (i = 1; i <= m; i++)
      {  GLPROW *row = P->row[i];
         /* determine ordinal number of x[k] in working LP */
         if ((k = map[i]) < 0)
            k = -k;
         if (k == 0)
         {  /* corresponding original variable was excluded */
            continue;
         }
         xassert(1 <= k && k <= n);
         if (row->stat == GLP_BS)
         {  /* x[k] is basic variable xB[ii] */
            ii = row->bind;
            xassert(1 <= ii && ii <= m);
            xassert(head[ii] == 0);
            head[ii] = k;
         }
         else
         {  /* x[k] is non-basic variable xN[jj] */
            jj++;
            head[m+jj] = k;
            flag[jj] = (row->stat == GLP_NU);
         }
      }
      /* scan columns of original LP */
      for (j = 1; j <= P->n; j++)
      {  GLPCOL *col = P->col[j];
         /* determine ordinal number of x[k] in working LP */
         if ((k = map[m+j]) < 0)
            k = -k;
         if (k == 0)
         {  /* corresponding original variable was excluded */
            continue;
         }
         xassert(1 <= k && k <= n);
         if (col->stat == GLP_BS)
         {  /* x[k] is basic variable xB[ii] */
            ii = col->bind;
            xassert(1 <= ii && ii <= m);
            xassert(head[ii] == 0);
            head[ii] = k;
         }
         else
         {  /* x[k] is non-basic variable xN[jj] */
            jj++;
            head[m+jj] = k;
            flag[jj] = (col->stat == GLP_NU);
         }
      }
      xassert(m+jj == n);
      /* acquire basis factorization */
      lp->valid = 1;
      lp->bfd = P->bfd;
      P->valid = 0;
      P->bfd = NULL;
      return;
}

/***********************************************************************
*  spx_store_basis - convert working LP basis to original LP basis
*
*  This routine converts the current working LP basis to corresponding
*  original LP basis. This operations includes determining and setting
*  statuses of all rows (auxiliary variables) and columns (structural
*  variables), and building the basis header.
*
*  The array map should contain information provided by the routine
*  spx_build_lp.
*
*  On exit the routine fills the array daeh. This array should have
*  1+lp->n locations (location [0] is not used) and contain the inverse
*  of the working basis header lp->head, i.e. head[k'] = k means that
*  daeh[k] = k'. */

void spx_store_basis(SPXLP *lp, glp_prob *P, const int map[],
      int daeh[/*1+n*/])
{     int m = lp->m;
      int n = lp->n;
      int *head = lp->head;
      char *flag = lp->flag;
      int i, j, k, kk;
      /* determine inverse of working basis header */
      for (kk = 1; kk <= n; kk++)
         daeh[head[kk]] = kk;
      /* set row statuses */
      xassert(P->m == m);
      for (i = 1; i <= m; i++)
      {  GLPROW *row = P->row[i];
         if ((k = map[i]) < 0)
            k = -k;
         if (k == 0)
         {  /* non-basic fixed auxiliary variable was excluded */
            xassert(row->type == GLP_FX);
            row->stat = GLP_NS;
            row->bind = 0;
         }
         else
         {  /* auxiliary variable corresponds to variable x[k] */
            kk = daeh[k];
            if (kk <= m)
            {  /* x[k] = xB[kk] */
               P->head[kk] = i;
               row->stat = GLP_BS;
               row->bind = kk;
            }
            else
            {  /* x[k] = xN[kk-m] */
               switch (row->type)
               {  case GLP_FR:
                     row->stat = GLP_NF;
                     break;
                  case GLP_LO:
                     row->stat = GLP_NL;
                     break;
                  case GLP_UP:
                     row->stat = GLP_NU;
                     break;
                  case GLP_DB:
                     row->stat = (flag[kk-m] ? GLP_NU : GLP_NL);
                     break;
                  case GLP_FX:
                     row->stat = GLP_NS;
                     break;
                  default:
                     xassert(row != row);
               }
               row->bind = 0;
            }
         }
      }
      /* set column statuses */
      for (j = 1; j <= P->n; j++)
      {  GLPCOL *col = P->col[j];
         if ((k = map[m+j]) < 0)
            k = -k;
         if (k == 0)
         {  /* non-basic fixed structural variable was excluded */
            xassert(col->type == GLP_FX);
            col->stat = GLP_NS;
            col->bind = 0;
         }
         else
         {  /* structural variable corresponds to variable x[k] */
            kk = daeh[k];
            if (kk <= m)
            {  /* x[k] = xB[kk] */
               P->head[kk] = m+j;
               col->stat = GLP_BS;
               col->bind = kk;
            }
            else
            {  /* x[k] = xN[kk-m] */
               switch (col->type)
               {  case GLP_FR:
                     col->stat = GLP_NF;
                     break;
                  case GLP_LO:
                     col->stat = GLP_NL;
                     break;
                  case GLP_UP:
                     col->stat = GLP_NU;
                     break;
                  case GLP_DB:
                     col->stat = (flag[kk-m] ? GLP_NU : GLP_NL);
                     break;
                  case GLP_FX:
                     col->stat = GLP_NS;
                     break;
                  default:
                     xassert(col != col);
               }
               col->bind = 0;
            }
         }
      }
      return;
}

/***********************************************************************
*  spx_store_sol - convert working LP solution to original LP solution
*
*  This routine converts the current basic solution of the working LP
*  (values of basic variables, simplex multipliers, reduced costs of
*  non-basic variables) to corresponding basic solution of the original
*  LP (values and reduced costs of auxiliary and structural variables).
*  This conversion includes unscaling all basic solution components,
*  computing reduced costs of excluded non-basic variables, recovering
*  unshifted values of basic variables, changing the signs of reduced
*  costs (if the original LP is maximization), and computing the value
*  of the objective function.
*
*  The flag shift should have the same value as it has in a call to the
*  routine spx_build_lp.
*
*  The array map should contain information provided by the routine
*  spx_build_lp.
*
*  The array daeh should contain information provided by the routine
*  spx_store_basis.
*
*  The arrays beta, pi, and d should contain basic solution components
*  for the working LP:
*
*  array locations beta[1], ..., beta[m] should contain values of basic
*  variables beta = (beta[i]);
*
*  array locations pi[1], ..., pi[m] should contain simplex multipliers
*  pi = (pi[i]);
*
*  array locations d[1], ..., d[n-m] should contain reduced costs of
*  non-basic variables d = (d[j]). */

void spx_store_sol(SPXLP *lp, glp_prob *P, int shift,
      const int map[], const int daeh[], const double beta[],
      const double pi[], const double d[])
{     int m = lp->m;
      char *flag = lp->flag;
      int i, j, k, kk;
      double dir;
      /* working LP is always minimization */
      switch (P->dir)
      {  case GLP_MIN:
            dir = +1.0;
            break;
         case GLP_MAX:
            dir = -1.0;
            break;
         default:
            xassert(P != P);
      }
      /* compute row solution components */
      xassert(P->m == m);
      for (i = 1; i <= m; i++)
      {  GLPROW *row = P->row[i];
         if ((k = map[i]) < 0)
            k = -k;
         if (k == 0)
         {  /* non-basic fixed auxiliary variable was excluded */
            xassert(row->type == GLP_FX);
            row->prim = row->lb;
            /* compute reduced cost d[k] = c[k] - A'[k] * pi as if x[k]
             * would be non-basic in working LP */
            row->dual = - dir * pi[i] * row->rii;
         }
         else
         {  /* auxiliary variable corresponds to variable x[k] */
            kk = daeh[k];
            if (kk <= m)
            {  /* x[k] = xB[kk] */
               row->prim = beta[kk] / row->rii;
               if (shift)
                  row->prim += (map[i] < 0 ? row->ub : row->lb);
               row->dual = 0.0;
            }
            else
            {  /* x[k] = xN[kk-m] */
               row->prim = (flag[kk-m] ? row->ub : row->lb);
               row->dual = (dir * d[kk-m]) * row->rii;
            }
         }
      }
      /* compute column solution components and objective value */
      P->obj_val = P->c0;
      for (j = 1; j <= P->n; j++)
      {  GLPCOL *col = P->col[j];
         if ((k = map[m+j]) < 0)
            k = -k;
         if (k == 0)
         {  /* non-basic fixed structural variable was excluded */
            GLPAIJ *aij;
            double dk;
            xassert(col->type == GLP_FX);
            col->prim = col->lb;
            /* compute reduced cost d[k] = c[k] - A'[k] * pi as if x[k]
             * would be non-basic in working LP */
            /* (note that sjj scale factor is cancelled) */
            dk = dir * col->coef;
            for (aij = col->ptr; aij != NULL; aij = aij->c_next)
               dk += (aij->row->rii * aij->val) * pi[aij->row->i];
            col->dual = dir * dk;
         }
         else
         {  /* structural variable corresponds to variable x[k] */
            kk = daeh[k];
            if (kk <= m)
            {  /* x[k] = xB[kk] */
               col->prim = beta[kk] * col->sjj;
               if (shift)
                  col->prim += (map[m+j] < 0 ? col->ub : col->lb);
               col->dual = 0.0;
            }
            else
            {  /* x[k] = xN[kk-m] */
               col->prim = (flag[kk-m] ? col->ub : col->lb);
               col->dual = (dir * d[kk-m]) / col->sjj;
            }
         }
         P->obj_val += col->coef * col->prim;
      }
      return;
}

/***********************************************************************
*  spx_free_lp - deallocate working LP arrays
*
*  This routine deallocates the memory used for arrays of the working
*  LP object. */

void spx_free_lp(SPXLP *lp)
{     tfree(lp->A_ptr);
      tfree(lp->A_ind);
      tfree(lp->A_val);
      tfree(lp->b);
      tfree(lp->c);
      tfree(lp->l);
      tfree(lp->u);
      tfree(lp->head);
      tfree(lp->flag);
      return;
}

/* eof */
