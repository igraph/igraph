/* mirgen.c (mixed integer rounding cuts generator) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2007-2018 Free Software Foundation, Inc.
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

#if 1 /* 29/II-2016 by Chris */
/*----------------------------------------------------------------------
Subject: Mir cut generation performance improvement
From: Chris Matrakidis <cmatraki@gmail.com>
To: Andrew Makhorin <mao@gnu.org>, help-glpk <help-glpk@gnu.org>

Andrew,

I noticed that mir cut generation takes considerable time on some large
problems (like rocII-4-11 from miplib). The attached patch makes two
improvements that considerably improve performance in such instances:
1. A lot of time was spent on generating a temporary vector in function
aggregate_row. It is a lot faster to reuse an existing vector.
2. A search for an element in the same function was done in row order,
where using the elements in the order they are in the column is more
efficient. This changes the generated cuts in some cases, but seems
neutral overall (0.3% less cuts in a test set of 64 miplib instances).

Best Regards,

Chris Matrakidis
----------------------------------------------------------------------*/
#endif

#include "env.h"
#include "prob.h"
#include "spv.h"

#define MIR_DEBUG 0

#define MAXAGGR 5
/* maximal number of rows that can be aggregated */

struct glp_mir
{     /* MIR cut generator working area */
      /*--------------------------------------------------------------*/
      /* global information valid for the root subproblem */
      int m;
      /* number of rows (in the root subproblem) */
      int n;
      /* number of columns */
      char *skip; /* char skip[1+m]; */
      /* skip[i], 1 <= i <= m, is a flag that means that row i should
         not be used because (1) it is not suitable, or (2) because it
         has been used in the aggregated constraint */
      char *isint; /* char isint[1+m+n]; */
      /* isint[k], 1 <= k <= m+n, is a flag that means that variable
         x[k] is integer (otherwise, continuous) */
      double *lb; /* double lb[1+m+n]; */
      /* lb[k], 1 <= k <= m+n, is lower bound of x[k]; -DBL_MAX means
         that x[k] has no lower bound */
      int *vlb; /* int vlb[1+m+n]; */
      /* vlb[k] = k', 1 <= k <= m+n, is the number of integer variable,
         which defines variable lower bound x[k] >= lb[k] * x[k']; zero
         means that x[k] has simple lower bound */
      double *ub; /* double ub[1+m+n]; */
      /* ub[k], 1 <= k <= m+n, is upper bound of x[k]; +DBL_MAX means
         that x[k] has no upper bound */
      int *vub; /* int vub[1+m+n]; */
      /* vub[k] = k', 1 <= k <= m+n, is the number of integer variable,
         which defines variable upper bound x[k] <= ub[k] * x[k']; zero
         means that x[k] has simple upper bound */
      /*--------------------------------------------------------------*/
      /* current (fractional) point to be separated */
      double *x; /* double x[1+m+n]; */
      /* x[k] is current value of auxiliary (1 <= k <= m) or structural
         (m+1 <= k <= m+n) variable */
      /*--------------------------------------------------------------*/
      /* aggregated constraint sum a[k] * x[k] = b, which is a linear
         combination of original constraints transformed to equalities
         by introducing auxiliary variables */
      int agg_cnt;
      /* number of rows (original constraints) used to build aggregated
         constraint, 1 <= agg_cnt <= MAXAGGR */
      int *agg_row; /* int agg_row[1+MAXAGGR]; */
      /* agg_row[k], 1 <= k <= agg_cnt, is the row number used to build
         aggregated constraint */
      SPV *agg_vec; /* SPV agg_vec[1:m+n]; */
      /* sparse vector of aggregated constraint coefficients, a[k] */
      double agg_rhs;
      /* right-hand side of the aggregated constraint, b */
      /*--------------------------------------------------------------*/
      /* bound substitution flags for modified constraint */
      char *subst; /* char subst[1+m+n]; */
      /* subst[k], 1 <= k <= m+n, is a bound substitution flag used for
         variable x[k]:
         '?' - x[k] is missing in modified constraint
         'L' - x[k] = (lower bound) + x'[k]
         'U' - x[k] = (upper bound) - x'[k] */
      /*--------------------------------------------------------------*/
      /* modified constraint sum a'[k] * x'[k] = b', where x'[k] >= 0,
         derived from aggregated constraint by substituting bounds;
         note that due to substitution of variable bounds there may be
         additional terms in the modified constraint */
      SPV *mod_vec; /* SPV mod_vec[1:m+n]; */
      /* sparse vector of modified constraint coefficients, a'[k] */
      double mod_rhs;
      /* right-hand side of the modified constraint, b' */
      /*--------------------------------------------------------------*/
      /* cutting plane sum alpha[k] * x[k] <= beta */
      SPV *cut_vec; /* SPV cut_vec[1:m+n]; */
      /* sparse vector of cutting plane coefficients, alpha[k] */
      double cut_rhs;
      /* right-hand size of the cutting plane, beta */
};

/***********************************************************************
*  NAME
*
*  glp_mir_init - create and initialize MIR cut generator
*
*  SYNOPSIS
*
*  glp_mir *glp_mir_init(glp_prob *P);
*
*  DESCRIPTION
*
*  This routine creates and initializes the MIR cut generator for the
*  specified problem object.
*
*  RETURNS
*
*  The routine returns a pointer to the MIR cut generator workspace. */

static void set_row_attrib(glp_prob *mip, glp_mir *mir)
{     /* set global row attributes */
      int m = mir->m;
      int k;
      for (k = 1; k <= m; k++)
      {  GLPROW *row = mip->row[k];
         mir->skip[k] = 0;
         mir->isint[k] = 0;
         switch (row->type)
         {  case GLP_FR:
               mir->lb[k] = -DBL_MAX, mir->ub[k] = +DBL_MAX; break;
            case GLP_LO:
               mir->lb[k] = row->lb, mir->ub[k] = +DBL_MAX; break;
            case GLP_UP:
               mir->lb[k] = -DBL_MAX, mir->ub[k] = row->ub; break;
            case GLP_DB:
               mir->lb[k] = row->lb, mir->ub[k] = row->ub; break;
            case GLP_FX:
               mir->lb[k] = mir->ub[k] = row->lb; break;
            default:
               xassert(row != row);
         }
         mir->vlb[k] = mir->vub[k] = 0;
      }
      return;
}

static void set_col_attrib(glp_prob *mip, glp_mir *mir)
{     /* set global column attributes */
      int m = mir->m;
      int n = mir->n;
      int k;
      for (k = m+1; k <= m+n; k++)
      {  GLPCOL *col = mip->col[k-m];
         switch (col->kind)
         {  case GLP_CV:
               mir->isint[k] = 0; break;
            case GLP_IV:
               mir->isint[k] = 1; break;
            default:
               xassert(col != col);
         }
         switch (col->type)
         {  case GLP_FR:
               mir->lb[k] = -DBL_MAX, mir->ub[k] = +DBL_MAX; break;
            case GLP_LO:
               mir->lb[k] = col->lb, mir->ub[k] = +DBL_MAX; break;
            case GLP_UP:
               mir->lb[k] = -DBL_MAX, mir->ub[k] = col->ub; break;
            case GLP_DB:
               mir->lb[k] = col->lb, mir->ub[k] = col->ub; break;
            case GLP_FX:
               mir->lb[k] = mir->ub[k] = col->lb; break;
            default:
               xassert(col != col);
         }
         mir->vlb[k] = mir->vub[k] = 0;
      }
      return;
}

static void set_var_bounds(glp_prob *mip, glp_mir *mir)
{     /* set variable bounds */
      int m = mir->m;
      GLPAIJ *aij;
      int i, k1, k2;
      double a1, a2;
      for (i = 1; i <= m; i++)
      {  /* we need the row to be '>= 0' or '<= 0' */
         if (!(mir->lb[i] == 0.0 && mir->ub[i] == +DBL_MAX ||
               mir->lb[i] == -DBL_MAX && mir->ub[i] == 0.0)) continue;
         /* take first term */
         aij = mip->row[i]->ptr;
         if (aij == NULL) continue;
         k1 = m + aij->col->j, a1 = aij->val;
         /* take second term */
         aij = aij->r_next;
         if (aij == NULL) continue;
         k2 = m + aij->col->j, a2 = aij->val;
         /* there must be only two terms */
         if (aij->r_next != NULL) continue;
         /* interchange terms, if needed */
         if (!mir->isint[k1] && mir->isint[k2])
            ;
         else if (mir->isint[k1] && !mir->isint[k2])
         {  k2 = k1, a2 = a1;
            k1 = m + aij->col->j, a1 = aij->val;
         }
         else
         {  /* both terms are either continuous or integer */
            continue;
         }
         /* x[k2] should be double-bounded */
         if (mir->lb[k2] == -DBL_MAX || mir->ub[k2] == +DBL_MAX ||
             mir->lb[k2] == mir->ub[k2]) continue;
         /* change signs, if necessary */
         if (mir->ub[i] == 0.0) a1 = - a1, a2 = - a2;
         /* now the row has the form a1 * x1 + a2 * x2 >= 0, where x1
            is continuous, x2 is integer */
         if (a1 > 0.0)
         {  /* x1 >= - (a2 / a1) * x2 */
            if (mir->vlb[k1] == 0)
            {  /* set variable lower bound for x1 */
               mir->lb[k1] = - a2 / a1;
               mir->vlb[k1] = k2;
               /* the row should not be used */
               mir->skip[i] = 1;
            }
         }
         else /* a1 < 0.0 */
         {  /* x1 <= - (a2 / a1) * x2 */
            if (mir->vub[k1] == 0)
            {  /* set variable upper bound for x1 */
               mir->ub[k1] = - a2 / a1;
               mir->vub[k1] = k2;
               /* the row should not be used */
               mir->skip[i] = 1;
            }
         }
      }
      return;
}

static void mark_useless_rows(glp_prob *mip, glp_mir *mir)
{     /* mark rows which should not be used */
      int m = mir->m;
      GLPAIJ *aij;
      int i, k, nv;
      for (i = 1; i <= m; i++)
      {  /* free rows should not be used */
         if (mir->lb[i] == -DBL_MAX && mir->ub[i] == +DBL_MAX)
         {  mir->skip[i] = 1;
            continue;
         }
         nv = 0;
         for (aij = mip->row[i]->ptr; aij != NULL; aij = aij->r_next)
         {  k = m + aij->col->j;
            /* rows with free variables should not be used */
            if (mir->lb[k] == -DBL_MAX && mir->ub[k] == +DBL_MAX)
            {  mir->skip[i] = 1;
               break;
            }
            /* rows with integer variables having infinite (lower or
               upper) bound should not be used */
            if (mir->isint[k] && mir->lb[k] == -DBL_MAX ||
                mir->isint[k] && mir->ub[k] == +DBL_MAX)
            {  mir->skip[i] = 1;
               break;
            }
            /* count non-fixed variables */
            if (!(mir->vlb[k] == 0 && mir->vub[k] == 0 &&
                  mir->lb[k] == mir->ub[k])) nv++;
         }
         /* rows with all variables fixed should not be used */
         if (nv == 0)
         {  mir->skip[i] = 1;
            continue;
         }
      }
      return;
}

glp_mir *glp_mir_init(glp_prob *mip)
{     /* create and initialize MIR cut generator */
      int m = mip->m;
      int n = mip->n;
      glp_mir *mir;
#if MIR_DEBUG
      xprintf("ios_mir_init: warning: debug mode enabled\n");
#endif
      /* allocate working area */
      mir = xmalloc(sizeof(glp_mir));
      mir->m = m;
      mir->n = n;
      mir->skip = xcalloc(1+m, sizeof(char));
      mir->isint = xcalloc(1+m+n, sizeof(char));
      mir->lb = xcalloc(1+m+n, sizeof(double));
      mir->vlb = xcalloc(1+m+n, sizeof(int));
      mir->ub = xcalloc(1+m+n, sizeof(double));
      mir->vub = xcalloc(1+m+n, sizeof(int));
      mir->x = xcalloc(1+m+n, sizeof(double));
      mir->agg_row = xcalloc(1+MAXAGGR, sizeof(int));
      mir->agg_vec = spv_create_vec(m+n);
      mir->subst = xcalloc(1+m+n, sizeof(char));
      mir->mod_vec = spv_create_vec(m+n);
      mir->cut_vec = spv_create_vec(m+n);
      /* set global row attributes */
      set_row_attrib(mip, mir);
      /* set global column attributes */
      set_col_attrib(mip, mir);
      /* set variable bounds */
      set_var_bounds(mip, mir);
      /* mark rows which should not be used */
      mark_useless_rows(mip, mir);
      return mir;
}

/***********************************************************************
*  NAME
*
*  glp_mir_gen - generate mixed integer rounding (MIR) cuts
*
*  SYNOPSIS
*
*  int glp_mir_gen(glp_prob *P, glp_mir *mir, glp_prob *pool);
*
*  DESCRIPTION
*
*  This routine attempts to generate mixed integer rounding (MIR) cuts
*  for current basic solution to the specified problem object.
*
*  The cutting plane inequalities generated by the routine are added to
*  the specified cut pool.
*
*  RETURNS
*
*  The routine returns the number of cuts that have been generated and
*  added to the cut pool. */

static void get_current_point(glp_prob *mip, glp_mir *mir)
{     /* obtain current point */
      int m = mir->m;
      int n = mir->n;
      int k;
      for (k = 1; k <= m; k++)
         mir->x[k] = mip->row[k]->prim;
      for (k = m+1; k <= m+n; k++)
         mir->x[k] = mip->col[k-m]->prim;
      return;
}

#if MIR_DEBUG
static void check_current_point(glp_mir *mir)
{     /* check current point */
      int m = mir->m;
      int n = mir->n;
      int k, kk;
      double lb, ub, eps;
      for (k = 1; k <= m+n; k++)
      {  /* determine lower bound */
         lb = mir->lb[k];
         kk = mir->vlb[k];
         if (kk != 0)
         {  xassert(lb != -DBL_MAX);
            xassert(!mir->isint[k]);
            xassert(mir->isint[kk]);
            lb *= mir->x[kk];
         }
         /* check lower bound */
         if (lb != -DBL_MAX)
         {  eps = 1e-6 * (1.0 + fabs(lb));
            xassert(mir->x[k] >= lb - eps);
         }
         /* determine upper bound */
         ub = mir->ub[k];
         kk = mir->vub[k];
         if (kk != 0)
         {  xassert(ub != +DBL_MAX);
            xassert(!mir->isint[k]);
            xassert(mir->isint[kk]);
            ub *= mir->x[kk];
         }
         /* check upper bound */
         if (ub != +DBL_MAX)
         {  eps = 1e-6 * (1.0 + fabs(ub));
            xassert(mir->x[k] <= ub + eps);
         }
      }
      return;
}
#endif

static void initial_agg_row(glp_prob *mip, glp_mir *mir, int i)
{     /* use original i-th row as initial aggregated constraint */
      int m = mir->m;
      GLPAIJ *aij;
      xassert(1 <= i && i <= m);
      xassert(!mir->skip[i]);
      /* mark i-th row in order not to use it in the same aggregated
         constraint */
      mir->skip[i] = 2;
      mir->agg_cnt = 1;
      mir->agg_row[1] = i;
      /* use x[i] - sum a[i,j] * x[m+j] = 0, where x[i] is auxiliary
         variable of row i, x[m+j] are structural variables */
      spv_clear_vec(mir->agg_vec);
      spv_set_vj(mir->agg_vec, i, 1.0);
      for (aij = mip->row[i]->ptr; aij != NULL; aij = aij->r_next)
         spv_set_vj(mir->agg_vec, m + aij->col->j, - aij->val);
      mir->agg_rhs = 0.0;
#if MIR_DEBUG
      spv_check_vec(mir->agg_vec);
#endif
      return;
}

#if MIR_DEBUG
static void check_agg_row(glp_mir *mir)
{     /* check aggregated constraint */
      int m = mir->m;
      int n = mir->n;
      int j, k;
      double r, big;
      /* compute the residual r = sum a[k] * x[k] - b and determine
         big = max(1, |a[k]|, |b|) */
      r = 0.0, big = 1.0;
      for (j = 1; j <= mir->agg_vec->nnz; j++)
      {  k = mir->agg_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         r += mir->agg_vec->val[j] * mir->x[k];
         if (big < fabs(mir->agg_vec->val[j]))
            big = fabs(mir->agg_vec->val[j]);
      }
      r -= mir->agg_rhs;
      if (big < fabs(mir->agg_rhs))
         big = fabs(mir->agg_rhs);
      /* the residual must be close to zero */
      xassert(fabs(r) <= 1e-6 * big);
      return;
}
#endif

static void subst_fixed_vars(glp_mir *mir)
{     /* substitute fixed variables into aggregated constraint */
      int m = mir->m;
      int n = mir->n;
      int j, k;
      for (j = 1; j <= mir->agg_vec->nnz; j++)
      {  k = mir->agg_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         if (mir->vlb[k] == 0 && mir->vub[k] == 0 &&
             mir->lb[k] == mir->ub[k])
         {  /* x[k] is fixed */
            mir->agg_rhs -= mir->agg_vec->val[j] * mir->lb[k];
            mir->agg_vec->val[j] = 0.0;
         }
      }
      /* remove terms corresponding to fixed variables */
      spv_clean_vec(mir->agg_vec, DBL_EPSILON);
#if MIR_DEBUG
      spv_check_vec(mir->agg_vec);
#endif
      return;
}

static void bound_subst_heur(glp_mir *mir)
{     /* bound substitution heuristic */
      int m = mir->m;
      int n = mir->n;
      int j, k, kk;
      double d1, d2;
      for (j = 1; j <= mir->agg_vec->nnz; j++)
      {  k = mir->agg_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         if (mir->isint[k]) continue; /* skip integer variable */
         /* compute distance from x[k] to its lower bound */
         kk = mir->vlb[k];
         if (kk == 0)
         {  if (mir->lb[k] == -DBL_MAX)
               d1 = DBL_MAX;
            else
               d1 = mir->x[k] - mir->lb[k];
         }
         else
         {  xassert(1 <= kk && kk <= m+n);
            xassert(mir->isint[kk]);
            xassert(mir->lb[k] != -DBL_MAX);
            d1 = mir->x[k] - mir->lb[k] * mir->x[kk];
         }
         /* compute distance from x[k] to its upper bound */
         kk = mir->vub[k];
         if (kk == 0)
         {  if (mir->vub[k] == +DBL_MAX)
               d2 = DBL_MAX;
            else
               d2 = mir->ub[k] - mir->x[k];
         }
         else
         {  xassert(1 <= kk && kk <= m+n);
            xassert(mir->isint[kk]);
            xassert(mir->ub[k] != +DBL_MAX);
            d2 = mir->ub[k] * mir->x[kk] - mir->x[k];
         }
         /* x[k] cannot be free */
         xassert(d1 != DBL_MAX || d2 != DBL_MAX);
         /* choose the bound which is closer to x[k] */
         xassert(mir->subst[k] == '?');
         if (d1 <= d2)
            mir->subst[k] = 'L';
         else
            mir->subst[k] = 'U';
      }
      return;
}

static void build_mod_row(glp_mir *mir)
{     /* substitute bounds and build modified constraint */
      int m = mir->m;
      int n = mir->n;
      int j, jj, k, kk;
      /* initially modified constraint is aggregated constraint */
      spv_copy_vec(mir->mod_vec, mir->agg_vec);
      mir->mod_rhs = mir->agg_rhs;
#if MIR_DEBUG
      spv_check_vec(mir->mod_vec);
#endif
      /* substitute bounds for continuous variables; note that due to
         substitution of variable bounds additional terms may appear in
         modified constraint */
      for (j = mir->mod_vec->nnz; j >= 1; j--)
      {  k = mir->mod_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         if (mir->isint[k]) continue; /* skip integer variable */
         if (mir->subst[k] == 'L')
         {  /* x[k] = (lower bound) + x'[k] */
            xassert(mir->lb[k] != -DBL_MAX);
            kk = mir->vlb[k];
            if (kk == 0)
            {  /* x[k] = lb[k] + x'[k] */
               mir->mod_rhs -= mir->mod_vec->val[j] * mir->lb[k];
            }
            else
            {  /* x[k] = lb[k] * x[kk] + x'[k] */
               xassert(mir->isint[kk]);
               jj = mir->mod_vec->pos[kk];
               if (jj == 0)
               {  spv_set_vj(mir->mod_vec, kk, 1.0);
                  jj = mir->mod_vec->pos[kk];
                  mir->mod_vec->val[jj] = 0.0;
               }
               mir->mod_vec->val[jj] +=
                  mir->mod_vec->val[j] * mir->lb[k];
            }
         }
         else if (mir->subst[k] == 'U')
         {  /* x[k] = (upper bound) - x'[k] */
            xassert(mir->ub[k] != +DBL_MAX);
            kk = mir->vub[k];
            if (kk == 0)
            {  /* x[k] = ub[k] - x'[k] */
               mir->mod_rhs -= mir->mod_vec->val[j] * mir->ub[k];
            }
            else
            {  /* x[k] = ub[k] * x[kk] - x'[k] */
               xassert(mir->isint[kk]);
               jj = mir->mod_vec->pos[kk];
               if (jj == 0)
               {  spv_set_vj(mir->mod_vec, kk, 1.0);
                  jj = mir->mod_vec->pos[kk];
                  mir->mod_vec->val[jj] = 0.0;
               }
               mir->mod_vec->val[jj] +=
                  mir->mod_vec->val[j] * mir->ub[k];
            }
            mir->mod_vec->val[j] = - mir->mod_vec->val[j];
         }
         else
            xassert(k != k);
      }
#if MIR_DEBUG
      spv_check_vec(mir->mod_vec);
#endif
      /* substitute bounds for integer variables */
      for (j = 1; j <= mir->mod_vec->nnz; j++)
      {  k = mir->mod_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         if (!mir->isint[k]) continue; /* skip continuous variable */
         xassert(mir->subst[k] == '?');
         xassert(mir->vlb[k] == 0 && mir->vub[k] == 0);
         xassert(mir->lb[k] != -DBL_MAX && mir->ub[k] != +DBL_MAX);
         if (fabs(mir->lb[k]) <= fabs(mir->ub[k]))
         {  /* x[k] = lb[k] + x'[k] */
            mir->subst[k] = 'L';
            mir->mod_rhs -= mir->mod_vec->val[j] * mir->lb[k];
         }
         else
         {  /* x[k] = ub[k] - x'[k] */
            mir->subst[k] = 'U';
            mir->mod_rhs -= mir->mod_vec->val[j] * mir->ub[k];
            mir->mod_vec->val[j] = - mir->mod_vec->val[j];
         }
      }
#if MIR_DEBUG
      spv_check_vec(mir->mod_vec);
#endif
      return;
}

#if MIR_DEBUG
static void check_mod_row(glp_mir *mir)
{     /* check modified constraint */
      int m = mir->m;
      int n = mir->n;
      int j, k, kk;
      double r, big, x;
      /* compute the residual r = sum a'[k] * x'[k] - b' and determine
         big = max(1, |a[k]|, |b|) */
      r = 0.0, big = 1.0;
      for (j = 1; j <= mir->mod_vec->nnz; j++)
      {  k = mir->mod_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         if (mir->subst[k] == 'L')
         {  /* x'[k] = x[k] - (lower bound) */
            xassert(mir->lb[k] != -DBL_MAX);
            kk = mir->vlb[k];
            if (kk == 0)
               x = mir->x[k] - mir->lb[k];
            else
               x = mir->x[k] - mir->lb[k] * mir->x[kk];
         }
         else if (mir->subst[k] == 'U')
         {  /* x'[k] = (upper bound) - x[k] */
            xassert(mir->ub[k] != +DBL_MAX);
            kk = mir->vub[k];
            if (kk == 0)
               x = mir->ub[k] - mir->x[k];
            else
               x = mir->ub[k] * mir->x[kk] - mir->x[k];
         }
         else
            xassert(k != k);
         r += mir->mod_vec->val[j] * x;
         if (big < fabs(mir->mod_vec->val[j]))
            big = fabs(mir->mod_vec->val[j]);
      }
      r -= mir->mod_rhs;
      if (big < fabs(mir->mod_rhs))
         big = fabs(mir->mod_rhs);
      /* the residual must be close to zero */
      xassert(fabs(r) <= 1e-6 * big);
      return;
}
#endif

/***********************************************************************
*  mir_ineq - construct MIR inequality
*
*  Given the single constraint mixed integer set
*
*                    |N|
*     X = {(x,s) in Z    x R  : sum   a[j] * x[j] <= b + s},
*                    +      +  j in N
*
*  this routine constructs the mixed integer rounding (MIR) inequality
*
*     sum   alpha[j] * x[j] <= beta + gamma * s,
*    j in N
*
*  which is valid for X.
*
*  If the MIR inequality has been successfully constructed, the routine
*  returns zero. Otherwise, if b is close to nearest integer, there may
*  be numeric difficulties due to big coefficients; so in this case the
*  routine returns non-zero. */

static int mir_ineq(const int n, const double a[], const double b,
      double alpha[], double *beta, double *gamma)
{     int j;
      double f, t;
      if (fabs(b - floor(b + .5)) < 0.01)
         return 1;
      f = b - floor(b);
      for (j = 1; j <= n; j++)
      {  t = (a[j] - floor(a[j])) - f;
         if (t <= 0.0)
            alpha[j] = floor(a[j]);
         else
            alpha[j] = floor(a[j]) + t / (1.0 - f);
      }
      *beta = floor(b);
      *gamma = 1.0 / (1.0 - f);
      return 0;
}

/***********************************************************************
*  cmir_ineq - construct c-MIR inequality
*
*  Given the mixed knapsack set
*
*      MK              |N|
*     X   = {(x,s) in Z    x R  : sum   a[j] * x[j] <= b + s,
*                      +      +  j in N
*
*             x[j] <= u[j]},
*
*  a subset C of variables to be complemented, and a divisor delta > 0,
*  this routine constructs the complemented MIR (c-MIR) inequality
*
*     sum   alpha[j] * x[j] <= beta + gamma * s,
*    j in N
*                      MK
*  which is valid for X  .
*
*  If the c-MIR inequality has been successfully constructed, the
*  routine returns zero. Otherwise, if there is a risk of numerical
*  difficulties due to big coefficients (see comments to the routine
*  mir_ineq), the routine cmir_ineq returns non-zero. */

static int cmir_ineq(const int n, const double a[], const double b,
      const double u[], const char cset[], const double delta,
      double alpha[], double *beta, double *gamma)
{     int j;
      double *aa, bb;
      aa = alpha, bb = b;
      for (j = 1; j <= n; j++)
      {  aa[j] = a[j] / delta;
         if (cset[j])
            aa[j] = - aa[j], bb -= a[j] * u[j];
      }
      bb /= delta;
      if (mir_ineq(n, aa, bb, alpha, beta, gamma)) return 1;
      for (j = 1; j <= n; j++)
      {  if (cset[j])
            alpha[j] = - alpha[j], *beta += alpha[j] * u[j];
      }
      *gamma /= delta;
      return 0;
}

/***********************************************************************
*  cmir_sep - c-MIR separation heuristic
*
*  Given the mixed knapsack set
*
*      MK              |N|
*     X   = {(x,s) in Z    x R  : sum   a[j] * x[j] <= b + s,
*                      +      +  j in N
*
*             x[j] <= u[j]}
*
*                           *   *
*  and a fractional point (x , s ), this routine tries to construct
*  c-MIR inequality
*
*     sum   alpha[j] * x[j] <= beta + gamma * s,
*    j in N
*                      MK
*  which is valid for X   and has (desirably maximal) violation at the
*  fractional point given. This is attained by choosing an appropriate
*  set C of variables to be complemented and a divisor delta > 0, which
*  together define corresponding c-MIR inequality.
*
*  If a violated c-MIR inequality has been successfully constructed,
*  the routine returns its violation:
*
*                       *                      *
*     sum   alpha[j] * x [j] - beta - gamma * s ,
*    j in N
*
*  which is positive. In case of failure the routine returns zero. */

struct vset { int j; double v; };

static int CDECL cmir_cmp(const void *p1, const void *p2)
{     const struct vset *v1 = p1, *v2 = p2;
      if (v1->v < v2->v) return -1;
      if (v1->v > v2->v) return +1;
      return 0;
}

static double cmir_sep(const int n, const double a[], const double b,
      const double u[], const double x[], const double s,
      double alpha[], double *beta, double *gamma)
{     int fail, j, k, nv, v;
      double delta, eps, d_try[1+3], r, r_best;
      char *cset;
      struct vset *vset;
      /* allocate working arrays */
      cset = xcalloc(1+n, sizeof(char));
      vset = xcalloc(1+n, sizeof(struct vset));
      /* choose initial C */
      for (j = 1; j <= n; j++)
         cset[j] = (char)(x[j] >= 0.5 * u[j]);
      /* choose initial delta */
      r_best = delta = 0.0;
      for (j = 1; j <= n; j++)
      {  xassert(a[j] != 0.0);
         /* if x[j] is close to its bounds, skip it */
         eps = 1e-9 * (1.0 + fabs(u[j]));
         if (x[j] < eps || x[j] > u[j] - eps) continue;
         /* try delta = |a[j]| to construct c-MIR inequality */
         fail = cmir_ineq(n, a, b, u, cset, fabs(a[j]), alpha, beta,
            gamma);
         if (fail) continue;
         /* compute violation */
         r = - (*beta) - (*gamma) * s;
         for (k = 1; k <= n; k++) r += alpha[k] * x[k];
         if (r_best < r) r_best = r, delta = fabs(a[j]);
      }
      if (r_best < 0.001) r_best = 0.0;
      if (r_best == 0.0) goto done;
      xassert(delta > 0.0);
      /* try to increase violation by dividing delta by 2, 4, and 8,
         respectively */
      d_try[1] = delta / 2.0;
      d_try[2] = delta / 4.0;
      d_try[3] = delta / 8.0;
      for (j = 1; j <= 3; j++)
      {  /* construct c-MIR inequality */
         fail = cmir_ineq(n, a, b, u, cset, d_try[j], alpha, beta,
            gamma);
         if (fail) continue;
         /* compute violation */
         r = - (*beta) - (*gamma) * s;
         for (k = 1; k <= n; k++) r += alpha[k] * x[k];
         if (r_best < r) r_best = r, delta = d_try[j];
      }
      /* build subset of variables lying strictly between their bounds
         and order it by nondecreasing values of |x[j] - u[j]/2| */
      nv = 0;
      for (j = 1; j <= n; j++)
      {  /* if x[j] is close to its bounds, skip it */
         eps = 1e-9 * (1.0 + fabs(u[j]));
         if (x[j] < eps || x[j] > u[j] - eps) continue;
         /* add x[j] to the subset */
         nv++;
         vset[nv].j = j;
         vset[nv].v = fabs(x[j] - 0.5 * u[j]);
      }
      qsort(&vset[1], nv, sizeof(struct vset), cmir_cmp);
      /* try to increase violation by successively complementing each
         variable in the subset */
      for (v = 1; v <= nv; v++)
      {  j = vset[v].j;
         /* replace x[j] by its complement or vice versa */
         cset[j] = (char)!cset[j];
         /* construct c-MIR inequality */
         fail = cmir_ineq(n, a, b, u, cset, delta, alpha, beta, gamma);
         /* restore the variable */
         cset[j] = (char)!cset[j];
         /* do not replace the variable in case of failure */
         if (fail) continue;
         /* compute violation */
         r = - (*beta) - (*gamma) * s;
         for (k = 1; k <= n; k++) r += alpha[k] * x[k];
         if (r_best < r) r_best = r, cset[j] = (char)!cset[j];
      }
      /* construct the best c-MIR inequality chosen */
      fail = cmir_ineq(n, a, b, u, cset, delta, alpha, beta, gamma);
      xassert(!fail);
done: /* free working arrays */
      xfree(cset);
      xfree(vset);
      /* return to the calling routine */
      return r_best;
}

static double generate(glp_mir *mir)
{     /* try to generate violated c-MIR cut for modified constraint */
      int m = mir->m;
      int n = mir->n;
      int j, k, kk, nint;
      double s, *u, *x, *alpha, r_best = 0.0, b, beta, gamma;
      spv_copy_vec(mir->cut_vec, mir->mod_vec);
      mir->cut_rhs = mir->mod_rhs;
      /* remove small terms, which can appear due to substitution of
         variable bounds */
      spv_clean_vec(mir->cut_vec, DBL_EPSILON);
#if MIR_DEBUG
      spv_check_vec(mir->cut_vec);
#endif
      /* remove positive continuous terms to obtain MK relaxation */
      for (j = 1; j <= mir->cut_vec->nnz; j++)
      {  k = mir->cut_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         if (!mir->isint[k] && mir->cut_vec->val[j] > 0.0)
            mir->cut_vec->val[j] = 0.0;
      }
      spv_clean_vec(mir->cut_vec, 0.0);
#if MIR_DEBUG
      spv_check_vec(mir->cut_vec);
#endif
      /* move integer terms to the beginning of the sparse vector and
         determine the number of integer variables */
      nint = 0;
      for (j = 1; j <= mir->cut_vec->nnz; j++)
      {  k = mir->cut_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         if (mir->isint[k])
         {  double temp;
            nint++;
            /* interchange elements [nint] and [j] */
            kk = mir->cut_vec->ind[nint];
            mir->cut_vec->pos[k] = nint;
            mir->cut_vec->pos[kk] = j;
            mir->cut_vec->ind[nint] = k;
            mir->cut_vec->ind[j] = kk;
            temp = mir->cut_vec->val[nint];
            mir->cut_vec->val[nint] = mir->cut_vec->val[j];
            mir->cut_vec->val[j] = temp;
         }
      }
#if MIR_DEBUG
      spv_check_vec(mir->cut_vec);
#endif
      /* if there is no integer variable, nothing to generate */
      if (nint == 0) goto done;
      /* allocate working arrays */
      u = xcalloc(1+nint, sizeof(double));
      x = xcalloc(1+nint, sizeof(double));
      alpha = xcalloc(1+nint, sizeof(double));
      /* determine u and x */
      for (j = 1; j <= nint; j++)
      {  k = mir->cut_vec->ind[j];
         xassert(m+1 <= k && k <= m+n);
         xassert(mir->isint[k]);
         u[j] = mir->ub[k] - mir->lb[k];
         xassert(u[j] >= 1.0);
         if (mir->subst[k] == 'L')
            x[j] = mir->x[k] - mir->lb[k];
         else if (mir->subst[k] == 'U')
            x[j] = mir->ub[k] - mir->x[k];
         else
            xassert(k != k);
#if 0 /* 06/III-2016; notorious bug reported many times */
         xassert(x[j] >= -0.001);
#else
         if (x[j] < -0.001)
         {  xprintf("glp_mir_gen: warning: x[%d] = %g\n", j, x[j]);
            r_best = 0.0;
            goto skip;
         }
#endif
         if (x[j] < 0.0) x[j] = 0.0;
      }
      /* compute s = - sum of continuous terms */
      s = 0.0;
      for (j = nint+1; j <= mir->cut_vec->nnz; j++)
      {  double x;
         k = mir->cut_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         /* must be continuous */
         xassert(!mir->isint[k]);
         if (mir->subst[k] == 'L')
         {  xassert(mir->lb[k] != -DBL_MAX);
            kk = mir->vlb[k];
            if (kk == 0)
               x = mir->x[k] - mir->lb[k];
            else
               x = mir->x[k] - mir->lb[k] * mir->x[kk];
         }
         else if (mir->subst[k] == 'U')
         {  xassert(mir->ub[k] != +DBL_MAX);
            kk = mir->vub[k];
            if (kk == 0)
               x = mir->ub[k] - mir->x[k];
            else
               x = mir->ub[k] * mir->x[kk] - mir->x[k];
         }
         else
            xassert(k != k);
#if 0 /* 06/III-2016; notorious bug reported many times */
         xassert(x >= -0.001);
#else
         if (x < -0.001)
         {  xprintf("glp_mir_gen: warning: x = %g\n", x);
            r_best = 0.0;
            goto skip;
         }
#endif
         if (x < 0.0) x = 0.0;
         s -= mir->cut_vec->val[j] * x;
      }
      xassert(s >= 0.0);
      /* apply heuristic to obtain most violated c-MIR inequality */
      b = mir->cut_rhs;
      r_best = cmir_sep(nint, mir->cut_vec->val, b, u, x, s, alpha,
         &beta, &gamma);
      if (r_best == 0.0) goto skip;
      xassert(r_best > 0.0);
      /* convert to raw cut */
      /* sum alpha[j] * x[j] <= beta + gamma * s */
      for (j = 1; j <= nint; j++)
         mir->cut_vec->val[j] = alpha[j];
      for (j = nint+1; j <= mir->cut_vec->nnz; j++)
      {  k = mir->cut_vec->ind[j];
         if (k <= m+n) mir->cut_vec->val[j] *= gamma;
      }
      mir->cut_rhs = beta;
#if MIR_DEBUG
      spv_check_vec(mir->cut_vec);
#endif
skip: /* free working arrays */
      xfree(u);
      xfree(x);
      xfree(alpha);
done: return r_best;
}

#if MIR_DEBUG
static void check_raw_cut(glp_mir *mir, double r_best)
{     /* check raw cut before back bound substitution */
      int m = mir->m;
      int n = mir->n;
      int j, k, kk;
      double r, big, x;
      /* compute the residual r = sum a[k] * x[k] - b and determine
         big = max(1, |a[k]|, |b|) */
      r = 0.0, big = 1.0;
      for (j = 1; j <= mir->cut_vec->nnz; j++)
      {  k = mir->cut_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         if (mir->subst[k] == 'L')
         {  xassert(mir->lb[k] != -DBL_MAX);
            kk = mir->vlb[k];
            if (kk == 0)
               x = mir->x[k] - mir->lb[k];
            else
               x = mir->x[k] - mir->lb[k] * mir->x[kk];
         }
         else if (mir->subst[k] == 'U')
         {  xassert(mir->ub[k] != +DBL_MAX);
            kk = mir->vub[k];
            if (kk == 0)
               x = mir->ub[k] - mir->x[k];
            else
               x = mir->ub[k] * mir->x[kk] - mir->x[k];
         }
         else
            xassert(k != k);
         r += mir->cut_vec->val[j] * x;
         if (big < fabs(mir->cut_vec->val[j]))
            big = fabs(mir->cut_vec->val[j]);
      }
      r -= mir->cut_rhs;
      if (big < fabs(mir->cut_rhs))
         big = fabs(mir->cut_rhs);
      /* the residual must be close to r_best */
      xassert(fabs(r - r_best) <= 1e-6 * big);
      return;
}
#endif

static void back_subst(glp_mir *mir)
{     /* back substitution of original bounds */
      int m = mir->m;
      int n = mir->n;
      int j, jj, k, kk;
      /* at first, restore bounds of integer variables (because on
         restoring variable bounds of continuous variables we need
         original, not shifted, bounds of integer variables) */
      for (j = 1; j <= mir->cut_vec->nnz; j++)
      {  k = mir->cut_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         if (!mir->isint[k]) continue; /* skip continuous */
         if (mir->subst[k] == 'L')
         {  /* x'[k] = x[k] - lb[k] */
            xassert(mir->lb[k] != -DBL_MAX);
            xassert(mir->vlb[k] == 0);
            mir->cut_rhs += mir->cut_vec->val[j] * mir->lb[k];
         }
         else if (mir->subst[k] == 'U')
         {  /* x'[k] = ub[k] - x[k] */
            xassert(mir->ub[k] != +DBL_MAX);
            xassert(mir->vub[k] == 0);
            mir->cut_rhs -= mir->cut_vec->val[j] * mir->ub[k];
            mir->cut_vec->val[j] = - mir->cut_vec->val[j];
         }
         else
            xassert(k != k);
      }
      /* now restore bounds of continuous variables */
      for (j = 1; j <= mir->cut_vec->nnz; j++)
      {  k = mir->cut_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         if (mir->isint[k]) continue; /* skip integer */
         if (mir->subst[k] == 'L')
         {  /* x'[k] = x[k] - (lower bound) */
            xassert(mir->lb[k] != -DBL_MAX);
            kk = mir->vlb[k];
            if (kk == 0)
            {  /* x'[k] = x[k] - lb[k] */
               mir->cut_rhs += mir->cut_vec->val[j] * mir->lb[k];
            }
            else
            {  /* x'[k] = x[k] - lb[k] * x[kk] */
               jj = mir->cut_vec->pos[kk];
#if 0
               xassert(jj != 0);
#else
               if (jj == 0)
               {  spv_set_vj(mir->cut_vec, kk, 1.0);
                  jj = mir->cut_vec->pos[kk];
                  xassert(jj != 0);
                  mir->cut_vec->val[jj] = 0.0;
               }
#endif
               mir->cut_vec->val[jj] -= mir->cut_vec->val[j] *
                  mir->lb[k];
            }
         }
         else if (mir->subst[k] == 'U')
         {  /* x'[k] = (upper bound) - x[k] */
            xassert(mir->ub[k] != +DBL_MAX);
            kk = mir->vub[k];
            if (kk == 0)
            {  /* x'[k] = ub[k] - x[k] */
               mir->cut_rhs -= mir->cut_vec->val[j] * mir->ub[k];
            }
            else
            {  /* x'[k] = ub[k] * x[kk] - x[k] */
               jj = mir->cut_vec->pos[kk];
               if (jj == 0)
               {  spv_set_vj(mir->cut_vec, kk, 1.0);
                  jj = mir->cut_vec->pos[kk];
                  xassert(jj != 0);
                  mir->cut_vec->val[jj] = 0.0;
               }
               mir->cut_vec->val[jj] += mir->cut_vec->val[j] *
                  mir->ub[k];
            }
            mir->cut_vec->val[j] = - mir->cut_vec->val[j];
         }
         else
            xassert(k != k);
      }
#if MIR_DEBUG
      spv_check_vec(mir->cut_vec);
#endif
      return;
}

#if MIR_DEBUG
static void check_cut_row(glp_mir *mir, double r_best)
{     /* check the cut after back bound substitution or elimination of
         auxiliary variables */
      int m = mir->m;
      int n = mir->n;
      int j, k;
      double r, big;
      /* compute the residual r = sum a[k] * x[k] - b and determine
         big = max(1, |a[k]|, |b|) */
      r = 0.0, big = 1.0;
      for (j = 1; j <= mir->cut_vec->nnz; j++)
      {  k = mir->cut_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         r += mir->cut_vec->val[j] * mir->x[k];
         if (big < fabs(mir->cut_vec->val[j]))
            big = fabs(mir->cut_vec->val[j]);
      }
      r -= mir->cut_rhs;
      if (big < fabs(mir->cut_rhs))
         big = fabs(mir->cut_rhs);
      /* the residual must be close to r_best */
      xassert(fabs(r - r_best) <= 1e-6 * big);
      return;
}
#endif

static void subst_aux_vars(glp_prob *mip, glp_mir *mir)
{     /* final substitution to eliminate auxiliary variables */
      int m = mir->m;
      int n = mir->n;
      GLPAIJ *aij;
      int j, k, kk, jj;
      for (j = mir->cut_vec->nnz; j >= 1; j--)
      {  k = mir->cut_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         if (k > m) continue; /* skip structurals */
         for (aij = mip->row[k]->ptr; aij != NULL; aij = aij->r_next)
         {  kk = m + aij->col->j; /* structural */
            jj = mir->cut_vec->pos[kk];
            if (jj == 0)
            {  spv_set_vj(mir->cut_vec, kk, 1.0);
               jj = mir->cut_vec->pos[kk];
               mir->cut_vec->val[jj] = 0.0;
            }
            mir->cut_vec->val[jj] += mir->cut_vec->val[j] * aij->val;
         }
         mir->cut_vec->val[j] = 0.0;
      }
      spv_clean_vec(mir->cut_vec, 0.0);
      return;
}

static void add_cut(glp_mir *mir, glp_prob *pool)
{     /* add constructed cut inequality to the cut pool */
      int m = mir->m;
      int n = mir->n;
      int j, k, len;
      int *ind = xcalloc(1+n, sizeof(int));
      double *val = xcalloc(1+n, sizeof(double));
      len = 0;
      for (j = mir->cut_vec->nnz; j >= 1; j--)
      {  k = mir->cut_vec->ind[j];
         xassert(m+1 <= k && k <= m+n);
         len++, ind[len] = k - m, val[len] = mir->cut_vec->val[j];
      }
#if 0
#if 0
      ios_add_cut_row(tree, pool, GLP_RF_MIR, len, ind, val, GLP_UP,
         mir->cut_rhs);
#else
      glp_ios_add_row(tree, NULL, GLP_RF_MIR, 0, len, ind, val, GLP_UP,
         mir->cut_rhs);
#endif
#else
      {  int i;
         i = glp_add_rows(pool, 1);
         glp_set_row_bnds(pool, i, GLP_UP, 0, mir->cut_rhs);
         glp_set_mat_row(pool, i, len, ind, val);
      }
#endif
      xfree(ind);
      xfree(val);
      return;
}

#if 0 /* 29/II-2016 by Chris */
static int aggregate_row(glp_prob *mip, glp_mir *mir)
#else
static int aggregate_row(glp_prob *mip, glp_mir *mir, SPV *v)
#endif
{     /* try to aggregate another row */
      int m = mir->m;
      int n = mir->n;
      GLPAIJ *aij;
#if 0 /* 29/II-2016 by Chris */
      SPV *v;
#endif
      int ii, j, jj, k, kk, kappa = 0, ret = 0;
      double d1, d2, d, d_max = 0.0;
      /* choose appropriate structural variable in the aggregated row
         to be substituted */
      for (j = 1; j <= mir->agg_vec->nnz; j++)
      {  k = mir->agg_vec->ind[j];
         xassert(1 <= k && k <= m+n);
         if (k <= m) continue; /* skip auxiliary var */
         if (mir->isint[k]) continue; /* skip integer var */
         if (fabs(mir->agg_vec->val[j]) < 0.001) continue;
         /* compute distance from x[k] to its lower bound */
         kk = mir->vlb[k];
         if (kk == 0)
         {  if (mir->lb[k] == -DBL_MAX)
               d1 = DBL_MAX;
            else
               d1 = mir->x[k] - mir->lb[k];
         }
         else
         {  xassert(1 <= kk && kk <= m+n);
            xassert(mir->isint[kk]);
            xassert(mir->lb[k] != -DBL_MAX);
            d1 = mir->x[k] - mir->lb[k] * mir->x[kk];
         }
         /* compute distance from x[k] to its upper bound */
         kk = mir->vub[k];
         if (kk == 0)
         {  if (mir->vub[k] == +DBL_MAX)
               d2 = DBL_MAX;
            else
               d2 = mir->ub[k] - mir->x[k];
         }
         else
         {  xassert(1 <= kk && kk <= m+n);
            xassert(mir->isint[kk]);
            xassert(mir->ub[k] != +DBL_MAX);
            d2 = mir->ub[k] * mir->x[kk] - mir->x[k];
         }
         /* x[k] cannot be free */
         xassert(d1 != DBL_MAX || d2 != DBL_MAX);
         /* d = min(d1, d2) */
         d = (d1 <= d2 ? d1 : d2);
         xassert(d != DBL_MAX);
         /* should not be close to corresponding bound */
         if (d < 0.001) continue;
         if (d_max < d) d_max = d, kappa = k;
      }
      if (kappa == 0)
      {  /* nothing chosen */
         ret = 1;
         goto done;
      }
      /* x[kappa] has been chosen */
      xassert(m+1 <= kappa && kappa <= m+n);
      xassert(!mir->isint[kappa]);
      /* find another row, which have not been used yet, to eliminate
         x[kappa] from the aggregated row */
#if 0 /* 29/II-2016 by Chris */
      for (ii = 1; ii <= m; ii++)
      {  if (mir->skip[ii]) continue;
         for (aij = mip->row[ii]->ptr; aij != NULL; aij = aij->r_next)
            if (aij->col->j == kappa - m) break;
         if (aij != NULL && fabs(aij->val) >= 0.001) break;
#else
      ii = 0;
      for (aij = mip->col[kappa - m]->ptr; aij != NULL;
         aij = aij->c_next)
      {  if (aij->row->i > m) continue;
         if (mir->skip[aij->row->i]) continue;
         if (fabs(aij->val) >= 0.001)
         {  ii = aij->row->i;
            break;
         }
#endif
      }
#if 0 /* 29/II-2016 by Chris */
      if (ii > m)
#else
      if (ii == 0)
#endif
      {  /* nothing found */
         ret = 2;
         goto done;
      }
      /* row ii has been found; include it in the aggregated list */
      mir->agg_cnt++;
      xassert(mir->agg_cnt <= MAXAGGR);
      mir->agg_row[mir->agg_cnt] = ii;
      mir->skip[ii] = 2;
      /* v := new row */
#if 0 /* 29/II-2016 by Chris */
      v = ios_create_vec(m+n);
#else
      spv_clear_vec(v);
#endif
      spv_set_vj(v, ii, 1.0);
      for (aij = mip->row[ii]->ptr; aij != NULL; aij = aij->r_next)
         spv_set_vj(v, m + aij->col->j, - aij->val);
#if MIR_DEBUG
      spv_check_vec(v);
#endif
      /* perform gaussian elimination to remove x[kappa] */
      j = mir->agg_vec->pos[kappa];
      xassert(j != 0);
      jj = v->pos[kappa];
      xassert(jj != 0);
      spv_linear_comb(mir->agg_vec,
         - mir->agg_vec->val[j] / v->val[jj], v);
#if 0 /* 29/II-2016 by Chris */
      ios_delete_vec(v);
#endif
      spv_set_vj(mir->agg_vec, kappa, 0.0);
#if MIR_DEBUG
      spv_check_vec(mir->agg_vec);
#endif
done: return ret;
}

int glp_mir_gen(glp_prob *mip, glp_mir *mir, glp_prob *pool)
{     /* main routine to generate MIR cuts */
      int m = mir->m;
      int n = mir->n;
      int i, nnn = 0;
      double r_best;
#if 1 /* 29/II-2016 by Chris */
      SPV *work;
#endif
      xassert(mip->m >= m);
      xassert(mip->n == n);
      /* obtain current point */
      get_current_point(mip, mir);
#if MIR_DEBUG
      /* check current point */
      check_current_point(mir);
#endif
      /* reset bound substitution flags */
      memset(&mir->subst[1], '?', m+n);
#if 1 /* 29/II-2016 by Chris */
      work = spv_create_vec(m+n);
#endif
      /* try to generate a set of violated MIR cuts */
      for (i = 1; i <= m; i++)
      {  if (mir->skip[i]) continue;
         /* use original i-th row as initial aggregated constraint */
         initial_agg_row(mip, mir, i);
loop:    ;
#if MIR_DEBUG
         /* check aggregated row */
         check_agg_row(mir);
#endif
         /* substitute fixed variables into aggregated constraint */
         subst_fixed_vars(mir);
#if MIR_DEBUG
         /* check aggregated row */
         check_agg_row(mir);
#endif
#if MIR_DEBUG
         /* check bound substitution flags */
         {  int k;
            for (k = 1; k <= m+n; k++)
               xassert(mir->subst[k] == '?');
         }
#endif
         /* apply bound substitution heuristic */
         bound_subst_heur(mir);
         /* substitute bounds and build modified constraint */
         build_mod_row(mir);
#if MIR_DEBUG
         /* check modified row */
         check_mod_row(mir);
#endif
         /* try to generate violated c-MIR cut for modified row */
         r_best = generate(mir);
         if (r_best > 0.0)
         {  /* success */
#if MIR_DEBUG
            /* check raw cut before back bound substitution */
            check_raw_cut(mir, r_best);
#endif
            /* back substitution of original bounds */
            back_subst(mir);
#if MIR_DEBUG
            /* check the cut after back bound substitution */
            check_cut_row(mir, r_best);
#endif
            /* final substitution to eliminate auxiliary variables */
            subst_aux_vars(mip, mir);
#if MIR_DEBUG
            /* check the cut after elimination of auxiliaries */
            check_cut_row(mir, r_best);
#endif
            /* add constructed cut inequality to the cut pool */
            add_cut(mir, pool), nnn++;
         }
         /* reset bound substitution flags */
         {  int j, k;
            for (j = 1; j <= mir->mod_vec->nnz; j++)
            {  k = mir->mod_vec->ind[j];
               xassert(1 <= k && k <= m+n);
               xassert(mir->subst[k] != '?');
               mir->subst[k] = '?';
            }
         }
         if (r_best == 0.0)
         {  /* failure */
            if (mir->agg_cnt < MAXAGGR)
            {  /* try to aggregate another row */
#if 0 /* 29/II-2016 by Chris */
               if (aggregate_row(mip, mir) == 0) goto loop;
#else
               if (aggregate_row(mip, mir, work) == 0) goto loop;
#endif
            }
         }
         /* unmark rows used in the aggregated constraint */
         {  int k, ii;
            for (k = 1; k <= mir->agg_cnt; k++)
            {  ii = mir->agg_row[k];
               xassert(1 <= ii && ii <= m);
               xassert(mir->skip[ii] == 2);
               mir->skip[ii] = 0;
            }
         }
      }
#if 1 /* 29/II-2016 by Chris */
      spv_delete_vec(work);
#endif
      return nnn;
}

/***********************************************************************
*  NAME
*
*  glp_mir_free - delete MIR cut generator workspace
*
*  SYNOPSIS
*
*  void glp_mir_free(glp_mir *mir);
*
*  DESCRIPTION
*
*  This routine deletes the MIR cut generator workspace and frees all
*  the memory allocated to it. */

void glp_mir_free(glp_mir *mir)
{     xfree(mir->skip);
      xfree(mir->isint);
      xfree(mir->lb);
      xfree(mir->vlb);
      xfree(mir->ub);
      xfree(mir->vub);
      xfree(mir->x);
      xfree(mir->agg_row);
      spv_delete_vec(mir->agg_vec);
      xfree(mir->subst);
      spv_delete_vec(mir->mod_vec);
      spv_delete_vec(mir->cut_vec);
      xfree(mir);
      return;
}

/* eof */
