/* glpios05.c (Gomory's mixed integer cut generator) */

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
*  NAME
*
*  ios_gmi_gen - generate Gomory's mixed integer cuts.
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_gmi_gen(glp_tree *tree, IOSPOOL *pool);
*
*  DESCRIPTION
*
*  The routine ios_gmi_gen generates Gomory's mixed integer cuts for
*  the current point and adds them to the cut pool. */

#define MAXCUTS 50
/* maximal number of cuts to be generated for one round */

struct worka
{     /* Gomory's cut generator working area */
      int *ind; /* int ind[1+n]; */
      double *val; /* double val[1+n]; */
      double *phi; /* double phi[1+m+n]; */
};

#define f(x) ((x) - floor(x))
/* compute fractional part of x */

static void gen_cut(glp_tree *tree, struct worka *worka, int j)
{     /* this routine tries to generate Gomory's mixed integer cut for
         specified structural variable x[m+j] of integer kind, which is
         basic and has fractional value in optimal solution to current
         LP relaxation */
      glp_prob *mip = tree->mip;
      int m = mip->m;
      int n = mip->n;
      int *ind = worka->ind;
      double *val = worka->val;
      double *phi = worka->phi;
      int i, k, len, kind, stat;
      double lb, ub, alfa, beta, ksi, phi1, rhs;
      /* compute row of the simplex tableau, which (row) corresponds
         to specified basic variable xB[i] = x[m+j]; see (23) */
      len = glp_eval_tab_row(mip, m+j, ind, val);
      /* determine beta[i], which a value of xB[i] in optimal solution
         to current LP relaxation; note that this value is the same as
         if it would be computed with formula (27); it is assumed that
         beta[i] is fractional enough */
      beta = mip->col[j]->prim;
      /* compute cut coefficients phi and right-hand side rho, which
         correspond to formula (30); dense format is used, because rows
         of the simplex tableau is usually dense */
      for (k = 1; k <= m+n; k++) phi[k] = 0.0;
      rhs = f(beta); /* initial value of rho; see (28), (32) */
      for (j = 1; j <= len; j++)
      {  /* determine original number of non-basic variable xN[j] */
         k = ind[j];
         xassert(1 <= k && k <= m+n);
         /* determine the kind, bounds and current status of xN[j] in
            optimal solution to LP relaxation */
         if (k <= m)
         {  /* auxiliary variable */
            GLPROW *row = mip->row[k];
            kind = GLP_CV;
            lb = row->lb;
            ub = row->ub;
            stat = row->stat;
         }
         else
         {  /* structural variable */
            GLPCOL *col = mip->col[k-m];
            kind = col->kind;
            lb = col->lb;
            ub = col->ub;
            stat = col->stat;
         }
         /* xN[j] cannot be basic */
         xassert(stat != GLP_BS);
         /* determine row coefficient ksi[i,j] at xN[j]; see (23) */
         ksi = val[j];
         /* if ksi[i,j] is too large in the magnitude, do not generate
            the cut */
         if (fabs(ksi) > 1e+05) goto fini;
         /* if ksi[i,j] is too small in the magnitude, skip it */
         if (fabs(ksi) < 1e-10) goto skip;
         /* compute row coefficient alfa[i,j] at y[j]; see (26) */
         switch (stat)
         {  case GLP_NF:
               /* xN[j] is free (unbounded) having non-zero ksi[i,j];
                  do not generate the cut */
               goto fini;
            case GLP_NL:
               /* xN[j] has active lower bound */
               alfa = - ksi;
               break;
            case GLP_NU:
               /* xN[j] has active upper bound */
               alfa = + ksi;
               break;
            case GLP_NS:
               /* xN[j] is fixed; skip it */
               goto skip;
            default:
               xassert(stat != stat);
         }
         /* compute cut coefficient phi'[j] at y[j]; see (21), (28) */
         switch (kind)
         {  case GLP_IV:
               /* y[j] is integer */
               if (fabs(alfa - floor(alfa + 0.5)) < 1e-10)
               {  /* alfa[i,j] is close to nearest integer; skip it */
                  goto skip;
               }
               else if (f(alfa) <= f(beta))
                  phi1 = f(alfa);
               else
                  phi1 = (f(beta) / (1.0 - f(beta))) * (1.0 - f(alfa));
               break;
            case GLP_CV:
               /* y[j] is continuous */
               if (alfa >= 0.0)
                  phi1 = + alfa;
               else
                  phi1 = (f(beta) / (1.0 - f(beta))) * (- alfa);
               break;
            default:
               xassert(kind != kind);
         }
         /* compute cut coefficient phi[j] at xN[j] and update right-
            hand side rho; see (31), (32) */
         switch (stat)
         {  case GLP_NL:
               /* xN[j] has active lower bound */
               phi[k] = + phi1;
               rhs += phi1 * lb;
               break;
            case GLP_NU:
               /* xN[j] has active upper bound */
               phi[k] = - phi1;
               rhs -= phi1 * ub;
               break;
            default:
               xassert(stat != stat);
         }
skip:    ;
      }
      /* now the cut has the form sum_k phi[k] * x[k] >= rho, where cut
         coefficients are stored in the array phi in dense format;
         x[1,...,m] are auxiliary variables, x[m+1,...,m+n] are struc-
         tural variables; see (30) */
      /* eliminate auxiliary variables in order to express the cut only
         through structural variables; see (33) */
      for (i = 1; i <= m; i++)
      {  GLPROW *row;
         GLPAIJ *aij;
         if (fabs(phi[i]) < 1e-10) continue;
         /* auxiliary variable x[i] has non-zero cut coefficient */
         row = mip->row[i];
         /* x[i] cannot be fixed */
         xassert(row->type != GLP_FX);
         /* substitute x[i] = sum_j a[i,j] * x[m+j] */
         for (aij = row->ptr; aij != NULL; aij = aij->r_next)
            phi[m+aij->col->j] += phi[i] * aij->val;
      }
      /* convert the final cut to sparse format and substitute fixed
         (structural) variables */
      len = 0;
      for (j = 1; j <= n; j++)
      {  GLPCOL *col;
         if (fabs(phi[m+j]) < 1e-10) continue;
         /* structural variable x[m+j] has non-zero cut coefficient */
         col = mip->col[j];
         if (col->type == GLP_FX)
         {  /* eliminate x[m+j] */
            rhs -= phi[m+j] * col->lb;
         }
         else
         {  len++;
            ind[len] = j;
            val[len] = phi[m+j];
         }
      }
      if (fabs(rhs) < 1e-12) rhs = 0.0;
      /* if the cut inequality seems to be badly scaled, reject it to
         avoid numeric difficulties */
      for (k = 1; k <= len; k++)
      {  if (fabs(val[k]) < 1e-03) goto fini;
         if (fabs(val[k]) > 1e+03) goto fini;
      }
      /* add the cut to the cut pool for further consideration */
#if 0
      ios_add_cut_row(tree, pool, GLP_RF_GMI, len, ind, val, GLP_LO,
         rhs);
#else
      glp_ios_add_row(tree, NULL, GLP_RF_GMI, 0, len, ind, val, GLP_LO,
         rhs);
#endif
fini: return;
}

struct var { int j; double f; };

static int fcmp(const void *p1, const void *p2)
{     const struct var *v1 = p1, *v2 = p2;
      if (v1->f > v2->f) return -1;
      if (v1->f < v2->f) return +1;
      return 0;
}

void ios_gmi_gen(glp_tree *tree)
{     /* main routine to generate Gomory's cuts */
      glp_prob *mip = tree->mip;
      int m = mip->m;
      int n = mip->n;
      struct var *var;
      int k, nv, j, size;
      struct worka _worka, *worka = &_worka;
      /* allocate working arrays */
      var = xcalloc(1+n, sizeof(struct var));
      worka->ind = xcalloc(1+n, sizeof(int));
      worka->val = xcalloc(1+n, sizeof(double));
      worka->phi = xcalloc(1+m+n, sizeof(double));
      /* build the list of integer structural variables, which are
         basic and have fractional value in optimal solution to current
         LP relaxation */
      nv = 0;
      for (j = 1; j <= n; j++)
      {  GLPCOL *col = mip->col[j];
         double frac;
         if (col->kind != GLP_IV) continue;
         if (col->type == GLP_FX) continue;
         if (col->stat != GLP_BS) continue;
         frac = f(col->prim);
         if (!(0.05 <= frac && frac <= 0.95)) continue;
         /* add variable to the list */
         nv++, var[nv].j = j, var[nv].f = frac;
      }
      /* order the list by descending fractionality */
      qsort(&var[1], nv, sizeof(struct var), fcmp);
      /* try to generate cuts by one for each variable in the list, but
         not more than MAXCUTS cuts */
      size = glp_ios_pool_size(tree);
      for (k = 1; k <= nv; k++)
      {  if (glp_ios_pool_size(tree) - size >= MAXCUTS) break;
         gen_cut(tree, worka, var[k].j);
      }
      /* free working arrays */
      xfree(var);
      xfree(worka->ind);
      xfree(worka->val);
      xfree(worka->phi);
      return;
}

/* eof */
