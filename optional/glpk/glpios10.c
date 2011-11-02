/* glpios10.c (feasibility pump heuristic) */

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
#include "glprng.h"

/***********************************************************************
*  NAME
*
*  ios_feas_pump - feasibility pump heuristic
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_feas_pump(glp_tree *T);
*
*  DESCRIPTION
*
*  The routine ios_feas_pump is a simple implementation of the Feasi-
*  bility Pump heuristic.
*
*  REFERENCES
*
*  M.Fischetti, F.Glover, and A.Lodi. "The feasibility pump." Math.
*  Program., Ser. A 104, pp. 91-104 (2005). */

struct VAR
{     /* binary variable */
      int j;
      /* ordinal number */
      int x;
      /* value in the rounded solution (0 or 1) */
      double d;
      /* sorting key */
};

static int fcmp(const void *x, const void *y)
{     /* comparison routine */
      const struct VAR *vx = x, *vy = y;
      if (vx->d > vy->d)
         return -1;
      else if (vx->d < vy->d)
         return +1;
      else
         return 0;
}

void ios_feas_pump(glp_tree *T)
{     glp_prob *P = T->mip;
      int n = P->n;
      glp_prob *lp = NULL;
      struct VAR *var = NULL;
      RNG *rand = NULL;
      GLPCOL *col;
      glp_smcp parm;
      int j, k, new_x, nfail, npass, nv, ret, stalling;
      double dist, tol;
      xassert(glp_get_status(P) == GLP_OPT);
      /* this heuristic is applied only once on the root level */
      if (!(T->curr->level == 0 && T->curr->solved == 1)) goto done;
      /* determine number of binary variables */
      nv = 0;
      for (j = 1; j <= n; j++)
      {  col = P->col[j];
         /* if x[j] is continuous, skip it */
         if (col->kind == GLP_CV) continue;
         /* if x[j] is fixed, skip it */
         if (col->type == GLP_FX) continue;
         /* x[j] is non-fixed integer */
         xassert(col->kind == GLP_IV);
         if (col->type == GLP_DB && col->lb == 0.0 && col->ub == 1.0)
         {  /* x[j] is binary */
            nv++;
         }
         else
         {  /* x[j] is general integer */
            if (T->parm->msg_lev >= GLP_MSG_ALL)
               xprintf("FPUMP heuristic cannot be applied due to genera"
                  "l integer variables\n");
            goto done;
         }
      }
      /* there must be at least one binary variable */
      if (nv == 0) goto done;
      if (T->parm->msg_lev >= GLP_MSG_ALL)
         xprintf("Applying FPUMP heuristic...\n");
      /* build the list of binary variables */
      var = xcalloc(1+nv, sizeof(struct VAR));
      k = 0;
      for (j = 1; j <= n; j++)
      {  col = P->col[j];
         if (col->kind == GLP_IV && col->type == GLP_DB)
            var[++k].j = j;
      }
      xassert(k == nv);
      /* create working problem object */
      lp = glp_create_prob();
more: /* copy the original problem object to keep it intact */
      glp_copy_prob(lp, P, GLP_OFF);
      /* we are interested to find an integer feasible solution, which
         is better than the best known one */
      if (P->mip_stat == GLP_FEAS)
      {  int *ind;
         double *val, bnd;
         /* add a row and make it identical to the objective row */
         glp_add_rows(lp, 1);
         ind = xcalloc(1+n, sizeof(int));
         val = xcalloc(1+n, sizeof(double));
         for (j = 1; j <= n; j++)
         {  ind[j] = j;
            val[j] = P->col[j]->coef;
         }
         glp_set_mat_row(lp, lp->m, n, ind, val);
         xfree(ind);
         xfree(val);
         /* introduce upper (minimization) or lower (maximization)
            bound to the original objective function; note that this
            additional constraint is not violated at the optimal point
            to LP relaxation */
#if 0 /* modified by xypron <xypron.glpk@gmx.de> */
         if (P->dir == GLP_MIN)
         {  bnd = P->mip_obj - 0.10 * (1.0 + fabs(P->mip_obj));
            if (bnd < P->obj_val) bnd = P->obj_val;
            glp_set_row_bnds(lp, lp->m, GLP_UP, 0.0, bnd - P->c0);
         }
         else if (P->dir == GLP_MAX)
         {  bnd = P->mip_obj + 0.10 * (1.0 + fabs(P->mip_obj));
            if (bnd > P->obj_val) bnd = P->obj_val;
            glp_set_row_bnds(lp, lp->m, GLP_LO, bnd - P->c0, 0.0);
         }
         else
            xassert(P != P);
#else
         bnd = 0.1 * P->obj_val + 0.9 * P->mip_obj;
         /* xprintf("bnd = %f\n", bnd); */
         if (P->dir == GLP_MIN)
            glp_set_row_bnds(lp, lp->m, GLP_UP, 0.0, bnd - P->c0);
         else if (P->dir == GLP_MAX)
            glp_set_row_bnds(lp, lp->m, GLP_LO, bnd - P->c0, 0.0);
         else
            xassert(P != P);
#endif
      }
      /* reset pass count */
      npass = 0;
      /* invalidate the rounded point */
      for (k = 1; k <= nv; k++)
         var[k].x = -1;
pass: /* next pass starts here */
      npass++;
      if (T->parm->msg_lev >= GLP_MSG_ALL)
         xprintf("Pass %d\n", npass);
      /* initialize minimal distance between the basic point and the
         rounded one obtained during this pass */
      dist = DBL_MAX;
      /* reset failure count (the number of succeeded iterations failed
         to improve the distance) */
      nfail = 0;
      /* if it is not the first pass, perturb the last rounded point
         rather than construct it from the basic solution */
      if (npass > 1)
      {  double rho, temp;
         if (rand == NULL)
            rand = rng_create_rand();
         for (k = 1; k <= nv; k++)
         {  j = var[k].j;
            col = lp->col[j];
            rho = rng_uniform(rand, -0.3, 0.7);
            if (rho < 0.0) rho = 0.0;
            temp = fabs((double)var[k].x - col->prim);
            if (temp + rho > 0.5) var[k].x = 1 - var[k].x;
         }
         goto skip;
      }
loop: /* innermost loop begins here */
      /* round basic solution (which is assumed primal feasible) */
      stalling = 1;
      for (k = 1; k <= nv; k++)
      {  col = lp->col[var[k].j];
         if (col->prim < 0.5)
         {  /* rounded value is 0 */
            new_x = 0;
         }
         else
         {  /* rounded value is 1 */
            new_x = 1;
         }
         if (var[k].x != new_x)
         {  stalling = 0;
            var[k].x = new_x;
         }
      }
      /* if the rounded point has not changed (stalling), choose and
         flip some its entries heuristically */
      if (stalling)
      {  /* compute d[j] = |x[j] - round(x[j])| */
         for (k = 1; k <= nv; k++)
         {  col = lp->col[var[k].j];
            var[k].d = fabs(col->prim - (double)var[k].x);
         }
         /* sort the list of binary variables by descending d[j] */
         qsort(&var[1], nv, sizeof(struct VAR), fcmp);
         /* choose and flip some rounded components */
         for (k = 1; k <= nv; k++)
         {  if (k >= 5 && var[k].d < 0.35 || k >= 10) break;
            var[k].x = 1 - var[k].x;
         }
      }
skip: /* check if the time limit has been exhausted */
      if (T->parm->tm_lim < INT_MAX &&
         (double)(T->parm->tm_lim - 1) <=
         1000.0 * xdifftime(xtime(), T->tm_beg)) goto done;
      /* build the objective, which is the distance between the current
         (basic) point and the rounded one */
      lp->dir = GLP_MIN;
      lp->c0 = 0.0;
      for (j = 1; j <= n; j++)
         lp->col[j]->coef = 0.0;
      for (k = 1; k <= nv; k++)
      {  j = var[k].j;
         if (var[k].x == 0)
            lp->col[j]->coef = +1.0;
         else
         {  lp->col[j]->coef = -1.0;
            lp->c0 += 1.0;
         }
      }
      /* minimize the distance with the simplex method */
      glp_init_smcp(&parm);
      if (T->parm->msg_lev <= GLP_MSG_ERR)
         parm.msg_lev = T->parm->msg_lev;
      else if (T->parm->msg_lev <= GLP_MSG_ALL)
      {  parm.msg_lev = GLP_MSG_ON;
         parm.out_dly = 10000;
      }
      ret = glp_simplex(lp, &parm);
      if (ret != 0)
      {  if (T->parm->msg_lev >= GLP_MSG_ERR)
            xprintf("Warning: glp_simplex returned %d\n", ret);
         goto done;
      }
      ret = glp_get_status(lp);
      if (ret != GLP_OPT)
      {  if (T->parm->msg_lev >= GLP_MSG_ERR)
            xprintf("Warning: glp_get_status returned %d\n", ret);
         goto done;
      }
      if (T->parm->msg_lev >= GLP_MSG_DBG)
         xprintf("delta = %g\n", lp->obj_val);
      /* check if the basic solution is integer feasible; note that it
         may be so even if the minimial distance is positive */
      tol = 0.3 * T->parm->tol_int;
      for (k = 1; k <= nv; k++)
      {  col = lp->col[var[k].j];
         if (tol < col->prim && col->prim < 1.0 - tol) break;
      }
      if (k > nv)
      {  /* okay; the basic solution seems to be integer feasible */
         double *x = xcalloc(1+n, sizeof(double));
         for (j = 1; j <= n; j++)
         {  x[j] = lp->col[j]->prim;
            if (P->col[j]->kind == GLP_IV) x[j] = floor(x[j] + 0.5);
         }
#if 1 /* modified by xypron <xypron.glpk@gmx.de> */
         /* reset direction and right-hand side of objective */
         lp->c0  = P->c0;
         lp->dir = P->dir;
         /* fix integer variables */
         for (k = 1; k <= nv; k++)
         {  lp->col[var[k].j]->lb   = x[var[k].j];
            lp->col[var[k].j]->ub   = x[var[k].j];
            lp->col[var[k].j]->type = GLP_FX;
         }
         /* copy original objective function */
         for (j = 1; j <= n; j++)
            lp->col[j]->coef = P->col[j]->coef;
         /* solve original LP and copy result */
         ret = glp_simplex(lp, &parm);
         if (ret != 0)
         {  if (T->parm->msg_lev >= GLP_MSG_ERR)
               xprintf("Warning: glp_simplex returned %d\n", ret);
            goto done;
         }
         ret = glp_get_status(lp);
         if (ret != GLP_OPT)
         {  if (T->parm->msg_lev >= GLP_MSG_ERR)
               xprintf("Warning: glp_get_status returned %d\n", ret);
            goto done;
         }
         for (j = 1; j <= n; j++)
            if (P->col[j]->kind != GLP_IV) x[j] = lp->col[j]->prim;
#endif
         ret = glp_ios_heur_sol(T, x);
         xfree(x);
         if (ret == 0)
         {  /* the integer solution is accepted */
            if (ios_is_hopeful(T, T->curr->bound))
            {  /* it is reasonable to apply the heuristic once again */
               goto more;
            }
            else
            {  /* the best known integer feasible solution just found
                  is close to optimal solution to LP relaxation */
               goto done;
            }
         }
      }
      /* the basic solution is fractional */
      if (dist == DBL_MAX ||
          lp->obj_val <= dist - 1e-6 * (1.0 + dist))
      {  /* the distance is reducing */
         nfail = 0, dist = lp->obj_val;
      }
      else
      {  /* improving the distance failed */
         nfail++;
      }
      if (nfail < 3) goto loop;
      if (npass < 5) goto pass;
done: /* delete working objects */
      if (lp != NULL) glp_delete_prob(lp);
      if (var != NULL) xfree(var);
      if (rand != NULL) rng_delete_rand(rand);
      return;
}

/* eof */
