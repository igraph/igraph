/* glpios09.c (branching heuristics) */

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

/***********************************************************************
*  NAME
*
*  ios_choose_var - select variable to branch on
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  int ios_choose_var(glp_tree *T, int *next);
*
*  The routine ios_choose_var chooses a variable from the candidate
*  list to branch on. Additionally the routine provides a flag stored
*  in the location next to suggests which of the child subproblems
*  should be solved next.
*
*  RETURNS
*
*  The routine ios_choose_var returns the ordinal number of the column
*  choosen. */

static int branch_first(glp_tree *T, int *next);
static int branch_last(glp_tree *T, int *next);
static int branch_mostf(glp_tree *T, int *next);
static int branch_drtom(glp_tree *T, int *next);

int ios_choose_var(glp_tree *T, int *next)
{     int j;
      if (T->parm->br_tech == GLP_BR_FFV)
      {  /* branch on first fractional variable */
         j = branch_first(T, next);
      }
      else if (T->parm->br_tech == GLP_BR_LFV)
      {  /* branch on last fractional variable */
         j = branch_last(T, next);
      }
      else if (T->parm->br_tech == GLP_BR_MFV)
      {  /* branch on most fractional variable */
         j = branch_mostf(T, next);
      }
      else if (T->parm->br_tech == GLP_BR_DTH)
      {  /* branch using the heuristic by Dreebeck and Tomlin */
         j = branch_drtom(T, next);
      }
      else if (T->parm->br_tech == GLP_BR_PCH)
      {  /* hybrid pseudocost heuristic */
         j = ios_pcost_branch(T, next);
      }
      else
         xassert(T != T);
      return j;
}

/***********************************************************************
*  branch_first - choose first branching variable
*
*  This routine looks up the list of structural variables and chooses
*  the first one, which is of integer kind and has fractional value in
*  optimal solution to the current LP relaxation.
*
*  This routine also selects the branch to be solved next where integer
*  infeasibility of the chosen variable is less than in other one. */

static int branch_first(glp_tree *T, int *_next)
{     int j, next;
      double beta;
      /* choose the column to branch on */
      for (j = 1; j <= T->n; j++)
         if (T->non_int[j]) break;
      xassert(1 <= j && j <= T->n);
      /* select the branch to be solved next */
      beta = glp_get_col_prim(T->mip, j);
      if (beta - floor(beta) < ceil(beta) - beta)
         next = GLP_DN_BRNCH;
      else
         next = GLP_UP_BRNCH;
      *_next = next;
      return j;
}

/***********************************************************************
*  branch_last - choose last branching variable
*
*  This routine looks up the list of structural variables and chooses
*  the last one, which is of integer kind and has fractional value in
*  optimal solution to the current LP relaxation.
*
*  This routine also selects the branch to be solved next where integer
*  infeasibility of the chosen variable is less than in other one. */

static int branch_last(glp_tree *T, int *_next)
{     int j, next;
      double beta;
      /* choose the column to branch on */
      for (j = T->n; j >= 1; j--)
         if (T->non_int[j]) break;
      xassert(1 <= j && j <= T->n);
      /* select the branch to be solved next */
      beta = glp_get_col_prim(T->mip, j);
      if (beta - floor(beta) < ceil(beta) - beta)
         next = GLP_DN_BRNCH;
      else
         next = GLP_UP_BRNCH;
      *_next = next;
      return j;
}

/***********************************************************************
*  branch_mostf - choose most fractional branching variable
*
*  This routine looks up the list of structural variables and chooses
*  that one, which is of integer kind and has most fractional value in
*  optimal solution to the current LP relaxation.
*
*  This routine also selects the branch to be solved next where integer
*  infeasibility of the chosen variable is less than in other one.
*
*  (Alexander Martin notices that "...most infeasible is as good as
*  random...".) */

static int branch_mostf(glp_tree *T, int *_next)
{     int j, jj, next;
      double beta, most, temp;
      /* choose the column to branch on */
      jj = 0, most = DBL_MAX;
      for (j = 1; j <= T->n; j++)
      {  if (T->non_int[j])
         {  beta = glp_get_col_prim(T->mip, j);
            temp = floor(beta) + 0.5;
            if (most > fabs(beta - temp))
            {  jj = j, most = fabs(beta - temp);
               if (beta < temp)
                  next = GLP_DN_BRNCH;
               else
                  next = GLP_UP_BRNCH;
            }
         }
      }
      *_next = next;
      return jj;
}

/***********************************************************************
*  branch_drtom - choose branching var using Driebeck-Tomlin heuristic
*
*  This routine chooses a structural variable, which is required to be
*  integral and has fractional value in optimal solution of the current
*  LP relaxation, using a heuristic proposed by Driebeck and Tomlin.
*
*  The routine also selects the branch to be solved next, again due to
*  Driebeck and Tomlin.
*
*  This routine is based on the heuristic proposed in:
*
*  Driebeck N.J. An algorithm for the solution of mixed-integer
*  programming problems, Management Science, 12: 576-87 (1966);
*
*  and improved in:
*
*  Tomlin J.A. Branch and bound methods for integer and non-convex
*  programming, in J.Abadie (ed.), Integer and Nonlinear Programming,
*  North-Holland, Amsterdam, pp. 437-50 (1970).
*
*  Must note that this heuristic is time-expensive, because computing
*  one-step degradation (see the routine below) requires one BTRAN for
*  each fractional-valued structural variable. */

static int branch_drtom(glp_tree *T, int *_next)
{     glp_prob *mip = T->mip;
      int m = mip->m;
      int n = mip->n;
      unsigned char *non_int = T->non_int;
      int j, jj, k, t, next, kase, len, stat, *ind;
      double x, dk, alfa, delta_j, delta_k, delta_z, dz_dn, dz_up,
         dd_dn, dd_up, degrad, *val;
      /* basic solution of LP relaxation must be optimal */
      xassert(glp_get_status(mip) == GLP_OPT);
      /* allocate working arrays */
      ind = xcalloc(1+n, sizeof(int));
      val = xcalloc(1+n, sizeof(double));
      /* nothing has been chosen so far */
      jj = 0, degrad = -1.0;
      /* walk through the list of columns (structural variables) */
      for (j = 1; j <= n; j++)
      {  /* if j-th column is not marked as fractional, skip it */
         if (!non_int[j]) continue;
         /* obtain (fractional) value of j-th column in basic solution
            of LP relaxation */
         x = glp_get_col_prim(mip, j);
         /* since the value of j-th column is fractional, the column is
            basic; compute corresponding row of the simplex table */
         len = glp_eval_tab_row(mip, m+j, ind, val);
         /* the following fragment computes a change in the objective
            function: delta Z = new Z - old Z, where old Z is the
            objective value in the current optimal basis, and new Z is
            the objective value in the adjacent basis, for two cases:
            1) if new upper bound ub' = floor(x[j]) is introduced for
               j-th column (down branch);
            2) if new lower bound lb' = ceil(x[j]) is introduced for
               j-th column (up branch);
            since in both cases the solution remaining dual feasible
            becomes primal infeasible, one implicit simplex iteration
            is performed to determine the change delta Z;
            it is obvious that new Z, which is never better than old Z,
            is a lower (minimization) or upper (maximization) bound of
            the objective function for down- and up-branches. */
         for (kase = -1; kase <= +1; kase += 2)
         {  /* if kase < 0, the new upper bound of x[j] is introduced;
               in this case x[j] should decrease in order to leave the
               basis and go to its new upper bound */
            /* if kase > 0, the new lower bound of x[j] is introduced;
               in this case x[j] should increase in order to leave the
               basis and go to its new lower bound */
            /* apply the dual ratio test in order to determine which
               auxiliary or structural variable should enter the basis
               to keep dual feasibility */
            k = glp_dual_rtest(mip, len, ind, val, kase, 1e-9);
            if (k != 0) k = ind[k];
            /* if no non-basic variable has been chosen, LP relaxation
               of corresponding branch being primal infeasible and dual
               unbounded has no primal feasible solution; in this case
               the change delta Z is formally set to infinity */
            if (k == 0)
            {  delta_z =
                  (T->mip->dir == GLP_MIN ? +DBL_MAX : -DBL_MAX);
               goto skip;
            }
            /* row of the simplex table that corresponds to non-basic
               variable x[k] choosen by the dual ratio test is:
                  x[j] = ... + alfa * x[k] + ...
               where alfa is the influence coefficient (an element of
               the simplex table row) */
            /* determine the coefficient alfa */
            for (t = 1; t <= len; t++) if (ind[t] == k) break;
            xassert(1 <= t && t <= len);
            alfa = val[t];
            /* since in the adjacent basis the variable x[j] becomes
               non-basic, knowing its value in the current basis we can
               determine its change delta x[j] = new x[j] - old x[j] */
            delta_j = (kase < 0 ? floor(x) : ceil(x)) - x;
            /* and knowing the coefficient alfa we can determine the
               corresponding change delta x[k] = new x[k] - old x[k],
               where old x[k] is a value of x[k] in the current basis,
               and new x[k] is a value of x[k] in the adjacent basis */
            delta_k = delta_j / alfa;
            /* Tomlin noticed that if the variable x[k] is of integer
               kind, its change cannot be less (eventually) than one in
               the magnitude */
            if (k > m && glp_get_col_kind(mip, k-m) != GLP_CV)
            {  /* x[k] is structural integer variable */
               if (fabs(delta_k - floor(delta_k + 0.5)) > 1e-3)
               {  if (delta_k > 0.0)
                     delta_k = ceil(delta_k);  /* +3.14 -> +4 */
                  else
                     delta_k = floor(delta_k); /* -3.14 -> -4 */
               }
            }
            /* now determine the status and reduced cost of x[k] in the
               current basis */
            if (k <= m)
            {  stat = glp_get_row_stat(mip, k);
               dk = glp_get_row_dual(mip, k);
            }
            else
            {  stat = glp_get_col_stat(mip, k-m);
               dk = glp_get_col_dual(mip, k-m);
            }
            /* if the current basis is dual degenerate, some reduced
               costs which are close to zero may have wrong sign due to
               round-off errors, so correct the sign of d[k] */
            switch (T->mip->dir)
            {  case GLP_MIN:
                  if (stat == GLP_NL && dk < 0.0 ||
                      stat == GLP_NU && dk > 0.0 ||
                      stat == GLP_NF) dk = 0.0;
                  break;
               case GLP_MAX:
                  if (stat == GLP_NL && dk > 0.0 ||
                      stat == GLP_NU && dk < 0.0 ||
                      stat == GLP_NF) dk = 0.0;
                  break;
               default:
                  xassert(T != T);
            }
            /* now knowing the change of x[k] and its reduced cost d[k]
               we can compute the corresponding change in the objective
               function delta Z = new Z - old Z = d[k] * delta x[k];
               note that due to Tomlin's modification new Z can be even
               worse than in the adjacent basis */
            delta_z = dk * delta_k;
skip:       /* new Z is never better than old Z, therefore the change
               delta Z is always non-negative (in case of minimization)
               or non-positive (in case of maximization) */
            switch (T->mip->dir)
            {  case GLP_MIN: xassert(delta_z >= 0.0); break;
               case GLP_MAX: xassert(delta_z <= 0.0); break;
               default: xassert(T != T);
            }
            /* save the change in the objective fnction for down- and
               up-branches, respectively */
            if (kase < 0) dz_dn = delta_z; else dz_up = delta_z;
         }
         /* thus, in down-branch no integer feasible solution can be
            better than Z + dz_dn, and in up-branch no integer feasible
            solution can be better than Z + dz_up, where Z is value of
            the objective function in the current basis */
         /* following the heuristic by Driebeck and Tomlin we choose a
            column (i.e. structural variable) which provides largest
            degradation of the objective function in some of branches;
            besides, we select the branch with smaller degradation to
            be solved next and keep other branch with larger degradation
            in the active list hoping to minimize the number of further
            backtrackings */
         if (degrad < fabs(dz_dn) || degrad < fabs(dz_up))
         {  jj = j;
            if (fabs(dz_dn) < fabs(dz_up))
            {  /* select down branch to be solved next */
               next = GLP_DN_BRNCH;
               degrad = fabs(dz_up);
            }
            else
            {  /* select up branch to be solved next */
               next = GLP_UP_BRNCH;
               degrad = fabs(dz_dn);
            }
            /* save the objective changes for printing */
            dd_dn = dz_dn, dd_up = dz_up;
            /* if down- or up-branch has no feasible solution, we does
               not need to consider other candidates (in principle, the
               corresponding branch could be pruned right now) */
            if (degrad == DBL_MAX) break;
         }
      }
      /* free working arrays */
      xfree(ind);
      xfree(val);
      /* something must be chosen */
      xassert(1 <= jj && jj <= n);
#if 1 /* 02/XI-2009 */
      if (degrad < 1e-6 * (1.0 + 0.001 * fabs(mip->obj_val)))
      {  jj = branch_mostf(T, &next);
         goto done;
      }
#endif
      if (T->parm->msg_lev >= GLP_MSG_DBG)
      {  xprintf("branch_drtom: column %d chosen to branch on\n", jj);
         if (fabs(dd_dn) == DBL_MAX)
            xprintf("branch_drtom: down-branch is infeasible\n");
         else
            xprintf("branch_drtom: down-branch bound is %.9e\n",
               glp_get_obj_val(mip) + dd_dn);
         if (fabs(dd_up) == DBL_MAX)
            xprintf("branch_drtom: up-branch   is infeasible\n");
         else
            xprintf("branch_drtom: up-branch   bound is %.9e\n",
               glp_get_obj_val(mip) + dd_up);
      }
done: *_next = next;
      return jj;
}

/**********************************************************************/

struct csa
{     /* common storage area */
      int *dn_cnt; /* int dn_cnt[1+n]; */
      /* dn_cnt[j] is the number of subproblems, whose LP relaxations
         have been solved and which are down-branches for variable x[j];
         dn_cnt[j] = 0 means the down pseudocost is uninitialized */
      double *dn_sum; /* double dn_sum[1+n]; */
      /* dn_sum[j] is the sum of per unit degradations of the objective
         over all dn_cnt[j] subproblems */
      int *up_cnt; /* int up_cnt[1+n]; */
      /* up_cnt[j] is the number of subproblems, whose LP relaxations
         have been solved and which are up-branches for variable x[j];
         up_cnt[j] = 0 means the up pseudocost is uninitialized */
      double *up_sum; /* double up_sum[1+n]; */
      /* up_sum[j] is the sum of per unit degradations of the objective
         over all up_cnt[j] subproblems */
};

void *ios_pcost_init(glp_tree *tree)
{     /* initialize working data used on pseudocost branching */
      struct csa *csa;
      int n = tree->n, j;
      csa = xmalloc(sizeof(struct csa));
      csa->dn_cnt = xcalloc(1+n, sizeof(int));
      csa->dn_sum = xcalloc(1+n, sizeof(double));
      csa->up_cnt = xcalloc(1+n, sizeof(int));
      csa->up_sum = xcalloc(1+n, sizeof(double));
      for (j = 1; j <= n; j++)
      {  csa->dn_cnt[j] = csa->up_cnt[j] = 0;
         csa->dn_sum[j] = csa->up_sum[j] = 0.0;
      }
      return csa;
}

static double eval_degrad(glp_prob *P, int j, double bnd)
{     /* compute degradation of the objective on fixing x[j] at given
         value with a limited number of dual simplex iterations */
      /* this routine fixes column x[j] at specified value bnd,
         solves resulting LP, and returns a lower bound to degradation
         of the objective, degrad >= 0 */
      glp_prob *lp;
      glp_smcp parm;
      int ret;
      double degrad;
      /* the current basis must be optimal */
      xassert(glp_get_status(P) == GLP_OPT);
      /* create a copy of P */
      lp = glp_create_prob();
      glp_copy_prob(lp, P, 0);
      /* fix column x[j] at specified value */
      glp_set_col_bnds(lp, j, GLP_FX, bnd, bnd);
      /* try to solve resulting LP */
      glp_init_smcp(&parm);
      parm.msg_lev = GLP_MSG_OFF;
      parm.meth = GLP_DUAL;
      parm.it_lim = 30;
      parm.out_dly = 1000;
      parm.meth = GLP_DUAL;
      ret = glp_simplex(lp, &parm);
      if (ret == 0 || ret == GLP_EITLIM)
      {  if (glp_get_prim_stat(lp) == GLP_NOFEAS)
         {  /* resulting LP has no primal feasible solution */
            degrad = DBL_MAX;
         }
         else if (glp_get_dual_stat(lp) == GLP_FEAS)
         {  /* resulting basis is optimal or at least dual feasible,
               so we have the correct lower bound to degradation */
            if (P->dir == GLP_MIN)
               degrad = lp->obj_val - P->obj_val;
            else if (P->dir == GLP_MAX)
               degrad = P->obj_val - lp->obj_val;
            else
               xassert(P != P);
            /* degradation cannot be negative by definition */
            /* note that the lower bound to degradation may be close
               to zero even if its exact value is zero due to round-off
               errors on computing the objective value */
            if (degrad < 1e-6 * (1.0 + 0.001 * fabs(P->obj_val)))
               degrad = 0.0;
         }
         else
         {  /* the final basis reported by the simplex solver is dual
               infeasible, so we cannot determine a non-trivial lower
               bound to degradation */
            degrad = 0.0;
         }
      }
      else
      {  /* the simplex solver failed */
         degrad = 0.0;
      }
      /* delete the copy of P */
      glp_delete_prob(lp);
      return degrad;
}

void ios_pcost_update(glp_tree *tree)
{     /* update history information for pseudocost branching */
      /* this routine is called every time when LP relaxation of the
         current subproblem has been solved to optimality with all lazy
         and cutting plane constraints included */
      int j;
      double dx, dz, psi;
      struct csa *csa = tree->pcost;
      xassert(csa != NULL);
      xassert(tree->curr != NULL);
      /* if the current subproblem is the root, skip updating */
      if (tree->curr->up == NULL) goto skip;
      /* determine branching variable x[j], which was used in the
         parent subproblem to create the current subproblem */
      j = tree->curr->up->br_var;
      xassert(1 <= j && j <= tree->n);
      /* determine the change dx[j] = new x[j] - old x[j],
         where new x[j] is a value of x[j] in optimal solution to LP
         relaxation of the current subproblem, old x[j] is a value of
         x[j] in optimal solution to LP relaxation of the parent
         subproblem */
      dx = tree->mip->col[j]->prim - tree->curr->up->br_val;
      xassert(dx != 0.0);
      /* determine corresponding change dz = new dz - old dz in the
         objective function value */
      dz = tree->mip->obj_val - tree->curr->up->lp_obj;
      /* determine per unit degradation of the objective function */
      psi = fabs(dz / dx);
      /* update history information */
      if (dx < 0.0)
      {  /* the current subproblem is down-branch */
         csa->dn_cnt[j]++;
         csa->dn_sum[j] += psi;
      }
      else /* dx > 0.0 */
      {  /* the current subproblem is up-branch */
         csa->up_cnt[j]++;
         csa->up_sum[j] += psi;
      }
skip: return;
}

void ios_pcost_free(glp_tree *tree)
{     /* free working area used on pseudocost branching */
      struct csa *csa = tree->pcost;
      xassert(csa != NULL);
      xfree(csa->dn_cnt);
      xfree(csa->dn_sum);
      xfree(csa->up_cnt);
      xfree(csa->up_sum);
      xfree(csa);
      tree->pcost = NULL;
      return;
}

static double eval_psi(glp_tree *T, int j, int brnch)
{     /* compute estimation of pseudocost of variable x[j] for down-
         or up-branch */
      struct csa *csa = T->pcost;
      double beta, degrad, psi;
      xassert(csa != NULL);
      xassert(1 <= j && j <= T->n);
      if (brnch == GLP_DN_BRNCH)
      {  /* down-branch */
         if (csa->dn_cnt[j] == 0)
         {  /* initialize down pseudocost */
            beta = T->mip->col[j]->prim;
            degrad = eval_degrad(T->mip, j, floor(beta));
            if (degrad == DBL_MAX)
            {  psi = DBL_MAX;
               goto done;
            }
            csa->dn_cnt[j] = 1;
            csa->dn_sum[j] = degrad / (beta - floor(beta));
         }
         psi = csa->dn_sum[j] / (double)csa->dn_cnt[j];
      }
      else if (brnch == GLP_UP_BRNCH)
      {  /* up-branch */
         if (csa->up_cnt[j] == 0)
         {  /* initialize up pseudocost */
            beta = T->mip->col[j]->prim;
            degrad = eval_degrad(T->mip, j, ceil(beta));
            if (degrad == DBL_MAX)
            {  psi = DBL_MAX;
               goto done;
            }
            csa->up_cnt[j] = 1;
            csa->up_sum[j] = degrad / (ceil(beta) - beta);
         }
         psi = csa->up_sum[j] / (double)csa->up_cnt[j];
      }
      else
         xassert(brnch != brnch);
done: return psi;
}

static void progress(glp_tree *T)
{     /* display progress of pseudocost initialization */
      struct csa *csa = T->pcost;
      int j, nv = 0, ni = 0;
      for (j = 1; j <= T->n; j++)
      {  if (glp_ios_can_branch(T, j))
         {  nv++;
            if (csa->dn_cnt[j] > 0 && csa->up_cnt[j] > 0) ni++;
         }
      }
      xprintf("Pseudocosts initialized for %d of %d variables\n",
         ni, nv);
      return;
}

int ios_pcost_branch(glp_tree *T, int *_next)
{     /* choose branching variable with pseudocost branching */
#if 0 /* 10/VI-2013 */
      glp_long t = xtime();
#else
      double t = xtime();
#endif
      int j, jjj, sel;
      double beta, psi, d1, d2, d, dmax;
      /* initialize the working arrays */
      if (T->pcost == NULL)
         T->pcost = ios_pcost_init(T);
      /* nothing has been chosen so far */
      jjj = 0, dmax = -1.0;
      /* go through the list of branching candidates */
      for (j = 1; j <= T->n; j++)
      {  if (!glp_ios_can_branch(T, j)) continue;
         /* determine primal value of x[j] in optimal solution to LP
            relaxation of the current subproblem */
         beta = T->mip->col[j]->prim;
         /* estimate pseudocost of x[j] for down-branch */
         psi = eval_psi(T, j, GLP_DN_BRNCH);
         if (psi == DBL_MAX)
         {  /* down-branch has no primal feasible solution */
            jjj = j, sel = GLP_DN_BRNCH;
            goto done;
         }
         /* estimate degradation of the objective for down-branch */
         d1 = psi * (beta - floor(beta));
         /* estimate pseudocost of x[j] for up-branch */
         psi = eval_psi(T, j, GLP_UP_BRNCH);
         if (psi == DBL_MAX)
         {  /* up-branch has no primal feasible solution */
            jjj = j, sel = GLP_UP_BRNCH;
            goto done;
         }
         /* estimate degradation of the objective for up-branch */
         d2 = psi * (ceil(beta) - beta);
         /* determine d = max(d1, d2) */
         d = (d1 > d2 ? d1 : d2);
         /* choose x[j] which provides maximal estimated degradation of
            the objective either in down- or up-branch */
         if (dmax < d)
         {  dmax = d;
            jjj = j;
            /* continue the search from a subproblem, where degradation
               is less than in other one */
            sel = (d1 <= d2 ? GLP_DN_BRNCH : GLP_UP_BRNCH);
         }
         /* display progress of pseudocost initialization */
         if (T->parm->msg_lev >= GLP_ON)
         {  if (xdifftime(xtime(), t) >= 10.0)
            {  progress(T);
               t = xtime();
            }
         }
      }
      if (dmax == 0.0)
      {  /* no degradation is indicated; choose a variable having most
            fractional value */
         jjj = branch_mostf(T, &sel);
      }
done: *_next = sel;
      return jjj;
}

/* eof */
