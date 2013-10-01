/* glpios03.c (branch-and-cut driver) */

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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wsometimes-uninitialized"
#pragma clang diagnostic ignored "-Wlogical-op-parentheses"
#endif

#include "glpios.h"

/***********************************************************************
*  show_progress - display current progress of the search
*
*  This routine displays some information about current progress of the
*  search.
*
*  The information includes:
*
*  the current number of iterations performed by the simplex solver;
*
*  the objective value for the best known integer feasible solution,
*  which is upper (minimization) or lower (maximization) global bound
*  for optimal solution of the original mip problem;
*
*  the best local bound for active nodes, which is lower (minimization)
*  or upper (maximization) global bound for optimal solution of the
*  original mip problem;
*
*  the relative mip gap, in percents;
*
*  the number of open (active) subproblems;
*
*  the number of completely explored subproblems, i.e. whose nodes have
*  been removed from the tree. */

static void show_progress(glp_tree *T, int bingo)
{     int p;
      double temp;
      char best_mip[50], best_bound[50], *rho, rel_gap[50];
      /* format the best known integer feasible solution */
      if (T->mip->mip_stat == GLP_FEAS)
         sprintf(best_mip, "%17.9e", T->mip->mip_obj);
      else
         sprintf(best_mip, "%17s", "not found yet");
      /* determine reference number of an active subproblem whose local
         bound is best */
      p = ios_best_node(T);
      /* format the best bound */
      if (p == 0)
         sprintf(best_bound, "%17s", "tree is empty");
      else
      {  temp = T->slot[p].node->bound;
         if (temp == -DBL_MAX)
            sprintf(best_bound, "%17s", "-inf");
         else if (temp == +DBL_MAX)
            sprintf(best_bound, "%17s", "+inf");
         else
            sprintf(best_bound, "%17.9e", temp);
      }
      /* choose the relation sign between global bounds */
      if (T->mip->dir == GLP_MIN)
         rho = ">=";
      else if (T->mip->dir == GLP_MAX)
         rho = "<=";
      else
         xassert(T != T);
      /* format the relative mip gap */
      temp = ios_relative_gap(T);
      if (temp == 0.0)
         sprintf(rel_gap, "  0.0%%");
      else if (temp < 0.001)
         sprintf(rel_gap, "< 0.1%%");
      else if (temp <= 9.999)
         sprintf(rel_gap, "%5.1f%%", 100.0 * temp);
      else
         sprintf(rel_gap, "%6s", "");
      /* display progress of the search */
      xprintf("+%6d: %s %s %s %s %s (%d; %d)\n",
         T->mip->it_cnt, bingo ? ">>>>>" : "mip =", best_mip, rho,
         best_bound, rel_gap, T->a_cnt, T->t_cnt - T->n_cnt);
      T->tm_lag = xtime();
      return;
}

/***********************************************************************
*  is_branch_hopeful - check if specified branch is hopeful
*
*  This routine checks if the specified subproblem can have an integer
*  optimal solution which is better than the best known one.
*
*  The check is based on comparison of the local objective bound stored
*  in the subproblem descriptor and the incumbent objective value which
*  is the global objective bound.
*
*  If there is a chance that the specified subproblem can have a better
*  integer optimal solution, the routine returns non-zero. Otherwise, if
*  the corresponding branch can pruned, zero is returned. */

static int is_branch_hopeful(glp_tree *T, int p)
{     xassert(1 <= p && p <= T->nslots);
      xassert(T->slot[p].node != NULL);
      return ios_is_hopeful(T, T->slot[p].node->bound);
}

/***********************************************************************
*  check_integrality - check integrality of basic solution
*
*  This routine checks if the basic solution of LP relaxation of the
*  current subproblem satisfies to integrality conditions, i.e. that all
*  variables of integer kind have integral primal values. (The solution
*  is assumed to be optimal.)
*
*  For each variable of integer kind the routine computes the following
*  quantity:
*
*     ii(x[j]) = min(x[j] - floor(x[j]), ceil(x[j]) - x[j]),         (1)
*
*  which is a measure of the integer infeasibility (non-integrality) of
*  x[j] (for example, ii(2.1) = 0.1, ii(3.7) = 0.3, ii(5.0) = 0). It is
*  understood that 0 <= ii(x[j]) <= 0.5, and variable x[j] is integer
*  feasible if ii(x[j]) = 0. However, due to floating-point arithmetic
*  the routine checks less restrictive condition:
*
*     ii(x[j]) <= tol_int,                                           (2)
*
*  where tol_int is a given tolerance (small positive number) and marks
*  each variable which does not satisfy to (2) as integer infeasible by
*  setting its fractionality flag.
*
*  In order to characterize integer infeasibility of the basic solution
*  in the whole the routine computes two parameters: ii_cnt, which is
*  the number of variables with the fractionality flag set, and ii_sum,
*  which is the sum of integer infeasibilities (1). */

static void check_integrality(glp_tree *T)
{     glp_prob *mip = T->mip;
      int j, type, ii_cnt = 0;
      double lb, ub, x, temp1, temp2, ii_sum = 0.0;
      /* walk through the set of columns (structural variables) */
      for (j = 1; j <= mip->n; j++)
      {  GLPCOL *col = mip->col[j];
         T->non_int[j] = 0;
         /* if the column is not integer, skip it */
         if (col->kind != GLP_IV) continue;
         /* if the column is non-basic, it is integer feasible */
         if (col->stat != GLP_BS) continue;
         /* obtain the type and bounds of the column */
         type = col->type, lb = col->lb, ub = col->ub;
         /* obtain value of the column in optimal basic solution */
         x = col->prim;
         /* if the column's primal value is close to the lower bound,
            the column is integer feasible within given tolerance */
         if (type == GLP_LO || type == GLP_DB || type == GLP_FX)
         {  temp1 = lb - T->parm->tol_int;
            temp2 = lb + T->parm->tol_int;
            if (temp1 <= x && x <= temp2) continue;
#if 0
            /* the lower bound must not be violated */
            xassert(x >= lb);
#else
            if (x < lb) continue;
#endif
         }
         /* if the column's primal value is close to the upper bound,
            the column is integer feasible within given tolerance */
         if (type == GLP_UP || type == GLP_DB || type == GLP_FX)
         {  temp1 = ub - T->parm->tol_int;
            temp2 = ub + T->parm->tol_int;
            if (temp1 <= x && x <= temp2) continue;
#if 0
            /* the upper bound must not be violated */
            xassert(x <= ub);
#else
            if (x > ub) continue;
#endif
         }
         /* if the column's primal value is close to nearest integer,
            the column is integer feasible within given tolerance */
         temp1 = floor(x + 0.5) - T->parm->tol_int;
         temp2 = floor(x + 0.5) + T->parm->tol_int;
         if (temp1 <= x && x <= temp2) continue;
         /* otherwise the column is integer infeasible */
         T->non_int[j] = 1;
         /* increase the number of fractional-valued columns */
         ii_cnt++;
         /* compute the sum of integer infeasibilities */
         temp1 = x - floor(x);
         temp2 = ceil(x) - x;
         xassert(temp1 > 0.0 && temp2 > 0.0);
         ii_sum += (temp1 <= temp2 ? temp1 : temp2);
      }
      /* store ii_cnt and ii_sum to the current problem descriptor */
      xassert(T->curr != NULL);
      T->curr->ii_cnt = ii_cnt;
      T->curr->ii_sum = ii_sum;
      /* and also display these parameters */
      if (T->parm->msg_lev >= GLP_MSG_DBG)
      {  if (ii_cnt == 0)
            xprintf("There are no fractional columns\n");
         else if (ii_cnt == 1)
            xprintf("There is one fractional column, integer infeasibil"
               "ity is %.3e\n", ii_sum);
         else
            xprintf("There are %d fractional columns, integer infeasibi"
               "lity is %.3e\n", ii_cnt, ii_sum);
      }
      return;
}

/***********************************************************************
*  record_solution - record better integer feasible solution
*
*  This routine records optimal basic solution of LP relaxation of the
*  current subproblem, which being integer feasible is better than the
*  best known integer feasible solution. */

static void record_solution(glp_tree *T)
{     glp_prob *mip = T->mip;
      int i, j;
      mip->mip_stat = GLP_FEAS;
      mip->mip_obj = mip->obj_val;
      for (i = 1; i <= mip->m; i++)
      {  GLPROW *row = mip->row[i];
         row->mipx = row->prim;
      }
      for (j = 1; j <= mip->n; j++)
      {  GLPCOL *col = mip->col[j];
         if (col->kind == GLP_CV)
            col->mipx = col->prim;
         else if (col->kind == GLP_IV)
         {  /* value of the integer column must be integral */
            col->mipx = floor(col->prim + 0.5);
         }
         else
            xassert(col != col);
      }
      T->sol_cnt++;
      return;
}

/***********************************************************************
*  fix_by_red_cost - fix non-basic integer columns by reduced costs
*
*  This routine fixes some non-basic integer columns if their reduced
*  costs indicate that increasing (decreasing) the column at least by
*  one involves the objective value becoming worse than the incumbent
*  objective value. */

static void fix_by_red_cost(glp_tree *T)
{     glp_prob *mip = T->mip;
      int j, stat, fixed = 0;
      double obj, lb, ub, dj;
      /* the global bound must exist */
      xassert(T->mip->mip_stat == GLP_FEAS);
      /* basic solution of LP relaxation must be optimal */
      xassert(mip->pbs_stat == GLP_FEAS && mip->dbs_stat == GLP_FEAS);
      /* determine the objective function value */
      obj = mip->obj_val;
      /* walk through the column list */
      for (j = 1; j <= mip->n; j++)
      {  GLPCOL *col = mip->col[j];
         /* if the column is not integer, skip it */
         if (col->kind != GLP_IV) continue;
         /* obtain bounds of j-th column */
         lb = col->lb, ub = col->ub;
         /* and determine its status and reduced cost */
         stat = col->stat, dj = col->dual;
         /* analyze the reduced cost */
         switch (mip->dir)
         {  case GLP_MIN:
               /* minimization */
               if (stat == GLP_NL)
               {  /* j-th column is non-basic on its lower bound */
                  if (dj < 0.0) dj = 0.0;
                  if (obj + dj >= mip->mip_obj)
                     glp_set_col_bnds(mip, j, GLP_FX, lb, lb), fixed++;
               }
               else if (stat == GLP_NU)
               {  /* j-th column is non-basic on its upper bound */
                  if (dj > 0.0) dj = 0.0;
                  if (obj - dj >= mip->mip_obj)
                     glp_set_col_bnds(mip, j, GLP_FX, ub, ub), fixed++;
               }
               break;
            case GLP_MAX:
               /* maximization */
               if (stat == GLP_NL)
               {  /* j-th column is non-basic on its lower bound */
                  if (dj > 0.0) dj = 0.0;
                  if (obj + dj <= mip->mip_obj)
                     glp_set_col_bnds(mip, j, GLP_FX, lb, lb), fixed++;
               }
               else if (stat == GLP_NU)
               {  /* j-th column is non-basic on its upper bound */
                  if (dj < 0.0) dj = 0.0;
                  if (obj - dj <= mip->mip_obj)
                     glp_set_col_bnds(mip, j, GLP_FX, ub, ub), fixed++;
               }
               break;
            default:
               xassert(T != T);
         }
      }
      if (T->parm->msg_lev >= GLP_MSG_DBG)
      {  if (fixed == 0)
            /* nothing to say */;
         else if (fixed == 1)
            xprintf("One column has been fixed by reduced cost\n");
         else
            xprintf("%d columns have been fixed by reduced costs\n",
               fixed);
      }
      /* fixing non-basic columns on their current bounds does not
         change the basic solution */
      xassert(mip->pbs_stat == GLP_FEAS && mip->dbs_stat == GLP_FEAS);
      return;
}

/***********************************************************************
*  branch_on - perform branching on specified variable
*
*  This routine performs branching on j-th column (structural variable)
*  of the current subproblem. The specified column must be of integer
*  kind and must have a fractional value in optimal basic solution of
*  LP relaxation of the current subproblem (i.e. only columns for which
*  the flag non_int[j] is set are valid candidates to branch on).
*
*  Let x be j-th structural variable, and beta be its primal fractional
*  value in the current basic solution. Branching on j-th variable is
*  dividing the current subproblem into two new subproblems, which are
*  identical to the current subproblem with the following exception: in
*  the first subproblem that begins the down-branch x has a new upper
*  bound x <= floor(beta), and in the second subproblem that begins the
*  up-branch x has a new lower bound x >= ceil(beta).
*
*  Depending on estimation of local bounds for down- and up-branches
*  this routine returns the following:
*
*  0 - both branches have been created;
*  1 - one branch is hopeless and has been pruned, so now the current
*      subproblem is other branch;
*  2 - both branches are hopeless and have been pruned; new subproblem
*      selection is needed to continue the search. */

static int branch_on(glp_tree *T, int j, int next)
{     glp_prob *mip = T->mip;
      IOSNPD *node;
      int m = mip->m;
      int n = mip->n;
      int type, dn_type, up_type, dn_bad, up_bad, p, ret, clone[1+2];
      double lb, ub, beta, new_ub, new_lb, dn_lp, up_lp, dn_bnd, up_bnd;
      /* determine bounds and value of x[j] in optimal solution to LP
         relaxation of the current subproblem */
      xassert(1 <= j && j <= n);
      type = mip->col[j]->type;
      lb = mip->col[j]->lb;
      ub = mip->col[j]->ub;
      beta = mip->col[j]->prim;
      /* determine new bounds of x[j] for down- and up-branches */
      new_ub = floor(beta);
      new_lb = ceil(beta);
      switch (type)
      {  case GLP_FR:
            dn_type = GLP_UP;
            up_type = GLP_LO;
            break;
         case GLP_LO:
            xassert(lb <= new_ub);
            dn_type = (lb == new_ub ? GLP_FX : GLP_DB);
            xassert(lb + 1.0 <= new_lb);
            up_type = GLP_LO;
            break;
         case GLP_UP:
            xassert(new_ub <= ub - 1.0);
            dn_type = GLP_UP;
            xassert(new_lb <= ub);
            up_type = (new_lb == ub ? GLP_FX : GLP_DB);
            break;
         case GLP_DB:
            xassert(lb <= new_ub && new_ub <= ub - 1.0);
            dn_type = (lb == new_ub ? GLP_FX : GLP_DB);
            xassert(lb + 1.0 <= new_lb && new_lb <= ub);
            up_type = (new_lb == ub ? GLP_FX : GLP_DB);
            break;
         default:
            xassert(type != type);
      }
      /* compute local bounds to LP relaxation for both branches */
      ios_eval_degrad(T, j, &dn_lp, &up_lp);
      /* and improve them by rounding */
      dn_bnd = ios_round_bound(T, dn_lp);
      up_bnd = ios_round_bound(T, up_lp);
      /* check local bounds for down- and up-branches */
      dn_bad = !ios_is_hopeful(T, dn_bnd);
      up_bad = !ios_is_hopeful(T, up_bnd);
      if (dn_bad && up_bad)
      {  if (T->parm->msg_lev >= GLP_MSG_DBG)
            xprintf("Both down- and up-branches are hopeless\n");
         ret = 2;
         goto done;
      }
      else if (up_bad)
      {  if (T->parm->msg_lev >= GLP_MSG_DBG)
            xprintf("Up-branch is hopeless\n");
         glp_set_col_bnds(mip, j, dn_type, lb, new_ub);
         T->curr->lp_obj = dn_lp;
         if (mip->dir == GLP_MIN)
         {  if (T->curr->bound < dn_bnd)
                T->curr->bound = dn_bnd;
         }
         else if (mip->dir == GLP_MAX)
         {  if (T->curr->bound > dn_bnd)
                T->curr->bound = dn_bnd;
         }
         else
            xassert(mip != mip);
         ret = 1;
         goto done;
      }
      else if (dn_bad)
      {  if (T->parm->msg_lev >= GLP_MSG_DBG)
            xprintf("Down-branch is hopeless\n");
         glp_set_col_bnds(mip, j, up_type, new_lb, ub);
         T->curr->lp_obj = up_lp;
         if (mip->dir == GLP_MIN)
         {  if (T->curr->bound < up_bnd)
                T->curr->bound = up_bnd;
         }
         else if (mip->dir == GLP_MAX)
         {  if (T->curr->bound > up_bnd)
                T->curr->bound = up_bnd;
         }
         else
            xassert(mip != mip);
         ret = 1;
         goto done;
      }
      /* both down- and up-branches seem to be hopeful */
      if (T->parm->msg_lev >= GLP_MSG_DBG)
         xprintf("Branching on column %d, primal value is %.9e\n",
            j, beta);
      /* determine the reference number of the current subproblem */
      xassert(T->curr != NULL);
      p = T->curr->p;
      T->curr->br_var = j;
      T->curr->br_val = beta;
      /* freeze the current subproblem */
      ios_freeze_node(T);
      /* create two clones of the current subproblem; the first clone
         begins the down-branch, the second one begins the up-branch */
      ios_clone_node(T, p, 2, clone);
      if (T->parm->msg_lev >= GLP_MSG_DBG)
         xprintf("Node %d begins down branch, node %d begins up branch "
            "\n", clone[1], clone[2]);
      /* set new upper bound of j-th column in the down-branch */
      node = T->slot[clone[1]].node;
      xassert(node != NULL);
      xassert(node->up != NULL);
      xassert(node->b_ptr == NULL);
      node->b_ptr = dmp_get_atom(T->pool, sizeof(IOSBND));
      node->b_ptr->k = m + j;
      node->b_ptr->type = (unsigned char)dn_type;
      node->b_ptr->lb = lb;
      node->b_ptr->ub = new_ub;
      node->b_ptr->next = NULL;
      node->lp_obj = dn_lp;
      if (mip->dir == GLP_MIN)
      {  if (node->bound < dn_bnd)
             node->bound = dn_bnd;
      }
      else if (mip->dir == GLP_MAX)
      {  if (node->bound > dn_bnd)
             node->bound = dn_bnd;
      }
      else
         xassert(mip != mip);
      /* set new lower bound of j-th column in the up-branch */
      node = T->slot[clone[2]].node;
      xassert(node != NULL);
      xassert(node->up != NULL);
      xassert(node->b_ptr == NULL);
      node->b_ptr = dmp_get_atom(T->pool, sizeof(IOSBND));
      node->b_ptr->k = m + j;
      node->b_ptr->type = (unsigned char)up_type;
      node->b_ptr->lb = new_lb;
      node->b_ptr->ub = ub;
      node->b_ptr->next = NULL;
      node->lp_obj = up_lp;
      if (mip->dir == GLP_MIN)
      {  if (node->bound < up_bnd)
             node->bound = up_bnd;
      }
      else if (mip->dir == GLP_MAX)
      {  if (node->bound > up_bnd)
             node->bound = up_bnd;
      }
      else
         xassert(mip != mip);
      /* suggest the subproblem to be solved next */
      xassert(T->child == 0);
      if (next == GLP_NO_BRNCH)
         T->child = 0;
      else if (next == GLP_DN_BRNCH)
         T->child = clone[1];
      else if (next == GLP_UP_BRNCH)
         T->child = clone[2];
      else
         xassert(next != next);
      ret = 0;
done: return ret;
}

/***********************************************************************
*  cleanup_the_tree - prune hopeless branches from the tree
*
*  This routine walks through the active list and checks the local
*  bound for every active subproblem. If the local bound indicates that
*  the subproblem cannot have integer optimal solution better than the
*  incumbent objective value, the routine deletes such subproblem that,
*  in turn, involves pruning the corresponding branch of the tree. */

static void cleanup_the_tree(glp_tree *T)
{     IOSNPD *node, *next_node;
      int count = 0;
      /* the global bound must exist */
      xassert(T->mip->mip_stat == GLP_FEAS);
      /* walk through the list of active subproblems */
      for (node = T->head; node != NULL; node = next_node)
      {  /* deleting some active problem node may involve deleting its
            parents recursively; however, all its parents being created
            *before* it are always *precede* it in the node list, so
            the next problem node is never affected by such deletion */
         next_node = node->next;
         /* if the branch is hopeless, prune it */
         if (!is_branch_hopeful(T, node->p))
            ios_delete_node(T, node->p), count++;
      }
      if (T->parm->msg_lev >= GLP_MSG_DBG)
      {  if (count == 1)
            xprintf("One hopeless branch has been pruned\n");
         else if (count > 1)
            xprintf("%d hopeless branches have been pruned\n", count);
      }
      return;
}

/**********************************************************************/

static void generate_cuts(glp_tree *T)
{     /* generate generic cuts with built-in generators */
      if (!(T->parm->mir_cuts == GLP_ON ||
            T->parm->gmi_cuts == GLP_ON ||
            T->parm->cov_cuts == GLP_ON ||
            T->parm->clq_cuts == GLP_ON)) goto done;
#if 1 /* 20/IX-2008 */
      {  int i, max_cuts, added_cuts;
         max_cuts = T->n;
         if (max_cuts < 1000) max_cuts = 1000;
         added_cuts = 0;
         for (i = T->orig_m+1; i <= T->mip->m; i++)
         {  if (T->mip->row[i]->origin == GLP_RF_CUT)
               added_cuts++;
         }
         /* xprintf("added_cuts = %d\n", added_cuts); */
         if (added_cuts >= max_cuts) goto done;
      }
#endif
      /* generate and add to POOL all cuts violated by x* */
      if (T->parm->gmi_cuts == GLP_ON)
      {  if (T->curr->changed < 5)
            ios_gmi_gen(T);
      }
      if (T->parm->mir_cuts == GLP_ON)
      {  xassert(T->mir_gen != NULL);
         ios_mir_gen(T, T->mir_gen);
      }
      if (T->parm->cov_cuts == GLP_ON)
      {  /* cover cuts works well along with mir cuts */
         /*if (T->round <= 5)*/
            ios_cov_gen(T);
      }
      if (T->parm->clq_cuts == GLP_ON)
      {  if (T->clq_gen != NULL)
         {  if (T->curr->level == 0 && T->curr->changed < 50 ||
                T->curr->level >  0 && T->curr->changed < 5)
               ios_clq_gen(T, T->clq_gen);
         }
      }
done: return;
}

/**********************************************************************/

static void remove_cuts(glp_tree *T)
{     /* remove inactive cuts (some valueable globally valid cut might
         be saved in the global cut pool) */
      int i, cnt = 0, *num = NULL;
      xassert(T->curr != NULL);
      for (i = T->orig_m+1; i <= T->mip->m; i++)
      {  if (T->mip->row[i]->origin == GLP_RF_CUT &&
             T->mip->row[i]->level == T->curr->level &&
             T->mip->row[i]->stat == GLP_BS)
         {  if (num == NULL)
               num = xcalloc(1+T->mip->m, sizeof(int));
            num[++cnt] = i;
         }
      }
      if (cnt > 0)
      {  glp_del_rows(T->mip, cnt, num);
#if 0
         xprintf("%d inactive cut(s) removed\n", cnt);
#endif
         xfree(num);
         xassert(glp_factorize(T->mip) == 0);
      }
      return;
}

/**********************************************************************/

static void display_cut_info(glp_tree *T)
{     glp_prob *mip = T->mip;
      int i, gmi = 0, mir = 0, cov = 0, clq = 0, app = 0;
      for (i = mip->m; i > 0; i--)
      {  GLPROW *row;
         row = mip->row[i];
         /* if (row->level < T->curr->level) break; */
         if (row->origin == GLP_RF_CUT)
         {  if (row->klass == GLP_RF_GMI)
               gmi++;
            else if (row->klass == GLP_RF_MIR)
               mir++;
            else if (row->klass == GLP_RF_COV)
               cov++;
            else if (row->klass == GLP_RF_CLQ)
               clq++;
            else
               app++;
         }
      }
      xassert(T->curr != NULL);
      if (gmi + mir + cov + clq + app > 0)
      {  xprintf("Cuts on level %d:", T->curr->level);
         if (gmi > 0) xprintf(" gmi = %d;", gmi);
         if (mir > 0) xprintf(" mir = %d;", mir);
         if (cov > 0) xprintf(" cov = %d;", cov);
         if (clq > 0) xprintf(" clq = %d;", clq);
         if (app > 0) xprintf(" app = %d;", app);
         xprintf("\n");
      }
      return;
}

/***********************************************************************
*  NAME
*
*  ios_driver - branch-and-cut driver
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  int ios_driver(glp_tree *T);
*
*  DESCRIPTION
*
*  The routine ios_driver is a branch-and-cut driver. It controls the
*  MIP solution process.
*
*  RETURNS
*
*  0  The MIP problem instance has been successfully solved. This code
*     does not necessarily mean that the solver has found optimal
*     solution. It only means that the solution process was successful.
*
*  GLP_EFAIL
*     The search was prematurely terminated due to the solver failure.
*
*  GLP_EMIPGAP
*     The search was prematurely terminated, because the relative mip
*     gap tolerance has been reached.
*
*  GLP_ETMLIM
*     The search was prematurely terminated, because the time limit has
*     been exceeded.
*
*  GLP_ESTOP
*     The search was prematurely terminated by application. */

int ios_driver(glp_tree *T)
{     int p, curr_p, p_stat, d_stat, ret;
#if 1 /* carry out to glp_tree */
      int pred_p = 0;
      /* if the current subproblem has been just created due to
         branching, pred_p is the reference number of its parent
         subproblem, otherwise pred_p is zero */
#endif
      glp_long ttt = T->tm_beg;
#if 0
      ((glp_iocp *)T->parm)->msg_lev = GLP_MSG_DBG;
#endif
      /* on entry to the B&B driver it is assumed that the active list
         contains the only active (i.e. root) subproblem, which is the
         original MIP problem to be solved */
loop: /* main loop starts here */
      /* at this point the current subproblem does not exist */
      xassert(T->curr == NULL);
      /* if the active list is empty, the search is finished */
      if (T->head == NULL)
      {  if (T->parm->msg_lev >= GLP_MSG_DBG)
            xprintf("Active list is empty!\n");
         xassert(dmp_in_use(T->pool).lo == 0);
         ret = 0;
         goto done;
      }
      /* select some active subproblem to continue the search */
      xassert(T->next_p == 0);
      /* let the application program select subproblem */
      if (T->parm->cb_func != NULL)
      {  xassert(T->reason == 0);
         T->reason = GLP_ISELECT;
         T->parm->cb_func(T, T->parm->cb_info);
         T->reason = 0;
         if (T->stop)
         {  ret = GLP_ESTOP;
            goto done;
         }
      }
      if (T->next_p != 0)
      {  /* the application program has selected something */
         ;
      }
      else if (T->a_cnt == 1)
      {  /* the only active subproblem exists, so select it */
         xassert(T->head->next == NULL);
         T->next_p = T->head->p;
      }
      else if (T->child != 0)
      {  /* select one of branching childs suggested by the branching
            heuristic */
         T->next_p = T->child;
      }
      else
      {  /* select active subproblem as specified by the backtracking
            technique option */
         T->next_p = ios_choose_node(T);
      }
      /* the active subproblem just selected becomes current */
      ios_revive_node(T, T->next_p);
      T->next_p = T->child = 0;
      /* invalidate pred_p, if it is not the reference number of the
         parent of the current subproblem */
      if (T->curr->up != NULL && T->curr->up->p != pred_p) pred_p = 0;
      /* determine the reference number of the current subproblem */
      p = T->curr->p;
      if (T->parm->msg_lev >= GLP_MSG_DBG)
      {  xprintf("-----------------------------------------------------"
            "-------------------\n");
         xprintf("Processing node %d at level %d\n", p, T->curr->level);
      }
      /* if it is the root subproblem, initialize cut generators */
      if (p == 1)
      {  if (T->parm->gmi_cuts == GLP_ON)
         {  if (T->parm->msg_lev >= GLP_MSG_ALL)
               xprintf("Gomory's cuts enabled\n");
         }
         if (T->parm->mir_cuts == GLP_ON)
         {  if (T->parm->msg_lev >= GLP_MSG_ALL)
               xprintf("MIR cuts enabled\n");
            xassert(T->mir_gen == NULL);
            T->mir_gen = ios_mir_init(T);
         }
         if (T->parm->cov_cuts == GLP_ON)
         {  if (T->parm->msg_lev >= GLP_MSG_ALL)
               xprintf("Cover cuts enabled\n");
         }
         if (T->parm->clq_cuts == GLP_ON)
         {  xassert(T->clq_gen == NULL);
            if (T->parm->msg_lev >= GLP_MSG_ALL)
               xprintf("Clique cuts enabled\n");
            T->clq_gen = ios_clq_init(T);
         }
      }
more: /* minor loop starts here */
      /* at this point the current subproblem needs either to be solved
         for the first time or re-optimized due to reformulation */
      /* display current progress of the search */
      if (T->parm->msg_lev >= GLP_MSG_DBG ||
          T->parm->msg_lev >= GLP_MSG_ON &&
        (double)(T->parm->out_frq - 1) <=
            1000.0 * xdifftime(xtime(), T->tm_lag))
         show_progress(T, 0);
      if (T->parm->msg_lev >= GLP_MSG_ALL &&
            xdifftime(xtime(), ttt) >= 60.0)
      {  glp_long total;
         glp_mem_usage(NULL, NULL, &total, NULL);
         xprintf("Time used: %.1f secs.  Memory used: %.1f Mb.\n",
            xdifftime(xtime(), T->tm_beg), xltod(total) / 1048576.0);
         ttt = xtime();
      }
      /* check the mip gap */
      if (T->parm->mip_gap > 0.0 &&
          ios_relative_gap(T) <= T->parm->mip_gap)
      {  if (T->parm->msg_lev >= GLP_MSG_DBG)
            xprintf("Relative gap tolerance reached; search terminated "
               "\n");
         ret = GLP_EMIPGAP;
         goto done;
      }
      /* check if the time limit has been exhausted */
      if (T->parm->tm_lim < INT_MAX &&
         (double)(T->parm->tm_lim - 1) <=
         1000.0 * xdifftime(xtime(), T->tm_beg))
      {  if (T->parm->msg_lev >= GLP_MSG_DBG)
            xprintf("Time limit exhausted; search terminated\n");
         ret = GLP_ETMLIM;
         goto done;
      }
      /* let the application program preprocess the subproblem */
      if (T->parm->cb_func != NULL)
      {  xassert(T->reason == 0);
         T->reason = GLP_IPREPRO;
         T->parm->cb_func(T, T->parm->cb_info);
         T->reason = 0;
         if (T->stop)
         {  ret = GLP_ESTOP;
            goto done;
         }
      }
      /* perform basic preprocessing */
      if (T->parm->pp_tech == GLP_PP_NONE)
         ;
      else if (T->parm->pp_tech == GLP_PP_ROOT)
      {  if (T->curr->level == 0)
         {  if (ios_preprocess_node(T, 100))
               goto fath;
         }
      }
      else if (T->parm->pp_tech == GLP_PP_ALL)
      {  if (ios_preprocess_node(T, T->curr->level == 0 ? 100 : 10))
            goto fath;
      }
      else
         xassert(T != T);
      /* preprocessing may improve the global bound */
      if (!is_branch_hopeful(T, p))
      {  xprintf("*** not tested yet ***\n");
         goto fath;
      }
      /* solve LP relaxation of the current subproblem */
      if (T->parm->msg_lev >= GLP_MSG_DBG)
         xprintf("Solving LP relaxation...\n");
      ret = ios_solve_node(T);
      if (!(ret == 0 || ret == GLP_EOBJLL || ret == GLP_EOBJUL))
      {  if (T->parm->msg_lev >= GLP_MSG_ERR)
            xprintf("ios_driver: unable to solve current LP relaxation;"
               " glp_simplex returned %d\n", ret);
         ret = GLP_EFAIL;
         goto done;
      }
      /* analyze status of the basic solution to LP relaxation found */
      p_stat = T->mip->pbs_stat;
      d_stat = T->mip->dbs_stat;
      if (p_stat == GLP_FEAS && d_stat == GLP_FEAS)
      {  /* LP relaxation has optimal solution */
         if (T->parm->msg_lev >= GLP_MSG_DBG)
            xprintf("Found optimal solution to LP relaxation\n");
      }
      else if (d_stat == GLP_NOFEAS)
      {  /* LP relaxation has no dual feasible solution */
         /* since the current subproblem cannot have a larger feasible
            region than its parent, there is something wrong */
         if (T->parm->msg_lev >= GLP_MSG_ERR)
            xprintf("ios_driver: current LP relaxation has no dual feas"
               "ible solution\n");
         ret = GLP_EFAIL;
         goto done;
      }
      else if (p_stat == GLP_INFEAS && d_stat == GLP_FEAS)
      {  /* LP relaxation has no primal solution which is better than
            the incumbent objective value */
         xassert(T->mip->mip_stat == GLP_FEAS);
         if (T->parm->msg_lev >= GLP_MSG_DBG)
            xprintf("LP relaxation has no solution better than incumben"
               "t objective value\n");
         /* prune the branch */
         goto fath;
      }
      else if (p_stat == GLP_NOFEAS)
      {  /* LP relaxation has no primal feasible solution */
         if (T->parm->msg_lev >= GLP_MSG_DBG)
            xprintf("LP relaxation has no feasible solution\n");
         /* prune the branch */
         goto fath;
      }
      else
      {  /* other cases cannot appear */
         xassert(T->mip != T->mip);
      }
      /* at this point basic solution to LP relaxation of the current
         subproblem is optimal */
      xassert(p_stat == GLP_FEAS && d_stat == GLP_FEAS);
      xassert(T->curr != NULL);
      T->curr->lp_obj = T->mip->obj_val;
      /* thus, it defines a local bound to integer optimal solution of
         the current subproblem */
      {  double bound = T->mip->obj_val;
         /* some local bound to the current subproblem could be already
            set before, so we should only improve it */
         bound = ios_round_bound(T, bound);
         if (T->mip->dir == GLP_MIN)
         {  if (T->curr->bound < bound)
               T->curr->bound = bound;
         }
         else if (T->mip->dir == GLP_MAX)
         {  if (T->curr->bound > bound)
               T->curr->bound = bound;
         }
         else
            xassert(T->mip != T->mip);
         if (T->parm->msg_lev >= GLP_MSG_DBG)
            xprintf("Local bound is %.9e\n", bound);
      }
      /* if the local bound indicates that integer optimal solution of
         the current subproblem cannot be better than the global bound,
         prune the branch */
      if (!is_branch_hopeful(T, p))
      {  if (T->parm->msg_lev >= GLP_MSG_DBG)
            xprintf("Current branch is hopeless and can be pruned\n");
         goto fath;
      }
      /* let the application program generate additional rows ("lazy"
         constraints) */
      xassert(T->reopt == 0);
      xassert(T->reinv == 0);
      if (T->parm->cb_func != NULL)
      {  xassert(T->reason == 0);
         T->reason = GLP_IROWGEN;
         T->parm->cb_func(T, T->parm->cb_info);
         T->reason = 0;
         if (T->stop)
         {  ret = GLP_ESTOP;
            goto done;
         }
         if (T->reopt)
         {  /* some rows were added; re-optimization is needed */
            T->reopt = T->reinv = 0;
            goto more;
         }
         if (T->reinv)
         {  /* no rows were added, however, some inactive rows were
               removed */
            T->reinv = 0;
            xassert(glp_factorize(T->mip) == 0);
         }
      }
      /* check if the basic solution is integer feasible */
      check_integrality(T);
      /* if the basic solution satisfies to all integrality conditions,
         it is a new, better integer feasible solution */
      if (T->curr->ii_cnt == 0)
      {  if (T->parm->msg_lev >= GLP_MSG_DBG)
            xprintf("New integer feasible solution found\n");
         if (T->parm->msg_lev >= GLP_MSG_ALL)
            display_cut_info(T);
         record_solution(T);
         if (T->parm->msg_lev >= GLP_MSG_ON)
            show_progress(T, 1);
         /* make the application program happy */
         if (T->parm->cb_func != NULL)
         {  xassert(T->reason == 0);
            T->reason = GLP_IBINGO;
            T->parm->cb_func(T, T->parm->cb_info);
            T->reason = 0;
            if (T->stop)
            {  ret = GLP_ESTOP;
               goto done;
            }
         }
         /* since the current subproblem has been fathomed, prune its
            branch */
         goto fath;
      }
      /* at this point basic solution to LP relaxation of the current
         subproblem is optimal, but integer infeasible */
      /* try to fix some non-basic structural variables of integer kind
         on their current bounds due to reduced costs */
      if (T->mip->mip_stat == GLP_FEAS)
         fix_by_red_cost(T);
      /* let the application program try to find some solution to the
         original MIP with a primal heuristic */
      if (T->parm->cb_func != NULL)
      {  xassert(T->reason == 0);
         T->reason = GLP_IHEUR;
         T->parm->cb_func(T, T->parm->cb_info);
         T->reason = 0;
         if (T->stop)
         {  ret = GLP_ESTOP;
            goto done;
         }
         /* check if the current branch became hopeless */
         if (!is_branch_hopeful(T, p))
         {  if (T->parm->msg_lev >= GLP_MSG_DBG)
               xprintf("Current branch became hopeless and can be prune"
                  "d\n");
            goto fath;
         }
      }
      /* try to find solution with the feasibility pump heuristic */
      if (T->parm->fp_heur)
      {  xassert(T->reason == 0);
         T->reason = GLP_IHEUR;
         ios_feas_pump(T);
         T->reason = 0;
         /* check if the current branch became hopeless */
         if (!is_branch_hopeful(T, p))
         {  if (T->parm->msg_lev >= GLP_MSG_DBG)
               xprintf("Current branch became hopeless and can be prune"
                  "d\n");
            goto fath;
         }
      }
      /* it's time to generate cutting planes */
      xassert(T->local != NULL);
      xassert(T->local->size == 0);
      /* let the application program generate some cuts; note that it
         can add cuts either to the local cut pool or directly to the
         current subproblem */
      if (T->parm->cb_func != NULL)
      {  xassert(T->reason == 0);
         T->reason = GLP_ICUTGEN;
         T->parm->cb_func(T, T->parm->cb_info);
         T->reason = 0;
         if (T->stop)
         {  ret = GLP_ESTOP;
            goto done;
         }
      }
      /* try to generate generic cuts with built-in generators
         (as suggested by Matteo Fischetti et al. the built-in cuts
         are not generated at each branching node; an intense attempt
         of generating new cuts is only made at the root node, and then
         a moderate effort is spent after each backtracking step) */
      if (T->curr->level == 0 || pred_p == 0)
      {  xassert(T->reason == 0);
         T->reason = GLP_ICUTGEN;
         generate_cuts(T);
         T->reason = 0;
      }
      /* if the local cut pool is not empty, select useful cuts and add
         them to the current subproblem */
      if (T->local->size > 0)
      {  xassert(T->reason == 0);
         T->reason = GLP_ICUTGEN;
         ios_process_cuts(T);
         T->reason = 0;
      }
      /* clear the local cut pool */
      ios_clear_pool(T, T->local);
      /* perform re-optimization, if necessary */
      if (T->reopt)
      {  T->reopt = 0;
         T->curr->changed++;
         goto more;
      }
      /* no cuts were generated; remove inactive cuts */
      remove_cuts(T);
      if (T->parm->msg_lev >= GLP_MSG_ALL && T->curr->level == 0)
         display_cut_info(T);
      /* update history information used on pseudocost branching */
      if (T->pcost != NULL) ios_pcost_update(T);
      /* it's time to perform branching */
      xassert(T->br_var == 0);
      xassert(T->br_sel == 0);
      /* let the application program choose variable to branch on */
      if (T->parm->cb_func != NULL)
      {  xassert(T->reason == 0);
         xassert(T->br_var == 0);
         xassert(T->br_sel == 0);
         T->reason = GLP_IBRANCH;
         T->parm->cb_func(T, T->parm->cb_info);
         T->reason = 0;
         if (T->stop)
         {  ret = GLP_ESTOP;
            goto done;
         }
      }
      /* if nothing has been chosen, choose some variable as specified
         by the branching technique option */
      if (T->br_var == 0)
         T->br_var = ios_choose_var(T, &T->br_sel);
      /* perform actual branching */
      curr_p = T->curr->p;
      ret = branch_on(T, T->br_var, T->br_sel);
      T->br_var = T->br_sel = 0;
      if (ret == 0)
      {  /* both branches have been created */
         pred_p = curr_p;
         goto loop;
      }
      else if (ret == 1)
      {  /* one branch is hopeless and has been pruned, so now the
            current subproblem is other branch */
         /* the current subproblem should be considered as a new one,
            since one bound of the branching variable was changed */
         T->curr->solved = T->curr->changed = 0;
         goto more;
      }
      else if (ret == 2)
      {  /* both branches are hopeless and have been pruned; new
            subproblem selection is needed to continue the search */
         goto fath;
      }
      else
         xassert(ret != ret);
fath: /* the current subproblem has been fathomed */
      if (T->parm->msg_lev >= GLP_MSG_DBG)
         xprintf("Node %d fathomed\n", p);
      /* freeze the current subproblem */
      ios_freeze_node(T);
      /* and prune the corresponding branch of the tree */
      ios_delete_node(T, p);
      /* if a new integer feasible solution has just been found, other
         branches may become hopeless and therefore must be pruned */
      if (T->mip->mip_stat == GLP_FEAS) cleanup_the_tree(T);
      /* new subproblem selection is needed due to backtracking */
      pred_p = 0;
      goto loop;
done: /* display progress of the search on exit from the solver */
      if (T->parm->msg_lev >= GLP_MSG_ON)
         show_progress(T, 0);
      if (T->mir_gen != NULL)
         ios_mir_term(T->mir_gen), T->mir_gen = NULL;
      if (T->clq_gen != NULL)
         ios_clq_term(T->clq_gen), T->clq_gen = NULL;
      /* return to the calling program */
      return ret;
}

/* eof */
