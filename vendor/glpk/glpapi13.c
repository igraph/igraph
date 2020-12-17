/* glpapi13.c (branch-and-bound interface routines) */

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
*  glp_ios_reason - determine reason for calling the callback routine
*
*  SYNOPSIS
*
*  glp_ios_reason(glp_tree *tree);
*
*  RETURNS
*
*  The routine glp_ios_reason returns a code, which indicates why the
*  user-defined callback routine is being called. */

int glp_ios_reason(glp_tree *tree)
{     return
         tree->reason;
}

/***********************************************************************
*  NAME
*
*  glp_ios_get_prob - access the problem object
*
*  SYNOPSIS
*
*  glp_prob *glp_ios_get_prob(glp_tree *tree);
*
*  DESCRIPTION
*
*  The routine glp_ios_get_prob can be called from the user-defined
*  callback routine to access the problem object, which is used by the
*  MIP solver. It is the original problem object passed to the routine
*  glp_intopt if the MIP presolver is not used; otherwise it is an
*  internal problem object built by the presolver. If the current
*  subproblem exists, LP segment of the problem object corresponds to
*  its LP relaxation.
*
*  RETURNS
*
*  The routine glp_ios_get_prob returns a pointer to the problem object
*  used by the MIP solver. */

glp_prob *glp_ios_get_prob(glp_tree *tree)
{     return
         tree->mip;
}

/***********************************************************************
*  NAME
*
*  glp_ios_tree_size - determine size of the branch-and-bound tree
*
*  SYNOPSIS
*
*  void glp_ios_tree_size(glp_tree *tree, int *a_cnt, int *n_cnt,
*     int *t_cnt);
*
*  DESCRIPTION
*
*  The routine glp_ios_tree_size stores the following three counts which
*  characterize the current size of the branch-and-bound tree:
*
*  a_cnt is the current number of active nodes, i.e. the current size of
*        the active list;
*
*  n_cnt is the current number of all (active and inactive) nodes;
*
*  t_cnt is the total number of nodes including those which have been
*        already removed from the tree. This count is increased whenever
*        a new node appears in the tree and never decreased.
*
*  If some of the parameters a_cnt, n_cnt, t_cnt is a null pointer, the
*  corresponding count is not stored. */

void glp_ios_tree_size(glp_tree *tree, int *a_cnt, int *n_cnt,
      int *t_cnt)
{     if (a_cnt != NULL) *a_cnt = tree->a_cnt;
      if (n_cnt != NULL) *n_cnt = tree->n_cnt;
      if (t_cnt != NULL) *t_cnt = tree->t_cnt;
      return;
}

/***********************************************************************
*  NAME
*
*  glp_ios_curr_node - determine current active subproblem
*
*  SYNOPSIS
*
*  int glp_ios_curr_node(glp_tree *tree);
*
*  RETURNS
*
*  The routine glp_ios_curr_node returns the reference number of the
*  current active subproblem. However, if the current subproblem does
*  not exist, the routine returns zero. */

int glp_ios_curr_node(glp_tree *tree)
{     IOSNPD *node;
      /* obtain pointer to the current subproblem */
      node = tree->curr;
      /* return its reference number */
      return node == NULL ? 0 : node->p;
}

/***********************************************************************
*  NAME
*
*  glp_ios_next_node - determine next active subproblem
*
*  SYNOPSIS
*
*  int glp_ios_next_node(glp_tree *tree, int p);
*
*  RETURNS
*
*  If the parameter p is zero, the routine glp_ios_next_node returns
*  the reference number of the first active subproblem. However, if the
*  tree is empty, zero is returned.
*
*  If the parameter p is not zero, it must specify the reference number
*  of some active subproblem, in which case the routine returns the
*  reference number of the next active subproblem. However, if there is
*  no next active subproblem in the list, zero is returned.
*
*  All subproblems in the active list are ordered chronologically, i.e.
*  subproblem A precedes subproblem B if A was created before B. */

int glp_ios_next_node(glp_tree *tree, int p)
{     IOSNPD *node;
      if (p == 0)
      {  /* obtain pointer to the first active subproblem */
         node = tree->head;
      }
      else
      {  /* obtain pointer to the specified subproblem */
         if (!(1 <= p && p <= tree->nslots))
err:        xerror("glp_ios_next_node: p = %d; invalid subproblem refer"
               "ence number\n", p);
         node = tree->slot[p].node;
         if (node == NULL) goto err;
         /* the specified subproblem must be active */
         if (node->count != 0)
            xerror("glp_ios_next_node: p = %d; subproblem not in the ac"
               "tive list\n", p);
         /* obtain pointer to the next active subproblem */
         node = node->next;
      }
      /* return the reference number */
      return node == NULL ? 0 : node->p;
}

/***********************************************************************
*  NAME
*
*  glp_ios_prev_node - determine previous active subproblem
*
*  SYNOPSIS
*
*  int glp_ios_prev_node(glp_tree *tree, int p);
*
*  RETURNS
*
*  If the parameter p is zero, the routine glp_ios_prev_node returns
*  the reference number of the last active subproblem. However, if the
*  tree is empty, zero is returned.
*
*  If the parameter p is not zero, it must specify the reference number
*  of some active subproblem, in which case the routine returns the
*  reference number of the previous active subproblem. However, if there
*  is no previous active subproblem in the list, zero is returned.
*
*  All subproblems in the active list are ordered chronologically, i.e.
*  subproblem A precedes subproblem B if A was created before B. */

int glp_ios_prev_node(glp_tree *tree, int p)
{     IOSNPD *node;
      if (p == 0)
      {  /* obtain pointer to the last active subproblem */
         node = tree->tail;
      }
      else
      {  /* obtain pointer to the specified subproblem */
         if (!(1 <= p && p <= tree->nslots))
err:        xerror("glp_ios_prev_node: p = %d; invalid subproblem refer"
               "ence number\n", p);
         node = tree->slot[p].node;
         if (node == NULL) goto err;
         /* the specified subproblem must be active */
         if (node->count != 0)
            xerror("glp_ios_prev_node: p = %d; subproblem not in the ac"
               "tive list\n", p);
         /* obtain pointer to the previous active subproblem */
         node = node->prev;
      }
      /* return the reference number */
      return node == NULL ? 0 : node->p;
}

/***********************************************************************
*  NAME
*
*  glp_ios_up_node - determine parent subproblem
*
*  SYNOPSIS
*
*  int glp_ios_up_node(glp_tree *tree, int p);
*
*  RETURNS
*
*  The parameter p must specify the reference number of some (active or
*  inactive) subproblem, in which case the routine iet_get_up_node
*  returns the reference number of its parent subproblem. However, if
*  the specified subproblem is the root of the tree and, therefore, has
*  no parent, the routine returns zero. */

int glp_ios_up_node(glp_tree *tree, int p)
{     IOSNPD *node;
      /* obtain pointer to the specified subproblem */
      if (!(1 <= p && p <= tree->nslots))
err:     xerror("glp_ios_up_node: p = %d; invalid subproblem reference "
            "number\n", p);
      node = tree->slot[p].node;
      if (node == NULL) goto err;
      /* obtain pointer to the parent subproblem */
      node = node->up;
      /* return the reference number */
      return node == NULL ? 0 : node->p;
}

/***********************************************************************
*  NAME
*
*  glp_ios_node_level - determine subproblem level
*
*  SYNOPSIS
*
*  int glp_ios_node_level(glp_tree *tree, int p);
*
*  RETURNS
*
*  The routine glp_ios_node_level returns the level of the subproblem,
*  whose reference number is p, in the branch-and-bound tree. (The root
*  subproblem has level 0, and the level of any other subproblem is the
*  level of its parent plus one.) */

int glp_ios_node_level(glp_tree *tree, int p)
{     IOSNPD *node;
      /* obtain pointer to the specified subproblem */
      if (!(1 <= p && p <= tree->nslots))
err:     xerror("glp_ios_node_level: p = %d; invalid subproblem referen"
            "ce number\n", p);
      node = tree->slot[p].node;
      if (node == NULL) goto err;
      /* return the node level */
      return node->level;
}

/***********************************************************************
*  NAME
*
*  glp_ios_node_bound - determine subproblem local bound
*
*  SYNOPSIS
*
*  double glp_ios_node_bound(glp_tree *tree, int p);
*
*  RETURNS
*
*  The routine glp_ios_node_bound returns the local bound for (active or
*  inactive) subproblem, whose reference number is p.
*
*  COMMENTS
*
*  The local bound for subproblem p is an lower (minimization) or upper
*  (maximization) bound for integer optimal solution to this subproblem
*  (not to the original problem). This bound is local in the sense that
*  only subproblems in the subtree rooted at node p cannot have better
*  integer feasible solutions.
*
*  On creating a subproblem (due to the branching step) its local bound
*  is inherited from its parent and then may get only stronger (never
*  weaker). For the root subproblem its local bound is initially set to
*  -DBL_MAX (minimization) or +DBL_MAX (maximization) and then improved
*  as the root LP relaxation has been solved.
*
*  Note that the local bound is not necessarily the optimal objective
*  value to corresponding LP relaxation; it may be stronger. */

double glp_ios_node_bound(glp_tree *tree, int p)
{     IOSNPD *node;
      /* obtain pointer to the specified subproblem */
      if (!(1 <= p && p <= tree->nslots))
err:     xerror("glp_ios_node_bound: p = %d; invalid subproblem referen"
            "ce number\n", p);
      node = tree->slot[p].node;
      if (node == NULL) goto err;
      /* return the node local bound */
      return node->bound;
}

/***********************************************************************
*  NAME
*
*  glp_ios_best_node - find active subproblem with best local bound
*
*  SYNOPSIS
*
*  int glp_ios_best_node(glp_tree *tree);
*
*  RETURNS
*
*  The routine glp_ios_best_node returns the reference number of the
*  active subproblem, whose local bound is best (i.e. smallest in case
*  of minimization or largest in case of maximization). However, if the
*  tree is empty, the routine returns zero.
*
*  COMMENTS
*
*  The best local bound is an lower (minimization) or upper
*  (maximization) bound for integer optimal solution to the original
*  MIP problem. */

int glp_ios_best_node(glp_tree *tree)
{     return
         ios_best_node(tree);
}

/***********************************************************************
*  NAME
*
*  glp_ios_mip_gap - compute relative MIP gap
*
*  SYNOPSIS
*
*  double glp_ios_mip_gap(glp_tree *tree);
*
*  DESCRIPTION
*
*  The routine glp_ios_mip_gap computes the relative MIP gap with the
*  following formula:
*
*     gap = |best_mip - best_bnd| / (|best_mip| + DBL_EPSILON),
*
*  where best_mip is the best integer feasible solution found so far,
*  best_bnd is the best (global) bound. If no integer feasible solution
*  has been found yet, gap is set to DBL_MAX.
*
*  RETURNS
*
*  The routine glp_ios_mip_gap returns the relative MIP gap. */

double glp_ios_mip_gap(glp_tree *tree)
{     return
         ios_relative_gap(tree);
}

/***********************************************************************
*  NAME
*
*  glp_ios_node_data - access subproblem application-specific data
*
*  SYNOPSIS
*
*  void *glp_ios_node_data(glp_tree *tree, int p);
*
*  DESCRIPTION
*
*  The routine glp_ios_node_data allows the application accessing a
*  memory block allocated for the subproblem (which may be active or
*  inactive), whose reference number is p.
*
*  The size of the block is defined by the control parameter cb_size
*  passed to the routine glp_intopt. The block is initialized by binary
*  zeros on creating corresponding subproblem, and its contents is kept
*  until the subproblem will be removed from the tree.
*
*  The application may use these memory blocks to store specific data
*  for each subproblem.
*
*  RETURNS
*
*  The routine glp_ios_node_data returns a pointer to the memory block
*  for the specified subproblem. Note that if cb_size = 0, the routine
*  returns a null pointer. */

void *glp_ios_node_data(glp_tree *tree, int p)
{     IOSNPD *node;
      /* obtain pointer to the specified subproblem */
      if (!(1 <= p && p <= tree->nslots))
err:     xerror("glp_ios_node_level: p = %d; invalid subproblem referen"
            "ce number\n", p);
      node = tree->slot[p].node;
      if (node == NULL) goto err;
      /* return pointer to the application-specific data */
      return node->data;
}

/***********************************************************************
*  NAME
*
*  glp_ios_row_attr - retrieve additional row attributes
*
*  SYNOPSIS
*
*  void glp_ios_row_attr(glp_tree *tree, int i, glp_attr *attr);
*
*  DESCRIPTION
*
*  The routine glp_ios_row_attr retrieves additional attributes of row
*  i and stores them in the structure glp_attr. */

void glp_ios_row_attr(glp_tree *tree, int i, glp_attr *attr)
{     GLPROW *row;
      if (!(1 <= i && i <= tree->mip->m))
         xerror("glp_ios_row_attr: i = %d; row number out of range\n",
            i);
      row = tree->mip->row[i];
      attr->level = row->level;
      attr->origin = row->origin;
      attr->klass = row->klass;
      return;
}

/**********************************************************************/

int glp_ios_pool_size(glp_tree *tree)
{     /* determine current size of the cut pool */
      if (tree->reason != GLP_ICUTGEN)
         xerror("glp_ios_pool_size: operation not allowed\n");
      xassert(tree->local != NULL);
      return tree->local->size;
}

/**********************************************************************/

int glp_ios_add_row(glp_tree *tree,
      const char *name, int klass, int flags, int len, const int ind[],
      const double val[], int type, double rhs)
{     /* add row (constraint) to the cut pool */
      int num;
      if (tree->reason != GLP_ICUTGEN)
         xerror("glp_ios_add_row: operation not allowed\n");
      xassert(tree->local != NULL);
      num = ios_add_row(tree, tree->local, name, klass, flags, len,
         ind, val, type, rhs);
      return num;
}

/**********************************************************************/

void glp_ios_del_row(glp_tree *tree, int i)
{     /* remove row (constraint) from the cut pool */
      if (tree->reason != GLP_ICUTGEN)
         xerror("glp_ios_del_row: operation not allowed\n");
      ios_del_row(tree, tree->local, i);
      return;
}

/**********************************************************************/

void glp_ios_clear_pool(glp_tree *tree)
{     /* remove all rows (constraints) from the cut pool */
      if (tree->reason != GLP_ICUTGEN)
         xerror("glp_ios_clear_pool: operation not allowed\n");
      ios_clear_pool(tree, tree->local);
      return;
}

/***********************************************************************
*  NAME
*
*  glp_ios_can_branch - check if can branch upon specified variable
*
*  SYNOPSIS
*
*  int glp_ios_can_branch(glp_tree *tree, int j);
*
*  RETURNS
*
*  If j-th variable (column) can be used to branch upon, the routine
*  glp_ios_can_branch returns non-zero, otherwise zero. */

int glp_ios_can_branch(glp_tree *tree, int j)
{     if (!(1 <= j && j <= tree->mip->n))
         xerror("glp_ios_can_branch: j = %d; column number out of range"
            "\n", j);
      return tree->non_int[j];
}

/***********************************************************************
*  NAME
*
*  glp_ios_branch_upon - choose variable to branch upon
*
*  SYNOPSIS
*
*  void glp_ios_branch_upon(glp_tree *tree, int j, int sel);
*
*  DESCRIPTION
*
*  The routine glp_ios_branch_upon can be called from the user-defined
*  callback routine in response to the reason GLP_IBRANCH to choose a
*  branching variable, whose ordinal number is j. Should note that only
*  variables, for which the routine glp_ios_can_branch returns non-zero,
*  can be used to branch upon.
*
*  The parameter sel is a flag that indicates which branch (subproblem)
*  should be selected next to continue the search:
*
*  GLP_DN_BRNCH - select down-branch;
*  GLP_UP_BRNCH - select up-branch;
*  GLP_NO_BRNCH - use general selection technique. */

void glp_ios_branch_upon(glp_tree *tree, int j, int sel)
{     if (!(1 <= j && j <= tree->mip->n))
         xerror("glp_ios_branch_upon: j = %d; column number out of rang"
            "e\n", j);
      if (!(sel == GLP_DN_BRNCH || sel == GLP_UP_BRNCH ||
            sel == GLP_NO_BRNCH))
         xerror("glp_ios_branch_upon: sel = %d: invalid branch selectio"
            "n flag\n", sel);
      if (!(tree->non_int[j]))
         xerror("glp_ios_branch_upon: j = %d; variable cannot be used t"
            "o branch upon\n", j);
      if (tree->br_var != 0)
         xerror("glp_ios_branch_upon: branching variable already chosen"
            "\n");
      tree->br_var = j;
      tree->br_sel = sel;
      return;
}

/***********************************************************************
*  NAME
*
*  glp_ios_select_node - select subproblem to continue the search
*
*  SYNOPSIS
*
*  void glp_ios_select_node(glp_tree *tree, int p);
*
*  DESCRIPTION
*
*  The routine glp_ios_select_node can be called from the user-defined
*  callback routine in response to the reason GLP_ISELECT to select an
*  active subproblem, whose reference number is p. The search will be
*  continued from the subproblem selected. */

void glp_ios_select_node(glp_tree *tree, int p)
{     IOSNPD *node;
      /* obtain pointer to the specified subproblem */
      if (!(1 <= p && p <= tree->nslots))
err:     xerror("glp_ios_select_node: p = %d; invalid subproblem refere"
            "nce number\n", p);
      node = tree->slot[p].node;
      if (node == NULL) goto err;
      /* the specified subproblem must be active */
      if (node->count != 0)
         xerror("glp_ios_select_node: p = %d; subproblem not in the act"
            "ive list\n", p);
      /* no subproblem must be selected yet */
      if (tree->next_p != 0)
         xerror("glp_ios_select_node: subproblem already selected\n");
      /* select the specified subproblem to continue the search */
      tree->next_p = p;
      return;
}

/***********************************************************************
*  NAME
*
*  glp_ios_heur_sol - provide solution found by heuristic
*
*  SYNOPSIS
*
*  int glp_ios_heur_sol(glp_tree *tree, const double x[]);
*
*  DESCRIPTION
*
*  The routine glp_ios_heur_sol can be called from the user-defined
*  callback routine in response to the reason GLP_IHEUR to provide an
*  integer feasible solution found by a primal heuristic.
*
*  Primal values of *all* variables (columns) found by the heuristic
*  should be placed in locations x[1], ..., x[n], where n is the number
*  of columns in the original problem object. Note that the routine
*  glp_ios_heur_sol *does not* check primal feasibility of the solution
*  provided.
*
*  Using the solution passed in the array x the routine computes value
*  of the objective function. If the objective value is better than the
*  best known integer feasible solution, the routine computes values of
*  auxiliary variables (rows) and stores all solution components in the
*  problem object.
*
*  RETURNS
*
*  If the provided solution is accepted, the routine glp_ios_heur_sol
*  returns zero. Otherwise, if the provided solution is rejected, the
*  routine returns non-zero. */

int glp_ios_heur_sol(glp_tree *tree, const double x[])
{     glp_prob *mip = tree->mip;
      int m = tree->orig_m;
      int n = tree->n;
      int i, j;
      double obj;
      xassert(mip->m >= m);
      xassert(mip->n == n);
      /* check values of integer variables and compute value of the
         objective function */
      obj = mip->c0;
      for (j = 1; j <= n; j++)
      {  GLPCOL *col = mip->col[j];
         if (col->kind == GLP_IV)
         {  /* provided value must be integral */
            if (x[j] != floor(x[j])) return 1;
         }
         obj += col->coef * x[j];
      }
      /* check if the provided solution is better than the best known
         integer feasible solution */
      if (mip->mip_stat == GLP_FEAS)
      {  switch (mip->dir)
         {  case GLP_MIN:
               if (obj >= tree->mip->mip_obj) return 1;
               break;
            case GLP_MAX:
               if (obj <= tree->mip->mip_obj) return 1;
               break;
            default:
               xassert(mip != mip);
         }
      }
      /* it is better; store it in the problem object */
      if (tree->parm->msg_lev >= GLP_MSG_ON)
         xprintf("Solution found by heuristic: %.12g\n", obj);
      mip->mip_stat = GLP_FEAS;
      mip->mip_obj = obj;
      for (j = 1; j <= n; j++)
         mip->col[j]->mipx = x[j];
      for (i = 1; i <= m; i++)
      {  GLPROW *row = mip->row[i];
         GLPAIJ *aij;
         row->mipx = 0.0;
         for (aij = row->ptr; aij != NULL; aij = aij->r_next)
            row->mipx += aij->val * aij->col->mipx;
      }
      return 0;
}

/***********************************************************************
*  NAME
*
*  glp_ios_terminate - terminate the solution process.
*
*  SYNOPSIS
*
*  void glp_ios_terminate(glp_tree *tree);
*
*  DESCRIPTION
*
*  The routine glp_ios_terminate sets a flag indicating that the MIP
*  solver should prematurely terminate the search. */

void glp_ios_terminate(glp_tree *tree)
{     if (tree->parm->msg_lev >= GLP_MSG_DBG)
         xprintf("The search is prematurely terminated due to applicati"
            "on request\n");
      tree->stop = 1;
      return;
}

/* eof */
