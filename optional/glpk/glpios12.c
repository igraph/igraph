/* glpios12.c (node selection heuristics) */

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
*  ios_choose_node - select subproblem to continue the search
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  int ios_choose_node(glp_tree *T);
*
*  DESCRIPTION
*
*  The routine ios_choose_node selects a subproblem from the active
*  list to continue the search. The choice depends on the backtracking
*  technique option.
*
*  RETURNS
*
*  The routine ios_choose_node return the reference number of the
*  subproblem selected. */

static int most_feas(glp_tree *T);
static int best_proj(glp_tree *T);
static int best_node(glp_tree *T);

int ios_choose_node(glp_tree *T)
{     int p;
      if (T->parm->bt_tech == GLP_BT_DFS)
      {  /* depth first search */
         xassert(T->tail != NULL);
         p = T->tail->p;
      }
      else if (T->parm->bt_tech == GLP_BT_BFS)
      {  /* breadth first search */
         xassert(T->head != NULL);
         p = T->head->p;
      }
      else if (T->parm->bt_tech == GLP_BT_BLB)
      {  /* select node with best local bound */
         p = best_node(T);
      }
      else if (T->parm->bt_tech == GLP_BT_BPH)
      {  if (T->mip->mip_stat == GLP_UNDEF)
         {  /* "most integer feasible" subproblem */
            p = most_feas(T);
         }
         else
         {  /* best projection heuristic */
            p = best_proj(T);
         }
      }
      else
         xassert(T != T);
      return p;
}

static int most_feas(glp_tree *T)
{     /* select subproblem whose parent has minimal sum of integer
         infeasibilities */
      IOSNPD *node;
      int p;
      double best;
      p = 0, best = DBL_MAX;
      for (node = T->head; node != NULL; node = node->next)
      {  xassert(node->up != NULL);
         if (best > node->up->ii_sum)
            p = node->p, best = node->up->ii_sum;
      }
      return p;
}

static int best_proj(glp_tree *T)
{     /* select subproblem using the best projection heuristic */
      IOSNPD *root, *node;
      int p;
      double best, deg, obj;
      /* the global bound must exist */
      xassert(T->mip->mip_stat == GLP_FEAS);
      /* obtain pointer to the root node, which must exist */
      root = T->slot[1].node;
      xassert(root != NULL);
      /* deg estimates degradation of the objective function per unit
         of the sum of integer infeasibilities */
      xassert(root->ii_sum > 0.0);
      deg = (T->mip->mip_obj - root->bound) / root->ii_sum;
      /* nothing has been selected so far */
      p = 0, best = DBL_MAX;
      /* walk through the list of active subproblems */
      for (node = T->head; node != NULL; node = node->next)
      {  xassert(node->up != NULL);
         /* obj estimates optimal objective value if the sum of integer
            infeasibilities were zero */
         obj = node->up->bound + deg * node->up->ii_sum;
         if (T->mip->dir == GLP_MAX) obj = - obj;
         /* select the subproblem which has the best estimated optimal
            objective value */
         if (best > obj) p = node->p, best = obj;
      }
      return p;
}

static int best_node(glp_tree *T)
{     /* select subproblem with best local bound */
      IOSNPD *node, *best = NULL;
      double bound, eps;
      switch (T->mip->dir)
      {  case GLP_MIN:
            bound = +DBL_MAX;
            for (node = T->head; node != NULL; node = node->next)
               if (bound > node->bound) bound = node->bound;
            xassert(bound != +DBL_MAX);
            eps = 0.001 * (1.0 + fabs(bound));
            for (node = T->head; node != NULL; node = node->next)
            {  if (node->bound <= bound + eps)
               {  xassert(node->up != NULL);
                  if (best == NULL ||
#if 1
                  best->up->ii_sum > node->up->ii_sum) best = node;
#else
                  best->lp_obj > node->lp_obj) best = node;
#endif
               }
            }
            break;
         case GLP_MAX:
            bound = -DBL_MAX;
            for (node = T->head; node != NULL; node = node->next)
               if (bound < node->bound) bound = node->bound;
            xassert(bound != -DBL_MAX);
            eps = 0.001 * (1.0 + fabs(bound));
            for (node = T->head; node != NULL; node = node->next)
            {  if (node->bound >= bound - eps)
               {  xassert(node->up != NULL);
                  if (best == NULL ||
#if 1
                  best->up->ii_sum > node->up->ii_sum) best = node;
#else
                  best->lp_obj < node->lp_obj) best = node;
#endif
               }
            }
            break;
         default:
            xassert(T != T);
      }
      xassert(best != NULL);
      return best->p;
}

/* eof */
