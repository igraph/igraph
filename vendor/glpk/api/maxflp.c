/* maxflp.c (convert maximum flow problem to LP) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2009-2016 Free Software Foundation, Inc.
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
#include "glpk.h"

/***********************************************************************
*  NAME
*
*  glp_maxflow_lp - convert maximum flow problem to LP
*
*  SYNOPSIS
*
*  void glp_maxflow_lp(glp_prob *lp, glp_graph *G, int names, int s,
*     int t, int a_cap);
*
*  DESCRIPTION
*
*  The routine glp_maxflow_lp builds an LP problem, which corresponds
*  to the maximum flow problem on the specified network G. */

void glp_maxflow_lp(glp_prob *lp, glp_graph *G, int names, int s,
      int t, int a_cap)
{     glp_vertex *v;
      glp_arc *a;
      int i, j, type, ind[1+2];
      double cap, val[1+2];
      if (!(names == GLP_ON || names == GLP_OFF))
         xerror("glp_maxflow_lp: names = %d; invalid parameter\n",
            names);
      if (!(1 <= s && s <= G->nv))
         xerror("glp_maxflow_lp: s = %d; source node number out of rang"
            "e\n", s);
      if (!(1 <= t && t <= G->nv))
         xerror("glp_maxflow_lp: t = %d: sink node number out of range "
            "\n", t);
      if (s == t)
         xerror("glp_maxflow_lp: s = t = %d; source and sink nodes must"
            " be distinct\n", s);
      if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
         xerror("glp_maxflow_lp: a_cap = %d; invalid offset\n", a_cap);
      glp_erase_prob(lp);
      if (names) glp_set_prob_name(lp, G->name);
      glp_set_obj_dir(lp, GLP_MAX);
      glp_add_rows(lp, G->nv);
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         if (names) glp_set_row_name(lp, i, v->name);
         if (i == s)
            type = GLP_LO;
         else if (i == t)
            type = GLP_UP;
         else
            type = GLP_FX;
         glp_set_row_bnds(lp, i, type, 0.0, 0.0);
      }
      if (G->na > 0) glp_add_cols(lp, G->na);
      for (i = 1, j = 0; i <= G->nv; i++)
      {  v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  j++;
            if (names)
            {  char name[50+1];
               sprintf(name, "x[%d,%d]", a->tail->i, a->head->i);
               xassert(strlen(name) < sizeof(name));
               glp_set_col_name(lp, j, name);
            }
            if (a->tail->i != a->head->i)
            {  ind[1] = a->tail->i, val[1] = +1.0;
               ind[2] = a->head->i, val[2] = -1.0;
               glp_set_mat_col(lp, j, 2, ind, val);
            }
            if (a_cap >= 0)
               memcpy(&cap, (char *)a->data + a_cap, sizeof(double));
            else
               cap = 1.0;
            if (cap == DBL_MAX)
               type = GLP_LO;
            else if (cap != 0.0)
               type = GLP_DB;
            else
               type = GLP_FX;
            glp_set_col_bnds(lp, j, type, 0.0, cap);
            if (a->tail->i == s)
               glp_set_obj_coef(lp, j, +1.0);
            else if (a->head->i == s)
               glp_set_obj_coef(lp, j, -1.0);
         }
      }
      xassert(j == G->na);
      return;
}

/* eof */
