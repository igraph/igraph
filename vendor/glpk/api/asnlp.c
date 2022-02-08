/* asnlp.c (convert assignment problem to LP) */

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
*  glp_asnprob_lp - convert assignment problem to LP
*
*  SYNOPSIS
*
*  int glp_asnprob_lp(glp_prob *P, int form, glp_graph *G, int names,
*     int v_set, int a_cost);
*
*  DESCRIPTION
*
*  The routine glp_asnprob_lp builds an LP problem, which corresponds
*  to the assignment problem on the specified graph G.
*
*  RETURNS
*
*  If the LP problem has been successfully built, the routine returns
*  zero, otherwise, non-zero. */

int glp_asnprob_lp(glp_prob *P, int form, glp_graph *G, int names,
      int v_set, int a_cost)
{     glp_vertex *v;
      glp_arc *a;
      int i, j, ret, ind[1+2];
      double cost, val[1+2];
      if (!(form == GLP_ASN_MIN || form == GLP_ASN_MAX ||
            form == GLP_ASN_MMP))
         xerror("glp_asnprob_lp: form = %d; invalid parameter\n",
            form);
      if (!(names == GLP_ON || names == GLP_OFF))
         xerror("glp_asnprob_lp: names = %d; invalid parameter\n",
            names);
      if (v_set >= 0 && v_set > G->v_size - (int)sizeof(int))
         xerror("glp_asnprob_lp: v_set = %d; invalid offset\n",
            v_set);
      if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
         xerror("glp_asnprob_lp: a_cost = %d; invalid offset\n",
            a_cost);
      ret = glp_check_asnprob(G, v_set);
      if (ret != 0) goto done;
      glp_erase_prob(P);
      if (names) glp_set_prob_name(P, G->name);
      glp_set_obj_dir(P, form == GLP_ASN_MIN ? GLP_MIN : GLP_MAX);
      if (G->nv > 0) glp_add_rows(P, G->nv);
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         if (names) glp_set_row_name(P, i, v->name);
         glp_set_row_bnds(P, i, form == GLP_ASN_MMP ? GLP_UP : GLP_FX,
            1.0, 1.0);
      }
      if (G->na > 0) glp_add_cols(P, G->na);
      for (i = 1, j = 0; i <= G->nv; i++)
      {  v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  j++;
            if (names)
            {  char name[50+1];
               sprintf(name, "x[%d,%d]", a->tail->i, a->head->i);
               xassert(strlen(name) < sizeof(name));
               glp_set_col_name(P, j, name);
            }
            ind[1] = a->tail->i, val[1] = +1.0;
            ind[2] = a->head->i, val[2] = +1.0;
            glp_set_mat_col(P, j, 2, ind, val);
            glp_set_col_bnds(P, j, GLP_DB, 0.0, 1.0);
            if (a_cost >= 0)
               memcpy(&cost, (char *)a->data + a_cost, sizeof(double));
            else
               cost = 1.0;
            glp_set_obj_coef(P, j, cost);
         }
      }
      xassert(j == G->na);
done: return ret;
}

/* eof */
