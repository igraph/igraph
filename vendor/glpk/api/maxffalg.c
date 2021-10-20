/* maxffalg.c (find maximal flow with Ford-Fulkerson algorithm) */

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
#include "ffalg.h"
#include "glpk.h"

int glp_maxflow_ffalg(glp_graph *G, int s, int t, int a_cap,
      double *sol, int a_x, int v_cut)
{     /* find maximal flow with Ford-Fulkerson algorithm */
      glp_vertex *v;
      glp_arc *a;
      int nv, na, i, k, flag, *tail, *head, *cap, *x, ret;
      char *cut;
      double temp;
      if (!(1 <= s && s <= G->nv))
         xerror("glp_maxflow_ffalg: s = %d; source node number out of r"
            "ange\n", s);
      if (!(1 <= t && t <= G->nv))
         xerror("glp_maxflow_ffalg: t = %d: sink node number out of ran"
            "ge\n", t);
      if (s == t)
         xerror("glp_maxflow_ffalg: s = t = %d; source and sink nodes m"
            "ust be distinct\n", s);
      if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
         xerror("glp_maxflow_ffalg: a_cap = %d; invalid offset\n",
            a_cap);
      if (v_cut >= 0 && v_cut > G->v_size - (int)sizeof(int))
         xerror("glp_maxflow_ffalg: v_cut = %d; invalid offset\n",
            v_cut);
      /* allocate working arrays */
      nv = G->nv;
      na = G->na;
      tail = xcalloc(1+na, sizeof(int));
      head = xcalloc(1+na, sizeof(int));
      cap = xcalloc(1+na, sizeof(int));
      x = xcalloc(1+na, sizeof(int));
      if (v_cut < 0)
         cut = NULL;
      else
         cut = xcalloc(1+nv, sizeof(char));
      /* copy the flow network */
      k = 0;
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  k++;
            tail[k] = a->tail->i;
            head[k] = a->head->i;
            if (tail[k] == head[k])
            {  ret = GLP_EDATA;
               goto done;
            }
            if (a_cap >= 0)
               memcpy(&temp, (char *)a->data + a_cap, sizeof(double));
            else
               temp = 1.0;
            if (!(0.0 <= temp && temp <= (double)INT_MAX &&
                  temp == floor(temp)))
            {  ret = GLP_EDATA;
               goto done;
            }
            cap[k] = (int)temp;
         }
      }
      xassert(k == na);
      /* find maximal flow in the flow network */
      ffalg(nv, na, tail, head, s, t, cap, x, cut);
      ret = 0;
      /* store solution components */
      /* (objective function = total flow through the network) */
      if (sol != NULL)
      {  temp = 0.0;
         for (k = 1; k <= na; k++)
         {  if (tail[k] == s)
               temp += (double)x[k];
            else if (head[k] == s)
               temp -= (double)x[k];
         }
         *sol = temp;
      }
      /* (arc flows) */
      if (a_x >= 0)
      {  k = 0;
         for (i = 1; i <= G->nv; i++)
         {  v = G->v[i];
            for (a = v->out; a != NULL; a = a->t_next)
            {  temp = (double)x[++k];
               memcpy((char *)a->data + a_x, &temp, sizeof(double));
            }
         }
      }
      /* (node flags) */
      if (v_cut >= 0)
      {  for (i = 1; i <= G->nv; i++)
         {  v = G->v[i];
            flag = cut[i];
            memcpy((char *)v->data + v_cut, &flag, sizeof(int));
         }
      }
done: /* free working arrays */
      xfree(tail);
      xfree(head);
      xfree(cap);
      xfree(x);
      if (cut != NULL) xfree(cut);
      return ret;
}

/* eof */
