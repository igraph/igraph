/* asnokalg.c (solve assignment problem with out-of-kilter alg.) */

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
#include "okalg.h"

int glp_asnprob_okalg(int form, glp_graph *G, int v_set, int a_cost,
      double *sol, int a_x)
{     /* solve assignment problem with out-of-kilter algorithm */
      glp_vertex *v;
      glp_arc *a;
      int nv, na, i, k, *tail, *head, *low, *cap, *cost, *x, *pi, ret;
      double temp;
      if (!(form == GLP_ASN_MIN || form == GLP_ASN_MAX ||
            form == GLP_ASN_MMP))
         xerror("glp_asnprob_okalg: form = %d; invalid parameter\n",
            form);
      if (v_set >= 0 && v_set > G->v_size - (int)sizeof(int))
         xerror("glp_asnprob_okalg: v_set = %d; invalid offset\n",
            v_set);
      if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
         xerror("glp_asnprob_okalg: a_cost = %d; invalid offset\n",
            a_cost);
      if (a_x >= 0 && a_x > G->a_size - (int)sizeof(int))
         xerror("glp_asnprob_okalg: a_x = %d; invalid offset\n", a_x);
      if (glp_check_asnprob(G, v_set))
         return GLP_EDATA;
      /* nv is the total number of nodes in the resulting network */
      nv = G->nv + 1;
      /* na is the total number of arcs in the resulting network */
      na = G->na + G->nv;
      /* allocate working arrays */
      tail = xcalloc(1+na, sizeof(int));
      head = xcalloc(1+na, sizeof(int));
      low = xcalloc(1+na, sizeof(int));
      cap = xcalloc(1+na, sizeof(int));
      cost = xcalloc(1+na, sizeof(int));
      x = xcalloc(1+na, sizeof(int));
      pi = xcalloc(1+nv, sizeof(int));
      /* construct the resulting network */
      k = 0;
      /* (original arcs) */
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  k++;
            tail[k] = a->tail->i;
            head[k] = a->head->i;
            low[k] = 0;
            cap[k] = 1;
            if (a_cost >= 0)
               memcpy(&temp, (char *)a->data + a_cost, sizeof(double));
            else
               temp = 1.0;
            if (!(fabs(temp) <= (double)INT_MAX && temp == floor(temp)))
            {  ret = GLP_EDATA;
               goto done;
            }
            cost[k] = (int)temp;
            if (form != GLP_ASN_MIN) cost[k] = - cost[k];
         }
      }
      /* (artificial arcs) */
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         k++;
         if (v->out == NULL)
            tail[k] = i, head[k] = nv;
         else if (v->in == NULL)
            tail[k] = nv, head[k] = i;
         else
            xassert(v != v);
         low[k] = (form == GLP_ASN_MMP ? 0 : 1);
         cap[k] = 1;
         cost[k] = 0;
      }
      xassert(k == na);
      /* find minimal-cost circulation in the resulting network */
      ret = okalg(nv, na, tail, head, low, cap, cost, x, pi);
      switch (ret)
      {  case 0:
            /* optimal circulation found */
            ret = 0;
            break;
         case 1:
            /* no feasible circulation exists */
            ret = GLP_ENOPFS;
            break;
         case 2:
            /* integer overflow occured */
            ret = GLP_ERANGE;
            goto done;
         case 3:
            /* optimality test failed (logic error) */
            ret = GLP_EFAIL;
            goto done;
         default:
            xassert(ret != ret);
      }
      /* store solution components */
      /* (objective function = the total cost) */
      if (sol != NULL)
      {  temp = 0.0;
         for (k = 1; k <= na; k++)
            temp += (double)cost[k] * (double)x[k];
         if (form != GLP_ASN_MIN) temp = - temp;
         *sol = temp;
      }
      /* (arc flows) */
      if (a_x >= 0)
      {  k = 0;
         for (i = 1; i <= G->nv; i++)
         {  v = G->v[i];
            for (a = v->out; a != NULL; a = a->t_next)
            {  k++;
               if (ret == 0)
                  xassert(x[k] == 0 || x[k] == 1);
               memcpy((char *)a->data + a_x, &x[k], sizeof(int));
            }
         }
      }
done: /* free working arrays */
      xfree(tail);
      xfree(head);
      xfree(low);
      xfree(cap);
      xfree(cost);
      xfree(x);
      xfree(pi);
      return ret;
}

/* eof */
