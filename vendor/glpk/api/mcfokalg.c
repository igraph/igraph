/* mcfokalg.c (find minimum-cost flow with out-of-kilter algorithm) */

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

int glp_mincost_okalg(glp_graph *G, int v_rhs, int a_low, int a_cap,
      int a_cost, double *sol, int a_x, int v_pi)
{     /* find minimum-cost flow with out-of-kilter algorithm */
      glp_vertex *v;
      glp_arc *a;
      int nv, na, i, k, s, t, *tail, *head, *low, *cap, *cost, *x, *pi,
         ret;
      double sum, temp;
      if (v_rhs >= 0 && v_rhs > G->v_size - (int)sizeof(double))
         xerror("glp_mincost_okalg: v_rhs = %d; invalid offset\n",
            v_rhs);
      if (a_low >= 0 && a_low > G->a_size - (int)sizeof(double))
         xerror("glp_mincost_okalg: a_low = %d; invalid offset\n",
            a_low);
      if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
         xerror("glp_mincost_okalg: a_cap = %d; invalid offset\n",
            a_cap);
      if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
         xerror("glp_mincost_okalg: a_cost = %d; invalid offset\n",
            a_cost);
      if (a_x >= 0 && a_x > G->a_size - (int)sizeof(double))
         xerror("glp_mincost_okalg: a_x = %d; invalid offset\n", a_x);
      if (v_pi >= 0 && v_pi > G->v_size - (int)sizeof(double))
         xerror("glp_mincost_okalg: v_pi = %d; invalid offset\n", v_pi);
      /* s is artificial source node */
      s = G->nv + 1;
      /* t is artificial sink node */
      t = s + 1;
      /* nv is the total number of nodes in the resulting network */
      nv = t;
      /* na is the total number of arcs in the resulting network */
      na = G->na + 1;
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         if (v_rhs >= 0)
            memcpy(&temp, (char *)v->data + v_rhs, sizeof(double));
         else
            temp = 0.0;
         if (temp != 0.0) na++;
      }
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
            if (tail[k] == head[k])
            {  ret = GLP_EDATA;
               goto done;
            }
            if (a_low >= 0)
               memcpy(&temp, (char *)a->data + a_low, sizeof(double));
            else
               temp = 0.0;
            if (!(0.0 <= temp && temp <= (double)INT_MAX &&
                  temp == floor(temp)))
            {  ret = GLP_EDATA;
               goto done;
            }
            low[k] = (int)temp;
            if (a_cap >= 0)
               memcpy(&temp, (char *)a->data + a_cap, sizeof(double));
            else
               temp = 1.0;
            if (!((double)low[k] <= temp && temp <= (double)INT_MAX &&
                  temp == floor(temp)))
            {  ret = GLP_EDATA;
               goto done;
            }
            cap[k] = (int)temp;
            if (a_cost >= 0)
               memcpy(&temp, (char *)a->data + a_cost, sizeof(double));
            else
               temp = 0.0;
            if (!(fabs(temp) <= (double)INT_MAX && temp == floor(temp)))
            {  ret = GLP_EDATA;
               goto done;
            }
            cost[k] = (int)temp;
         }
      }
      /* (artificial arcs) */
      sum = 0.0;
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         if (v_rhs >= 0)
            memcpy(&temp, (char *)v->data + v_rhs, sizeof(double));
         else
            temp = 0.0;
         if (!(fabs(temp) <= (double)INT_MAX && temp == floor(temp)))
         {  ret = GLP_EDATA;
            goto done;
         }
         if (temp > 0.0)
         {  /* artificial arc from s to original source i */
            k++;
            tail[k] = s;
            head[k] = i;
            low[k] = cap[k] = (int)(+temp); /* supply */
            cost[k] = 0;
            sum += (double)temp;
         }
         else if (temp < 0.0)
         {  /* artificial arc from original sink i to t */
            k++;
            tail[k] = i;
            head[k] = t;
            low[k] = cap[k] = (int)(-temp); /* demand */
            cost[k] = 0;
         }
      }
      /* (feedback arc from t to s) */
      k++;
      xassert(k == na);
      tail[k] = t;
      head[k] = s;
      if (sum > (double)INT_MAX)
      {  ret = GLP_EDATA;
         goto done;
      }
      low[k] = cap[k] = (int)sum; /* total supply/demand */
      cost[k] = 0;
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
      /* (node potentials = Lagrange multipliers) */
      if (v_pi >= 0)
      {  for (i = 1; i <= G->nv; i++)
         {  v = G->v[i];
            temp = - (double)pi[i];
            memcpy((char *)v->data + v_pi, &temp, sizeof(double));
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
