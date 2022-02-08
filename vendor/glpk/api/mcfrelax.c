/* mcfrelax.c (find minimum-cost flow with RELAX-IV) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2013-2016 Free Software Foundation, Inc.
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
#include "relax4.h"

static int overflow(int u, int v)
{     /* check for integer overflow on computing u + v */
      if (u > 0 && v > 0 && u + v < 0) return 1;
      if (u < 0 && v < 0 && u + v > 0) return 1;
      return 0;
}

int glp_mincost_relax4(glp_graph *G, int v_rhs, int a_low, int a_cap,
      int a_cost, int crash, double *sol, int a_x, int a_rc)
{     /* find minimum-cost flow with Bertsekas-Tseng relaxation method
         (RELAX-IV) */
      glp_vertex *v;
      glp_arc *a;
      struct relax4_csa csa;
      int i, k, large, n, na, ret;
      double cap, cost, low, rc, rhs, sum, x;
      if (v_rhs >= 0 && v_rhs > G->v_size - (int)sizeof(double))
         xerror("glp_mincost_relax4: v_rhs = %d; invalid offset\n",
            v_rhs);
      if (a_low >= 0 && a_low > G->a_size - (int)sizeof(double))
         xerror("glp_mincost_relax4: a_low = %d; invalid offset\n",
            a_low);
      if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
         xerror("glp_mincost_relax4: a_cap = %d; invalid offset\n",
            a_cap);
      if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
         xerror("glp_mincost_relax4: a_cost = %d; invalid offset\n",
            a_cost);
      if (a_x >= 0 && a_x > G->a_size - (int)sizeof(double))
         xerror("glp_mincost_relax4: a_x = %d; invalid offset\n",
            a_x);
      if (a_rc >= 0 && a_rc > G->a_size - (int)sizeof(double))
         xerror("glp_mincost_relax4: a_rc = %d; invalid offset\n",
            a_rc);
      csa.n = n = G->nv; /* number of nodes */
      csa.na = na = G->na; /* number of arcs */
      csa.large = large = INT_MAX / 4;
      csa.repeat = 0;
      csa.crash = crash;
      /* allocate working arrays */
      csa.startn = xcalloc(1+na, sizeof(int));
      csa.endn = xcalloc(1+na, sizeof(int));
      csa.fou = xcalloc(1+n, sizeof(int));
      csa.nxtou = xcalloc(1+na, sizeof(int));
      csa.fin = xcalloc(1+n, sizeof(int));
      csa.nxtin = xcalloc(1+na, sizeof(int));
      csa.rc = xcalloc(1+na, sizeof(int));
      csa.u = xcalloc(1+na, sizeof(int));
      csa.dfct = xcalloc(1+n, sizeof(int));
      csa.x = xcalloc(1+na, sizeof(int));
      csa.label = xcalloc(1+n, sizeof(int));
      csa.prdcsr = xcalloc(1+n, sizeof(int));
      csa.save = xcalloc(1+na, sizeof(int));
      csa.tfstou = xcalloc(1+n, sizeof(int));
      csa.tnxtou = xcalloc(1+na, sizeof(int));
      csa.tfstin = xcalloc(1+n, sizeof(int));
      csa.tnxtin = xcalloc(1+na, sizeof(int));
      csa.nxtqueue = xcalloc(1+n, sizeof(int));
      csa.scan = xcalloc(1+n, sizeof(char));
      csa.mark = xcalloc(1+n, sizeof(char));
      if (crash)
      {  csa.extend_arc = xcalloc(1+n, sizeof(int));
         csa.sb_level = xcalloc(1+n, sizeof(int));
         csa.sb_arc = xcalloc(1+n, sizeof(int));
      }
      else
      {  csa.extend_arc = NULL;
         csa.sb_level = NULL;
         csa.sb_arc = NULL;
      }
      /* scan nodes */
      for (i = 1; i <= n; i++)
      {  v = G->v[i];
         /* get supply at i-th node */
         if (v_rhs >= 0)
            memcpy(&rhs, (char *)v->data + v_rhs, sizeof(double));
         else
            rhs = 0.0;
         if (!(fabs(rhs) <= (double)large && rhs == floor(rhs)))
         {  ret = GLP_EDATA;
            goto done;
         }
         /* set demand at i-th node */
         csa.dfct[i] = -(int)rhs;
      }
      /* scan arcs */
      k = 0;
      for (i = 1; i <= n; i++)
      {  v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  k++;
            /* set endpoints of k-th arc */
            if (a->tail->i == a->head->i)
            {  /* self-loops not allowed */
               ret = GLP_EDATA;
               goto done;
            }
            csa.startn[k] = a->tail->i;
            csa.endn[k] = a->head->i;
            /* set per-unit cost for k-th arc flow */
            if (a_cost >= 0)
               memcpy(&cost, (char *)a->data + a_cost, sizeof(double));
            else
               cost = 0.0;
            if (!(fabs(cost) <= (double)large && cost == floor(cost)))
            {  ret = GLP_EDATA;
               goto done;
            }
            csa.rc[k] = (int)cost;
            /* get lower bound for k-th arc flow */
            if (a_low >= 0)
               memcpy(&low, (char *)a->data + a_low, sizeof(double));
            else
               low = 0.0;
            if (!(0.0 <= low && low <= (double)large &&
                  low == floor(low)))
            {  ret = GLP_EDATA;
               goto done;
            }
            /* get upper bound for k-th arc flow */
            if (a_cap >= 0)
               memcpy(&cap, (char *)a->data + a_cap, sizeof(double));
            else
               cap = 1.0;
            if (!(low <= cap && cap <= (double)large &&
                  cap == floor(cap)))
            {  ret = GLP_EDATA;
               goto done;
            }
            /* substitute x = x' + low, where 0 <= x' <= cap - low */
            csa.u[k] = (int)(cap - low);
            /* correct demands at endpoints of k-th arc */
            if (overflow(csa.dfct[a->tail->i], +low))
            {  ret = GLP_ERANGE;
               goto done;
            }
#if 0 /* 29/IX-2017 */
            csa.dfct[a->tail->i] += low;
#else
            csa.dfct[a->tail->i] += (int)low;
#endif
            if (overflow(csa.dfct[a->head->i], -low))
            {  ret = GLP_ERANGE;
               goto done;
            }
#if 0 /* 29/IX-2017 */
            csa.dfct[a->head->i] -= low;
#else
            csa.dfct[a->head->i] -= (int)low;
#endif
         }
      }
      /* construct linked list for network topology */
      relax4_inidat(&csa);
      /* find minimum-cost flow */
      ret = relax4(&csa);
      if (ret != 0)
      {  /* problem is found to be infeasible */
         xassert(1 <= ret && ret <= 8);
         ret = GLP_ENOPFS;
         goto done;
      }
      /* store solution */
      sum = 0.0;
      k = 0;
      for (i = 1; i <= n; i++)
      {  v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  k++;
            /* get lower bound for k-th arc flow */
            if (a_low >= 0)
               memcpy(&low, (char *)a->data + a_low, sizeof(double));
            else
               low = 0.0;
            /* store original flow x = x' + low thru k-th arc */
            x = (double)csa.x[k] + low;
            if (a_x >= 0)
               memcpy((char *)a->data + a_x, &x, sizeof(double));
            /* store reduced cost for k-th arc flow */
            rc = (double)csa.rc[k];
            if (a_rc >= 0)
               memcpy((char *)a->data + a_rc, &rc, sizeof(double));
            /* get per-unit cost for k-th arc flow */
            if (a_cost >= 0)
               memcpy(&cost, (char *)a->data + a_cost, sizeof(double));
            else
               cost = 0.0;
            /* compute the total cost */
            sum += cost * x;
         }
      }
      /* store the total cost */
      if (sol != NULL)
         *sol = sum;
done: /* free working arrays */
      xfree(csa.startn);
      xfree(csa.endn);
      xfree(csa.fou);
      xfree(csa.nxtou);
      xfree(csa.fin);
      xfree(csa.nxtin);
      xfree(csa.rc);
      xfree(csa.u);
      xfree(csa.dfct);
      xfree(csa.x);
      xfree(csa.label);
      xfree(csa.prdcsr);
      xfree(csa.save);
      xfree(csa.tfstou);
      xfree(csa.tnxtou);
      xfree(csa.tfstin);
      xfree(csa.tnxtin);
      xfree(csa.nxtqueue);
      xfree(csa.scan);
      xfree(csa.mark);
      if (crash)
      {  xfree(csa.extend_arc);
         xfree(csa.sb_level);
         xfree(csa.sb_arc);
      }
      return ret;
}

/* eof */
