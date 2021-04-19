/* glpapi17.c (flow network problems) */

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

#include "glpapi.h"
#include "glpnet.h"

/***********************************************************************
*  NAME
*
*  glp_mincost_lp - convert minimum cost flow problem to LP
*
*  SYNOPSIS
*
*  void glp_mincost_lp(glp_prob *lp, glp_graph *G, int names,
*     int v_rhs, int a_low, int a_cap, int a_cost);
*
*  DESCRIPTION
*
*  The routine glp_mincost_lp builds an LP problem, which corresponds
*  to the minimum cost flow problem on the specified network G. */

void glp_mincost_lp(glp_prob *lp, glp_graph *G, int names, int v_rhs,
      int a_low, int a_cap, int a_cost)
{     glp_vertex *v;
      glp_arc *a;
      int i, j, type, ind[1+2];
      double rhs, low, cap, cost, val[1+2];
      if (!(names == GLP_ON || names == GLP_OFF))
         xerror("glp_mincost_lp: names = %d; invalid parameter\n",
            names);
      if (v_rhs >= 0 && v_rhs > G->v_size - (int)sizeof(double))
         xerror("glp_mincost_lp: v_rhs = %d; invalid offset\n", v_rhs);
      if (a_low >= 0 && a_low > G->a_size - (int)sizeof(double))
         xerror("glp_mincost_lp: a_low = %d; invalid offset\n", a_low);
      if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
         xerror("glp_mincost_lp: a_cap = %d; invalid offset\n", a_cap);
      if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
         xerror("glp_mincost_lp: a_cost = %d; invalid offset\n", a_cost)
            ;
      glp_erase_prob(lp);
      if (names) glp_set_prob_name(lp, G->name);
      if (G->nv > 0) glp_add_rows(lp, G->nv);
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         if (names) glp_set_row_name(lp, i, v->name);
         if (v_rhs >= 0)
            memcpy(&rhs, (char *)v->data + v_rhs, sizeof(double));
         else
            rhs = 0.0;
         glp_set_row_bnds(lp, i, GLP_FX, rhs, rhs);
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
            if (a_low >= 0)
               memcpy(&low, (char *)a->data + a_low, sizeof(double));
            else
               low = 0.0;
            if (a_cap >= 0)
               memcpy(&cap, (char *)a->data + a_cap, sizeof(double));
            else
               cap = 1.0;
            if (cap == DBL_MAX)
               type = GLP_LO;
            else if (low != cap)
               type = GLP_DB;
            else
               type = GLP_FX;
            glp_set_col_bnds(lp, j, type, low, cap);
            if (a_cost >= 0)
               memcpy(&cost, (char *)a->data + a_cost, sizeof(double));
            else
               cost = 0.0;
            glp_set_obj_coef(lp, j, cost);
         }
      }
      xassert(j == G->na);
      return;
}

/**********************************************************************/

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

/***********************************************************************
*  NAME
*
*  glp_check_asnprob - check correctness of assignment problem data
*
*  SYNOPSIS
*
*  int glp_check_asnprob(glp_graph *G, int v_set);
*
*  RETURNS
*
*  If the specified assignment problem data are correct, the routine
*  glp_check_asnprob returns zero, otherwise, non-zero. */

int glp_check_asnprob(glp_graph *G, int v_set)
{     glp_vertex *v;
      int i, k, ret = 0;
      if (v_set >= 0 && v_set > G->v_size - (int)sizeof(int))
         xerror("glp_check_asnprob: v_set = %d; invalid offset\n",
            v_set);
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         if (v_set >= 0)
         {  memcpy(&k, (char *)v->data + v_set, sizeof(int));
            if (k == 0)
            {  if (v->in != NULL)
               {  ret = 1;
                  break;
               }
            }
            else if (k == 1)
            {  if (v->out != NULL)
               {  ret = 2;
                  break;
               }
            }
            else
            {  ret = 3;
               break;
            }
         }
         else
         {  if (v->in != NULL && v->out != NULL)
            {  ret = 4;
               break;
            }
         }
      }
      return ret;
}

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

/**********************************************************************/

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

/***********************************************************************
*  NAME
*
*  glp_asnprob_hall - find bipartite matching of maximum cardinality
*
*  SYNOPSIS
*
*  int glp_asnprob_hall(glp_graph *G, int v_set, int a_x);
*
*  DESCRIPTION
*
*  The routine glp_asnprob_hall finds a matching of maximal cardinality
*  in the specified bipartite graph G. It uses a version of the Fortran
*  routine MC21A developed by I.S.Duff [1], which implements Hall's
*  algorithm [2].
*
*  RETURNS
*
*  The routine glp_asnprob_hall returns the cardinality of the matching
*  found. However, if the specified graph is incorrect (as detected by
*  the routine glp_check_asnprob), the routine returns negative value.
*
*  REFERENCES
*
*  1. I.S.Duff, Algorithm 575: Permutations for zero-free diagonal, ACM
*     Trans. on Math. Softw. 7 (1981), 387-390.
*
*  2. M.Hall, "An Algorithm for distinct representatives," Amer. Math.
*     Monthly 63 (1956), 716-717. */

int glp_asnprob_hall(glp_graph *G, int v_set, int a_x)
{     glp_vertex *v;
      glp_arc *a;
      int card, i, k, loc, n, n1, n2, xij;
      int *num, *icn, *ip, *lenr, *iperm, *pr, *arp, *cv, *out;
      if (v_set >= 0 && v_set > G->v_size - (int)sizeof(int))
         xerror("glp_asnprob_hall: v_set = %d; invalid offset\n",
            v_set);
      if (a_x >= 0 && a_x > G->a_size - (int)sizeof(int))
         xerror("glp_asnprob_hall: a_x = %d; invalid offset\n", a_x);
      if (glp_check_asnprob(G, v_set))
         return -1;
      /* determine the number of vertices in sets R and S and renumber
         vertices in S which correspond to columns of the matrix; skip
         all isolated vertices */
      num = xcalloc(1+G->nv, sizeof(int));
      n1 = n2 = 0;
      for (i = 1; i <= G->nv; i++)
      {  v = G->v[i];
         if (v->in == NULL && v->out != NULL)
            n1++, num[i] = 0; /* vertex in R */
         else if (v->in != NULL && v->out == NULL)
            n2++, num[i] = n2; /* vertex in S */
         else
         {  xassert(v->in == NULL && v->out == NULL);
            num[i] = -1; /* isolated vertex */
         }
      }
      /* the matrix must be square, thus, if it has more columns than
         rows, extra rows will be just empty, and vice versa */
      n = (n1 >= n2 ? n1 : n2);
      /* allocate working arrays */
      icn = xcalloc(1+G->na, sizeof(int));
      ip = xcalloc(1+n, sizeof(int));
      lenr = xcalloc(1+n, sizeof(int));
      iperm = xcalloc(1+n, sizeof(int));
      pr = xcalloc(1+n, sizeof(int));
      arp = xcalloc(1+n, sizeof(int));
      cv = xcalloc(1+n, sizeof(int));
      out = xcalloc(1+n, sizeof(int));
      /* build the adjacency matrix of the bipartite graph in row-wise
         format (rows are vertices in R, columns are vertices in S) */
      k = 0, loc = 1;
      for (i = 1; i <= G->nv; i++)
      {  if (num[i] != 0) continue;
         /* vertex i in R */
         ip[++k] = loc;
         v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  xassert(num[a->head->i] != 0);
            icn[loc++] = num[a->head->i];
         }
         lenr[k] = loc - ip[k];
      }
      xassert(loc-1 == G->na);
      /* make all extra rows empty (all extra columns are empty due to
         the row-wise format used) */
      for (k++; k <= n; k++)
         ip[k] = loc, lenr[k] = 0;
      /* find a row permutation that maximizes the number of non-zeros
         on the main diagonal */
      card = mc21a(n, icn, ip, lenr, iperm, pr, arp, cv, out);
#if 1 /* 18/II-2010 */
      /* FIXED: if card = n, arp remains clobbered on exit */
      for (i = 1; i <= n; i++)
         arp[i] = 0;
      for (i = 1; i <= card; i++)
      {  k = iperm[i];
         xassert(1 <= k && k <= n);
         xassert(arp[k] == 0);
         arp[k] = i;
      }
#endif
      /* store solution, if necessary */
      if (a_x < 0) goto skip;
      k = 0;
      for (i = 1; i <= G->nv; i++)
      {  if (num[i] != 0) continue;
         /* vertex i in R */
         k++;
         v = G->v[i];
         for (a = v->out; a != NULL; a = a->t_next)
         {  /* arp[k] is the number of matched column or zero */
            if (arp[k] == num[a->head->i])
            {  xassert(arp[k] != 0);
               xij = 1;
            }
            else
               xij = 0;
            memcpy((char *)a->data + a_x, &xij, sizeof(int));
         }
      }
skip: /* free working arrays */
      xfree(num);
      xfree(icn);
      xfree(ip);
      xfree(lenr);
      xfree(iperm);
      xfree(pr);
      xfree(arp);
      xfree(cv);
      xfree(out);
      return card;
}

/***********************************************************************
*  NAME
*
*  glp_cpp - solve critical path problem
*
*  SYNOPSIS
*
*  double glp_cpp(glp_graph *G, int v_t, int v_es, int v_ls);
*
*  DESCRIPTION
*
*  The routine glp_cpp solves the critical path problem represented in
*  the form of the project network.
*
*  The parameter G is a pointer to the graph object, which specifies
*  the project network. This graph must be acyclic. Multiple arcs are
*  allowed being considered as single arcs.
*
*  The parameter v_t specifies an offset of the field of type double
*  in the vertex data block, which contains time t[i] >= 0 needed to
*  perform corresponding job j. If v_t < 0, it is assumed that t[i] = 1
*  for all jobs.
*
*  The parameter v_es specifies an offset of the field of type double
*  in the vertex data block, to which the routine stores earliest start
*  time for corresponding job. If v_es < 0, this time is not stored.
*
*  The parameter v_ls specifies an offset of the field of type double
*  in the vertex data block, to which the routine stores latest start
*  time for corresponding job. If v_ls < 0, this time is not stored.
*
*  RETURNS
*
*  The routine glp_cpp returns the minimal project duration, that is,
*  minimal time needed to perform all jobs in the project. */

static void sorting(glp_graph *G, int list[]);

double glp_cpp(glp_graph *G, int v_t, int v_es, int v_ls)
{     glp_vertex *v;
      glp_arc *a;
      int i, j, k, nv, *list;
      double temp, total, *t, *es, *ls;
      if (v_t >= 0 && v_t > G->v_size - (int)sizeof(double))
         xerror("glp_cpp: v_t = %d; invalid offset\n", v_t);
      if (v_es >= 0 && v_es > G->v_size - (int)sizeof(double))
         xerror("glp_cpp: v_es = %d; invalid offset\n", v_es);
      if (v_ls >= 0 && v_ls > G->v_size - (int)sizeof(double))
         xerror("glp_cpp: v_ls = %d; invalid offset\n", v_ls);
      nv = G->nv;
      if (nv == 0)
      {  total = 0.0;
         goto done;
      }
      /* allocate working arrays */
      t = xcalloc(1+nv, sizeof(double));
      es = xcalloc(1+nv, sizeof(double));
      ls = xcalloc(1+nv, sizeof(double));
      list = xcalloc(1+nv, sizeof(int));
      /* retrieve job times */
      for (i = 1; i <= nv; i++)
      {  v = G->v[i];
         if (v_t >= 0)
         {  memcpy(&t[i], (char *)v->data + v_t, sizeof(double));
            if (t[i] < 0.0)
               xerror("glp_cpp: t[%d] = %g; invalid time\n", i, t[i]);
         }
         else
            t[i] = 1.0;
      }
      /* perform topological sorting to determine the list of nodes
         (jobs) such that if list[k] = i and list[kk] = j and there
         exists arc (i->j), then k < kk */
      sorting(G, list);
      /* FORWARD PASS */
      /* determine earliest start times */
      for (k = 1; k <= nv; k++)
      {  j = list[k];
         es[j] = 0.0;
         for (a = G->v[j]->in; a != NULL; a = a->h_next)
         {  i = a->tail->i;
            /* there exists arc (i->j) in the project network */
            temp = es[i] + t[i];
            if (es[j] < temp) es[j] = temp;
         }
      }
      /* determine the minimal project duration */
      total = 0.0;
      for (i = 1; i <= nv; i++)
      {  temp = es[i] + t[i];
         if (total < temp) total = temp;
      }
      /* BACKWARD PASS */
      /* determine latest start times */
      for (k = nv; k >= 1; k--)
      {  i = list[k];
         ls[i] = total - t[i];
         for (a = G->v[i]->out; a != NULL; a = a->t_next)
         {  j = a->head->i;
            /* there exists arc (i->j) in the project network */
            temp = ls[j] - t[i];
            if (ls[i] > temp) ls[i] = temp;
         }
         /* avoid possible round-off errors */
         if (ls[i] < es[i]) ls[i] = es[i];
      }
      /* store results, if necessary */
      if (v_es >= 0)
      {  for (i = 1; i <= nv; i++)
         {  v = G->v[i];
            memcpy((char *)v->data + v_es, &es[i], sizeof(double));
         }
      }
      if (v_ls >= 0)
      {  for (i = 1; i <= nv; i++)
         {  v = G->v[i];
            memcpy((char *)v->data + v_ls, &ls[i], sizeof(double));
         }
      }
      /* free working arrays */
      xfree(t);
      xfree(es);
      xfree(ls);
      xfree(list);
done: return total;
}

static void sorting(glp_graph *G, int list[])
{     /* perform topological sorting to determine the list of nodes
         (jobs) such that if list[k] = i and list[kk] = j and there
         exists arc (i->j), then k < kk */
      int i, k, nv, v_size, *num;
      void **save;
      nv = G->nv;
      v_size = G->v_size;
      save = xcalloc(1+nv, sizeof(void *));
      num = xcalloc(1+nv, sizeof(int));
      G->v_size = sizeof(int);
      for (i = 1; i <= nv; i++)
      {  save[i] = G->v[i]->data;
         G->v[i]->data = &num[i];
         list[i] = 0;
      }
      if (glp_top_sort(G, 0) != 0)
         xerror("glp_cpp: project network is not acyclic\n");
      G->v_size = v_size;
      for (i = 1; i <= nv; i++)
      {  G->v[i]->data = save[i];
         k = num[i];
         xassert(1 <= k && k <= nv);
         xassert(list[k] == 0);
         list[k] = i;
      }
      xfree(save);
      xfree(num);
      return;
}

/* eof */
