/* cpp.c (solve critical path problem) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2010-2016 Free Software Foundation, Inc.
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
