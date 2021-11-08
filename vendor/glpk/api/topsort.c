/* topsort.c (topological sorting of acyclic digraph) */

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
*  glp_top_sort - topological sorting of acyclic digraph
*
*  SYNOPSIS
*
*  int glp_top_sort(glp_graph *G, int v_num);
*
*  DESCRIPTION
*
*  The routine glp_top_sort performs topological sorting of vertices of
*  the specified acyclic digraph.
*
*  The parameter v_num specifies an offset of the field of type int in
*  the vertex data block, to which the routine stores the vertex number
*  assigned. If v_num < 0, vertex numbers are not stored.
*
*  The vertices are numbered from 1 to n, where n is the total number
*  of vertices in the graph. The vertex numbering has the property that
*  for every arc (i->j) in the graph the condition num(i) < num(j)
*  holds. Special case num(i) = 0 means that vertex i is not assigned a
*  number, because the graph is *not* acyclic.
*
*  RETURNS
*
*  If the graph is acyclic and therefore all the vertices have been
*  assigned numbers, the routine glp_top_sort returns zero. Otherwise,
*  if the graph is not acyclic, the routine returns the number of
*  vertices which have not been numbered, i.e. for which num(i) = 0. */

static int top_sort(glp_graph *G, int num[])
{     glp_arc *a;
      int i, j, cnt, top, *stack, *indeg;
      /* allocate working arrays */
      indeg = xcalloc(1+G->nv, sizeof(int));
      stack = xcalloc(1+G->nv, sizeof(int));
      /* determine initial indegree of each vertex; push into the stack
         the vertices having zero indegree */
      top = 0;
      for (i = 1; i <= G->nv; i++)
      {  num[i] = indeg[i] = 0;
         for (a = G->v[i]->in; a != NULL; a = a->h_next)
            indeg[i]++;
         if (indeg[i] == 0)
            stack[++top] = i;
      }
      /* assign numbers to vertices in the sorted order */
      cnt = 0;
      while (top > 0)
      {  /* pull vertex i from the stack */
         i = stack[top--];
         /* it has zero indegree in the current graph */
         xassert(indeg[i] == 0);
         /* so assign it a next number */
         xassert(num[i] == 0);
         num[i] = ++cnt;
         /* remove vertex i from the current graph, update indegree of
            its adjacent vertices, and push into the stack new vertices
            whose indegree becomes zero */
         for (a = G->v[i]->out; a != NULL; a = a->t_next)
         {  j = a->head->i;
            /* there exists arc (i->j) in the graph */
            xassert(indeg[j] > 0);
            indeg[j]--;
            if (indeg[j] == 0)
               stack[++top] = j;
         }
      }
      /* free working arrays */
      xfree(indeg);
      xfree(stack);
      return G->nv - cnt;
}

int glp_top_sort(glp_graph *G, int v_num)
{     glp_vertex *v;
      int i, cnt, *num;
      if (v_num >= 0 && v_num > G->v_size - (int)sizeof(int))
         xerror("glp_top_sort: v_num = %d; invalid offset\n", v_num);
      if (G->nv == 0)
      {  cnt = 0;
         goto done;
      }
      num = xcalloc(1+G->nv, sizeof(int));
      cnt = top_sort(G, num);
      if (v_num >= 0)
      {  for (i = 1; i <= G->nv; i++)
         {  v = G->v[i];
            memcpy((char *)v->data + v_num, &num[i], sizeof(int));
         }
      }
      xfree(num);
done: return cnt;
}

/* eof */
