/* keller.c (cover edges by cliques, Kellerman's heuristic) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2009-2013 Free Software Foundation, Inc.
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

#include "glpk.h"
#include "env.h"
#include "keller.h"

/***********************************************************************
*  NAME
*
*  kellerman - cover edges by cliques with Kellerman's heuristic
*
*  SYNOPSIS
*
*  #include "keller.h"
*  int kellerman(int n, int (*func)(void *info, int i, int ind[]),
*     void *info, glp_graph *H);
*
*  DESCRIPTION
*
*  The routine kellerman implements Kellerman's heuristic algorithm
*  to find a minimal set of cliques which cover all edges of specified
*  graph G = (V, E).
*
*  The parameter n specifies the number of vertices |V|, n >= 0.
*
*  Formal routine func specifies the set of edges E in the following
*  way. Running the routine kellerman calls the routine func and passes
*  to it parameter i, which is the number of some vertex, 1 <= i <= n.
*  In response the routine func should store numbers of all vertices
*  adjacent to vertex i to locations ind[1], ind[2], ..., ind[len] and
*  return the value of len, which is the number of adjacent vertices,
*  0 <= len <= n. Self-loops are allowed, but ignored. Multiple edges
*  are not allowed.
*
*  The parameter info is a transit pointer (magic cookie) passed to the
*  formal routine func as its first parameter.
*
*  The result provided by the routine kellerman is the bipartite graph
*  H = (V union C, F), which defines the covering found. (The program
*  object of type glp_graph specified by the parameter H should be
*  previously created with the routine glp_create_graph. On entry the
*  routine kellerman erases the content of this object with the routine
*  glp_erase_graph.) Vertices of first part V correspond to vertices of
*  the graph G and have the same ordinal numbers 1, 2, ..., n. Vertices
*  of second part C correspond to cliques and have ordinal numbers
*  n+1, n+2, ..., n+k, where k is the total number of cliques in the
*  edge covering found. Every edge f in F in the program object H is
*  represented as arc f = (i->j), where i in V and j in C, which means
*  that vertex i of the graph G is in clique C[j], 1 <= j <= k. (Thus,
*  if two vertices of the graph G are in the same clique, these vertices
*  are adjacent in G, and corresponding edge is covered by that clique.)
*
*  RETURNS
*
*  The routine Kellerman returns k, the total number of cliques in the
*  edge covering found.
*
*  REFERENCE
*
*  For more details see: glpk/doc/notes/keller.pdf (in Russian). */

struct set
{     /* set of vertices */
      int size;
      /* size (cardinality) of the set, 0 <= card <= n */
      int *list; /* int list[1+n]; */
      /* the set contains vertices list[1,...,size] */
      int *pos; /* int pos[1+n]; */
      /* pos[i] > 0 means that vertex i is in the set and
       * list[pos[i]] = i; pos[i] = 0 means that vertex i is not in
       * the set */
};

int kellerman(int n, int (*func)(void *info, int i, int ind[]),
      void *info, void /* glp_graph */ *H_)
{     glp_graph *H = H_;
      struct set W_, *W = &W_, V_, *V = &V_;
      glp_arc *a;
      int i, j, k, m, t, len, card, best;
      xassert(n >= 0);
      /* H := (V, 0; 0), where V is the set of vertices of graph G */
      glp_erase_graph(H, H->v_size, H->a_size);
      glp_add_vertices(H, n);
      /* W := 0 */
      W->size = 0;
      W->list = xcalloc(1+n, sizeof(int));
      W->pos = xcalloc(1+n, sizeof(int));
      memset(&W->pos[1], 0, sizeof(int) * n);
      /* V := 0 */
      V->size = 0;
      V->list = xcalloc(1+n, sizeof(int));
      V->pos = xcalloc(1+n, sizeof(int));
      memset(&V->pos[1], 0, sizeof(int) * n);
      /* main loop */
      for (i = 1; i <= n; i++)
      {  /* W must be empty */
         xassert(W->size == 0);
         /* W := { j : i > j and (i,j) in E } */
         len = func(info, i, W->list);
         xassert(0 <= len && len <= n);
         for (t = 1; t <= len; t++)
         {  j = W->list[t];
            xassert(1 <= j && j <= n);
            if (j >= i) continue;
            xassert(W->pos[j] == 0);
            W->list[++W->size] = j, W->pos[j] = W->size;
         }
         /* on i-th iteration we need to cover edges (i,j) for all
          * j in W */
         /* if W is empty, it is a special case */
         if (W->size == 0)
         {  /* set k := k + 1 and create new clique C[k] = { i } */
            k = glp_add_vertices(H, 1) - n;
            glp_add_arc(H, i, n + k);
            continue;
         }
         /* try to include vertex i into existing cliques */
         /* V must be empty */
         xassert(V->size == 0);
         /* k is the number of cliques found so far */
         k = H->nv - n;
         for (m = 1; m <= k; m++)
         {  /* do while V != W; since here V is within W, we can use
             * equivalent condition: do while |V| < |W| */
            if (V->size == W->size) break;
            /* check if C[m] is within W */
            for (a = H->v[n + m]->in; a != NULL; a = a->h_next)
            {  j = a->tail->i;
               if (W->pos[j] == 0) break;
            }
            if (a != NULL) continue;
            /* C[m] is within W, expand clique C[m] with vertex i */
            /* C[m] := C[m] union {i} */
            glp_add_arc(H, i, n + m);
            /* V is a set of vertices whose incident edges are already
             * covered by existing cliques */
            /* V := V union C[m] */
            for (a = H->v[n + m]->in; a != NULL; a = a->h_next)
            {  j = a->tail->i;
               if (V->pos[j] == 0)
                  V->list[++V->size] = j, V->pos[j] = V->size;
            }
         }
         /* remove from set W the vertices whose incident edges are
          * already covered by existing cliques */
         /* W := W \ V, V := 0 */
         for (t = 1; t <= V->size; t++)
         {  j = V->list[t], V->pos[j] = 0;
            if (W->pos[j] != 0)
            {  /* remove vertex j from W */
               if (W->pos[j] != W->size)
               {  int jj = W->list[W->size];
                  W->list[W->pos[j]] = jj;
                  W->pos[jj] = W->pos[j];
               }
               W->size--, W->pos[j] = 0;
            }
         }
         V->size = 0;
         /* now set W contains only vertices whose incident edges are
          * still not covered by existing cliques; create new cliques
          * to cover remaining edges until set W becomes empty */
         while (W->size > 0)
         {  /* find clique C[m], 1 <= m <= k, which shares maximal
             * number of vertices with W; to break ties choose clique
             * having smallest number m */
            m = 0, best = -1;
            k = H->nv - n;
            for (t = 1; t <= k; t++)
            {  /* compute cardinality of intersection of W and C[t] */
               card = 0;
               for (a = H->v[n + t]->in; a != NULL; a = a->h_next)
               {  j = a->tail->i;
                  if (W->pos[j] != 0) card++;
               }
               if (best < card)
                  m = t, best = card;
            }
            xassert(m > 0);
            /* set k := k + 1 and create new clique:
             * C[k] := (W intersect C[m]) union { i }, which covers all
             * edges incident to vertices from (W intersect C[m]) */
            k = glp_add_vertices(H, 1) - n;
            for (a = H->v[n + m]->in; a != NULL; a = a->h_next)
            {  j = a->tail->i;
               if (W->pos[j] != 0)
               {  /* vertex j is in both W and C[m]; include it in new
                   * clique C[k] */
                  glp_add_arc(H, j, n + k);
                  /* remove vertex j from W, since edge (i,j) will be
                   * covered by new clique C[k] */
                  if (W->pos[j] != W->size)
                  {  int jj = W->list[W->size];
                     W->list[W->pos[j]] = jj;
                     W->pos[jj] = W->pos[j];
                  }
                  W->size--, W->pos[j] = 0;
               }
            }
            /* include vertex i to new clique C[k] to cover edges (i,j)
             * incident to all vertices j just removed from W */
            glp_add_arc(H, i, n + k);
         }
      }
      /* free working arrays */
      xfree(W->list);
      xfree(W->pos);
      xfree(V->list);
      xfree(V->pos);
      /* return the number of cliques in the edge covering found */
      return H->nv - n;
}

/* eof */
