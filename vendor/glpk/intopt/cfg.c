/* cfg.c (conflict graph) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2012-2013 Free Software Foundation, Inc.
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

#include "cfg.h"
#include "env.h"

/***********************************************************************
*  cfg_create_graph - create conflict graph
*
*  This routine creates the conflict graph, which initially is empty,
*  and returns a pointer to the graph descriptor.
*
*  The parameter n specifies the number of *all* variables in MIP, for
*  which the conflict graph will be built.
*
*  The parameter nv_max specifies maximal number of vertices in the
*  conflict graph. It should be the double number of binary variables
*  in corresponding MIP. */

CFG *cfg_create_graph(int n, int nv_max)
{     CFG *G;
      xassert(n >= 0);
      xassert(0 <= nv_max && nv_max <= n + n);
      G = talloc(1, CFG);
      G->n = n;
      G->pos = talloc(1+n, int);
      memset(&G->pos[1], 0, n * sizeof(int));
      G->neg = talloc(1+n, int);
      memset(&G->neg[1], 0, n * sizeof(int));
      G->pool = dmp_create_pool();
      G->nv_max = nv_max;
      G->nv = 0;
      G->ref = talloc(1+nv_max, int);
      G->vptr = talloc(1+nv_max, CFGVLE *);
      G->cptr = talloc(1+nv_max, CFGCLE *);
      return G;
}

/***********************************************************************
*  cfg_add_clique - add clique to conflict graph
*
*  This routine adds a clique to the conflict graph.
*
*  The parameter size specifies the clique size, size >= 2. Note that
*  any edge can be considered as a clique of size 2.
*
*  The array ind specifies vertices constituting the clique in elements
*  ind[k], 1 <= k <= size:
*
*  ind[k] = +j means a vertex of the conflict graph that corresponds to
*  original binary variable x[j], 1 <= j <= n.
*
*  ind[k] = -j means a vertex of the conflict graph that corresponds to
*  complement of original binary variable x[j], 1 <= j <= n.
*
*  Note that if both vertices for x[j] and (1 - x[j]) have appeared in
*  the conflict graph, the routine automatically adds an edge incident
*  to these vertices. */

static void add_edge(CFG *G, int v, int w)
{     /* add clique of size 2 */
      DMP *pool = G->pool;
      int nv = G->nv;
      CFGVLE **vptr = G->vptr;
      CFGVLE *vle;
      xassert(1 <= v && v <= nv);
      xassert(1 <= w && w <= nv);
      xassert(v != w);
      vle = dmp_talloc(pool, CFGVLE);
      vle->v = w;
      vle->next = vptr[v];
      vptr[v] = vle;
      vle = dmp_talloc(pool, CFGVLE);
      vle->v = v;
      vle->next = vptr[w];
      vptr[w] = vle;
      return;
}

void cfg_add_clique(CFG *G, int size, const int ind[])
{     int n = G->n;
      int *pos = G->pos;
      int *neg = G->neg;
      DMP *pool = G->pool;
      int nv_max = G->nv_max;
      int *ref = G->ref;
      CFGVLE **vptr = G->vptr;
      CFGCLE **cptr = G->cptr;
      int j, k, v;
      xassert(2 <= size && size <= nv_max);
      /* add new vertices to the conflict graph */
      for (k = 1; k <= size; k++)
      {  j = ind[k];
         if (j > 0)
         {  /* vertex corresponds to x[j] */
            xassert(1 <= j && j <= n);
            if (pos[j] == 0)
            {  /* no such vertex exists; add it */
               v = pos[j] = ++(G->nv);
               xassert(v <= nv_max);
               ref[v] = j;
               vptr[v] = NULL;
               cptr[v] = NULL;
               if (neg[j] != 0)
               {  /* now both vertices for x[j] and (1 - x[j]) exist */
                  add_edge(G, v, neg[j]);
               }
            }
         }
         else
         {  /* vertex corresponds to (1 - x[j]) */
            j = -j;
            xassert(1 <= j && j <= n);
            if (neg[j] == 0)
            {  /* no such vertex exists; add it */
               v = neg[j] = ++(G->nv);
               xassert(v <= nv_max);
               ref[v] = j;
               vptr[v] = NULL;
               cptr[v] = NULL;
               if (pos[j] != 0)
               {  /* now both vertices for x[j] and (1 - x[j]) exist */
                  add_edge(G, v, pos[j]);
               }
            }
         }
      }
      /* add specified clique to the conflict graph */
      if (size == 2)
         add_edge(G,
            ind[1] > 0 ? pos[+ind[1]] : neg[-ind[1]],
            ind[2] > 0 ? pos[+ind[2]] : neg[-ind[2]]);
      else
      {  CFGVLE *vp, *vle;
         CFGCLE *cle;
         /* build list of clique vertices */
         vp = NULL;
         for (k = 1; k <= size; k++)
         {  vle = dmp_talloc(pool, CFGVLE);
            vle->v = ind[k] > 0 ? pos[+ind[k]] : neg[-ind[k]];
            vle->next = vp;
            vp = vle;
         }
         /* attach the clique to all its vertices */
         for (k = 1; k <= size; k++)
         {  cle = dmp_talloc(pool, CFGCLE);
            cle->vptr = vp;
            v = ind[k] > 0 ? pos[+ind[k]] : neg[-ind[k]];
            cle->next = cptr[v];
            cptr[v] = cle;
         }
      }
      return;
}

/***********************************************************************
*  cfg_get_adjacent - get vertices adjacent to specified vertex
*
*  This routine stores numbers of all vertices adjacent to specified
*  vertex v of the conflict graph in locations ind[1], ..., ind[len],
*  and returns len, 1 <= len <= nv-1, where nv is the total number of
*  vertices in the conflict graph.
*
*  Note that the conflict graph defined by this routine has neither
*  self-loops nor multiple edges. */

int cfg_get_adjacent(CFG *G, int v, int ind[])
{     int nv = G->nv;
      int *ref = G->ref;
      CFGVLE **vptr = G->vptr;
      CFGCLE **cptr = G->cptr;
      CFGVLE *vle;
      CFGCLE *cle;
      int k, w, len;
      xassert(1 <= v && v <= nv);
      len = 0;
      /* walk thru the list of adjacent vertices */
      for (vle = vptr[v]; vle != NULL; vle = vle->next)
      {  w = vle->v;
         xassert(1 <= w && w <= nv);
         xassert(w != v);
         if (ref[w] > 0)
         {  ind[++len] = w;
            ref[w] = -ref[w];
         }
      }
      /* walk thru the list of incident cliques */
      for (cle = cptr[v]; cle != NULL; cle = cle->next)
      {  /* walk thru the list of clique vertices */
         for (vle = cle->vptr; vle != NULL; vle = vle->next)
         {  w = vle->v;
            xassert(1 <= w && w <= nv);
            if (w != v && ref[w] > 0)
            {  ind[++len] = w;
               ref[w] = -ref[w];
            }
         }
      }
      xassert(1 <= len && len < nv);
      /* unmark vertices included in the resultant adjacency list */
      for (k = 1; k <= len; k++)
      {  w = ind[k];
         ref[w] = -ref[w];
      }
      return len;
}

/***********************************************************************
*  cfg_expand_clique - expand specified clique to maximal clique
*
*  Given some clique in the conflict graph this routine expands it to
*  a maximal clique by including in it new vertices.
*
*  On entry vertex indices constituting the initial clique should be
*  stored in locations c_ind[1], ..., c_ind[c_len], where c_len is the
*  initial clique size. On exit the routine stores new vertex indices
*  to locations c_ind[c_len+1], ..., c_ind[c_len'], where c_len' is the
*  size of the maximal clique found, and returns c_len'.
*
*  ALGORITHM
*
*  Let G = (V, E) be a graph, C within V be a current clique to be
*  expanded, and D within V \ C be a subset of vertices adjacent to all
*  vertices from C. On every iteration the routine chooses some vertex
*  v in D, includes it into C, and removes from D the vertex v as well
*  as all vertices not adjacent to v. Initially C is empty and D = V.
*  Iterations repeat until D becomes an empty set. Obviously, the final
*  set C is a maximal clique in G.
*
*  Now let C0 be an initial clique, and we want C0 to be a subset of
*  the final maximal clique C. To provide this condition the routine
*  starts constructing C by choosing only such vertices v in D, which
*  are in C0, until all vertices from C0 have been included in C. May
*  note that if on some iteration C0 \ C is non-empty (i.e. if not all
*  vertices from C0 have been included in C), C0 \ C is a subset of D,
*  because C0 is a clique. */

static int intersection(int d_len, int d_ind[], int d_pos[], int len,
      const int ind[])
{     /* compute intersection D := D inter W, where W is some specified
       * set of vertices */
      int k, t, v, new_len;
      /* walk thru vertices in W and mark vertices in D */
      for (t = 1; t <= len; t++)
      {  /* v in W */
         v = ind[t];
         /* determine position of v in D */
         k = d_pos[v];
         if (k != 0)
         {  /* v in D */
            xassert(d_ind[k] == v);
            /* mark v to keep it in D */
            d_ind[k] = -v;
         }
      }
      /* remove all unmarked vertices from D */
      new_len = 0;
      for (k = 1; k <= d_len; k++)
      {  /* v in D */
         v = d_ind[k];
         if (v < 0)
         {  /* v is marked; keep it */
            v = -v;
            new_len++;
            d_ind[new_len] = v;
            d_pos[v] = new_len;
         }
         else
         {  /* v is not marked; remove it */
            d_pos[v] = 0;
         }
      }
      return new_len;
}

int cfg_expand_clique(CFG *G, int c_len, int c_ind[])
{     int nv = G->nv;
      int d_len, *d_ind, *d_pos, len, *ind;
      int k, v;
      xassert(0 <= c_len && c_len <= nv);
      /* allocate working arrays */
      d_ind = talloc(1+nv, int);
      d_pos = talloc(1+nv, int);
      ind = talloc(1+nv, int);
      /* initialize C := 0, D := V */
      d_len = nv;
      for (k = 1; k <= nv; k++)
         d_ind[k] = d_pos[k] = k;
      /* expand C by vertices of specified initial clique C0 */
      for (k = 1; k <= c_len; k++)
      {  /* v in C0 */
         v = c_ind[k];
         xassert(1 <= v && v <= nv);
         /* since C0 is clique, v should be in D */
         xassert(d_pos[v] != 0);
         /* W := set of vertices adjacent to v */
         len = cfg_get_adjacent(G, v, ind);
         /* D := D inter W */
         d_len = intersection(d_len, d_ind, d_pos, len, ind);
         /* since v not in W, now v should be not in D */
         xassert(d_pos[v] == 0);
      }
      /* expand C by some other vertices until D is empty */
      while (d_len > 0)
      {  /* v in D */
         v = d_ind[1];
         xassert(1 <= v && v <= nv);
         /* note that v is adjacent to all vertices in C (by design),
          * so add v to C */
         c_ind[++c_len] = v;
         /* W := set of vertices adjacent to v */
         len = cfg_get_adjacent(G, v, ind);
         /* D := D inter W */
         d_len = intersection(d_len, d_ind, d_pos, len, ind);
         /* since v not in W, now v should be not in D */
         xassert(d_pos[v] == 0);
      }
      /* free working arrays */
      tfree(d_ind);
      tfree(d_pos);
      tfree(ind);
      /* bring maximal clique to calling routine */
      return c_len;
}

/***********************************************************************
*  cfg_check_clique - check clique in conflict graph
*
*  This routine checks that vertices of the conflict graph specified
*  in locations c_ind[1], ..., c_ind[c_len] constitute a clique.
*
*  NOTE: for testing/debugging only. */

void cfg_check_clique(CFG *G, int c_len, const int c_ind[])
{     int nv = G->nv;
      int k, kk, v, w, len, *ind;
      char *flag;
      ind = talloc(1+nv, int);
      flag = talloc(1+nv, char);
      memset(&flag[1], 0, nv);
      /* walk thru clique vertices */
      xassert(c_len >= 0);
      for (k = 1; k <= c_len; k++)
      {  /* get clique vertex v */
         v = c_ind[k];
         xassert(1 <= v && v <= nv);
         /* get vertices adjacent to vertex v */
         len = cfg_get_adjacent(G, v, ind);
         for (kk = 1; kk <= len; kk++)
         {  w = ind[kk];
            xassert(1 <= w && w <= nv);
            xassert(w != v);
            flag[w] = 1;
         }
         /* check that all clique vertices other than v are adjacent
            to v */
         for (kk = 1; kk <= c_len; kk++)
         {  w = c_ind[kk];
            xassert(1 <= w && w <= nv);
            if (w != v)
               xassert(flag[w]);
         }
         /* reset vertex flags */
         for (kk = 1; kk <= len; kk++)
            flag[ind[kk]] = 0;
      }
      tfree(ind);
      tfree(flag);
      return;
}

/***********************************************************************
*  cfg_delete_graph - delete conflict graph
*
*  This routine deletes the conflict graph by freeing all the memory
*  allocated to this program object. */

void cfg_delete_graph(CFG *G)
{     tfree(G->pos);
      tfree(G->neg);
      dmp_delete_pool(G->pool);
      tfree(G->ref);
      tfree(G->vptr);
      tfree(G->cptr);
      tfree(G);
      return;
}

/* eof */
