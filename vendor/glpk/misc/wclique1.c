/* wclique1.c (maximum weight clique, greedy heuristic) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2012-2018 Free Software Foundation, Inc.
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
#include "wclique1.h"

/***********************************************************************
*  NAME
*
*  wclique1 - find maximum weight clique with greedy heuristic
*
*  SYNOPSIS
*
*  #include "wclique1.h"
*  int wclique1(int n, const double w[],
*     int (*func)(void *info, int i, int ind[]), void *info, int c[]);
*
*  DESCRIPTION
*
*  The routine wclique1 implements a sequential greedy heuristic to
*  find maximum weight clique in a given (undirected) graph G = (V, E).
*
*  The parameter n specifies the number of vertices |V| in the graph,
*  n >= 0.
*
*  The array w specifies vertex weights in locations w[i], i = 1,...,n.
*  All weights must be non-negative.
*
*  The formal routine func specifies the graph. For a given vertex i,
*  1 <= i <= n, it stores indices of all vertices adjacent to vertex i
*  in locations ind[1], ..., ind[deg], where deg is the degree of
*  vertex i, 0 <= deg < n, returned on exit. Note that self-loops and
*  multiple edges are not allowed.
*
*  The parameter info is a cookie passed to the routine func.
*
*  On exit the routine wclique1 stores vertex indices included in
*  the clique found to locations c[1], ..., c[size], where size is the
*  clique size returned by the routine, 0 <= size <= n.
*
*  RETURNS
*
*  The routine wclique1 returns the size of the clique found. */

struct vertex { int i; double cw; };

static int CDECL fcmp(const void *xx, const void *yy)
{     const struct vertex *x = xx, *y = yy;
      if (x->cw > y->cw) return -1;
      if (x->cw < y->cw) return +1;
      return 0;
}

int wclique1(int n, const double w[],
      int (*func)(void *info, int i, int ind[]), void *info, int c[])
{     struct vertex *v_list;
      int deg, c_size, d_size, i, j, k, kk, l, *ind, *c_list, *d_list,
         size = 0;
      double c_wght, d_wght, *sw, best = 0.0;
      char *d_flag, *skip;
      /* perform sanity checks */
      xassert(n >= 0);
      for (i = 1; i <= n; i++)
         xassert(w[i] >= 0.0);
      /* if the graph is empty, nothing to do */
      if (n == 0) goto done;
      /* allocate working arrays */
      ind = xcalloc(1+n, sizeof(int));
      v_list = xcalloc(1+n, sizeof(struct vertex));
      c_list = xcalloc(1+n, sizeof(int));
      d_list = xcalloc(1+n, sizeof(int));
      d_flag = xcalloc(1+n, sizeof(char));
      skip = xcalloc(1+n, sizeof(char));
      sw = xcalloc(1+n, sizeof(double));
      /* build the vertex list */
      for (i = 1; i <= n; i++)
      {  v_list[i].i = i;
         /* compute the cumulative weight of each vertex i, which is
          * cw[i] = w[i] + sum{j : (i,j) in E} w[j] */
         v_list[i].cw = w[i];
         deg = func(info, i, ind);
         xassert(0 <= deg && deg < n);
         for (k = 1; k <= deg; k++)
         {  j = ind[k];
            xassert(1 <= j && j <= n && j != i);
            v_list[i].cw += w[j];
         }
      }
      /* sort the vertex list to access vertices in descending order of
       * cumulative weights */
      qsort(&v_list[1], n, sizeof(struct vertex), fcmp);
      /* initially all vertices are unmarked */
      memset(&skip[1], 0, sizeof(char) * n);
      /* clear flags of all vertices */
      memset(&d_flag[1], 0, sizeof(char) * n);
      /* look through all vertices of the graph */
      for (l = 1; l <= n; l++)
      {  /* take vertex i */
         i = v_list[l].i;
         /* if this vertex was already included in one of previosuly
          * constructed cliques, skip it */
         if (skip[i]) continue;
         /* use vertex i as the initial clique vertex */
         c_size = 1;    /* size of current clique */
         c_list[1] = i; /* list of vertices in current clique */
         c_wght = w[i]; /* weight of current clique */
         /* determine the candidate set D = { j : (i,j) in E } */
         d_size = func(info, i, d_list);
         xassert(0 <= d_size && d_size < n);
         d_wght = 0.0;  /* weight of set D */
         for (k = 1; k <= d_size; k++)
         {  j = d_list[k];
            xassert(1 <= j && j <= n && j != i);
            xassert(!d_flag[j]);
            d_flag[j] = 1;
            d_wght += w[j];
         }
         /* check an upper bound to the final clique weight */
         if (c_wght + d_wght < best + 1e-5 * (1.0 + fabs(best)))
         {  /* skip constructing the current clique */
            goto next;
         }
         /* compute the summary weight of each vertex i in D, which is
          * sw[i] = w[i] + sum{j in D and (i,j) in E} w[j] */
         for (k = 1; k <= d_size; k++)
         {  i = d_list[k];
            sw[i] = w[i];
            /* consider vertices adjacent to vertex i */
            deg = func(info, i, ind);
            xassert(0 <= deg && deg < n);
            for (kk = 1; kk <= deg; kk++)
            {  j = ind[kk];
               xassert(1 <= j && j <= n && j != i);
               if (d_flag[j]) sw[i] += w[j];
            }
         }
         /* grow the current clique by adding vertices from D */
         while (d_size > 0)
         {  /* check an upper bound to the final clique weight */
            if (c_wght + d_wght < best + 1e-5 * (1.0 + fabs(best)))
            {  /* skip constructing the current clique */
               goto next;
            }
            /* choose vertex i in D having maximal summary weight */
            i = d_list[1];
            for (k = 2; k <= d_size; k++)
            {  j = d_list[k];
               if (sw[i] < sw[j]) i = j;
            }
            /* include vertex i in the current clique */
            c_size++;
            c_list[c_size] = i;
            c_wght += w[i];
            /* remove all vertices not adjacent to vertex i, including
             * vertex i itself, from the candidate set D */
            deg = func(info, i, ind);
            xassert(0 <= deg && deg < n);
            for (k = 1; k <= deg; k++)
            {  j = ind[k];
               xassert(1 <= j && j <= n && j != i);
               /* vertex j is adjacent to vertex i */
               if (d_flag[j])
               {  xassert(d_flag[j] == 1);
                  /* mark vertex j to keep it in D */
                  d_flag[j] = 2;
               }
            }
            kk = d_size, d_size = 0;
            for (k = 1; k <= kk; k++)
            {  j = d_list[k];
               if (d_flag[j] == 1)
               {  /* remove vertex j from D */
                  d_flag[j] = 0;
                  d_wght -= w[j];
               }
               else if (d_flag[j] == 2)
               {  /* keep vertex j in D */
                  d_list[++d_size] = j;
                  d_flag[j] = 1;
               }
               else
                  xassert(d_flag != d_flag);
            }
         }
         /* the current clique has been completely constructed */
         if (best < c_wght)
         {  best = c_wght;
            size = c_size;
            xassert(1 <= size && size <= n);
            memcpy(&c[1], &c_list[1], size * sizeof(int));
         }
next:    /* mark the current clique vertices in order not to use them
          * as initial vertices anymore */
         for (k = 1; k <= c_size; k++)
            skip[c_list[k]] = 1;
         /* set D can be non-empty, so clean up vertex flags */
         for (k = 1; k <= d_size; k++)
            d_flag[d_list[k]] = 0;
      }
      /* free working arrays */
      xfree(ind);
      xfree(v_list);
      xfree(c_list);
      xfree(d_list);
      xfree(d_flag);
      xfree(skip);
      xfree(sw);
done: /* return to the calling program */
      return size;
}

/**********************************************************************/

#ifdef GLP_TEST
#include "glpk.h"
#include "rng.h"

typedef struct { double w; } v_data;

#define weight(v) (((v_data *)((v)->data))->w)

glp_graph *G;

char *flag;

int func(void *info, int i, int ind[])
{     glp_arc *e;
      int j, k, deg = 0;
      xassert(info == NULL);
      xassert(1 <= i && i <= G->nv);
      /* look through incoming arcs */
      for (e = G->v[i]->in; e != NULL; e = e->h_next)
      {  j = e->tail->i; /* j->i */
         if (j != i && !flag[j]) ind[++deg] = j, flag[j] = 1;
      }
      /* look through outgoing arcs */
      for (e = G->v[i]->out; e != NULL; e = e->t_next)
      {  j = e->head->i; /* i->j */
         if (j != i && !flag[j]) ind[++deg] = j, flag[j] = 1;
      }
      /* clear the flag array */
      xassert(deg < G->nv);
      for (k = 1; k <= deg; k++) flag[ind[k]] = 0;
      return deg;
}

int main(int argc, char *argv[])
{     RNG *rand;
      int i, k, kk, size, *c, *ind, deg;
      double *w, sum, t;
      /* read graph in DIMACS format */
      G = glp_create_graph(sizeof(v_data), 0);
      xassert(argc == 2);
      xassert(glp_read_ccdata(G, offsetof(v_data, w), argv[1]) == 0);
      /* print the number of connected components */
      xprintf("nc = %d\n", glp_weak_comp(G, -1));
      /* assign random weights unformly distributed in [1,100] */
      w = xcalloc(1+G->nv, sizeof(double));
      rand = rng_create_rand();
      for (i = 1; i <= G->nv; i++)
#if 0
         w[i] = weight(G->v[i]) = 1.0;
#else
         w[i] = weight(G->v[i]) = rng_unif_rand(rand, 100) + 1;
#endif
      /* write graph in DIMACS format */
      xassert(glp_write_ccdata(G, offsetof(v_data, w), "graph") == 0);
      /* find maximum weight clique */
      c = xcalloc(1+G->nv, sizeof(int));
      flag = xcalloc(1+G->nv, sizeof(char));
      memset(&flag[1], 0, G->nv);
      t = xtime();
      size = wclique1(G->nv, w, func, NULL, c);
      xprintf("Time used: %.1f s\n", xdifftime(xtime(), t));
      /* check the clique found */
      ind = xcalloc(1+G->nv, sizeof(int));
      for (k = 1; k <= size; k++)
      {  i = c[k];
         deg = func(NULL, i, ind);
         for (kk = 1; kk <= size; kk++)
            flag[c[kk]] = 1;
         flag[i] = 0;
         for (kk = 1; kk <= deg; kk++)
            flag[ind[kk]] = 0;
         for (kk = 1; kk <= size; kk++)
            xassert(flag[c[kk]] == 0);
      }
      /* compute the clique weight */
      sum = 0.0;
      for (i = 1; i <= size; i++)
         sum += w[c[i]];
      xprintf("size = %d; sum = %g\n", size, sum);
      return 0;
}
#endif

/* eof */
