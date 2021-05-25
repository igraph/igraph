/* glpapi18.c (maximum clique problem) */

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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "glpapi.h"
#include "glpnet.h"

static void set_edge(int nv, unsigned char a[], int i, int j)
{     int k;
      xassert(1 <= j && j < i && i <= nv);
      k = ((i - 1) * (i - 2)) / 2 + (j - 1);
      a[k / CHAR_BIT] |=
         (unsigned char)(1 << ((CHAR_BIT - 1) - k % CHAR_BIT));
      return;
}

int glp_wclique_exact(glp_graph *G, int v_wgt, double *sol, int v_set)
{     /* find maximum weight clique with exact algorithm */
      glp_arc *e;
      int i, j, k, len, x, *w, *ind, ret = 0;
      unsigned char *a;
      double s, t;
      if (v_wgt >= 0 && v_wgt > G->v_size - (int)sizeof(double))
         xerror("glp_wclique_exact: v_wgt = %d; invalid parameter\n",
            v_wgt);
      if (v_set >= 0 && v_set > G->v_size - (int)sizeof(int))
         xerror("glp_wclique_exact: v_set = %d; invalid parameter\n",
            v_set);
      if (G->nv == 0)
      {  /* empty graph has only empty clique */
         if (sol != NULL) *sol = 0.0;
         return 0;
      }
      /* allocate working arrays */
      w = xcalloc(1+G->nv, sizeof(int));
      ind = xcalloc(1+G->nv, sizeof(int));
      len = G->nv; /* # vertices */
      len = len * (len - 1) / 2; /* # entries in lower triangle */
      len = (len + (CHAR_BIT - 1)) / CHAR_BIT; /* # bytes needed */
      a = xcalloc(len, sizeof(char));
      memset(a, 0, len * sizeof(char));
      /* determine vertex weights */
      s = 0.0;
      for (i = 1; i <= G->nv; i++)
      {  if (v_wgt >= 0)
         {  memcpy(&t, (char *)G->v[i]->data + v_wgt, sizeof(double));
            if (!(0.0 <= t && t <= (double)INT_MAX && t == floor(t)))
            {  ret = GLP_EDATA;
               goto done;
            }
            w[i] = (int)t;
         }
         else
            w[i] = 1;
         s += (double)w[i];
      }
      if (s > (double)INT_MAX)
      {  ret = GLP_EDATA;
         goto done;
      }
      /* build the adjacency matrix */
      for (i = 1; i <= G->nv; i++)
      {  for (e = G->v[i]->in; e != NULL; e = e->h_next)
         {  j = e->tail->i;
            /* there exists edge (j,i) in the graph */
            if (i > j) set_edge(G->nv, a, i, j);
         }
         for (e = G->v[i]->out; e != NULL; e = e->t_next)
         {  j = e->head->i;
            /* there exists edge (i,j) in the graph */
            if (i > j) set_edge(G->nv, a, i, j);
         }
      }
      /* find maximum weight clique in the graph */
      len = wclique(G->nv, w, a, ind);
      /* compute the clique weight */
      s = 0.0;
      for (k = 1; k <= len; k++)
      {  i = ind[k];
         xassert(1 <= i && i <= G->nv);
         s += (double)w[i];
      }
      if (sol != NULL) *sol = s;
      /* mark vertices included in the clique */
      if (v_set >= 0)
      {  x = 0;
         for (i = 1; i <= G->nv; i++)
            memcpy((char *)G->v[i]->data + v_set, &x, sizeof(int));
         x = 1;
         for (k = 1; k <= len; k++)
         {  i = ind[k];
            memcpy((char *)G->v[i]->data + v_set, &x, sizeof(int));
         }
      }
done: /* free working arrays */
      xfree(w);
      xfree(ind);
      xfree(a);
      return ret;
}

/* eof */
