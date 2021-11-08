/* asnhall.c (find bipartite matching of maximum cardinality) */

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
#include "mc21a.h"

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

/* eof */
