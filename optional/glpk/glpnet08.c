/* glpnet08.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Two subroutines sub() and wclique() below are intended to find a
*  maximum weight clique in a given undirected graph. These subroutines
*  are slightly modified version of the program WCLIQUE developed by
*  Patric Ostergard <http://www.tcs.hut.fi/~pat/wclique.html> and based
*  on ideas from the article "P. R. J. Ostergard, A new algorithm for
*  the maximum-weight clique problem, submitted for publication", which
*  in turn is a generalization of the algorithm for unweighted graphs
*  presented in "P. R. J. Ostergard, A fast algorithm for the maximum
*  clique problem, submitted for publication".
*
*  USED WITH PERMISSION OF THE AUTHOR OF THE ORIGINAL CODE.
*
*  Changes were made by Andrew Makhorin <mao@gnu.org>.
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

#include "glpenv.h"
#include "glpnet.h"

/***********************************************************************
*  NAME
*
*  wclique - find maximum weight clique with Ostergard's algorithm
*
*  SYNOPSIS
*
*  int wclique(int n, const int w[], const unsigned char a[],
*     int ind[]);
*
*  DESCRIPTION
*
*  The routine wclique finds a maximum weight clique in an undirected
*  graph with Ostergard's algorithm.
*
*  INPUT PARAMETERS
*
*  n is the number of vertices, n > 0.
*
*  w[i], i = 1,...,n, is a weight of vertex i.
*
*  a[*] is the strict (without main diagonal) lower triangle of the
*  graph adjacency matrix in packed format.
*
*  OUTPUT PARAMETER
*
*  ind[k], k = 1,...,size, is the number of a vertex included in the
*  clique found, 1 <= ind[k] <= n, where size is the number of vertices
*  in the clique returned on exit.
*
*  RETURNS
*
*  The routine returns the clique size, i.e. the number of vertices in
*  the clique. */

struct csa
{     /* common storage area */
      int n;
      /* number of vertices */
      const int *wt; /* int wt[0:n-1]; */
      /* weights */
      const unsigned char *a;
      /* adjacency matrix (packed lower triangle without main diag.) */
      int record;
      /* weight of best clique */
      int rec_level;
      /* number of vertices in best clique */
      int *rec; /* int rec[0:n-1]; */
      /* best clique so far */
      int *clique; /* int clique[0:n-1]; */
      /* table for pruning */
      int *set; /* int set[0:n-1]; */
      /* current clique */
};

#define n         (csa->n)
#define wt        (csa->wt)
#define a         (csa->a)
#define record    (csa->record)
#define rec_level (csa->rec_level)
#define rec       (csa->rec)
#define clique    (csa->clique)
#define set       (csa->set)

#if 0
static int is_edge(struct csa *csa, int i, int j)
{     /* if there is arc (i,j), the routine returns true; otherwise
         false; 0 <= i, j < n */
      int k;
      xassert(0 <= i && i < n);
      xassert(0 <= j && j < n);
      if (i == j) return 0;
      if (i < j) k = i, i = j, j = k;
      k = (i * (i - 1)) / 2 + j;
      return a[k / CHAR_BIT] &
         (unsigned char)(1 << ((CHAR_BIT - 1) - k % CHAR_BIT));
}
#else
#define is_edge(csa, i, j) ((i) == (j) ? 0 : \
      (i) > (j) ? is_edge1(i, j) : is_edge1(j, i))
#define is_edge1(i, j) is_edge2(((i) * ((i) - 1)) / 2 + (j))
#define is_edge2(k) (a[(k) / CHAR_BIT] & \
      (unsigned char)(1 << ((CHAR_BIT - 1) - (k) % CHAR_BIT)))
#endif

static void sub(struct csa *csa, int ct, int table[], int level,
      int weight, int l_weight)
{     int i, j, k, curr_weight, left_weight, *p1, *p2, *newtable;
      newtable = xcalloc(n, sizeof(int));
      if (ct <= 0)
      {  /* 0 or 1 elements left; include these */
         if (ct == 0)
         {  set[level++] = table[0];
            weight += l_weight;
         }
         if (weight > record)
         {  record = weight;
            rec_level = level;
            for (i = 0; i < level; i++) rec[i] = set[i];
         }
         goto done;
      }
      for (i = ct; i >= 0; i--)
      {  if ((level == 0) && (i < ct)) goto done;
         k = table[i];
         if ((level > 0) && (clique[k] <= (record - weight)))
            goto done; /* prune */
         set[level] = k;
         curr_weight = weight + wt[k];
         l_weight -= wt[k];
         if (l_weight <= (record - curr_weight))
            goto done; /* prune */
         p1 = newtable;
         p2 = table;
         left_weight = 0;
         while (p2 < table + i)
         {  j = *p2++;
            if (is_edge(csa, j, k))
            {  *p1++ = j;
               left_weight += wt[j];
            }
         }
         if (left_weight <= (record - curr_weight)) continue;
         sub(csa, p1 - newtable - 1, newtable, level + 1, curr_weight,
            left_weight);
      }
done: xfree(newtable);
      return;
}

int wclique(int _n, const int w[], const unsigned char _a[], int ind[])
{     struct csa _csa, *csa = &_csa;
      int i, j, p, max_wt, max_nwt, wth, *used, *nwt, *pos;
      glp_long timer;
      n = _n;
      xassert(n > 0);
      wt = &w[1];
      a = _a;
      record = 0;
      rec_level = 0;
      rec = &ind[1];
      clique = xcalloc(n, sizeof(int));
      set = xcalloc(n, sizeof(int));
      used = xcalloc(n, sizeof(int));
      nwt = xcalloc(n, sizeof(int));
      pos = xcalloc(n, sizeof(int));
      /* start timer */
      timer = xtime();
      /* order vertices */
      for (i = 0; i < n; i++)
      {  nwt[i] = 0;
         for (j = 0; j < n; j++)
            if (is_edge(csa, i, j)) nwt[i] += wt[j];
      }
      for (i = 0; i < n; i++)
         used[i] = 0;
      for (i = n-1; i >= 0; i--)
      {  max_wt = -1;
         max_nwt = -1;
         for (j = 0; j < n; j++)
         {  if ((!used[j]) && ((wt[j] > max_wt) || (wt[j] == max_wt
               && nwt[j] > max_nwt)))
            {  max_wt = wt[j];
               max_nwt = nwt[j];
               p = j;
            }
         }
         pos[i] = p;
         used[p] = 1;
         for (j = 0; j < n; j++)
            if ((!used[j]) && (j != p) && (is_edge(csa, p, j)))
               nwt[j] -= wt[p];
      }
      /* main routine */
      wth = 0;
      for (i = 0; i < n; i++)
      {  wth += wt[pos[i]];
         sub(csa, i, pos, 0, 0, wth);
         clique[pos[i]] = record;
         if (xdifftime(xtime(), timer) >= 5.0 - 0.001)
         {  /* print current record and reset timer */
            xprintf("level = %d (%d); best = %d\n", i+1, n, record);
            timer = xtime();
         }
      }
      xfree(clique);
      xfree(set);
      xfree(used);
      xfree(nwt);
      xfree(pos);
      /* return the solution found */
      for (i = 1; i <= rec_level; i++) ind[i]++;
      return rec_level;
}

#undef n
#undef wt
#undef a
#undef record
#undef rec_level
#undef rec
#undef clique
#undef set

/* eof */
