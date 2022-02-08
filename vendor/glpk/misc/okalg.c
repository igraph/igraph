/* okalg.c (out-of-kilter algorithm) */

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

#include "env.h"
#include "okalg.h"

/***********************************************************************
*  NAME
*
*  okalg - out-of-kilter algorithm
*
*  SYNOPSIS
*
*  #include "okalg.h"
*  int okalg(int nv, int na, const int tail[], const int head[],
*     const int low[], const int cap[], const int cost[], int x[],
*     int pi[]);
*
*  DESCRIPTION
*
*  The routine okalg implements the out-of-kilter algorithm to find a
*  minimal-cost circulation in the specified flow network.
*
*  INPUT PARAMETERS
*
*  nv is the number of nodes, nv >= 0.
*
*  na is the number of arcs, na >= 0.
*
*  tail[a], a = 1,...,na, is the index of tail node of arc a.
*
*  head[a], a = 1,...,na, is the index of head node of arc a.
*
*  low[a], a = 1,...,na, is an lower bound to the flow through arc a.
*
*  cap[a], a = 1,...,na, is an upper bound to the flow through arc a,
*  which is the capacity of the arc.
*
*  cost[a], a = 1,...,na, is a per-unit cost of the flow through arc a.
*
*  NOTES
*
*  1. Multiple arcs are allowed, but self-loops are not allowed.
*
*  2. It is required that 0 <= low[a] <= cap[a] for all arcs.
*
*  3. Arc costs may have any sign.
*
*  OUTPUT PARAMETERS
*
*  x[a], a = 1,...,na, is optimal value of the flow through arc a.
*
*  pi[i], i = 1,...,nv, is Lagrange multiplier for flow conservation
*  equality constraint corresponding to node i (the node potential).
*
*  RETURNS
*
*  0  optimal circulation found;
*
*  1  there is no feasible circulation;
*
*  2  integer overflow occured;
*
*  3  optimality test failed (logic error).
*
*  REFERENCES
*
*  L.R.Ford, Jr., and D.R.Fulkerson, "Flows in Networks," The RAND
*  Corp., Report R-375-PR (August 1962), Chap. III "Minimal Cost Flow
*  Problems," pp.113-26. */

static int overflow(int u, int v)
{     /* check for integer overflow on computing u + v */
      if (u > 0 && v > 0 && u + v < 0) return 1;
      if (u < 0 && v < 0 && u + v > 0) return 1;
      return 0;
}

int okalg(int nv, int na, const int tail[], const int head[],
      const int low[], const int cap[], const int cost[], int x[],
      int pi[])
{     int a, aok, delta, i, j, k, lambda, pos1, pos2, s, t, temp, ret,
         *ptr, *arc, *link, *list;
      /* sanity checks */
      xassert(nv >= 0);
      xassert(na >= 0);
      for (a = 1; a <= na; a++)
      {  i = tail[a], j = head[a];
         xassert(1 <= i && i <= nv);
         xassert(1 <= j && j <= nv);
         xassert(i != j);
         xassert(0 <= low[a] && low[a] <= cap[a]);
      }
      /* allocate working arrays */
      ptr = xcalloc(1+nv+1, sizeof(int));
      arc = xcalloc(1+na+na, sizeof(int));
      link = xcalloc(1+nv, sizeof(int));
      list = xcalloc(1+nv, sizeof(int));
      /* ptr[i] := (degree of node i) */
      for (i = 1; i <= nv; i++)
         ptr[i] = 0;
      for (a = 1; a <= na; a++)
      {  ptr[tail[a]]++;
         ptr[head[a]]++;
      }
      /* initialize arc pointers */
      ptr[1]++;
      for (i = 1; i < nv; i++)
         ptr[i+1] += ptr[i];
      ptr[nv+1] = ptr[nv];
      /* build arc lists */
      for (a = 1; a <= na; a++)
      {  arc[--ptr[tail[a]]] = a;
         arc[--ptr[head[a]]] = a;
      }
      xassert(ptr[1] == 1);
      xassert(ptr[nv+1] == na+na+1);
      /* now the indices of arcs incident to node i are stored in
       * locations arc[ptr[i]], arc[ptr[i]+1], ..., arc[ptr[i+1]-1] */
      /* initialize arc flows and node potentials */
      for (a = 1; a <= na; a++)
         x[a] = 0;
      for (i = 1; i <= nv; i++)
         pi[i] = 0;
loop: /* main loop starts here */
      /* find out-of-kilter arc */
      aok = 0;
      for (a = 1; a <= na; a++)
      {  i = tail[a], j = head[a];
         if (overflow(cost[a], pi[i] - pi[j]))
         {  ret = 2;
            goto done;
         }
         lambda = cost[a] + (pi[i] - pi[j]);
         if (x[a] < low[a] || (lambda < 0 && x[a] < cap[a]))
         {  /* arc a = i->j is out of kilter, and we need to increase
             * the flow through this arc */
            aok = a, s = j, t = i;
            break;
         }
         if (x[a] > cap[a] || (lambda > 0 && x[a] > low[a]))
         {  /* arc a = i->j is out of kilter, and we need to decrease
             * the flow through this arc */
            aok = a, s = i, t = j;
            break;
         }
      }
      if (aok == 0)
      {  /* all arcs are in kilter */
         /* check for feasibility */
         for (a = 1; a <= na; a++)
         {  if (!(low[a] <= x[a] && x[a] <= cap[a]))
            {  ret = 3;
               goto done;
            }
         }
         for (i = 1; i <= nv; i++)
         {  temp = 0;
            for (k = ptr[i]; k < ptr[i+1]; k++)
            {  a = arc[k];
               if (tail[a] == i)
               {  /* a is outgoing arc */
                  temp += x[a];
               }
               else if (head[a] == i)
               {  /* a is incoming arc */
                  temp -= x[a];
               }
               else
                  xassert(a != a);
            }
            if (temp != 0)
            {  ret = 3;
               goto done;
            }
         }
         /* check for optimality */
         for (a = 1; a <= na; a++)
         {  i = tail[a], j = head[a];
            lambda = cost[a] + (pi[i] - pi[j]);
            if ((lambda > 0 && x[a] != low[a]) ||
                (lambda < 0 && x[a] != cap[a]))
            {  ret = 3;
               goto done;
            }
         }
         /* current circulation is optimal */
         ret = 0;
         goto done;
      }
      /* now we need to find a cycle (t, a, s, ..., t), which allows
       * increasing the flow along it, where a is the out-of-kilter arc
       * just found */
      /* link[i] = 0 means that node i is not labelled yet;
       * link[i] = a means that arc a immediately precedes node i */
      /* initially only node s is labelled */
      for (i = 1; i <= nv; i++)
         link[i] = 0;
      link[s] = aok, list[1] = s, pos1 = pos2 = 1;
      /* breadth first search */
      while (pos1 <= pos2)
      {  /* dequeue node i */
         i = list[pos1++];
         /* consider all arcs incident to node i */
         for (k = ptr[i]; k < ptr[i+1]; k++)
         {  a = arc[k];
            if (tail[a] == i)
            {  /* a = i->j is a forward arc from s to t */
               j = head[a];
               /* if node j has been labelled, skip the arc */
               if (link[j] != 0) continue;
               /* if the arc does not allow increasing the flow through
                * it, skip the arc */
               if (x[a] >= cap[a]) continue;
               if (overflow(cost[a], pi[i] - pi[j]))
               {  ret = 2;
                  goto done;
               }
               lambda = cost[a] + (pi[i] - pi[j]);
               if (lambda > 0 && x[a] >= low[a]) continue;
            }
            else if (head[a] == i)
            {  /* a = i<-j is a backward arc from s to t */
               j = tail[a];
               /* if node j has been labelled, skip the arc */
               if (link[j] != 0) continue;
               /* if the arc does not allow decreasing the flow through
                * it, skip the arc */
               if (x[a] <= low[a]) continue;
               if (overflow(cost[a], pi[j] - pi[i]))
               {  ret = 2;
                  goto done;
               }
               lambda = cost[a] + (pi[j] - pi[i]);
               if (lambda < 0 && x[a] <= cap[a]) continue;
            }
            else
               xassert(a != a);
            /* label node j and enqueue it */
            link[j] = a, list[++pos2] = j;
            /* check for breakthrough */
            if (j == t) goto brkt;
         }
      }
      /* NONBREAKTHROUGH */
      /* consider all arcs, whose one endpoint is labelled and other is
       * not, and determine maximal change of node potentials */
      delta = 0;
      for (a = 1; a <= na; a++)
      {  i = tail[a], j = head[a];
         if (link[i] != 0 && link[j] == 0)
         {  /* a = i->j, where node i is labelled, node j is not */
            if (overflow(cost[a], pi[i] - pi[j]))
            {  ret = 2;
               goto done;
            }
            lambda = cost[a] + (pi[i] - pi[j]);
            if (x[a] <= cap[a] && lambda > 0)
               if (delta == 0 || delta > + lambda) delta = + lambda;
         }
         else if (link[i] == 0 && link[j] != 0)
         {  /* a = j<-i, where node j is labelled, node i is not */
            if (overflow(cost[a], pi[i] - pi[j]))
            {  ret = 2;
               goto done;
            }
            lambda = cost[a] + (pi[i] - pi[j]);
            if (x[a] >= low[a] && lambda < 0)
               if (delta == 0 || delta > - lambda) delta = - lambda;
         }
      }
      if (delta == 0)
      {  /* there is no feasible circulation */
         ret = 1;
         goto done;
      }
      /* increase potentials of all unlabelled nodes */
      for (i = 1; i <= nv; i++)
      {  if (link[i] == 0)
         {  if (overflow(pi[i], delta))
            {  ret = 2;
               goto done;
            }
            pi[i] += delta;
         }
      }
      goto loop;
brkt: /* BREAKTHROUGH */
      /* walk through arcs of the cycle (t, a, s, ..., t) found in the
       * reverse order and determine maximal change of the flow */
      delta = 0;
      for (j = t;; j = i)
      {  /* arc a immediately precedes node j in the cycle */
         a = link[j];
         if (head[a] == j)
         {  /* a = i->j is a forward arc of the cycle */
            i = tail[a];
            lambda = cost[a] + (pi[i] - pi[j]);
            if (lambda > 0 && x[a] < low[a])
            {  /* x[a] may be increased until its lower bound */
               temp = low[a] - x[a];
            }
            else if (lambda <= 0 && x[a] < cap[a])
            {  /* x[a] may be increased until its upper bound */
               temp = cap[a] - x[a];
            }
            else
               xassert(a != a);
         }
         else if (tail[a] == j)
         {  /* a = i<-j is a backward arc of the cycle */
            i = head[a];
            lambda = cost[a] + (pi[j] - pi[i]);
            if (lambda < 0 && x[a] > cap[a])
            {  /* x[a] may be decreased until its upper bound */
               temp = x[a] - cap[a];
            }
            else if (lambda >= 0 && x[a] > low[a])
            {  /* x[a] may be decreased until its lower bound */
               temp = x[a] - low[a];
            }
            else
               xassert(a != a);
         }
         else
            xassert(a != a);
         if (delta == 0 || delta > temp) delta = temp;
         /* check for end of the cycle */
         if (i == t) break;
      }
      xassert(delta > 0);
      /* increase the flow along the cycle */
      for (j = t;; j = i)
      {  /* arc a immediately precedes node j in the cycle */
         a = link[j];
         if (head[a] == j)
         {  /* a = i->j is a forward arc of the cycle */
            i = tail[a];
            /* overflow cannot occur */
            x[a] += delta;
         }
         else if (tail[a] == j)
         {  /* a = i<-j is a backward arc of the cycle */
            i = head[a];
            /* overflow cannot occur */
            x[a] -= delta;
         }
         else
            xassert(a != a);
         /* check for end of the cycle */
         if (i == t) break;
      }
      goto loop;
done: /* free working arrays */
      xfree(ptr);
      xfree(arc);
      xfree(link);
      xfree(list);
      return ret;
}

/* eof */
