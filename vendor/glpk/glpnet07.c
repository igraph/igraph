/* glpnet07.c (Ford-Fulkerson algorithm) */

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

#include "glpenv.h"
#include "glpnet.h"

/***********************************************************************
*  NAME
*
*  ffalg - Ford-Fulkerson algorithm
*
*  SYNOPSIS
*
*  #include "glpnet.h"
*  void ffalg(int nv, int na, const int tail[], const int head[],
*     int s, int t, const int cap[], int x[], char cut[]);
*
*  DESCRIPTION
*
*  The routine ffalg implements the Ford-Fulkerson algorithm to find a
*  maximal flow in the specified flow network.
*
*  INPUT PARAMETERS
*
*  nv is the number of nodes, nv >= 2.
*
*  na is the number of arcs, na >= 0.
*
*  tail[a], a = 1,...,na, is the index of tail node of arc a.
*
*  head[a], a = 1,...,na, is the index of head node of arc a.
*
*  s is the source node index, 1 <= s <= nv.
*
*  t is the sink node index, 1 <= t <= nv, t != s.
*
*  cap[a], a = 1,...,na, is the capacity of arc a, cap[a] >= 0.
*
*  NOTE: Multiple arcs are allowed, but self-loops are not allowed.
*
*  OUTPUT PARAMETERS
*
*  x[a], a = 1,...,na, is optimal value of the flow through arc a.
*
*  cut[i], i = 1,...,nv, is 1 if node i is labelled, and 0 otherwise.
*  The set of arcs, whose one endpoint is labelled and other is not,
*  defines the minimal cut corresponding to the maximal flow found.
*  If the parameter cut is NULL, the cut information are not stored.
*
*  REFERENCES
*
*  L.R.Ford, Jr., and D.R.Fulkerson, "Flows in Networks," The RAND
*  Corp., Report R-375-PR (August 1962), Chap. I "Static Maximal Flow,"
*  pp.30-33. */

void ffalg(int nv, int na, const int tail[], const int head[],
      int s, int t, const int cap[], int x[], char cut[])
{     int a, delta, i, j, k, pos1, pos2, temp,
         *ptr, *arc, *link, *list;
      /* sanity checks */
      xassert(nv >= 2);
      xassert(na >= 0);
      xassert(1 <= s && s <= nv);
      xassert(1 <= t && t <= nv);
      xassert(s != t);
      for (a = 1; a <= na; a++)
      {  i = tail[a], j = head[a];
         xassert(1 <= i && i <= nv);
         xassert(1 <= j && j <= nv);
         xassert(i != j);
         xassert(cap[a] >= 0);
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
         locations arc[ptr[i]], arc[ptr[i]+1], ..., arc[ptr[i+1]-1] */
      /* initialize arc flows */
      for (a = 1; a <= na; a++)
         x[a] = 0;
loop: /* main loop starts here */
      /* build augmenting tree rooted at s */
      /* link[i] = 0 means that node i is not labelled yet;
         link[i] = a means that arc a immediately precedes node i */
      /* initially node s is labelled as the root */
      for (i = 1; i <= nv; i++)
         link[i] = 0;
      link[s] = -1, list[1] = s, pos1 = pos2 = 1;
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
                  it, skip the arc */
               if (x[a] == cap[a]) continue;
            }
            else if (head[a] == i)
            {  /* a = i<-j is a backward arc from s to t */
               j = tail[a];
               /* if node j has been labelled, skip the arc */
               if (link[j] != 0) continue;
               /* if the arc does not allow decreasing the flow through
                  it, skip the arc */
               if (x[a] == 0) continue;
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
      /* no augmenting path exists; current flow is maximal */
      /* store minimal cut information, if necessary */
      if (cut != NULL)
      {  for (i = 1; i <= nv; i++)
            cut[i] = (char)(link[i] != 0);
      }
      goto done;
brkt: /* BREAKTHROUGH */
      /* walk through arcs of the augmenting path (s, ..., t) found in
         the reverse order and determine maximal change of the flow */
      delta = 0;
      for (j = t; j != s; j = i)
      {  /* arc a immediately precedes node j in the path */
         a = link[j];
         if (head[a] == j)
         {  /* a = i->j is a forward arc of the cycle */
            i = tail[a];
            /* x[a] may be increased until its upper bound */
            temp = cap[a] - x[a];
         }
         else if (tail[a] == j)
         {  /* a = i<-j is a backward arc of the cycle */
            i = head[a];
            /* x[a] may be decreased until its lower bound */
            temp = x[a];
         }
         else
            xassert(a != a);
         if (delta == 0 || delta > temp) delta = temp;
      }
      xassert(delta > 0);
      /* increase the flow along the path */
      for (j = t; j != s; j = i)
      {  /* arc a immediately precedes node j in the path */
         a = link[j];
         if (head[a] == j)
         {  /* a = i->j is a forward arc of the cycle */
            i = tail[a];
            x[a] += delta;
         }
         else if (tail[a] == j)
         {  /* a = i<-j is a backward arc of the cycle */
            i = head[a];
            x[a] -= delta;
         }
         else
            xassert(a != a);
      }
      goto loop;
done: /* free working arrays */
      xfree(ptr);
      xfree(arc);
      xfree(link);
      xfree(list);
      return;
}

/* eof */
