/* glpqmd.c (quotient minimum degree algorithm) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  THIS CODE IS THE RESULT OF TRANSLATION OF THE FORTRAN SUBROUTINES
*  GENQMD, QMDRCH, QMDQT, QMDUPD, AND QMDMRG FROM THE BOOK:
*
*  ALAN GEORGE, JOSEPH W-H LIU. COMPUTER SOLUTION OF LARGE SPARSE
*  POSITIVE DEFINITE SYSTEMS. PRENTICE-HALL, 1981.
*
*  THE TRANSLATION HAS BEEN DONE WITH THE PERMISSION OF THE AUTHORS
*  OF THE ORIGINAL FORTRAN SUBROUTINES: ALAN GEORGE AND JOSEPH LIU,
*  UNIVERSITY OF WATERLOO, WATERLOO, ONTARIO, CANADA.
*
*  The translation was made by Andrew Makhorin <mao@gnu.org>.
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

#include "glpqmd.h"

/***********************************************************************
*  NAME
*
*  genqmd - GENeral Quotient Minimum Degree algorithm
*
*  SYNOPSIS
*
*  #include "glpqmd.h"
*  void genqmd(int *neqns, int xadj[], int adjncy[], int perm[],
*     int invp[], int deg[], int marker[], int rchset[], int nbrhd[],
*     int qsize[], int qlink[], int *nofsub);
*
*  PURPOSE
*
*  This routine implements the minimum degree algorithm. It makes use
*  of the implicit representation of the elimination graph by quotient
*  graphs, and the notion of indistinguishable nodes.
*
*  CAUTION
*
*  The adjancy vector adjncy will be destroyed.
*
*  INPUT PARAMETERS
*
*  neqns  - number of equations;
*  (xadj, adjncy) -
*           the adjancy structure.
*
*  OUTPUT PARAMETERS
*
*  perm   - the minimum degree ordering;
*  invp   - the inverse of perm.
*
*  WORKING PARAMETERS
*
*  deg    - the degree vector. deg[i] is negative means node i has been
*           numbered;
*  marker - a marker vector, where marker[i] is negative means node i
*           has been merged with another nodeand thus can be ignored;
*  rchset - vector used for the reachable set;
*  nbrhd  - vector used for neighborhood set;
*  qsize  - vector used to store the size of indistinguishable
*           supernodes;
*  qlink  - vector used to store indistinguishable nodes, i, qlink[i],
*           qlink[qlink[i]], ... are the members of the supernode
*           represented by i.
*
*  PROGRAM SUBROUTINES
*
*  qmdrch, qmdqt, qmdupd.
***********************************************************************/

void genqmd(int *_neqns, int xadj[], int adjncy[], int perm[],
      int invp[], int deg[], int marker[], int rchset[], int nbrhd[],
      int qsize[], int qlink[], int *_nofsub)
{     int inode, ip, irch, j, mindeg, ndeg, nhdsze, node, np, num,
         nump1, nxnode, rchsze, search, thresh;
#     define neqns  (*_neqns)
#     define nofsub (*_nofsub)
      /* Initialize degree vector and other working variables. */
      mindeg = neqns;
      nofsub = 0;
      for (node = 1; node <= neqns; node++)
      {  perm[node] = node;
         invp[node] = node;
         marker[node] = 0;
         qsize[node] = 1;
         qlink[node] = 0;
         ndeg = xadj[node+1] - xadj[node];
         deg[node] = ndeg;
         if (ndeg < mindeg) mindeg = ndeg;
      }
      num = 0;
      /* Perform threshold search to get a node of min degree.
         Variable search point to where search should start. */
s200: search = 1;
      thresh = mindeg;
      mindeg = neqns;
s300: nump1 = num + 1;
      if (nump1 > search) search = nump1;
      for (j = search; j <= neqns; j++)
      {  node = perm[j];
         if (marker[node] >= 0)
         {  ndeg = deg[node];
            if (ndeg <= thresh) goto s500;
            if (ndeg < mindeg) mindeg = ndeg;
         }
      }
      goto s200;
      /* Node has minimum degree. Find its reachable sets by calling
         qmdrch. */
s500: search = j;
      nofsub += deg[node];
      marker[node] = 1;
      qmdrch(&node, xadj, adjncy, deg, marker, &rchsze, rchset, &nhdsze,
         nbrhd);
      /* Eliminate all nodes indistinguishable from node. They are given
         by node, qlink[node], ... . */
      nxnode = node;
s600: num++;
      np = invp[nxnode];
      ip = perm[num];
      perm[np] = ip;
      invp[ip] = np;
      perm[num] = nxnode;
      invp[nxnode] = num;
      deg[nxnode] = -1;
      nxnode = qlink[nxnode];
      if (nxnode > 0) goto s600;
      if (rchsze > 0)
      {  /* Update the degrees of the nodes in the reachable set and
            identify indistinguishable nodes. */
         qmdupd(xadj, adjncy, &rchsze, rchset, deg, qsize, qlink,
            marker, &rchset[rchsze+1], &nbrhd[nhdsze+1]);
         /* Reset marker value of nodes in reach set. Update threshold
            value for cyclic search. Also call qmdqt to form new
            quotient graph. */
         marker[node] = 0;
         for (irch = 1; irch <= rchsze; irch++)
         {  inode = rchset[irch];
            if (marker[inode] >= 0)
            {  marker[inode] = 0;
               ndeg = deg[inode];
               if (ndeg < mindeg) mindeg = ndeg;
               if (ndeg <= thresh)
               {  mindeg = thresh;
                  thresh = ndeg;
                  search = invp[inode];
               }
            }
         }
         if (nhdsze > 0)
            qmdqt(&node, xadj, adjncy, marker, &rchsze, rchset, nbrhd);
      }
      if (num < neqns) goto s300;
      return;
#     undef neqns
#     undef nofsub
}

/***********************************************************************
*  NAME
*
*  qmdrch - Quotient MD ReaCHable set
*
*  SYNOPSIS
*
*  #include "glpqmd.h"
*  void qmdrch(int *root, int xadj[], int adjncy[], int deg[],
*     int marker[], int *rchsze, int rchset[], int *nhdsze,
*     int nbrhd[]);
*
*  PURPOSE
*
*  This subroutine determines the reachable set of a node through a
*  given subset. The adjancy structure is assumed to be stored in a
*  quotient graph format.
* 
*  INPUT PARAMETERS
*
*  root   - the given node not in the subset;
*  (xadj, adjncy) -
*           the adjancy structure pair;
*  deg    - the degree vector. deg[i] < 0 means the node belongs to the
*           given subset.
*
*  OUTPUT PARAMETERS
*
*  (rchsze, rchset) -
*           the reachable set;
*  (nhdsze, nbrhd) -
*           the neighborhood set.
*
*  UPDATED PARAMETERS
*
*  marker - the marker vector for reach and nbrhd sets. > 0 means the
*           node is in reach set. < 0 means the node has been merged
*           with others in the quotient or it is in nbrhd set.
***********************************************************************/

void qmdrch(int *_root, int xadj[], int adjncy[], int deg[],
      int marker[], int *_rchsze, int rchset[], int *_nhdsze,
      int nbrhd[])
{     int i, istop, istrt, j, jstop, jstrt, nabor, node;
#     define root   (*_root)
#     define rchsze (*_rchsze)
#     define nhdsze (*_nhdsze)
      /* Loop through the neighbors of root in the quotient graph. */
      nhdsze = 0;
      rchsze = 0;
      istrt = xadj[root];
      istop = xadj[root+1] - 1;
      if (istop < istrt) return;
      for (i = istrt; i <= istop; i++)
      {  nabor = adjncy[i];
         if (nabor == 0) return;
         if (marker[nabor] == 0)
         {  if (deg[nabor] >= 0)
            {  /* Include nabor into the reachable set. */
               rchsze++;
               rchset[rchsze] = nabor;
               marker[nabor] = 1;
               goto s600;
            }
            /* nabor has been eliminated. Find nodes reachable from
               it. */
            marker[nabor] = -1;
            nhdsze++;
            nbrhd[nhdsze] = nabor;
s300:       jstrt = xadj[nabor];
            jstop = xadj[nabor+1] - 1;
            for (j = jstrt; j <= jstop; j++)
            {  node = adjncy[j];
               nabor = - node;
               if (node < 0) goto s300;
               if (node == 0) goto s600;
               if (marker[node] == 0)
               {  rchsze++;
                  rchset[rchsze] = node;
                  marker[node] = 1;
               }
            }
         }
s600:    ;
      }
      return;
#     undef root
#     undef rchsze
#     undef nhdsze
}

/***********************************************************************
*  NAME
*
*  qmdqt - Quotient MD Quotient graph Transformation
*
*  SYNOPSIS
*
*  #include "glpqmd.h"
*  void qmdqt(int *root, int xadj[], int adjncy[], int marker[],
*     int *rchsze, int rchset[], int nbrhd[]);
*
*  PURPOSE
*
*  This subroutine performs the quotient graph transformation after a
*  node has been eliminated.
*
*  INPUT PARAMETERS
*
*  root   - the node just eliminated. It becomes the representative of
*           the new supernode;
*  (xadj, adjncy) -
*           the adjancy structure;
*  (rchsze, rchset) -
*           the reachable set of root in the old quotient graph;
*  nbrhd  - the neighborhood set which will be merged with root to form
*           the new supernode;
*  marker - the marker vector.
*
*  UPDATED PARAMETERS
*
*  adjncy - becomes the adjncy of the quotient graph.
***********************************************************************/

void qmdqt(int *_root, int xadj[], int adjncy[], int marker[],
      int *_rchsze, int rchset[], int nbrhd[])
{     int inhd, irch, j, jstop, jstrt, link, nabor, node;
#     define root   (*_root)
#     define rchsze (*_rchsze)
      irch = 0;
      inhd = 0;
      node = root;
s100: jstrt = xadj[node];
      jstop = xadj[node+1] - 2;
      if (jstop >= jstrt)
      {  /* Place reach nodes into the adjacent list of node. */
         for (j = jstrt; j <= jstop; j++)
         {  irch++;
            adjncy[j] = rchset[irch];
            if (irch >= rchsze) goto s400;
         }
      }
      /* Link to other space provided by the nbrhd set. */
      link = adjncy[jstop+1];
      node = - link;
      if (link >= 0)
      {  inhd++;
         node = nbrhd[inhd];
         adjncy[jstop+1] = - node;
      }
      goto s100;
      /* All reachable nodes have been saved. End the adjacent list.
         Add root to the neighborhood list of each node in the reach
         set. */
s400: adjncy[j+1] = 0;
      for (irch = 1; irch <= rchsze; irch++)
      {  node = rchset[irch];
         if (marker[node] >= 0)
         {  jstrt = xadj[node];
            jstop = xadj[node+1] - 1;
            for (j = jstrt; j <= jstop; j++)
            {  nabor = adjncy[j];
               if (marker[nabor] < 0)
               {  adjncy[j] = root;
                  goto s600;
               }
            }
         }
s600:    ;
      }
      return;
#     undef root
#     undef rchsze
}

/***********************************************************************
*  NAME
*
*  qmdupd - Quotient MD UPDate
*
*  SYNOPSIS
*
*  #include "glpqmd.h"
*  void qmdupd(int xadj[], int adjncy[], int *nlist, int list[],
*     int deg[], int qsize[], int qlink[], int marker[], int rchset[],
*     int nbrhd[]);
*
*  PURPOSE
*
*  This routine performs degree update for a set of nodes in the minimum
*  degree algorithm.
*
*  INPUT PARAMETERS
*
*  (xadj, adjncy) -
*           the adjancy structure;
*  (nlist, list) -
*           the list of nodes whose degree has to be updated.
*
*  UPDATED PARAMETERS
*
*  deg    - the degree vector;
*  qsize  - size of indistinguishable supernodes;
*  qlink  - linked list for indistinguishable nodes;
*  marker - used to mark those nodes in reach/nbrhd sets.
*
*  WORKING PARAMETERS
*
*  rchset - the reachable set;
*  nbrhd  - the neighborhood set.
*
*  PROGRAM SUBROUTINES
*
*  qmdmrg.
***********************************************************************/

void qmdupd(int xadj[], int adjncy[], int *_nlist, int list[],
      int deg[], int qsize[], int qlink[], int marker[], int rchset[],
      int nbrhd[])
{     int deg0, deg1, il, inhd, inode, irch, j, jstop, jstrt, mark,
         nabor, nhdsze, node, rchsze;
#     define nlist  (*_nlist)
      /* Find all eliminated supernodes that are adjacent to some nodes
         in the given list. Put them into (nhdsze, nbrhd). deg0 contains
         the number of nodes in the list. */
      if (nlist <= 0) return;
      deg0 = 0;
      nhdsze = 0;
      for (il = 1; il <= nlist; il++)
      {  node = list[il];
         deg0 += qsize[node];
         jstrt = xadj[node];
         jstop = xadj[node+1] - 1;
         for (j = jstrt; j <= jstop; j++)
         {  nabor = adjncy[j];
            if (marker[nabor] == 0 && deg[nabor] < 0)
            {  marker[nabor] = -1;
               nhdsze++;
               nbrhd[nhdsze] = nabor;
            }
         }
      }
      /* Merge indistinguishable nodes in the list by calling the
         subroutine qmdmrg. */
      if (nhdsze > 0)
         qmdmrg(xadj, adjncy, deg, qsize, qlink, marker, &deg0, &nhdsze,
            nbrhd, rchset, &nbrhd[nhdsze+1]);
      /* Find the new degrees of the nodes that have not been merged. */
      for (il = 1; il <= nlist; il++)
      {  node = list[il];
         mark = marker[node];
         if (mark == 0 || mark == 1)
         {  marker[node] = 2;
            qmdrch(&node, xadj, adjncy, deg, marker, &rchsze, rchset,
               &nhdsze, nbrhd);
            deg1 = deg0;
            if (rchsze > 0)
            {  for (irch = 1; irch <= rchsze; irch++)
               {  inode = rchset[irch];
                  deg1 += qsize[inode];
                  marker[inode] = 0;
               }
            }
            deg[node] = deg1 - 1;
            if (nhdsze > 0)
            {  for (inhd = 1; inhd <= nhdsze; inhd++)
               {  inode = nbrhd[inhd];
                  marker[inode] = 0;
               }
            }
         }
      }
      return;
#     undef nlist
}

/***********************************************************************
*  NAME
*
*  qmdmrg - Quotient MD MeRGe
*
*  SYNOPSIS
*
*  #include "qmdmrg.h"
*  void qmdmrg(int xadj[], int adjncy[], int deg[], int qsize[],
*     int qlink[], int marker[], int *deg0, int *nhdsze, int nbrhd[],
*     int rchset[], int ovrlp[]);
*
*  PURPOSE
*
*  This routine merges indistinguishable nodes in the minimum degree
*  ordering algorithm. It also computes the new degrees of these new
*  supernodes.
*
*  INPUT PARAMETERS
*
*  (xadj, adjncy) -
*           the adjancy structure;
*  deg0   - the number of nodes in the given set;
*  (nhdsze, nbrhd) -
*           the set of eliminated supernodes adjacent to some nodes in
*           the set.
*
*  UPDATED PARAMETERS
*
*  deg    - the degree vector;
*  qsize  - size of indistinguishable nodes;
*  qlink  - linked list for indistinguishable nodes;
*  marker - the given set is given by those nodes with marker value set
*           to 1. Those nodes with degree updated will have marker value
*           set to 2.
*
*  WORKING PARAMETERS
*
*  rchset - the reachable set;
*  ovrlp  - temp vector to store the intersection of two reachable sets.
***********************************************************************/

void qmdmrg(int xadj[], int adjncy[], int deg[], int qsize[],
      int qlink[], int marker[], int *_deg0, int *_nhdsze, int nbrhd[],
      int rchset[], int ovrlp[])
{     int deg1, head, inhd, iov, irch, j, jstop, jstrt, link, lnode,
         mark, mrgsze, nabor, node, novrlp, rchsze, root;
#     define deg0   (*_deg0)
#     define nhdsze (*_nhdsze)
      /* Initialization. */
      if (nhdsze <= 0) return;
      for (inhd = 1; inhd <= nhdsze; inhd++)
      {  root = nbrhd[inhd];
         marker[root] = 0;
      }
      /* Loop through each eliminated supernode in the set
         (nhdsze, nbrhd). */
      for (inhd = 1; inhd <= nhdsze; inhd++)
      {  root = nbrhd[inhd];
         marker[root] = -1;
         rchsze = 0;
         novrlp = 0;
         deg1 = 0;
s200:    jstrt = xadj[root];
         jstop = xadj[root+1] - 1;
         /* Determine the reachable set and its intersection with the
            input reachable set. */
         for (j = jstrt; j <= jstop; j++)
         {  nabor = adjncy[j];
            root = - nabor;
            if (nabor < 0) goto s200;
            if (nabor == 0) break;
            mark = marker[nabor];
            if (mark == 0)
            {  rchsze++;
               rchset[rchsze] = nabor;
               deg1 += qsize[nabor];
               marker[nabor] = 1;
            }
            else if (mark == 1)
            {  novrlp++;
               ovrlp[novrlp] = nabor;
               marker[nabor] = 2;
            }
         }
         /* From the overlapped set, determine the nodes that can be
            merged together. */
         head = 0;
         mrgsze = 0;
         for (iov = 1; iov <= novrlp; iov++)
         {  node = ovrlp[iov];
            jstrt = xadj[node];
            jstop = xadj[node+1] - 1;
            for (j = jstrt; j <= jstop; j++)
            {  nabor = adjncy[j];
               if (marker[nabor] == 0)
               {  marker[node] = 1;
                  goto s1100;
               }
            }
            /* Node belongs to the new merged supernode. Update the
               vectors qlink and qsize. */
            mrgsze += qsize[node];
            marker[node] = -1;
            lnode = node;
s900:       link = qlink[lnode];
            if (link > 0)
            {  lnode = link;
               goto s900;
            }
            qlink[lnode] = head;
            head = node;
s1100:      ;
         }
         if (head > 0)
         {  qsize[head] = mrgsze;
            deg[head] = deg0 + deg1 - 1;
            marker[head] = 2;
         }
         /* Reset marker values. */
         root = nbrhd[inhd];
         marker[root] = 0;
         if (rchsze > 0)
         {  for (irch = 1; irch <= rchsze; irch++)
            {  node = rchset[irch];
               marker[node] = 0;
            }
         }
      }
      return;
#     undef deg0
#     undef nhdsze
}

/* eof */
