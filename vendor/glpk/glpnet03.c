/* glpnet03.c (Klingman's network problem generator) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  This code is the result of translation of the Fortran program NETGEN
*  developed by Dr. Darwin Klingman, which is publically available from
*  NETLIB at <http://www.netlib.org/lp/generators>.
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

#include "glpapi.h"

/***********************************************************************
*  NAME
*
*  glp_netgen - Klingman's network problem generator
*
*  SYNOPSIS
*
*  int glp_netgen(glp_graph *G, int v_rhs, int a_cap, int a_cost,
*     const int parm[1+15]);
*
*  DESCRIPTION
*
*  The routine glp_netgen is a network problem generator developed by
*  Dr. Darwin Klingman. It can create capacitated and uncapacitated
*  minimum cost flow (or transshipment), transportation, and assignment
*  problems.
*
*  The parameter G specifies the graph object, to which the generated
*  problem data have to be stored. Note that on entry the graph object
*  is erased with the routine glp_erase_graph.
*
*  The parameter v_rhs specifies an offset of the field of type double
*  in the vertex data block, to which the routine stores the supply or
*  demand value. If v_rhs < 0, the value is not stored.
*
*  The parameter a_cap specifies an offset of the field of type double
*  in the arc data block, to which the routine stores the arc capacity.
*  If a_cap < 0, the capacity is not stored.
*
*  The parameter a_cost specifies an offset of the field of type double
*  in the arc data block, to which the routine stores the per-unit cost
*  if the arc flow. If a_cost < 0, the cost is not stored.
*
*  The array parm contains description of the network to be generated:
*
*  parm[0]           not used
*  parm[1]  (iseed)  8-digit positive random number seed
*  parm[2]  (nprob)  8-digit problem id number
*  parm[3]  (nodes)  total number of nodes
*  parm[4]  (nsorc)  total number of source nodes (including
*                    transshipment nodes)
*  parm[5]  (nsink)  total number of sink nodes (including
*                    transshipment nodes)
*  parm[6]  (iarcs)  number of arcs
*  parm[7]  (mincst) minimum cost for arcs
*  parm[8]  (maxcst) maximum cost for arcs
*  parm[9]  (itsup)  total supply
*  parm[10] (ntsorc) number of transshipment source nodes
*  parm[11] (ntsink) number of transshipment sink nodes
*  parm[12] (iphic)  percentage of skeleton arcs to be given
*                    the maximum cost
*  parm[13] (ipcap)  percentage of arcs to be capacitated
*  parm[14] (mincap) minimum upper bound for capacitated arcs
*  parm[15] (maxcap) maximum upper bound for capacitated arcs
*
*  The routine generates a transportation problem if:
*
*     nsorc + nsink = nodes, ntsorc = 0, and ntsink = 0.
*
*  The routine generates an assignment problem if the requirements for
*  a transportation problem are met and:
*
*     nsorc = nsink and itsup = nsorc.
*
*  RETURNS
*
*  If the instance was successfully generated, the routine glp_netgen
*  returns zero; otherwise, if specified parameters are inconsistent,
*  the routine returns a non-zero error code.
*
*  REFERENCES
*
*  D.Klingman, A.Napier, and J.Stutz. NETGEN: A program for generating
*  large scale capacitated assignment, transportation, and minimum cost
*  flow networks. Management Science 20 (1974), 814-20. */

struct csa
{     /* common storage area */
      glp_graph *G;
      int v_rhs, a_cap, a_cost;
      int nodes, iarcs, mincst, maxcst, itsup, nsorc, nsink, nonsor,
         nfsink, narcs, nsort, nftsor, ipcap, mincap, maxcap, ktl,
         nodlft, *ipred, *ihead, *itail, *iflag, *isup, *lsinks, mult,
         modul, i15, i16, jran;
};

#define G      (csa->G)
#define v_rhs  (csa->v_rhs)
#define a_cap  (csa->a_cap)
#define a_cost (csa->a_cost)
#define nodes  (csa->nodes)
#define iarcs  (csa->iarcs)
#define mincst (csa->mincst)
#define maxcst (csa->maxcst)
#define itsup  (csa->itsup)
#define nsorc  (csa->nsorc)
#define nsink  (csa->nsink)
#define nonsor (csa->nonsor)
#define nfsink (csa->nfsink)
#define narcs  (csa->narcs)
#define nsort  (csa->nsort)
#define nftsor (csa->nftsor)
#define ipcap  (csa->ipcap)
#define mincap (csa->mincap)
#define maxcap (csa->maxcap)
#define ktl    (csa->ktl)
#define nodlft (csa->nodlft)
#if 0
/* spent a day to find out this bug */
#define ist    (csa->ist)
#else
#define ist    (ipred[0])
#endif
#define ipred  (csa->ipred)
#define ihead  (csa->ihead)
#define itail  (csa->itail)
#define iflag  (csa->iflag)
#define isup   (csa->isup)
#define lsinks (csa->lsinks)
#define mult   (csa->mult)
#define modul  (csa->modul)
#define i15    (csa->i15)
#define i16    (csa->i16)
#define jran   (csa->jran)

static void cresup(struct csa *csa);
static void chain(struct csa *csa, int lpick, int lsorc);
static void chnarc(struct csa *csa, int lsorc);
static void sort(struct csa *csa);
static void pickj(struct csa *csa, int it);
static void assign(struct csa *csa);
static void setran(struct csa *csa, int iseed);
static int iran(struct csa *csa, int ilow, int ihigh);

int glp_netgen(glp_graph *G_, int _v_rhs, int _a_cap, int _a_cost,
      const int parm[1+15])
{     struct csa _csa, *csa = &_csa;
      int iseed, nprob, ntsorc, ntsink, iphic, i, nskel, nltr, ltsink,
         ntrans, npsink, nftr, npsorc, ntravl, ntrrem, lsorc, lpick,
         nsksr, nsrchn, j, item, l, ks, k, ksp, li, n, ii, it, ih, icap,
         jcap, icost, jcost, ret;
      G = G_;
      v_rhs = _v_rhs;
      a_cap = _a_cap;
      a_cost = _a_cost;
      if (G != NULL)
      {  if (v_rhs >= 0 && v_rhs > G->v_size - (int)sizeof(double))
            xerror("glp_netgen: v_rhs = %d; invalid offset\n", v_rhs);
         if (a_cap >= 0 && a_cap > G->a_size - (int)sizeof(double))
            xerror("glp_netgen: a_cap = %d; invalid offset\n", a_cap);
         if (a_cost >= 0 && a_cost > G->a_size - (int)sizeof(double))
            xerror("glp_netgen: a_cost = %d; invalid offset\n", a_cost);
      }
      /* Input the user's random number seed and fix it if
         non-positive. */
      iseed = parm[1];
      nprob = parm[2];
      if (iseed <= 0) iseed = 13502460;
      setran(csa, iseed);
      /* Input the user's problem characteristics. */
      nodes = parm[3];
      nsorc = parm[4];
      nsink = parm[5];
      iarcs = parm[6];
      mincst = parm[7];
      maxcst = parm[8];
      itsup = parm[9];
      ntsorc = parm[10];
      ntsink = parm[11];
      iphic = parm[12];
      ipcap = parm[13];
      mincap = parm[14];
      maxcap = parm[15];
      /* Check the size of the problem. */
      if (!(10 <= nodes && nodes <= 100000))
      {  ret = 1;
         goto done;
      }
      /* Check user supplied parameters for consistency. */
      if (!(nsorc >= 0 && nsink >= 0 && nsorc + nsink <= nodes))
      {  ret = 2;
         goto done;
      }
      if (iarcs < 0)
      {  ret = 3;
         goto done;
      }
      if (mincst > maxcst)
      {  ret = 4;
         goto done;
      }
      if (itsup < 0)
      {  ret = 5;
         goto done;
      }
      if (!(0 <= ntsorc && ntsorc <= nsorc))
      {  ret = 6;
         goto done;
      }
      if (!(0 <= ntsink && ntsink <= nsink))
      {  ret = 7;
         goto done;
      }
      if (!(0 <= iphic && iphic <= 100))
      {  ret = 8;
         goto done;
      }
      if (!(0 <= ipcap && ipcap <= 100))
      {  ret = 9;
         goto done;
      }
      if (mincap > maxcap)
      {  ret = 10;
         goto done;
      }
      /* Initailize the graph object. */
      if (G != NULL)
      {  glp_erase_graph(G, G->v_size, G->a_size);
         glp_add_vertices(G, nodes);
         if (v_rhs >= 0)
         {  double zero = 0.0;
            for (i = 1; i <= nodes; i++)
            {  glp_vertex *v = G->v[i];
               memcpy((char *)v->data + v_rhs, &zero, sizeof(double));
            }
         }
      }
      /* Allocate working arrays. */
      ipred = xcalloc(1+nodes, sizeof(int));
      ihead = xcalloc(1+nodes, sizeof(int));
      itail = xcalloc(1+nodes, sizeof(int));
      iflag = xcalloc(1+nodes, sizeof(int));
      isup = xcalloc(1+nodes, sizeof(int));
      lsinks = xcalloc(1+nodes, sizeof(int));
      /* Print the problem documentation records. */
      if (G == NULL)
      {  xprintf("BEGIN\n");
         xprintf("NETGEN PROBLEM%8d%10s%10d NODES AND%10d ARCS\n",
            nprob, "", nodes, iarcs);
         xprintf("USER:%11d%11d%11d%11d%11d%11d\nDATA:%11d%11d%11d%11d%"
            "11d%11d\n", iseed, nsorc, nsink, mincst,
            maxcst, itsup, ntsorc, ntsink, iphic, ipcap,
            mincap, maxcap);
      }
      else
         glp_set_graph_name(G, "NETGEN");
      /* Set various constants used in the program. */
      narcs = 0;
      nskel = 0;
      nltr = nodes - nsink;
      ltsink = nltr + ntsink;
      ntrans = nltr - nsorc;
      nfsink = nltr + 1;
      nonsor = nodes - nsorc + ntsorc;
      npsink = nsink - ntsink;
      nodlft = nodes - nsink + ntsink;
      nftr = nsorc + 1;
      nftsor = nsorc - ntsorc + 1;
      npsorc = nsorc - ntsorc;
      /* Randomly distribute the supply among the source nodes. */
      if (npsorc + npsink == nodes && npsorc == npsink &&
          itsup == nsorc)
      {  assign(csa);
         nskel = nsorc;
         goto L390;
      }
      cresup(csa);
      /* Print the supply records. */
      if (G == NULL)
      {  xprintf("SUPPLY\n");
         for (i = 1; i <= nsorc; i++)
            xprintf("%6s%6d%18s%10d\n", "", i, "", isup[i]);
         xprintf("ARCS\n");
      }
      else
      {  if (v_rhs >= 0)
         {  for (i = 1; i <= nsorc; i++)
            {  double temp = (double)isup[i];
               glp_vertex *v = G->v[i];
               memcpy((char *)v->data + v_rhs, &temp, sizeof(double));
            }
         }
      }
      /* Make the sources point to themselves in ipred array. */
      for (i = 1; i <= nsorc; i++)
         ipred[i] = i;
      if (ntrans == 0) goto L170;
      /* Chain the transshipment nodes together in the ipred array. */
      ist = nftr;
      ipred[nltr] = 0;
      for (i = nftr; i < nltr; i++)
         ipred[i] = i+1;
      /* Form even length chains for 60 percent of the transshipments.*/
      ntravl = 6 * ntrans / 10;
      ntrrem = ntrans - ntravl;
L140: lsorc = 1;
      while (ntravl != 0)
      {  lpick = iran(csa, 1, ntravl + ntrrem);
         ntravl--;
         chain(csa, lpick, lsorc);
         if (lsorc == nsorc) goto L140;
         lsorc++;
      }
      /* Add the remaining transshipments to the chains. */
      while (ntrrem != 0)
      {
         lpick = iran(csa, 1, ntrrem);
         ntrrem--;
         lsorc = iran(csa, 1, nsorc);
         chain(csa, lpick, lsorc);
      }
L170: /* Set all demands equal to zero. */
      for (i = nfsink; i <= nodes; i++)
         ipred[i] = 0;
      /* The following loop takes one chain at a time (through the use
         of logic contained in the loop and calls to other routines) and
         creates the remaining network arcs. */
      for (lsorc = 1; lsorc <= nsorc; lsorc++)
      {  chnarc(csa, lsorc);
         for (i = nfsink; i <= nodes; i++)
            iflag[i] = 0;
         /* Choose the number of sinks to be hooked up to the current
            chain. */
         if (ntrans != 0)
            nsksr = (nsort * 2 * nsink) / ntrans;
         else
            nsksr = nsink / nsorc + 1;
         if (nsksr < 2) nsksr = 2;
         if (nsksr > nsink) nsksr = nsink;
         nsrchn = nsort;
         /* Randomly pick nsksr sinks and put their names in lsinks. */
         ktl = nsink;
         for (j = 1; j <= nsksr; j++)
         {  item = iran(csa, 1, ktl);
            ktl--;
            for (l = nfsink; l <= nodes; l++)
            {  if (iflag[l] != 1)
               {  item--;
                  if (item == 0) goto L230;
               }
            }
            break;
L230:       lsinks[j] = l;
            iflag[l] = 1;
         }
         /* If last source chain, add all sinks with zero demand to
            lsinks list. */
         if (lsorc == nsorc)
         {  for (j = nfsink; j <= nodes; j++)
            {  if (ipred[j] == 0 && iflag[j] != 1)
               {  nsksr++;
                  lsinks[nsksr] = j;
                  iflag[j] = 1;
               }
            }
         }
         /* Create demands for group of sinks in lsinks. */
         ks = isup[lsorc] / nsksr;
         k = ipred[lsorc];
         for (i = 1; i <= nsksr; i++)
         {  nsort++;
            ksp = iran(csa, 1, ks);
            j = iran(csa, 1, nsksr);
            itail[nsort] = k;
            li = lsinks[i];
            ihead[nsort] = li;
            ipred[li] += ksp;
            li = lsinks[j];
            ipred[li] += ks - ksp;
            n = iran(csa, 1, nsrchn);
            k = lsorc;
            for (ii = 1; ii <= n; ii++)
               k = ipred[k];
         }
         li = lsinks[1];
         ipred[li] += isup[lsorc] - ks * nsksr;
         nskel += nsort;
         /* Sort the arcs in the chain from source lsorc using itail as
            sort key. */
         sort(csa);
         /* Print this part of skeleton and create the arcs for these
            nodes. */
         i = 1;
         itail[nsort+1] = 0;
L300:    for (j = nftsor; j <= nodes; j++)
            iflag[j] = 0;
         ktl = nonsor - 1;
         it = itail[i];
         iflag[it] = 1;
L320:    ih = ihead[i];
         iflag[ih] = 1;
         narcs++;
         ktl--;
         /* Determine if this skeleton arc should be capacitated. */
         icap = itsup;
         jcap = iran(csa, 1, 100);
         if (jcap <= ipcap)
         {  icap = isup[lsorc];
            if (mincap > icap) icap = mincap;
         }
         /* Determine if this skeleton arc should have the maximum
            cost. */
         icost = maxcst;
         jcost = iran(csa, 1, 100);
         if (jcost > iphic)
            icost = iran(csa, mincst, maxcst);
         if (G == NULL)
            xprintf("%6s%6d%6d%2s%10d%10d\n", "", it, ih, "", icost,
               icap);
         else
         {  glp_arc *a = glp_add_arc(G, it, ih);
            if (a_cap >= 0)
            {  double temp = (double)icap;
               memcpy((char *)a->data + a_cap, &temp, sizeof(double));
            }
            if (a_cost >= 0)
            {  double temp = (double)icost;
               memcpy((char *)a->data + a_cost, &temp, sizeof(double));
            }
         }
         i++;
         if (itail[i] == it) goto L320;
         pickj(csa, it);
         if (i <= nsort) goto L300;
      }
      /* Create arcs from the transshipment sinks. */
      if (ntsink != 0)
      {  for (i = nfsink; i <= ltsink; i++)
         {  for (j = nftsor; j <= nodes; j++)
               iflag[j] = 0;
            ktl = nonsor - 1;
            iflag[i] = 1;
            pickj(csa, i);
         }
      }
L390: /* Print the demand records and end record. */
      if (G == NULL)
      {  xprintf("DEMAND\n");
         for (i = nfsink; i <= nodes; i++)
            xprintf("%6s%6d%18s%10d\n", "", i, "", ipred[i]);
         xprintf("END\n");
      }
      else
      {  if (v_rhs >= 0)
         {  for (i = nfsink; i <= nodes; i++)
            {  double temp = - (double)ipred[i];
               glp_vertex *v = G->v[i];
               memcpy((char *)v->data + v_rhs, &temp, sizeof(double));
            }
         }
      }
      /* Free working arrays. */
      xfree(ipred);
      xfree(ihead);
      xfree(itail);
      xfree(iflag);
      xfree(isup);
      xfree(lsinks);
      /* The instance has been successfully generated. */
      ret = 0;
done: return ret;
}

/***********************************************************************
*  The routine cresup randomly distributes the total supply among the
*  source nodes. */

static void cresup(struct csa *csa)
{     int i, j, ks, ksp;
      xassert(itsup > nsorc);
      ks = itsup / nsorc;
      for (i = 1; i <= nsorc; i++)
         isup[i] = 0;
      for (i = 1; i <= nsorc; i++)
      {  ksp = iran(csa, 1, ks);
         j = iran(csa, 1, nsorc);
         isup[i] += ksp;
         isup[j] += ks - ksp;
      }
      j = iran(csa, 1, nsorc);
      isup[j] += itsup - ks * nsorc;
      return;
}

/***********************************************************************
*  The routine chain adds node lpick to the end of the chain with source
*  node lsorc. */

static void chain(struct csa *csa, int lpick, int lsorc)
{     int i, j, k, l, m;
      k = 0;
      m = ist;
      for (i = 1; i <= lpick; i++)
      {  l = k;
         k = m;
         m = ipred[k];
      }
      ipred[l] = m;
      j = ipred[lsorc];
      ipred[k] = j;
      ipred[lsorc] = k;
      return;
}

/***********************************************************************
*  The routine chnarc puts the arcs in the chain from source lsorc into
*  the ihead and itail arrays for sorting. */

static void chnarc(struct csa *csa, int lsorc)
{     int ito, ifrom;
      nsort = 0;
      ito = ipred[lsorc];
L10:  if (ito == lsorc) return;
      nsort++;
      ifrom = ipred[ito];
      ihead[nsort] = ito;
      itail[nsort] = ifrom;
      ito = ifrom;
      goto L10;
}

/***********************************************************************
*  The routine sort sorts the nsort arcs in the ihead and itail arrays.
*  ihead is used as the sort key (i.e. forward star sort order). */

static void sort(struct csa *csa)
{     int i, j, k, l, m, n, it;
      n = nsort;
      m = n;
L10:  m /= 2;
      if (m == 0) return;
      k = n - m;
      j = 1;
L20:  i = j;
L30:  l = i + m;
      if (itail[i] <= itail[l]) goto L40;
      it = itail[i];
      itail[i] = itail[l];
      itail[l] = it;
      it = ihead[i];
      ihead[i] = ihead[l];
      ihead[l] = it;
      i -= m;
      if (i >= 1) goto L30;
L40:  j++;
      if (j <= k) goto L20;
      goto L10;
}

/***********************************************************************
*  The routine pickj creates a random number of arcs out of node 'it'.
*  Various parameters are dynamically adjusted in an attempt to ensure
*  that the generated network has the correct number of arcs. */

static void pickj(struct csa *csa, int it)
{     int j, k, l, nn, nupbnd, icap, jcap, icost;
      if ((nodlft - 1) * 2 > iarcs - narcs - 1)
      {  nodlft--;
         return;
      }
      if ((iarcs - narcs + nonsor - ktl - 1) / nodlft - nonsor + 1 >= 0)
         k = nonsor;
      else
      {  nupbnd = (iarcs - narcs - nodlft) / nodlft * 2;
L40:     k = iran(csa, 1, nupbnd);
         if (nodlft == 1) k = iarcs - narcs;
         if ((nodlft - 1) * (nonsor - 1) < iarcs - narcs - k) goto L40;
      }
      nodlft--;
      for (j = 1; j <= k; j++)
      {  nn = iran(csa, 1, ktl);
         ktl--;
         for (l = nftsor; l <= nodes; l++)
         {  if (iflag[l] != 1)
            {  nn--;
               if (nn == 0) goto L70;
            }
         }
         return;
L70:     iflag[l] = 1;
         icap = itsup;
         jcap = iran(csa, 1, 100);
         if (jcap <= ipcap)
            icap = iran(csa, mincap, maxcap);
         icost = iran(csa, mincst, maxcst);
         if (G == NULL)
            xprintf("%6s%6d%6d%2s%10d%10d\n", "", it, l, "", icost,
               icap);
         else
         {  glp_arc *a = glp_add_arc(G, it, l);
            if (a_cap >= 0)
            {  double temp = (double)icap;
               memcpy((char *)a->data + a_cap, &temp, sizeof(double));
            }
            if (a_cost >= 0)
            {  double temp = (double)icost;
               memcpy((char *)a->data + a_cost, &temp, sizeof(double));
            }
         }
         narcs++;
      }
      return;
}

/***********************************************************************
*  The routine assign generate assignment problems. It defines the unit
*  supplies, builds a skeleton, then calls pickj to create the arcs. */

static void assign(struct csa *csa)
{     int i, it, nn, l, ll, icost;
      if (G == NULL)
         xprintf("SUPPLY\n");
      for (i = 1; i <= nsorc; i++)
      {  isup[i] = 1;
         iflag[i] = 0;
         if (G == NULL)
            xprintf("%6s%6d%18s%10d\n", "", i, "", isup[i]);
         else
         {  if (v_rhs >= 0)
            {  double temp = (double)isup[i];
               glp_vertex *v = G->v[i];
               memcpy((char *)v->data + v_rhs, &temp, sizeof(double));
            }
         }
      }
      if (G == NULL)
         xprintf("ARCS\n");
      for (i = nfsink; i <= nodes; i++)
         ipred[i] = 1;
      for (it = 1; it <= nsorc; it++)
      {  for (i = nfsink; i <= nodes; i++)
            iflag[i] = 0;
         ktl = nsink - 1;
         nn = iran(csa, 1, nsink - it + 1);
         for (l = 1; l <= nsorc; l++)
         {  if (iflag[l] != 1)
            {  nn--;
               if (nn == 0) break;
            }
         }
         narcs++;
         ll = nsorc + l;
         icost = iran(csa, mincst, maxcst);
         if (G == NULL)
            xprintf("%6s%6d%6d%2s%10d%10d\n", "", it, ll, "", icost,
               isup[1]);
         else
         {  glp_arc *a = glp_add_arc(G, it, ll);
            if (a_cap >= 0)
            {  double temp = (double)isup[1];
               memcpy((char *)a->data + a_cap, &temp, sizeof(double));
            }
            if (a_cost >= 0)
            {  double temp = (double)icost;
               memcpy((char *)a->data + a_cost, &temp, sizeof(double));
            }
         }
         iflag[l] = 1;
         iflag[ll] = 1;
         pickj(csa, it);
      }
      return;
}

/***********************************************************************
*  Portable congruential (uniform) random number generator:
*
*     next_value = ((7**5) * previous_value) modulo ((2**31)-1)
*
*  This generator consists of three routines:
*
*  (1) setran - initializes constants and seed
*  (2) iran   - generates an integer random number
*  (3) rran   - generates a real random number
*
*  The generator requires a machine with at least 32 bits of precision.
*  The seed (iseed) must be in the range [1,(2**31)-1]. */

static void setran(struct csa *csa, int iseed)
{     xassert(iseed >= 1);
      mult = 16807;
      modul = 2147483647;
      i15 = 1 << 15;
      i16 = 1 << 16;
      jran = iseed;
      return;
}

/***********************************************************************
*  The routine iran generates an integer random number between ilow and
*  ihigh. If ilow > ihigh then iran returns ihigh. */

static int iran(struct csa *csa, int ilow, int ihigh)
{     int ixhi, ixlo, ixalo, leftlo, ixahi, ifulhi, irtlo, iover,
         irthi, j;
      ixhi = jran / i16;
      ixlo = jran - ixhi * i16;
      ixalo = ixlo * mult;
      leftlo = ixalo / i16;
      ixahi = ixhi * mult;
      ifulhi = ixahi + leftlo;
      irtlo = ixalo - leftlo * i16;
      iover = ifulhi / i15;
      irthi = ifulhi - iover * i15;
      jran = ((irtlo - modul) + irthi * i16) + iover;
      if (jran < 0) jran += modul;
      j = ihigh - ilow + 1;
      if (j > 0)
         return jran % j + ilow;
      else
         return ihigh;
}

/**********************************************************************/

#if 0
static int scan(char card[80+1], int pos, int len)
{     char buf[10+1];
      memcpy(buf, &card[pos-1], len);
      buf[len] = '\0';
      return atoi(buf);
}

int main(void)
{     int parm[1+15];
      char card[80+1];
      xassert(fgets(card, sizeof(card), stdin) == card);
      parm[1] = scan(card, 1, 8);
      parm[2] = scan(card, 9, 8);
      xassert(fgets(card, sizeof(card), stdin) == card);
      parm[3] = scan(card, 1, 5);
      parm[4] = scan(card, 6, 5);
      parm[5] = scan(card, 11, 5);
      parm[6] = scan(card, 16, 5);
      parm[7] = scan(card, 21, 5);
      parm[8] = scan(card, 26, 5);
      parm[9] = scan(card, 31, 10);
      parm[10] = scan(card, 41, 5);
      parm[11] = scan(card, 46, 5);
      parm[12] = scan(card, 51, 5);
      parm[13] = scan(card, 56, 5);
      parm[14] = scan(card, 61, 10);
      parm[15] = scan(card, 71, 10);
      glp_netgen(NULL, 0, 0, 0, parm);
      return 0;
}
#endif

/* eof */
