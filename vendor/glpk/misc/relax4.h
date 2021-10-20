/* relax4.h */

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

#ifndef RELAX4_H
#define RELAX4_H

struct relax4_csa
{     /* common storage area */
      /* input parameters --------------------------------------------*/
      int n;
      /* number of nodes */
      int na;
      /* number of arcs */
      int large;
      /* very large int to represent infinity */
      int repeat;
      /* true if initialization is to be skipped (false otherwise) */
      int crash;
      /* 0 if default initialization is used
       * 1 if auction initialization is used */
      int *startn; /* int startn[1+na]; */
      /* startn[j] = starting node for arc j, j = 1,...,na */
      int *endn; /* int endn[1+na] */
      /* endn[j] = ending node for arc j, j = 1,...,na */
      int *fou; /* int fou[1+n]; */
      /* fou[i] = first arc out of node i, i = 1,...,n */
      int *nxtou; /* int nxtou[1+na]; */
      /* nxtou[j] = next arc out of the starting node of arc j,
       * j = 1,...,na */
      int *fin; /* int fin[1+n]; */
      /* fin[i] = first arc into node i, i = 1,...,n */
      int *nxtin; /* int nxtin[1+na]; */
      /* nxtin[j] = next arc into the ending node of arc j,
       * j = 1,...,na */
      /* updated parameters ------------------------------------------*/
      int *rc; /* int rc[1+na]; */
      /* rc[j] = reduced cost of arc j, j = 1,...,na */
      int *u; /* int u[1+na]; */
      /* u[j] = capacity of arc j on input
       * and (capacity of arc j) - x(j) on output, j = 1,...,na */
      int *dfct; /* int dfct[1+n]; */
      /* dfct[i] = demand at node i on input
       * and zero on output, i = 1,...,n */
      /* output parameters -------------------------------------------*/
      int *x; /* int x[1+na]; */
      /* x[j] = flow on arc j, j = 1,...,na */
      int nmultinode;
      /* number of multinode relaxation iterations in RELAX4 */
      int iter;
      /* number of relaxation iterations in RELAX4 */
      int num_augm;
      /* number of flow augmentation steps in RELAX4 */
      int num_ascnt;
      /* number of multinode ascent steps in RELAX4 */
      int nsp;
      /* number of auction/shortest path iterations */
      /* working parameters ------------------------------------------*/
      int *label; /* int label, tempin, p[1+n]; */
      int *prdcsr; /* int prdcsr, tempou, price[1+n]; */
      int *save; /* int save[1+na]; */
      int *tfstou; /* int tfstou, fpushf[1+n]; */
      int *tnxtou; /* int tnxtou, nxtpushf[1+na]; */
      int *tfstin; /* int tfstin, fpushb[1+n]; */
      int *tnxtin; /* int tnxtin, nxtpushb[1+na]; */
      int *nxtqueue; /* int nxtqueue[1+n]; */
      char *scan; /* bool scan[1+n]; */
      char *mark; /* bool mark, path_id[1+n]; */
      /* working parameters used by routine auction only -------------*/
      int *extend_arc; /* int extend_arc[1+n]; */
      int *sb_level; /* int sb_level[1+n]; */
      int *sb_arc; /* int sb_arc[1+n]; */
};

#define relax4 _glp_relax4
int relax4(struct relax4_csa *csa);

#define relax4_inidat _glp_relax4_inidat
void relax4_inidat(struct relax4_csa *csa);

#endif

/* eof */
