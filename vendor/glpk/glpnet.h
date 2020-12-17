/* glpnet.h (graph and network algorithms) */

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

#ifndef GLPNET_H
#define GLPNET_H

#define mc21a _glp_mc21a
int mc21a(int n, const int icn[], const int ip[], const int lenr[],
      int iperm[], int pr[], int arp[], int cv[], int out[]);
/* permutations for zero-free diagonal */

#define mc13d _glp_mc13d
int mc13d(int n, const int icn[], const int ip[], const int lenr[],
      int ior[], int ib[], int lowl[], int numb[], int prev[]);
/* permutations to block triangular form */

#define okalg _glp_okalg
int okalg(int nv, int na, const int tail[], const int head[],
      const int low[], const int cap[], const int cost[], int x[],
      int pi[]);
/* out-of-kilter algorithm */

#define ffalg _glp_ffalg
void ffalg(int nv, int na, const int tail[], const int head[],
      int s, int t, const int cap[], int x[], char cut[]);
/* Ford-Fulkerson algorithm */

#define wclique _glp_wclique
int wclique(int n, const int w[], const unsigned char a[], int ind[]);
/* find maximum weight clique with Ostergard's algorithm */

#define kellerman _glp_kellerman
int kellerman(int n, int (*func)(void *info, int i, int ind[]),
      void *info, void /* glp_graph */ *H);
/* cover edges by cliques with Kellerman's heuristic */

#endif

/* eof */
