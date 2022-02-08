/* cfg.h (conflict graph) */

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

#ifndef CFG_H
#define CFG_H

#include "dmp.h"

/***********************************************************************
*  The structure CFG describes the conflict graph.
*
*  Conflict graph is an undirected graph G = (V, E), where V is a set
*  of vertices, E <= V x V is a set of edges. Each vertex v in V of the
*  conflict graph corresponds to a binary variable z[v], which is
*  either an original binary variable x[j] or its complement 1 - x[j].
*  Edge (v,w) in E means that z[v] and z[w] cannot take the value 1 at
*  the same time, i.e. it defines an inequality z[v] + z[w] <= 1, which
*  is assumed to be valid for original MIP.
*
*  Since the conflict graph may be dense, it is stored as an union of
*  its cliques rather than explicitly. */

#if 0 /* 08/III-2016 */
typedef struct CFG CFG;
#else
typedef struct glp_cfg CFG;
#endif
typedef struct CFGVLE CFGVLE;
typedef struct CFGCLE CFGCLE;

#if 0 /* 08/III-2016 */
struct CFG
#else
struct glp_cfg
#endif
{     /* conflict graph descriptor */
      int n;
      /* number of *all* variables (columns) in corresponding MIP */
      int *pos; /* int pos[1+n]; */
      /* pos[0] is not used;
       * pos[j] = v, 1 <= j <= n, means that vertex v corresponds to
       * original binary variable x[j], and pos[j] = 0 means that the
       * conflict graph has no such vertex */
      int *neg; /* int neg[1+n]; */
      /* neg[0] is not used;
       * neg[j] = v, 1 <= j <= n, means that vertex v corresponds to
       * complement of original binary variable x[j], and neg[j] = 0
       * means that the conflict graph has no such vertex */
      DMP *pool;
      /* memory pool to allocate elements of the conflict graph */
      int nv_max;
      /* maximal number of vertices in the conflict graph */
      int nv;
      /* current number of vertices in the conflict graph */
      int *ref; /* int ref[1+nv_max]; */
      /* ref[v] = j, 1 <= v <= nv, means that vertex v corresponds
       * either to original binary variable x[j] or to its complement,
       * i.e. either pos[j] = v or neg[j] = v */
      CFGVLE **vptr; /* CFGVLE *vptr[1+nv_max]; */
      /* vptr[v], 1 <= v <= nv, is an initial pointer to the list of
       * vertices adjacent to vertex v */
      CFGCLE **cptr; /* CFGCLE *cptr[1+nv_max]; */
      /* cptr[v], 1 <= v <= nv, is an initial pointer to the list of
       * cliques that contain vertex v */
};

struct CFGVLE
{     /* vertex list element */
      int v;
      /* vertex number, 1 <= v <= nv */
      CFGVLE *next;
      /* pointer to next vertex list element */
};

struct CFGCLE
{     /* clique list element */
      CFGVLE *vptr;
      /* initial pointer to the list of clique vertices */
      CFGCLE *next;
      /* pointer to next clique list element */
};

#define cfg_create_graph _glp_cfg_create_graph
CFG *cfg_create_graph(int n, int nv_max);
/* create conflict graph */

#define cfg_add_clique _glp_cfg_add_clique
void cfg_add_clique(CFG *G, int size, const int ind[]);
/* add clique to conflict graph */

#define cfg_get_adjacent _glp_cfg_get_adjacent
int cfg_get_adjacent(CFG *G, int v, int ind[]);
/* get vertices adjacent to specified vertex */

#define cfg_expand_clique _glp_cfg_expand_clique
int cfg_expand_clique(CFG *G, int c_len, int c_ind[]);
/* expand specified clique to maximal clique */

#define cfg_check_clique _glp_cfg_check_clique
void cfg_check_clique(CFG *G, int c_len, const int c_ind[]);
/* check clique in conflict graph */

#define cfg_delete_graph _glp_cfg_delete_graph
void cfg_delete_graph(CFG *G);
/* delete conflict graph */

#define cfg_build_graph _glp_cfg_build_graph
CFG *cfg_build_graph(void /* glp_prob */ *P);
/* build conflict graph */

#define cfg_find_clique _glp_cfg_find_clique
int cfg_find_clique(void /* glp_prob */ *P, CFG *G, int ind[],
      double *sum);
/* find maximum weight clique in conflict graph */

#endif

/* eof */
