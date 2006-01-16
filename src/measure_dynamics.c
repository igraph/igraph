/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "igraph.h"
#include "memory.h"

#include <math.h>

int igraph_measure_dynamics_idage(const igraph_t *graph, igraph_matrix_t *akl,
				  igraph_matrix_t *sd,
				  const igraph_vector_t *st, integer_t pagebins,
				  integer_t pmaxind, bool_t lsd) {

  long int agebins=pagebins;
  long int maxind=pmaxind;
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth;
  
  int *indegree;
  igraph_matrix_t ntkl, ch, normfact, notnull;
  igraph_vector_t neis;
  
  long int node;
  long int i, j, k;
  long int edges=0;

  binwidth = no_of_nodes/agebins+1;

  igraph_vector_init(&neis, 0);
  indegree=Calloc(no_of_nodes, int);
  igraph_matrix_resize(akl, maxind+1, agebins);
  igraph_matrix_null(akl);
  if (lsd) {
    igraph_matrix_resize(sd, maxind+1, agebins);
    igraph_matrix_null(sd);
  }
  igraph_matrix_init(&ntkl, maxind+1, agebins+1);
  igraph_matrix_init(&ch, maxind+1, agebins+1);
  igraph_matrix_init(&normfact, maxind+1, agebins);
  igraph_matrix_init(&notnull, maxind+1, agebins);

  for (node=0; node<no_of_nodes; node++) {
    
    /* inspect the edges */
   
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      double xk=VECTOR(*st)[node]/MATRIX(ntkl, xidx, yidx);
      double oldm=MATRIX(*akl, xidx, yidx);
      MATRIX(notnull, xidx, yidx) += 1;
      MATRIX(*akl, xidx, yidx) += (xk-oldm)/MATRIX(notnull, xidx, yidx);
      if (lsd) {
	MATRIX(*sd, xidx, yidx) += (xk-oldm)*(xk-MATRIX(*akl, xidx, yidx));
      }

      indegree[to] ++;
      MATRIX(ntkl, xidx, yidx)--;
      if (MATRIX(ntkl,xidx, yidx)==0) {
	MATRIX(normfact, xidx, yidx) += (edges-MATRIX(ch, xidx, yidx)+1);
	MATRIX(ch, xidx, yidx)=edges;
      }
      MATRIX(ntkl, xidx+1, yidx)++;
      if (MATRIX(ntkl, xidx+1, yidx)==1) {
	MATRIX(ch, xidx+1, yidx)=edges;
      }
      edges++;
    }

    /* new node, aging */
    MATRIX(ntkl, 0, 0)++;
    if (MATRIX(ntkl, 0, 0)==1) {
      MATRIX(ch, 0, 0)=edges;
    }
    for (k=1; node-binwidth*k+1 >=1; k++) {
      long int shnode=node-binwidth*k;
      long int deg=indegree[shnode];
      MATRIX(ntkl, deg, k-1)--;
      if (MATRIX(ntkl, deg, k-1)==0) {
	MATRIX(normfact, deg, k-1) += (edges-MATRIX(ch, deg, k-1)+1);
	MATRIX(ch, deg, k-1)=edges;
      }
      MATRIX(ntkl, deg, k)++;
      if (MATRIX(ntkl, deg, k)==1) {
	MATRIX(ch, deg, k)=edges;
      }
    }
  }

  /* Ok, measurement done, update change */
  for (i=0; i<maxind+1; i++) {
    for (j=0; j<agebins; j++) {
      if (MATRIX(ntkl, i, j) != 0) {
	MATRIX(normfact, i, j) += (edges-MATRIX(ch, i, j)+1);
      }
      MATRIX(*akl, i, j) *= MATRIX(notnull, i, j) / MATRIX(normfact, i, j);
      if (lsd) {
	MATRIX(*sd, i, j) +=
	  MATRIX(*akl, i, j)*MATRIX(*akl,i,j)*
	  (MATRIX(normfact,i,j)/MATRIX(notnull,i,j)-1);
	if (MATRIX(normfact,i,j) > 0) {
	  MATRIX(*sd, i, j) =
	    sqrt(MATRIX(*sd, i, j)/(MATRIX(normfact,i,j)-1));
	  MATRIX(*sd, i, j) =
	    2 * MATRIX(*sd,i,j)/sqrt(MATRIX(normfact,i,j));
	}
      }
    }
  }
  
  Free(indegree);
  igraph_matrix_destroy(&ntkl);
  igraph_matrix_destroy(&ch);
  igraph_matrix_destroy(&normfact);
  igraph_matrix_destroy(&notnull);
  igraph_vector_destroy(&neis);

  return 0;
}

int igraph_measure_dynamics_idage_st(const igraph_t *graph, igraph_vector_t *res,
				     const igraph_matrix_t *akl) {

  long int agebins=igraph_matrix_ncol(akl);
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth;
  
  int *indegree;
  igraph_vector_t neis;
  
  long int node;
  long int i, k;

  igraph_vector_init(&neis, 0);
  
  indegree=Calloc(no_of_nodes, int);
  binwidth=no_of_nodes/agebins+1;
  
  igraph_vector_resize(res, no_of_nodes);
  igraph_vector_null(res);
  VECTOR(*res)[0]=MATRIX(*akl, 0, 0);

  for (node=1; node<no_of_nodes; node++) {
    
    /* new node, aging */
    VECTOR(*res)[node] = VECTOR(*res)[node-1] + MATRIX(*akl, 0, 0);
    for (k=1; node-binwidth*k+1 >= 1; k++) {
      long int shnode=node-binwidth*k;
      long int deg=indegree[shnode];
      VECTOR(*res)[node] += -MATRIX(*akl, deg, k-1)+MATRIX(*akl, deg, k);
    }
    
    /* inspect the outgoing edges */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      indegree[to] ++;
      
      VECTOR(*res)[node] +=
	-MATRIX(*akl, xidx, yidx) + MATRIX(*akl, xidx+1, yidx);
    }
  }
  
  igraph_vector_destroy(&neis);
  Free(indegree);
  
  return 0;
}

int igraph_measure_dynamics_idage_debug(const igraph_t *graph, igraph_matrix_t *akl,
					igraph_matrix_t *sd,
					const igraph_vector_t *st, integer_t pagebins,
					integer_t pmaxind, bool_t lsd,
					igraph_vector_t *estimates, 
					integer_t est_ind, integer_t est_age) {
  
  long int agebins=pagebins;
  long int maxind=pmaxind;
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth;
  
  int *indegree;
  igraph_matrix_t ntkl, ch, normfact, notnull;
  igraph_vector_t neis;
  
  long int node;
  long int i, j, k;
  long int edges=0;

  binwidth = no_of_nodes/agebins+1;

  igraph_vector_init(&neis, 0);
  indegree=Calloc(no_of_nodes, int);
  igraph_matrix_resize(akl, maxind+1, agebins);
  igraph_matrix_null(akl);
  if (lsd) {
    igraph_matrix_resize(sd, maxind+1, agebins);
    igraph_matrix_null(sd);
  }

  igraph_vector_clear(estimates);

  igraph_matrix_init(&ntkl, maxind+1, agebins+1);
  igraph_matrix_init(&ch, maxind+1, agebins+1);
  igraph_matrix_init(&normfact, maxind+1, agebins);
  igraph_matrix_init(&notnull, maxind+1, agebins);

  for (node=0; node<no_of_nodes; node++) {
    
    /* inspect the edges */
   
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      double xk=VECTOR(*st)[node]/MATRIX(ntkl, xidx, yidx);
      double oldm=MATRIX(*akl, xidx, yidx);
      MATRIX(notnull, xidx, yidx) += 1;
      MATRIX(*akl, xidx, yidx) += (xk-oldm)/MATRIX(notnull, xidx, yidx);
      if (lsd) {
	MATRIX(*sd, xidx, yidx) += (xk-oldm)*(xk-MATRIX(*akl, xidx, yidx));
      }

      if (xidx==est_ind && yidx==est_age) {
	igraph_vector_push_back(estimates, xk);
      }

      indegree[to] ++;
      MATRIX(ntkl, xidx, yidx)--;
      if (MATRIX(ntkl,xidx, yidx)==0) {
	MATRIX(normfact, xidx, yidx) += (edges-MATRIX(ch, xidx, yidx)+1);
	MATRIX(ch, xidx, yidx)=edges;
      }
      MATRIX(ntkl, xidx+1, yidx)++;
      if (MATRIX(ntkl, xidx+1, yidx)==1) {
	MATRIX(ch, xidx+1, yidx)=edges;
      }
      edges++;
    }

    /* new node, aging */
    MATRIX(ntkl, 0, 0)++;
    if (MATRIX(ntkl, 0, 0)==1) {
      MATRIX(ch, 0, 0)=edges;
    }
    for (k=1; node-binwidth*k+1 >=1; k++) {
      long int shnode=node-binwidth*k;
      long int deg=indegree[shnode];
      MATRIX(ntkl, deg, k-1)--;
      if (MATRIX(ntkl, deg, k-1)==0) {
	MATRIX(normfact, deg, k-1) += (edges-MATRIX(ch, deg, k-1)+1);
	MATRIX(ch, deg, k-1)=edges;
      }
      MATRIX(ntkl, deg, k)++;
      if (MATRIX(ntkl, deg, k)==1) {
	MATRIX(ch, deg, k)=edges;
      }
    }
  }

  /* Ok, measurement done, update change */
  for (i=0; i<maxind+1; i++) {
    for (j=0; j<agebins; j++) {
      if (MATRIX(ntkl, i, j) != 0) {
	MATRIX(normfact, i, j) += (edges-MATRIX(ch, i, j)+1);
      }
      MATRIX(*akl, i, j) *= MATRIX(notnull, i, j) / MATRIX(normfact, i, j);
      if (lsd) {
	MATRIX(*sd, i, j) +=
	  MATRIX(*akl, i, j)*MATRIX(*akl,i,j)*
	  (MATRIX(normfact,i,j)/MATRIX(notnull,i,j)-1);
	if (MATRIX(normfact,i,j) > 0) {
	  MATRIX(*sd, i, j) =
	    sqrt(MATRIX(*sd, i, j)/(MATRIX(normfact,i,j)-1));
	  MATRIX(*sd, i, j) =
	    2 * MATRIX(*sd,i,j)/sqrt(MATRIX(normfact,i,j));
	}
      }
    }
  }

  for (i=0; i<igraph_vector_size(estimates); i++) {
    VECTOR(*estimates)[i] /= MATRIX(normfact, (long int)est_ind, 
				    (long int)est_age);
  }
/*   igraph_vector_push_back(estimates, MATRIX(normfact, est_ind, est_age)); */
  
  Free(indegree);
  igraph_matrix_destroy(&ntkl);
  igraph_matrix_destroy(&ch);
  igraph_matrix_destroy(&normfact);
  igraph_matrix_destroy(&notnull);
  igraph_vector_destroy(&neis);

  return 0;
}
