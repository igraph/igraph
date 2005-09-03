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

int igraph_measure_dynamics_idage(igraph_t *graph, matrix_t *akl,
				  matrix_t *sd,
				  vector_t *st, integer_t pagebins,
				  integer_t pmaxind, bool_t lsd) {

  long int agebins=pagebins;
  long int maxind=pmaxind;
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth;
  
  int *indegree;
  matrix_t ntkl, ch, normfact, notnull;
  vector_t neis;
  
  long int node;
  long int i, j, k;
  long int edges=0;

/*   int mptr=1; int m_x=3, m_y=10; */
  
  binwidth = no_of_nodes/agebins+1;

  vector_init(&neis, 0);
  indegree=Calloc(no_of_nodes, int);
  matrix_resize(akl, maxind+1, agebins);
  matrix_null(akl);
  if (lsd) {
    matrix_resize(sd, maxind+1, agebins);
    matrix_null(sd);
  }
/*   SET_VECTOR_ELT(result, 2, NEW_NUMERIC(10000)); */
/*   memset(REAL(VECTOR_ELT(result,2)), 0, sizeof(double)*10000); */
  matrix_init(&ntkl, maxind+1, agebins+1);
  matrix_init(&ch, maxind+1, agebins+1);
  matrix_init(&normfact, maxind+1, agebins);
  matrix_init(&notnull, maxind+1, agebins);

  for (node=0; node<no_of_nodes; node++) {
    
    /* inspect the edges */
   
    igraph_neighbors(graph, &neis, node, 1);
    for (i=0; i<vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
/*       if (xidx+1 >= maxind+1 || yidx >=agebins) {  */
/* 	printf("GEBAASZ: %li %li %li %li %li %li\n", node, to, xidx+1, maxind+1, yidx, agebins);  */
/*       } */
      
      double xk=VECTOR(*st)[node]/MATRIX(ntkl, xidx, yidx);
      double oldm=MATRIX(*akl, xidx, yidx);
      MATRIX(notnull, xidx, yidx) += 1;
      MATRIX(*akl, xidx, yidx) += (xk-oldm)/MATRIX(notnull, xidx, yidx);
      if (lsd) {
	MATRIX(*sd, xidx, yidx) += (xk-oldm)*(xk-MATRIX(*akl, xidx, yidx));
      }

/*       if (xidx==m_x && yidx==m_y) { */
/* 	REAL(VECTOR_ELT(result, 2))[mptr++] = xk; */
/*       } */

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
/*       if (deg >= maxind+1 || k >=agebins) {  */
/* 	printf("GEBAASZ 2: %li %li %li %li %li\n", node, */
/* 	       deg, maxind+1, k, agebins); */
/*       } */
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

/*   REAL(VECTOR_ELT(result,2))[0]=RMATRIX(normfact, m_x, m_y); */
/*   for (i=1; i<mptr; i++) { */
/*     REAL(VECTOR_ELT(result, 2))[i] /= RMATRIX(normfact, m_x, m_y); */
/*   } */
  
  Free(indegree);
  matrix_destroy(&ntkl);
  matrix_destroy(&ch);
  matrix_destroy(&normfact);
  matrix_destroy(&notnull);
  vector_destroy(&neis);

  return 0;
}

int igraph_measure_dynamics_idage_st(igraph_t *graph, vector_t *res,
				     matrix_t *akl) {

  long int maxind=matrix_nrow(akl);
  long int agebins=matrix_ncol(akl);
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth;
  
  int *indegree;
  vector_t neis;
  
  long int node;
  long int i, j, k;

  vector_init(&neis, 0);
  
  indegree=Calloc(no_of_nodes, int);
  binwidth=no_of_nodes/agebins+1;
  
  vector_resize(res, no_of_nodes);
  vector_null(res);
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
    igraph_neighbors(graph, &neis, node, 1);
    for (i=0; i<vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      indegree[to] ++;
      
      VECTOR(*res)[node] +=
	-MATRIX(*akl, xidx, yidx) + MATRIX(*akl, xidx+1, yidx);
    }
  }
  
  vector_destroy(&neis);
  Free(indegree);
  
  return 0;
}
