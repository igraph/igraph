/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2003, 2004  Gabor Csardi <csardi@rmki.kfki.hu>
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

SEXP REST_create_lattice(SEXP dimvector, SEXP pnrnodes, 
			 SEXP pcirc, SEXP pdir) {

  SEXP result;
  
  long int dims;
  long int no_of_nodes;
  int circular;
  int directed;

  vector_t edges;

  long int *coords, *weights;

  long int i, j;

  long int resp=0;

  int carry, pos;

  long int nei;

  dims=GET_LENGTH(dimvector);
  no_of_nodes=R(pnrnodes);
  circular=LOGICAL(pcirc)[0];
  directed=LOGICAL(pdir)[0];

  /* init coords & weights */

  coords=(long int *) R_alloc(dims, sizeof(long int));
  weights=(long int *) R_alloc(dims, sizeof(long int));
  memset(coords, 0, dims*sizeof(long int));
  weights[0]=1;
  for (i=1; i<dims; i++) {
    weights[i]=weights[i-1]*REAL(dimvector)[i-1];
  }
  
  vector_init(&edges, 0);
  vector_reserve(&edges, no_of_nodes*dims + 
		 directed * no_of_nodes*dims);
  
  for (i=0; i<no_of_nodes; i++) {
    for (j=0; j<dims; j++) {
      if (circular || coords[j] != REAL(dimvector)[j]-1) {
	long int new_nei;
	if (coords[j] != REAL(dimvector)[j]-1) {
	  new_nei = i + weights[j] + 1;
	} else {
	  new_nei = i - (REAL(dimvector)[j]-1) * weights[j] + 1;
	}
	if (new_nei != i+1 && (REAL(dimvector)[j] != 2 || coords[j] != 1)) {
	  vector_push_back(&edges, i+1);
	  vector_push_back(&edges, new_nei);
	}
      } /* if circular || coords[j] */
      if (directed && (circular || coords[j] != 0)) {
	long int new_nei;
	if (coords[j]!=0) {
	  new_nei=i-weights[j]+1;
	} else {
	  new_nei=i+(REAL(dimvector)[j]-1) * weights[j]+1;
	}
	if (REAL(dimvector)[j] != 2 || coords[j] != 0) {
	  vector_push_back(&edges, i+1);
	  vector_push_back(&edges, new_nei);
	}
      } /* if circular || coords[0] */
    } /* for j<dims */
    
    /* increase coords */
    carry=1;
    pos=0;
    
    while (carry==1 && pos != dims) {
      if (coords[pos] != REAL(dimvector)[pos]-1) {
	coords[pos]++;
	carry=0;
      } else {
	coords[pos]=0;
	pos++;
      }
    }
    
    R_CheckUserInterrupt();
  } /* for i<no_of_nodes */

  /* Convert vector to SEXP */
  
  i=vector_size(&edges);
  PROTECT(result=NEW_NUMERIC(i));
  for (j=0; j<i; j++) {
    REAL(result)[j] = vector_e(&edges, j);
  }

  /* clean up */
  
  vector_destroy(&edges);
  UNPROTECT(1);
  return result;
}

SEXP REST_connect_neighborhood(SEXP neis, SEXP pradius, SEXP pmutual) {

  SEXP result;
  
  long int no_of_nodes;
  long int radius;
  dqueue_t q;
  long int *already_visited;
  vector_t add;
  long int i,j;
  int mutual;

  no_of_nodes=GET_LENGTH(neis);
  radius=R(pradius);
  mutual=LOGICAL(pmutual)[0];

  already_visited=(long int*) R_alloc(no_of_nodes, sizeof(long int));
  memset(already_visited, 0, no_of_nodes*sizeof(long int));

  dqueue_init(&q, 100);
  vector_init(&add, 0);
  vector_reserve(&add, radius*no_of_nodes);
  
  for (i=1; i<=no_of_nodes; i++) {
    dqueue_push(&q, i);
    dqueue_push(&q, 0);
    already_visited[i-1]=i;
    
    while (!dqueue_empty(&q)) {
      long int actnode=dqueue_pop(&q);
      long int actdist=dqueue_pop(&q);
      if (actdist >= 2 && (actnode > i || mutual)) {
	vector_push_back(&add, i);
	vector_push_back(&add, actnode);
      }
      
      if (actdist+1 <= radius) {
	for (j=0; j<GET_LENGTH(VECTOR_ELT(neis, actnode-1)); j++) {
	  long int neighbor=REAL(VECTOR_ELT(neis, actnode-1))[j];
	  if (already_visited[neighbor-1] == i) { continue; }
	  already_visited[neighbor-1] = i;
	  dqueue_push(&q, neighbor);
	  dqueue_push(&q, actdist+1);
	}
      }
    } /* while !dqueue_empty(q) */
    R_CheckUserInterrupt();
  } /* for i<=no_of_nodes */

  dqueue_destroy(&q);
  
  /* Copy vector to result */
  j=vector_size(&add);
  PROTECT(result=NEW_NUMERIC(j));
  for (i=0; i<j; i++) {
    REAL(result)[i]=vector_e(&add, i);
  }

  /* Clean and return */
  vector_destroy(&add);
  UNPROTECT(1);
  return result;  
}
