/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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

SEXP REST_betweenness (SEXP interface, SEXP graph, SEXP pdirected) {

  REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph);

  SEXP result;
  
  long int no_of_nodes;
  dqueue_t q;
  long int *distance;
  long int *nrgeo;
  double *tmpscore;
  stack_t stack;
  long int source;
  long int j;

  SEXP tmp, modein, modeout;

  no_of_nodes=R(ptrtable.vcount(interface, graph));

  if (LOGICAL(pdirected)[0]) {
    PROTECT(modeout=ScalarString(CREATE_STRING_VECTOR("out")));
    PROTECT(modein =ScalarString(CREATE_STRING_VECTOR("in")));
  } else {
    PROTECT(modeout=ScalarString(CREATE_STRING_VECTOR("all")));
    PROTECT(modein =ScalarString(CREATE_STRING_VECTOR("all")));
  }        
  
  distance=(long int*) R_alloc(no_of_nodes, sizeof(long int));
  nrgeo=(long int*) R_alloc(no_of_nodes, sizeof(long int));
  tmpscore=(double*) R_alloc(no_of_nodes, sizeof(double));
  /* memsetting later */

  PROTECT(result=NEW_NUMERIC(no_of_nodes));
  memset(REAL(result), 0, no_of_nodes*sizeof(double));

  dqueue_init(&q, 100);
  stack_init(&stack, no_of_nodes);

  /* here we go */
  
  for (source=1; source<=no_of_nodes; source++) {
    
    memset(distance, 0, no_of_nodes*sizeof(long int));
    memset(nrgeo, 0, no_of_nodes*sizeof(long int));
    memset(tmpscore, 0, no_of_nodes*sizeof(double));
    stack_clear(&stack); /* it should be empty anyway... */
    
    dqueue_push(&q, source);
    nrgeo[source-1]=1;
    distance[source-1]=0;
    
    while (!dqueue_empty(&q)) {
      long int actnode=dqueue_pop(&q);
      
      tmp=ptrtable.neighbors(interface, graph, actnode, modeout);
      for (j=0; j<GET_LENGTH(tmp); j++) {
	long int neighbor=REAL(tmp)[j];
	if (nrgeo[neighbor-1] != 0) {
	  /* we've already seen this node, another shortest path? */
	  if (distance[neighbor-1]==distance[actnode-1]+1) {
	    nrgeo[neighbor-1]++;
	  }
	} else {
	  /* we haven't seen this node yet */
	  nrgeo[neighbor-1]++;
	  distance[neighbor-1]=distance[actnode-1]+1;
	  dqueue_push(&q, neighbor);
	  stack_push(&stack, neighbor);
	}
      }
    } /* while !dqueue_empty */
    
    /* Ok, we've the distance of each node and also the number of 
       shortest paths to them. Now we do an inverse search, starting
       with the farthest nodes. */
    while (!stack_empty(&stack)) {
      long int actnode=stack_pop(&stack);
      long int friends=0;
      if (distance[actnode-1]<=1) { continue; } /* skip source node */
      
      /* search for the neighbors on the geodesics */
      tmp=ptrtable.neighbors(interface, graph, actnode, modein);
      for (j=0; j<GET_LENGTH(tmp); j++) {
	long int neighbor=REAL(tmp)[j];
	if (distance[neighbor-1]==distance[actnode-1]-1) { friends++; }
      }

      /* set the temporary score of the friends */
      for (j=0; j<GET_LENGTH(tmp); j++) {
	long int neighbor=REAL(tmp)[j];
	if (distance[neighbor-1]==distance[actnode-1]-1) { 
	  tmpscore[neighbor-1] += (tmpscore[actnode-1]+1)/friends;
	}
      }
    }
    
    /* Ok, we've the scores for this source */
    for (j=0; j<no_of_nodes; j++) {
      REAL(result)[j] += tmpscore[j];
    }
    R_CheckUserInterrupt();

  } /* for source <= no_of_nodes */
  
  
  /* clean and return */
  dqueue_destroy(&q);
  stack_destroy(&stack);
  UNPROTECT(3);
  return result;
}

SEXP REST_edge_betweenness (SEXP interface, SEXP graph, SEXP pdirected) {

  REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph);

  SEXP result;
  SEXP dim;
  
  long int no_of_nodes;
  long int no_of_edges;
  dqueue_t q;
  long int *distance;
  long int *nrgeo;
  double *tmpscore;
  stack_t stack;
  long int source;
  long int j;

  SEXP tmp, modein, modeout;

  no_of_nodes=R(ptrtable.vcount(interface, graph));
  no_of_edges=R(ptrtable.ecount(interface, graph));

  if (LOGICAL(pdirected)[0]) {
    PROTECT(modeout=ScalarString(CREATE_STRING_VECTOR("out")));
    PROTECT(modein =ScalarString(CREATE_STRING_VECTOR("in")));
  } else {
    PROTECT(modeout=ScalarString(CREATE_STRING_VECTOR("all")));
    PROTECT(modein =ScalarString(CREATE_STRING_VECTOR("all")));
  }        
  
  distance=(long int*) R_alloc(no_of_nodes, sizeof(long int));
  nrgeo=(long int*) R_alloc(no_of_nodes, sizeof(long int));
  tmpscore=(double*) R_alloc(no_of_nodes, sizeof(double));
  /* others memsetting later */

  PROTECT(result=NEW_NUMERIC(no_of_nodes*no_of_nodes));
  memset(REAL(result), 0, no_of_nodes*no_of_nodes*sizeof(double));
  PROTECT(dim=NEW_INTEGER(2));
  INTEGER(dim)[0]=no_of_nodes;
  INTEGER(dim)[1]=no_of_nodes;
  SET_DIM(result, dim);

  dqueue_init(&q, 100);
  stack_init(&stack, no_of_nodes);

  /* here we go */
  
  for (source=1; source<=no_of_nodes; source++) {
    
    memset(distance, 0, no_of_nodes*sizeof(long int));
    memset(nrgeo, 0, no_of_nodes*sizeof(long int));
    memset(tmpscore, 0, no_of_nodes*sizeof(double));
    stack_clear(&stack); /* it should be empty anyway... */
    
    dqueue_push(&q, source);
    nrgeo[source-1]=1;
    distance[source-1]=0;
    
    while (!dqueue_empty(&q)) {
      long int actnode=dqueue_pop(&q);
      
      tmp=ptrtable.neighbors(interface, graph, actnode, modeout);
      for (j=0; j<GET_LENGTH(tmp); j++) {
	long int neighbor=REAL(tmp)[j];
	if (nrgeo[neighbor-1] != 0) {
	  /* we've already seen this node, another shortest path? */
	  if (distance[neighbor-1]==distance[actnode-1]+1) {
	    nrgeo[neighbor-1]++;
	  }
	} else {
	  /* we haven't seen this node yet */
	  nrgeo[neighbor-1]++;
	  distance[neighbor-1]=distance[actnode-1]+1;
	  dqueue_push(&q, neighbor);
	  stack_push(&stack, neighbor);
	}
      }
    } /* while !dqueue_empty */
    
    /* Ok, we've the distance of each node and also the number of 
       shortest paths to them. Now we do an inverse search, starting
       with the farthest nodes. */
    while (!stack_empty(&stack)) {
      long int actnode=stack_pop(&stack);
      long int friends=0;
      if (distance[actnode-1]<1) { continue; } /* skip source node */
      
      /* search for the neighbors on the geodesics */
      tmp=ptrtable.neighbors(interface, graph, actnode, modein);
      for (j=0; j<GET_LENGTH(tmp); j++) {
	long int neighbor=REAL(tmp)[j];
	if (distance[neighbor-1]==distance[actnode-1]-1) { friends++; }
      }

      /* set the temporary score of the friends */
      for (j=0; j<GET_LENGTH(tmp); j++) {
	long int neighbor=REAL(tmp)[j];
	if (distance[neighbor-1]==distance[actnode-1]-1) { 	  
	  tmpscore[neighbor-1] += (tmpscore[actnode-1]+1)/friends;
	  RMATRIX(result,neighbor, actnode) += (tmpscore[actnode-1]+1)/friends;
	  if (! LOGICAL(pdirected)[0]) {
	    RMATRIX(result, actnode, neighbor) += 
	      (tmpscore[actnode-1]+1)/friends;
	  }
	}
      }
    }
    
    /* Ok, we've the scores for this source */
    R_CheckUserInterrupt();

  } /* for source <= no_of_nodes */
  
  
  /* clean and return */
  dqueue_destroy(&q);
  stack_destroy(&stack);
  UNPROTECT(4);
  return result;

}
