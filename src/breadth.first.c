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

SEXP REST_diameter(SEXP interface, SEXP graph, SEXP pdirected, SEXP punconn) {

  REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph);

  SEXP result;
  long int no_of_nodes;
  long int i, j;
  long int dia=0;
  long int *already_added;
  long int nodes_reached;
  long int unconn;

  dqueue_t q;
  SEXP tmp;
  SEXP dirmode;

  no_of_nodes=R(ptrtable.vcount(interface, graph));
  unconn=LOGICAL(punconn)[0];
  
  if (LOGICAL(pdirected)[0]) {
    PROTECT(dirmode=ScalarString(CREATE_STRING_VECTOR("out")));
  } else {
    PROTECT(dirmode=ScalarString(CREATE_STRING_VECTOR("all")));
  }

  already_added=(long int*) R_alloc(no_of_nodes, sizeof(long int));
  memset(already_added, 0, no_of_nodes*sizeof(long int));

  dqueue_init(&q, 100);
  
  for (i=0; i<no_of_nodes; i++) {
    nodes_reached=1;
    dqueue_push(&q, i);
    dqueue_push(&q, 0);
    already_added[i]=i+1;
    
    while (!dqueue_empty(&q)) {
      long int actnode=dqueue_pop(&q);
      long int actdist=dqueue_pop(&q);
      if (actdist>dia) { dia=actdist; }
      
      tmp=ptrtable.neighbors(interface,graph, actnode+1, dirmode);
      for (j=0; j<GET_LENGTH(tmp); j++) {
	long int neighbor=REAL(tmp)[j]-1;
	if (already_added[neighbor] == i+1) { continue; }
	already_added[neighbor]=i+1;
	nodes_reached++;
	dqueue_push(&q, neighbor);
	dqueue_push(&q, actdist+1);
      }
    } /* while !dqueue_empty */
    
    /* not connected, return largest possible */
    if (nodes_reached != no_of_nodes && !unconn) {
      dia=no_of_nodes;
      break;
    }
    R_CheckUserInterrupt();
  } /* for i<no_of_nodes */
  
  PROTECT(result=NEW_NUMERIC(1));
  R(result)=dia;
  
  /* clean */
  dqueue_destroy(&q);
  UNPROTECT(2);
  return result;
}

SEXP REST_closeness(SEXP interface, SEXP graph, SEXP nodes, SEXP pmode) {
  
  REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph);

  SEXP result;
  
  long int no_of_nodes;
  long int *already_counted;
  long int i, j;
  long int nodes_reached; 

  dqueue_t q;
  
  long int nodes_to_calc;
  SEXP tmp;

  no_of_nodes=R(ptrtable.vcount(interface, graph));
  nodes_to_calc=GET_LENGTH(nodes);

  PROTECT(result=NEW_NUMERIC(nodes_to_calc));
  memset(REAL(result), 0, nodes_to_calc*sizeof(double));
  
  already_counted=(long int*)R_alloc(no_of_nodes, sizeof(long int));
  memset(already_counted, 0, no_of_nodes*sizeof(long int));

  dqueue_init(&q, 100);

  for (i=0; i<nodes_to_calc; i++) {
    dqueue_push(&q, REAL(nodes)[i]-1);
    dqueue_push(&q, 0);
    nodes_reached=1;
    already_counted[(long int)REAL(nodes)[i]-1]=i+1;
    
    while (!dqueue_empty(&q)) {
      long int act=dqueue_pop(&q);
      long int actdist=dqueue_pop(&q);
      REAL(result)[i] += actdist;

      tmp=ptrtable.neighbors(interface,graph, act+1, pmode);
      for (j=0; j<GET_LENGTH(tmp); j++) {
	long int neighbor=REAL(tmp)[j]-1;
	if (already_counted[neighbor] == i+1) { continue; }
	already_counted[neighbor] = i+1;
	nodes_reached++;
	dqueue_push(&q, neighbor);
	dqueue_push(&q, actdist+1);
      }
    }
    REAL(result)[i] += (no_of_nodes * (no_of_nodes-nodes_reached));
    REAL(result)[i] = (no_of_nodes-1) / REAL(result)[i];
    R_CheckUserInterrupt();
  }
  
  /* Clean */
  dqueue_destroy(&q);
  UNPROTECT(1);
  return result;
}

SEXP REST_clusters(SEXP interface, SEXP graph) {

  REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph);

  SEXP result;
  SEXP names;
  
  long int no_of_nodes; 
  char *already_added;
  long int first_node, act_cluster_size=0, no_of_clusters=1;
  
  dqueue_t q;
  dqueue_t cluster_sizes;

  long int i, j;
  SEXP tmp;
  SEXP mode;

  no_of_nodes=R(ptrtable.vcount(interface, graph));
  already_added=(char*) R_alloc(no_of_nodes, sizeof(char));
  memset(already_added, 0, no_of_nodes*sizeof(char));
  PROTECT(mode=NEW_CHARACTER(1));
  SET_STRING_ELT(mode, 0, CREATE_STRING_VECTOR("all"));

  /* Memory for result */
  
  PROTECT(result=NEW_LIST(2));
  PROTECT(names=NEW_LIST(2));
  SET_VECTOR_ELT(names, 0, CREATE_STRING_VECTOR("membership"));
  SET_VECTOR_ELT(names, 1, CREATE_STRING_VECTOR("csize"));
  SET_NAMES(result, names);
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(no_of_nodes));

  dqueue_init(&q, no_of_nodes > 100000 ? 10000 : no_of_nodes/10);
  dqueue_init(&cluster_sizes, 10);

  /* The algorithm */

  for (first_node=0; first_node < no_of_nodes; ++first_node) {
    if (already_added[first_node]==1) continue;

    already_added[first_node]=1;
    act_cluster_size=1;
    REAL(VECTOR_ELT(result, 0))[first_node]=no_of_clusters;
    dqueue_push(&q, first_node);
    
    while ( !dqueue_empty(&q) ) {
      long int act_node=dqueue_pop(&q);
      tmp=ptrtable.neighbors(interface,graph, act_node+1, mode);
      for (i=0; i<GET_LENGTH(tmp); i++) {
	long int neighbor=REAL(tmp)[i]-1;
	if (already_added[neighbor]==1) { continue; }
	dqueue_push(&q, neighbor);
	already_added[neighbor]=1;
	act_cluster_size++;
	REAL(VECTOR_ELT(result, 0))[neighbor]=no_of_clusters;
      }
    }
    no_of_clusters++;
    dqueue_push(&cluster_sizes, act_cluster_size);
    R_CheckUserInterrupt();
  }
  
  /* Save the result */

  dqueue_destroy(&q);
  j=dqueue_size(&cluster_sizes);
  SET_VECTOR_ELT(result, 1, NEW_NUMERIC(j));
  
  for (i=0; i<j; i++) {
    REAL(VECTOR_ELT(result, 1))[i]=dqueue_pop(&cluster_sizes);
  }

  /* Cleaning up */
  
  dqueue_destroy(&cluster_sizes);
  
  UNPROTECT(3);
  return result;
}

SEXP REST_strong_components(SEXP interface, SEXP graph) {

  REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph);

  SEXP result;
  SEXP names;

  long int no_of_nodes;
  long int *out;
  long int *next_nei;
  
  long int i, t=0;

  dqueue_t q, cluster_sizes;
  
  indheap_t h;

  long int no_of_clusters=1;
  long int act_cluster_size;
  SEXP modein, modeout;
  SEXP tmp;

  no_of_nodes=R(ptrtable.vcount(interface, graph));
  PROTECT(modein =ScalarString(CREATE_STRING_VECTOR("in")));
  PROTECT(modeout=ScalarString(CREATE_STRING_VECTOR("out")));  

  /* The result */

  PROTECT(result=NEW_LIST(2));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(no_of_nodes));
  PROTECT(names=NEW_LIST(2));
  SET_VECTOR_ELT(names, 0, CREATE_STRING_VECTOR("membership"));
  SET_VECTOR_ELT(names, 1, CREATE_STRING_VECTOR("csize"));
  SET_NAMES(result, names);
  
  out=(long int*) Calloc(no_of_nodes, long int);
  next_nei=(long int*) R_alloc(no_of_nodes, sizeof(long int));
  memset(next_nei, 0, no_of_nodes*sizeof(long int));

  dqueue_init(&q, 100);

  for (i=0; i<no_of_nodes; i++) {
    if (next_nei[i] > GET_LENGTH(ptrtable.neighbors(interface,graph, i+1, modeout))) { continue; }
    
    dqueue_push(&q, i);
    while (!dqueue_empty(&q)) {
      long int act_node=dqueue_back(&q);
      tmp=ptrtable.neighbors(interface,graph, act_node+1, modeout);
      if (next_nei[act_node]==0) {
	/* this is the first time we've met this vertex */
	next_nei[act_node]++;
      } else if (next_nei[act_node] <=
		 GET_LENGTH(tmp)) {
	/* we've already met this vertex but it has more children */
	long int neighbor=REAL(tmp)[next_nei[act_node]-1]-1;
	if (next_nei[neighbor] == 0) {
	  dqueue_push(&q, neighbor);
	}
	next_nei[act_node]++;
      } else {
	/* we've met this vertex and it has no more children */
	t++;
	out[act_node]=t;
	dqueue_pop_back(&q);
      }
    } /* while q */    
    R_CheckUserInterrupt();
  }  /* for */

  /* OK, we've the 'out' values for the nodes, let's use them in
     descreasing order with the help of a heap */

  indheap_init_array(&h, out, no_of_nodes);
  Free(out); /* gain some memory... */
  memset(next_nei, 0, no_of_nodes*sizeof(long int)); /* mark already
							added vertices */
  dqueue_init(&cluster_sizes, 100);
  
  while (!indheap_empty(&h)) {
    long int grandfather=indheap_max_index(&h)-1;
    indheap_delete_max(&h);
    if (next_nei[grandfather] != 0) { continue; }
    next_nei[grandfather]=1;
    act_cluster_size=1;
    REAL(VECTOR_ELT(result, 0))[grandfather]=no_of_clusters;
    dqueue_push(&q, grandfather);
    
    while (!dqueue_empty(&q)) {
      long int act_node=dqueue_pop_back(&q);
      tmp=ptrtable.neighbors(interface,graph, act_node+1, modein);
      for (i=0; i<GET_LENGTH(tmp); i++) {
	long int neighbor=REAL(tmp)[i]-1;
	if (next_nei[neighbor] != 0) { continue; }
	dqueue_push(&q, neighbor);
	next_nei[neighbor]=1;
	act_cluster_size++;
	REAL(VECTOR_ELT(result, 0))[neighbor]=no_of_clusters;
      }
    }
    no_of_clusters++;
    dqueue_push(&cluster_sizes, act_cluster_size);
    R_CheckUserInterrupt();
  }
  
  /* Save the result */
  
  dqueue_destroy(&q);
  indheap_destroy(&h);
  t=dqueue_size(&cluster_sizes);
  SET_VECTOR_ELT(result, 1, NEW_NUMERIC(t));
  
  for (i=0; i<t; i++) {
    REAL(VECTOR_ELT(result, 1))[i]=dqueue_pop(&cluster_sizes);
  }
  
  /* Clean up, return */
  
  dqueue_destroy(&cluster_sizes);
  UNPROTECT(4);
  return result;
}

SEXP REST_shortest_paths(SEXP interface, SEXP graph, SEXP from, SEXP pmode) {
  
  REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph);

  SEXP result;
  SEXP dim;
  
  long int no_of_nodes;
  long int no_of_from;
  long int *already_counted;
  
  dqueue_t q;

  long int i, j;
  SEXP tmp;

  no_of_nodes=R(ptrtable.vcount(interface, graph));
  no_of_from=GET_LENGTH(from);
  
  already_counted=(long int*) R_alloc(no_of_nodes, sizeof(long int));
  memset(already_counted, 0, no_of_nodes*sizeof(long int));

  PROTECT(result=NEW_NUMERIC(no_of_from*no_of_nodes));
  PROTECT(dim=NEW_INTEGER(2));
  INTEGER(dim)[0]=no_of_nodes;
  INTEGER(dim)[1]=no_of_from;
  SET_DIM(result, dim);
  memset(REAL(result), 0, no_of_nodes*no_of_from*sizeof(double));

  dqueue_init(&q, 100);

  for (i=0; i<no_of_from; i++) {
    long int reached=1;
    dqueue_push(&q, REAL(from)[i]);
    dqueue_push(&q, 0);
    already_counted[ (long int) REAL(from)[i]-1 ] = i+1;
    
    while (!dqueue_empty(&q)) {
      long int act=dqueue_pop(&q);
      long int actdist=dqueue_pop(&q);
      RMATRIX(result, act, i+1)=actdist;
      
      tmp=ptrtable.neighbors(interface,graph, act, pmode);
      for (j=0; j<GET_LENGTH(tmp); j++) {
	long int neighbor=REAL(tmp)[j];
	if (already_counted[neighbor-1] == i+1) { continue; }
	already_counted[neighbor-1] = i+1;
	reached++;
	dqueue_push(&q, neighbor);
	dqueue_push(&q, actdist+1);
      }
    }
    /* Plus the unreachable nodes */
    j=1;
    while (reached < no_of_nodes) {
      if (RMATRIX(result, j, i+1) == 0 && j != REAL(from)[i]) {
	RMATRIX(result, j, i+1)=no_of_nodes;
	reached++;
      }
      j++;
    }
    R_CheckUserInterrupt();
  }

  /* Clean */
  dqueue_destroy(&q);
  UNPROTECT(2);
  return result;
}

SEXP REST_subcomponent(SEXP interface, SEXP graph, SEXP pvertex, SEXP pmode) {
  
  REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph);

  SEXP result;
  
  long int no_of_nodes;
  long int vertex;
  dqueue_t q;
  dqueue_t comp;
  char *already_added;
  long int i,j;
  SEXP tmp;
  
  no_of_nodes=R(ptrtable.vcount(interface, graph));
  vertex=R(pvertex);

  already_added=(char*) R_alloc(no_of_nodes, sizeof(char));
  memset(already_added, 0, sizeof(char)*no_of_nodes);  

  dqueue_init(&q, 100);
  dqueue_init(&comp, 100);
  
  dqueue_push(&q, vertex);
  dqueue_push(&comp, vertex);
  already_added[vertex-1]=1;
  
  while (!dqueue_empty(&q)) {
    long int actnode=dqueue_pop(&q);
    
    tmp=ptrtable.neighbors(interface,graph, actnode, pmode);
    for (i=0; i<GET_LENGTH(tmp); i++) {
      long int neighbor=REAL(tmp)[i];
      
      if (already_added[neighbor-1]) { continue; }
      already_added[neighbor-1]=1;
      dqueue_push(&comp, neighbor);
      dqueue_push(&q, neighbor);
    }
    R_CheckUserInterrupt();
  }

  dqueue_destroy(&q);

  j=dqueue_size(&comp);
  PROTECT(result=NEW_NUMERIC(j));
  for (i=0; i<j; i++) {
    REAL(result)[i]=dqueue_pop(&comp);
  }
  
  /* Clean */
  
  dqueue_destroy(&comp);
  UNPROTECT(1);
  return result;    
}
