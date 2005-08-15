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

SEXP REST_minimum_spanning_tree_unweighted(SEXP interface, 
					   SEXP graph) {
  
  REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph);
  
  SEXP result;
  long int no_of_nodes;
  char *already_added;
  
  dqueue_t q;
  dqueue_t edges;
  SEXP mode, tmp;
  long int i, j;

  no_of_nodes=R(ptrtable.vcount(interface, graph));
  already_added=(char*) R_alloc(no_of_nodes, sizeof(char));
  memset(already_added, 0, no_of_nodes*sizeof(char));

  PROTECT(mode=NEW_CHARACTER(1));
  SET_STRING_ELT(mode, 0, CREATE_STRING_VECTOR("all"));

  dqueue_init(&q, 100);
  dqueue_init(&edges, (no_of_nodes-1)*2);
  
  for (i=1; i<=no_of_nodes; i++) {
    if (already_added[i-1]>0) { continue; }

    already_added[i-1]=1;
    dqueue_push(&q, i);
    while (! dqueue_empty(&q)) {
      long int act_node=dqueue_pop(&q);
      tmp=ptrtable.neighbors(interface, graph, act_node, mode);
      for (j=0; j<GET_LENGTH(tmp); j++) {
	long int neighbor=REAL(tmp)[j];
	if (already_added[neighbor-1]==0) {
	  already_added[neighbor-1]=1;
	  dqueue_push(&q, neighbor);
	  dqueue_push(&edges, act_node);
	  dqueue_push(&edges, neighbor);
	}
      }
    }
    R_CheckUserInterrupt();
  }
  
  /* Save the result */
  
  dqueue_destroy(&q);
  j=dqueue_size(&edges);
  PROTECT(result=NEW_NUMERIC(j));
  for (i=0; i<j; i++) {
    REAL(result)[i] = dqueue_pop(&edges);
  }
    
  dqueue_destroy(&edges);
  UNPROTECT(2);
  return result;
}

SEXP REST_minimum_spanning_tree_prim(SEXP interface, 
				     SEXP graph) {

  REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph);

  SEXP result;
  long int no_of_nodes;
  char *already_added;

  d_indheap_t heap;
  dqueue_t edges;

  SEXP mode, weight, tmp, tmp2;
  long int i, j;
  
  no_of_nodes=R(ptrtable.vcount(interface, graph));
  already_added=(char*) R_alloc(no_of_nodes, sizeof(char));
  memset(already_added, 0, no_of_nodes*sizeof(char));
  
  PROTECT(mode=ScalarString(CREATE_STRING_VECTOR("all")));
  PROTECT(weight=ScalarString(CREATE_STRING_VECTOR("weight")));
  
  d_indheap_init(&heap, 0);
  dqueue_init(&edges, (no_of_nodes-1)*2);
  
  for (i=1; i<=no_of_nodes; i++) {
    if (already_added[i-1]>0) { continue; }

    already_added[i-1]=1;
    tmp=ptrtable.neighbors(interface, graph, i, mode);
    tmp2=ptrtable.get_edge_attribute(interface, graph, weight, 
				     ScalarReal((double)i),
				     NULL_USER_OBJECT);
    /* add all edges of the first vertex */
    for (j=0; j<GET_LENGTH(tmp); j++) {
      long int neighbor=REAL(tmp)[j];
      if (already_added[neighbor-1] == 0) {
	d_indheap_push(&heap, -REAL(tmp2)[j], i, neighbor);
      }
    }

    while(! d_indheap_empty(&heap)) {
      /* Get minimal edge */
      long int from, to;
      d_indheap_max_index(&heap, &from, &to);

      /* Erase it */
      d_indheap_delete_max(&heap);
      
      /* Does it point to a visited node? */
      if (already_added[to-1]==0) {
	already_added[to-1]=1;
	dqueue_push(&edges, from);
	dqueue_push(&edges, to);
	tmp=ptrtable.neighbors(interface, graph, to, mode);
	tmp2=ptrtable.get_edge_attribute(interface, graph, weight, 
					 ScalarReal((double)to), 
					 NULL_USER_OBJECT);
	/* add all outgoing edges */
	for (j=0; j<GET_LENGTH(tmp); j++) {
	  long int neighbor=REAL(tmp)[j];
	  if (already_added[neighbor-1] == 0) {
	    d_indheap_push(&heap, -REAL(tmp2)[j], to, neighbor);
	  }
	} /* for */
      } /* if !already_added */
    } /* while in the same component */
  } /* for all nodes */

  d_indheap_destroy(&heap);
  j=dqueue_size(&edges);
  PROTECT(result=NEW_NUMERIC(j));
  for (i=0; i<j; i++) {
    REAL(result)[i] = dqueue_pop(&edges);
  }  

  dqueue_destroy(&edges);
  UNPROTECT(3);
  return result;
}
