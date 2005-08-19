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

/*
 */

SEXP REST_edgelist_delete_edges_sorted(SEXP edges, SEXP deleted) {

  SEXP result;
  long int no_of_edges;
  long int edges_to_delete;

  long int eidx=0, didx=0;
  
  dqueue_t q;
  long int i, j;

  no_of_edges=GET_LENGTH(edges)/2;
  edges_to_delete=GET_LENGTH(deleted)/2;
  
  dqueue_init(&q, no_of_edges-edges_to_delete);
  
  while (didx < edges_to_delete) {
    if (REAL(edges)[eidx] < REAL(deleted)[didx] ||
	(REAL(edges)[eidx]==REAL(deleted)[didx] && 
	 REAL(edges)[no_of_edges+eidx] < REAL(deleted)[edges_to_delete+didx])){
      dqueue_push(&q, eidx+1);
      eidx++;
    } else if (REAL(edges)[eidx] > REAL(deleted)[didx] ||
	       (REAL(edges)[eidx]==REAL(deleted)[didx] && 
		REAL(edges)[no_of_edges+eidx] > 
		REAL(deleted)[edges_to_delete+didx])) {
      dqueue_destroy(&q);
      error("No such edge to remove");
    } else {
      eidx++;
      didx++;      
    }
  }
  
  i=dqueue_size(&q);
  PROTECT(result=NEW_NUMERIC(i+no_of_edges-eidx));
  for (j=0; j<i; j++) {
    REAL(result)[j] = dqueue_pop(&q);
  }
  for (j=i; eidx < no_of_edges; eidx++) {
    REAL(result)[j++] = eidx+1;
  }

  dqueue_destroy(&q);
  UNPROTECT(1);
  return result;
}

/*
 */

SEXP REST_edgelist_neighbors(SEXP edges, SEXP pvert, SEXP pmode) {
  
  SEXP result;
  
  long int no_of_edges;
  double v;
  int mode;

  dqueue_t q;
  long int i, j;

  no_of_edges=GET_LENGTH(edges)/2;
  v=R(pvert);
  mode=R(pmode);
  
  dqueue_init(&q, 10);

  if (mode & 1) {
    for (i=0; i<no_of_edges; i++) {
      if (REAL(edges)[i] == v) {
	dqueue_push(&q, (long int)REAL(edges)[no_of_edges+i]);
      }
    }
  } 
  if (mode & 2) {
    for (i=0; i<no_of_edges; i++) {
      if (REAL(edges)[no_of_edges+i] == v) {
	dqueue_push(&q, (long int)REAL(edges)[i]);
      }
    }
  }
  
  i=dqueue_size(&q);
  PROTECT(result=NEW_NUMERIC(i));
  for (j=0; j<i; j++) {
    REAL(result)[j] = dqueue_pop(&q);
  }
  
  dqueue_destroy(&q);
  UNPROTECT(1);
  return result;
}

/**
 */

SEXP REST_edgelist_delete_vertices(SEXP edges, SEXP pnonodes, SEXP vids) {
  
  SEXP result;
  SEXP dim;
  
  long int no_of_nodes;
  long int remaining_nodes=0;
  long int *index;
  long int i, j;
  long int no_of_edges;

  dqueue_t q;

  no_of_nodes=R(pnonodes);
  no_of_edges=GET_LENGTH(edges)/2;
  
  /* Create index */
  index=(long int*)R_alloc(no_of_nodes, sizeof(long int));
  memset(index, 0, sizeof(long int)*no_of_nodes);
  for (i=0; i<GET_LENGTH(vids); i++) {
    index[(long int)REAL(vids)[i]-1] = 1;
  }
  j=1;
  for (i=0; i<no_of_nodes; i++) {
    if (index[i]==0) {
      index[i]=j++;
      remaining_nodes++;
    } else {
      index[i]=0;
    }
    R_CheckUserInterrupt();
  }
  
  PROTECT(result=NEW_LIST(3));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(1));
  REAL(VECTOR_ELT(result, 0))[0] = remaining_nodes;
  
  dqueue_init(&q, no_of_edges * (remaining_nodes / no_of_nodes +1));

  for (i=0; i<no_of_edges; i++) {
    long int from=REAL(edges)[i];
    long int to=REAL(edges)[no_of_edges+i];
    if (index[from-1]!=0 && index[to-1]!=0) {
      dqueue_push(&q, i);
    }
  }

  i=dqueue_size(&q);
  SET_VECTOR_ELT(result, 1, NEW_NUMERIC(i*2));
  SET_VECTOR_ELT(result, 2, NEW_NUMERIC(i));
  PROTECT(dim=NEW_INTEGER(2));
  INTEGER(dim)[0] = i;
  INTEGER(dim)[1] = 2;
  SET_DIM(VECTOR_ELT(result,1), dim);
  for (j=0; j<i; j++) {
    long int id=dqueue_pop(&q);
    REAL(VECTOR_ELT(result, 2))[j]=id+1;
    REAL(VECTOR_ELT(result, 1))[j]=index[ (long int) REAL(edges)[id]-1 ];
    REAL(VECTOR_ELT(result, 1))[i+j]=
      index[ (long int) REAL(edges)[no_of_edges+id]-1 ];    
  }

  dqueue_destroy(&q);
  UNPROTECT(2);
  return result;
}

/**
 */

SEXP REST_edgelist_iterator(SEXP interface, SEXP graph, SEXP ptype) {
  
  SEXP result;
  SEXP names, attrnames;
  SEXP graph_type, graph_id;
  
  char * type;
  int prot=0;

  graph_type=REST_i_get_list_element(graph, "gal");
  graph_id=REST_i_get_list_element(graph_type, "id");
  graph_type=REST_i_get_list_element(graph_type, "type");

  type = CHAR(STRING_ELT(ptype, 0));
  if (!strncmp(type, "vid", 3)) {

    PROTECT(result=NEW_LIST(2));
    SET_VECTOR_ELT(result, 0, NEW_LIST(3));
    SET_VECTOR_ELT(result, 1, NEW_NUMERIC(1));
    SET_VECTOR_ELT(VECTOR_ELT(result, 0), 0, duplicate(graph_type));
    SET_VECTOR_ELT(VECTOR_ELT(result, 0), 1, duplicate(ptype));
    SET_VECTOR_ELT(VECTOR_ELT(result, 0), 2, duplicate(graph_id));
    R(VECTOR_ELT(result, 1))=1.0;

    PROTECT(attrnames=NEW_CHARACTER(3));
    SET_STRING_ELT(attrnames, 0, CREATE_STRING_VECTOR("graph.type"));
    SET_STRING_ELT(attrnames, 1, CREATE_STRING_VECTOR("type"));
    SET_STRING_ELT(attrnames, 2, CREATE_STRING_VECTOR("id"));
    SET_NAMES(VECTOR_ELT(result, 0), attrnames);
    
    prot+=2;
  } else if (!strncmp(type, "eid", 3)) {

    PROTECT(result=NEW_LIST(2));
    SET_VECTOR_ELT(result, 0, NEW_LIST(3));
    SET_VECTOR_ELT(result, 1, NEW_NUMERIC(1));
    SET_VECTOR_ELT(VECTOR_ELT(result, 0), 0, duplicate(graph_type));
    SET_VECTOR_ELT(VECTOR_ELT(result, 0), 1, duplicate(ptype));
    SET_VECTOR_ELT(VECTOR_ELT(result, 0), 2, duplicate(graph_id));
    R(VECTOR_ELT(result, 1))=1.0;

    PROTECT(attrnames=NEW_CHARACTER(3));
    SET_STRING_ELT(attrnames, 0, CREATE_STRING_VECTOR("graph.type"));
    SET_STRING_ELT(attrnames, 1, CREATE_STRING_VECTOR("type"));
    SET_STRING_ELT(attrnames, 2, CREATE_STRING_VECTOR("id"));
    SET_NAMES(VECTOR_ELT(result, 0), attrnames);
    
    prot+=2;
    
  } else {
    error("Unknown iterator type");
  }
  
  /* set names attribute for result & class*/
  PROTECT(names=NEW_CHARACTER(2));
  SET_STRING_ELT(names, 0, CREATE_STRING_VECTOR("attr"));
  SET_STRING_ELT(names, 1, CREATE_STRING_VECTOR("data"));
  SET_NAMES(result, names);  				 
  SET_CLASS(result, ScalarString(CREATE_STRING_VECTOR("igraph.iterator")));
  
  UNPROTECT(prot+1);
  return result;
}

/**
 */

SEXP REST_edgelist_next(SEXP interface, SEXP graph, SEXP it) {

  SEXP result, data, ptype;
  char *type;
  int prot=0;
  
  ptype=REST_i_get_list_element(REST_i_get_list_element(it, "attr"), "type");
  type=CHAR(STRING_ELT(ptype, 0));
  if (!strncmp(type, "vid", 3)) {
    long int no_of_nodes=
      R(REST_i_get_list_element(REST_i_get_list_element(graph, "gal"), "n"));
    PROTECT(result=duplicate(it));      
    data=REST_i_get_list_element(result, "data");
    R(data) += 1;
    if (R(data) > no_of_nodes+1) {
      error("No more elements");
    }
    prot++;
  } else if (!strncmp(type, "eid", 3)) {
    long int no_of_edges=
      I(GET_DIM(REST_i_get_list_element
		(REST_i_get_list_element(graph, "data"), "data")));
    PROTECT(result=duplicate(it));      
    data=REST_i_get_list_element(result, "data");
    R(data) += 1;
    if (R(data) > no_of_edges+1) {
      error("No more elements");
    }
    prot++;
  } else {
    error("Unknown iterator type");
  }

  UNPROTECT(prot);
  return result;
}
    

    
  
