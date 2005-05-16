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

SEXP REST_i_adjacencylist_vcount(SEXP interface, SEXP graph) {
  SEXP data=REST_i_get_list_element(graph, "data");
  SEXP out=REST_i_get_list_element(data, "out");
  SEXP result;
  
  PROTECT(result=NEW_NUMERIC(1));
  R(result)=GET_LENGTH(out);
  UNPROTECT(1);
  return result;  
}

SEXP REST_i_adjacencylist_neighbors(SEXP interface, SEXP graph, 
				    long int vertex, SEXP mode) {
  
  SEXP result;
  SEXP dir;
  SEXP out, inc;

  dir=REST_i_get_list_element(graph, "gal");
  dir=REST_i_get_list_element(dir, "directed");

  if (LOGICAL(dir)[0]) {
    if (!strcmp(STRING_VALUE(mode), "in")) {
      inc=REST_i_get_list_element(graph, "data");
      inc=VECTOR_ELT(REST_i_get_list_element(inc, "inc"), vertex-1);
      PROTECT(result=NEW_NUMERIC(GET_LENGTH(inc)));
      memcpy(REAL(result), REAL(inc), sizeof(double) * GET_LENGTH(result));
    } else if (!strcmp(STRING_VALUE(mode), "out")) {
      out=REST_i_get_list_element(graph, "data");
      out=VECTOR_ELT(REST_i_get_list_element(out, "out"), vertex-1);
      PROTECT(result=NEW_NUMERIC(GET_LENGTH(out)));
      memcpy(REAL(result), REAL(out), sizeof(double) * GET_LENGTH(result));
    } else {
      long int len;
      inc=REST_i_get_list_element(graph, "data");
      inc=VECTOR_ELT(REST_i_get_list_element(inc, "inc"), vertex-1);
      out=REST_i_get_list_element(graph, "data");
      out=VECTOR_ELT(REST_i_get_list_element(out, "out"), vertex-1);
      len=GET_LENGTH(inc)+GET_LENGTH(out);
      PROTECT(result=NEW_NUMERIC(len));
      memcpy(REAL(result), REAL(inc), sizeof(double)*GET_LENGTH(inc));
      memcpy(REAL(result) + GET_LENGTH(inc), REAL(out), 
	     sizeof(double)*GET_LENGTH(out));
    }
  } else {
    out=REST_i_get_list_element(graph, "data");
    out=VECTOR_ELT(REST_i_get_list_element(out, "out"), vertex-1);
    PROTECT(result=NEW_NUMERIC(GET_LENGTH(out)));
    memcpy(REAL(result), REAL(out), sizeof(double) * GET_LENGTH(result));
  }
  
  UNPROTECT(1);
  return result;
}    



/**
 */

SEXP REST_add_edges_adjacencylist(SEXP data, SEXP edges, 
				  SEXP pdirected) {
  
  SEXP out, inc=0;
  SEXP result;
  SEXP names;
  SEXP newout, newinc=0;

  long int *length;
  int directed;
  long int no_of_nodes;
  long int no_of_edges;
  long int i;
  
  no_of_nodes=GET_LENGTH(VECTOR_ELT(data, 0));
  directed=LOGICAL(pdirected)[0];
  out=VECTOR_ELT(data, 0);  
  if (directed) {
    inc=VECTOR_ELT(data, 1);
  }
  
  length=(long int*) R_alloc(no_of_nodes, sizeof(long int));

  /* create result */

  if (directed) {
    PROTECT(result=NEW_LIST(2));
    SET_VECTOR_ELT(result, 0, NEW_LIST(no_of_nodes));
    SET_VECTOR_ELT(result, 1, NEW_LIST(no_of_nodes));
    PROTECT(names=NEW_LIST(2));
    SET_VECTOR_ELT(names, 0, CREATE_STRING_VECTOR("out"));
    SET_VECTOR_ELT(names, 1, CREATE_STRING_VECTOR("inc"));
    SET_NAMES(result, names);
    newout=VECTOR_ELT(result, 0);
    newinc=VECTOR_ELT(result, 1);
  } else { 
    PROTECT(result=NEW_LIST(1));
    SET_VECTOR_ELT(result, 0, NEW_LIST(no_of_nodes));
    PROTECT(names=NEW_LIST(1));
    SET_VECTOR_ELT(names, 0, CREATE_STRING_VECTOR("out"));
    SET_NAMES(result, names);
    newout=VECTOR_ELT(result, 0);
  }
  
  /*********************/
  /* I. out only first */
  /*********************/

  for (i=0; i<no_of_nodes; i++) {
    length[i]=GET_LENGTH(VECTOR_ELT(out, i));
    R_CheckUserInterrupt();
  }
  
  for (i=0; i<GET_LENGTH(edges); i+=2) {
    long int from=REAL(edges)[i];
    long int to=REAL(edges)[i+1];
    length[from-1] ++;
    if (!directed) { length[to-1]++; }
    R_CheckUserInterrupt();
  }
  
  /* allocate memory, copy old elements */
  for (i=0; i<no_of_nodes; i++) {
    SET_VECTOR_ELT(newout, i, NEW_NUMERIC(length[i]));
    memcpy(REAL(VECTOR_ELT(newout, i)), REAL(VECTOR_ELT(out, i)),
	   GET_LENGTH(VECTOR_ELT(out, i))*sizeof(double));
    R_CheckUserInterrupt();
  }

  /* add the new edges */
  for (i=0; i<GET_LENGTH(edges); i+=2) {
    long int from=REAL(edges)[i];
    long int to=REAL(edges)[i+1];
    REAL(VECTOR_ELT(newout, from-1))[ --length[from-1] ] = to;
    if (!directed) {     
      REAL(VECTOR_ELT(newout, to-1))[ --length[to-1] ] = from;
    }
    R_CheckUserInterrupt();
  }

  /*******************************/
  /* II. inc for directed graphs */
  /*******************************/
  
  if (directed) {
    for (i=0; i<no_of_nodes; i++) {
      length[i]=GET_LENGTH(VECTOR_ELT(inc, i));
      R_CheckUserInterrupt();
    }
  
    for (i=0; i<GET_LENGTH(edges); i+=2) {
      long int to=REAL(edges)[i+1];
      length[to-1]++;
      R_CheckUserInterrupt();
    }
    
    /* allocate memory, copy old elements */
    for (i=0; i<no_of_nodes; i++) {
      SET_VECTOR_ELT(newinc, i, NEW_NUMERIC(length[i]));
      memcpy(REAL(VECTOR_ELT(newinc, i)), REAL(VECTOR_ELT(inc, i)),
	     GET_LENGTH(VECTOR_ELT(inc, i))*sizeof(double));
      R_CheckUserInterrupt();
    }
    
    /* add the new edges */
    for (i=0; i<GET_LENGTH(edges); i+=2) {
      long int from=REAL(edges)[i];
      long int to=REAL(edges)[i+1];
      REAL(VECTOR_ELT(newinc, to-1))[ --length[to-1] ] = from;
      R_CheckUserInterrupt();
    }
  } /* if directed */
  
  /* That's it */
  UNPROTECT(2);
  return result;
}

/**
 */
SEXP REST_neighborlist_delete_vertices(SEXP neilist, SEXP vids) {
  
  SEXP result;

  long int no_of_nodes;
  long int remaining_nodes=0;
  long int *index;
  long int i, j, k, l;

  no_of_nodes=GET_LENGTH(neilist);

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
  
  PROTECT(result=NEW_LIST(remaining_nodes));
  
  k=0;
  for (i=0; i<no_of_nodes; i++) {
    if (index[i]==0) { continue; }  /* deleted node */

    l=0;
    for (j=0; j<GET_LENGTH(VECTOR_ELT(neilist, i)); j++) {
      if (index[(long int) REAL(VECTOR_ELT(neilist, i))[j]-1] !=0 ) { l++; }
    }    
    
    SET_VECTOR_ELT(result, k, NEW_NUMERIC(l));

    for (j=0; j<GET_LENGTH(VECTOR_ELT(neilist, i)); j++) {
      long int idx=index[(long int) REAL(VECTOR_ELT(neilist, i))[j]-1];
      if (idx!=0) {
	REAL(VECTOR_ELT(result, k))[ --l ] = idx;
      }
    }
    k++;        
    R_CheckUserInterrupt();
  }
  
  /* Clean up */
  UNPROTECT(1);
  return result;
}
