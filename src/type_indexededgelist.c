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

#define EDGE(i) (REAL(el)[ (long int) REAL(index)[(i)-1] -1 ])

SEXP REST_indexededgelist_create_startindex(SEXP el, SEXP index, SEXP nodes) {
    
  SEXP result, SEXP ;
  long int no_of_nodes;
  long int no_of_edges;
  long int i, j, idx;
  
  no_of_nodes=R(nodes);
  no_of_edges=GET_LENGTH(el);
  
  /* result */
  
  PROTECT(result=NEW_NUMERIC(no_of_nodes+1));
  
  /* create the index */

  idx=0;
  for (i=1; i<=EDGE(1); i++) {
     idx++; REAL(result)[idx-1]=1;
  }
  for (i=1; i<=no_of_edges; i++) {
    long int n=EDGE(i) - EDGE((long int)REAL(result)[idx-1]);
    for (j=0; j<n; j++) {
      idx++; REAL(result)[idx-1]=i;
    }
  }
  j=EDGE((long int)REAL(result)[idx-1]);
  for (i=0; i<=no_of_nodes-j; i++) {
    idx++; REAL(result)[idx-1]=no_of_edges+1;
  }

  /* clean */

  UNPROTECT(1);
  return result;
}

SEXP REST_indexededgelist_neighbors(SEXP interface, SEXP graph, SEXP v,
				    SEXP pmode) {
  
  SEXP result;
  long int node;
  int mode;
  long int length=0, idx=0;
  long int no_of_edges;
  SEXP data, el, is, ii, os, oi;
  long int i;

  node=R(v);
  mode=R(pmode);
  data=REST_i_get_list_element(graph, "data");
  el=REST_i_get_list_element(data, "el");
  no_of_edges=GET_LENGTH(el)/2;

  /* Calculate needed space first & allocate it*/

  if (mode & 1) {
    os=REST_i_get_list_element(data, "os");
    oi=REST_i_get_list_element(data, "oi");
    length += (REAL(os)[node-1+1] - REAL(os)[node-1]);
  }
  if (mode & 2) {
    is=REST_i_get_list_element(data, "is");
    ii=REST_i_get_list_element(data, "ii");
    length += (REAL(is)[node-1+1] - REAL(is)[node-1]);
  }
  
  PROTECT(result=NEW_NUMERIC(length));
  
  if (mode & 1) {
    for (i=REAL(os)[node-1]; i<REAL(os)[node+1-1]; i++) {
      REAL(result)[idx++] = REAL(el)[ (long int)REAL(oi)[i-1]+no_of_edges-1 ];
    }
  }
  if (mode & 2) {
    for (i=REAL(is)[node-1]; i<REAL(is)[node+1-1]; i++) {
      REAL(result)[idx++] = REAL(el)[ (long int)REAL(ii)[i-1]-1 ];
    }
  }

  UNPROTECT(1);
  return result;
}

/**
 * Modifies the graph!!!! Dirty hack.....
 */

SEXP REST_indexededgelist_delete_edges(SEXP interface, SEXP graph, 
				       SEXP edges, SEXP pdirected) {

  REST_i_ptrtable_t ptrtable=REST_i_getptrtable(graph);
  
  SEXP result;
  SEXP data, el, is, ii, os, oi;
  int directed;
  long int no_of_edges;
  long int edges_to_delete;
  long int i;

  directed=LOGICAL(pdirected)[0];
  edges_to_delete=GET_LENGTH(edges)/2;
  data=REST_i_get_list_element(graph, "data");
  el=REST_i_get_list_element(data, "el");
  os=REST_i_get_list_element(data, "os");
  is=REST_i_get_list_element(data, "is");
  oi=REST_i_get_list_element(data, "oi");
  ii=REST_i_get_list_element(data, "ii");
  no_of_edges=GET_LENGTH(el)/2;
    
  /* result */

  PROTECT(result=NEW_NUMERIC(edges_to_delete));
  
  for (i=0; i<edges_to_delete; i++) {
    long int from=REAL(edges)[2*i];
    long int to  =REAL(edges)[2*i+1];
    long int j=REAL(os)[from-1];
    long int d=-1;
    while (d==-1 && j<REAL(os)[from+1-1]) {
      long int idx=REAL(oi)[j-1];
      if ( REAL(el)[idx+no_of_edges-1] == to) {
	d=idx;
      }
      j++;
    }
    if (d!=-1) {
      REAL(el)[d-1]=0;
      REAL(el)[d+no_of_edges-1]=0;
    }
    if (! directed && d==-1) {
      j=REAL(is)[from-1];
      while(d==-1 && j<REAL(is)[from+1-1]) {
	long int idx=REAL(ii)[j-1];
	if( REAL(el)[idx-1] == to) {
	  d=idx;
	}
	j++;
      }
      if (d!=-1) {
	REAL(el)[d-1]=0;
	REAL(el)[d+no_of_edges-1]=0;
      }
    }
    if (d==-1) {
      UNPROTECT(1);
      error("No such edge to delete");
    }
    REAL(result)[i]=d;    
  }
  
  /* clean  */
  
  UNPROTECT(1);
  return result;
}

/*
 * Modifies the graph!!!
 */

SEXP REST_indexededgelist_delete_vertices(SEXP interface, SEXP graph,
					  SEXP vids, SEXP pdirected) {

  SEXP result;
  dqueue_t q;
  int directed;
  SEXP data, el, is, ii, os, oi;
  long int no_of_edges;
  long int vertices_to_delete;
  long int i, j;

  directed=LOGICAL(pdirected)[0];
  vertices_to_delete=GET_LENGTH(vids);
  data=REST_i_get_list_element(graph, "data");
  el=REST_i_get_list_element(data, "el");
  os=REST_i_get_list_element(data, "os");
  is=REST_i_get_list_element(data, "is");
  oi=REST_i_get_list_element(data, "oi");
  ii=REST_i_get_list_element(data, "ii");
  no_of_edges=GET_LENGTH(el)/2;

  dqueue_init(&q, 100);

  for (i=0; i<vertices_to_delete; i++) {
    long int vid=REAL(vids)[i];
    for (j=REAL(os)[vid-1]; j<REAL(os)[vid+1-1]; j++) {
      long int idx=REAL(oi)[j-1];
      REAL(el)[idx-1]=0;
      dqueue_push(&q, idx);
    }
    if (!directed) {
      for (j=REAL(is)[vid-1]; j<REAL(is)[vid+1-1]; j++) {
	long int idx=REAL(ii)[j-1];
	REAL(el)[idx-1]=0;
	dqueue_push(&q, idx);
      }
    }
  }
  
  i=dqueue_size(&q);
  PROTECT(result=NEW_NUMERIC(i));
  for (j=0; j<i; j++) {
    REAL(result)[j]=dqueue_pop(&q);
  }
  
  UNPROTECT(1);
  return result;  
}

/**
 */

SEXP REST_indexededgelist_degree(SEXP interface, SEXP graph, 
				 SEXP vids, SEXP pmode, SEXP ploops) {
  
  SEXP result;
  SEXP data, os, is;
  long int nodes_to_calc;
  long int i, j;
  long int mode;
  int loops;
  
  nodes_to_calc=GET_LENGTH(vids);
  mode=R(pmode);
  loops=LOGICAL(ploops)[0];
  
  data=REST_i_get_list_element(graph, "data");
  os=REST_i_get_list_element(data, "os");
  is=REST_i_get_list_element(data, "is");
  
  PROTECT(result=NEW_NUMERIC(nodes_to_calc));
  memset(REAL(result), 0, nodes_to_calc*sizeof(double));

  if (loops) {
    if (mode & 1) {
      for (i=0; i<nodes_to_calc; i++) {
	long int vid=REAL(vids)[i];
	REAL(result)[i] += (REAL(os)[vid+1-1]-REAL(os)[vid-1]);
      }
    }
    if (mode & 2) {
      for (i=0; i<nodes_to_calc; i++) {
	long int vid=REAL(vids)[i];
	REAL(result)[i] += (REAL(is)[vid+1-1]-REAL(is)[vid-1]);
      }
    }
  } else { /* no loops */
    SEXP oi=REST_i_get_list_element(data, "oi");
    SEXP ii=REST_i_get_list_element(data, "ii");
    SEXP el=REST_i_get_list_element(data, "el");
    long int no_of_edges=GET_LENGTH(el)/2;
    if (mode & 1) {
      for (i=0; i<nodes_to_calc; i++) {
	long int vid=REAL(vids)[i];
	REAL(result)[i] += (REAL(os)[vid+1-1]-REAL(os)[vid-1]);
	for (j=REAL(os)[vid-1]; j<REAL(os)[vid+1-1]; j++) {
	  if (REAL(el)[ (long int)REAL(oi)[j-1]+no_of_edges-1 ]==vid) {
	    REAL(result)[i] -= 1;
	  }
	}
      }
    }
    if (mode & 2) {
      for (i=0; i<nodes_to_calc; i++) {
	long int vid=REAL(vids)[i];
	REAL(result)[i] += (REAL(is)[vid+1-1]-REAL(is)[vid-1]);
	for (j=REAL(is)[vid-1]; j<REAL(is)[vid+1-1]; j++) {
	  if (REAL(el)[ (long int)REAL(ii)[j-1]-1 ]==vid) {
	    REAL(result)[i] -= 1;
	  }
	}
      }
    } 
  }  /* loops */

  UNPROTECT(1);
  return result;
}
   
/**
 */

SEXP REST_i_indexededgelist_neighbors(SEXP interface, SEXP graph, 
				      long int vertex, SEXP pmode){
  
  double mode;
  SEXP dir;
  const char *str;

  dir=REST_i_get_list_element(graph, "gal");
  dir=REST_i_get_list_element(dir, "directed");
  
  if (LOGICAL(dir)[0]) {
    str=CHAR(STRING_ELT(pmode, 0));  
    if (!strcmp(str, "out")) {
      mode=1;
    } else if (!strcmp(str, "in")) {
      mode=2;
    } else if (!strcmp(str, "all") || !strcmp(str,"total")) {
      mode=3;
    }
  } else {
    mode=3;
  }

  return REST_indexededgelist_neighbors(interface, graph, 
					ScalarReal((double)vertex), 
					ScalarReal(mode));
}
