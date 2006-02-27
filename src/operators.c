/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include "error.h"
#include "memory.h"

int igraph_disjoint_union(igraph_t *res, igraph_t *left, igraph_t *right) {

  long int no_of_nodes_left=igraph_vcount(left);
  long int no_of_nodes_right=igraph_vcount(right);
  long int no_of_edges_left=igraph_ecount(left);
  long int no_of_edges_right=igraph_ecount(right);
  igraph_vector_t edges;
  bool_t directed_left=igraph_is_directed(left);
  integer_t from, to;
  long int i;
  
  if (directed_left != igraph_is_directed(right)) {
    IGRAPH_ERROR("Cannot union directed and undirected graphs",
		 IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, 
				     2*(no_of_edges_left+no_of_edges_right)));
  for (i=0; i<no_of_edges_left; i++) {
    igraph_edge(left, i, &from, &to);
    igraph_vector_push_back(&edges, from);
    igraph_vector_push_back(&edges, to);
  }
  for (i=0; i<no_of_edges_right; i++) {
    igraph_edge(right, i, &from, &to);
    igraph_vector_push_back(&edges, from+no_of_nodes_left);
    igraph_vector_push_back(&edges, to+no_of_nodes_left);
  }
  
  IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes_left+no_of_nodes_right, 
			     directed_left));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

int igraph_disjoint_union_many(igraph_t *res, igraph_vector_ptr_t *graphs) {
  long int no_of_graphs=igraph_vector_ptr_size(graphs);
  bool_t directed=1;
  igraph_vector_t edges;
  long int no_of_edges=0;
  long int shift=0;
  igraph_t *graph;
  long int i, j;
  integer_t from, to;
  
  if (no_of_graphs != 0) {
    graph=VECTOR(*graphs)[0];
    directed=igraph_is_directed(graph);
    for (i=0; i<no_of_graphs; i++) {      
      graph=VECTOR(*graphs)[i];
      no_of_edges += igraph_ecount(graph);
      if (directed != igraph_is_directed(graph)) {
	IGRAPH_ERROR("Cannot union directed and undirected graphs", 
		     IGRAPH_EINVAL);
      }
    }
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, 2*no_of_edges));
  
  for (i=0; i<no_of_graphs; i++) {
    long int ec;
    graph=VECTOR(*graphs)[i];    
    ec=igraph_ecount(graph);
    for (j=0; j<ec; j++) {
      igraph_edge(graph, j, &from, &to);
      igraph_vector_push_back(&edges, from+shift);
      igraph_vector_push_back(&edges, to+shift);
    }
    shift += igraph_vcount(graph);
  }
  
  IGRAPH_CHECK(igraph_create(res, &edges, shift, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

int igraph_intersection(igraph_t *res, igraph_t *left, igraph_t *right) {
  
  long int no_of_nodes_left=igraph_vcount(left);
  long int no_of_nodes_right=igraph_vcount(right);
  long int no_of_nodes;
  long int smaller_nodes;
  bool_t directed=igraph_is_directed(left);
  igraph_vector_t edges;
  igraph_vector_t nei1, nei2;
  long int i,j1,j2,n1,n2;
  integer_t v1, v2;

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&nei1, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&nei2, 0);

  no_of_nodes= no_of_nodes_left > no_of_nodes_right ? 
    no_of_nodes_left : no_of_nodes_right;
  smaller_nodes= no_of_nodes_left < no_of_nodes_right ? 
    no_of_nodes_left : no_of_nodes_right;

  if (directed != igraph_is_directed(right)) {
    IGRAPH_ERROR("Cannot intersect directed and undirected graph",
		 IGRAPH_EINVAL);
  }

  for (i=0; i<smaller_nodes; i++) {
    IGRAPH_CHECK(igraph_neighbors(left, &nei1, i, IGRAPH_OUT));
    IGRAPH_CHECK(igraph_neighbors(right, &nei2, i, IGRAPH_OUT));
    igraph_vector_sort(&nei1);
    igraph_vector_sort(&nei2);
    if (!directed) {
      igraph_vector_filter_smaller(&nei1, i);
      igraph_vector_filter_smaller(&nei2, i);
    }
    n1=igraph_vector_size(&nei1);
    n2=igraph_vector_size(&nei2);
    j1=j2=0;
    while (j1<n1 && j2<n2) {
      v1=VECTOR(nei1)[j1]; v2=VECTOR(nei2)[j2];
      if (v1 < v2) {
	j1++;
      } else if (v1 > v2) {
	j2++; 
      } else {
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, v1));
	j1++; 
	j2++;
      }
    }
  }

  IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges);
  igraph_vector_destroy(&nei1);
  igraph_vector_destroy(&nei2);
  IGRAPH_FINALLY_CLEAN(3);
  return 0;
}

void igraph_i_intersection_many_free(igraph_vector_ptr_t *v) {
  long int i, n=igraph_vector_ptr_size(v);
  for (i=0; i<n; i++) { 
    if (VECTOR(*v)[i] != 0) {
      igraph_vector_destroy(VECTOR(*v)[i]);
    }
  }
  igraph_vector_ptr_destroy(v);
}

int igraph_intersection_many(igraph_t *res, igraph_vector_ptr_t *graphs) {

  long int no_of_graphs=igraph_vector_ptr_size(graphs);
  long int no_of_nodes=0;
  bool_t directed=1;
  igraph_vector_t edges;
  igraph_vector_ptr_t neivects;
  igraph_vector_t neiptr;
  long int i, j;
  
  /* Check directedness */
  if (no_of_graphs != 0) {
    directed=igraph_is_directed(VECTOR(*graphs)[0]);
    no_of_nodes=igraph_vcount(VECTOR(*graphs)[0]);
  }
  for (i=1; i<no_of_graphs; i++) {
    if (directed != igraph_is_directed(VECTOR(*graphs)[i])) {
      IGRAPH_ERROR("Cannot intersect directed and undirected graphs",
		   IGRAPH_EINVAL);
    }
  }

  /* Calculate number of nodes */
  for (i=0; i<no_of_graphs; i++) {
    long int n=igraph_vcount(VECTOR(*graphs)[i]);
    if (n < no_of_nodes) {
      no_of_nodes=n;
    }
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&neiptr, no_of_graphs);

  /* Init neighbor vectors */
  if (no_of_graphs != 0) {
    IGRAPH_CHECK(igraph_vector_ptr_init(&neivects, no_of_graphs));
    IGRAPH_FINALLY(igraph_i_intersection_many_free, &neivects);
  }
  for (i=0; i<no_of_graphs; i++) {
    VECTOR(neivects)[i]=Calloc(1, igraph_vector_t);
    if (VECTOR(neivects)[i]==0) { 
      IGRAPH_ERROR("Cannot intersect graphs", IGRAPH_ENOMEM);
    }
    IGRAPH_CHECK(igraph_vector_init(VECTOR(neivects)[i], 0));
  }

  /* Main part */
  for (i=0; i<no_of_nodes; i++) {
    bool_t l;
    
    /* get neighbors */
    for (j=0; j<no_of_graphs; j++) {
      IGRAPH_CHECK(igraph_neighbors(VECTOR(*graphs)[j], VECTOR(neivects)[j], i,
				    IGRAPH_OUT));
      igraph_vector_sort(VECTOR(neivects)[j]);
      if (!directed) {
	igraph_vector_filter_smaller(VECTOR(neivects)[j], i);
      }
    }
    igraph_vector_null(&neiptr);

    /* check if there are more edges */
    l=1; j=0;
    while (l && j<no_of_graphs) {
      l = VECTOR(neiptr)[j] < igraph_vector_size(VECTOR(neivects)[j]);
      j++;
    }

    while (l) {
      
      /* get largest head element and decide whether all head elements
	 are the same */
      bool_t k=1;
      long int head=VECTOR(*(igraph_vector_t*)VECTOR(neivects)[0])
	[(long int)VECTOR(neiptr)[0]];
      long int maxhead=head;
      for (j=1; j<no_of_graphs; j++) {
	long int h=VECTOR(*(igraph_vector_t*)VECTOR(neivects)[j])
	  [(long int)VECTOR(neiptr)[j]];
	k= (k && head==h);
	if (h > maxhead) {
	  maxhead=h;
	}
      }
      
      /* add edge if common, remove head elements */
      if (k) {
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, head));
	for (j=0; j<no_of_graphs; j++) {
	  VECTOR(neiptr)[j] += 1;
	}
      } else {
	for (j=0; j<no_of_graphs; j++) {
	  while (VECTOR(neiptr)[j] < igraph_vector_size(VECTOR(neivects)[j]) &&
		 VECTOR(*(igraph_vector_t*)VECTOR(neivects)[j])
		 [ (long int) VECTOR(neiptr)[j]] < maxhead){
	    VECTOR(neiptr)[j] += 1;
	  }
	}
      }

      /* update l, check if there are more edges */
      l=1; j=0;
      while (l && j<no_of_graphs) {
	l = VECTOR(neiptr)[j] < igraph_vector_size(VECTOR(neivects)[j]);
	j++;
      }
      
    }
  }
      
  if (no_of_graphs != 0) {
    for (i=0; i<no_of_graphs; i++) {
      igraph_vector_destroy(VECTOR(neivects)[i]);
    }
    igraph_vector_ptr_destroy(&neivects);
    IGRAPH_FINALLY_CLEAN(1);
  }

  IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges);
  igraph_vector_destroy(&neiptr);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}
