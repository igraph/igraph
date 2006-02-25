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

int igraph_union(igraph_t *res, igraph_t *left, igraph_t *right) {

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

int igraph_union_many(igraph_t *res, igraph_vector_ptr_t *graphs) {
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
