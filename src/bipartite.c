/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2008  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph.h"
#include "attributes.h"

/**
 * \section about_bipartite Bipartite networks in igraph
 *
 * <para>
 * A bipartite network contains two kinds of vertices and connections
 * are only possible between two vertices of different kind. There are
 * many natural examples, e.g. movies and actors as vertices and a
 * movie is connected to all participating actors, etc.
 * 
 * </para><para>
 * igraph does not have direct support for bipartite networks, at
 * least not at the C language level. In other words the igraph_t
 * structure does not contain information about the vertex types. 
 * The C functions for bipartite networks usually have an additional
 * input argument to graph, called \c types, a boolean vector giving
 * the vertex types.
 * 
 * </para><para>
 * Most functions creating bipartite networks are able to create this
 * extra vector, you just need to supply an initialized boolean vector
 * to them.</para>
 */

int igraph_i_bipartite_projection(const igraph_t *graph,
				  const igraph_vector_bool_t *types,
				  igraph_t *proj,
				  int which) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int i, j, k, remaining_nodes=0;
  igraph_vector_t vertex_perm, vertex_index;
  igraph_vector_t edges;
  igraph_adjlist_t adjlist;
  igraph_vector_t *neis1, *neis2;
  long int neilen1, neilen2;
  igraph_vector_long_t added;

  if (which < 0) { return 0; }
  
  IGRAPH_VECTOR_INIT_FINALLY(&vertex_perm, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&vertex_perm, no_of_nodes));
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&vertex_index, no_of_nodes);
  IGRAPH_CHECK(igraph_vector_long_init(&added, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &added);
  IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
  
  for (i=0; i<no_of_nodes; i++) {
    if (VECTOR(*types)[i] == which) {
      VECTOR(vertex_index)[i] = ++remaining_nodes;
      igraph_vector_push_back(&vertex_perm, i);
    }
  }

  for (i=0; i<no_of_nodes; i++) {
    if (VECTOR(*types)[i] == which) {
      long int new_i=VECTOR(vertex_index)[i]-1;
      neis1=igraph_adjlist_get(&adjlist, i);
      neilen1=igraph_vector_size(neis1);
      for (j=0; j<neilen1; j++) {
	long int nei=VECTOR(*neis1)[j];
	neis2=igraph_adjlist_get(&adjlist, nei);
	neilen2=igraph_vector_size(neis2);
	for (k=0; k<neilen2; k++) {
	  long int nei2=VECTOR(*neis2)[k], new_nei2;
	  if (nei2 <= i) { continue; }
	  if (VECTOR(added)[nei2] == i+1) { continue; }
	  VECTOR(added)[nei2] = i+1;
	  new_nei2=VECTOR(vertex_index)[nei2]-1;
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, new_i));
	  IGRAPH_CHECK(igraph_vector_push_back(&edges, new_nei2));
	}
      }
    }
  }

  igraph_adjlist_destroy(&adjlist);
  igraph_vector_long_destroy(&added);
  igraph_vector_destroy(&vertex_index);
  IGRAPH_FINALLY_CLEAN(2);
  
  IGRAPH_CHECK(igraph_create(proj, &edges, remaining_nodes, 
			     /*directed=*/ 0));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_FINALLY(igraph_destroy, proj);
  
  IGRAPH_I_ATTRIBUTE_DESTROY(proj);
  IGRAPH_I_ATTRIBUTE_COPY(proj, graph, 1, 0, 0);
  /*  For this we need the new attribute handling interface first... */
  /*   IGRAPH_CHECK(igraph_i_attribute_permute_vertices(graph, proj, &vertex_perm)); */
  igraph_vector_destroy(&vertex_perm);
  IGRAPH_FINALLY_CLEAN(2);
  
  return 0;
}

/**
 * \function igraph_bipartite_projection
 * Create one or both projections of a bipartite (two-mode) network
 * 
 * Creates one or both projections of a bipartite graph. 
 * \param graph The bipartite input graph. Directedness of the edges
 *   is ignored.
 * \param types Boolean vector giving the vertex types of the graph.
 * \param proj1 Pointer to an uninitialized graph object, the first
 *   projection will be created here. It a null pointer, then it is
 *   ignored, see also the \p probe1 argument.
 * \param proj2 Pointer to an uninitialized graph object, the second
 *   projection is created here, if it is not a null pointer. See also
 *   the \p probe1 argument.
 * \return Error code.
 * 
 * Time complexity: O(|V|*d^2+|E|), |V| is the number of vertices, |E|
 * is the number of edges, d is the average (total) degree of the
 * graphs.
 */

int igraph_bipartite_projection(const igraph_t *graph, 
				const igraph_vector_bool_t *types,
				igraph_t *proj1,
				igraph_t *proj2,
				igraph_integer_t probe1) {
  
  long int no_of_nodes=igraph_vcount(graph);

  /* t1 is -1 if proj1 is omitted, it is 0 if it belongs to type zero,
     it is 1 if it belongs to type one. The same for t2 */
  int t1, t2;
  
  if (igraph_vector_bool_size(types) != no_of_nodes) {
    IGRAPH_ERROR("Invalid bipartite type vector size", IGRAPH_EINVAL);
  }
  
  if (probe1 >= no_of_nodes) {
    IGRAPH_ERROR("No such vertex to probe", IGRAPH_EINVAL);
  }
  
  if (probe1 >= 0 && !proj1) {
    IGRAPH_ERROR("`probe1' given, but `proj1' is a null pointer", IGRAPH_EINVAL);
  }
  
  if (probe1 >=0) {
    t1=VECTOR(*types)[(long int)probe1];
    if (proj2) {
      t2=1-t1;
    } else {
      t2=-1;
    }
  } else {
    t1 = proj1 ? 0 : -1;
    t2 = proj2 ? 1 : -1;
  }

  IGRAPH_CHECK(igraph_i_bipartite_projection(graph, types, proj1, t1));
  IGRAPH_FINALLY(igraph_destroy, proj1);
  IGRAPH_CHECK(igraph_i_bipartite_projection(graph, types, proj2, t2));
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}
  
  
  
      
    
