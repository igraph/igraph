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
  

/** 
 * \function igraph_full_bipartite
 * Create a full bipartite network
 * 
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
 * to them.
 * 
 * \param graph Pointer to an igraph_t object, the graph will be
 *   created here.
 * \param types Pointer to a boolean vector. If not a null pointer,
 *   then the vertex types will be stored here.
 * \param n1 Integer, the number of vertices of the first kind.
 * \param n2 Integer, the number of vertices of the second kind.
 * \param directed Boolean, whether to create a directed graph.
 * \param mode A constant that gives the type of connections for
 *   directed graphs. If \c IGRAPH_OUT, then edges point from vertices
 *   of the first kind to vertices of the second kind; if \c
 *   IGRAPH_IN, then the opposite direction is realized; if \c
 *   IGRAPH_ALL, then mutual edges will be created.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 * 
 * \sa \ref igraph_full() for non-bipartite full graphs.
 */

int igraph_full_bipartite(igraph_t *graph, 
			  igraph_vector_bool_t *types,
			  igraph_integer_t n1, igraph_integer_t n2,
			  igraph_bool_t directed, 
			  igraph_neimode_t mode) {
  
  long int nn1=n1, nn2=n2;
  long int no_of_nodes=nn1+nn2;
  igraph_vector_t edges;
  long int no_of_edges;
  long int ptr=0;
  long int i, j;
  
  if (!directed) {
    no_of_edges=nn1 * nn2;
  } else if (mode==IGRAPH_OUT || mode==IGRAPH_IN) {
    no_of_edges=nn1 * nn2;
  } else { /* mode==IGRAPH_ALL */
    no_of_edges=nn1 * nn2 * 2;
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges*2);

  if (!directed || mode==IGRAPH_OUT) {

    for (i=0; i<nn1; i++) {
      for (j=0; j<nn2; j++) {
	VECTOR(edges)[ptr++]=i;
	VECTOR(edges)[ptr++]=nn1+j;
      }
    }
    
  } else if (mode==IGRAPH_IN) {

    for (i=0; i<nn1; i++) {
      for (j=0; j<nn2; j++) {
	VECTOR(edges)[ptr++]=nn1+j;
	VECTOR(edges)[ptr++]=i;
      }
    }
    
  } else {
    
    for (i=0; i<nn1; i++) {
      for (j=0; j<nn2; j++) {
	VECTOR(edges)[ptr++]=i;
	VECTOR(edges)[ptr++]=nn1+j;
	VECTOR(edges)[ptr++]=nn1+j;
	VECTOR(edges)[ptr++]=i;
      }
    }
  }
  
  IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_FINALLY(igraph_destroy, graph);

  if (types) {
    IGRAPH_CHECK(igraph_vector_bool_resize(types, no_of_nodes));
    igraph_vector_bool_null(types);
    for (i=nn1; i<no_of_nodes; i++) {
      VECTOR(*types)[i] = 1;
    }
  }
  
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/**
 * \function igraph_create_bipartite
 * Create a bipartite graph
 * 
 * This is a simple wrapper function to create a bipartite graph. It
 * does a little more than \ref igraph_create(), e.g. it checks that
 * the graph is indeed bipartite with respect to the given \p types
 * vector. If there is an edge connecting two vertices of the same
 * kind, then an error is reported.
 * \param graph Pointer to an uninitlized graph object, the result is
 *   created here.
 * \param types Boolean vector giving the vertex types. The length of
 *   the vector defines the number of vertices in the graph.
 * \param edges Vector giving the edges of the graph. The highest
 *   vertex id in this vector must be smaller than the length of the
 *   \p types vector.
 * \param directed Boolean scalar, whether to create a directed
 *   graph.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 */
 
int igraph_create_bipartite(igraph_t *graph, const igraph_vector_bool_t *types,
			    const igraph_vector_t *edges, 
			    igraph_bool_t directed) {

  long int no_of_nodes=igraph_vector_bool_size(types);
  long int no_of_edges=igraph_vector_size(edges);
  igraph_real_t min_edge=0, max_edge=0;
  igraph_bool_t min_type=0, max_type=0;
  long int i;

  if (no_of_edges % 2 != 0) {
    IGRAPH_ERROR("Invalid (odd) edges vector", IGRAPH_EINVEVECTOR);
  }
  no_of_edges /= 2;
  
  if (no_of_edges != 0) {
    igraph_vector_minmax(edges, &min_edge, &max_edge);
  }
  if (min_edge < 0 || max_edge >= no_of_nodes) {
    IGRAPH_ERROR("Invalid (negative) vertex id", IGRAPH_EINVVID);
  }

  /* Check types vector */
  if (no_of_nodes != 0) {
    igraph_vector_bool_minmax(types, &min_type, &max_type);
    if (min_type < 0 || max_type > 1) {
      IGRAPH_WARNING("Non-binary type vector when creating a bipartite graph");
    }
  }

  /* Check bipartiteness */
  for (i=0; i<no_of_edges*2; i+=2) {
    long int from=VECTOR(*edges)[i];
    long int to=VECTOR(*edges)[i+1];
    long int t1=VECTOR(*types)[from];
    long int t2=VECTOR(*types)[to];
    if ( (t1 && t2) || (!t1 && !t2) ) {
      IGRAPH_ERROR("Invalid edges, not a bipartite graph", IGRAPH_EINVAL);
    }
  }
  
  IGRAPH_CHECK(igraph_empty(graph, no_of_nodes, directed));
  IGRAPH_FINALLY(igraph_destroy, graph);
  IGRAPH_CHECK(igraph_add_edges(graph, edges, 0));
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}
  
/**
 * \function igraph_incidence
 * Create a bipartite graph from an incidence matrix
 *
 * A bipartite (or two-mode) graph contains two types of vertices and
 * edges always connect vertices of different types. An incidence
 * matrix is an nxm matrix, n and m are the number of vertices of the
 * two types, respectively. Nonzero elements in the matrix denote
 * edges between the two corresponding vertices.
 * 
 * </para><para>
 * Note that this function can operate in two modes, depending on the
 * \p multiple argument. If it is FALSE (i.e. 0), then a single edge is
 * created for every non-zero element in the incidence matrix. If \p
 * multiple is TRUE (i.e. 1), then the matrix elements are rounded up
 * to the closest non-negative integer to get the number of edges to
 * create between a pair of vertices.
 * 
 * </para><para>
 * This function does not create multiple edges if \p multiple is
 * FALSE, but might create some if it is TRUE.
 * 
 * \param graph Pointer to an uninitialized graph object.
 * \param types Pointer to an initialized boolean vector, or a null
 *   pointer. If not a null pointer, then the vertex types are stored
 *   here. It is resized as needed.
 * \param incidence The incidence matrix.
 * \param directed Gives whether to create an undirected or a directed
 *   graph.
 * \param mode Specifies the direction of the edges in a directed
 *   graph. If \c IGRAPH_OUT, then edges point from vertices
 *   of the first kind (corresponding to rows) to vertices of the 
 *   second kind (corresponding to columns); if \c
 *   IGRAPH_IN, then the opposite direction is realized; if \c
 *   IGRAPH_ALL, then mutual edges will be created.
 * \param multiple How to interpret the incidence matrix elements. See
 *   details below.
 * \return Error code.
 * 
 * Time complexity: O(n*m), the size of the incidence matrix.
 */

int igraph_incidence(igraph_t *graph, igraph_vector_bool_t *types,
		     const igraph_matrix_t *incidence, 
		     igraph_bool_t directed,
		     igraph_neimode_t mode, igraph_bool_t multiple) {
  
  long int n1=igraph_matrix_nrow(incidence);
  long int n2=igraph_matrix_ncol(incidence);
  long int no_of_nodes=n1+n2;
  igraph_vector_t edges;
  long int i, j, k;

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

  if (multiple) { 

    for (i=0; i<n1; i++) {
      for (j=0; j<n2; j++) { 
	long int elem=MATRIX(*incidence, i, j);
	long int from, to;
	
	if (!elem) { continue; }
	
	if (mode == IGRAPH_IN) {
	  from=n1+j;
	  to=i;
	} else {
	  from=i;
	  to=n1+j;
	}
	
	if (mode != IGRAPH_ALL || !directed) {
	  for (k=0; k<elem; k++) {
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, from));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
	  }
	} else {
	  for (k=0; k<elem; k++) {
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, from));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, from));
	  }
	}	  
      }
    }

  } else {
    
    for (i=0; i<n1; i++) {
      for (j=0; j<n2; j++) { 
	long int from, to;
	
	if (MATRIX(*incidence, i, j) != 0) {
	  if (mode == IGRAPH_IN) {
	    from=n1+j;
	    to=i;
	  } else {
	    from=i;
	    to=n1+j;
	  }
	  if (mode != IGRAPH_ALL || !directed) {
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, from));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
	  } else {
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, from));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
	    IGRAPH_CHECK(igraph_vector_push_back(&edges, from));
	  }	    
	}
      }
    }
    
  }
  
  IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges); 
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_FINALLY(igraph_destroy, graph);
  
  if (types) {
    IGRAPH_CHECK(igraph_vector_bool_resize(types, no_of_nodes));
    igraph_vector_bool_null(types);
    for (i=n1; i<no_of_nodes; i++) {
      VECTOR(*types)[i] = 1;
    }
  }
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/** 
 * \function igraph_get_incidence
 * Convert a bipartite graph into an incidence matrix
 * 
 * \param graph The input graph, edge directions are ignored.
 * \param types Boolean vector containing the vertex types.
 * \param res Pointer to an initialized matrix, the result is stored
 *   here. An element of the matrix gives the number of edges
 *   (irrespectively of their direction) between the two corresponding
 *   vertices. 
 * \return Error code.
 * 
 * Time complexity: O(n*m), n and m are number of vertices of the two
 * different kind.
 * 
 * \sa \ref igraph_incidence() for the opposite operation.
 */

int igraph_get_incidence(const igraph_t *graph,
			 const igraph_vector_bool_t *types,
			 igraph_matrix_t *res) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  long int n1=0, n2=0, i;

  if (igraph_vector_bool_size(types) != no_of_nodes) {
    IGRAPH_ERROR("Invalid vertex type vector for bipartite graph", 
		 IGRAPH_EINVAL);
  }
  
  for (i=0; i<no_of_nodes; i++) {
    n1 += VECTOR(*types)[i] == 0 ? 1 : 0;
  }
  n2 = no_of_nodes-n1;
  
  IGRAPH_CHECK(igraph_matrix_resize(res, n1, n2));
  igraph_matrix_null(res);
  for (i=0; i<no_of_edges; i++) {
    long int from=IGRAPH_FROM(graph, i);
    long int to=IGRAPH_TO(graph, i);
    if (! VECTOR(*types)[from]) {
      MATRIX(*res, from, to-n1) += 1;
    } else {
      MATRIX(*res, to, from-n1) += 1;
    }
  }
  
  return 0;
}
    
