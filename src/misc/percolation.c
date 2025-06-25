/*
    IGraph library.
    Copyright (C) 2025  The igraph development team <igraph@igraph.org>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
 
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
 
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "igraph_components.h"
#include "igraph_constants.h"
#include "igraph_conversion.h"
#include "igraph_error.h"
#include "igraph_iterators.h"
#include "igraph_structural.h"
#include "igraph_random.h"
#include "igraph_interface.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_vector_list.h"
#include <stdio.h>

/**
 * \function igraph_i_edge_shuffle
 * \brief In-place Fisher-Yates shuffle for edge lists.
 * \param edges List of edges, will be permuted.
 */
static igraph_error_t igraph_i_edge_shuffle(igraph_vector_int_t *edges) {
    igraph_integer_t k;
    igraph_integer_t edge_count = igraph_vector_int_size(edges);
    igraph_integer_t dummy;
    if (edge_count % 2 == 1) IGRAPH_ERROR("Invalid edge list, odd number of elements", IGRAPH_EINVAL);
    edge_count >>= 1;

    RNG_BEGIN();
    while (edge_count > 1) {
        edge_count --;
        k = RNG_INTEGER(0, edge_count);
        dummy = VECTOR(*edges)[edge_count*2];
        VECTOR(*edges)[edge_count*2] = VECTOR(*edges)[k*2];
        VECTOR(*edges)[k*2] = dummy;

        dummy = VECTOR(*edges)[edge_count*2+1];
        VECTOR(*edges)[edge_count*2] = VECTOR(*edges)[k*2+1];
        VECTOR(*edges)[k*2+1] = dummy;
    }
    RNG_END();
    return IGRAPH_SUCCESS;
}
/**
 * \function igraph_i_percolate_edge
 * \brief Percolates a single edge.
 * \param links Vector representing parents.
 * \param sizes sizes[i] is the number of children of links[i]
 * \param biggest The biggest value in sizes, is updated if a bigger cluster is created.
 * \param a A vertex incident to the edge.
 * \param b The other vertex incident to the edge.
 */

static igraph_error_t igraph_i_percolate_edge(igraph_vector_int_t *links,
					      igraph_vector_int_t *sizes,
					      igraph_integer_t* biggest,
					      igraph_integer_t a,
					      igraph_integer_t b) {
    // find head of each tree
    // TODO: Path compression
    while (VECTOR(*links)[a] != a) {
        a = VECTOR(*links)[a];
    }
    while (VECTOR(*links)[b] != b) {
        b = VECTOR(*links)[b];
    }
    
    // if they are already connected, exit early.
    if (a == b) {
      return IGRAPH_SUCCESS;
    }
    // make a child of b
    VECTOR(*links)[a] = b;
    VECTOR(*sizes)[b] += VECTOR(*sizes)[a];
    // if made new biggest component, update biggest.
    if (VECTOR(*sizes)[b] >= *biggest) *biggest = VECTOR(*sizes)[b];
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_i_edge_list_percolation
 * \brief Gives the size of the largest connected component as edges are added.
 * \param edges Vector of edges, will be added in order.
 * \param output output[i] is the size of the largest connected component after edges[i] is added
 */

static igraph_error_t igraph_i_edge_list_percolation(const igraph_vector_int_t *edges,
						     igraph_vector_int_t* output) {
    igraph_integer_t biggest = 0;

    igraph_integer_t lower, upper;

    igraph_vector_int_minmax(edges, &lower, &upper);

    if (lower < 0) IGRAPH_ERROR("Negative number in edge list.", IGRAPH_EINVVID);

    igraph_vector_int_t sizes;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&sizes, upper+1);

    igraph_vector_int_t links;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&links, upper+1);

    for (igraph_integer_t i = 0; i < upper+1; i++) {
        VECTOR(sizes)[i] = 1;
        VECTOR(links)[i] = i;
    }

    igraph_integer_t edge_count = igraph_vector_int_size(edges);
    if (edge_count % 2 == 1) IGRAPH_ERROR("Invalid edge list, odd number of elements", IGRAPH_EINVAL);
    edge_count >>= 1;
    IGRAPH_CHECK(igraph_vector_int_resize(output, edge_count));

    for (igraph_integer_t i = 0; i < edge_count; i++) {
        igraph_i_percolate_edge(&links, &sizes, &biggest, VECTOR(*edges)[2*i], VECTOR(*edges)[2*i+1]);
        VECTOR(*output)[i] = biggest;
    }
    igraph_vector_int_destroy(&sizes);
    igraph_vector_int_destroy(&links);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_bond_percolation
 * \brief calculates the bond-percolation curve of a graph.
 * \param graph The input graph
 * \param output output[i] is the size of the largest component after adding i+1 edges, will be created.
 */

IGRAPH_EXPORT igraph_error_t igraph_bond_percolation(const igraph_t *graph,
						     igraph_vector_int_t * output,
						     igraph_vector_int_t * edge_indices) {
  if (edge_indices == NULL) {
    // generate random vertex list
    igraph_integer_t size = igraph_ecount(graph);
    igraph_vector_int_t new_edges;
    IGRAPH_CHECK(igraph_vector_int_init_range(&new_edges, 0, size));

    RNG_BEGIN();
    IGRAPH_CHECK(igraph_vector_int_shuffle(&new_edges));
    RNG_END();

    IGRAPH_FINALLY(igraph_vector_int_destroy, &new_edges);
    igraph_bond_percolation(graph, output, &new_edges);
    
    igraph_vector_int_destroy(&new_edges);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
  }

  igraph_vector_int_t edgelist;
  IGRAPH_CHECK(igraph_vector_int_init(&edgelist, 2*igraph_vector_int_size(edge_indices)));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &edgelist);
  
  for (igraph_integer_t i = 0; i < igraph_vector_int_size(edge_indices); i++) {
    IGRAPH_CHECK(igraph_edge(graph,
			      VECTOR(*edge_indices)[i],
			     &VECTOR(edgelist)[2*i],
			     &VECTOR(edgelist)[2*i+1]));
  }
  
  igraph_i_edge_list_percolation(&edgelist, output);
    
  igraph_vector_int_destroy(&edgelist);
  IGRAPH_FINALLY_CLEAN(1);
  
  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_percolate_site(const igraph_t *graph,
                                              igraph_vector_int_t *links,
                                              igraph_vector_int_t *sizes,
					      igraph_integer_t *biggest,
                                              igraph_integer_t vertex) {
  if (VECTOR(*sizes)[vertex] != 0) {
    IGRAPH_ERROR("Vertex already added.", IGRAPH_EINVAL);

  }
  
  VECTOR(*sizes)[vertex] = 1;
  
  igraph_vector_int_t neighbors;
  IGRAPH_CHECK(igraph_vector_int_init(&neighbors, 0));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &neighbors);
  
  IGRAPH_CHECK(igraph_neighbors(graph, &neighbors, vertex, IGRAPH_IN));
  
  for (igraph_integer_t i = 0; i < igraph_vector_int_size(&neighbors); i++) {
    if (VECTOR(*sizes)[VECTOR(neighbors)[i]] == 0) continue;
    IGRAPH_CHECK(igraph_i_percolate_edge(links, sizes, biggest, vertex, VECTOR(neighbors)[i]));
  }
  
  igraph_vector_int_destroy(&neighbors);
  IGRAPH_FINALLY_CLEAN(1);
  
  return IGRAPH_SUCCESS;
}

IGRAPH_EXPORT igraph_error_t igraph_site_percolation(const igraph_t *graph,
						     igraph_vector_int_t *output,
						     igraph_vector_int_t *vertices
						     ) {
  // if vertex list is null, generate a random one.
  if (vertices == NULL) {
    // generate random vertex list
    igraph_integer_t size = igraph_vcount(graph);
    igraph_vector_int_t new_vertices;
    IGRAPH_CHECK(igraph_vector_int_init_range(&new_vertices, 0, size));

    RNG_BEGIN();
    IGRAPH_CHECK(igraph_vector_int_shuffle(&new_vertices));
    RNG_END();

    IGRAPH_FINALLY(igraph_vector_int_destroy, &new_vertices);
    igraph_site_percolation(graph, output, &new_vertices);
    
    igraph_vector_int_destroy(&new_vertices);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
  }

  // Initialize variables
  igraph_integer_t size = igraph_vcount(graph);
  if (size != igraph_vector_int_size(vertices)) IGRAPH_ERROR("Graph size and vertex size do not match", IGRAPH_EINVAL);
  igraph_integer_t biggest = 1;
  
  igraph_vector_int_t sizes;
  IGRAPH_CHECK(igraph_vector_int_init(&sizes, size));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &sizes);

  igraph_vector_int_t links;
  IGRAPH_CHECK(igraph_vector_int_init(&links, size));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &links);

  IGRAPH_CHECK(igraph_vector_int_resize(output, size));
  
  for (igraph_integer_t i = 0; i < size; i++) {
    VECTOR(sizes)[i] = 0;
    VECTOR(links)[i] = i;
  }
  
  for (igraph_integer_t i = 0; i < size; i++) {
    IGRAPH_CHECK(igraph_i_percolate_site(graph, &links, &sizes, &biggest, VECTOR(*vertices)[i]));
    VECTOR(*output)[i] = biggest;
  }

  for (igraph_integer_t i = 0; i < size; i++) {
    if (VECTOR(sizes)[i] == 0) IGRAPH_ERROR("Vertex list is missing vertices from graph", IGRAPH_EINVAL);
  }
  // cleanup
  igraph_vector_int_destroy(&links);
  igraph_vector_int_destroy(&sizes);
  IGRAPH_FINALLY_CLEAN(2);
  
  return IGRAPH_SUCCESS;
}
