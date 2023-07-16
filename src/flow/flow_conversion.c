/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2023  The igraph development team

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

#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_interface.h"

#include "flow/flow_internal.h"

/**
 * \function igraph_i_split_vertices
 * \brief Splits each vertex in the graph into an input and an output vertex.
 *
 * This function implements a transformation that allows us to calculate the
 * vertex connectivity of a directed graph either for a specific s-t pair or
 * for all s-t pairs using flows. The transformation splits each vertex into
 * an input vertex and an output vertex. All inbound edges of the original
 * vertex are rewired to point to the input vertex, and all outbound edges of
 * the original vertex are rewired to original from the output vertex, while
 * adding a single directed edge from the input vertex to the output vertex.
 *
 * </para><para>
 * s-t vertex connectivities can then be calculated on this modified graph by
 * setting the capacity of each edge to 1, \em except for the following edges:
 * the edges incident of the input half of the source vertex and the edges
 * incident on the output half of the target vertex. The max flow on this
 * modified graph will be equal to the s-t vertex connectivity of the original
 * graph.
 *
 * </para><para>
 * This function prepares the graph only but does not supply a capacity vector;
 * it is the responsibility of the caller to provide the capacities.
 *
 * </para><para>
 * If the original graph had \em n vertices, he function guarantees that the
 * first \em n vertices of the result graph will correspond to the \em output
 * halves of the vertices and the remaining \em n vertices will correspond to
 * the \em input halves, in the same order as in the original graph.
 *
 * \param graph the input graph
 * \param result an uninitialized graph object; the result will be returned here
 */
igraph_error_t igraph_i_split_vertices(const igraph_t* graph, igraph_t* result) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t i;
    igraph_vector_int_t edges;

    if (!igraph_is_directed(graph)) {
        IGRAPH_ERROR("Input graph must be directed.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, 2 * (no_of_edges + no_of_nodes)));
    IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));
    IGRAPH_CHECK(igraph_vector_int_resize(&edges, 2 * (no_of_edges + no_of_nodes)));

    for (i = 0; i < 2 * no_of_edges; i += 2) {
        igraph_integer_t to = VECTOR(edges)[i + 1];
        VECTOR(edges)[i + 1] = no_of_nodes + to;
    }

    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(edges)[2 * (no_of_edges + i)] = no_of_nodes + i;
        VECTOR(edges)[2 * (no_of_edges + i) + 1] = i;
    }

    IGRAPH_CHECK(igraph_create(result, &edges, 2 * no_of_nodes, IGRAPH_DIRECTED));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
