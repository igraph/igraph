/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

#include "igraph_structural.h"
#include "igraph_interface.h"
#include "igraph_conversion.h"

/* HELPER for igraph_i_rich_club_density_sequence()
 *
 * Given an ordered list of vertex IDs, builds a list where the index of each vertex ID
 * corresponds to its position in this returned list.
 *
 * Index: ordering  -> Index: vertex ID
 * Value: vertex ID    Value: ordering
 *
 * order[id] = placement of vertex ID in vertex_order
 */
static igraph_error_t igraph_i_get_vertex_order_map(const igraph_vector_int_t *vertex_order,
                                                    igraph_vector_int_t *map) {
    igraph_integer_t size = igraph_vector_int_size(vertex_order);
    IGRAPH_CHECK(igraph_vector_int_resize(map, size)); // resize map

    for (igraph_integer_t i = 0; i < size; i++) {
        VECTOR(*map)[VECTOR(*vertex_order)[i]] = i;
    }
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_rich_club_density_sequence(const igraph_t *graph,
                                                 const igraph_vector_int_t *vertex_order,
                                                 igraph_bool_t directed,
                                                 igraph_bool_t loops,
                                                 igraph_vector_t *res) {
    /* TODO:
     * implement directed vs. undirected functionality
     */
    igraph_integer_t numVertices = igraph_vcount(graph);

    // Error: vertex_order wrong size
    if (igraph_vector_int_size(vertex_order) != numVertices) {
        IGRAPH_ERROR("Invalid vertex order length.", IGRAPH_EINVAL);
    }

    igraph_integer_t numEdges = igraph_ecount(graph);
    igraph_vector_int_t edges;
    igraph_vector_int_t edgesRemainingAfter;

    IGRAPH_CHECK(igraph_vector_int_init(&edges, 0));
    IGRAPH_CHECK(igraph_vector_int_init(&edgesRemainingAfter, numVertices));

    IGRAPH_CHECK(igraph_vector_resize(res, numVertices)); // resize res

    IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0)); // get list of edges

    // get map of vertex ordering -> orderOf[id] = index of vertexID in vertex_order
    igraph_vector_int_t orderOf;
    IGRAPH_CHECK(igraph_vector_int_init(&orderOf, numVertices));
    igraph_i_get_vertex_order_map(vertex_order, &orderOf);

    for (igraph_integer_t i = 0; i < numEdges; i++) {   // O(E)
        igraph_integer_t v1 = VECTOR(edges)[2 * i];     // endpoints
        igraph_integer_t v2 = VECTOR(edges)[2 * i + 1];
        igraph_integer_t orderV1 = VECTOR(orderOf)[v1]; // ordering of endpoints
        igraph_integer_t orderV2 = VECTOR(orderOf)[v2];

        igraph_integer_t edgeRemovalIndex = (orderV1 < orderV2 ? orderV1 : orderV2);

        VECTOR(edgesRemainingAfter)[edgeRemovalIndex]++; // add to edge remaining count
    }

    // derive the final edgesRemainingAfter vector via accumulating from the end
    igraph_integer_t total = 0;
    for (igraph_integer_t i = numVertices - 1; i >= 0; i--) { // O(V)
        total += VECTOR(edgesRemainingAfter)[i];
        VECTOR(edgesRemainingAfter)[i] = total;
    }

    // density calculation
    for (igraph_integer_t i = 0; i < numVertices; i++) { // O(V)
        // (numVertices - i) = the number of vertices left on this loop
        igraph_real_t totalPossibleEdges = (numVertices - i) * ((numVertices - i) - 1) / 2;
        VECTOR(*res)[i] = VECTOR(edgesRemainingAfter)[i] / totalPossibleEdges;
    }

    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&orderOf);
    igraph_vector_int_destroy(&edgesRemainingAfter);

    return IGRAPH_SUCCESS;

    // Complexity: O(E + V)
}

igraph_error_t igraph_rich_club_coefficient(const igraph_t *graph,
                                            const igraph_vector_t *weights,
                                            igraph_neimode_t mode,
                                            igraph_vector_t *res) {
    /* Computes the rich-club coefficient as a function of degree. */

    return IGRAPH_SUCCESS;
}
