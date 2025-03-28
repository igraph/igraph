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

static igraph_real_t igraph_i_total_possible_edges(igraph_integer_t numVertices,
                                                   igraph_bool_t directed,
                                                   igraph_bool_t loops){
    igraph_real_t nv = (igraph_real_t) numVertices; // int -> real
    if (!loops) {
        return (directed ? nv * (nv - 1) : nv * (nv - 1) / 2);
    } else {
        return (directed ? nv * nv : nv * (nv + 1) / 2);
    }
}

igraph_error_t igraph_rich_club_density_sequence(const igraph_t *graph,
                                                 const igraph_vector_int_t *vertex_order,
                                                 igraph_bool_t directed,
                                                 igraph_bool_t loops,
                                                 igraph_vector_t *weights,
                                                 igraph_vector_t *res) {
    igraph_integer_t numVertices = igraph_vcount(graph);

    // Error: vertex_order wrong size
    if (igraph_vector_int_size(vertex_order) != numVertices) {
        IGRAPH_ERROR("Invalid vertex order length.", IGRAPH_EINVAL);
    }

    igraph_integer_t numEdges = igraph_ecount(graph);
    igraph_vector_int_t edges;
    igraph_vector_int_t edgesRemainingAfter;
    igraph_vector_int_t orderOf;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edgesRemainingAfter, numVertices);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&orderOf, numVertices);

    IGRAPH_CHECK(igraph_vector_resize(res, numVertices)); // resize res
    IGRAPH_CHECK(igraph_i_get_vertex_order_map(vertex_order, &orderOf)); // get map of vertex order

    igraph_bool_t warningIssued = 0;

    for (igraph_integer_t eid = 0; eid < numEdges; eid++) {
        igraph_integer_t v1 = IGRAPH_FROM(graph, eid);
        igraph_integer_t v2 = IGRAPH_TO(graph, eid);
        if (v1 == v2 && !loops && !warningIssued) { // if loops = false and there's a loop
            IGRAPH_WARNING("Self-loop encountered while `loops = false`. "
                           "Proceeding as if loops were not possible (density computed accordingly)");
            warningIssued = 1;
        }
        igraph_integer_t orderV1 = VECTOR(orderOf)[v1]; // order of endpoints
        igraph_integer_t orderV2 = VECTOR(orderOf)[v2];

        igraph_integer_t edgeRemovalIndex = (orderV1 < orderV2 ? orderV1 : orderV2);

        // add to edge remaining count
        VECTOR(edgesRemainingAfter)[edgeRemovalIndex] += (weights ? VECTOR(*weights)[eid] : 1);
    }

    // derive the final edgesRemainingAfter vector via accumulating from the end
    igraph_integer_t total = 0;
    for (igraph_integer_t i = numVertices - 1; i >= 0; i--) {
        total += VECTOR(edgesRemainingAfter)[i];
        VECTOR(edgesRemainingAfter)[i] = total;
    }

    // density calculation
    for (igraph_integer_t i = 0; i < numVertices; i++) {
        // (numVertices - i) = the number of vertices left on this loop
        igraph_real_t totalPossibleEdges =
                      igraph_i_total_possible_edges(numVertices - i, directed, loops);
        VECTOR(*res)[i] = VECTOR(edgesRemainingAfter)[i] / totalPossibleEdges;
    }

    igraph_vector_int_destroy(&orderOf);
    igraph_vector_int_destroy(&edgesRemainingAfter);
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(3);

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
