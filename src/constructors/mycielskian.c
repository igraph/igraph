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

#include "igraph_constructors.h"
#include "igraph_conversion.h"

#include "igraph_interface.h"

#include "core/interruption.h"
#include "math/safe_intop.h"

igraph_error_t igraph_mycielskian(igraph_t *res, const igraph_t *graph, igraph_integer_t k) {
    if (k < 0) {
        IGRAPH_ERROR("Number of iterations (k) must be non-negative", IGRAPH_EINVAL);
    }
    if (graph == NULL || res == NULL) {
        IGRAPH_ERROR("Invalid input/output graph", IGRAPH_EINVAL);
    }
    if (igraph_is_directed(graph)) {
        IGRAPH_ERROR("Mycielski's construction is not defined for directed graphs", IGRAPH_EINVAL);
    }
    
    igraph_integer_t vcount = igraph_vcount(graph);
    igraph_integer_t ecount = igraph_ecount(graph);

    if (k == 0) { // if k is 0, return the original graph
        IGRAPH_CHECK(igraph_copy(res, graph));
        return IGRAPH_SUCCESS;
    }

    igraph_integer_t new_vcount = vcount;
    igraph_integer_t new_ecount = ecount;

    for (igraph_integer_t i = 0; i < k; i++) {
        new_ecount = 3 * new_ecount + new_vcount; // the number of edges after each iteration
        new_vcount = new_vcount * 2 + 1; // the new number of vertices after each iteration
    }

    IGRAPH_CHECK(igraph_empty(res, new_vcount, IGRAPH_UNDIRECTED));

    igraph_vector_int_t edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, false)); // copy the edges from the original graph to the new vector
    IGRAPH_CHECK(igraph_vector_int_resize(&edges, new_ecount * 2));

    igraph_integer_t edge_index = 2 * ecount;  // Current last edge index in edge vector
    igraph_integer_t offset = vcount;          // Tracks where new vertices start

    for (igraph_integer_t i = 0; i < k; i++) {
        igraph_integer_t prev_vcount = offset;  // Number of vertices before this step
        igraph_integer_t w = offset * 2;        // The new 'w' node index
        igraph_integer_t last_edge_index = edge_index;  // Mark where edges before this step end

        // For each edge before this step, add two new edges
        for (igraph_integer_t j = 0; j < last_edge_index; j += 2) {
            igraph_integer_t v1 = VECTOR(edges)[j];
            igraph_integer_t v2 = VECTOR(edges)[j + 1];

            VECTOR(edges)[edge_index++] = v1;
            VECTOR(edges)[edge_index++] = offset + v2;

            VECTOR(edges)[edge_index++] = v2;
            VECTOR(edges)[edge_index++] = offset + v1;
        }

        // Add edges connecting each `ui` to `w` (forming a star)
        for (igraph_integer_t j = prev_vcount; j < w; j++) {
            VECTOR(edges)[edge_index++] = j;
            VECTOR(edges)[edge_index++] = w;
        }

        // Update offset for next step
        offset = offset * 2 + 1;
    }

    // Add all edges in one go
    IGRAPH_CHECK(igraph_add_edges(res, &edges, 0));
    igraph_vector_int_destroy(&edges);

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_mycielski_graph(igraph_t *graph, igraph_integer_t k) {
    igraph_t g;
    // igraph_small(&g, 1, 0, -1);

    if (k <= 0) {
        IGRAPH_ERROR("Number of iterations (k) must be positive", IGRAPH_EINVAL);
    }
    if (graph == NULL) {
        IGRAPH_ERROR("Invalid Output graph", IGRAPH_EINVAL);
    }
    if (k == 1) {
        igraph_small(&g, 1, IGRAPH_UNDIRECTED, -1); // single vertex
        igraph_copy(graph, &g);
        igraph_destroy(&g);
        return IGRAPH_SUCCESS;
    }
    igraph_small(&g, 1, IGRAPH_UNDIRECTED, 0, 1, -1); // a path
    igraph_mycielskian(graph, &g, k - 2);
    igraph_destroy(&g);

    return IGRAPH_SUCCESS;
}
