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
        IGRAPH_ERROR("The number of Mycielski iterations must not be negative.", IGRAPH_EINVAL);
    }
    
    igraph_integer_t vcount = igraph_vcount(graph);
    igraph_integer_t ecount = igraph_ecount(graph);

    igraph_integer_t new_vcount = vcount;
    igraph_integer_t new_ecount = ecount;

    if (vcount == 0) { // empty graph
        if (k <= 1) { // 0-> return null graph, 1-> return a single vertex
            IGRAPH_CHECK(igraph_empty(res, k, IGRAPH_UNDIRECTED));
            return IGRAPH_SUCCESS;
        }
        igraph_t g;
        // create a path 0---1
        IGRAPH_CHECK(igraph_ring(&g, 2, IGRAPH_UNDIRECTED, 0, 0));
        IGRAPH_FINALLY(igraph_destroy, &g);

        igraph_mycielskian(res, &g, k - 2);

        igraph_destroy(&g);
        return IGRAPH_SUCCESS;
    }

    if (vcount == 1) { // single vertex, assuming no self loop
        if (k == 0) { // 0-> return single vertex
            IGRAPH_CHECK(igraph_empty(res, 1, IGRAPH_UNDIRECTED));
            return IGRAPH_SUCCESS;
        }
        igraph_t g;
        // create a path 0---1
        IGRAPH_CHECK(igraph_ring(&g, 2, IGRAPH_UNDIRECTED, 0, 0));
        IGRAPH_FINALLY(igraph_destroy, &g);

        igraph_mycielskian(res, &g, k - 1);

        igraph_destroy(&g);
        return IGRAPH_SUCCESS;
    }

    for (igraph_integer_t i = 0; i < k; i++) {
        IGRAPH_SAFE_MULT(new_ecount, 3, &new_ecount);
        IGRAPH_SAFE_ADD(new_ecount, new_vcount, &new_ecount); // new edges = 3 * old edges + old vertices

        IGRAPH_SAFE_MULT(new_vcount, 2, &new_vcount);
        IGRAPH_SAFE_ADD(new_vcount, 1, &new_vcount); // new vertices = 2 * old vertices + 1
    }

    if (igraph_is_directed(graph))
        IGRAPH_CHECK(igraph_empty(res, new_vcount, IGRAPH_DIRECTED));
    else
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

    if (k <= 0) {
        IGRAPH_ERROR("The Mycielski graph order must be a positive integer.", IGRAPH_EINVAL);
    }
    if (k == 1) {
        IGRAPH_CHECK(igraph_empty(graph, 1, IGRAPH_UNDIRECTED));
        return IGRAPH_SUCCESS;
    }

    // create a path 0---1
    IGRAPH_CHECK(igraph_ring(&g, 2, IGRAPH_UNDIRECTED, 0, 0));
    IGRAPH_FINALLY(igraph_destroy, &g);

    igraph_mycielskian(graph, &g, k - 2);

    igraph_destroy(&g);

    return IGRAPH_SUCCESS;
}
