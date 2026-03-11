/*
   igraph library.
   Copyright (C) 2024  The igraph development team

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

#include <igraph.h>
#include <cstdlib>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    igraph_t graph;
    igraph_vector_int_t edges;

    igraph_set_warning_handler(igraph_warning_handler_ignore);

    if (Size % 2 == 0 || Size > 512+1 || Size < 1) {
        return 0;
    }

    igraph_vector_int_init(&edges, Size-1);
    for (size_t i=0; i < Size-1; ++i) {
        VECTOR(edges)[i] = Data[i+1];
    }

    /* Directed */
    if (igraph_create(&graph, &edges, Data[0], IGRAPH_DIRECTED) == IGRAPH_SUCCESS) {
        igraph_bool_t has_loop, has_multi, has_mutual, is_simple, is_complete;
        igraph_bool_t is_connected, is_dag, is_forest, is_tree, has_eulerian_path, has_eulerian_cycles;
        igraph_int_t ecount, vcount;
        igraph_real_t r;

        vcount = igraph_vcount(&graph);
        ecount = igraph_ecount(&graph);

        igraph_has_multiple(&graph, &has_multi);
        igraph_has_loop(&graph, &has_loop);
        igraph_has_mutual(&graph, &has_mutual, false);
        igraph_invalidate_cache(&graph);

        igraph_is_simple(&graph, &is_simple, IGRAPH_DIRECTED);
        igraph_invalidate_cache(&graph);

        IGRAPH_ASSERT((has_multi || has_loop) == !is_simple);

        igraph_is_complete(&graph, &is_complete);
        igraph_invalidate_cache(&graph);

        IGRAPH_ASSERT(!is_complete || ecount >= vcount*(vcount-1));

        igraph_is_connected(&graph, &is_connected, IGRAPH_STRONG);
        igraph_invalidate_cache(&graph);

        IGRAPH_ASSERT(!is_connected || ecount >= vcount || vcount <= 1);
        IGRAPH_ASSERT(!is_complete || is_connected || vcount == 0);

        igraph_is_dag(&graph, &is_dag);
        igraph_invalidate_cache(&graph);

        IGRAPH_ASSERT(!is_complete || !is_dag || vcount <= 1);
        IGRAPH_ASSERT(!is_connected || !is_dag || vcount == 1);

        igraph_is_forest(&graph, &is_forest, NULL, IGRAPH_OUT);
        igraph_invalidate_cache(&graph);

        IGRAPH_ASSERT(!is_forest || is_dag);
        IGRAPH_ASSERT(!is_forest || !is_connected || vcount == 1);

        igraph_is_tree(&graph, &is_tree, NULL, IGRAPH_OUT);
        igraph_invalidate_cache(&graph);

        IGRAPH_ASSERT(!is_tree || is_forest);
        IGRAPH_ASSERT(!is_tree || ecount >= vcount-1);

        igraph_is_eulerian(&graph, &has_eulerian_path, &has_eulerian_cycles);
        igraph_invalidate_cache(&graph);

        IGRAPH_ASSERT(!has_eulerian_cycles || has_eulerian_path);

        igraph_density(&graph, NULL, &r, false);

        IGRAPH_ASSERT(!is_complete || r >= 1 || vcount <= 1);

        igraph_density(&graph, NULL, &r, true);

        igraph_destroy(&graph);
    }

    igraph_vector_int_destroy(&edges);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;  // Non-zero return values are reserved for future use.
}
