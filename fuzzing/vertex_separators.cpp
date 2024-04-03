/*
   IGraph library.
   Copyright (C) 2021-2024  The igraph development team

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

    /* We work with small, up-to 16-vertex graphs, as the algorithms
     * tested here can be slow. Each byte is interpreted as an edge.
     * A simple graph can have at most 120 edges, but we allow up
     * to 240, as the fuzzer usually generates multigraphs. */

    if (Size > 240) {
        return 0;
    }

    igraph_vector_int_init(&edges, 2*Size);
    size_t j = 0;
    for (size_t i=0; i < Size; ++i) {
        VECTOR(edges)[j++] = Data[i] / 16;
        VECTOR(edges)[j++] = Data[i] % 16;
    }

    if (igraph_create(&graph, &edges, 0, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS) {
        igraph_vector_int_list_t ivl;
        igraph_vector_int_t iv1, iv2;
        igraph_bool_t is_separator, is_minimal_separator;

        igraph_vector_int_list_init(&ivl, 0);
        igraph_vector_int_init(&iv1, 0);
        igraph_vector_int_init(&iv2, 0);

        igraph_all_minimal_st_separators(&graph, &ivl);

        // Check that all returned sets are indeed separators.
        for (igraph_integer_t i=0; i < igraph_vector_int_list_size(&ivl); i++) {
            igraph_is_separator(
                &graph,
                igraph_vss_vector(igraph_vector_int_list_get_ptr(&ivl, i)),
                &is_separator);
            IGRAPH_ASSERT(is_separator);
        }

        igraph_minimum_size_separators(&graph, &ivl);

        // Simplification is necessary for cohesive_blocks() and
        // enables a straightforward check for complete graphs below.
        igraph_simplify(&graph, true, true, NULL);

        const igraph_integer_t vcount = igraph_vcount(&graph);
        const igraph_integer_t ecount = igraph_ecount(&graph);

        if (ecount != vcount*(vcount-1)/2) {
            // minimum_size_separators() returns all size n-1 subsets
            // of the complete graph K_n. is_minimal_separator() does not
            // consider these to be separators. Therefore we skip complete
            // graphs. For non-complete graphs we check that all results
            // are minimal separators.
            for (igraph_integer_t i=0; i < igraph_vector_int_list_size(&ivl); i++) {
                igraph_is_minimal_separator(
                    &graph,
                    igraph_vss_vector(igraph_vector_int_list_get_ptr(&ivl, i)),
                    &is_minimal_separator);
                IGRAPH_ASSERT(is_minimal_separator);
            }
        }

        {
            igraph_t g;
            igraph_cohesive_blocks(&graph, &ivl, &iv1, &iv2, &g);
            igraph_destroy(&g);
        }

        igraph_vector_int_destroy(&iv2);
        igraph_vector_int_destroy(&iv1);
        igraph_vector_int_list_destroy(&ivl);
        igraph_destroy(&graph);
    }

    igraph_vector_int_destroy(&edges);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;  // Non-zero return values are reserved for future use.
}
