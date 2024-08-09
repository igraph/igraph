/*
   IGraph library.
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

    igraph_rng_seed(igraph_rng_default(), 137);

    /* Directed */
    if (igraph_create(&graph, &edges, Data[0], IGRAPH_DIRECTED) == IGRAPH_SUCCESS) {
        igraph_vector_int_list_t ivl1;
        igraph_vector_t v1, v2;
        igraph_vector_int_t iv1, iv2;

        igraph_vector_int_list_init(&ivl1, 0);
        igraph_vector_init(&v1, 0);
        igraph_vector_init(&v2, 0);
        igraph_vector_int_init(&iv1, 0);
        igraph_vector_int_init(&iv2, 0);

        igraph_minimum_cycle_basis(&graph, &ivl1, -1, true, true, NULL);

        igraph_vector_resize(&v2, 3);
        igraph_vector_null(&v2);
        igraph_motifs_randesu(&graph, &v1, 3, &v2);

        igraph_list_triangles(&graph, &iv1);

        igraph_count_reachable(&graph, &iv1, IGRAPH_OUT);

        if (igraph_vcount(&graph) >= 2) {
            igraph_get_all_simple_paths(&graph,&iv1,  0, igraph_vss_1(1), 5, IGRAPH_ALL);
        }

        if (igraph_vcount(&graph) >= 1) {
            igraph_t subg;

            igraph_random_walk(&graph, NULL, &iv1, &iv2, 0, IGRAPH_ALL, igraph_ecount(&graph), IGRAPH_RANDOM_WALK_STUCK_RETURN);

            igraph_induced_subgraph(&graph, &subg, igraph_vss_vector(&iv1), IGRAPH_SUBGRAPH_AUTO);
            igraph_destroy(&subg);

            igraph_subgraph_from_edges(&graph, &subg, igraph_ess_vector(&iv2), true);
            igraph_destroy(&subg);

            igraph_reverse_edges(&graph, igraph_ess_vector(&iv2));
        }

        igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_COLLAPSE, NULL);

        igraph_vector_resize(&v2, 4);
        igraph_vector_null(&v2);
        igraph_motifs_randesu(&graph, &v1, 4, &v2);

        if (igraph_vcount(&graph) >= 1) {
            igraph_random_walk(&graph, NULL, &iv1, &iv2, 0, IGRAPH_ALL, igraph_ecount(&graph), IGRAPH_RANDOM_WALK_STUCK_RETURN);
        }

        igraph_vector_int_destroy(&iv2);
        igraph_vector_int_destroy(&iv1);
        igraph_vector_destroy(&v2);
        igraph_vector_destroy(&v1);
        igraph_vector_int_list_destroy(&ivl1);

        igraph_destroy(&graph);
    }

    igraph_vector_int_destroy(&edges);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;  // Non-zero return values are reserved for future use.
}
