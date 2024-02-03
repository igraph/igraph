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

inline void check_err(igraph_error_t err) {
    if (err != IGRAPH_SUCCESS)
        abort();
}

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    igraph_t graph;
    igraph_vector_int_t edges;

    igraph_set_warning_handler(igraph_warning_handler_ignore);

    if (Size % 2 == 0 || Size > 512+1 || Size < 1) {
        return 0;
    }

    check_err(igraph_vector_int_init(&edges, Size-1));
    for (size_t i=0; i < Size-1; ++i) {
        VECTOR(edges)[i] = Data[i+1];
    }

    /* Undirected */
    if (igraph_create(&graph, &edges, Data[0], IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS) {
        igraph_vector_int_list_t ivl1, ivl2;
        igraph_vector_int_t iv1, iv2, iv3, iv4, iv5;
        igraph_matrix_t m;
        igraph_integer_t i, i2;
        igraph_real_t r;
        igraph_t g;

        check_err(igraph_vector_int_list_init(&ivl1, 0));
        check_err(igraph_vector_int_list_init(&ivl2, 0));
        check_err(igraph_vector_int_init(&iv1, 0));
        check_err(igraph_vector_int_init(&iv2, 0));
        check_err(igraph_vector_int_init(&iv3, 0));
        check_err(igraph_vector_int_init(&iv4, 0));
        check_err(igraph_vector_int_init(&iv5, 0));
        check_err(igraph_matrix_init(&m, 0, 0));

        igraph_biconnected_components(&graph, &i, NULL, &ivl1, &ivl2, &iv1);
        igraph_maximum_cardinality_search(&graph, &iv1, &iv2);
        igraph_coreness(&graph, &iv1, IGRAPH_ALL);
        igraph_girth(&graph, &r, &iv1);
        igraph_bridges(&graph, &iv1);
        igraph_assortativity_degree(&graph, &r, IGRAPH_UNDIRECTED);
        igraph_count_multiple(&graph, &iv1, igraph_ess_all(IGRAPH_EDGEORDER_FROM));

        // These algorithms require a starting vertex,
        // so we require the graph to have at least one vertex.
        if (igraph_vcount(&graph) >=1) {
            igraph_distances(&graph, &m, igraph_vss_1(0), igraph_vss_all(), IGRAPH_ALL);
            igraph_get_shortest_paths(&graph, &ivl1, &ivl2, 0, igraph_vss_all(), IGRAPH_ALL, &iv1, &iv2);
            igraph_get_all_shortest_paths(&graph, &ivl1, &ivl2, &iv1, 0, igraph_vss_all(), IGRAPH_ALL);
            igraph_pseudo_diameter(&graph, &r, 0, &i, &i2, false, true);
            igraph_bfs(&graph, 0, NULL, IGRAPH_ALL, true, NULL, &iv1, &iv2, &iv3, &iv4, NULL, &iv5, NULL, NULL);
            igraph_dfs(&graph, 0, IGRAPH_ALL, true, &iv1, &iv2, &iv3, &iv4, NULL, NULL, NULL);
        }

        igraph_connected_components(&graph, &iv1, &iv2, &i, IGRAPH_WEAK);
        igraph_minimum_spanning_tree_unweighted(&graph, &g);
        // Only when there is precisely one connected component:
        if (i == 1 && igraph_vcount(&g) >= 2) {
            igraph_to_prufer(&g, &iv1);

            igraph_t t;
            igraph_from_prufer(&t, &iv1);
            igraph_destroy(&t);
        }
        igraph_destroy(&g);

        igraph_simplify(&graph, true, true, NULL);

        igraph_matrix_destroy(&m);
        igraph_vector_int_destroy(&iv5);
        igraph_vector_int_destroy(&iv4);
        igraph_vector_int_destroy(&iv3);
        igraph_vector_int_destroy(&iv2);
        igraph_vector_int_destroy(&iv1);
        igraph_vector_int_list_destroy(&ivl2);
        igraph_vector_int_list_destroy(&ivl1);

        igraph_destroy(&graph);
    }

    igraph_vector_int_destroy(&edges);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;  // Non-zero return values are reserved for future use.
}
