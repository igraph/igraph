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
    igraph_vector_t weights;

    igraph_set_warning_handler(igraph_warning_handler_ignore);

    if (Size % 3 == 0 || Size % 3 == 2 || Size > (3 * 256) + 1 || Size < 1) {
        return 0;
    }

    igraph_vector_int_init(&edges, ((Size-1) / 3) * 2);
    igraph_vector_init(&weights, (Size-1) / 3);
    for (size_t i=0; i < ((Size-1) / 3); ++i) {
        VECTOR(edges)[i * 2] = Data[i * 3 + 1];
        VECTOR(edges)[i * 2 + 1] = Data[i * 3 + 2];
        // We keep the weights strictly positive, as this is required by some algorithms.
        VECTOR(weights)[i] = ((double) Data[i * 3 + 3] + 1.0) / 105.0;
    }

    // Turn on attribute handling.
    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_rng_seed(igraph_rng_default(), 42);

    if (igraph_create(&graph, &edges, Data[0], IGRAPH_DIRECTED) == IGRAPH_SUCCESS) {
        igraph_real_t r, r2;
        igraph_vector_int_list_t ivl1;
        igraph_vector_t v1, v2;
        igraph_vector_int_t iv1, iv2;
        igraph_matrix_t m;

        igraph_vector_int_list_init(&ivl1, 0);
        igraph_vector_init(&v1, 0);
        igraph_vector_init(&v2, 0);
        igraph_vector_int_init(&iv1, 0);
        igraph_vector_int_init(&iv2, 0);
        igraph_matrix_init(&m, 0, 0);

        /* Directed */

        igraph_get_adjacency(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, &weights, IGRAPH_NO_LOOPS);
        igraph_get_laplacian(&graph, &m, IGRAPH_OUT, IGRAPH_LAPLACIAN_LEFT, &weights);
        igraph_get_stochastic(&graph, &m, false, &weights);
        igraph_joint_degree_matrix(&graph, &weights, &m, -1, -1);
        igraph_joint_degree_distribution(&graph, &weights, &m, IGRAPH_IN, IGRAPH_OUT, true, true, -1, -1);
        igraph_degree_correlation_vector(&graph, &weights, &v1, IGRAPH_OUT, IGRAPH_IN, true);
        igraph_avg_nearest_neighbor_degree(&graph, igraph_vss_all(), IGRAPH_OUT, IGRAPH_IN, &v1, &v2, &weights);

        igraph_pseudo_diameter(&graph, &weights, &r, -1, NULL, NULL, true, true);
        igraph_diameter(&graph, &weights, &r, NULL, NULL, NULL, NULL, true, false);
        igraph_average_path_length(&graph, &weights, &r, &r2, true, true);
        igraph_radius(&graph, &weights, &r, IGRAPH_OUT);

        igraph_feedback_arc_set(&graph, &iv1, &weights, IGRAPH_FAS_APPROX_EADES);

        if (igraph_vcount(&graph) >= 1) {
            igraph_random_walk(&graph, &weights, &iv1, &iv2, 0, IGRAPH_ALL, igraph_ecount(&graph), IGRAPH_RANDOM_WALK_STUCK_RETURN);
            igraph_distances_dijkstra(&graph, &m, igraph_vss_1(0), igraph_vss_all(), &weights, IGRAPH_OUT);
            igraph_distances_bellman_ford(&graph, &m, igraph_vss_1(0), igraph_vss_all(), &weights, IGRAPH_OUT);
            igraph_widest_path_widths_dijkstra(&graph, &m, igraph_vss_1(0), igraph_vss_all(), &weights, IGRAPH_IN);
        }

        if (igraph_vcount(&graph) >= 2) {
            igraph_get_k_shortest_paths(&graph, &weights, &ivl1, NULL, 5, 0, 1, IGRAPH_IN);
        }

        /* Undirected */

        {
            igraph_attribute_combination_t comb;

            SETEANV(&graph, "weight", &weights);
            igraph_attribute_combination(&comb,
                                         "weight", IGRAPH_ATTRIBUTE_COMBINE_SUM,
                                         IGRAPH_NO_MORE_ATTRIBUTES);
            igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_COLLAPSE, &comb);
            igraph_attribute_combination_destroy(&comb);
            EANV(&graph, "weight", &weights);
        }

        igraph_get_adjacency(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, &weights, IGRAPH_LOOPS_ONCE);
        igraph_get_laplacian(&graph, &m, IGRAPH_OUT, IGRAPH_LAPLACIAN_UNNORMALIZED, &weights);
        igraph_get_stochastic(&graph, &m, true, &weights);
        igraph_joint_degree_matrix(&graph, &weights, &m, -1, -1);
        igraph_joint_degree_distribution(&graph, &weights, &m, IGRAPH_ALL, IGRAPH_ALL, false, true, -1, -1);
        igraph_degree_correlation_vector(&graph, &weights, &v1, IGRAPH_ALL, IGRAPH_ALL, false);
        igraph_avg_nearest_neighbor_degree(&graph, igraph_vss_all(), IGRAPH_ALL, IGRAPH_ALL, &v1, &v2, &weights);

        igraph_pseudo_diameter(&graph, &weights, &r, -1, NULL, NULL, false, false);
        igraph_diameter(&graph, &weights, &r, NULL, NULL, NULL, NULL, false, true);
        igraph_average_path_length(&graph, &weights, &r, &r, true, true);
        igraph_radius(&graph, &weights, &r, IGRAPH_OUT);

        if (igraph_vcount(&graph) >= 1) {
            igraph_random_walk(&graph, &weights, &iv1, &iv2, 0, IGRAPH_ALL, igraph_ecount(&graph), IGRAPH_RANDOM_WALK_STUCK_RETURN);
            igraph_distances_dijkstra(&graph, &m, igraph_vss_1(0), igraph_vss_all(), &weights, IGRAPH_ALL);
            igraph_distances_bellman_ford(&graph, &m, igraph_vss_1(0), igraph_vss_all(), &weights, IGRAPH_ALL);
            igraph_widest_path_widths_dijkstra(&graph, &m, igraph_vss_1(0), igraph_vss_all(), &weights, IGRAPH_ALL);
        }

        igraph_minimum_spanning_tree(&graph, &iv1, &weights, IGRAPH_MST_PRIM);
        igraph_minimum_spanning_tree(&graph, &iv1, &weights, IGRAPH_MST_KRUSKAL);

        igraph_spanner(&graph, &iv1, 2.34, &weights);

        igraph_matrix_destroy(&m);
        igraph_vector_int_destroy(&iv2);
        igraph_vector_int_destroy(&iv1);
        igraph_vector_destroy(&v2);
        igraph_vector_destroy(&v1);
        igraph_vector_int_list_destroy(&ivl1);

        igraph_destroy(&graph);
    }

    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&weights);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;  // Non-zero return values are reserved for future use.
}
