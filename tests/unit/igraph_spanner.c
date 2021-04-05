/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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


#include <igraph.h>
#include "test_utilities.inc"
#include <stdlib.h>

igraph_bool_t _test_spanner (igraph_t *graph, igraph_t *spanner, double stretch, igraph_vector_t *weights, igraph_vector_t *spanner_weights) {
    /* check the spanner and compare it to the original graph in several topics:
        1) Compare the number of nodes.
        2) If the weights of the edges are equal between the graphs (doesn't work on multigraph).
        3) the stretch factor of the spanner.
        returns True or false
    */
   long int spanner_no_of_nodes = igraph_vcount(spanner);
   long int no_of_nodes = igraph_vcount(graph);
   long int no_of_edges = igraph_ecount(graph);
   long int no_of_edges_spanner = igraph_ecount(spanner);
   igraph_integer_t from, to, edge_spanner, edge;
   igraph_eit_t graph_edges;
   igraph_matrix_t res_spanner, res_graph;

    // compare number of nodes
    IGRAPH_ASSERT(spanner_no_of_nodes == no_of_nodes);
    // compare number of edges
    IGRAPH_ASSERT(no_of_edges > no_of_edges_spanner);

    igraph_matrix_init(&res_spanner, 0, 0);
    igraph_matrix_init(&res_graph, 0, 0);

    igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &graph_edges);
    // compare the edges and their weight
    for (IGRAPH_EIT_RESET(graph_edges); !IGRAPH_EIT_END(graph_edges); IGRAPH_EIT_NEXT(graph_edges)) {
        edge = IGRAPH_EIT_GET(graph_edges);
        igraph_edge(graph, edge, &from, &to);
        igraph_get_eid(spanner, &edge_spanner, from, to, IGRAPH_UNDIRECTED, 0);
        if (edge != -1) {
            double weight_graph = VECTOR(*weights)[edge];
            double weight_spanner = VECTOR(*spanner_weights)[edge_spanner];
            IGRAPH_ASSERT(weight_graph == weight_graph);
        }
    }

    // Validate the stretch factor
    igraph_shortest_paths_dijkstra(graph, &res_graph, igraph_vss_all(), igraph_vss_all(), weights, IGRAPH_ALL);
    igraph_shortest_paths_dijkstra(spanner, &res_spanner, igraph_vss_all(), igraph_vss_all(), spanner_weights, IGRAPH_ALL);
    for (int x = 0; x < no_of_nodes; x++) {
        for (int y = 0; y < no_of_nodes; y++) {
            if (x == y) {
                continue;
            }
            IGRAPH_ASSERT(MATRIX(res_spanner, x, y) < MATRIX(res_graph, x, y) * stretch);
        }
    }

    return IGRAPH_SUCCESS;
}

int main () {
    igraph_t graph, spanner;
    igraph_vector_t weights, spanner_weight;
    igraph_bool_t check_spanner;
    igraph_rng_t rng;
    long int no_of_edges, no_of_edges_spanner;
    IGRAPH_CHECK(igraph_rng_init(&rng, &igraph_rngtype_mt19937));
    IGRAPH_CHECK(igraph_rng_seed(&rng, time(0)));

    // trevial spanner with stretch of one
    igraph_full(&graph, 20, IGRAPH_UNDIRECTED, 0);
    unsigned long seed = (unsigned long) time(0);
    igraph_spanner(&graph, &spanner, 1, NULL, NULL, &seed);
    no_of_edges = igraph_ecount(&graph);
    no_of_edges_spanner = igraph_ecount(&spanner);
    IGRAPH_ASSERT(no_of_edges_spanner == no_of_edges);
    igraph_destroy(&spanner);

    //Test spanner for random weighted complete graph
    no_of_edges = igraph_ecount(&graph);
    igraph_vector_init(&weights, no_of_edges);
    for (int i = 0; i < no_of_edges; i++) {
        double generated_number = igraph_rng_get_unif(&rng, 1, 100);
        VECTOR(weights)[i] = generated_number;
    }
    igraph_spanner(&graph, &spanner, 10, &weights, &spanner_weight, NULL);
    check_spanner = _test_spanner(&graph, &spanner, 10, &weights, &spanner_weight);
    igraph_vector_destroy(&spanner_weight);
    igraph_vector_destroy(&weights);
    igraph_destroy(&spanner);

    // Test spanner for unweighted complete graph
    igraph_spanner(&graph, &spanner, 5, NULL, NULL, &seed);
    int spanner_no_of_edges = igraph_ecount(&spanner);
    no_of_edges = igraph_ecount(&graph);
    igraph_vector_init(&weights, no_of_edges);
    igraph_vector_fill(&weights, 1);
    igraph_vector_init(&spanner_weight, spanner_no_of_edges);
    igraph_vector_fill(&spanner_weight, 1);
    _test_spanner(&graph, &spanner, 5, &weights, &spanner_weight);
    igraph_vector_destroy(&spanner_weight);
    igraph_vector_destroy(&weights);
    igraph_destroy(&spanner);
    igraph_destroy(&graph);

    // Random erdos renyi graph
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, 200, 0.9, IGRAPH_UNDIRECTED, 0);
    no_of_edges = igraph_ecount(&graph);
    igraph_vector_init(&weights, no_of_edges);
    for (int i = 0; i < no_of_edges; i++) {
        double generated_number = igraph_rng_get_unif(&rng, 1, 100);
        VECTOR(weights)[i] = generated_number;
    }
    igraph_spanner(&graph, &spanner, 7, &weights, &spanner_weight, NULL);
    _test_spanner(&graph, &spanner, 5, &weights, &spanner_weight);
    igraph_vector_destroy(&spanner_weight);
    igraph_vector_destroy(&weights);
    igraph_destroy(&spanner);
    igraph_destroy(&graph);

    // Singlton graph
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_spanner(&graph, &spanner, 2, NULL, NULL, NULL);
    int no_of_nodes = igraph_vcount(&spanner);
    IGRAPH_ASSERT(no_of_nodes == 1);
    igraph_destroy(&spanner);
    igraph_destroy(&graph);

    igraph_rng_destroy(&rng);
    return IGRAPH_SUCCESS;
}