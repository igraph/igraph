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

void test_spanner(igraph_t *graph, igraph_t *spanner, double stretch, igraph_vector_t *weights, igraph_vector_t *spanner_weights) {
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

    // compare the edges and their weight
    igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &graph_edges);
    for (IGRAPH_EIT_RESET(graph_edges); !IGRAPH_EIT_END(graph_edges); IGRAPH_EIT_NEXT(graph_edges)) {
        edge = IGRAPH_EIT_GET(graph_edges);
        igraph_edge(graph, edge, &from, &to);
        igraph_get_eid(spanner, &edge_spanner, from, to, IGRAPH_UNDIRECTED, 0);
        if (edge_spanner != -1) {
            double weight_graph = weights ? VECTOR(*weights)[edge] : 1;
            double weight_spanner = VECTOR(*spanner_weights)[edge_spanner];
            IGRAPH_ASSERT(weight_graph == weight_spanner);
        }
    }

    // Validate the stretch factor
    igraph_matrix_init(&res_spanner, 0, 0);
    igraph_matrix_init(&res_graph, 0, 0);
    igraph_shortest_paths_dijkstra(graph, &res_graph, igraph_vss_all(), igraph_vss_all(), weights, IGRAPH_ALL);
    igraph_shortest_paths_dijkstra(spanner, &res_spanner, igraph_vss_all(), igraph_vss_all(), spanner_weights, IGRAPH_ALL);
    for (int x = 0; x < no_of_nodes; x++) {
        for (int y = 0; y < no_of_nodes; y++) {
            if (x == y) {
                continue;
            }
            IGRAPH_ASSERT(MATRIX(res_spanner, x, y) <= MATRIX(res_graph, x, y) * stretch);
        }
    }
    igraph_matrix_destroy(&res_graph);
    igraph_matrix_destroy(&res_spanner);
}

int main () {
    igraph_t graph, spanner;
    igraph_vector_t weights, spanner_weight;
    long int no_of_nodes, no_of_edges, no_of_edges_spanner;

    /* Seed the RNG to make the test output predictable */
    igraph_rng_seed(igraph_rng_default(), 42);

    /* Trivial spanner with stretch of one */
    printf("Trivial case with stretch of one\n");
    igraph_full(&graph, 20, IGRAPH_UNDIRECTED, 0);
    igraph_spanner(&graph, &spanner, 1, NULL, NULL);
    no_of_edges = igraph_ecount(&graph);
    no_of_edges_spanner = igraph_ecount(&spanner);
    IGRAPH_ASSERT(no_of_edges_spanner == no_of_edges);
    igraph_destroy(&spanner);

    /* Test spanner for random weighted complete graph */
    printf("Weighted complete graph\n");
    no_of_edges = igraph_ecount(&graph);
    igraph_vector_init(&weights, no_of_edges);
    for (int i = 0; i < no_of_edges; i++) {
        double generated_number = igraph_rng_get_unif(igraph_rng_default(), 1, 100);
        VECTOR(weights)[i] = generated_number;
    }
    igraph_spanner(&graph, &spanner, 10, &weights, &spanner_weight);
    test_spanner(&graph, &spanner, 10, &weights, &spanner_weight);
    igraph_vector_destroy(&spanner_weight);
    igraph_vector_destroy(&weights);
    igraph_destroy(&spanner);

    /* Test spanner for unweighted complete graph */
    printf("Unweighted complete graph\n");
    igraph_spanner(&graph, &spanner, 5, NULL, NULL);
    int spanner_no_of_edges = igraph_ecount(&spanner);
    no_of_edges = igraph_ecount(&graph);
    igraph_vector_init(&weights, no_of_edges);
    igraph_vector_fill(&weights, 1);
    igraph_vector_init(&spanner_weight, spanner_no_of_edges);
    igraph_vector_fill(&spanner_weight, 1);
    test_spanner(&graph, &spanner, 5, &weights, &spanner_weight);
    igraph_vector_destroy(&spanner_weight);
    igraph_vector_destroy(&weights);
    igraph_destroy(&spanner);
    igraph_destroy(&graph);

    /* Random Erdos-Renyi graph */
    printf("Random Erdos-Renyi graph\n");
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, 200, 0.25, IGRAPH_UNDIRECTED, 0);
    no_of_edges = igraph_ecount(&graph);
    igraph_vector_init(&weights, no_of_edges);
    for (int i = 0; i < no_of_edges; i++) {
        double generated_number = igraph_rng_get_unif(igraph_rng_default(), 1, 100);
        VECTOR(weights)[i] = generated_number;
    }
    igraph_spanner(&graph, &spanner, 7, &weights, &spanner_weight);
    test_spanner(&graph, &spanner, 5, &weights, &spanner_weight);
    igraph_vector_destroy(&spanner_weight);
    igraph_vector_destroy(&weights);
    igraph_destroy(&spanner);
    igraph_destroy(&graph);

    /* Geometric random graph */
    printf("Geometric random graph, unweighted, but requesting output weights\n");
    igraph_grg_game(&graph, 100, 0.2, /* torus = */ 0, 0, 0);
    igraph_spanner(&graph, &spanner, 7, 0, &spanner_weight);
    test_spanner(&graph, &spanner, 5, 0, &spanner_weight);
    igraph_vector_destroy(&spanner_weight);
    igraph_destroy(&spanner);
    igraph_destroy(&graph);

    /* Singleton graph */
    printf("Singleton graph\n");
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_spanner(&graph, &spanner, 2, NULL, NULL);
    no_of_nodes = igraph_vcount(&spanner);
    IGRAPH_ASSERT(no_of_nodes == 1);
    igraph_destroy(&spanner);
    igraph_destroy(&graph);

    /* Null graph */
    printf("Null graph\n");
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_spanner(&graph, &spanner, 2, NULL, NULL);
    no_of_nodes = igraph_vcount(&spanner);
    IGRAPH_ASSERT(no_of_nodes == 0);
    igraph_destroy(&spanner);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    /* Error conditions */
    igraph_set_error_handler(igraph_error_handler_ignore);

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, 200, 0.9, IGRAPH_UNDIRECTED, 0);
    no_of_edges = igraph_ecount(&graph);
    igraph_vector_init(&weights, no_of_edges);
    for (int i = 0; i < no_of_edges; i++) {
        double generated_number = igraph_rng_get_unif(igraph_rng_default(), 1, 100);
        VECTOR(weights)[i] = generated_number;
    }

    printf("Negative weight\n");
    VECTOR(weights)[10] = -42;
    IGRAPH_ASSERT(igraph_spanner(&graph, &spanner, 7, &weights, &spanner_weight) == IGRAPH_EINVAL);
    VECTOR(weights)[10] = 42;

    printf("NaN weight\n");
    VECTOR(weights)[10] = IGRAPH_NAN;
    IGRAPH_ASSERT(igraph_spanner(&graph, &spanner, 7, &weights, &spanner_weight) == IGRAPH_EINVAL);
    VECTOR(weights)[10] = 42;

    printf("Invalid spanning factor\n");
    IGRAPH_ASSERT(igraph_spanner(&graph, &spanner, 0.5, &weights, &spanner_weight) == IGRAPH_EINVAL);

    printf("Invalid weight vector length\n");
    igraph_vector_resize(&weights, no_of_edges - 1);
    IGRAPH_ASSERT(igraph_spanner(&graph, &spanner, 7, &weights, &spanner_weight) == IGRAPH_EINVAL);

    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    return IGRAPH_SUCCESS;
}
