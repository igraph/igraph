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
#include "test_utilities.h"
#include <stdlib.h>

void test_spanner(igraph_t *graph, igraph_vector_int_t *spanner, double stretch, igraph_vector_t *weights) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_t spanner_graph;
    igraph_vector_t spanner_weights;
    igraph_matrix_t res_spanner, res_graph;

    /* create the spanner graph with igraph_subgraph_from_edges() as recommended in the docs,
       then compare it with the original graph and validate the stretch factor */
    if (weights) {
        igraph_cattribute_EAN_setv(graph, "weight", weights);
    }
    igraph_subgraph_from_edges(graph, &spanner_graph, igraph_ess_vector(spanner), 0);
    igraph_vector_init(&spanner_weights, igraph_vector_int_size(spanner));
    if (weights){
        igraph_cattribute_EANV(&spanner_graph, "weight", igraph_ess_all(IGRAPH_EDGEORDER_ID), &spanner_weights);
    } else {
        igraph_vector_fill(&spanner_weights, 1);
    }

    // compare number of nodes
    IGRAPH_ASSERT(igraph_vcount(&spanner_graph) == no_of_nodes);

    // compare number of edges. We expect the number of edges to decrease in
    // all cases but the trivial ones
    if (stretch > 1 && no_of_edges > 1) {
        IGRAPH_ASSERT(igraph_ecount(&spanner_graph) < no_of_edges);
    }

    // print the number of nodes, the number of original edges and the new edge
    // count. This is not validated (there is no expected output) but it helps
    // to gauge whether the algorithm is not simply keeping most of the edges
    printf(
        "%" IGRAPH_PRId ", %" IGRAPH_PRId " --> %" IGRAPH_PRId "\n",
        no_of_nodes, no_of_edges, igraph_ecount(&spanner_graph)
    );

    // Validate the stretch factor
    igraph_matrix_init(&res_spanner, 0, 0);
    igraph_matrix_init(&res_graph, 0, 0);
    igraph_distances_dijkstra(graph, &res_graph, igraph_vss_all(), igraph_vss_all(), weights, IGRAPH_ALL);
    igraph_distances_dijkstra(&spanner_graph, &res_spanner, igraph_vss_all(), igraph_vss_all(), &spanner_weights, IGRAPH_ALL);
    for (igraph_integer_t x = 0; x < no_of_nodes; x++) {
        for (igraph_integer_t y = 0; y < no_of_nodes; y++) {
            if (x == y) {
                continue;
            }
            IGRAPH_ASSERT(MATRIX(res_spanner, x, y) <= MATRIX(res_graph, x, y) * stretch);
        }
    }
    igraph_matrix_destroy(&res_graph);
    igraph_matrix_destroy(&res_spanner);

    // Clean up
    igraph_vector_destroy(&spanner_weights);
    igraph_destroy(&spanner_graph);
}

int main(void) {
    igraph_t graph;
    igraph_vector_t weights;
    igraph_vector_int_t spanner;
    igraph_integer_t no_of_edges;

    /* Initialize attribute handler; we will use edge attributes in test_spanner() */
    igraph_set_attribute_table(&igraph_cattribute_table);

    /* Seed the RNG to make the test output predictable */
    igraph_rng_seed(igraph_rng_default(), 42);

    /* Create the output vector -- this will be re-used several times */
    igraph_vector_int_init(&spanner, 0);

    /* Trivial spanner with stretch of one */
    printf("Complete graph with stretch of one\n");
    igraph_full(&graph, 20, IGRAPH_UNDIRECTED, 0);
    igraph_spanner(&graph, &spanner, 1, NULL);
    test_spanner(&graph, &spanner, 1, NULL);
    igraph_destroy(&graph);

    /* Test spanner for random weighted complete graph */
    printf("Weighted complete graph with random weights and stretch = 10\n");
    igraph_full(&graph, 20, IGRAPH_UNDIRECTED, 0);
    no_of_edges = igraph_ecount(&graph);
    igraph_vector_init(&weights, no_of_edges);
    for (int i = 0; i < no_of_edges; i++) {
        double generated_number = igraph_rng_get_unif(igraph_rng_default(), 1, 100);
        VECTOR(weights)[i] = generated_number;
    }
    igraph_spanner(&graph, &spanner, 10, &weights);
    test_spanner(&graph, &spanner, 10, &weights);
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    /* Test spanner for unweighted complete graph */
    printf("Unweighted complete graph with stretch = 5\n");
    igraph_full(&graph, 20, IGRAPH_UNDIRECTED, 0);
    igraph_spanner(&graph, &spanner, 5, NULL);
    test_spanner(&graph, &spanner, 5, NULL);
    igraph_destroy(&graph);

    /* Random Erdos-Renyi graph */
    printf("Random Erdos-Renyi graph\n");
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, 200, 0.25, IGRAPH_UNDIRECTED, 0);
    no_of_edges = igraph_ecount(&graph);
    igraph_vector_init(&weights, no_of_edges);
    for (igraph_integer_t i = 0; i < no_of_edges; i++) {
        double generated_number = igraph_rng_get_unif(igraph_rng_default(), 1, 100);
        VECTOR(weights)[i] = generated_number;
    }
    igraph_spanner(&graph, &spanner, 7, &weights);
    test_spanner(&graph, &spanner, 7, &weights);
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    /* Geometric random graph */
    printf("Geometric random graph, unweighted\n");
    igraph_grg_game(&graph, 100, 0.2, /* torus = */ 0, 0, 0);
    igraph_spanner(&graph, &spanner, 7, 0);
    test_spanner(&graph, &spanner, 7, 0);
    igraph_destroy(&graph);

    /* Singleton graph */
    printf("Singleton graph\n");
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_spanner(&graph, &spanner, 2, 0);
    test_spanner(&graph, &spanner, 2, 0);
    igraph_destroy(&graph);

    /* Null graph */
    printf("Null graph\n");
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_spanner(&graph, &spanner, 2, 0);
    test_spanner(&graph, &spanner, 2, 0);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    /* Error conditions */
    igraph_set_error_handler(igraph_error_handler_ignore);

    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, 200, 0.9, IGRAPH_UNDIRECTED, 0);
    no_of_edges = igraph_ecount(&graph);
    igraph_vector_init(&weights, no_of_edges);
    for (igraph_integer_t i = 0; i < no_of_edges; i++) {
        double generated_number = igraph_rng_get_unif(igraph_rng_default(), 1, 100);
        VECTOR(weights)[i] = generated_number;
    }

    printf("Negative weight\n");
    VECTOR(weights)[10] = -42;
    IGRAPH_ASSERT(igraph_spanner(&graph, &spanner, 7, &weights) == IGRAPH_EINVAL);
    VECTOR(weights)[10] = 42;

    printf("NaN weight\n");
    VECTOR(weights)[10] = IGRAPH_NAN;
    IGRAPH_ASSERT(igraph_spanner(&graph, &spanner, 7, &weights) == IGRAPH_EINVAL);
    VECTOR(weights)[10] = 42;

    printf("Invalid spanning factor\n");
    IGRAPH_ASSERT(igraph_spanner(&graph, &spanner, 0.5, &weights) == IGRAPH_EINVAL);

    printf("Invalid weight vector length\n");
    igraph_vector_resize(&weights, no_of_edges - 1);
    IGRAPH_ASSERT(igraph_spanner(&graph, &spanner, 7, &weights) == IGRAPH_EINVAL);

    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    igraph_vector_int_destroy(&spanner);

    return IGRAPH_SUCCESS;
}
