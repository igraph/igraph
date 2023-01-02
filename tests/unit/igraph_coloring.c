#include <igraph.h>

#include "test_utilities.h"

void verify_coloring(igraph_t *graph, igraph_vector_int_t *colors) {
    /* Verify that the colouring is valid, i.e. no two adjacent vertices have the same colour. */
    igraph_integer_t i;
    /* Store the edge count to avoid the overhead from igraph_ecount in the for loop. */
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    for (i = 0; i < no_of_edges; ++i) {
        if (IGRAPH_FROM(graph, i) == IGRAPH_TO(graph, i)) {
            continue;
        }
        if ( VECTOR(*colors)[ IGRAPH_FROM(graph, i) ] == VECTOR(*colors)[ IGRAPH_TO(graph, i) ]  ) {
            IGRAPH_FATALF("Inconsistent coloring! Vertices %" IGRAPH_PRId " and %" IGRAPH_PRId " are adjacent but have the same color.\n",
                          IGRAPH_FROM(graph, i), IGRAPH_TO(graph, i));
        }
    }
    igraph_integer_t vertex_count = igraph_vcount(graph);
    /*verify all assigned colors are valid*/
    for (i = 0 ; i < vertex_count ; i++) {
        igraph_integer_t color = VECTOR(*colors)[i];
        if (color < 0 || color >= vertex_count) {
            IGRAPH_FATALF("The vertex %" IGRAPH_PRId " has invalid color %" IGRAPH_PRId "\n", i, color);
        }
    }
}

void run_tests(igraph_t *graph, igraph_vector_int_t *colors) {
    igraph_vertex_coloring_greedy(graph, colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS);
    verify_coloring(graph, colors);

    igraph_vector_int_fill(colors, 0);
    igraph_vertex_coloring_greedy(graph, colors, IGRAPH_COLORING_GREEDY_DSATUR);
    verify_coloring(graph, colors);
}

void test_empty_graph(void) {
    igraph_t graph;
    igraph_vector_int_t colors;
    igraph_vector_int_init(&colors, 0);

    //run it on empty graph
    igraph_empty(&graph, 0, IGRAPH_DIRECTED);
    run_tests(&graph, &colors);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&colors);
}

void test_graph_1_vertex(void) {
    igraph_t graph;
    igraph_vector_int_t colors;
    igraph_vector_int_init(&colors, 0);

    //run it on empty graph
    igraph_empty(&graph, 0, IGRAPH_DIRECTED);
    run_tests(&graph, &colors);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&colors);
}

void test_graph_large(void) {
    igraph_t graph;
    igraph_vector_int_t colors;
    igraph_vector_int_init(&colors, 0);

    //simple large undirected graphx
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 1000, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    run_tests(&graph, &colors);

    //graph with loops and parallel edges
    igraph_rewire_edges(&graph, 1.0, true, true);
    run_tests(&graph, &colors);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&colors);
}

int main(void) {
    igraph_rng_seed(igraph_rng_default(), 42);

    test_empty_graph();
    test_graph_1_vertex();
    test_graph_large();

    VERIFY_FINALLY_STACK();

    return 0;
}

