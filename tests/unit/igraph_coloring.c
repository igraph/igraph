#include <igraph.h>

#include "test_utilities.h"

void verify_coloring(igraph_t *graph, igraph_vector_int_t *colors) {
    /* Verify that the colouring is valid, i.e. no two adjacent vertices have the same colour. */
    igraph_integer_t i;
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t vertex_count = igraph_vcount(graph);

    for (i = 0; i < no_of_edges; ++i) {
        if (IGRAPH_FROM(graph, i) == IGRAPH_TO(graph, i)) {
            continue;
        }
        if ( VECTOR(*colors)[ IGRAPH_FROM(graph, i) ] == VECTOR(*colors)[ IGRAPH_TO(graph, i) ]  ) {
            IGRAPH_FATALF("Inconsistent coloring! Vertices %" IGRAPH_PRId " and %" IGRAPH_PRId " are adjacent but have the same color.\n",
                          IGRAPH_FROM(graph, i), IGRAPH_TO(graph, i));
        }
    }

    for (i = 0; i < vertex_count; i++) {
        igraph_integer_t color = VECTOR(*colors)[i];
        if (color < 0 || color >= vertex_count) {
            IGRAPH_FATALF("The vertex %" IGRAPH_PRId " has invalid color %" IGRAPH_PRId "\n", i, color);
        }
    }
}

void print_result(
    const char* result_name, igraph_vector_int_t *colors,
    igraph_bool_t print_color_vector
) {
    printf("%s", result_name);
    if (!igraph_vector_int_empty(colors)) {
        printf("Colors used %" IGRAPH_PRId "\n", igraph_vector_int_max(colors) + 1);
    }
    if (print_color_vector) {
        print_vector_int(colors);
    }
}

void run_tests(igraph_t *graph, igraph_vector_int_t *colors, igraph_bool_t print_color_vector) {
    igraph_vector_int_fill(colors, 0);
    igraph_vertex_coloring_greedy(graph, colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS);
    verify_coloring(graph, colors);
    print_result("testing greedy coloring with COLORED_NEIGHBORS heuristic\n", colors, print_color_vector);

    igraph_vector_int_fill(colors, 0);
    igraph_vertex_coloring_greedy(graph, colors, IGRAPH_COLORING_GREEDY_DSATUR);
    verify_coloring(graph, colors);
    print_result("testing greedy coloring with DSATUR heuristic\n", colors, print_color_vector);
    printf("\n");
}

void test_empty_graph(void) {
    igraph_t graph;
    igraph_vector_int_t colors;
    igraph_vector_int_init(&colors, 0);

    /* run it on empty graph */
    igraph_empty(&graph, 0, IGRAPH_DIRECTED);
    printf("Testing empty graph\n");
    run_tests(&graph, &colors, true);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&colors);
}

void test_graph_1_vertex(void) {
    igraph_t graph;
    igraph_vector_int_t colors;
    igraph_vector_int_init(&colors, 0);

    /* run it on empty graph */
    igraph_empty(&graph, 1, IGRAPH_DIRECTED);
    printf("testing 1 vertice graph\n");
    run_tests(&graph, &colors, true);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&colors);
}

void test_graph_large(void) {
    igraph_t graph;
    igraph_vector_int_t colors;
    igraph_vector_int_init(&colors, 0);

    /* simple large undirected graph */
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 1000, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    printf("testing large simple graph\n");
    run_tests(&graph, &colors, false);

    /* graph with loops and parallel edges */
    igraph_rewire_edges(&graph, 1.0, true, true);
    printf("testing large multiedge looped graph\n");
    run_tests(&graph, &colors, false);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&colors);
}

void test_graph_small(void) {
    igraph_t graph;
    igraph_vector_int_t colors;
    igraph_vector_int_init(&colors, 0);

    igraph_small(&graph, 8, IGRAPH_UNDIRECTED,
                 0, 3, 0, 4, 1, 4, 2, 4, 1, 5, 2, 5, 3, 5, 0, 6, 1, 6, 2, 6, 3, 6, 0, 7, 1, 7, 2, 7, 4, 7, 5, 7, -1);
    printf("testing small simple graph\n");
    run_tests(&graph, &colors, true);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&colors);
}

void test_graph_lcf(void) {
    igraph_t graph;
    igraph_vector_int_t colors;
    igraph_vector_int_init(&colors, 0);

    igraph_lcf(&graph, /* n= */ 8, /* shifts= */ 2, /* repeats= */ 8, 0);
    printf("testing small lcf graph\n");
    run_tests(&graph, &colors, true);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&colors);
}

void test_isolated_vertex(void) {
    igraph_t graph;
    igraph_vector_int_t colors;
    igraph_vector_int_init(&colors, 0);

    igraph_small(&graph, 4, IGRAPH_UNDIRECTED, 0, 1, -1);
    printf("testing graph with isolated vertex\n");
    run_tests(&graph, &colors, true);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&colors);
}

void test_wheel_graph(void) {
    igraph_t graph_odd, graph_even;
    igraph_vector_int_t colors;
    igraph_vector_int_init(&colors, 0);

    igraph_wheel(&graph_odd, 11, IGRAPH_WHEEL_UNDIRECTED, 0);
    printf("testing wheel graph with odd number of vertices\n");
    run_tests(&graph_odd, &colors, true);

    igraph_wheel(&graph_even, 12, IGRAPH_WHEEL_UNDIRECTED, 0);
    printf("testing wheel graph with even number of vertices\n");
    run_tests(&graph_even, &colors, true);

    igraph_destroy(&graph_odd);
    igraph_destroy(&graph_even);
    igraph_vector_int_destroy(&colors);
}

void test_bipartite_graph(void) {
    igraph_t graph_1, graph_2;
    igraph_vector_int_t colors;
    igraph_vector_int_init(&colors, 0);

    igraph_bipartite_game_gnm(&graph_1, 0, 1000, 1000, 1000000, false, IGRAPH_ALL);
    printf("testing complete bipartite graph\n");
    run_tests(&graph_1, &colors, false);

    igraph_bipartite_game_gnm(&graph_2, 0, 1000, 1000, 1000000 - 100, false, IGRAPH_ALL);
    printf("testing large bipartite graph\n");
    run_tests(&graph_2, &colors, false);

    igraph_destroy(&graph_1);
    igraph_destroy(&graph_2);
    igraph_vector_int_destroy(&colors);
}

int main(void) {
    igraph_rng_seed(igraph_rng_default(), 42);

    test_empty_graph();
    test_graph_1_vertex();
    test_graph_large();
    test_graph_small();
    test_graph_lcf();
    test_isolated_vertex();
    test_wheel_graph();
    test_bipartite_graph();

    VERIFY_FINALLY_STACK();

    return 0;
}
