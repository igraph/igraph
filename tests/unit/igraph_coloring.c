/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

/* Verify that the colouring is valid, i.e. no two adjacent vertices have the same colour. */
void verify_coloring(igraph_t *graph, const igraph_vector_int_t *colors) {
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    for (igraph_integer_t i = 0; i < no_of_edges; ++i) {
        if (IGRAPH_FROM(graph, i) == IGRAPH_TO(graph, i)) {
            continue;
        }
        if ( VECTOR(*colors)[ IGRAPH_FROM(graph, i) ] == VECTOR(*colors)[ IGRAPH_TO(graph, i) ]  ) {
            IGRAPH_FATALF("Inconsistent coloring! Vertices %" IGRAPH_PRId " and %" IGRAPH_PRId " are adjacent but have the same color.\n",
                          IGRAPH_FROM(graph, i), IGRAPH_TO(graph, i));
        }
    }

    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_integer_t color = VECTOR(*colors)[i];
        if (color < 0 || color >= no_of_nodes) {
            IGRAPH_FATALF("The vertex %" IGRAPH_PRId " has invalid color %" IGRAPH_PRId "\n", i, color);
        }
    }
}

void print_result(
    const char* result_name, const igraph_vector_int_t *colors,
    igraph_bool_t print_color_vector
) {
    printf("%s", result_name);
    if (!igraph_vector_int_empty(colors)) {
        printf("Number of colors used: %" IGRAPH_PRId "\n", igraph_vector_int_max(colors) + 1);
    }
    if (print_color_vector) {
        print_vector_int(colors);
    }
}

void run_tests(igraph_t *graph, igraph_bool_t print_color_vector) {
    igraph_vector_int_t colors;
    igraph_vector_int_init(&colors, 0);

    igraph_vector_int_fill(&colors, -1);
    igraph_vertex_coloring_greedy(graph, &colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS);
    verify_coloring(graph, &colors);
    print_result("testing greedy coloring with COLORED_NEIGHBORS heuristic\n", &colors, print_color_vector);

    igraph_vector_int_fill(&colors, -1);
    igraph_vertex_coloring_greedy(graph, &colors, IGRAPH_COLORING_GREEDY_DSATUR);
    verify_coloring(graph, &colors);
    print_result("testing greedy coloring with DSATUR heuristic\n", &colors, print_color_vector);
    printf("\n");

    igraph_vector_int_destroy(&colors);
}

void test_null_graph(void) {
    igraph_t graph;

    /* run it on the null graph */
    igraph_empty(&graph, 0, IGRAPH_DIRECTED);
    printf("Testing null graph\n");
    run_tests(&graph, true);

    igraph_destroy(&graph);
}

void test_graph_1_vertex(void) {
    igraph_t graph;

    igraph_empty(&graph, 1, IGRAPH_DIRECTED);
    printf("Testing singleton graph\n");
    run_tests(&graph, true);

    igraph_destroy(&graph);
}

void test_graph_large(void) {
    igraph_t graph;

    /* simple large undirected graph */
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 1000, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    printf("Testing large simple graph\n");
    run_tests(&graph, false);

    /* graph with loops and parallel edges */
    igraph_rewire_edges(&graph, 1.0, true, true);
    printf("Testing large graph with loops and multi-edges\n");
    run_tests(&graph, false);

    igraph_destroy(&graph);
}

// This graph has chromatic number of 3.
void test_graph_small(void) {
    igraph_t graph;

    igraph_small(&graph, 8, IGRAPH_UNDIRECTED,
                 0, 3, 0, 4, 1, 4, 2, 4, 1, 5, 2, 5, 3, 5, 0, 6, 1, 6, 2, 6, 3, 6, 0, 7, 1, 7, 2, 7, 4, 7, 5, 7,
                 -1);
    printf("Testing small simple graph\n");
    run_tests(&graph, true);

    igraph_destroy(&graph);
}

// Same as the graph above with multi-edges and self-loops added
void test_graph_small_multi(void) {
    igraph_t graph;

    igraph_small(&graph, 8, IGRAPH_UNDIRECTED,
                 0, 3, 0, 4, 1, 4, 2, 4, 1, 5, 2, 5, 3, 5, 0, 6, 1, 6, 2, 6, 3, 6, 0, 7, 1, 7, 2, 7, 4, 7, 5, 7,
                 0, 4, 0, 4, 3, 5, 3, 5, 3, 5, 0, 0, 0, 0, 1, 1,
                 -1);
    printf("Testing small multigraph\n");
    run_tests(&graph, true);

    igraph_destroy(&graph);
}

// Note: This graph has chromatic number 4.
// Currently implemented greedy heuristics don't produce an exact result.
void test_graph_lcf(void) {
    igraph_t graph;

    igraph_lcf(&graph, /* n= */ 8, /* shifts= */ 2, /* repeats= */ 8, 0);
    printf("Testing small LCF graph\n");
    run_tests(&graph, true);

    igraph_destroy(&graph);
}

// Isolated vertices must be assigned color 0.
void test_isolated_vertices(void) {
    igraph_t graph;

    igraph_small(&graph, 4, IGRAPH_UNDIRECTED, 0, 1, -1);
    printf("Testing graph with isolated vertices\n");
    run_tests(&graph, true);

    igraph_destroy(&graph);
}

// Isolated vertices must be assigned color 0.
void test_isolated_vertices_with_loops(void) {
    igraph_t graph;

    igraph_small(&graph, 3, IGRAPH_UNDIRECTED,
                 0,0, 2,2, 2,2,
                 -1);
    printf("Testing graph with isolated vertices and self-loops\n");
    run_tests(&graph, true);

    igraph_destroy(&graph);
}

// Wheel graphs on odd/even number of vertices have chromatic number of 3/4.
// DSatur is expected to find an exact minimum coloring.
void test_wheel_graph(void) {
    igraph_t graph_odd, graph_even;

    igraph_wheel(&graph_odd, 11, IGRAPH_WHEEL_UNDIRECTED, 0);
    printf("Testing wheel graph with odd number of vertices\n");
    run_tests(&graph_odd, true);

    igraph_wheel(&graph_even, 12, IGRAPH_WHEEL_UNDIRECTED, 0);
    printf("Testing wheel graph with even number of vertices\n");
    run_tests(&graph_even, true);

    igraph_destroy(&graph_odd);
    igraph_destroy(&graph_even);
}

// Bipartite graphs have a chromatic number of 2 (by definition).
// DSatur is expected to find an exact minimum coloring.
void test_bipartite_graph(void) {
    igraph_t graph_1, graph_2;

    igraph_bipartite_game_gnm(&graph_1, 0, 100, 100, 10000, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    printf("Testing complete bipartite graph\n");
    run_tests(&graph_1, false);

    igraph_bipartite_game_gnm(&graph_2, 0, 100, 100, 10000 - 100, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    printf("Testing large bipartite graph\n");
    run_tests(&graph_2, false);

    igraph_destroy(&graph_1);
    igraph_destroy(&graph_2);
}

int main(void) {
    igraph_rng_seed(igraph_rng_default(), 42);

    test_null_graph();
    test_graph_1_vertex();
    test_graph_large();
    test_graph_small();
    test_graph_small_multi();
    test_graph_lcf();
    test_isolated_vertices();
    test_isolated_vertices_with_loops();
    test_wheel_graph();
    test_bipartite_graph();

    VERIFY_FINALLY_STACK();

    return 0;
}
