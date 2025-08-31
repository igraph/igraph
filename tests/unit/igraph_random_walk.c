/*
   igraph library.
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

void random_vertex_walk(igraph_t *graph, igraph_vector_t *weights,
    igraph_int_t start, igraph_neimode_t mode,
    igraph_int_t steps, igraph_random_walk_stuck_t stuck) {
    igraph_vector_int_t vertices;

    igraph_vector_int_init(&vertices, 0);
    igraph_random_walk(graph, weights, &vertices, NULL, start, mode, steps, stuck);
    print_vector_int(&vertices);
    igraph_vector_int_destroy(&vertices);
}

void random_walk(igraph_t *graph, igraph_vector_t *weights,
    igraph_int_t start, igraph_neimode_t mode,
    igraph_int_t steps, igraph_random_walk_stuck_t stuck) {
    igraph_vector_int_t vertices;
    igraph_vector_int_t edges;

    igraph_vector_int_init(&vertices, 0);
    igraph_vector_int_init(&edges, 0);
    igraph_random_walk(graph, weights, &vertices, &edges, start, mode, steps, stuck);
    print_vector_int(&vertices);
    igraph_vector_int_destroy(&vertices);
    igraph_vector_int_destroy(&edges);
}

int main(void) {
    igraph_t g_1, g_line, g_full, g_loop, g_de_bruijn;
    igraph_vector_int_t vertices, edges;
    igraph_vector_t g_line_weights, g_de_bruijn_weights, error_weights;
    igraph_int_t ec, i;

    igraph_vector_int_init(&vertices, 0);
    igraph_vector_int_init(&edges, 0);
    igraph_vector_init(&g_line_weights, 0);
    igraph_vector_init(&g_de_bruijn_weights, 0);
    igraph_vector_init(&error_weights, 0);
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_small(&g_1, 1, IGRAPH_UNDIRECTED, -1);
    igraph_small(&g_line, 5, IGRAPH_DIRECTED, 0,1, 1,2, 2,3, 3,4, -1);
    igraph_full(&g_full, 20, IGRAPH_UNDIRECTED, 0);
    igraph_small(&g_loop, 1, IGRAPH_DIRECTED, 0,0, -1);

    /* This directed graph has loop edges.
       It also has multi-edges when considered as undirected. */
    igraph_de_bruijn(&g_de_bruijn, 3, 2);

    /* Initialize de bruijn graph weights */
    igraph_rng_seed(igraph_rng_default(), 42);
    ec = igraph_ecount(&g_de_bruijn);
    igraph_vector_resize(&g_de_bruijn_weights, ec);
    for (i = 0; i < ec; ++i) {
        VECTOR(g_de_bruijn_weights)[i] = igraph_rng_get_unif01(igraph_rng_default());
    }

    /* Initialize g_line graph weights */
    igraph_rng_seed(igraph_rng_default(), 42);
    ec = igraph_ecount(&g_line);
    igraph_vector_resize(&g_line_weights, ec);
    for (i = 0; i < ec; ++i) {
        VECTOR(g_line_weights)[i] = igraph_rng_get_unif01(igraph_rng_default());
    }


    /* 1. only vertices required, edges = NULL
          unweighted -> igraph_i_random_walk_adjlist
          weighted   -> igraph_i_random_walk_inclist */
    printf("Only vertices required, edges = NULL:\n");

    printf("Singleton graph with zero steps:\n");
    random_vertex_walk(&g_1, NULL, 0, IGRAPH_OUT, 0, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Singleton graph with one step:\n");
    random_vertex_walk(&g_line, NULL, 0, IGRAPH_OUT, 1, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Line graph:\n");
    random_vertex_walk(&g_line, NULL, 0, IGRAPH_OUT, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Line graph, 10 steps, returns on stuck:\n");
    random_vertex_walk(&g_line, NULL, 0, IGRAPH_OUT, 10, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Line graph backward:\n");
    random_vertex_walk(&g_line, NULL, 4, IGRAPH_IN, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Loop at vertex 0:\n");
    random_vertex_walk(&g_loop, NULL, 0, IGRAPH_OUT, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    igraph_rng_seed(igraph_rng_default(), 42);
    printf("Checking an actual random walk with seed 42:\n");
    random_vertex_walk(&g_full, NULL, 4, IGRAPH_IN, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    /* weighted, directed */
    igraph_random_walk(&g_de_bruijn, &g_de_bruijn_weights, &vertices, NULL, 0, IGRAPH_OUT, 1000, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_int_size(&vertices) == 1001);

    /* weighted, undirecetd */
    igraph_random_walk(&g_de_bruijn, &g_de_bruijn_weights, &vertices, NULL, 0, IGRAPH_ALL, 1000, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_int_size(&vertices) == 1001);

    /* weighted, directed line graph, 4 edges, 10 steps, returns on stuck */
    igraph_random_walk(&g_line, &g_line_weights, &vertices, NULL, 4, IGRAPH_IN, 10, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_int_size(&vertices) == 5);


    /* 2. both vertices and edges required -> igraph_i_random_walk_inclist */
    printf("\nBoth vertices and edges required:\n");

    printf("Singleton graph with zero steps:\n");
    random_walk(&g_1, NULL, 0, IGRAPH_OUT, 0, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Singleton graph with one step:\n");
    random_walk(&g_line, NULL, 0, IGRAPH_OUT, 1, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Line graph:\n");
    random_walk(&g_line, NULL, 0, IGRAPH_OUT, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Line graph, 10 steps, returns on stuck:\n");
    random_walk(&g_line, NULL, 0, IGRAPH_OUT, 10, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Line graph backward:\n");
    random_walk(&g_line, NULL, 4, IGRAPH_IN, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Loop at vertex 0:\n");
    random_walk(&g_loop, NULL, 0, IGRAPH_OUT, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    igraph_rng_seed(igraph_rng_default(), 42);
    printf("Checking an actual random walk with seed 42:\n");
    random_walk(&g_full, NULL, 4, IGRAPH_IN, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    /* weighted, directed */
    igraph_random_walk(&g_de_bruijn, &g_de_bruijn_weights, &vertices, &edges, 0, IGRAPH_OUT, 1000, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_int_size(&vertices) == 1001);
    IGRAPH_ASSERT(igraph_vector_int_size(&edges) == 1000);

    /* weighted, undirecetd */
    igraph_random_walk(&g_de_bruijn, &g_de_bruijn_weights, &vertices, &edges, 0, IGRAPH_ALL, 1000, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_int_size(&vertices) == 1001);
    IGRAPH_ASSERT(igraph_vector_int_size(&edges) == 1000);

    /* weighted, directed line graph, 4 edges, 10 steps, returns on stuck */
    igraph_random_walk(&g_line, &g_line_weights, &vertices, &edges, 4, IGRAPH_IN, 10, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_int_size(&vertices) == 5);
    IGRAPH_ASSERT(igraph_vector_int_size(&edges) == 4);


    /* 3. only edges required, vertices = NULL -> igraph_i_random_walk_inclist */
    printf("\nOnly edges required, vertices = NULL:\n");

    /* weighted, directed */
    igraph_random_walk(&g_de_bruijn, &g_de_bruijn_weights, NULL, &edges, 0, IGRAPH_OUT, 1000, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_int_size(&edges) == 1000);

    /* weighted, undirecetd */
    igraph_random_walk(&g_de_bruijn, &g_de_bruijn_weights, NULL, &edges, 0, IGRAPH_ALL, 1000, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_int_size(&edges) == 1000);


    igraph_set_error_handler(igraph_error_handler_ignore);
    printf("Checking error handling:\n");

    /* 1. only vertices required, edges = NULL
          unweighted -> igraph_i_random_walk_adjlist
          weighted   -> igraph_i_random_walk_inclist */
    printf("\nOnly vertices required, edges = NULL:\n");

    printf("unweighted directed line graph, 10 steps, errors on stuck.\n");
    IGRAPH_ASSERT(igraph_random_walk(&g_line, NULL, &vertices, NULL, 0, IGRAPH_OUT, 10, IGRAPH_RANDOM_WALK_STUCK_ERROR) == IGRAPH_ERWSTUCK);

    printf("Vertex out of range.\n");
    IGRAPH_ASSERT(igraph_random_walk(&g_1, NULL, &vertices, NULL, 10, IGRAPH_OUT, 0, IGRAPH_RANDOM_WALK_STUCK_RETURN) == IGRAPH_EINVAL);

    printf("Negative number of steps.\n");
    IGRAPH_ASSERT(igraph_random_walk(&g_1, NULL, &vertices, NULL, 0, IGRAPH_OUT, -10, IGRAPH_RANDOM_WALK_STUCK_RETURN) == IGRAPH_EINVAL);

    /* weighted, directed line graph, 4 edges, 10 steps, errors on stuck */
    IGRAPH_ASSERT(igraph_random_walk(&g_line, &g_line_weights, &vertices, NULL, 4, IGRAPH_IN, 10, IGRAPH_RANDOM_WALK_STUCK_ERROR) == IGRAPH_ERWSTUCK);

    /* weighted, directed line graph, negative weight value for edge-0 */
    ec = igraph_ecount(&g_line);
    igraph_vector_resize(&error_weights, ec);
    VECTOR(error_weights)[0] = -10;
    IGRAPH_ASSERT(igraph_random_walk(&g_line, &error_weights, &vertices, NULL, 0, IGRAPH_OUT, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN) == IGRAPH_EINVAL);

    /* weighted, directed line graph, invalid weight vector length (!= ec) */
    igraph_vector_resize(&error_weights, 2 * ec);
    IGRAPH_ASSERT(igraph_random_walk(&g_line, &error_weights, &vertices, NULL, 0, IGRAPH_OUT, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN) == IGRAPH_EINVAL);


    /* 2. both vertices and edges required -> igraph_i_random_walk_inclist */
    printf("\nBoth vertices and edges required:\n");

    printf("unweighted directed line graph, 10 steps, errors on stuck.\n");
    IGRAPH_ASSERT(igraph_random_walk(&g_line, NULL, &vertices, &edges, 0, IGRAPH_OUT, 10, IGRAPH_RANDOM_WALK_STUCK_ERROR) == IGRAPH_ERWSTUCK);

    /* weighted, directed line graph, 4 edges, 10 steps, errors on stuck */
    IGRAPH_ASSERT(igraph_random_walk(&g_line, &g_line_weights, &vertices, &edges, 4, IGRAPH_IN, 10, IGRAPH_RANDOM_WALK_STUCK_ERROR) == IGRAPH_ERWSTUCK);


    igraph_destroy(&g_1);
    igraph_destroy(&g_line);
    igraph_destroy(&g_loop);
    igraph_destroy(&g_full);
    igraph_destroy(&g_de_bruijn);

    igraph_vector_destroy(&g_de_bruijn_weights);
    igraph_vector_destroy(&g_line_weights);
    igraph_vector_destroy(&error_weights);
    igraph_vector_int_destroy(&vertices);
    igraph_vector_int_destroy(&edges);

    VERIFY_FINALLY_STACK();

    return 0;
}
