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

void random_walk(igraph_t *graph, igraph_vector_t* weights,
    igraph_integer_t start, igraph_neimode_t mode,
    igraph_integer_t steps, igraph_random_walk_stuck_t stuck) {
    igraph_vector_int_t walk;

    igraph_vector_int_init(&walk, 0);
    igraph_random_walk(graph, weights, &walk, start, mode, steps, stuck);
    print_vector_int(&walk);
    igraph_vector_int_destroy(&walk);
}

int main() {
    igraph_t g_1, g_line, g_full, g_loop, graph;
    igraph_vector_int_t walk;
    igraph_vector_t weights;
    igraph_integer_t ec, i;

    igraph_vector_int_init(&walk, 0);
    igraph_vector_init(&weights, 0);

    igraph_small(&g_1, 1, IGRAPH_UNDIRECTED, -1);
    igraph_small(&g_line, 5, IGRAPH_DIRECTED, 0,1, 1,2, 2,3, 3,4, -1);
    igraph_full(&g_full, 20, IGRAPH_UNDIRECTED, 0);
    igraph_small(&g_loop, 1, IGRAPH_DIRECTED, 0,0, -1);

    printf("Singleton graph with zero steps:\n");
    random_walk(&g_1, NULL, 0, IGRAPH_OUT, 0, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Singleton graph with one step:\n");
    random_walk(&g_1, NULL, 0, IGRAPH_OUT, 1, IGRAPH_RANDOM_WALK_STUCK_RETURN);

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

    /* This directed graph has loop edges.
       It also has multi-edges when considered as undirected. */
    igraph_de_bruijn(&graph, 3, 2);
    ec = igraph_ecount(&graph);

    igraph_vector_resize(&weights, ec);
    for (i = 0; i < ec; ++i) {
        VECTOR(weights)[i] = igraph_rng_get_unif01(igraph_rng_default());
    }

    /* weighted, directed */
    igraph_random_walk(&graph, &weights, &walk, 0, IGRAPH_OUT, 1000, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_int_size(&walk) == 1000);

    /* weighted, undirecetd */
    igraph_random_walk(&graph, &weights, &walk, 0, IGRAPH_ALL, 1000, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_int_size(&walk) == 1000);

    ec = igraph_ecount(&g_line);
    igraph_vector_resize(&weights, ec);
    for (i = 0; i < ec; ++i) {
        VECTOR(weights)[i] = igraph_rng_get_unif01(igraph_rng_default());
    }

    /* weighted, directed line graph, 4 edges, 10 steps, returns on stuck */
    igraph_random_walk(&g_line, &weights, &walk, 4, IGRAPH_IN, 10, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_int_size(&walk) == ec);
    
    igraph_set_error_handler(igraph_error_handler_ignore);
    printf("Checking error handling:\n");
    printf("Line graph, 10 steps, errors on stuck.\n");
    IGRAPH_ASSERT(igraph_random_walk(&g_line, NULL, &walk, 0, IGRAPH_OUT, 10, IGRAPH_RANDOM_WALK_STUCK_ERROR) == IGRAPH_ERWSTUCK);

    printf("Vertex out of range.\n");
    IGRAPH_ASSERT(igraph_random_walk(&g_1, NULL, &walk, 10, IGRAPH_OUT, 0, IGRAPH_RANDOM_WALK_STUCK_RETURN) == IGRAPH_EINVAL);

    printf("Negative number of steps.\n");
    IGRAPH_ASSERT(igraph_random_walk(&g_1, NULL, &walk, 0, IGRAPH_OUT, -10, IGRAPH_RANDOM_WALK_STUCK_RETURN) == IGRAPH_EINVAL);

    /* weighted, directed line graph, 4 edges, 10 steps, errors on stuck */
    IGRAPH_ASSERT(igraph_random_walk(&g_line, &weights, &walk, 4, IGRAPH_IN, 10, IGRAPH_RANDOM_WALK_STUCK_ERROR) == IGRAPH_ERWSTUCK);

    /* weighted, directed line graph, negative weight value for edge-0 */
    VECTOR(weights)[0] = -10;
    IGRAPH_ASSERT(igraph_random_walk(&g_line, &weights, &walk, 0, IGRAPH_OUT, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN) == IGRAPH_EINVAL);

    /* weighted, directed line graph, invalid weight vector length (!= ec) */
    igraph_vector_resize(&weights, 2 * ec);
    for (i = 0; i < ec; ++i) {
        VECTOR(weights)[i] = igraph_rng_get_unif01(igraph_rng_default());
    }
    IGRAPH_ASSERT(igraph_random_walk(&g_line, &weights, &walk, 0, IGRAPH_OUT, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN) == IGRAPH_EINVAL);

    igraph_destroy(&g_1);
    igraph_destroy(&g_line);
    igraph_destroy(&g_loop);
    igraph_destroy(&g_full);
    igraph_destroy(&graph);
    
    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&walk);

    VERIFY_FINALLY_STACK();

    return 0;
}
