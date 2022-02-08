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

void random_walk(igraph_t *graph, igraph_integer_t start, igraph_neimode_t mode,
    igraph_integer_t steps, igraph_random_walk_stuck_t stuck) {
    igraph_vector_int_t walk;

    igraph_vector_int_init(&walk, 0);
    igraph_random_walk(graph, &walk, start, mode, steps, stuck);
    print_vector_int(&walk);
    igraph_vector_int_destroy(&walk);
}

int main() {
    igraph_t g_1, g_line, g_full, g_loop;
    igraph_vector_int_t walk;

    igraph_vector_int_init(&walk, 0);

    igraph_small(&g_1, 1, IGRAPH_UNDIRECTED, -1);
    igraph_small(&g_line, 5, IGRAPH_DIRECTED, 0,1, 1,2, 2,3, 3,4, -1);
    igraph_full(&g_full, 20, IGRAPH_UNDIRECTED, 0);
    igraph_small(&g_loop, 1, IGRAPH_DIRECTED, 0,0, -1);

    printf("Singleton graph with zero steps:\n");
    random_walk(&g_1, 0, IGRAPH_OUT, 0, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Singleton graph with one step:\n");
    random_walk(&g_1, 0, IGRAPH_OUT, 1, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Line graph:\n");
    random_walk(&g_line, 0, IGRAPH_OUT, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Line graph, 10 steps, returns on stuck:\n");
    random_walk(&g_line, 0, IGRAPH_OUT, 10, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Line graph backward:\n");
    random_walk(&g_line, 4, IGRAPH_IN, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    printf("Loop at vertex 0:\n");
    random_walk(&g_loop, 0, IGRAPH_OUT, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    igraph_rng_seed(igraph_rng_default(), 42);
    printf("Checking an actual random walk with seed 42:\n");
    random_walk(&g_full, 4, IGRAPH_IN, 4, IGRAPH_RANDOM_WALK_STUCK_RETURN);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);
    printf("Checking error handling:\n");
    printf("Line graph, 10 steps, errors on stuck.\n");
    IGRAPH_ASSERT(igraph_random_walk(&g_line, &walk, 0, IGRAPH_OUT, 10, IGRAPH_RANDOM_WALK_STUCK_ERROR) == IGRAPH_ERWSTUCK);

    printf("Vertex out of range.\n");
    IGRAPH_ASSERT(igraph_random_walk(&g_1, &walk, 10, IGRAPH_OUT, 0, IGRAPH_RANDOM_WALK_STUCK_RETURN) == IGRAPH_EINVAL);

    printf("Negative number of steps.\n");
    IGRAPH_ASSERT(igraph_random_walk(&g_1, &walk, 0, IGRAPH_OUT, -10, IGRAPH_RANDOM_WALK_STUCK_RETURN) == IGRAPH_EINVAL);

    igraph_destroy(&g_1);
    igraph_destroy(&g_line);
    igraph_destroy(&g_loop);
    igraph_destroy(&g_full);

    igraph_vector_int_destroy(&walk);

    return 0;
}
