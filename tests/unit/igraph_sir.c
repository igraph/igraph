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

void print_sir(igraph_sir_t *sir) {
    printf("times: ");
    print_vector(&sir->times);
    printf("susceptibles: ");
    print_vector_int(&sir->no_s);
    printf("infected: ");
    print_vector_int(&sir->no_i);
    printf("recovered: ");
    print_vector_int(&sir->no_r);
}


void print_result(igraph_t *g, igraph_real_t beta, igraph_real_t gamma, igraph_integer_t no_sim) {
    igraph_vector_ptr_t result;
    igraph_vector_ptr_init(&result, 0);
    IGRAPH_ASSERT(igraph_sir(g, beta, gamma, no_sim, &result) == IGRAPH_SUCCESS);
    for (int i = 0; i < igraph_vector_ptr_size(&result); i++) {
        print_sir(VECTOR(result)[i]);
        igraph_sir_destroy(VECTOR(result)[i]);
    }
    igraph_vector_ptr_destroy_all(&result);
    printf("\n");
}

int main() {
    igraph_t g_empty, g_lm, g_line, g_1, g_2, g_full;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_small(&g_empty, 0, 0, -1);
    igraph_small(&g_lm, 6, 0, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_2, 2, 0, -1);
    igraph_small(&g_line, 5, 0, 0,1, 1,2, 2,3, 3,4, -1);
    igraph_full(&g_full, 5, 0, IGRAPH_NO_LOOPS);

    printf("Only one person, low recovery rate:\n");
    print_result(&g_1, 0.1, 0.0001, 2);

    printf("Two people, not connected, only one infection expected:\n");
    print_result(&g_2, 1, 1, 2);

    printf("Line:\n");
    print_result(&g_line, 1, 1, 2);

    printf("Line, low infection rate, few infections expected:\n");
    print_result(&g_line, 0.0001, 1, 2);

    printf("Full graph, more infections expected than line with same rates:\n");
    print_result(&g_full, 1, 1, 2);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    IGRAPH_ASSERT(igraph_sir(&g_lm, 1, 1, 1, NULL) == IGRAPH_EINVAL);
    IGRAPH_ASSERT(igraph_sir(&g_empty, 1, 1, 1, NULL) == IGRAPH_EINVAL);

    igraph_destroy(&g_empty);
    igraph_destroy(&g_lm);
    igraph_destroy(&g_line);
    igraph_destroy(&g_1);
    igraph_destroy(&g_2);
    igraph_destroy(&g_full);

    VERIFY_FINALLY_STACK();
    return 0;
}
