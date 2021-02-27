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

int main() {
    igraph_t g_0, g_1, g_simple, g_ml;
    igraph_real_t result;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_simple, 6, 1, 0,1, 0,2, 1,2, 1,3, 2,3, 3,4, -1);
    igraph_small(&g_ml, 6, 0, 0,1, 0,2, 1,1, 1,2, 1,3, 2,1, 2,3, 3,4, 3,4, -1);

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("No vertices, transitivity zero:\n");
    IGRAPH_ASSERT(igraph_transitivity_avglocal_undirected(&g_0, &result, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_SUCCESS);
    print_real(stdout, result, "%g");
    printf("\n");

    printf("No vertices, transitivity NaN:\n");
    IGRAPH_ASSERT(igraph_transitivity_avglocal_undirected(&g_0, &result, IGRAPH_TRANSITIVITY_NAN) == IGRAPH_SUCCESS);
    print_real(stdout, result, "%g");
    printf("\n");

    printf("One vertex, transitivity zero:\n");
    IGRAPH_ASSERT(igraph_transitivity_avglocal_undirected(&g_1, &result, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_SUCCESS);
    print_real(stdout, result, "%g");
    printf("\n");

    printf("One vertex, transitivity NaN:\n");
    IGRAPH_ASSERT(igraph_transitivity_avglocal_undirected(&g_1, &result, IGRAPH_TRANSITIVITY_NAN) == IGRAPH_SUCCESS);
    print_real(stdout, result, "%g");
    printf("\n");

    printf("Simple graph:\n");
    IGRAPH_ASSERT(igraph_transitivity_avglocal_undirected(&g_simple, &result, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_SUCCESS);
    print_real(stdout, result, "%g");
    printf("\n");

    printf("Multigraph:\n");
    IGRAPH_ASSERT(igraph_transitivity_avglocal_undirected(&g_ml, &result, IGRAPH_TRANSITIVITY_ZERO) == IGRAPH_SUCCESS);
    print_real(stdout, result, "%g");
    printf("\n");

    printf("Multigraph, TRANSITIVITY_NAN:\n");
    IGRAPH_ASSERT(igraph_transitivity_avglocal_undirected(&g_ml, &result, IGRAPH_TRANSITIVITY_NAN) == IGRAPH_SUCCESS);
    print_real(stdout, result, "%g");
    printf("\n");

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_simple);
    igraph_destroy(&g_ml);

    VERIFY_FINALLY_STACK();
    return 0;
}
