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
    igraph_t g;
    igraph_matrix_t result;

    igraph_rng_seed(igraph_rng_default(), 42);

    printf("No vertices:\n");
    igraph_small(&g, 0, 0, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_random_3d(&g, &result) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_destroy(&g);
    igraph_matrix_destroy(&result);

    printf("One vertex:\n");
    igraph_small(&g, 1, 0, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_random_3d(&g, &result) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_destroy(&g);
    igraph_matrix_destroy(&result);

    printf("10 vertices:\n");
    igraph_small(&g, 10, 0, 0,1, 0,1, 0,1, 2,2, 2,2, 2,2, 3,4, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_random_3d(&g, &result) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_destroy(&g);
    igraph_matrix_destroy(&result);

    VERIFY_FINALLY_STACK();
    return 0;
}
