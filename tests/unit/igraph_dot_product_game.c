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

int main(void) {
    igraph_t g;
    igraph_matrix_t vecs;

    igraph_rng_seed(igraph_rng_default(), 42);

    printf("No vertices:\n");
    {
        matrix_init_real_row_major(&vecs, 0, 0, NULL);
        igraph_dot_product_game(&g, &vecs, 1);
        print_graph_canon(&g);
        igraph_matrix_destroy(&vecs);
        igraph_destroy(&g);
    }

    printf("5 vertices, 0-1, 0-4, 1-4 are connected:\n");
    {
        igraph_real_t elems[] = {
            1, 0, 0, -10, 100,
            0, 0.01, 0, -100, 100,
            0, 0, 0, 0, 0,
            1, 1, 0, 0, 0
        };
        matrix_init_real_row_major(&vecs, 4, 5, elems);
        igraph_dot_product_game(&g, &vecs, 1);
        print_graph_canon(&g);
        igraph_matrix_destroy(&vecs);
        igraph_destroy(&g);
    }

    VERIFY_FINALLY_STACK();
    return 0;
}
