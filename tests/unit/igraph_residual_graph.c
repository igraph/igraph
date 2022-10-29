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

    igraph_t g, residual, expected_residual;
    igraph_vector_t capacity, residual_capacity, flow, expected_residual_capacity;
    igraph_bool_t iso;

    igraph_vector_init(&residual_capacity, 0);

    igraph_small(&g, 6, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 1, 2, 1, 3, 2, 4, 3, 4, 3, 5, 4, 5, -1);
    igraph_small(&expected_residual, 6, IGRAPH_DIRECTED,
         0, 1, 1, 2, 1, 3, 2, 4, 3, 5, 4, 5, -1);
    igraph_vector_init_int_end(&capacity, -1, 4, 2, 2, 3, 4, 1, 2, 5, -1);
    igraph_vector_init_int_end(&flow, -1, 3, 2, 1, 2, 3, 1, 1, 4, -1);
    igraph_vector_init_int_end(&expected_residual_capacity, -1, 1, 1, 1, 1, 1, 1, -1);

    igraph_residual_graph(&g, &capacity, &residual, &residual_capacity, &flow);

    /*    tests    */

    IGRAPH_ASSERT(!igraph_isomorphic(&residual, &expected_residual, &iso));

    IGRAPH_ASSERT(iso);

    IGRAPH_ASSERT(igraph_vector_all_e(&expected_residual_capacity, &residual_capacity));

    /*    cleanup    */

    igraph_vector_destroy(&capacity);
    igraph_vector_destroy(&residual_capacity);
    igraph_vector_destroy(&flow);
    igraph_vector_destroy(&expected_residual_capacity);

    igraph_destroy(&g);
    igraph_destroy(&residual);
    igraph_destroy(&expected_residual);

    VERIFY_FINALLY_STACK();
    return 0;
}
