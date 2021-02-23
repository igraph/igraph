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

void print_and_destroy(igraph_t *g, igraph_integer_t center, igraph_vector_t *order, igraph_error_type_t error) {
    igraph_matrix_t result;
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_star(g, &result, center, order) == error);
    if (error == IGRAPH_SUCCESS) {
        matrix_chop(&result, 1e-13);
        print_matrix(&result);
    }

    igraph_matrix_destroy(&result);
    igraph_destroy(g);
    if(order) {
        igraph_vector_destroy(order);
    }

}

int main() {
    igraph_t g;
    igraph_vector_t order;

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Star of 8 points and a center:\n");
    igraph_small(&g, 9, 0, -1);
    print_and_destroy(&g, 0, NULL, IGRAPH_SUCCESS);

    printf("Star of 8 points and a center in reverse:\n");
    igraph_small(&g, 9, 0, -1);
    igraph_vector_init_int(&order, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    print_and_destroy(&g, 0, &order, IGRAPH_SUCCESS);

    printf("Checking if negative center fails nicely.\n");
    igraph_small(&g, 9, 0, -1);
    print_and_destroy(&g, -10, NULL, IGRAPH_EINVAL);

    printf("Checking if order out of range fails nicely.\n");
    igraph_small(&g, 9, 0, -1);
    igraph_vector_init_int(&order, 9, -1, -1, -1, 10, 10, 10, 2, 1, 0);
    print_and_destroy(&g, 0, &order, IGRAPH_EINVAL);

    VERIFY_FINALLY_STACK();
    return 0;
}
