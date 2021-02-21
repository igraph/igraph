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
    igraph_vector_t data, result;

    igraph_set_error_handler(igraph_error_handler_ignore);
    igraph_vector_init_int(&result, 0);

    printf("No values, binwidth 0 should fail.\n");
    igraph_vector_init_int(&data, 0);
    IGRAPH_ASSERT(igraph_running_mean(&data, &result, /*binwidth*/ 0) == IGRAPH_EINVAL);
    igraph_vector_destroy(&data);

    printf("No values, binwidth 1 should fail.\n");
    igraph_vector_init_int(&data, 0);
    IGRAPH_ASSERT(igraph_running_mean(&data, &result, /*binwidth*/ 1) == IGRAPH_EINVAL);
    igraph_vector_destroy(&data);

    printf("One value, binwidth 1:\n");
    igraph_vector_init_int(&data, 1, 1);
    IGRAPH_ASSERT(igraph_running_mean(&data, &result, /*binwidth*/ 1) == IGRAPH_SUCCESS);
    print_vector(&result);
    igraph_vector_destroy(&data);

    printf("1, 2, 3, 4, 5, binwidth 1:\n");
    igraph_vector_init_int(&data, 5, 1, 2, 3, 4, 5);
    IGRAPH_ASSERT(igraph_running_mean(&data, &result, /*binwidth*/ 1) == IGRAPH_SUCCESS);
    print_vector(&result);
    igraph_vector_destroy(&data);

    printf("1, 2, 3, 4, 5, binwidth 2:\n");
    igraph_vector_init_int(&data, 5, 1, 2, 3, 4, 5);
    IGRAPH_ASSERT(igraph_running_mean(&data, &result, /*binwidth*/ 2) == IGRAPH_SUCCESS);
    print_vector(&result);
    igraph_vector_destroy(&data);

    igraph_vector_destroy(&result);

    VERIFY_FINALLY_STACK();
    return 0;
}
