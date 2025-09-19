/* igraph library.
   Copyright (C) 2010-2024  The igraph development team <igraph@igraph.org>

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
    igraph_vector_int_list_t res, res_all;
    igraph_int_t n;

    igraph_small(&g, 6, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 5,
                 0, 3, 3, 4, 4, 5,
                 3, 2, 3, 5,
                 -1);

    igraph_vector_int_list_init(&res, 0);

    printf("TEST MINLEN\n\n");
    for (igraph_int_t i = 0; i <= 5; i++) {
        igraph_get_all_simple_paths(&g, &res, 0, igraph_vss_1(5), IGRAPH_ALL, i, -1, IGRAPH_UNLIMITED);

        printf("Paths for minlen=%" IGRAPH_PRId ":\n", i);
        print_vector_int_list(&res);
    }

    printf("\nTEST MAXLEN\n\n");
    for (igraph_int_t i = 0; i <= 5; i++) {
        igraph_get_all_simple_paths(&g, &res, 0, igraph_vss_1(5), IGRAPH_ALL, -1, i, IGRAPH_UNLIMITED);

        printf("Paths for maxlen=%" IGRAPH_PRId ":\n", i);
        print_vector_int_list(&res);
    }

    igraph_vector_int_list_init(&res_all, 0);

    igraph_get_all_simple_paths(&g, &res_all, 0, igraph_vss_1(5), IGRAPH_ALL, -1, -1, IGRAPH_UNLIMITED);

    n = igraph_vector_int_list_size(&res);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&res_all) == n);
    for (igraph_int_t i = 0; i < n; i++) {
        IGRAPH_ASSERT(igraph_vector_int_all_e(
                          igraph_vector_int_list_get_ptr(&res, i),
                          igraph_vector_int_list_get_ptr(&res_all, i)));
    }

    igraph_vector_int_list_destroy(&res_all);
    igraph_vector_int_list_destroy(&res);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
