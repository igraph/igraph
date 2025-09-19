/*
   igraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdlib.h>

int main(void) {

    igraph_t g;
    igraph_vector_int_list_t result;
    igraph_int_t n;
    igraph_int_t alpha;
    const int params[] = {4, -1, 2, 2, 0, 0, -1, -1};

    /* Initialize the library. */
    igraph_setup();

    igraph_vector_int_list_init(&result, 0);

    igraph_kary_tree(&g, 5, 2, IGRAPH_TREE_OUT);
    for (size_t j = 0; j < sizeof(params) / (2 * sizeof(params[0])); j++) {
        if (params[2 * j + 1] != 0) {
            igraph_independent_vertex_sets(&g, &result, params[2 * j], params[2 * j + 1], IGRAPH_UNLIMITED);
        } else {
            igraph_largest_independent_vertex_sets(&g, &result);
        }
        n = igraph_vector_int_list_size(&result);
        printf("%" IGRAPH_PRId " independent sets found\n", n);
        for (igraph_int_t i = 0; i < n; i++) {
            igraph_vector_int_print(igraph_vector_int_list_get_ptr(&result, i));
        }
    }
    igraph_destroy(&g);

    igraph_kary_tree(&g, 10, 2, IGRAPH_TREE_OUT);
    igraph_maximal_independent_vertex_sets(&g, &result, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED);
    n = igraph_vector_int_list_size(&result);
    printf("%" IGRAPH_PRId " maximal independent sets found\n", n);
    for (igraph_int_t i = 0; i < n; i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&result, i));
    }

    igraph_maximal_independent_vertex_sets(&g, &result, 4, 5, IGRAPH_UNLIMITED);
    n = igraph_vector_int_list_size(&result);
    printf("%" IGRAPH_PRId " maximal independent sets between sizes 4 and 5 found\n", n);
    for (igraph_int_t i = 0; i < n; i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&result, i));
    }

    igraph_independence_number(&g, &alpha);
    printf("alpha=%" IGRAPH_PRId "\n", alpha);

    igraph_destroy(&g);

    igraph_vector_int_list_destroy(&result);

    return 0;
}
