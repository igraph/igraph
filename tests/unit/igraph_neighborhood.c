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

void print_and_destroy(igraph_vector_ptr_t *result) {
    int i;
    igraph_vector_t *v;
    for (i = 0; i < igraph_vector_ptr_size(result); i++) {
        v = VECTOR(*result)[i];
        print_vector(v);
        igraph_vector_destroy(v);
        igraph_free(v);
    }
}

int main() {
    igraph_t g_empty, g_lm;
    igraph_vector_ptr_t result;
    igraph_vs_t vids;

    igraph_vector_ptr_init(&result, 0);
    igraph_vs_all(&vids);

    igraph_small(&g_empty, 0, 0, -1);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);

    printf("No vertices:\n");
    IGRAPH_ASSERT(igraph_neighborhood(&g_empty, &result, vids, /*order*/ 1,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ 0) == IGRAPH_SUCCESS);
    print_and_destroy(&result);

    printf("Directed graph with loops and multi-edges, order 0:\n");
    IGRAPH_ASSERT(igraph_neighborhood(&g_lm, &result, vids, /*order*/ 0,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ 0) == IGRAPH_SUCCESS);
    print_and_destroy(&result);

    printf("Directed graph with loops and multi-edges, order 1, ignoring direction:\n");
    IGRAPH_ASSERT(igraph_neighborhood(&g_lm, &result, vids, /*order*/ 1,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ 0) == IGRAPH_SUCCESS);
    print_and_destroy(&result);

    printf("Directed graph with loops and multi-edges, order 1, only checking IGRAPH_IN:\n");
    IGRAPH_ASSERT(igraph_neighborhood(&g_lm, &result, vids, /*order*/ 1,
                  /*mode*/ IGRAPH_IN, /*mindist*/ 0) == IGRAPH_SUCCESS);
    print_and_destroy(&result);

    printf("Directed graph with loops and multi-edges, order 10, ignoring direction:\n");
    IGRAPH_ASSERT(igraph_neighborhood(&g_lm, &result, vids, /*order*/ 10,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ 0) == IGRAPH_SUCCESS);
    print_and_destroy(&result);

    printf("Directed graph with loops and multi-edges, order 2, mindist 2, IGRAPH_OUT:\n");
    IGRAPH_ASSERT(igraph_neighborhood(&g_lm, &result, vids, /*order*/ 2,
                  /*mode*/ IGRAPH_OUT, /*mindist*/ 2) == IGRAPH_SUCCESS);
    print_and_destroy(&result);

    printf("Directed graph with loops and multi-edges, order 4, mindist 4, IGRAPH_ALL:\n");
    IGRAPH_ASSERT(igraph_neighborhood(&g_lm, &result, vids, /*order*/ 4,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ 4) == IGRAPH_SUCCESS);
    print_and_destroy(&result);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Negative order.\n");
    IGRAPH_ASSERT(igraph_neighborhood(&g_lm, &result, vids, /*order*/ -4,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ 4) == IGRAPH_EINVAL);

    printf("Negative mindist.\n");
    IGRAPH_ASSERT(igraph_neighborhood(&g_lm, &result, vids, /*order*/ 4,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ -4) == IGRAPH_EINVAL);

    igraph_vector_ptr_destroy(&result);
    igraph_destroy(&g_empty);
    igraph_destroy(&g_lm);

    VERIFY_FINALLY_STACK();
    return 0;
}
