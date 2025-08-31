/*
   igraph library.
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

void print_and_clear(igraph_graph_list_t *result) {
    igraph_int_t i;
    igraph_t *g;
    for (i = 0; i < igraph_graph_list_size(result); i++) {
        g = igraph_graph_list_get_ptr(result, i);
        if (igraph_vcount(g) == 0) {
            printf("null graph\n");
        } else if (igraph_ecount(g) == 0 && igraph_vcount(g) == 1) {
            printf("One vertex, no edges\n");
        } else {
            print_graph_canon(g);
        }
    }
    printf("\n");
    igraph_graph_list_clear(result);
}

int main(void) {
    igraph_t g_empty, g_lm;
    igraph_graph_list_t result;
    igraph_vs_t vids;

    igraph_graph_list_init(&result, 0);
    igraph_vs_all(&vids);

    igraph_small(&g_empty, 0, 0, -1);
    igraph_small(&g_lm, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);

    printf("No vertices:\n");
    IGRAPH_ASSERT(igraph_neighborhood_graphs(&g_empty, &result, vids, /*order*/ 1,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ 0) == IGRAPH_SUCCESS);
    print_and_clear(&result);

    printf("Directed graph with loops and multi-edges, order 0:\n");
    IGRAPH_ASSERT(igraph_neighborhood_graphs(&g_lm, &result, vids, /*order*/ 0,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ 0) == IGRAPH_SUCCESS);
    print_and_clear(&result);

    printf("Directed graph with loops and multi-edges, order 1, ignoring direction:\n");
    IGRAPH_ASSERT(igraph_neighborhood_graphs(&g_lm, &result, vids, /*order*/ 1,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ 0) == IGRAPH_SUCCESS);
    print_and_clear(&result);

    printf("Directed graph with loops and multi-edges, order 1, only checking IGRAPH_IN:\n");
    IGRAPH_ASSERT(igraph_neighborhood_graphs(&g_lm, &result, vids, /*order*/ 1,
                  /*mode*/ IGRAPH_IN, /*mindist*/ 0) == IGRAPH_SUCCESS);
    print_and_clear(&result);

    printf("Directed graph with loops and multi-edges, order 10, ignoring direction:\n");
    IGRAPH_ASSERT(igraph_neighborhood_graphs(&g_lm, &result, vids, /*order*/ 10,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ 0) == IGRAPH_SUCCESS);
    print_and_clear(&result);

    printf("Directed graph with loops and multi-edges, order 2, mindist 2, IGRAPH_OUT:\n");
    IGRAPH_ASSERT(igraph_neighborhood_graphs(&g_lm, &result, vids, /*order*/ 2,
                  /*mode*/ IGRAPH_OUT, /*mindist*/ 2) == IGRAPH_SUCCESS);
    print_and_clear(&result);

    printf("Directed graph with loops and multi-edges, order 4, mindist 4, IGRAPH_ALL:\n");
    IGRAPH_ASSERT(igraph_neighborhood_graphs(&g_lm, &result, vids, /*order*/ 4,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ 4) == IGRAPH_SUCCESS);
    print_and_clear(&result);

    printf("Directed graph with loops and multi-edges, infinite order, IGRAPH_OUT:\n");
    igraph_neighborhood_graphs(&g_lm, &result, vids, /*order*/ -1,
                               /*mode*/ IGRAPH_OUT, /*mindist*/ 0);
    print_and_clear(&result);

    printf("Directed graph with loops and multi-edges, infinite order, mindist 2, IGRAPH_OUT:\n");
    igraph_neighborhood_graphs(&g_lm, &result, vids, /*order*/ -1,
                               /*mode*/ IGRAPH_OUT, /*mindist*/ 2);
    print_and_clear(&result);

    printf("Directed graph with loops and multi-edges, infinite order, mindist 2, IGRAPH_IN:\n");
    igraph_neighborhood_graphs(&g_lm, &result, vids, /*order*/ -1,
                               /*mode*/ IGRAPH_IN, /*mindist*/ 2);
    print_and_clear(&result);

    VERIFY_FINALLY_STACK();

    printf("Negative mindist.\n");
    CHECK_ERROR(igraph_neighborhood_graphs(&g_lm, &result, vids, /*order*/ 4,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ -4), IGRAPH_EINVAL);

    igraph_graph_list_destroy(&result);
    igraph_destroy(&g_empty);
    igraph_destroy(&g_lm);

    VERIFY_FINALLY_STACK();
    return 0;
}
