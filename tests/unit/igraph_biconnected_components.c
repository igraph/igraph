/*
   igraph library.
   Copyright (C) 2006-2024  The igraph development team <igraph@igraph.org>

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
    igraph_vector_int_list_t result;
    igraph_vector_int_t ap;
    igraph_int_t no;

    igraph_vector_int_list_init(&result, 0);
    igraph_vector_int_init(&ap, 0);
    igraph_small(&g, 10 /* extra isolated vertex */, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,3, 3,0,
                 2,4, 4,5, 5,2,
                 5,6,
                 7,8,
                 -1);

    printf("Vertices in biconnected components:\n");
    igraph_biconnected_components(&g, &no, NULL, NULL, /* components */ &result, NULL);
    IGRAPH_ASSERT(no == igraph_vector_int_list_size(&result));
    print_vector_int_list(&result);

    printf("Edges in biconnected components:\n");
    igraph_biconnected_components(&g, &no, NULL, /* component_edges */ &result, NULL, NULL);
    IGRAPH_ASSERT(no == igraph_vector_int_list_size(&result));
    print_vector_int_list(&result);

    printf("Edges in biconnected component spanning trees:\n");
    igraph_biconnected_components(&g, &no, /* tree_edges */ &result, NULL, NULL, &ap);
    IGRAPH_ASSERT(no == igraph_vector_int_list_size(&result));
    print_vector_int_list(&result);

    printf("Articulation points:\n");
    print_vector_int(&ap);

    igraph_destroy(&g);
    igraph_vector_int_destroy(&ap);
    igraph_vector_int_list_destroy(&result);

    VERIFY_FINALLY_STACK();

    return 0;
}
