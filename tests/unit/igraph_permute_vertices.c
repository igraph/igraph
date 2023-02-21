/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

#include "test_utilities.h"

int main(void) {
    igraph_t graph, pgraph;
    igraph_vector_int_t perm;

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_small(&graph, 6, IGRAPH_DIRECTED,
                 0,1, 1,0, 2,1, 2,2, 5,4, 5,4,
                 -1);

    printf("Graph before permuting vertices:\n");
    print_graph(&graph);

    igraph_vector_int_init_int(&perm, 6,
                               0, 1, 5, 3, 4, 2);

    igraph_permute_vertices(&graph, &pgraph, &perm);

    printf("\nGraph after permuting vertices:\n");
    print_graph(&pgraph);

    igraph_destroy(&pgraph);

    /* Now we test invalid input. */

    /* Out-of-range value in permutation */
    VECTOR(perm)[0] = 7;
    CHECK_ERROR(igraph_permute_vertices(&graph, &pgraph, &perm), IGRAPH_EINVAL);

    /* Duplicate value in permutation */
    VECTOR(perm)[0] = 1;
    CHECK_ERROR(igraph_permute_vertices(&graph, &pgraph, &perm), IGRAPH_EINVAL);

    /* Invalid permutation vector length */
    VECTOR(perm)[0] = 0;
    VECTOR(perm)[2] = 2;
    igraph_vector_int_resize(&perm, 5);
    CHECK_ERROR(igraph_permute_vertices(&graph, &pgraph, &perm), IGRAPH_EINVAL);

    igraph_vector_int_destroy(&perm);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
