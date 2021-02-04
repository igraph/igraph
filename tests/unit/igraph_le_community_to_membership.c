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

void print_and_destroy(igraph_vector_t *membership, igraph_vector_t *csize, igraph_matrix_t *mat) {
    printf("Membership: ");
    igraph_vector_print(membership);
    printf("Csize: ");
    igraph_vector_print(csize);
    printf("\n");
    igraph_vector_destroy(membership);
    igraph_matrix_destroy(mat);
}

int main() {
    igraph_matrix_t merges;
    igraph_vector_t membership;
    igraph_vector_t csize;

    igraph_vector_init_int(&csize, 0);

    printf("One member:\n");
    igraph_vector_init_int(&membership, 1, 0);
    igraph_matrix_init(&merges, 0, 2);
    IGRAPH_ASSERT(igraph_le_community_to_membership(&merges, /*steps*/ 0, &membership, &csize) == IGRAPH_SUCCESS);
    print_and_destroy(&membership, &csize, &merges);

    printf("Five singleton clusters, one merge:\n");
    igraph_vector_init_int(&membership, 5, 0, 1, 2, 3, 4);
    {
        int elem[] = {1, 3};
        matrix_init_int_row_major(&merges, 1, 2, elem);
    }
    IGRAPH_ASSERT(igraph_le_community_to_membership(&merges, /*steps*/ 1, &membership, &csize) == IGRAPH_SUCCESS);
    print_and_destroy(&membership, &csize, &merges);

    printf("Six clusters, two merges:\n");
    igraph_vector_init_int(&membership, 12, 0, 0, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5);
    {
        int elem[] = {0,3, 2,5};
        matrix_init_int_row_major(&merges, 2, 2, elem);
    }
    IGRAPH_ASSERT(igraph_le_community_to_membership(&merges, /*steps*/ 2, &membership, &csize) == IGRAPH_SUCCESS);
    print_and_destroy(&membership, &csize, &merges);


    printf("These should fail nicely:\n");
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Five singleton clusters, two merges, but second merge refers to singleton which is already merged.\n");
    igraph_vector_init_int(&membership, 5, 0, 1, 2, 3, 4);
    {
        int elem[] = {1,3, 1,4,};
        matrix_init_int_row_major(&merges, 2, 2, elem);
    }
    IGRAPH_ASSERT(igraph_le_community_to_membership(&merges, /*steps*/ 2, &membership, &csize) == IGRAPH_EINVAL);
    igraph_vector_destroy(&membership);
    igraph_matrix_destroy(&merges);

    printf("Negative cluster index.\n");
    igraph_vector_init_int(&membership, 5, -1, 0, 1, 2, 3);
    {
        int elem[] = {1,2, 3,4,};
        matrix_init_int_row_major(&merges, 2, 2, elem);
    }
    IGRAPH_ASSERT(igraph_le_community_to_membership(&merges, /*steps*/ 2, &membership, &csize) == IGRAPH_EINVAL);
    igraph_vector_destroy(&membership);
    igraph_matrix_destroy(&merges);

    printf("Skip a cluster index.\n");
    igraph_vector_init_int(&membership, 5, 0, 0, 2, 3, 4);
    {
        int elem[] = {1,2, 3,4,};
        matrix_init_int_row_major(&merges, 2, 2, elem);
    }
    IGRAPH_ASSERT(igraph_le_community_to_membership(&merges, /*steps*/ 2, &membership, &csize) == IGRAPH_EINVAL);
    igraph_vector_destroy(&membership);
    igraph_matrix_destroy(&merges);

    printf("Too many steps.\n");
    igraph_vector_init_int(&membership, 5, 0, 1, 2, 3, 4);
    {
        int elem[] = {1,2, 3,4,};
        matrix_init_int_row_major(&merges, 2, 2, elem);
    }
    IGRAPH_ASSERT(igraph_le_community_to_membership(&merges, /*steps*/ 20, &membership, &csize) == IGRAPH_EINVAL);
    igraph_vector_destroy(&membership);
    igraph_matrix_destroy(&merges);

    igraph_vector_destroy(&csize);
    VERIFY_FINALLY_STACK();
    return 0;
}
