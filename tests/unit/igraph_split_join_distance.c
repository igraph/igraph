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
    igraph_vector_t comm1, comm2;
    igraph_integer_t distance12, distance21;

    printf("No vertices:\n");
    igraph_vector_init_int(&comm1, 0);
    igraph_vector_init_int(&comm2, 0);
    IGRAPH_ASSERT(igraph_split_join_distance(&comm1, &comm2, &distance12, &distance21) == IGRAPH_SUCCESS);
    printf("%" IGRAPH_PRId ", %" IGRAPH_PRId "\n", distance12, distance21);
    igraph_vector_destroy(&comm1);
    igraph_vector_destroy(&comm2);

    printf("Comparing 5 separate vertices and one 5-element cluster:\n");
    igraph_vector_init_int(&comm1, 5, 0, 1, 2, 3, 4);
    igraph_vector_init_int(&comm2, 5, 0, 0, 0, 0, 0);
    IGRAPH_ASSERT(igraph_split_join_distance(&comm1, &comm2, &distance12, &distance21) == IGRAPH_SUCCESS);
    printf("%" IGRAPH_PRId ", %" IGRAPH_PRId "\n", distance12, distance21);
    igraph_vector_destroy(&comm1);
    igraph_vector_destroy(&comm2);

    printf("Comparing ((6, 1), (2,4), (3,5,0)) with ((2), (6,0,3), (4,5), (1)):\n");
    igraph_vector_init_int(&comm1, 7, 2, 0, 1, 2, 1, 2, 0);
    igraph_vector_init_int(&comm2, 7, 1, 3, 0, 1, 2, 2, 1);
    IGRAPH_ASSERT(igraph_split_join_distance(&comm1, &comm2, &distance12, &distance21) == IGRAPH_SUCCESS);
    printf("%" IGRAPH_PRId ", %" IGRAPH_PRId "\n", distance12, distance21);
    igraph_vector_destroy(&comm1);
    igraph_vector_destroy(&comm2);

    printf("Comparing ((0,1), (), 2) with ((0), (), (1,2))\n");
    igraph_vector_init_int(&comm1, 3, 0, 0, 2);
    igraph_vector_init_int(&comm2, 3, 0, 2, 2);
    IGRAPH_ASSERT(igraph_split_join_distance(&comm1, &comm2, &distance12, &distance21) == IGRAPH_SUCCESS);
    printf("%" IGRAPH_PRId ", %" IGRAPH_PRId "\n", distance12, distance21);
    igraph_vector_destroy(&comm1);
    igraph_vector_destroy(&comm2);

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("\nExpected to fail nicely:\n\n");

    printf("Differently sized clusterings\n");
    igraph_vector_init_int(&comm1, 3, 0, 1, 2);
    igraph_vector_init_int(&comm2, 5, 0, 0, 0, 0, 0);
    IGRAPH_ASSERT(igraph_split_join_distance(&comm1, &comm2, &distance12, &distance21) == IGRAPH_EINVAL);
    igraph_vector_destroy(&comm1);
    igraph_vector_destroy(&comm2);

    printf("Member index too high\n");
    igraph_vector_init_int(&comm1, 3, 0, 1, 2);
    igraph_vector_init_int(&comm2, 3, 9, 0, 0);
    IGRAPH_ASSERT(igraph_split_join_distance(&comm1, &comm2, &distance12, &distance21) == IGRAPH_EINVAL);
    igraph_vector_destroy(&comm1);
    igraph_vector_destroy(&comm2);


    VERIFY_FINALLY_STACK();
    return 0;
}
