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

char *names[] = {
    [IGRAPH_COMMCMP_VI] = "VI",
    [IGRAPH_COMMCMP_NMI] = "NMI",
    [IGRAPH_COMMCMP_SPLIT_JOIN] = "Split-join",
    [IGRAPH_COMMCMP_RAND] = "Rand",
    [IGRAPH_COMMCMP_ADJUSTED_RAND] = "Adjusted Rand",
};


void compare_and_print(igraph_vector_t *comm1, igraph_vector_t *comm2, igraph_community_comparison_t t, igraph_error_type_t e) {
    igraph_real_t result;
    printf("%s result: ", names[t]);
    IGRAPH_ASSERT(igraph_compare_communities(comm1, comm2, &result, t) == e);
    if (e == IGRAPH_EINVAL) {
        printf("failed as expected\n");
    } else {
        print_real(stdout, result, "%g");
        printf("\n");
    }
}


int main() {
    igraph_vector_t comm1, comm2;
    igraph_real_t result;

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Only one member, both partitions equal to whole set:\n");
    igraph_vector_init_int(&comm1, 1, 0);
    igraph_vector_init_int(&comm2, 1, 0);


    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_VI, IGRAPH_SUCCESS);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_RAND, IGRAPH_EINVAL);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_ADJUSTED_RAND, IGRAPH_EINVAL);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_NMI, IGRAPH_SUCCESS);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_SPLIT_JOIN, IGRAPH_SUCCESS);

    igraph_vector_destroy(&comm1);
    igraph_vector_destroy(&comm2);

    printf("\nEmpty sets:\n");
    igraph_vector_init(&comm1, 0);
    igraph_vector_init(&comm2, 0);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_VI, IGRAPH_SUCCESS);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_RAND, IGRAPH_EINVAL);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_ADJUSTED_RAND, IGRAPH_EINVAL);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_NMI, IGRAPH_SUCCESS);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_SPLIT_JOIN, IGRAPH_SUCCESS);

    igraph_vector_destroy(&comm1);
    igraph_vector_destroy(&comm2);

    printf("\nTwo equal, differenly labeled partitions:\n");
    igraph_vector_init_int(&comm1, 2, 0, 1);
    igraph_vector_init_int(&comm2, 2, 1, 0);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_VI, IGRAPH_SUCCESS);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_RAND, IGRAPH_SUCCESS);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_ADJUSTED_RAND, IGRAPH_SUCCESS);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_NMI, IGRAPH_SUCCESS);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_SPLIT_JOIN, IGRAPH_SUCCESS);

    igraph_vector_destroy(&comm1);
    igraph_vector_destroy(&comm2);

    printf("\nTwo different partitions: ((5,1), (8,3,4), (0,6,7,2,9)) and ((5,8), (1,3,4,0), (6,7,2,9))\n");
    igraph_vector_init_int(&comm1, 10, 2, 0, 2, 1, 1, 0, 2, 2, 1, 2);
    igraph_vector_init_int(&comm2, 10, 1, 1, 2, 1, 1, 0, 2, 2, 0, 2);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_VI, IGRAPH_SUCCESS);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_RAND, IGRAPH_SUCCESS);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_ADJUSTED_RAND, IGRAPH_SUCCESS);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_NMI, IGRAPH_SUCCESS);
    compare_and_print(&comm1, &comm2, IGRAPH_COMMCMP_SPLIT_JOIN, IGRAPH_SUCCESS);

    igraph_vector_destroy(&comm1);
    igraph_vector_destroy(&comm2);

    VERIFY_FINALLY_STACK();
    return 0;
}
