/*
   igraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <igraph.h>
#include "test_utilities.h"

int main(void) {
    igraph_vector_int_t membership;
    igraph_vector_int_list_t groups;
    igraph_vector_int_t group1, group2;

    igraph_vector_int_init(&membership, 0);
    igraph_vector_int_list_init(&groups, 0);

    printf("Test 1: Two groups, no singletons\n");
    {
        igraph_vector_int_init(&group1, 3);
        igraph_vector_int_init(&group2, 2);

        VECTOR(group1)[0] = 0;
        VECTOR(group1)[1] = 1;
        VECTOR(group1)[2] = 2;

        VECTOR(group2)[0] = 3;
        VECTOR(group2)[1] = 4;

        igraph_vector_int_list_clear(&groups);
        igraph_vector_int_list_push_back(&groups, &group1);
        igraph_vector_int_list_push_back(&groups, &group2);

        igraph_groups_to_membership(5, &groups, &membership);
        printf("Groups: "); print_vector_int_list(&groups);
        printf("Membership: "); print_vector_int(&membership);
        printf("\n");

        // Verify: vertices 0,1,2 should be in group 0, vertices 3,4 in group 1
        IGRAPH_ASSERT(VECTOR(membership)[0] == 0);
        IGRAPH_ASSERT(VECTOR(membership)[1] == 0);
        IGRAPH_ASSERT(VECTOR(membership)[2] == 0);
        IGRAPH_ASSERT(VECTOR(membership)[3] == 1);
        IGRAPH_ASSERT(VECTOR(membership)[4] == 1);
    }

    printf("Test 2: Groups with singletons\n");
    {
        igraph_vector_int_init(&group1, 2);
        VECTOR(group1)[0] = 0;
        VECTOR(group1)[1] = 1;

        igraph_vector_int_list_clear(&groups);
        igraph_vector_int_list_push_back(&groups, &group1);

        igraph_groups_to_membership(4, &groups, &membership);
        printf("Groups: "); print_vector_int_list(&groups);
        printf("Membership: "); print_vector_int(&membership);
        printf("\n");

        // Verify: 0,1 in group 0, 2,3 should be singletons (group 1, 2)
        IGRAPH_ASSERT(VECTOR(membership)[0] == 0);
        IGRAPH_ASSERT(VECTOR(membership)[1] == 0);
        // Vertices 2 and 3 should be in different singleton groups
        IGRAPH_ASSERT(VECTOR(membership)[2] >= 1);
        IGRAPH_ASSERT(VECTOR(membership)[3] >= 1);
        IGRAPH_ASSERT(VECTOR(membership)[2] != VECTOR(membership)[3]);
    }

    printf("Test 3: Empty groups (all singletons)\n");
    {
        igraph_vector_int_list_clear(&groups);
        igraph_groups_to_membership(3, &groups, &membership);
        printf("Groups: []\n");
        printf("Membership: "); print_vector_int(&membership);
        printf("\n");

        // All vertices should be in different singleton groups
        IGRAPH_ASSERT(VECTOR(membership)[0] != VECTOR(membership)[1]);
        IGRAPH_ASSERT(VECTOR(membership)[1] != VECTOR(membership)[2]);
        IGRAPH_ASSERT(VECTOR(membership)[0] != VECTOR(membership)[2]);
    }

    printf("Test 4: Zero vertices\n");
    {
        igraph_vector_int_list_clear(&groups);
        igraph_groups_to_membership(0, &groups, &membership);
        printf("Membership: "); print_vector_int(&membership);
        printf("\n");

        IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);
    }

    igraph_vector_int_destroy(&membership);
    igraph_vector_int_list_destroy(&groups);

    VERIFY_FINALLY_STACK();
    return 0;
}
