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

    igraph_vector_int_init(&membership, 0);
    igraph_vector_int_list_init(&groups, 0);

    /* Test 1: Two groups, no singletons */
    {
        igraph_vector_int_t group0, group1;

        igraph_vector_int_init(&group0, 3);
        igraph_vector_int_init(&group1, 2);

        VECTOR(group0)[0] = 0;
        VECTOR(group0)[1] = 1;
        VECTOR(group0)[2] = 2;
        VECTOR(group1)[0] = 3;
        VECTOR(group1)[1] = 4;

        igraph_vector_int_list_clear(&groups);
        igraph_vector_int_list_push_back(&groups, &group0);
        igraph_vector_int_list_push_back(&groups, &group1);

        IGRAPH_ASSERT(igraph_groups_to_membership(5, &groups, &membership) == IGRAPH_SUCCESS);

        IGRAPH_ASSERT(VECTOR(membership)[0] == 0);
        IGRAPH_ASSERT(VECTOR(membership)[1] == 0);
        IGRAPH_ASSERT(VECTOR(membership)[2] == 0);
        IGRAPH_ASSERT(VECTOR(membership)[3] == 1);
        IGRAPH_ASSERT(VECTOR(membership)[4] == 1);
    }

    /* Test 2: Groups with singletons */
    {
        igraph_vector_int_t group0;

        igraph_vector_int_init(&group0, 2);

        VECTOR(group0)[0] = 0;
        VECTOR(group0)[1] = 1;

        igraph_vector_int_list_clear(&groups);
        igraph_vector_int_list_push_back(&groups, &group0);

        IGRAPH_ASSERT(igraph_groups_to_membership(4, &groups, &membership) == IGRAPH_SUCCESS);

        IGRAPH_ASSERT(VECTOR(membership)[0] == 0);
        IGRAPH_ASSERT(VECTOR(membership)[1] == 0);
        IGRAPH_ASSERT(VECTOR(membership)[2] == 1);
        IGRAPH_ASSERT(VECTOR(membership)[3] == 2);
    }

    /* Test 3: Empty groups, all singletons */
    {
        igraph_vector_int_list_clear(&groups);

        IGRAPH_ASSERT(igraph_groups_to_membership(3, &groups, &membership) == IGRAPH_SUCCESS);

        IGRAPH_ASSERT(VECTOR(membership)[0] == 0);
        IGRAPH_ASSERT(VECTOR(membership)[1] == 1);
        IGRAPH_ASSERT(VECTOR(membership)[2] == 2);
    }

    /* Test 4: Zero vertices */
    {
        igraph_vector_int_list_clear(&groups);

        IGRAPH_ASSERT(igraph_groups_to_membership(0, &groups, &membership) == IGRAPH_SUCCESS);
        IGRAPH_ASSERT(igraph_vector_int_size(&membership) == 0);
    }

    /* Test 5: Multiple groups with singletons */
    {
        igraph_vector_int_t group0, group1, group2;

        igraph_vector_int_init(&group0, 2);
        igraph_vector_int_init(&group1, 2);
        igraph_vector_int_init(&group2, 1);

        VECTOR(group0)[0] = 0;
        VECTOR(group0)[1] = 1;
        VECTOR(group1)[0] = 3;
        VECTOR(group1)[1] = 4;
        VECTOR(group2)[0] = 6;

        igraph_vector_int_list_clear(&groups);
        igraph_vector_int_list_push_back(&groups, &group0);
        igraph_vector_int_list_push_back(&groups, &group1);
        igraph_vector_int_list_push_back(&groups, &group2);

        IGRAPH_ASSERT(igraph_groups_to_membership(8, &groups, &membership) == IGRAPH_SUCCESS);

        IGRAPH_ASSERT(VECTOR(membership)[0] == 0);
        IGRAPH_ASSERT(VECTOR(membership)[1] == 0);
        IGRAPH_ASSERT(VECTOR(membership)[2] == 3);  /* singleton */
        IGRAPH_ASSERT(VECTOR(membership)[3] == 1);
        IGRAPH_ASSERT(VECTOR(membership)[4] == 1);
        IGRAPH_ASSERT(VECTOR(membership)[5] == 4);  /* singleton */
        IGRAPH_ASSERT(VECTOR(membership)[6] == 2);
        IGRAPH_ASSERT(VECTOR(membership)[7] == 5);  /* singleton */
    }

    VERIFY_FINALLY_STACK();

    /* Test error conditions */
    /* Test 6: Duplicate vertex across groups */
    {
        igraph_vector_int_t group0, group1;

        igraph_vector_int_init(&group0, 2);
        igraph_vector_int_init(&group1, 2);

        VECTOR(group0)[0] = 0;
        VECTOR(group0)[1] = 1;
        VECTOR(group1)[0] = 1;  /* Vertex 1 is in both groups */
        VECTOR(group1)[1] = 2;

        igraph_vector_int_list_clear(&groups);
        igraph_vector_int_list_push_back(&groups, &group0);
        igraph_vector_int_list_push_back(&groups, &group1);

        CHECK_ERROR(igraph_groups_to_membership(3, &groups, &membership), IGRAPH_EINVAL);
    }

    igraph_vector_int_destroy(&membership);
    igraph_vector_int_list_destroy(&groups);

    VERIFY_FINALLY_STACK();
    return 0;
}
