/*
   igraph library.
   Copyright (C) 2023-2024  The igraph development team <igraph@igraph.org>

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

void print_and_destroy(igraph_t *g,
                      igraph_real_t value,
                      igraph_vector_int_list_t *partitions,
                      igraph_vector_int_list_t *cuts) {
    igraph_int_t i, e, m, n = igraph_vector_int_list_size(partitions);
    printf("Found %" IGRAPH_PRId " cuts, value: %g\n", n, value);
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *vec = igraph_vector_int_list_get_ptr(partitions, i);
        igraph_vector_int_t *vec2 = cuts ? igraph_vector_int_list_get_ptr(cuts, i) : 0;
        printf("Partition %" IGRAPH_PRId ": ", i);
        igraph_vector_int_print(vec);
        if (vec2) {
            printf("Cut %" IGRAPH_PRId ":\n", i);
            m = igraph_vector_int_size(vec2);
            for (e = 0; e < m; e++) {
                igraph_int_t from = IGRAPH_FROM(g, VECTOR(*vec2)[e]), to = IGRAPH_TO(g, VECTOR(*vec2)[e]);
                if (igraph_is_directed(g)) {
                    printf("  %" IGRAPH_PRId " -> %" IGRAPH_PRId "\n", from, to);
                } else {
                    printf("  %" IGRAPH_PRId " -- %" IGRAPH_PRId "\n", from, to);
                }
            }
        }
    }

    igraph_vector_int_list_destroy(partitions);
    if (cuts) {
        igraph_vector_int_list_destroy(cuts);
    }
    printf("\n");
}

int main(void) {

    igraph_t g;
    igraph_vector_int_list_t partitions;
    igraph_vector_int_list_t cuts;
    igraph_real_t value;

    printf("Graph 1:\n");
    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4,
                 -1);

    igraph_vector_int_list_init(&partitions, 0);
    igraph_vector_int_list_init(&cuts, 0);
    igraph_all_st_mincuts(&g, &value, &cuts, &partitions,
                          /*source=*/ 0, /*target=*/ 4,
                          /*capacity=*/ NULL);

    print_and_destroy(&g, value, &partitions, &cuts);
    igraph_destroy(&g);

    /* ---------------------------------------------------------------- */

    printf("Graph 2:\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0, 1, 1, 2, 1, 3, 2, 4, 3, 4, 4, 5, -1);
    igraph_vector_int_list_init(&partitions, 0);
    igraph_vector_int_list_init(&cuts, 0);
    igraph_all_st_mincuts(&g, &value, &cuts, &partitions,
                          /*source=*/ 0, /*target=*/ 5, /*capacity=*/ NULL);

    print_and_destroy(&g, value, &partitions, &cuts);
    igraph_destroy(&g);

    /* ---------------------------------------------------------------- */

    printf("Graph 3:\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0, 1, 1, 2, 1, 3, 2, 4, 3, 4, 4, 5, -1);
    igraph_vector_int_list_init(&partitions, 0);
    igraph_vector_int_list_init(&cuts, 0);
    igraph_all_st_mincuts(&g, &value, &cuts, &partitions,
                          /*source=*/ 0, /*target=*/ 4, /*capacity=*/ NULL);

    print_and_destroy(&g, value, &partitions, &cuts);
    igraph_destroy(&g);

    /* ---------------------------------------------------------------- */

    printf("Graph 4:\n");
    igraph_small(&g, 9, IGRAPH_DIRECTED, 0, 1, 0, 2, 1, 3, 2, 3,
                 1, 4, 4, 2, 1, 5, 5, 2, 1, 6, 6, 2, 1, 7, 7, 2, 1, 8, 8, 2,
                 -1);
    igraph_vector_int_list_init(&partitions, 0);
    igraph_vector_int_list_init(&cuts, 0);
    igraph_all_st_mincuts(&g, &value, &cuts, &partitions,
                          /*source=*/ 0, /*target=*/ 3, /*capacity=*/ NULL);

    print_and_destroy(&g, value, &partitions, &cuts);
    igraph_destroy(&g);

    /* ---------------------------------------------------------------- */

    printf("Graph 5:\n");
    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 1, 0, 2, 0, 2, 1, 3, 2,
                 -1);

    igraph_vector_int_list_init(&partitions, 0);
    igraph_vector_int_list_init(&cuts, 0);
    igraph_all_st_mincuts(&g, &value, &cuts, &partitions,
                          /*source=*/ 2, /*target=*/ 0, /*capacity=*/ NULL);

    print_and_destroy(&g, value, &partitions, &cuts);
    igraph_destroy(&g);

    /* ---------------------------------------------------------------- */

    printf("Graph 6:\n");
    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 1, 0, 2, 0, 2, 1, 2, 3,
                 -1);

    igraph_vector_int_list_init(&partitions, 0);
    igraph_vector_int_list_init(&cuts, 0);
    igraph_all_st_mincuts(&g, &value, &cuts, &partitions,
                          /*source=*/ 2, /*target=*/ 0, /*capacity=*/ NULL);

    print_and_destroy(&g, value, &partitions, &cuts);
    igraph_destroy(&g);

    /* ---------------------------------------------------------------- */

    printf("Graph 7:\n");
    igraph_small(&g, 9, IGRAPH_DIRECTED,
                 0, 4,  0, 7,  1, 6,  2, 1,  3, 8,  4, 0,  4, 2,
                 4, 5,  5, 0,  5, 3,  6, 7,  7, 8,
                 -1);

    igraph_vector_int_list_init(&partitions, 0);
    igraph_vector_int_list_init(&cuts, 0);
    igraph_all_st_mincuts(&g, &value, &cuts, &partitions,
                          /*source=*/ 0, /*target=*/ 8, /*capacity=*/ NULL);

    print_and_destroy(&g, value, &partitions, &cuts);
    igraph_destroy(&g);

    /* ---------------------------------------------------------------- */

    printf("Graph 8:\n");
    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 0, 2, 3, 0, 1, 2, 5, 1, 2, 4, 3, 5, 5, 4,
                 -1);

    igraph_vector_int_list_init(&partitions, 0);
    igraph_vector_int_list_init(&cuts, 0);
    igraph_all_st_mincuts(&g, &value, &cuts, &partitions,
                          /*source=*/ 3, /*target=*/ 4, /*capacity=*/ NULL);

    print_and_destroy(&g, value, &partitions, &cuts);
    igraph_destroy(&g);

    /* ---------------------------------------------------------------- */

    printf("Graph 9:\n");
    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 0, 1, 2, 0, 3, 0, 4, 0, 0, 5, 3, 1, 1, 4, 3, 2, 2, 5, 4, 3, 5, 3, 5, 4,
                 -1);

    igraph_vector_int_list_init(&partitions, 0);
    igraph_vector_int_list_init(&cuts, 0);
    igraph_all_st_mincuts(&g, &value, &cuts, &partitions,
                          /*source=*/ 3, /*target=*/ 0, /*capacity=*/ NULL);

    print_and_destroy(&g, value, &partitions, &cuts);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
