/*
   igraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

int main(void) {

    igraph_t g;
    igraph_vector_int_list_t partitions;
    igraph_vector_int_list_t cuts;
    igraph_real_t value;

    /* Initialize the library. */
    igraph_setup();

    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4,
                 -1);

    igraph_vector_int_list_init(&partitions, 0);
    igraph_vector_int_list_init(&cuts, 0);
    igraph_all_st_mincuts(&g, &value, &cuts, &partitions,
                          /*source=*/ 0, /*target=*/ 4,
                          /*capacity=*/ 0);

    igraph_int_t i, e, m, n = igraph_vector_int_list_size(&partitions);
    printf("Found %" IGRAPH_PRId " cuts, value: %g\n", n, value);
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *vec = igraph_vector_int_list_get_ptr(&partitions, i);
        igraph_vector_int_t *vec2 = igraph_vector_int_list_get_ptr(&cuts, i);
        printf("Partition %" IGRAPH_PRId ": ", i);
        igraph_vector_int_print(vec);
        if (vec2) {
            printf("Cut %" IGRAPH_PRId ":\n", i);
            m = igraph_vector_int_size(vec2);
            for (e = 0; e < m; e++) {
                igraph_int_t from = IGRAPH_FROM(&g, VECTOR(*vec2)[e]), to = IGRAPH_TO(&g, VECTOR(*vec2)[e]);
                printf("  %" IGRAPH_PRId " -> %" IGRAPH_PRId "\n", from, to);
            }
        }
    }

    igraph_vector_int_list_destroy(&partitions);
    igraph_vector_int_list_destroy(&cuts);
    printf("\n");
    igraph_destroy(&g);

    return 0;
}
