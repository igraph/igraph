/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.h"

int main(void) {
    igraph_t g_turan;
    igraph_t g_multipartite;
    igraph_vector_int_t partitions;
    igraph_vector_int_t types_turan;
    igraph_vector_int_t types_multipartite;
    igraph_bool_t res;

    printf("Empty graph with zero vertices:");
    igraph_vector_int_init(&types_turan, 0);
    igraph_turan(&g_turan, &types_turan, 0, 10);

    print_graph_canon(&g_turan);

    printf("\nPartition type:\n");
    igraph_vector_int_print(&types_turan);

    igraph_vector_int_destroy(&types_turan);
    igraph_destroy(&g_turan);


    printf("\nTuran graph with 10 vertices and 1 partition:");
    igraph_vector_int_init(&types_turan, 0);
    igraph_turan(&g_turan, &types_turan, 10, 1);

    print_graph_canon(&g_turan);

    printf("\nPartition type:\n");
    igraph_vector_int_print(&types_turan);

    igraph_vector_int_destroy(&types_turan);
    igraph_destroy(&g_turan);


    printf("\nTuran graph with 4 vertices and 6 partitions\n");
    printf("gives 4 vertices and 1 partition:");
    igraph_vector_int_init(&types_turan, 0);
    igraph_turan(&g_turan, &types_turan, 4, 6);

    print_graph_canon(&g_turan);

    printf("\nPartition type:\n");
    igraph_vector_int_print(&types_turan);

    igraph_vector_int_destroy(&types_turan);
    igraph_destroy(&g_turan);

    igraph_vector_int_init(&partitions, 4);
    igraph_vector_int_init(&types_multipartite, 0);
    igraph_vector_int_init(&types_turan, 0);

    VECTOR(partitions)[0] = 3;
    VECTOR(partitions)[1] = 3;
    VECTOR(partitions)[2] = 3;
    VECTOR(partitions)[3] = 4;

    igraph_full_multipartite(&g_multipartite, &types_multipartite, &partitions,
                        IGRAPH_UNDIRECTED, IGRAPH_ALL);
    igraph_turan(&g_turan, &types_turan, 13, 4);

    printf("\nTuran graph with 13 vertices and 4 partitions:");
    print_graph_canon(&g_turan);

    printf("\nPartition type:\n");
    igraph_vector_int_print(&types_turan);

    igraph_isomorphic(&g_multipartite, &g_turan, &res);
    if (res) {
        printf("\nTuran(13,4) is isomorphic to full multipartite(4,3,3,3)\n");
    }

    igraph_vector_int_destroy(&partitions);
    igraph_vector_int_destroy(&types_multipartite);
    igraph_vector_int_destroy(&types_turan);
    igraph_destroy(&g_turan);
    igraph_destroy(&g_multipartite);


    igraph_vector_int_init(&partitions, 3);
    igraph_vector_int_init(&types_multipartite, 0);
    igraph_vector_int_init(&types_turan, 0);

    VECTOR(partitions)[0] = 3;
    VECTOR(partitions)[1] = 3;
    VECTOR(partitions)[2] = 2;

    igraph_full_multipartite(&g_multipartite, &types_multipartite, &partitions,
                        IGRAPH_UNDIRECTED, IGRAPH_ALL);
    igraph_turan(&g_turan, &types_turan, 8, 3);

    printf("\nTuran graph with 8 vertices and 3 partitions:");
    print_graph_canon(&g_turan);

    printf("\nPartition type:\n");
    igraph_vector_int_print(&types_turan);

    igraph_isomorphic(&g_multipartite, &g_turan, &res);
    if (res) {
        printf("\nTuran(8,3) is isomorphic to full multipartite(3,3,2)\n");
    }

    igraph_vector_int_destroy(&partitions);
    igraph_vector_int_destroy(&types_multipartite);
    igraph_vector_int_destroy(&types_turan);
    igraph_destroy(&g_turan);
    igraph_destroy(&g_multipartite);


    igraph_vector_int_init(&partitions, 3);
    igraph_vector_int_init(&types_multipartite, 0);
    igraph_vector_int_init(&types_turan, 0);

    VECTOR(partitions)[0] = 2;
    VECTOR(partitions)[1] = 2;
    VECTOR(partitions)[2] = 2;

    igraph_full_multipartite(&g_multipartite, &types_multipartite, &partitions,
                        IGRAPH_UNDIRECTED, IGRAPH_ALL);
    igraph_turan(&g_turan, &types_turan, 6, 3);

    printf("\nTuran graph with 6 vertices and 3 partitions:");
    print_graph_canon(&g_turan);

    printf("\nPartition type:\n");
    igraph_vector_int_print(&types_turan);

    igraph_isomorphic(&g_multipartite, &g_turan, &res);
    if (res) {
        printf("\nTuran(6,3) is isomorphic to full multipartite(2,2,2)\n");
    }

    igraph_vector_int_destroy(&partitions);
    igraph_vector_int_destroy(&types_multipartite);
    igraph_vector_int_destroy(&types_turan);
    igraph_destroy(&g_turan);
    igraph_destroy(&g_multipartite);

    VERIFY_FINALLY_STACK();

    return IGRAPH_SUCCESS;
}
