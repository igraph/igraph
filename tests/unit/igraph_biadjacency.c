/*
   igraph library.
   Copyright (C) 2022-2024  The igraph development team <igraph@igraph.org>

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

void print_and_destroy(igraph_matrix_t *biadjmat, igraph_bool_t directed, igraph_neimode_t mode, igraph_bool_t multiple) {
    igraph_t g;
    igraph_vector_bool_t types;
    igraph_vector_bool_init(&types, 0);

    igraph_biadjacency(&g, &types, biadjmat, directed, mode, multiple);

    print_graph_canon(&g);
    printf("types: ");

    igraph_vector_bool_print(&types);
    igraph_matrix_destroy(biadjmat);
    igraph_destroy(&g);
    igraph_vector_bool_destroy(&types);
}

int main(void) {
    igraph_t g;
    igraph_vector_bool_t types;
    igraph_vector_bool_init(&types, 0);
    igraph_matrix_t biadjmat;

    {
        printf("Bipartite adjacency matrix with no rows and no columns:\n");
        igraph_matrix_init(&biadjmat, 0, 0);
        print_and_destroy(&biadjmat, IGRAPH_DIRECTED, IGRAPH_ALL, false);
    }

    {
        printf("Bipartite adjacency matrix no rows and some columns:\n");
        igraph_matrix_init(&biadjmat, 0, 5);
        print_and_destroy(&biadjmat, IGRAPH_DIRECTED, IGRAPH_ALL, false);
    }

    {
        printf("\nBipartite adjacency matrix for two vertices:\n");
        int elem[] = {5};
        matrix_init_int_row_major(&biadjmat, 1, 1, elem);
        print_and_destroy(&biadjmat, IGRAPH_DIRECTED, IGRAPH_ALL, false);
    }

    {
        printf("\nBipartite adjacency matrix for two vertices, multiple = true:\n");
        int elem[] = {5};
        matrix_init_int_row_major(&biadjmat, 1, 1, elem);
        print_and_destroy(&biadjmat, IGRAPH_DIRECTED, IGRAPH_ALL, true);
    }

    {
        printf("\nBipartite adjacency matrix for five vertices:\n");
        int elem[] = {0, 1, 2, 3, 4, 5};
        matrix_init_int_row_major(&biadjmat, 2, 3, elem);
        print_and_destroy(&biadjmat, IGRAPH_DIRECTED, IGRAPH_ALL, true);
    }

    {
        printf("\nSame graph, IGRAPH_OUT:\n");
        int elem[] = {0, 1, 2, 3, 4, 5};
        matrix_init_int_row_major(&biadjmat, 2, 3, elem);
        print_and_destroy(&biadjmat, IGRAPH_DIRECTED, IGRAPH_OUT, true);
    }

    {
        printf("\nSame graph, IGRAPH_IN:\n");
        int elem[] = {0, 1, 2, 3, 4, 5};
        matrix_init_int_row_major(&biadjmat, 2, 3, elem);
        print_and_destroy(&biadjmat, IGRAPH_DIRECTED, IGRAPH_IN, true);
    }

    {
        printf("\nSame graph, undirected:\n");
        int elem[] = {0, 1, 2, 3, 4, 5};
        matrix_init_int_row_major(&biadjmat, 2, 3, elem);
        print_and_destroy(&biadjmat, IGRAPH_UNDIRECTED, IGRAPH_OUT, true);
    }

    {
        printf("\nNon-integer elements, multiple=false:\n");
        igraph_real_t elem[] = {0.0, 1.2, -1.8,
                                5.0, 0.0,  2.2};
        matrix_init_real_row_major(&biadjmat, 2, 3, elem);
        print_and_destroy(&biadjmat, IGRAPH_UNDIRECTED, IGRAPH_ALL, false);
    }

    {
        printf("\nNon-integer elements, multiple=true:\n");
        igraph_real_t elem[] = {0.0, 1.2, 1.8,
                                5.0, 0.0, 2.2};
        matrix_init_real_row_major(&biadjmat, 2, 3, elem);
        print_and_destroy(&biadjmat, IGRAPH_UNDIRECTED, IGRAPH_ALL, true);
    }

    VERIFY_FINALLY_STACK();

    {
        printf("\nCheck error for negative element.\n");
        int elem[] = {-5};
        matrix_init_int_row_major(&biadjmat, 1, 1, elem);
        CHECK_ERROR(igraph_biadjacency(&g, &types, &biadjmat, IGRAPH_DIRECTED,
                    IGRAPH_ALL, true), IGRAPH_EINVAL);
        igraph_matrix_destroy(&biadjmat);
    }

    igraph_vector_bool_destroy(&types);

    VERIFY_FINALLY_STACK();
    return 0;
}
