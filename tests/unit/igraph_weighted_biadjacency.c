/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

void print_and_destroy(igraph_matrix_t *biadjmat, igraph_bool_t directed, igraph_neimode_t mode) {
    igraph_t g;
    igraph_vector_bool_t types;
    igraph_vector_t weights;

    igraph_vector_bool_init(&types, 0);
    igraph_vector_init(&weights, 0);

    igraph_weighted_biadjacency(&g, &types, &weights, biadjmat, directed, mode);

    print_weighted_graph_canon(&g, &weights);
    printf("types: ");
    igraph_vector_bool_print(&types);

    igraph_destroy(&g);
    igraph_vector_destroy(&weights);
    igraph_vector_bool_destroy(&types);

    igraph_matrix_destroy(biadjmat);
}

int main(void) {
    igraph_matrix_t biadjmat;

    {
        printf("Bipartite adjacency matrix with no rows and no columns:\n");
        igraph_matrix_init(&biadjmat, 0, 0);
        print_and_destroy(&biadjmat, IGRAPH_DIRECTED, IGRAPH_ALL);
    }

    {
        printf("Bipartite adjacency matrix with some rows and no columns:\n");
        igraph_matrix_init(&biadjmat, 3, 0);
        print_and_destroy(&biadjmat, IGRAPH_DIRECTED, IGRAPH_ALL);
    }

    {
        printf("\nBipartite adjacency matrix for two vertices:\n");
        igraph_real_t elem[] = {1.25};
        matrix_init_real_row_major(&biadjmat, 1, 1, elem);
        print_and_destroy(&biadjmat, IGRAPH_DIRECTED, IGRAPH_ALL);
    }

    {
        printf("\nBipartite adjacency matrix for five vertices:\n");
        igraph_real_t elem[] = {0.0, -4.5, 2.3,
                                -0.1, 0.0, 0.0};
        matrix_init_real_row_major(&biadjmat, 2, 3, elem);
        print_and_destroy(&biadjmat, IGRAPH_DIRECTED, IGRAPH_ALL);
    }

    {
        printf("\nSame graph, IGRAPH_OUT:\n");
        igraph_real_t elem[] = {0.0, -4.5, 2.3,
                                -0.1, 0.0, 0.0};
        matrix_init_real_row_major(&biadjmat, 2, 3, elem);
        print_and_destroy(&biadjmat, IGRAPH_DIRECTED, IGRAPH_OUT);
    }

    {
        printf("\nSame graph, IGRAPH_IN:\n");
        igraph_real_t elem[] = {0.0, -4.5, 2.3,
                                -0.1, 0.0, 0.0};
        matrix_init_real_row_major(&biadjmat, 2, 3, elem);
        print_and_destroy(&biadjmat, IGRAPH_DIRECTED, IGRAPH_IN);
    }

    {
        printf("\nSame graph, undirected:\n");
        igraph_real_t elem[] = {0.0, -4.5, 2.3,
                                -0.1, 0.0, 0.0};
        matrix_init_real_row_major(&biadjmat, 2, 3, elem);
        print_and_destroy(&biadjmat, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    }

    {
        printf("\nInfinity and NaN:\n");
        igraph_real_t elem[] = {0.0, -4.5,
                                IGRAPH_INFINITY, -IGRAPH_INFINITY,
                                IGRAPH_NAN, 0.0};
        matrix_init_real_row_major(&biadjmat, 3, 2, elem);
        print_and_destroy(&biadjmat, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    }

    VERIFY_FINALLY_STACK();
    return 0;
}
