/* IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

void print_destroy(igraph_matrix_t *adjmatrix, igraph_adjacency_t mode) {
    igraph_t g;

    igraph_adjacency(&g, adjmatrix, mode);
    print_graph_canon(&g);
    igraph_matrix_destroy(adjmatrix);
    igraph_destroy(&g);
}

int main() {
    igraph_matrix_t adjmatrix;
    igraph_t g;

    printf("\n0x0 matrix:\n");
    matrix_init_int_row_major(&adjmatrix, 0, 0, NULL);
    print_destroy(&adjmatrix, IGRAPH_ADJ_DIRECTED);

    printf("\n1x1 matrix:\n");
    {
        int e[] = {1};
        matrix_init_int_row_major(&adjmatrix, 1, 1, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_DIRECTED);
    }
    printf("\n3x3 matrix, IGRAP_ADJ_DIRECTED:\n");
    {
        int e[] = {1, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_DIRECTED);
    }
    printf("\n3x3 matrix, IGRAP_ADJ_UNDIRECTED:\n");
    {
        int e[] = {1, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_UNDIRECTED);
    }
    printf("\n3x3 matrix, IGRAP_ADJ_MAX:\n");
    {
        int e[] = {1, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_MAX);
    }
    printf("\n3x3 matrix, IGRAP_ADJ_MIN:\n");
    {
        int e[] = {1, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_MIN);
    }
    printf("\n3x3 matrix, IGRAP_ADJ_PLUS:\n");
    {
        int e[] = {1, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_PLUS);
    }
    printf("\n3x3 matrix, IGRAP_ADJ_UPPER:\n");
    {
        int e[] = {1, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_UPPER);
    }
    printf("\n3x3 matrix, IGRAP_ADJ_LOWER:\n");
    {
        int e[] = {1, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_LOWER);
    }
    VERIFY_FINALLY_STACK();

    printf("\nCheck handling of non-square matrix error.\n");
    {
        int e[] = {1, 2, 0};
        matrix_init_int_row_major(&adjmatrix, 3, 1, e);
        CHECK_ERROR(igraph_adjacency(&g, &adjmatrix, IGRAPH_ADJ_DIRECTED), IGRAPH_NONSQUARE);
        igraph_matrix_destroy(&adjmatrix);
    }
    printf("\nCheck handling of negative number of edges error.\n");
    {
        int e[] = {1, 2, 0, -3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        CHECK_ERROR(igraph_adjacency(&g, &adjmatrix, IGRAPH_ADJ_DIRECTED), IGRAPH_EINVAL);
        igraph_matrix_destroy(&adjmatrix);
    }

    return 0;
}
