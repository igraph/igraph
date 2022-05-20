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

void print_destroy(igraph_matrix_t *adjmatrix, igraph_adjacency_t mode, igraph_loops_t loops) {
    igraph_t g;

    igraph_adjacency(&g, adjmatrix, mode, loops);
    print_graph_canon(&g);
    igraph_matrix_destroy(adjmatrix);
    igraph_destroy(&g);
}

int main() {
    igraph_matrix_t adjmatrix;
    igraph_t g;

    printf("\n0x0 matrix:\n");
    matrix_init_int_row_major(&adjmatrix, 0, 0, NULL);
    print_destroy(&adjmatrix, IGRAPH_ADJ_DIRECTED, IGRAPH_LOOPS_ONCE);

    printf("\n1x1 matrix, no loops:\n");
    {
        int e[] = {1};
        matrix_init_int_row_major(&adjmatrix, 1, 1, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_DIRECTED, IGRAPH_NO_LOOPS);
    }
    printf("\n1x1 matrix, loops once:\n");
    {
        int e[] = {1};
        matrix_init_int_row_major(&adjmatrix, 1, 1, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_DIRECTED, IGRAPH_LOOPS_ONCE);
    }
    printf("\n1x1 matrix, loops twice:\n");
    {
        int e[] = {2};
        matrix_init_int_row_major(&adjmatrix, 1, 1, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_DIRECTED, IGRAPH_LOOPS_TWICE);
    }

    printf("\n3x3 matrix, IGRAPH_ADJ_DIRECTED, no loops:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_DIRECTED, IGRAPH_NO_LOOPS);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_DIRECTED, loops once:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_DIRECTED, IGRAPH_LOOPS_ONCE);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_DIRECTED, loops twice:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_DIRECTED, IGRAPH_LOOPS_TWICE);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_UNDIRECTED, no loops:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_UNDIRECTED, IGRAPH_NO_LOOPS);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_UNDIRECTED, loops once:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_UNDIRECTED, IGRAPH_LOOPS_ONCE);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_UNDIRECTED, loops twice:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_UNDIRECTED, IGRAPH_LOOPS_TWICE);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_MAX, no loops:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_MAX, IGRAPH_NO_LOOPS);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_MAX, loops once:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_MAX, IGRAPH_LOOPS_ONCE);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_MAX, loops twice:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_MAX, IGRAPH_LOOPS_TWICE);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_MIN, no loops:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 5, 0, 4, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_MIN, IGRAPH_NO_LOOPS);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_MIN, loops once:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_MIN, IGRAPH_LOOPS_ONCE);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_MIN, loops twice:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_MIN, IGRAPH_LOOPS_TWICE);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_PLUS, no loops:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_PLUS, IGRAPH_NO_LOOPS);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_PLUS, loops once:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_PLUS, IGRAPH_LOOPS_ONCE);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_PLUS, loops twice:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_PLUS, IGRAPH_LOOPS_TWICE);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_UPPER, no loops:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_UPPER, IGRAPH_NO_LOOPS);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_UPPER, loops once:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_UPPER, IGRAPH_LOOPS_ONCE);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_UPPER, loops twice:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_UPPER, IGRAPH_LOOPS_TWICE);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_LOWER, no loops:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_LOWER, IGRAPH_NO_LOOPS);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_LOWER, loops once:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_LOWER, IGRAPH_LOOPS_ONCE);
    }
    printf("\n3x3 matrix, IGRAPH_ADJ_LOWER, loops twice:\n");
    {
        int e[] = {4, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        print_destroy(&adjmatrix, IGRAPH_ADJ_LOWER, IGRAPH_LOOPS_TWICE);
    }

    VERIFY_FINALLY_STACK();

    printf("\nCheck handling of non-square matrix error.\n");
    {
        int e[] = {1, 2, 0};
        matrix_init_int_row_major(&adjmatrix, 3, 1, e);
        CHECK_ERROR(igraph_adjacency(&g, &adjmatrix, IGRAPH_ADJ_DIRECTED, IGRAPH_NO_LOOPS), IGRAPH_NONSQUARE);
        igraph_matrix_destroy(&adjmatrix);
    }
    printf("\nCheck handling of negative number of edges error.\n");
    {
        int e[] = {1, 2, 0, -3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        CHECK_ERROR(igraph_adjacency(&g, &adjmatrix, IGRAPH_ADJ_DIRECTED, IGRAPH_NO_LOOPS), IGRAPH_EINVAL);
        igraph_matrix_destroy(&adjmatrix);
    }
    printf("\nCheck handling of odd number in diagonal.\n");
    {
        int e[] = {1, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        CHECK_ERROR(igraph_adjacency(&g, &adjmatrix, IGRAPH_ADJ_DIRECTED, IGRAPH_LOOPS_TWICE), IGRAPH_EINVAL);
        igraph_matrix_destroy(&adjmatrix);
    }
    printf("\nCheck handling of invalid adjacency mode.\n");
    {
        int e[] = {0, 2, 0, 3, 0, 4, 0, 5, 6};
        matrix_init_int_row_major(&adjmatrix, 3, 3, e);
        CHECK_ERROR(igraph_adjacency(&g, &adjmatrix, 42, IGRAPH_LOOPS_TWICE), IGRAPH_EINVAL);
        igraph_matrix_destroy(&adjmatrix);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
