/*
   igraph library.
   Copyright (C) 2026 The igraph development team <igraph@igraph.org>

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
#include <stdlib.h>

#include "test_utilities.h"

int main(void) {

    igraph_t g;
    igraph_matrix_t coords;
    igraph_vector_int_t dimvector;

    igraph_matrix_init(&coords, 1, 1);

    /* Create 4-dimensional lattice dimvector for igraph_layout_square testing */
    igraph_vector_int_init(&dimvector, 4);
    VECTOR(dimvector)[0] = 4;
    VECTOR(dimvector)[1] = 3;
    VECTOR(dimvector)[2] = 2;
    VECTOR(dimvector)[3] = 2;

    /* 4-dimensional lattice when graph is initialized */
    igraph_empty(&g, 4*3*2*2, IGRAPH_UNDIRECTED);
    igraph_layout_square(&g, &coords, &dimvector);
    igraph_matrix_print(&coords);
    printf("=======\n");

    /* 4-dimensional lattice when graph is NULL */
    igraph_layout_square(NULL, &coords, &dimvector);
    igraph_matrix_print(&coords);
    printf("=======\n");

    igraph_destroy(&g);
    igraph_empty(&g, 4*3*2*1, IGRAPH_UNDIRECTED);
    printf("Check unequal number of nodes.\n");
    CHECK_ERROR(igraph_layout_square(&g, &coords, &dimvector), IGRAPH_EINVAL);
    igraph_destroy(&g);
    igraph_vector_int_destroy(&dimvector);

    igraph_vector_int_init(&dimvector, 2);
    VECTOR(dimvector)[0] = 4;
    VECTOR(dimvector)[1] = -1;
    printf("Check negative dimension.\n");
    CHECK_ERROR(igraph_layout_square(NULL, &coords, &dimvector), IGRAPH_EINVAL);

    igraph_matrix_destroy(&coords);
    igraph_vector_int_destroy(&dimvector);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
