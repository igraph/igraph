/*
   igraph library.
   Copyright (C) 2026 The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#include <igraph.h>

int main(void) {

    igraph_matrix_t res_m;
    igraph_vector_int_t dims;
    igraph_t graph;

    /* Initialize the library. */
    igraph_setup();

    /*
    Initialize matrix and set it to 1 row and 1 col -
    igraph_layout_square function will resize it anyway
    */
    igraph_matrix_init(&res_m, 1, 1);
    /* Pick dimensions for the lattice
    dims vector of length 2 is rectangular shape */
    igraph_vector_int_init(&dims, 2);
    VECTOR(dims)[0] = 4;
    VECTOR(dims)[1] = 3;


    IGRAPH_ASSERT(igraph_layout_triangular(NULL, &res_m, &dims) == IGRAPH_SUCCESS);
    printf("Triangular lattice layout coordinates:\n");
    igraph_matrix_print(&res_m);

    printf("--------------------------\n");
    igraph_triangular_lattice(
        &graph,
        &dims,
        IGRAPH_UNDIRECTED,
        /* mutual= */ false
    );

    IGRAPH_ASSERT(igraph_layout_triangular(&graph, &res_m, &dims) == IGRAPH_SUCCESS);
    printf("Triangular lattice layout coordinates:\n");
    igraph_matrix_print(&res_m);

    printf("triangle shape ----------------\n");
    igraph_vector_int_destroy(&dims);
    igraph_vector_int_init(&dims, 1);
    VECTOR(dims)[0] = 4;
    IGRAPH_ASSERT(igraph_layout_triangular(NULL, &res_m, &dims) == IGRAPH_SUCCESS);
    igraph_matrix_print(&res_m);

    igraph_matrix_destroy(&res_m);
    igraph_vector_int_destroy(&dims);
    igraph_destroy(&graph);
    return 0;
}

