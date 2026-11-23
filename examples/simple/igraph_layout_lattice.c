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

int main(void) {
    igraph_matrix_t res_m;
    igraph_vector_int_t dimvector;
    igraph_t graph;
    igraph_vector_bool_t periodic;

    /* Initialize the library. */
    igraph_setup();

    /*
    Initialize matrix and set it to 1 row and 1 col -
    igraph_layout_square function will resize it anyway
    */
    igraph_matrix_init(&res_m, 1, 1);

    /* Pick dimensions for the lattice */
    igraph_vector_int_init(&dimvector, 4);
    VECTOR(dimvector)[0] = 4;
    VECTOR(dimvector)[1] = 3;
    VECTOR(dimvector)[2] = 2;
    VECTOR(dimvector)[3] = 2;

    igraph_layout_square(NULL, &res_m, &dimvector);

    IGRAPH_ASSERT(igraph_layout_square(NULL, &res_m, &dimvector) == IGRAPH_SUCCESS);
    printf("Square lattice layout coordinates:\n");
    igraph_matrix_print(&res_m);


    /* Get layout coordinates on an initialized square lattice graph */ 
    igraph_vector_bool_init(&periodic, 4);
    igraph_vector_bool_fill(&periodic, false);
    igraph_square_lattice(
        &graph, &dimvector, /* nei= */ 0, IGRAPH_UNDIRECTED,
        /* mutual= */ false, &periodic
    );

    IGRAPH_ASSERT(igraph_layout_square(&graph, &res_m, &dimvector) == IGRAPH_SUCCESS);
    printf("Square lattice layout coordinates:\n");
    igraph_matrix_print(&res_m);

    igraph_matrix_destroy(&res_m);
    igraph_vector_int_destroy(&dimvector);
    igraph_vector_bool_destroy(&periodic);
    igraph_destroy(&graph);
    return 0;
}

