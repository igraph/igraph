/*
   igraph library.
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

int main(void) {
    igraph_t graph;
    igraph_vector_int_t dimvector;

    /* empty graph */
    IGRAPH_CHECK(igraph_vector_int_init(&dimvector, 2));
    VECTOR(dimvector)[0] = 3;

    IGRAPH_CHECK(igraph_hexagonal_lattice(&graph, &dimvector, true, false));
    printf("Empty graph:\n");
    print_graph_canon(&graph);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&dimvector);

    /* triangular triangular lattice with a single vertex and no edges*/
    IGRAPH_CHECK(igraph_vector_int_init(&dimvector, 1));
    VECTOR(dimvector)[0] = 1;

    IGRAPH_CHECK(igraph_hexagonal_lattice(&graph, &dimvector, true, false));
    printf("Triangular hexagonal lattice, single hexagon:\n");
    print_graph_canon(&graph);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&dimvector);

    /* triangular hexagonal lattice */
    IGRAPH_CHECK(igraph_vector_int_init(&dimvector, 1));
    VECTOR(dimvector)[0] = 5;

    IGRAPH_CHECK(igraph_hexagonal_lattice(&graph, &dimvector, true, false));
    printf("Triangular hexagonal lattice:\n");
    print_graph_canon(&graph);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&dimvector);

    /* rectangular hexagonal lattice */
    IGRAPH_CHECK(igraph_vector_int_init(&dimvector, 2));
    VECTOR(dimvector)[0] = 4;
    VECTOR(dimvector)[1] = 5;

    IGRAPH_CHECK(igraph_hexagonal_lattice(&graph, &dimvector, true, true));
    printf("Rectangular hexagonal lattice:\n");
    print_graph_canon(&graph);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&dimvector);

    /* hexagonal hexagonal lattice */
    IGRAPH_CHECK(igraph_vector_int_init(&dimvector, 3));
    VECTOR(dimvector)[0] = 3;
    VECTOR(dimvector)[1] = 4;
    VECTOR(dimvector)[2] = 5;

    IGRAPH_CHECK(igraph_hexagonal_lattice(&graph, &dimvector, false, true));
    printf("Hexagonal hexagonal lattice:\n");
    print_graph_canon(&graph);

    igraph_destroy(&graph);

    /* Erroneous calls */
    VECTOR(dimvector)[0] = -3;
    CHECK_ERROR(igraph_hexagonal_lattice(&graph, &dimvector, true, true), IGRAPH_EINVAL);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&dimvector);

    IGRAPH_CHECK(igraph_vector_int_init(&dimvector, 4));
    VECTOR(dimvector)[0] = 3;
    VECTOR(dimvector)[1] = 4;
    VECTOR(dimvector)[2] = 5;
    VECTOR(dimvector)[3] = 5;
    CHECK_ERROR(igraph_hexagonal_lattice(&graph, &dimvector, true, true), IGRAPH_EINVAL);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&dimvector);

    VERIFY_FINALLY_STACK();
    return 0;
}
