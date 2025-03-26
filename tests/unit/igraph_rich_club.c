/*
   IGraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

#include <stdio.h>
#include <igraph.h>
#include "test_utilities.h"

int main(void) {
    igraph_t graph;
    igraph_vector_t res;
    igraph_vector_int_t vertex_order;
    const igraph_integer_t N_VERTICES = 6;

    /* GRAPH FORMAT
       1-4 0-3-2-5
     */
    IGRAPH_CHECK(igraph_small(&graph, N_VERTICES, IGRAPH_UNDIRECTED,
                              1, 4, 2, 3, 0, 3, 2, 5, -1));

    igraph_vector_init(&res, 0);
    igraph_vector_int_init(&vertex_order, N_VERTICES);

    // vertex removal order: 0, 1, 2, 3, 4, 5
    for (int i = 0; i < N_VERTICES; i++){
        VECTOR(vertex_order)[i] = i;
    }

    IGRAPH_CHECK(igraph_rich_club_density_sequence(&graph, &vertex_order, 0, 0, &res));

    printf("Test 3 (disjoint graph):\n");
    for (int i = 0; i < igraph_vector_size(&res); i++){
        printf("%.4f ", VECTOR(res)[i]);
    }
    printf("\n");

    /* EXPECTED OUTPUT
       0.2667 0.3000 0.3333 0.0000 0.0000 nan
    */
    igraph_vector_int_destroy(&vertex_order);
    igraph_vector_destroy(&res);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
    return 0;
}
