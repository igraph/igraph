/*
   IGraph library.
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

#include "igraph_structural.h"
#include "igraph_interface.h"

igraph_error_t igraph_i_rich_club_density_sequence(const igraph_t *graph,
                                                  const igraph_vector_int_t *vertex_order,
                                                  igraph_bool_t directed,
                                                  igraph_bool_t loops,
                                                  igraph_vector_t *res) {
    /* TODO:

     * should we check that the vertex_order given is a permutation of the
       vertices of the graph? Or is it valid if vertex_order is only a subset of the vertices?

     * implement directed vs. undirected functionality

    */
    igraph_integer_t numVertices = igraph_vcount(graph);

    /* Error: vertex_order wrong size */
    if (igraph_vector_int_size(vertex_order) != numVertices) {
        IGRAPH_ERROR("Invalid vertex order length.", IGRAPH_EINVAL);
    }

    /* Make sure res is the right size (better way to do this?) */
    IGRAPH_CHECK(igraph_vector_resize(res, numVertices));

    printf("Checkpoint 1\n");

    /* create copy of the graph */
    igraph_t graphCopy;
    IGRAPH_CHECK(igraph_copy(&graphCopy, graph));

    /* loop through the number of vertices in the graph */
    for (int i = 0; i < numVertices; i++) {
        printf("Checkpoint 2a\n");

        /* compute the density */
        igraph_real_t density;
        IGRAPH_CHECK(igraph_density(&graphCopy, &density, loops));
        printf("Checkpoint 2b\n");

        /* store density in result */
        VECTOR(*res)[i] = density;
        printf("Checkpoint 2c\n");

        /* remove a vertex */
        IGRAPH_CHECK(igraph_delete_vertices(&graphCopy, igraph_vss_1(VECTOR(*vertex_order)[i])));
        printf("Checkpoint 2d\n");
    }

    /* destroy */
    igraph_destroy(&graphCopy);

    // TODO: add static back in

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_rich_club_coefficient(const igraph_t *graph,
                                            const igraph_vector_t *weights,
                                            igraph_neimode_t mode,
                                            igraph_vector_t *res) {
    /* Computes the rich-club coefficient as a function of degree. */

    return IGRAPH_SUCCESS;
}
