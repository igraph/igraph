/*
   igraph library.
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

#include <igraph.h>

#include "test_utilities.h"

#include <math.h>

void check_align(const igraph_t *graph) {
    const igraph_int_t vcount = igraph_vcount(graph);
    igraph_matrix_t layout;
    igraph_bool_t connected;

    igraph_matrix_init(&layout, 0, 0);

    /* Test with a 2D layout. */

    igraph_layout_fruchterman_reingold(
        graph, &layout,
        /* use_seed */ false, /* niter */ 1000,
        /* start_temp */ sqrt(vcount),
        IGRAPH_LAYOUT_NOGRID,
        /* weights */ NULL,
        NULL, NULL, NULL, NULL);

    igraph_layout_align(graph, &layout);

    IGRAPH_ASSERT(igraph_matrix_nrow(&layout) == vcount);
    IGRAPH_ASSERT(igraph_matrix_ncol(&layout) == 2);
    IGRAPH_ASSERT(igraph_vector_is_all_finite(&layout.data));

    /* Test with a 3D layout. */

    igraph_layout_fruchterman_reingold_3d(
        graph, &layout,
        /* use_seed */ false, /* niter */ 1000,
        /* start_temp */ sqrt(vcount),
        /* weights */ NULL,
        NULL, NULL, NULL, NULL, NULL, NULL);

    igraph_layout_align(graph, &layout);

    IGRAPH_ASSERT(igraph_matrix_nrow(&layout) == vcount);
    IGRAPH_ASSERT(igraph_matrix_ncol(&layout) == 3);
    IGRAPH_ASSERT(igraph_vector_is_all_finite(&layout.data));

    /* Test connected graphs having at least 4 vertices with a 4D layout.
     * igraph_layout_mds() only supports 2D layouts with disconnected graphs,
     * and requires at least as many vertices as the dimension. */

    igraph_is_connected(graph, &connected, IGRAPH_WEAK);

    if (connected && vcount >= 4) {
        igraph_layout_mds(graph, &layout, NULL, 4);

        igraph_layout_align(graph, &layout);

        IGRAPH_ASSERT(igraph_matrix_nrow(&layout) == vcount);
        IGRAPH_ASSERT(igraph_matrix_ncol(&layout) == 4);
        IGRAPH_ASSERT(igraph_vector_is_all_finite(&layout.data));
    }

    /* Test with a 1D layout. */

    igraph_matrix_resize(&layout, vcount, 1);
    igraph_layout_align(graph, &layout);

    igraph_matrix_destroy(&layout);

    VERIFY_FINALLY_STACK();
}

int main(void) {
    igraph_t graph;

    printf("Empty\n");
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    check_align(&graph);
    igraph_destroy(&graph);

    printf("Singleton\n");
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    check_align(&graph);
    igraph_destroy(&graph);

    printf("Singleton with loop\n");
    igraph_full(&graph, 1, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    check_align(&graph);
    igraph_destroy(&graph);

    printf("Two isolated vertices\n");
    igraph_empty(&graph, 2, IGRAPH_UNDIRECTED);
    check_align(&graph);
    igraph_destroy(&graph);

    printf("P_2\n");
    igraph_full(&graph, 2, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    check_align(&graph);
    igraph_destroy(&graph);

    printf("K_2 with loops\n");
    igraph_full(&graph, 2, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    check_align(&graph);
    igraph_destroy(&graph);

    printf("C_3\n");
    igraph_full(&graph, 3, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    check_align(&graph);
    igraph_destroy(&graph);

    printf("K_4\n");
    igraph_full(&graph, 4, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    check_align(&graph);
    igraph_destroy(&graph);

    printf("K_5\n");
    igraph_full(&graph, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    check_align(&graph);
    igraph_destroy(&graph);

    printf("C_4\n");
    igraph_ring(&graph, 4, IGRAPH_UNDIRECTED, false, /* circular */ true);
    check_align(&graph);
    igraph_destroy(&graph);

    printf("P_4\n");
    igraph_ring(&graph, 4, IGRAPH_UNDIRECTED, false, /* circular */ false);
    check_align(&graph);
    igraph_destroy(&graph);

    printf("C_7\n");
    igraph_ring(&graph, 7, IGRAPH_UNDIRECTED, false, /* circular */ true);
    check_align(&graph);
    igraph_destroy(&graph);

    printf("Frucht\n");
    igraph_famous(&graph, "Frucht");
    check_align(&graph);
    igraph_destroy(&graph);

    printf("Petersen\n");
    igraph_famous(&graph, "Petersen");
    check_align(&graph);
    igraph_destroy(&graph);

    printf("Groetzsch\n");
    igraph_famous(&graph, "Groetzsch");
    check_align(&graph);
    igraph_destroy(&graph);

    printf("Kautz(3,2)\n");
    igraph_kautz(&graph, 3, 2);
    check_align(&graph);
    igraph_destroy(&graph);

    printf("Karate club\n");
    igraph_famous(&graph, "Zachary");
    check_align(&graph);
    igraph_destroy(&graph);

    {
        igraph_matrix_t layout;

        /* Wrong row count */
        igraph_matrix_init(&layout, 5, 2);
        igraph_full(&graph, 3, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
        CHECK_ERROR(igraph_layout_align(&graph, &layout), IGRAPH_EINVAL);

        /* Zero-dimensional layout */
        igraph_matrix_resize(&layout, igraph_vcount(&graph), 0);
        CHECK_ERROR(igraph_layout_align(&graph, &layout), IGRAPH_EINVAL);

        igraph_destroy(&graph);
        igraph_matrix_destroy(&layout);
    }

    return 0;
}
