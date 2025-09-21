/* igraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

void call_and_print(igraph_t *graph, igraph_vs_t vids, igraph_vector_t *weights) {
    igraph_vector_t result;
    igraph_vector_init(&result, 0);
    IGRAPH_ASSERT(igraph_constraint(graph, &result, vids, weights) == IGRAPH_SUCCESS);
    print_vector(&result);
    igraph_vector_destroy(&result);
    printf("\n");
}


int main(void) {
    igraph_t g_0, g_1, g_4_full, g_4_full_split, g_4_full_loop, g_4_line, g_hole;
    igraph_vector_t weights_full, weights_full_split, weights_line, result;

    igraph_vector_init_real(&weights_line, 6, 1., 0., 0., 1., 0., 1.);
    igraph_vector_init_real(&weights_full, 6, 1., 1., 1., 1., 1., 1.);
    igraph_vector_init_real(&weights_full_split, 9, .5, .5, 1., 1., 1., 1./3., 1./3., 1./3., 1.);
    igraph_vector_init(&result, 0);
    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_4_full, 4, IGRAPH_UNDIRECTED, 0,1, 0,2, 0,3, 1,2, 1,3, 2,3, -1);
    igraph_small(&g_4_full_loop, 4, IGRAPH_UNDIRECTED, 0,0, 0,1, 0,2, 0,3, 1,2, 1,3, 2,3, -1);
    igraph_small(&g_4_full_split, 4, IGRAPH_UNDIRECTED, 0,1, 0,1, 0,2, 0,3, 1,2,
            1,3, 1,3, 1,3, 2,3, -1);
    igraph_small(&g_4_line, 4, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, -1);
    igraph_small(&g_hole, 9, IGRAPH_UNDIRECTED, 0,1, 0,2, 0,3, 1,2, 1,3, 2,3, 3,4, 4,5,
            5,6, 5,7, 5,8, 6,7, 6,8, 7,8, -1);

    printf("No vertices:\n");
    call_and_print(&g_0, igraph_vss_none(), NULL);

    printf("One vertex:\n");
    call_and_print(&g_1, igraph_vss_1(0), NULL);

    printf("Line graph, 4 vertices:\n");
    call_and_print(&g_4_line, igraph_vss_all(), NULL);

    printf("Full graph, 4 vertices, weights make it a line:\n");
    call_and_print(&g_4_full, igraph_vss_all(), &weights_line);

    /*
       each p_ij and p_ji is equal to 1/3, because each node has 3 connections with
       equal weight and p_ij is z_ij/sum_j(z_ij) where z are the edge weights.
       so C_ij = (1/3 + (1/3)^2) for each j, which means C_i = (1/3 + 2/9)^2 *3
       = 0.925925925925926
       */
    printf("Full graph, 4 vertices, all same weights:\n");
    call_and_print(&g_4_full, igraph_vss_all(), &weights_full);

    printf("Full graph, 4 vertices, all same weights, but split over multiple edges:\n");
    call_and_print(&g_4_full_split, igraph_vss_all(), &weights_full_split);

    printf("Full graph, 4 vertices, no weights:\n");
    call_and_print(&g_4_full, igraph_vss_all(), NULL);

    printf("Full graph, 4 vertices, no weights, with loop:\n");
    call_and_print(&g_4_full_loop, igraph_vss_all(), NULL);

    /*
       for node 0, each p_ij is again equal to 1/3, but p_ji equals 1/4 for j = 3,
       because that's connected to 4 other nodes, so for example for C_01 the contribution
       is
       C_01 = (p_01 + p_02p_21 + p_03p_31)^2 = (1/3 + 1/9 + 1/12)^2
       and you end up with:
       ((1/3 + 2/9)^2 ) + ((1/3 + 1/9 + 1/12)^2 * 2) = 0.8657407407407408
       */
    printf("Hole in middle of two clusters:\n");
    call_and_print(&g_hole, igraph_vss_all(), NULL);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_4_full);
    igraph_destroy(&g_4_full_loop);
    igraph_destroy(&g_4_full_split);
    igraph_destroy(&g_4_line);
    igraph_destroy(&g_hole);
    igraph_vector_destroy(&weights_line);
    igraph_vector_destroy(&weights_full);
    igraph_vector_destroy(&weights_full_split);
    igraph_vector_destroy(&result);

    VERIFY_FINALLY_STACK();
    return 0;
}
