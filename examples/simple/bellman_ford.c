/*
   igraph library.
   Copyright (C) 2008-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

int main(void) {

    igraph_t g;
    igraph_vector_t weights;
    igraph_real_t weights_data_0[] = { 0, 2, 1, 0, 5, 2, 1, 1, 0, 2, 2, 8, 1, 1, 3, 1, 1, 4, 2, 1 };
    igraph_real_t weights_data_1[] = { 6, 7, 8, -4, -2, -3, 9, 2, 7 };
    igraph_real_t weights_data_2[] = { 6, 7, 2, -4, -2, -3, 9, 2, 7 };
    igraph_matrix_t res;

    /* Initialize the library. */
    igraph_setup();

    /* Graph with only positive weights */
    igraph_small(&g, 10, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 0, 3,    1, 2, 1, 4, 1, 5,
                 2, 3, 2, 6,         3, 2, 3, 6,
                 4, 5, 4, 7,         5, 6, 5, 8, 5, 9,
                 7, 5, 7, 8,         8, 9,
                 5, 2,
                 2, 1,
                 -1);

    weights = igraph_vector_view(weights_data_0,
                       sizeof(weights_data_0) / sizeof(weights_data_0[0]));

    igraph_matrix_init(&res, 0, 0);
    igraph_distances_bellman_ford(&g, &res, igraph_vss_all(), igraph_vss_all(),
                                  &weights, IGRAPH_OUT);
    igraph_matrix_print(&res);

    igraph_matrix_destroy(&res);
    igraph_destroy(&g);

    printf("\n");

    /***************************************/

    /* Graph with negative weights */
    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 0, 1, 0, 3, 1, 3, 1, 4, 2, 1, 3, 2, 3, 4, 4, 0, 4, 2, -1);

    weights = igraph_vector_view(weights_data_1,
                       sizeof(weights_data_1) / sizeof(weights_data_1[0]));

    igraph_matrix_init(&res, 0, 0);
    igraph_distances_bellman_ford(&g, &res, igraph_vss_all(),
                                  igraph_vss_all(), &weights, IGRAPH_OUT);
    igraph_matrix_print(&res);

    /***************************************/

    /* Same graph with negative cycle */
    igraph_set_error_handler(igraph_error_handler_ignore);
    weights = igraph_vector_view(weights_data_2,
                       sizeof(weights_data_2) / sizeof(weights_data_2[0]));
    if (igraph_distances_bellman_ford(&g, &res, igraph_vss_all(),
                                      igraph_vss_all(),
                                      &weights, IGRAPH_OUT) != IGRAPH_ENEGCYCLE) {
        return 1;
    }

    igraph_matrix_destroy(&res);
    igraph_destroy(&g);

    if (!IGRAPH_FINALLY_STACK_EMPTY) {
        return 1;
    }

    return 0;
}
