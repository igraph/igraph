/*
   IGraph library.
   Copyright (C) 2008-2022  The igraph development team <igraph@igraph.org>

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

    igraph_t graph;
    igraph_vector_t weights;
    igraph_real_t weights_data[] = { 0, 2, 1, 0, 5, 2, 1, 1, 0, 2, 2, 8, 1, 1, 3, 1, 1, 4, 2, 1 };
    igraph_matrix_t res;
    igraph_real_t cutoff;

    igraph_small(&graph, 10, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 0, 3,    1, 2, 1, 4, 1, 5,
                 2, 3, 2, 6,          3, 2, 3, 6,
                 4, 5, 4, 7,          5, 6, 5, 8, 5, 9,
                 7, 5, 7, 8,          8, 9,
                 5, 2,
                 2, 1,
                 -1);

    igraph_matrix_init(&res, 0, 0);

    printf("Unweighted distances:\n\n");

    igraph_distances(&graph, &res, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT);
    igraph_matrix_print(&res);

    cutoff = 3; /* distances longer than this will be returned as infinity */
    printf("\nUnweighted distances with a cutoff of %g:\n\n", cutoff);
    igraph_distances_cutoff(&graph, &res, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT, cutoff);
    igraph_matrix_print(&res);

    printf("\nWeighted distances:\n\n");

    igraph_vector_view(&weights, weights_data,
                       sizeof(weights_data) / sizeof(weights_data[0]));

    igraph_distances_dijkstra(&graph, &res, igraph_vss_all(), igraph_vss_all(),
                              &weights, IGRAPH_OUT);
    igraph_matrix_print(&res);

    cutoff = 8; /* distances longer than this will be returned as infinity */
    printf("\nWeighted distances with a cutoff of %g:\n\n", cutoff);
    igraph_distances_dijkstra_cutoff(&graph, &res, igraph_vss_all(), igraph_vss_all(),
                              &weights, IGRAPH_OUT, cutoff);
    igraph_matrix_print(&res);

    igraph_matrix_destroy(&res);
    igraph_destroy(&graph);

    return 0;
}
