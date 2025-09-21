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

int main(void) {
    igraph_t graph;
    igraph_vector_t weights;
    igraph_matrix_int_t merges;

    /* Initialize the library. */
    igraph_setup();

    igraph_matrix_int_init(&merges, 0, 0);
    igraph_vector_init_int(&weights, 8, 10, 10, 1, 1, 1, 1, 1, 1);

    igraph_small(&graph, 6, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,3, 2,4, 2,5, 3,4, 3,5, 4,5, -1);
    igraph_community_fastgreedy(&graph, &weights, &merges,
                                /*modularity*/ NULL,
                                /*membership=*/ NULL);
    igraph_matrix_int_print(&merges);
    igraph_destroy(&graph);
    igraph_vector_destroy(&weights);
    igraph_matrix_int_destroy(&merges);

    return 0;
}
