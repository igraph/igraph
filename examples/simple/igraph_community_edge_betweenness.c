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

#include <igraph.h>

int main(void) {
    igraph_t graph;
    igraph_vector_int_t edges;
    igraph_vector_t eb, weights;
    igraph_matrix_int_t merges;
    igraph_vector_int_t bridges;

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 0,1, 0,1, 0,1, -1);
    igraph_vector_init_int(&weights, 3, 1, 2, 3);
    igraph_vector_init(&eb, 0);
    igraph_vector_int_init(&edges, 0);
    igraph_matrix_int_init(&merges, 0, 0);
    igraph_vector_int_init(&bridges, 0);
    igraph_community_edge_betweenness(&graph, &edges, &eb, &merges,
                                      &bridges, /*modularity*/ NULL,
                                      /*membership*/ NULL,
                                      IGRAPH_UNDIRECTED,
                                      &weights);
    printf("edges:\n");
    igraph_vector_int_print(&edges);
    printf("edge betweenness:\n");
    igraph_vector_print(&eb);
    printf("merges:\n");
    igraph_matrix_int_print(&merges);
    printf("bridges:\n");
    igraph_vector_int_print(&bridges);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&eb);
    igraph_vector_destroy(&weights);
    igraph_matrix_int_destroy(&merges);
    igraph_vector_int_destroy(&bridges);

    return 0;
}
