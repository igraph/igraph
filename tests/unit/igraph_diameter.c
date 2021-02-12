/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.inc"

int main() {

    igraph_t g;
    igraph_real_t result;
    igraph_integer_t from, to;
    igraph_vector_t path_edge, path_vertex, weights_vec;
    igraph_real_t weights[] = { 1, 2, 3, 4, 5, 1, 1, 1, 1};


    igraph_barabasi_game(&g, 30, /*power=*/ 1, 30, 0, 0, /*A=*/ 1,
                         IGRAPH_DIRECTED, IGRAPH_BARABASI_BAG,
                         /*start_from=*/ 0);
    igraph_diameter(&g, &result, 0, 0, 0, 0, IGRAPH_UNDIRECTED, 1);

    printf("Diameter: %li\n", (long int) result);

    igraph_destroy(&g);

    printf("diameter of ring and the path in terms of edges and vertices \n");

    igraph_ring(&g, 10, IGRAPH_DIRECTED, 0, 0);
    igraph_vector_init(&path_vertex, 0);
    igraph_vector_init(&path_edge, 0);
    igraph_vector_view(&weights_vec, weights, sizeof(weights) / sizeof(igraph_real_t));
    igraph_diameter(&g, &result, &from, &to, &path_vertex, &path_edge, IGRAPH_DIRECTED, 1);
    printf("diameter: %li, from %li to %li\n", (long int) result,
           (long int) from, (long int) to);
    print_vector_round(&path_vertex); 
    print_vector_round(&path_edge); 
    igraph_vector_destroy(&path_vertex);
    igraph_vector_destroy(&path_edge);

    printf("diameter of ring and the path in terms of edges with weights \n");

    igraph_vector_init(&path_edge, 0);
    igraph_diameter_dijkstra(&g, &weights_vec, &result, &from, &to, 0, &path_edge, IGRAPH_DIRECTED, 1);
    printf("diameter: %li, from %li to %li\n", (long int) result,
           (long int) from, (long int) to);
    print_vector_round(&path_edge);
    igraph_vector_destroy(&path_edge);

    printf("diameter of ring and the path in terms of vertices with weights \n");

    igraph_vector_init(&path_vertex, 0);
    igraph_diameter_dijkstra(&g, &weights_vec, &result, &from, &to, &path_vertex, 0, IGRAPH_DIRECTED, 1);
    printf("diameter: %li, from %li to %li\n", (long int) result,
           (long int) from, (long int) to);
    print_vector_round(&path_vertex);
    igraph_vector_destroy(&path_vertex);


    igraph_destroy(&g);

    return 0;
}
