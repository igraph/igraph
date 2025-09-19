/*
   igraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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
#include <math.h>
#include <stdlib.h>

#include "test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_matrix_t coords;
    igraph_vector_int_t edgelist, edgelist2;
    igraph_matrix_list_t routing;
    igraph_vector_int_t layers;

    igraph_matrix_init(&coords, 0, 0);
    igraph_matrix_list_init(&routing, 0);

    /* Layout on simple graph with predefined layers */
    printf("Simple graph with predefined layers\n");
    igraph_vector_int_init_int_end(&layers, -1, 0, 1, 1, 2, 3, 3, 4, 4, 5, -1);
    igraph_vector_int_init_int_end(&edgelist, -1,
                               0, 1, 0, 2, 0, 3, 1, 2, 2, 2, 1, 4, 2, 5, 4, 6, 5, 7, 6, 8, 7, 8,
                               3, 8, 8, 1, 8, 2, -1);
    igraph_create(&g, &edgelist, 0, IGRAPH_DIRECTED);

    igraph_layout_sugiyama(&g, &coords, 0, &layers,
                           /* hgap = */ 1,
                           /* vgap = */ 1,
                           /* maxiter = */ 100,
                           /* weights = */ 0);
    print_matrix(&coords);
    printf("\n");

    /* Same, but this time also return the routing information */
    printf("Simple graph with routing information\n");
    igraph_layout_sugiyama(&g, &coords, &routing, &layers,
                           /* hgap = */ 1,
                           /* vgap = */ 1,
                           /* maxiter = */ 100,
                           /* weights = */ 0);
    print_matrix(&coords);
    print_matrix_list(&routing);
    printf("\n");

    igraph_vector_int_destroy(&layers);

    /* Same, but with automatic layering */
    printf("Simple graph with automatic layers\n");
    igraph_layout_sugiyama(&g, &coords, &routing, 0,
                           /* hgap = */ 1,
                           /* vgap = */ 1,
                           /* maxiter = */ 100,
                           /* weights = */ 0);
    print_matrix(&coords);
    print_matrix_list(&routing);
    printf("\n");

    /* Layering with gaps in it */
    printf("Simple graph with gaps in layers\n");
    igraph_vector_int_init_int_end(&layers, -1, 0, 2, 2, 4, 6, 6, 12, 12, 15, -1);
    igraph_layout_sugiyama(&g, &coords, &routing, &layers,
                           /* hgap = */ 1,
                           /* vgap = */ 1,
                           /* maxiter = */ 100,
                           /* weights = */ 0);
    print_matrix(&coords);
    print_matrix_list(&routing);
    printf("\n");

    igraph_destroy(&g);
    igraph_vector_int_destroy(&edgelist);
    igraph_vector_int_destroy(&layers);

    /* Check edge directions and the order of control points within edges */
    igraph_vector_int_init_int_end(&edgelist, -1, 0, 1, 1, 2, 2, 3, 0, 3, 3, 0, -1);
    igraph_vector_int_init_int_end(&layers, -1, 0, 1, 2, 3, -1);

    printf("Path graph with two shortcuts, directed\n");
    igraph_create(&g, &edgelist, 0, IGRAPH_DIRECTED);
    igraph_vector_int_init(&edgelist2, 0);
    igraph_get_edgelist(&g, &edgelist2, 0);
    print_vector_int(&edgelist2);
    igraph_vector_int_destroy(&edgelist2);
    igraph_layout_sugiyama(&g, &coords, &routing, &layers,
                           /* hgap = */ 1,
                           /* vgap = */ 1,
                           /* maxiter = */ 100,
                           /* weights = */ 0);
    print_matrix(&coords);
    print_matrix_list(&routing);
    igraph_destroy(&g);
    printf("\n");

    printf("Path graph with two shortcuts, undirected\n");
    igraph_create(&g, &edgelist, 0, IGRAPH_UNDIRECTED);
    igraph_vector_int_init(&edgelist2, 0);
    igraph_get_edgelist(&g, &edgelist2, 0);
    print_vector_int(&edgelist2);
    igraph_vector_int_destroy(&edgelist2);
    igraph_layout_sugiyama(&g, &coords, &routing, &layers,
                           /* hgap = */ 1,
                           /* vgap = */ 1,
                           /* maxiter = */ 100,
                           /* weights = */ 0);
    print_matrix(&coords);
    print_matrix_list(&routing);
    igraph_destroy(&g);
    printf("\n");

    printf("Path graph with two shortcuts, undirected, opposite layer order\n");
    igraph_vector_int_reverse(&layers);
    igraph_create(&g, &edgelist, 0, IGRAPH_UNDIRECTED);
    igraph_layout_sugiyama(&g, &coords, &routing, &layers,
                           /* hgap = */ 1,
                           /* vgap = */ 1,
                           /* maxiter = */ 100,
                           /* weights = */ 0);
    print_matrix(&coords);
    print_matrix_list(&routing);
    igraph_destroy(&g);
    printf("\n");

    igraph_vector_int_destroy(&layers);
    igraph_vector_int_destroy(&edgelist);

    igraph_matrix_list_destroy(&routing);
    igraph_matrix_destroy(&coords);

    VERIFY_FINALLY_STACK();

    return 0;
}
