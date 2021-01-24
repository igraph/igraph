/* -*- mode: C -*-  */
/*
   IGraph library.
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

#include "test_utilities.inc"

int main() {
    igraph_t g, extd_g;
    igraph_matrix_t coords;
    igraph_vector_t edgelist, extd_edgelist, extd_to_orig_eids;
    igraph_vector_t layers;

    igraph_matrix_init(&coords, 0, 0);
    igraph_vector_init(&extd_to_orig_eids, 0);

    /* Layout on simple graph with predefined layers */
    igraph_vector_init_int_end(&layers, -1, 0, 1, 1, 2, 3, 3, 4, 4, 5, -1);
    igraph_vector_init_int_end(&edgelist, -1,
                               0, 1, 0, 2, 0, 3, 1, 2, 2, 2, 1, 4, 2, 5, 4, 6, 5, 7, 6, 8, 7, 8,
                               3, 8, 8, 1, 8, 2, -1);
    igraph_create(&g, &edgelist, 0, 1);

    igraph_layout_sugiyama(&g, &coords, 0, 0, &layers,
                           /* hgap = */ 1,
                           /* vgap = */ 1,
                           /* maxiter = */ 100,
                           /* weights = */ 0);
    igraph_matrix_print(&coords);
    printf("===\n");

    /* Same, but this time also return the extended graph */
    igraph_layout_sugiyama(&g, &coords, &extd_g, &extd_to_orig_eids, &layers,
                           /* hgap = */ 1,
                           /* vgap = */ 1,
                           /* maxiter = */ 100,
                           /* weights = */ 0);
    igraph_matrix_print(&coords);
    printf("===\n");
    igraph_vector_init(&extd_edgelist, 0);
    igraph_get_edgelist(&extd_g, &extd_edgelist, 0);
    igraph_vector_print(&extd_edgelist);
    igraph_vector_destroy(&extd_edgelist);
    igraph_destroy(&extd_g);
    printf("===\n");
    igraph_vector_print(&extd_to_orig_eids);
    printf("===\n");

    igraph_vector_destroy(&layers);

    /* Same, but with automatic layering */
    igraph_layout_sugiyama(&g, &coords, 0, 0, 0,
                           /* hgap = */ 1,
                           /* vgap = */ 1,
                           /* maxiter = */ 100,
                           /* weights = */ 0);
    igraph_matrix_print(&coords);
    printf("===\n");

    /* Layering with gaps in it */
    igraph_vector_init_int_end(&layers, -1, 0, 2, 2, 4, 6, 6, 12, 12, 15, -1);
    igraph_layout_sugiyama(&g, &coords, 0, 0, &layers,
                           /* hgap = */ 1,
                           /* vgap = */ 1,
                           /* maxiter = */ 100,
                           /* weights = */ 0);
    igraph_matrix_print(&coords);
    igraph_vector_destroy(&layers);
    printf("===\n");

    igraph_vector_destroy(&edgelist);
    igraph_matrix_destroy(&coords);
    igraph_vector_destroy(&extd_to_orig_eids);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
