/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
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

/*

    R
    library(igraph)
    g <- graph.tree(10, 3, mode="out")
    A <- get.adjacency(g)
    svd(A + .5 * degree(g) * diag(vcount(g)))

*/

int main() {

    igraph_t graph;
    igraph_matrix_t U, V;
    igraph_arpack_options_t options;
    igraph_vector_t cvec;

    igraph_tree(&graph, /*n=*/ 14, /*children=*/ 4, IGRAPH_TREE_OUT);

    igraph_matrix_init(&U, 0, 0);
    igraph_matrix_init(&V, 0, 0);
    igraph_arpack_options_init(&options);

    igraph_vector_init(&cvec, 0);
    igraph_degree(&graph, &cvec, igraph_vss_all(), IGRAPH_ALL,
                  IGRAPH_LOOPS);
    igraph_vector_scale(&cvec, .5);

    igraph_adjacency_spectral_embedding(&graph, 4, /*weights=*/ 0,
                                        IGRAPH_EIGEN_LA,
                                        /*scaled=*/ 0, &U, &V, /*D=*/ 0,
                                        &cvec, &options);

    /* eigenvectors are in the columns of U and V; make sure that the
     * first row contains positive values */
    print_matrix_first_row_positive(&U, "%8.4f");
    printf("--\n");
    print_matrix_first_row_positive(&V, "%8.4f");

    igraph_vector_destroy(&cvec);
    igraph_matrix_destroy(&V);
    igraph_matrix_destroy(&U);

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
