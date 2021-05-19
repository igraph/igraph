/* -*- mode: C -*-  */
/* vim:set sw=4 ts=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
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

/* This is a test for bug #1002140, reported by Luiz Fernando
   Bittencourt: https://bugs.launchpad.net/igraph/+bug/1002140 */

#include <igraph.h>

#include "test_utilities.inc"

int main() {

    int k;
    for (k = 0; k < 20; k++) {
        igraph_t g;
        igraph_matrix_t merges;
        igraph_vector_t membership;
        igraph_arpack_options_t options;
        double modularity;
        igraph_vector_t history;
        FILE *DLFile = fopen("input.dl", "r");

        igraph_read_graph_dl(&g, DLFile, /*directed=*/ 0);
        fclose(DLFile);

        igraph_matrix_init(&merges, 0, 0);
        igraph_vector_init(&membership, 0);
        igraph_vector_init(&history, 0);
        igraph_arpack_options_init(&options);

        igraph_community_leading_eigenvector(&g, /*weights=*/ 0, &merges,
                                             &membership, igraph_vcount(&g),
                                             &options, &modularity,
                                             /*start=*/ 0, /*eigenvalues=*/ 0,
                                             /*eigenvectors=*/ 0, &history,
                                             /*callback=*/ 0,
                                             /*callback_extra=*/ 0);

        igraph_vector_destroy(&history);
        igraph_vector_destroy(&membership);
        igraph_matrix_destroy(&merges);
        igraph_destroy(&g);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
