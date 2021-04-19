/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge MA, 02139 USA

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

int main() {
    igraph_t g;
    igraph_matrix_t merges;
    igraph_vector_t modularity;
    long int no_of_nodes;
    long int i;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_small(&g, 5, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4,
                 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 6, 9, 7, 8, 7, 9, 8, 9, 0, 5, -1);
    igraph_vector_init(&modularity, 0);
    igraph_matrix_init(&merges, 0, 0);

    igraph_community_walktrap(&g, 0 /* no weights */,
                              4 /* steps */,
                              &merges, &modularity,
                              /* membership=*/ 0);

    no_of_nodes = igraph_vcount(&g);
    printf("Merges:\n");
    for (i = 0; i < igraph_matrix_nrow(&merges); i++) {
        printf("%2.1li + %2.li -> %2.li (modularity %4.2f)\n",
               (long int)MATRIX(merges, i, 0),
               (long int)MATRIX(merges, i, 1),
               no_of_nodes + i,
               VECTOR(modularity)[i]);
    }

    igraph_destroy(&g);

    /* isolated vertices */
    igraph_small(&g, 5, IGRAPH_UNDIRECTED, -1);
    if (igraph_community_walktrap(&g, 0 /* no weights */, 4 /* steps */, &merges,
                                  &modularity, /* membership = */ 0)) {
        return 1;
    }
    if (igraph_vector_min(&modularity) != 0 || igraph_vector_max(&modularity) != 0) {
        return 2;
    }
    igraph_destroy(&g);

    igraph_matrix_destroy(&merges);
    igraph_vector_destroy(&modularity);
    return 0;
}
