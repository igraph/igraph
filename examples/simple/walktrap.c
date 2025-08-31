/*
   igraph library.
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

int main(void) {
    igraph_t g;
    igraph_matrix_int_t merges;
    igraph_vector_t modularity;
    igraph_int_t no_of_nodes;
    igraph_int_t i;

    /* Initialize the library. */
    igraph_setup();

    igraph_small(&g, 5, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4,
                 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 6, 9, 7, 8, 7, 9, 8, 9, 0, 5, 4, 9, -1);
    igraph_vector_init(&modularity, 0);
    igraph_matrix_int_init(&merges, 0, 0);

    igraph_community_walktrap(&g,
                              NULL /* no weights */,
                              4 /* steps */,
                              &merges, &modularity,
                              /* membership=*/ NULL);

    no_of_nodes = igraph_vcount(&g);
    printf("Merges:\n");
    for (i = 0; i < igraph_matrix_int_nrow(&merges); i++) {
        printf("%2.1" IGRAPH_PRId " + %2." IGRAPH_PRId " -> %2." IGRAPH_PRId " (modularity %4.2f)\n",
               MATRIX(merges, i, 0), MATRIX(merges, i, 1),
               no_of_nodes + i, VECTOR(modularity)[i]);
    }

    igraph_destroy(&g);

    igraph_matrix_int_destroy(&merges);
    igraph_vector_destroy(&modularity);

    return 0;
}
