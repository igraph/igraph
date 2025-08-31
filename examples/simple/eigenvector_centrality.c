/*
   igraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include <igraph.h>

int main(void) {

    igraph_t graph;
    igraph_vector_t vector, weights;
    igraph_real_t value;

    /* Initialize the library. */
    igraph_setup();

    /* Create a star graph, with vertex 0 at the center, and associated edge weights. */
    igraph_star(&graph, 10, IGRAPH_STAR_UNDIRECTED, 0);
    igraph_vector_init_range(&weights, 1, igraph_ecount(&graph)+1);

    /* Initialize the vector where the result will be stored. */
    igraph_vector_init(&vector, 0);

    /* Compute eigenvector centrality. */
    igraph_eigenvector_centrality(&graph, &vector, &value, IGRAPH_OUT,
                                  &weights, /*options=*/ NULL);

    /* Print results. */
    printf("eigenvalue: %g\n", value);
    printf("eigenvector:\n");
    igraph_vector_print(&vector);

    /* Free allocated data structures. */
    igraph_vector_destroy(&vector);
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    return 0;
}
