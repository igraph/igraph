/*
   igraph library.
   Copyright (C) 2006-2024  The igraph development team <igraph@igraph.org>

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
    igraph_real_t data[4][4] = { {   0, 1.2, 2.3,   0 },
                                 { 2.0,   0,   0, 1.0 },
                                 {   0,   0, 1.5,   0 },
                                 {   0, 1.0,   0,   0 } };

    /* C arrays use row-major storage, while igraph's matrix uses column-major.
     * The matrix 'mat' will be the transpose of 'data'. */
    const igraph_matrix_t mat =
        igraph_matrix_view(*data, sizeof(data[0]) / sizeof(data[0][0]),
                                  sizeof(data) / sizeof(data[0]));
    igraph_vector_t weights;
    igraph_vector_int_t edges;
    igraph_int_t n;

    /* Initialize the library. */
    igraph_setup();

    /* Initialize vector into which weights will be written. */
    igraph_vector_init(&weights, 0);

    igraph_weighted_adjacency(&graph, &mat, IGRAPH_ADJ_DIRECTED, &weights, IGRAPH_LOOPS_ONCE);

    /* When igraph_weighted_adjacency() returns, 'weights' will typically have
     * more capacity allocated than what it uses. We may optionally free any
     * unused capacity to save memory, although in most applications this
     * is not necessary. */
    igraph_vector_resize_min(&weights);

    /* Get the edge list of the graph and output it, along with the weights. */

    igraph_vector_int_init(&edges, 0);
    igraph_get_edgelist(&graph, &edges, 0);
    n = igraph_ecount(&graph);

    for (igraph_int_t i = 0; i < n; i++) {
        printf("%" IGRAPH_PRId " --> %" IGRAPH_PRId ": %g\n",
               VECTOR(edges)[2*i], VECTOR(edges)[2*i + 1], VECTOR(weights)[i]);
    }

    /* Free all allocated storage. */
    igraph_vector_int_destroy(&edges);
    igraph_destroy(&graph);
    igraph_vector_destroy(&weights);

    return 0;
}
