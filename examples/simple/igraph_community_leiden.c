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

int main(void) {
    igraph_t graph;
    igraph_vector_int_t membership;
    igraph_vector_int_t degree;
    igraph_vector_t weights;
    igraph_integer_t nb_clusters, i;
    igraph_real_t quality;

    /* Set default seed to get reproducible results */
    igraph_rng_seed(igraph_rng_default(), 0);

    /* Simple unweighted graph */
    igraph_small(&graph, 10, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4,
                 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 6, 9, 7, 8, 7, 9, 8, 9,
                 0, 5, -1);

    /* Perform Leiden algorithm using CPM for 1 iteration */
    igraph_vector_int_init(&membership, igraph_vcount(&graph));
    igraph_community_leiden(&graph, NULL, NULL, 0.05, 0.01, 0, 1, &membership, &nb_clusters, &quality);

    printf("Leiden found %" IGRAPH_PRId " clusters using CPM (resolution parameter 0.05), quality is %.4f.\n", nb_clusters, quality);
    printf("Membership: ");
    igraph_vector_int_print(&membership);
    printf("\n");

    /* Start from existing membership for 10 iterations to improve it further */
    igraph_community_leiden(&graph, NULL, NULL, 0.05, 0.01, 1, 10, &membership, &nb_clusters, &quality);

    printf("Iterated Leiden, using CPM (resolution parameter 0.05), quality is %.4f.\n", quality);
    printf("Membership: ");
    igraph_vector_int_print(&membership);
    printf("\n");

    /* Initialize degree vector to use for optimizing modularity */
    igraph_vector_int_init(&degree, igraph_vcount(&graph));
    igraph_degree(&graph, &degree, igraph_vss_all(), IGRAPH_ALL, 1);
    igraph_vector_init(&weights, igraph_vector_int_size(&degree));
    for (i = 0; i < igraph_vector_int_size(&degree); i++) {
        VECTOR(weights)[i] = VECTOR(degree)[i];
    }

    /* Perform Leiden algorithm using modularity until stable iteration */
    igraph_community_leiden(&graph, NULL, &weights, 1.0 / (2 * igraph_ecount(&graph)), 0.01, 0, -1, &membership, &nb_clusters, &quality);

    printf("Leiden found %" IGRAPH_PRId " clusters using modularity, quality is %.4f.\n", nb_clusters, quality);
    printf("Membership: ");
    igraph_vector_int_print(&membership);
    printf("\n");

    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&degree);
    igraph_vector_int_destroy(&membership);
    igraph_destroy(&graph);

    return 0;
}
