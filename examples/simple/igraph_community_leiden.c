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

int main() {
    igraph_t graph;
    igraph_vector_t membership, degree;
    igraph_integer_t nb_clusters;
    igraph_real_t quality;

    /* Set default seed to get reproducible results */
    igraph_rng_seed(igraph_rng_default(), 0);

    /* Simple unweighted graph */
    igraph_small(&graph, 10, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4,
                 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 6, 9, 7, 8, 7, 9, 8, 9,
                 0, 5, -1);

    /* Perform Leiden algorithm using CPM */
    igraph_vector_init(&membership, igraph_vcount(&graph));
    igraph_community_leiden(&graph, NULL, NULL, 0.05, 0.01, 0, &membership, &nb_clusters, &quality);

    printf("Leiden found %" IGRAPH_PRId " clusters using CPM (resolution parameter 0.05), quality is %.4f.\n", nb_clusters, quality);
    printf("Membership: ");
    igraph_vector_print(&membership);
    printf("\n");

    /* Start from existing membership to improve it further */
    igraph_community_leiden(&graph, NULL, NULL, 0.05, 0.01, 1, &membership, &nb_clusters, &quality);

    printf("Iterated Leiden, using CPM (resolution parameter 0.05), quality is %.4f.\n", quality);
    printf("Membership: ");
    igraph_vector_print(&membership);
    printf("\n");

    /* Initialize degree vector to use for optimizing modularity */
    igraph_vector_init(&degree, igraph_vcount(&graph));
    igraph_degree(&graph, &degree, igraph_vss_all(), IGRAPH_ALL, 1);

    /* Perform Leiden algorithm using modularity */
    igraph_community_leiden(&graph, NULL, &degree, 1.0 / (2 * igraph_ecount(&graph)), 0.01, 0, &membership, &nb_clusters, &quality);

    printf("Leiden found %" IGRAPH_PRId " clusters using modularity, quality is %.4f.\n", nb_clusters, quality);
    printf("Membership: ");
    igraph_vector_print(&membership);
    printf("\n");

    igraph_vector_destroy(&degree);
    igraph_vector_destroy(&membership);
    igraph_destroy(&graph);

    return 0;
}

