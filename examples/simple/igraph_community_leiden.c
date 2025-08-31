/*
   igraph library.
   Copyright (C) 2007-2024  The igraph development team <igraph@igraph.org>

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
#include <stdio.h>

int main(void) {
    igraph_t graph;
    igraph_vector_int_t membership;
    igraph_vector_t vertex_weights;
    igraph_vector_t vertex_out_weights, vertex_in_weights;
    igraph_int_t nb_clusters;
    igraph_real_t quality;

    /* Initialize the library. */
    igraph_setup();

    /* Set default seed to get reproducible results */
    igraph_rng_seed(igraph_rng_default(), 0);

    /* UNDIRECTED EXAMPLE */

    /* Simple unweighted undirected graph */
    igraph_small(&graph, 10, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4,
                 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 6, 9, 7, 8, 7, 9, 8, 9,
                 0, 5, -1);

    /* Perform Leiden algorithm using CPM for 1 iteration */
    igraph_vector_int_init(&membership, igraph_vcount(&graph));
    igraph_community_leiden(&graph, NULL, NULL, NULL,
                            /* resolution */ 0.05,
                            /* beta */ 0.01,
                            /* start */ false,
                            /* iterations */ 1,
                            &membership, &nb_clusters, &quality);

    printf("Leiden found %" IGRAPH_PRId " clusters using CPM (resolution parameter 0.05), quality is %.4f.\n", nb_clusters, quality);
    printf("Membership: ");
    igraph_vector_int_print(&membership);
    printf("\n");

    /* Start from existing membership for 10 iterations to improve it further */
    igraph_community_leiden(&graph, NULL, NULL, NULL,
                            /* resolution */ 0.05,
                            /* beta */ 0.01,
                            /* start */ true,
                            /* iterations */ 10,
                            &membership, &nb_clusters, &quality);

    printf("Iterated Leiden, using CPM (resolution parameter 0.05), quality is %.4f.\n", quality);
    printf("Membership: ");
    igraph_vector_int_print(&membership);
    printf("\n");

    /* Use degrees as vertex weights for optimizing modularity */
    igraph_vector_init(&vertex_weights, igraph_vcount(&graph));
    igraph_strength(&graph, &vertex_weights, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS, NULL);

    /* Perform Leiden algorithm using modularity until stable iteration */
    igraph_community_leiden(&graph, NULL, &vertex_weights, NULL,
                            /* resolution */ 1.0 / (2 * igraph_ecount(&graph)),
                            /* beta */ 0.01,
                            /* start */ false,
                            /* iterations */ -1,
                            &membership, &nb_clusters, &quality);

    printf("Leiden found %" IGRAPH_PRId " clusters using modularity, quality is %.4f.\n", nb_clusters, quality);
    printf("Membership: ");
    igraph_vector_int_print(&membership);
    printf("\n");

    igraph_vector_destroy(&vertex_weights);
    igraph_vector_int_destroy(&membership);
    igraph_destroy(&graph);

    /* DIRECTED EXAMPLE */

    /* Simple unweighted directed graph */
    igraph_small(&graph, 6, IGRAPH_DIRECTED,
                 0, 1, 0, 3, 0, 5, 1, 3, 1, 4, 2, 3, 4, 1, 5, 2, 5, 4,
                 -1);

    igraph_vector_int_init(&membership, igraph_vcount(&graph));
    igraph_vector_init(&vertex_out_weights, igraph_vcount(&graph));
    igraph_vector_init(&vertex_in_weights, igraph_vcount(&graph));
    igraph_strength(&graph, &vertex_out_weights, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS, NULL);
    igraph_strength(&graph, &vertex_in_weights, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS, NULL);

    /* Perform Leiden algorithm using modularity for two iterations,
     * which are usually sufficient to achieve stability.
     * Note that in the directed case, we divide by the edge count m, not 2*m. */
    igraph_community_leiden(&graph, NULL, &vertex_out_weights, &vertex_in_weights,
                            /* resolution */ 1.0 / igraph_ecount(&graph),
                            /* beta */ 0.01,
                            /* start */ false,
                            /* iterations */ 2,
                            &membership, &nb_clusters, &quality);

    printf("Leiden found %" IGRAPH_PRId " clusters using modularity, quality is %.4f.\n", nb_clusters, quality);
    printf("Membership: ");
    igraph_vector_int_print(&membership);
    printf("\n");

    igraph_vector_destroy(&vertex_in_weights);
    igraph_vector_destroy(&vertex_out_weights);
    igraph_vector_int_destroy(&membership);
    igraph_destroy(&graph);

    return 0;
}
