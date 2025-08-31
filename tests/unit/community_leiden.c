/*
   igraph library.
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

#include "test_utilities.h"

#define TOL (1e-15)

void run_leiden_CPM(const igraph_t *graph, const igraph_vector_t *edge_weights, const igraph_real_t resolution) {

    igraph_vector_int_t membership;
    igraph_int_t nb_clusters = igraph_vcount(graph);
    igraph_real_t quality, quality2;

    /* Initialize with singleton partition. */
    igraph_vector_int_init(&membership, igraph_vcount(graph));

    /* Use same seed as for the simplified interface below, to ensure the same result. */
    igraph_rng_seed(igraph_rng_default(), 123);
    igraph_community_leiden(graph,
                            edge_weights, NULL, NULL,
                            resolution, 0.01, false, 2, &membership,
                            &nb_clusters,
                            &quality);

    /* Handle negative zeros. */
    if (fabs(quality) < TOL) quality = 0.0;

    printf("Leiden found %" IGRAPH_PRId " clusters using CPM (resolution parameter=%.2f), quality is %.5f.\n", nb_clusters, resolution, quality);

    printf("Membership: ");
    igraph_vector_int_print(&membership);
    printf("\n");

    quality2 = quality;

    /* Use same seed as for the generic interface above, to ensure the same result. */
    igraph_rng_seed(igraph_rng_default(), 123);
    igraph_community_leiden_simple(graph,
                                   edge_weights,
                                   IGRAPH_LEIDEN_OBJECTIVE_CPM,
                                   resolution, 0.01, false, 2, &membership,
                                   NULL,
                                   &quality);
    if (fabs(quality) < TOL) quality = 0.0;
    IGRAPH_ASSERT((isnan(quality) && isnan(quality2)) ||
                  igraph_almost_equals(quality, quality2, TOL));

    igraph_vector_int_destroy(&membership);
}

void run_leiden_modularity(igraph_t *graph, igraph_vector_t *edge_weights) {

    const igraph_bool_t directed = igraph_is_directed(graph);
    igraph_vector_int_t membership;
    igraph_vector_t out_strength, in_strength;
    igraph_int_t nb_clusters = igraph_vcount(graph);
    igraph_real_t quality, quality2;
    const igraph_real_t directed_multiplier = directed ? 1.0 : 2.0;
    igraph_real_t m;

    igraph_vector_init(&out_strength, igraph_vcount(graph));
    igraph_strength(graph, &out_strength, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS, edge_weights);

    if (directed) {
        igraph_vector_init(&in_strength, igraph_vcount(graph));
        igraph_strength(graph, &in_strength, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS, edge_weights);
    }

    m = edge_weights ? igraph_vector_sum(edge_weights) : igraph_ecount(graph);

    /* Initialize with singleton partition. */
    igraph_vector_int_init(&membership, igraph_vcount(graph));

    /* Use same seed as for the simplified interface below, to ensure the same result. */
    igraph_rng_seed(igraph_rng_default(), 123);
    igraph_community_leiden(graph,
                            edge_weights, &out_strength, directed ? &in_strength : NULL,
                            1.0 / (directed_multiplier * m), 0.01, false, 2, &membership, &nb_clusters,
                            &quality);

    igraph_modularity(graph, &membership, edge_weights, 1.0, IGRAPH_DIRECTED, &quality2);
    if (isnan(quality)) {
        printf("Leiden found %" IGRAPH_PRId " clusters using modularity, %s, quality is nan.\n", nb_clusters, directed ? "directed" : "undirected");
        IGRAPH_ASSERT(isnan(quality2));
    } else {
        /* It is necessary to not only use igraph_almost_equals(), but also check
         * if the values are both very close to zero due to roundoff errors. */

        if (fabs(quality) < TOL) quality = 0.0;
        if (fabs(quality2) < TOL) quality2 = 0.0;

        printf("Leiden found %" IGRAPH_PRId " clusters using modularity, %s, quality is %.5f.\n", nb_clusters, directed ? "directed" : "undirected", quality);
        IGRAPH_ASSERT(igraph_almost_equals(quality, quality2, TOL));
    }

    printf("Membership: ");
    igraph_vector_int_print(&membership);
    printf("\n");

    /* Use same seed as for the generic interface above, to ensure the same result. */
    igraph_rng_seed(igraph_rng_default(), 123);
    igraph_community_leiden_simple(graph,
                                   edge_weights,
                                   IGRAPH_LEIDEN_OBJECTIVE_MODULARITY,
                                   1.0, 0.01, false, 2, &membership, NULL, &quality);
    if (fabs(quality) < TOL) quality = 0.0;
    IGRAPH_ASSERT((isnan(quality) && isnan(quality2)) ||
                  igraph_almost_equals(quality, quality2, TOL));

    igraph_vector_int_destroy(&membership);
    if (directed) igraph_vector_destroy(&in_strength);
    igraph_vector_destroy(&out_strength);
}

int main(void) {
    igraph_t graph;
    igraph_vector_t weights;

    igraph_vector_init(&weights, 0);

    /* Set default seed to get reproducible results */
    igraph_rng_seed(igraph_rng_default(), 0);

    /* Simple unweighted graph */
    igraph_small(&graph, 10, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4,
                 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 6, 9, 7, 8, 7, 9, 8, 9,
                 0, 5, -1);
    run_leiden_modularity(&graph, NULL);

    /* Same simple graph, with uniform edge weights */
    igraph_vector_resize(&weights, igraph_ecount(&graph));
    igraph_vector_fill(&weights, 2);
    run_leiden_modularity(&graph, &weights);

    /* Same simple graph, but directed with reciprocal edges */
    igraph_to_directed(&graph, IGRAPH_TO_DIRECTED_MUTUAL);
    run_leiden_modularity(&graph, NULL);

    igraph_destroy(&graph);

    /* Tiny directed graph; optimal community structure is different if
     * ignoring edge directions. */
    igraph_small(&graph, 4, IGRAPH_DIRECTED, 0, 2, 0, 3, 1, 2, 3, 1, 3, 2, -1);
    run_leiden_modularity(&graph, NULL);
    igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_EACH, NULL);
    run_leiden_modularity(&graph, NULL);
    igraph_destroy(&graph);

    /* Larger directed graph; optimal community structure is different if
     * ignoring edge directions. */
    igraph_small(
        &graph, 10, IGRAPH_DIRECTED,
        0, 3, 0, 4, 1, 0, 1, 4, 2, 1, 3, 0, 3, 2, 4, 0, 4, 3, 4, 7, 5, 0, 5,
        1, 5, 3, 5, 6, 5, 8, 5, 9, 7, 0, 8, 2, 8, 3, 9, 1, 9, 3, 9, 8,
        -1);
    run_leiden_modularity(&graph, NULL);
    igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_EACH, NULL);
    run_leiden_modularity(&graph, NULL);
    igraph_destroy(&graph);

    /* Simple nonuniform weighted graph, with and without weights */
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 3, 2, 4, 2, 5, 3, 4, 3, 5, 4, 5, -1);
    igraph_vector_resize(&weights, 8);
    igraph_vector_fill(&weights, 1);
    VECTOR(weights)[0] = 10;
    VECTOR(weights)[1] = 10;
    run_leiden_modularity(&graph, NULL);
    run_leiden_modularity(&graph, &weights);
    igraph_destroy(&graph);

    /* Zachary Karate club */
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 0,  1,  0,  2,  0,  3,  0,  4,  0,  5,
                 0,  6,  0,  7,  0,  8,  0, 10,  0, 11,
                 0, 12,  0, 13,  0, 17,  0, 19,  0, 21,
                 0, 31,  1,  2,  1,  3,  1,  7,  1, 13,
                 1, 17,  1, 19,  1, 21,  1, 30,  2,  3,
                 2,  7,  2,  8,  2,  9,  2, 13,  2, 27,
                 2, 28,  2, 32,  3,  7,  3, 12,  3, 13,
                 4,  6,  4, 10,  5,  6,  5, 10,  5, 16,
                 6, 16,  8, 30,  8, 32,  8, 33,  9, 33,
                 13, 33, 14, 32, 14, 33, 15, 32, 15, 33,
                 18, 32, 18, 33, 19, 33, 20, 32, 20, 33,
                 22, 32, 22, 33, 23, 25, 23, 27, 23, 29,
                 23, 32, 23, 33, 24, 25, 24, 27, 24, 31,
                 25, 31, 26, 29, 26, 33, 27, 33, 28, 31,
                 28, 33, 29, 32, 29, 33, 30, 32, 30, 33,
                 31, 32, 31, 33, 32, 33,
                 -1);
    run_leiden_modularity(&graph, NULL);
    run_leiden_CPM(&graph, NULL, 0.06);
    igraph_destroy(&graph);

    /* Simple disconnected graph with isolates */
    igraph_small(&graph, 9, IGRAPH_UNDIRECTED,
                 0,  1,  0,  2,  0,  3,  1,  2,  1,  3,  2,  3,
                 4,  5,  4,  6,  4,  7,  5,  6,  5,  7,  6,  7,
                 -1);
    run_leiden_modularity(&graph, NULL);
    igraph_destroy(&graph);

    /* Disjoint union of two rings */
    igraph_small(&graph, 20, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 0, 9,
                 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 10, 19, -1);
    run_leiden_modularity(&graph, NULL);
    run_leiden_CPM(&graph, NULL, 0.05);
    igraph_destroy(&graph);

    /* Completely empty graph */
    igraph_small(&graph, 10, IGRAPH_UNDIRECTED, -1);
    run_leiden_modularity(&graph, NULL);
    igraph_destroy(&graph);


    /* Set default seed to get reproducible results */
    igraph_rng_seed(igraph_rng_default(), 0);

    /* Ring graph without loop edges */
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,3, 3,4, 4,5, 5,0, -1);
    run_leiden_CPM(&graph, NULL, 0.4);
    igraph_destroy(&graph);

    /* Set default seed to get reproducible results */
    igraph_rng_seed(igraph_rng_default(), 0);

    /* Ring graph with loop edges */
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,3, 3,4, 4,5, 5,0,
                 0,0, 1,1, 2,2, 3,3, 4,4, 5,5,
                 -1);
    run_leiden_CPM(&graph, NULL, 0.4);
    igraph_destroy(&graph);

    /* Regression test -- graph with two vertices and two edges */
    igraph_small(&graph, 2, IGRAPH_UNDIRECTED, 0, 0, 1, 1, -1);
    run_leiden_modularity(&graph, NULL);
    igraph_destroy(&graph);

    /* The next two tests need an empty weight vector. */
    igraph_vector_clear(&weights);

    /* Null graph */
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    run_leiden_modularity(&graph, &weights);
    igraph_destroy(&graph);

    /* Edgeless graph */
    igraph_empty(&graph, 5, IGRAPH_UNDIRECTED);
    run_leiden_modularity(&graph, &weights);
    igraph_destroy(&graph);

    /* Check that the input is validated properly. */

    /* Small test graph. */
    igraph_small(&graph, 4, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,0, 0,3, 3,3,
                 -1);
    igraph_vector_range(&weights, 1, igraph_ecount(&graph) + 1);

    /* Omitting membership should raise no error. */
    igraph_community_leiden_simple(&graph, &weights, IGRAPH_LEIDEN_OBJECTIVE_MODULARITY,
                                   1.0, 0.01, false, 1, NULL, NULL, NULL);

    /* Negative weight. */
    VECTOR(weights)[0] = -1;
    CHECK_ERROR(igraph_community_leiden_simple(&graph, &weights, IGRAPH_LEIDEN_OBJECTIVE_MODULARITY,
                                               1.0, 0.01, false, 1, NULL, NULL, NULL), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_community_leiden_simple(&graph, &weights, IGRAPH_LEIDEN_OBJECTIVE_ER,
                                               1.0, 0.01, false, 1, NULL, NULL, NULL), IGRAPH_EINVAL);

    /* NaN weight. */
    VECTOR(weights)[0] = IGRAPH_NAN;
    CHECK_ERROR(igraph_community_leiden_simple(&graph, &weights, IGRAPH_LEIDEN_OBJECTIVE_CPM,
                                               1.0, 0.01, false, 1, NULL, NULL, NULL), IGRAPH_EINVAL);

    /* Invalid weight vector length. */
    igraph_vector_range(&weights, 1, igraph_ecount(&graph) + 2);
    CHECK_ERROR(igraph_community_leiden_simple(&graph, &weights, IGRAPH_LEIDEN_OBJECTIVE_MODULARITY,
                                               1.0, 0.01, false, 1, NULL, NULL, NULL), IGRAPH_EINVAL);

    igraph_destroy(&graph);

    igraph_vector_destroy(&weights);

    VERIFY_FINALLY_STACK();

    return 0;
}
