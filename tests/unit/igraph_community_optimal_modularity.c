/* IGraph library.
   Copyright (C) 2010-2024  The igraph development team <igraph@igraph.org>

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

void prepare_weights_vector(igraph_vector_t *weights, const igraph_t *graph) {
    const igraph_integer_t m =igraph_ecount(graph);
    igraph_edge_betweenness(graph, weights, true, NULL);
    for (igraph_integer_t i=0; i < m; i++) {
        VECTOR(*weights)[i] = 1 / VECTOR(*weights)[i];
    }
}

/* Only for undirected! */
void verify_with_leiden(const igraph_t *graph, const igraph_vector_t *weights,
                        igraph_real_t resolution, igraph_real_t modularity) {

    igraph_vector_int_t leiden_membership;
    igraph_vector_t vertex_weights;
    igraph_real_t Q, maxQ = 0.0;
    igraph_real_t m = weights ? igraph_vector_sum(weights) : igraph_ecount(graph);

    igraph_vector_init(&vertex_weights, 0);
    igraph_vector_int_init(&leiden_membership, 0);

    igraph_strength(graph, &vertex_weights, igraph_vss_all(), IGRAPH_ALL, /*loops*/ true, weights);

    for (int i=0; i < 10; i++) {
        igraph_community_leiden(graph, weights, &vertex_weights,
                                resolution / (2*m),
                                0.01, false, 2, &leiden_membership, NULL, &Q);
        if (Q > maxQ) {
            maxQ = Q;
        }
    }

    if (igraph_cmp_epsilon(modularity, maxQ, 1e-15) < 0) {
        printf("Partitioning doesn't achieve maximum modulairty. "
               "With Leiden: Q=%g; with optimal modularity: Q=%g\n",
               maxQ, modularity);
        IGRAPH_ASSERT(igraph_cmp_epsilon(modularity, maxQ, 1e-15) >= 0);
    }

    igraph_vector_int_destroy(&leiden_membership);
    igraph_vector_destroy(&vertex_weights);
}

int main(void) {
    igraph_t graph;
    igraph_vector_t weights;
    igraph_vector_int_t membership;
    igraph_real_t modularity;
    igraph_error_handler_t *handler;
    igraph_error_t errcode;

    igraph_vector_init(&weights, 0);
    igraph_vector_int_init(&membership, 0);

    igraph_famous(&graph, "zachary");

    /* Zachary karate club, unweighted */
    handler = igraph_set_error_handler(&igraph_error_handler_printignore);
    errcode = igraph_community_optimal_modularity(&graph, NULL, 1, &modularity, &membership);
    if (errcode == IGRAPH_UNIMPLEMENTED) {
        igraph_vector_int_destroy(&membership);
        igraph_vector_destroy(&weights);
        return 77; /* skip test */
    } else {
        IGRAPH_ASSERT(errcode == IGRAPH_SUCCESS);
    }
    igraph_set_error_handler(handler);

    verify_with_leiden(&graph, NULL, 1, modularity);
    IGRAPH_ASSERT(igraph_almost_equals(modularity, 0.4197896120973046, 1e-15));

    /* Zachary karate club, weighted */
    prepare_weights_vector(&weights, &graph);
    igraph_community_optimal_modularity(&graph, &weights, 1, &modularity, &membership);
    verify_with_leiden(&graph, &weights, 1, modularity);
    IGRAPH_ASSERT(igraph_almost_equals(modularity, 0.6138050900169181, 1e-15));

    /* Zachary karate club, unweighted, resolution=2 */
    igraph_community_optimal_modularity(&graph, NULL, 2, &modularity, &membership);
    verify_with_leiden(&graph, NULL, 2, modularity);
    IGRAPH_ASSERT(igraph_almost_equals(modularity, 0.16452991452991447, 1e-15));

    /* Zachary karate club, weighted, resolution=2 */
    igraph_community_optimal_modularity(&graph, &weights, 2, &modularity, &membership);
    verify_with_leiden(&graph, &weights, 2, modularity);
    IGRAPH_ASSERT(igraph_almost_equals(modularity, 0.4006998679817418, 1e-15));

    /* Zachary karate club, unweighted, resolution=0.5 */
    igraph_community_optimal_modularity(&graph, NULL, 0.5, &modularity, &membership);
    verify_with_leiden(&graph, NULL, 0.5, modularity);
    IGRAPH_ASSERT(igraph_almost_equals(modularity, 0.6217948717948718, 1e-15));

    /* Zachary karate club, weighted, resolution=0.5 */
    igraph_community_optimal_modularity(&graph, &weights, 0.5, &modularity, &membership);
    verify_with_leiden(&graph, &weights, 0.5, modularity);
    IGRAPH_ASSERT(igraph_almost_equals(modularity, 0.752255118181729, 1e-15));

    igraph_destroy(&graph);

    /* simple graph with self-loops and multi-edges, unweighted */
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED,
                 0, 1, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 0, 0, 2, 2, 2, 2, -1);

    igraph_community_optimal_modularity(&graph, NULL, 1, &modularity, &membership);
    verify_with_leiden(&graph, NULL, 1, modularity);
    IGRAPH_ASSERT(igraph_almost_equals(modularity, 0.36, 1e-15));

    /* simple graph with self-loops and multi-edges, weighted */
    igraph_vector_resize(&weights, igraph_ecount(&graph));
    igraph_vector_fill(&weights, 1.0);
    VECTOR(weights)[2] = 5.0;
    VECTOR(weights)[4] = 0.1;

    igraph_community_optimal_modularity(&graph, &weights, 1, &modularity, &membership);
    verify_with_leiden(&graph, &weights, 1, modularity);
    IGRAPH_ASSERT(igraph_almost_equals(modularity, 0.266855078375386, 1e-15));

    igraph_destroy(&graph);

    /* Directed graph, chosen so that:
     *  - the directed/undirected solutions differ
     *  - there is a single maximum on both cases */
    igraph_small(&graph, 5, IGRAPH_DIRECTED,
                 1, 0, 1, 2, 1, 4, 2, 0, 2, 3, 3, 1, 3, 4, 4, 0, -1);
    igraph_community_optimal_modularity(&graph, NULL, 1, &modularity, &membership);
    IGRAPH_ASSERT(igraph_almost_equals(modularity, 0.09375, 1e-15));

    igraph_destroy(&graph);

    igraph_vector_int_destroy(&membership);
    igraph_vector_destroy(&weights);

    return 0;
}
