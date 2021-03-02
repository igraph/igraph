/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2021  The igraph development team <igraph@igraph.org>

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
#include <math.h>
#include <float.h>

#include "test_utilities.inc"

int is_almost_one(igraph_real_t x) {
    /* 2^5 = 32 is 5 binary digits  of tolerance */
    if (fabs(x - 1) > 32*DBL_EPSILON) {
        printf("Expected value to be 1, but actually got %15g.", x);
        return 0;
    }
    return 1;
}

int main() {

    igraph_t g;
    igraph_vector_t res, reset, weights;
    igraph_arpack_options_t arpack_options;
    igraph_real_t value;
    int err;

    /* The ARPACK method uses a random perturbation to the in-degrees
       to set the starting vector for ARPACK. */
    igraph_rng_seed(igraph_rng_default(), 137);

    igraph_arpack_options_init(&arpack_options);

    /* Test graphs taken from http://www.iprcom.com/papers/pagerank/ */

    printf("\nTest graph 1\n");
    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,0, 3,2, 0,2,
                 -1);

    igraph_vector_init(&res, 0);
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 0, 0.85, 0, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 0, 0.85, 0, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_vector_destroy(&res);

    igraph_destroy(&g);

    printf("\nTest graph 2\n");
    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 0,1, 0,2, 0,3, 1,0, 2,0, 3,0,
                 3,4, 3,5, 3,6, 3,7, 4,0, 5,0,
                 6,0, 7,0,
                 -1);

    igraph_vector_init(&res, 0);
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 0, 0.85, 0, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 0, 0.85, 0, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_vector_destroy(&res);
    igraph_destroy(&g);

    /* Undirected star graph */

    printf("\nUndirected star\n");

    igraph_star(&g, 11, IGRAPH_STAR_UNDIRECTED, 0);
    igraph_vector_init(&res, 0);
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 0, 0.85, 0, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 0, 0.85, 0, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    /* Check twice more for consistency */
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 0, 0.85, 0, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 0, 0.85, 0, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 0, 0.85, 0, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 0, 0.85, 0, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    /* Check personalized PageRank */

    printf("\nPersonalized PageRank\n");

    igraph_personalized_pagerank_vs(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                                    igraph_vss_all(), 0, 0.5,
                                    igraph_vss_1(1), 0, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_personalized_pagerank_vs(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                                    igraph_vss_all(), 0, 0.5,
                                    igraph_vss_1(1), 0, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));


    /* Errors */

    igraph_set_error_handler(igraph_error_handler_ignore);
    igraph_vector_init(&reset, 2);
    err = igraph_personalized_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, 0,
                                       igraph_vss_all(), 0, 0.85, &reset, 0,
                                       &arpack_options);
    IGRAPH_ASSERT(err == IGRAPH_EINVAL);

    err = igraph_personalized_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, 0,
                                       igraph_vss_all(), 0, 0.85, &reset, 0, 0);
    IGRAPH_ASSERT(err == IGRAPH_EINVAL);

    igraph_vector_resize(&reset, 10);
    igraph_vector_fill(&reset, 0);
    err = igraph_personalized_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK,
                                       &res, 0, igraph_vss_all(), 0, 0.85,
                                       &reset, 0, &arpack_options);
    IGRAPH_ASSERT(err == IGRAPH_EINVAL);

    err = igraph_personalized_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK,
                                       &res, 0, igraph_vss_all(), 0, 0.85,
                                       &reset, 0, 0);
    IGRAPH_ASSERT(err == IGRAPH_EINVAL);

    igraph_vector_destroy(&reset);
    igraph_destroy(&g);
    igraph_set_error_handler(igraph_error_handler_abort);


    /* Special cases: check for empty graph */
    /* The ARPACK method has a special code path for this case */

    printf("\nEdgeless graph\n");

    igraph_empty(&g, 10, IGRAPH_UNDIRECTED);
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, 0, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, 0, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_destroy(&g);

    /* Special cases: check for empty graph, personalized */
    /* The ARPACK method has a special code path for this case */

    printf("\nEdgeless graph, personalized PageRank\n");

    igraph_empty(&g, 4, IGRAPH_UNDIRECTED);
    igraph_vector_init_seq(&reset, 1, 4);

    igraph_personalized_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, &reset, 0, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_personalized_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, &reset, 0, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    /* We also test the personalized case on a non-empty disconnected case,
     * which does not use the special code path for the ARPACK version.
     * The result must be the same for ARPACK and PRPACK. */

    printf("\nOne edge, two isolated vertices, personalized\n");
    igraph_add_edge(&g, 0, 1);

    igraph_personalized_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, &reset, 0, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_personalized_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, &reset, 0, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_vector_destroy(&reset);
    igraph_destroy(&g);

    /* Special cases: check for full graph, zero weights */
    /* The ARPACK method has a special code path for this case */

    printf("\nComplete graph, zero weights\n");

    igraph_full(&g, 10, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_init(&weights, 45);
    igraph_vector_fill(&weights, 0);
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, &weights, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, &weights, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* Another test case for PageRank (bug #792352) */

    printf("\nTest graph 3\n");

    igraph_small(&g, 9, IGRAPH_DIRECTED, 0, 5, 1, 5, 2, 0, 3, 1, 5, 4, 5, 7, 6, 0, 8, 0, 8, 1, -1);
    igraph_vector_init(&weights, 9);
    VECTOR(weights)[0] = 4;
    VECTOR(weights)[1] = 5;
    VECTOR(weights)[2] = 5;
    VECTOR(weights)[3] = 4;
    VECTOR(weights)[4] = 4;
    VECTOR(weights)[5] = 4;
    VECTOR(weights)[6] = 3;
    VECTOR(weights)[7] = 4;
    VECTOR(weights)[8] = 4;
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, &weights, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, &weights, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* Multigraph */

    printf("\nSmall undirected multigraph (unweighted)\n");

    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 1,2,
                 -1);

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, NULL, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, NULL, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    printf("\nSmall undirected multigraph (unit weights)\n");

    igraph_vector_init(&weights, 3);
    igraph_vector_fill(&weights, 1.0);

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, &weights, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, &weights, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    igraph_vector_destroy(&res);

    /* Graph with more than 127 vertices. PRPACK uses a different method above this size. */

    {
        igraph_vector_t edges_to_delete;
        igraph_vector_t res_arpack, res_prpack;
        igraph_vector_t weights;
        long int i, n;

        printf("\nLarge test graph, unweighted\n");

        /* 243 vertices, 729 edges */
        igraph_de_bruijn(&g, 3, 5);

        /* We delete some edges to break the symmetry of the graph.
         * Otherwise all vertices would have the same PageRank. */
        igraph_vector_init_seq(&edges_to_delete, 0, 37);
        igraph_delete_edges(&g, igraph_ess_vector(&edges_to_delete));
        igraph_vector_destroy(&edges_to_delete);

        /* Note: This test graph is not connected and has self-loops. */

        igraph_vector_init(&res_arpack, 0);
        igraph_vector_init(&res_prpack, 0);

        igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res_arpack, &value,
                        igraph_vss_all(), 1, 0.85, NULL, &arpack_options);
        IGRAPH_ASSERT(is_almost_one(value));

        igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res_prpack, &value,
                        igraph_vss_all(), 1, 0.85, NULL, NULL);
        IGRAPH_ASSERT(is_almost_one(value));

        n = igraph_vector_size(&res_arpack);
        for (i=0; i < n; ++i) {
            igraph_real_t ar = VECTOR(res_arpack)[i];
            igraph_real_t pr = VECTOR(res_prpack)[i];
            if (fabs(ar - pr) > 1e-12) {
                printf("Unexpected difference between ARPACK and PRPACK results for vertex %ld:\n"
                       "ARPACK: %g\n"
                       "PRPACK: %g\n"
                       "Difference: %g\n",
                       i, ar, pr, fabs(ar - pr));
            }
        }

        printf("\nLarge test graph, weighted\n");

        igraph_vector_init_seq(&weights, igraph_ecount(&g) + 1, 2*igraph_ecount(&g));

        igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res_arpack, &value,
                        igraph_vss_all(), 1, 0.85, &weights, &arpack_options);
        IGRAPH_ASSERT(is_almost_one(value));

        igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res_prpack, &value,
                        igraph_vss_all(), 1, 0.85, &weights, NULL);
        IGRAPH_ASSERT(is_almost_one(value));

        n = igraph_vector_size(&res_arpack);
        for (i=0; i < n; ++i) {
            igraph_real_t ar = VECTOR(res_arpack)[i];
            igraph_real_t pr = VECTOR(res_prpack)[i];
            if (fabs(ar - pr) > 1e-12) {
                printf("Unexpected difference between ARPACK and PRPACK results for vertex %ld:\n"
                       "ARPACK: %g\n"
                       "PRPACK: %g\n"
                       "Difference: %g\n",
                       i, ar, pr, fabs(ar - pr));
            }
        }

        igraph_vector_destroy(&weights);

        igraph_vector_destroy(&res_arpack);
        igraph_vector_destroy(&res_prpack);

        igraph_destroy(&g);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
