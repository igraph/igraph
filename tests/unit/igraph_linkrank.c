/* -*- mode: C -*-  */
/*
   IGraph library.
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
#include <math.h>
#include <float.h>

#include "test_utilities.h"

int is_almost_one(igraph_real_t x) {
    /* 2^5 = 32 is 5 binary digits  of tolerance */
    if (fabs(x - 1) > 32*DBL_EPSILON) {
        printf("Expected value to be 1, but actually got %15g.", x);
        return 0;
    }
    return 1;
}

int main(void) {

    igraph_t g;
    igraph_vector_t res, reset, weights;
    igraph_arpack_options_t arpack_options;
    igraph_real_t value;

    igraph_rng_seed(igraph_rng_default(), 137);
    igraph_arpack_options_init(&arpack_options);

    /* Test graph taken from PageRank tests - simple directed graph */
    printf("LinkRank Test 1 - Simple directed graph\n");
    igraph_small(&g, 4, IGRAPH_DIRECTED, 0, 1, 0, 2, 1, 0, 1, 3, 2, 3, 3, 2, -1);

    igraph_vector_init(&res, 0);

    /* Test unweighted LinkRank for all edges */
    igraph_linkrank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_ess_all(IGRAPH_EDGEORDER_ID), 1, 0.85, 0, 0);
    printf("PRPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    igraph_linkrank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_ess_all(IGRAPH_EDGEORDER_ID), 1, 0.85, 0, &arpack_options);
    printf("ARPACK: "); print_vector(&res);
    IGRAPH_ASSERT(is_almost_one(value));

    /* Test weighted LinkRank */
    printf("\nLinkRank Test 2 - Weighted edges\n");
    igraph_vector_init(&weights, 6);
    VECTOR(weights)[0] = 2.0;
    VECTOR(weights)[1] = 1.0;
    VECTOR(weights)[2] = 3.0;
    VECTOR(weights)[3] = 1.0;
    VECTOR(weights)[4] = 2.0;
    VECTOR(weights)[5] = 1.0;

    igraph_linkrank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_ess_all(IGRAPH_EDGEORDER_ID), 1, 0.85, &weights, 0);
    printf("PRPACK weighted: "); print_vector(&res);

    /* Test personalized LinkRank */
    printf("\nLinkRank Test 3 - Personalized LinkRank\n");
    igraph_vector_init(&reset, 4);
    VECTOR(reset)[0] = 1.0;
    VECTOR(reset)[1] = 0.0;
    VECTOR(reset)[2] = 0.0;
    VECTOR(reset)[3] = 0.0;

    igraph_personalized_linkrank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                                 igraph_ess_all(IGRAPH_EDGEORDER_ID), 1, 0.85, &reset, 0, 0);
    printf("PRPACK personalized: "); print_vector(&res);

    /* Test personalized LinkRank with vertex selector */
    printf("\nLinkRank Test 4 - Personalized LinkRank with vertex selector\n");
    igraph_personalized_linkrank_vs(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                                    igraph_ess_all(IGRAPH_EDGEORDER_ID), 1, 0.85, 
                                    igraph_vss_1(0), 0, 0);
    printf("PRPACK personalized vs: "); print_vector(&res);

    /* Test edge subset */
    printf("\nLinkRank Test 5 - Edge subset\n");
    igraph_es_t es;
    igraph_vector_int_t edge_vec;
    igraph_vector_int_init(&edge_vec, 2);
    VECTOR(edge_vec)[0] = 0;
    VECTOR(edge_vec)[1] = 3;
    igraph_es_vector(&es, &edge_vec);
    igraph_linkrank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    es, 1, 0.85, 0, 0);
    printf("PRPACK subset (edges 0,3): "); print_vector(&res);
    igraph_es_destroy(&es);
    igraph_vector_int_destroy(&edge_vec);

    /* Test star graph (undirected) */
    printf("\nLinkRank Test 6 - Undirected star\n");
    igraph_destroy(&g);
    igraph_star(&g, 5, IGRAPH_STAR_UNDIRECTED, 0);

    igraph_linkrank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_ess_all(IGRAPH_EDGEORDER_ID), 0, 0.85, 0, 0);
    printf("PRPACK star: "); print_vector(&res);

    /* Cleanup */
    igraph_vector_destroy(&res);
    igraph_vector_destroy(&reset);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}