/* -*- mode: C -*-  */
/*
  Test suite for stochastic imitation via uniform selection.
  Copyright (C) 2011 Minh Van Nguyen <nguyenminh2@gmail.com>

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
#include <stdio.h>

/* test parameters structure */
typedef struct {
    igraph_t *graph;
    igraph_integer_t vertex;
    igraph_imitate_algorithm_t algo;
    igraph_vector_t *quantities;
    igraph_vector_t *strategies;
    igraph_vector_t *known_strats;
    igraph_neimode_t mode;
    int retval;
} strategy_test_t;

/* Error tests. That is, we expect error codes to be returned from such tests.
 */
int error_tests() {
    igraph_t g, h;
    igraph_vector_t quant, strat;
    int i, n, ret;
    strategy_test_t *test;

    /* nonempty graph */
    igraph_small(&g, /*n vertices*/ 0, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, -1);
    igraph_empty(&h, 0, 0);         /* empty graph */
    igraph_vector_init(&quant, 1);  /* quantities vector */
    igraph_vector_init(&strat, 2);  /* strategies vector */

    /* test parameters */
    /*graph--vertex--algo--quantities--strategies--known_strats--mode--retval*/
    /* null pointer for graph */
    strategy_test_t null_graph = {NULL, 0, IGRAPH_IMITATE_BLIND, NULL, NULL, NULL, IGRAPH_ALL, IGRAPH_EINVAL};
    /* null pointer for quantities vector */
    strategy_test_t null_quant = {&g, 0, IGRAPH_IMITATE_BLIND, NULL, NULL, NULL, IGRAPH_ALL, IGRAPH_EINVAL};
    /* null pointer for strategies vector */
    strategy_test_t null_strat = {&g, 0, IGRAPH_IMITATE_BLIND, &quant, NULL, NULL, IGRAPH_ALL, IGRAPH_EINVAL};
    /* empty graph */
    strategy_test_t empty_graph = {&h, 0, IGRAPH_IMITATE_BLIND, &quant, &strat, NULL, IGRAPH_ALL, IGRAPH_EINVAL};
    /* length of quantities vector different from number of vertices */
    strategy_test_t qdiff_length = {&g, 0, IGRAPH_IMITATE_BLIND, &quant, &strat, NULL, IGRAPH_ALL, IGRAPH_EINVAL};
    /* length of strategies vector different from number of vertices */
    strategy_test_t sdiff_length = {&g, 0, IGRAPH_IMITATE_BLIND, &quant, &strat, NULL, IGRAPH_ALL, IGRAPH_EINVAL};
    strategy_test_t unknown_algo = {&g, 0, -1, &quant, &strat, NULL, IGRAPH_ALL, IGRAPH_EINVAL};
    strategy_test_t *all_checks[] = {/* 1 */ &null_graph,
                                             /* 2 */ &null_quant,
                                             /* 3 */ &null_strat,
                                             /* 4 */ &empty_graph,
                                             /* 5 */ &qdiff_length,
                                             /* 6 */ &sdiff_length,
                                             /* 7 */ &unknown_algo
                                    };
    /* Run the error tests. We expect error to be raised for each test. */
    igraph_set_error_handler(igraph_error_handler_ignore);
    n = 7;
    i = 0;
    while (i < n) {
        test = all_checks[i];
        ret = igraph_stochastic_imitation(test->graph, test->vertex, test->algo,
                                          test->quantities, test->strategies,
                                          test->mode);
        if (ret != test->retval) {
            printf("Error test no. %d failed.\n", (int)(i + 1));
            return IGRAPH_FAILURE;
        }
        i++;
    }
    /* clean up */
    igraph_destroy(&g);
    igraph_destroy(&h);
    igraph_vector_destroy(&quant);
    igraph_vector_destroy(&strat);

    return IGRAPH_SUCCESS;
}

/* Updating the strategy of an isolated vertex. In this case, the strategies
 * vector should not change at all.
 */
int isolated_vertex_test() {
    igraph_t g;
    igraph_vector_t quant, strat, v;
    int i, ret;

    /* graph with one isolated vertex */
    igraph_small(&g, /*n vertices*/ 0, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, -1);
    igraph_add_vertices(&g, 1, 0);  /* new vertex 3 is isolated */
    /* quantities vector: all vertices have the same fitness */
    igraph_vector_init_real(&quant, 4, 0.25, 0.25, 0.25, 0.25);
    /* strategies vector: 0 means aggressive strategy; 1 means passive */
    igraph_vector_init_real(&strat, 4, 1.0, 0.0, 1.0, 0.0);
    /* make a copy of the original strategies vector for comparison later on */
    igraph_vector_copy(&v, &strat);
    /* Now update strategy of vertex 3. Since this vertex is isolated, no */
    /* strategy update would take place. The resulting strategies vector */
    /* would be the same as it was originally. */
    ret = igraph_stochastic_imitation(/*graph*/ &g,
            /*vertex*/ 3,
            /*algorithm*/ IGRAPH_IMITATE_BLIND,
            /*quantities*/ &quant,
            /*strategies*/ &strat,
            /*mode*/ IGRAPH_ALL);
    if (ret) {
        printf("Isolated vertex test failed.\n");
        return IGRAPH_FAILURE;
    }
    for (i = 0; i < igraph_vector_size(&strat); i++) {
        if (VECTOR(strat)[i] != VECTOR(v)[i]) {
            printf("Isolated vertex test failed.\n");
            return IGRAPH_FAILURE;
        }
    }
    /* clean up */
    igraph_destroy(&g);
    igraph_vector_destroy(&quant);
    igraph_vector_destroy(&strat);
    igraph_vector_destroy(&v);

    return IGRAPH_SUCCESS;
}

/* A game on the Petersen graph. This graph has 10 vertices and 15 edges. The
 * Petersen graph is initialized with a default quantities vector and a
 * default strategies vector. Some vertices are chosen for strategy revision,
 * each one via a different stochastic imitation rule.
 */
int petersen_game_test() {
    igraph_t g;
    igraph_bool_t success;
    igraph_vector_t quant, strat, stratcopy, *knownstrats;
    igraph_vector_t known0, known2, known4;
    int i, k, n, nvert, ret;
    strategy_test_t *test;

    /* the Petersen graph */
    igraph_small(&g, /*n vertices*/ 0, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4, 0, 5, 1, 2, 1, 6, 2, 3, 2, 7, 3, 4, 3, 8, 4, 9,
                 5, 7, 5, 8, 6, 8, 6, 9, 7, 9, -1);
    nvert = igraph_vcount(&g);
    /* Strategies vector, one strategy for each vertex. Thus vec[i] is the */
    /* strategy of vertex i. The strategy space is: {0, 1, 2, 3}. */
    /* Each strategy should be an integer. */
    igraph_vector_init_real(&strat, nvert,
                            1.0, 1.0, 2.0, 2.0, 0.0,
                            0.0, 0.0, 1.0, 2.0, 3.0);
    /* Quantities vector, one quantity per vertex. Thus vec[i] is the */
    /* quantity for vertex i. */
    igraph_vector_init_real(&quant, nvert,
                            0.3, 1.1, 0.5, 1.0, 0.9,
                            0.8, 0.4, 0.1, 0.7, 0.7);
    /* parameter settings and known results */
    igraph_vector_init_real(&known0, 2, 0.0, 1.0);
    igraph_vector_init_real(&known2, 2, 1.0, 2.0);
    igraph_vector_init_real(&known4, 2, 0.0, 2.0);
    /*graph--vertex--algo--quantities--strategies--known_strats--mode--retval*/
    strategy_test_t blind0 = {&g, 0, IGRAPH_IMITATE_BLIND, &quant, NULL, &known0, IGRAPH_ALL, IGRAPH_SUCCESS};
    strategy_test_t augmented4 = {&g, 4, IGRAPH_IMITATE_AUGMENTED, &quant, NULL, &known4, IGRAPH_ALL, IGRAPH_SUCCESS};
    strategy_test_t contracted2 = {&g, 2, IGRAPH_IMITATE_CONTRACTED, &quant, NULL, &known2, IGRAPH_ALL, IGRAPH_SUCCESS};
    strategy_test_t *all_checks[] = {/* 1 */ &blind0,
                                             /* 2 */ &augmented4,
                                             /* 3 */ &contracted2
                                    };
    /* run the tests */
    n = 3;
    i = 0;
    while (i < n) {
        test = all_checks[i];
        igraph_vector_copy(&stratcopy, &strat);
        ret = igraph_stochastic_imitation(test->graph, test->vertex, test->algo,
                                          test->quantities, &stratcopy,
                                          test->mode);
        if (ret) {
            printf("Stochastic imitation failed for vertex %d.\n",
                   (int)test->vertex);
            return IGRAPH_FAILURE;
        }
        /* If the updated strategy for the vertex matches one of the known */
        /* strategies, then success. Default to failure. */
        success = 0;
        knownstrats = test->known_strats;
        for (k = 0; k < igraph_vector_size(knownstrats); k++) {
            if (VECTOR(*knownstrats)[k] == VECTOR(stratcopy)[test->vertex]) {
                success = 1;
                break;
            }
        }
        if (!success) {
            printf("Stochastic imitation failed for vertex %d.\n",
                   (int)test->vertex);
            return IGRAPH_FAILURE;
        }
        igraph_vector_destroy(&stratcopy);
        i++;
    }
    /* clean up */
    igraph_destroy(&g);
    igraph_vector_destroy(&known0);
    igraph_vector_destroy(&known2);
    igraph_vector_destroy(&known4);
    igraph_vector_destroy(&quant);
    igraph_vector_destroy(&strat);

    return IGRAPH_SUCCESS;
}

int main() {
    int ret;

    igraph_rng_seed(igraph_rng_default(), 547612);

    ret = error_tests();
    if (ret) {
        return ret;
    }
    ret = isolated_vertex_test();
    if (ret) {
        return ret;
    }
    ret = petersen_game_test();
    if (ret) {
        return ret;
    }

    return IGRAPH_SUCCESS;
}
