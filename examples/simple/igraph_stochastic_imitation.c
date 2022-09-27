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
    igraph_vector_int_t *strategies;
    igraph_vector_int_t *known_strats;
    igraph_neimode_t mode;
    igraph_integer_t retval;
} strategy_test_t;

/* Updating the strategy of an isolated vertex. In this case, the strategies
 * vector should not change at all.
 */
igraph_error_t isolated_vertex_test(void) {
    igraph_t g;
    igraph_vector_t quant;
    igraph_vector_int_t strat, v;
    igraph_integer_t i;
    igraph_error_t ret;

    /* graph with one isolated vertex */
    igraph_small(&g, /*n vertices*/ 0, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, -1);
    igraph_add_vertices(&g, 1, 0);  /* new vertex 3 is isolated */
    /* quantities vector: all vertices have the same fitness */
    igraph_vector_init_real(&quant, 4, 0.25, 0.25, 0.25, 0.25);
    /* strategies vector: 0 means aggressive strategy; 1 means passive */
    igraph_vector_int_init_int(&strat, 4, 1, 0, 1, 0);
    /* make a copy of the original strategies vector for comparison later on */
    igraph_vector_int_init_copy(&v, &strat);
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
    for (i = 0; i < igraph_vector_int_size(&strat); i++) {
        if (VECTOR(strat)[i] != VECTOR(v)[i]) {
            printf("Isolated vertex test failed.\n");
            return IGRAPH_FAILURE;
        }
    }
    /* clean up */
    igraph_destroy(&g);
    igraph_vector_destroy(&quant);
    igraph_vector_int_destroy(&strat);
    igraph_vector_int_destroy(&v);

    return IGRAPH_SUCCESS;
}

/* A game on the Petersen graph. This graph has 10 vertices and 15 edges. The
 * Petersen graph is initialized with a default quantities vector and a
 * default strategies vector. Some vertices are chosen for strategy revision,
 * each one via a different stochastic imitation rule.
 */
igraph_error_t petersen_game_test(void) {
    igraph_t g;
    igraph_bool_t success;
    igraph_vector_t quant;
    igraph_vector_int_t strat, stratcopy, *knownstrats;
    igraph_vector_int_t known0, known2, known4;
    igraph_integer_t i, k, n;
    igraph_error_t ret;
    int nvert;
    strategy_test_t *test;

    /* the Petersen graph */
    igraph_small(&g, /*n vertices*/ 0, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4, 0, 5, 1, 2, 1, 6, 2, 3, 2, 7, 3, 4, 3, 8, 4, 9,
                 5, 7, 5, 8, 6, 8, 6, 9, 7, 9, -1);
    nvert = igraph_vcount(&g);
    /* Strategies vector, one strategy for each vertex. Thus vec[i] is the */
    /* strategy of vertex i. The strategy space is: {0, 1, 2, 3}. */
    /* Each strategy should be an integer. */
    igraph_vector_int_init_int(&strat, nvert, 1, 1, 2, 2, 0, 0, 0, 1, 2, 3);
    /* Quantities vector, one quantity per vertex. Thus vec[i] is the */
    /* quantity for vertex i. */
    igraph_vector_init_real(&quant, nvert,
                            0.3, 1.1, 0.5, 1.0, 0.9,
                            0.8, 0.4, 0.1, 0.7, 0.7);
    /* parameter settings and known results */
    igraph_vector_int_init_int(&known0, 2, 0, 1);
    igraph_vector_int_init_int(&known2, 2, 1, 2);
    igraph_vector_int_init_int(&known4, 2, 0, 2);
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
        igraph_vector_int_init_copy(&stratcopy, &strat);
        ret = igraph_stochastic_imitation(test->graph, test->vertex, test->algo,
                                          test->quantities, &stratcopy,
                                          test->mode);
        if (ret) {
            printf("Stochastic imitation failed for vertex %" IGRAPH_PRId ".\n", test->vertex);
            return IGRAPH_FAILURE;
        }
        /* If the updated strategy for the vertex matches one of the known */
        /* strategies, then success. Default to failure. */
        success = 0;
        knownstrats = test->known_strats;
        for (k = 0; k < igraph_vector_int_size(knownstrats); k++) {
            if (VECTOR(*knownstrats)[k] == VECTOR(stratcopy)[test->vertex]) {
                success = 1;
                break;
            }
        }
        if (!success) {
            printf("Stochastic imitation failed for vertex %" IGRAPH_PRId ".\n", test->vertex);
            return IGRAPH_FAILURE;
        }
        igraph_vector_int_destroy(&stratcopy);
        i++;
    }
    /* clean up */
    igraph_destroy(&g);
    igraph_vector_int_destroy(&known0);
    igraph_vector_int_destroy(&known2);
    igraph_vector_int_destroy(&known4);
    igraph_vector_destroy(&quant);
    igraph_vector_int_destroy(&strat);

    return IGRAPH_SUCCESS;
}

int main(void) {
    igraph_error_t ret;

    igraph_rng_seed(igraph_rng_default(), 547612);

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
