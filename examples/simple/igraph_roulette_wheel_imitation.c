/* -*- mode: C -*-  */
/*
  Test suite for stochastic imitation via roulette wheel selection.
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

#define R_INTEGER(a,b) (igraph_rng_get_integer(igraph_rng_default(), (a), (b)))

/* test parameters structure */
typedef struct {
    igraph_t *graph;
    igraph_integer_t vertex;
    igraph_bool_t islocal;
    igraph_vector_t *quantities;
    igraph_vector_int_t *strategies;
    igraph_vector_int_t *known_strats;
    igraph_neimode_t mode;
    igraph_error_t retval;
} strategy_test_t;

/* A game on a graph with 5 vertices and 7 edges. Use roulette wheel selection
 * to update strategies. This example also illustrates how a choice of
 * perspective (whether local or global) could affect the range of
 * possible strategies a vertex could adopt.
 */
igraph_error_t roulette_test(void) {
    igraph_t g;
    igraph_bool_t success;
    igraph_vector_t quant;
    igraph_vector_int_t *known, strat, stratcopy;
    igraph_vector_int_t known0, known1, known2, known3, known4, known5;
    igraph_integer_t i, k, n, nvert;
    igraph_error_t ret;
    strategy_test_t *test;

    /* the game network */
    igraph_small(&g, /*nvert=*/ 0, IGRAPH_UNDIRECTED,
                 0, 3, 0, 4, 1, 2, 1, 4, 1, 5, 2, 3, 2, 4, 3, 4, -1);
    nvert = igraph_vcount(&g);
    /* strategies vector; the strategy space is {0, 1, 2, 3} */
    /* V[i] is strategy of vertex i */
    igraph_vector_int_init_int(&strat, nvert, 1, 0, 1, 2, 0, 3);
    /* quantities vector; V[i] is quantity of vertex i */
    igraph_vector_init_real(&quant, nvert, 0.56, 0.13, 0.26, 0.73, 0.67, 0.82);
    /* possible strategies each vertex can adopt */
    igraph_vector_int_init_int(&known0, /*n=*/ 3, 0, 1, 2);       /* local */
    igraph_vector_int_init_int(&known1, /*n=*/ 3, 0, 1, 3);       /* local */
    igraph_vector_int_init_int(&known2, /*n=*/ 3, 0, 1, 2);       /* local */
    igraph_vector_int_init_int(&known3, /*n=*/ 3, 0, 1, 2);       /* local */
    igraph_vector_int_init_int(&known4, /*n=*/ 3, 0, 1, 2);       /* local */
    igraph_vector_int_init_int(&known5, /*n=*/ 4, 0, 1, 2, 3);    /* global */

    /* test parameters */
    /*graph--vert--islocal--quantities--strategies--known_strats--mode-retval*/
    strategy_test_t game0 = {&g, 0, 1, &quant, NULL, &known0, IGRAPH_ALL, IGRAPH_SUCCESS};
    strategy_test_t game1 = {&g, 1, 1, &quant, NULL, &known1, IGRAPH_ALL, IGRAPH_SUCCESS};
    strategy_test_t game2 = {&g, 2, 1, &quant, NULL, &known2, IGRAPH_ALL, IGRAPH_SUCCESS};
    strategy_test_t game3 = {&g, 3, 1, &quant, NULL, &known3, IGRAPH_ALL, IGRAPH_SUCCESS};
    strategy_test_t game4 = {&g, 4, 1, &quant, NULL, &known4, IGRAPH_ALL, IGRAPH_SUCCESS};
    strategy_test_t game5 = {&g, 5, 0, &quant, NULL, &known5, IGRAPH_ALL, IGRAPH_SUCCESS};
    strategy_test_t *all_checks[] = {/* 1 */ &game0,
                                             /* 2 */ &game1,
                                             /* 3 */ &game2,
                                             /* 4 */ &game3,
                                             /* 5 */ &game4,
                                             /* 6 */ &game5
                                    };

    /* play game */
    n = 6;
    i = 0;
    while (i < n) {
        test = all_checks[i];
        igraph_vector_int_init_copy(&stratcopy, &strat);
        ret = igraph_roulette_wheel_imitation(test->graph, test->vertex,
                                              test->islocal, test->quantities,
                                              &stratcopy, test->mode);
        if (ret != test->retval) {
            printf("Test no. %" IGRAPH_PRId " failed.\n", i + 1);
            return IGRAPH_FAILURE;
        }
        /* If the revised strategy s matches one of the candidate strategies, */
        /* then success. If s doesn't match any of the possible strategies, then */
        /* failure. Default to failure. */
        success = 0;
        known = test->known_strats;
        for (k = 0; k < igraph_vector_int_size(known); k++) {
            if (VECTOR(*known)[k] == VECTOR(stratcopy)[test->vertex]) {
                success = 1;
                break;
            }
        }
        if (!success) {
            printf("Roulette wheel imitation failed for vertex %" IGRAPH_PRId ".\n", test->vertex);
            return IGRAPH_FAILURE;
        }
        igraph_vector_int_destroy(&stratcopy);
        i++;
    }
    /* game finished; pack up */
    igraph_destroy(&g);
    igraph_vector_int_destroy(&known0);
    igraph_vector_int_destroy(&known1);
    igraph_vector_int_destroy(&known2);
    igraph_vector_int_destroy(&known3);
    igraph_vector_int_destroy(&known4);
    igraph_vector_int_destroy(&known5);
    igraph_vector_destroy(&quant);
    igraph_vector_int_destroy(&strat);

    return IGRAPH_SUCCESS;
}

/* It is possible for a vertex to retain its current strategy. This can
 * happen both in the local and global perspectives.
 */
igraph_error_t retain_strategy_test(void) {
    igraph_t g;
    igraph_integer_t max, min, v;
    igraph_vector_t quant;
    igraph_vector_int_t strat, stratcp;
    igraph_integer_t i, ntry, nvert;

    /* the game network */
    igraph_small(&g, /*nvert=*/ 0, IGRAPH_UNDIRECTED,
                 0, 3, 0, 4, 1, 2, 1, 4, 1, 5, 2, 3, 2, 4, 3, 4, -1);
    nvert = igraph_vcount(&g);
    /* strategies vector; the strategy space is {0, 1, 2, 3} */
    /* V[i] is strategy of vertex i */
    igraph_vector_int_init_int(&strat, nvert, 1, 0, 1, 2, 0, 3);
    /* quantities vector; V[i] is quantity of vertex i */
    igraph_vector_init_real(&quant, nvert, 0.56, 0.13, 0.26, 0.73, 0.67, 0.82);

    /* random vertex */
    min = 0;
    max = 5;
    igraph_rng_seed(igraph_rng_default(), 42); /* make tests deterministic */
    v = R_INTEGER(min, max);  /* min <= v <= max */
    /* Ensure that it is possible for v to retain its current strategy. We */
    /* will try to do this at most ntry times. As there are at most 6 vertices */
    /* to choose from, it shouldn't take long before we encounter a strategy */
    /* revision round where v retains its current strategy. */
    /* With local perspective. */
    i = 0;
    ntry = 100;
    igraph_vector_int_init(&stratcp, 0);
    do {
        i++;
        if (i > ntry) {
            return IGRAPH_FAILURE;    /* ideally this should never happen */
        }
        igraph_vector_int_destroy(&stratcp);
        igraph_vector_int_init_copy(&stratcp, &strat);
        igraph_roulette_wheel_imitation(&g, v, /*is local?*/ 1, &quant, &stratcp,
                                        IGRAPH_ALL);
    } while (VECTOR(stratcp)[v] != VECTOR(strat)[v]);
    /* If we get to this point, we know that there was an update round */
    /* i <= ntry as a result of which v retains its current strategy. */
    /* Now try again, but this time with the global perspective. */
    i = 0;
    do {
        i++;
        if (i > ntry) {
            return IGRAPH_FAILURE;    /* ideally this should never happen */
        }
        igraph_vector_int_destroy(&stratcp);
        igraph_vector_int_init_copy(&stratcp, &strat);
        igraph_roulette_wheel_imitation(&g, v, /*is local?*/ 0, &quant, &stratcp,
                                        IGRAPH_ALL);
    } while (VECTOR(stratcp)[v] != VECTOR(strat)[v]);
    /* nothing further to do, but housekeeping */
    igraph_destroy(&g);
    igraph_vector_destroy(&quant);
    igraph_vector_int_destroy(&strat);
    igraph_vector_int_destroy(&stratcp);

    return IGRAPH_SUCCESS;
}

int main(void) {
    igraph_error_t ret;

    igraph_rng_seed(igraph_rng_default(), 3241);

    ret = roulette_test();
    if (ret) {
        return IGRAPH_FAILURE;
    }
    ret = retain_strategy_test();
    if (ret) {
        return IGRAPH_FAILURE;
    }

    return IGRAPH_SUCCESS;
}
