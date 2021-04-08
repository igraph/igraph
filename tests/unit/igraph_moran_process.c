/* -*- mode: C -*-  */
/*
  Test suite for the Moran process in a network setting.
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

#include "test_utilities.inc"

/* test parameters structure */
typedef struct {
    igraph_t *graph;
    igraph_vector_t *weights;
    igraph_vector_t *quantities;
    igraph_vector_t *strategies;
    igraph_neimode_t mode;
    int retval;
} strategy_test_t;

/* Error tests, i.e. we expect errors to be raised for each test.
 */
int error_tests() {
    igraph_t g, gzero, h;
    igraph_vector_t quant, quantnvert, quantzero;
    igraph_vector_t strat, stratnvert, stratzero;
    igraph_vector_t wgt, wgtnedge, wgtzero;
    int i, n, nvert, ret;
    strategy_test_t *test;

    igraph_empty(&h, 0, 0);  /* empty graph */
    /* nonempty graph */
    igraph_small(&g, /* n= */ 0, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, -1);
    nvert = igraph_vcount(&g);
    /* weights vectors */
    igraph_vector_init(&wgt, 0);
    igraph_vector_init(&wgtnedge, igraph_ecount(&g));
    /* quantities vectors */
    igraph_vector_init(&quant, 1);
    igraph_vector_init_real(&quantnvert, nvert, 0.1, 0.2, 0.3);
    /* strategies vectors */
    igraph_vector_init(&strat, 2);
    igraph_vector_init_real(&stratnvert, nvert, 0.0, 1.0, 2.0);

    igraph_small(&gzero, /* n= */ 0, IGRAPH_UNDIRECTED,
                 0, 3, 0, 4, 1, 2, 1, 4, 1, 5, 2, 3, 2, 4, 3, 4, -1);
    nvert = igraph_vcount(&gzero);
    igraph_vector_init(&quantzero, nvert);                /* vector of zeros */
    igraph_vector_init(&stratzero, nvert);                /* vector of zeros */
    igraph_vector_init(&wgtzero, igraph_ecount(&gzero));  /* vector of zeros */
    /* igraph_vector_init_real(&stratzero, nvert, 1.0, 0.0, 1.0, 2.0, 0.0, 3.0); */

    /* test parameters */
    /*------graph--weights--quantities--strategies--mode--retval------*/
    /* null pointer for graph */
    strategy_test_t null_graph = {NULL, NULL, NULL, NULL, IGRAPH_ALL, IGRAPH_EINVAL};
    /* null pointer for weights vector */
    strategy_test_t null_wgt = {&g, NULL, &quantnvert, &stratnvert, IGRAPH_ALL, IGRAPH_EINVAL};
    /* null pointer for quantities vector */
    strategy_test_t null_quant = {&g, &wgt, NULL, NULL, IGRAPH_ALL, IGRAPH_EINVAL};
    /* null pointer for strategies vector */
    strategy_test_t null_strat = {&g, &wgt, &quant, NULL, IGRAPH_ALL, IGRAPH_EINVAL};
    /* empty graph */
    strategy_test_t empty_graph = {&h, &wgt, &quant, &strat, IGRAPH_ALL, IGRAPH_EINVAL};
    /* length of quantities vector different from number of vertices */
    strategy_test_t qdiff_length = {&g, &wgtnedge, &quant, &strat, IGRAPH_ALL, IGRAPH_EINVAL};
    /* length of strategies vector different from number of vertices */
    strategy_test_t sdiff_length = {&g, &wgtnedge, &quantnvert, &strat, IGRAPH_ALL, IGRAPH_EINVAL};
    /* length of weights vector different from number of edges */
    strategy_test_t wdiff_length = {&g, &wgt, &quantnvert, &stratnvert, IGRAPH_ALL, IGRAPH_EINVAL};
    /* weights vector contains all zeros */
    strategy_test_t zero_wgt = {&g, &wgtnedge, &quantnvert, &stratnvert, IGRAPH_ALL, IGRAPH_EINVAL};
    /* quantities vector contains all zeros */
    strategy_test_t zero_quant = {&gzero, &wgtzero, &quantzero, &stratzero, IGRAPH_ALL, IGRAPH_EINVAL};
    strategy_test_t *all_checks[] = {/* 1 */  &null_graph,
                                              /* 2 */  &null_quant,
                                              /* 3 */  &null_strat,
                                              /* 4 */  &null_wgt,
                                              /* 5 */  &empty_graph,
                                              /* 6 */  &qdiff_length,
                                              /* 7 */  &sdiff_length,
                                              /* 8 */  &wdiff_length,
                                              /* 9 */  &zero_quant,
                                              /* 10 */ &zero_wgt
                                    };

    /* Run the error tests. We expect error to be raised for each test. */
    igraph_set_error_handler(igraph_error_handler_ignore);
    n = 10;
    i = 0;
    while (i < n) {
        test = all_checks[i];
        ret = igraph_moran_process(test->graph, test->weights, test->quantities,
                                   test->strategies, test->mode);
        if (ret != test->retval) {
            printf("Error test no. %d failed.\n", (int)(i + 1));
            return IGRAPH_FAILURE;
        }
        i++;
    }
    /* clean up */
    igraph_destroy(&g);
    igraph_destroy(&gzero);
    igraph_destroy(&h);
    igraph_vector_destroy(&quant);
    igraph_vector_destroy(&quantnvert);
    igraph_vector_destroy(&quantzero);
    igraph_vector_destroy(&strat);
    igraph_vector_destroy(&stratnvert);
    igraph_vector_destroy(&stratzero);
    igraph_vector_destroy(&wgt);
    igraph_vector_destroy(&wgtnedge);
    igraph_vector_destroy(&wgtzero);

    return IGRAPH_SUCCESS;
}

/* One iteration of the Moran process on a simple digraph.
 */
int moran_one_test() {
    igraph_t g;
    igraph_integer_t u = -1;  /* vertex chosen for reproduction */
    igraph_integer_t v = -1;  /* clone of u */
    igraph_integer_t nedge, nvert;
    igraph_real_t q = 0.0;
    igraph_vector_t quant, quantcp;
    igraph_vector_t strat, stratcp;
    igraph_vector_t wgt;
    long int i;

    /* graph representing the game network; quantities and strategies vectors */
    igraph_small(&g, /*nvert*/ 0, IGRAPH_DIRECTED,
                 0, 1, 0, 4, 1, 2, 1, 4, 2, 1, 3, 2, 4, 0, 4, 3, -1);
    nvert = igraph_vcount(&g);
    nedge = igraph_ecount(&g);
    igraph_vector_init_real(&quant, nvert, 0.77, 0.83, 0.64, 0.81, 0.05);
    igraph_vector_init_real(&strat, nvert, 2.0, 0.0, 0.0, 1.0, 2.0);
    /* Set the edge weights. Here we assume the following correspondence */
    /* between edge IDs and directed edges: */
    /* edge 0: 0 -> 1 */
    /* edge 1: 0 -> 4 */
    /* edge 2: 1 -> 2 */
    /* edge 3: 1 -> 4 */
    /* edge 4: 2 -> 1 */
    /* edge 5: 3 -> 2 */
    /* edge 6: 4 -> 0 */
    /* edge 7: 4 -> 3 */
    igraph_vector_init_real(&wgt, nedge, 1.9, 0.8, 6.2, 2.4, 1.1, 5.2, 7.3, 8.8);

    /* play game */
    igraph_vector_copy(&quantcp, &quant);
    igraph_vector_copy(&stratcp, &strat);
    igraph_moran_process(&g, &wgt, &quantcp, &stratcp, IGRAPH_OUT);

    /* Determine which vertex was chosen for death. The original quantities */
    /* vector contain unique values, i.e. no duplicates. Thus we compare the */
    /* updated quantities with the original one. */
    for (i = 0; i < igraph_vector_size(&quant); i++) {
        if (VECTOR(quant)[i] != VECTOR(quantcp)[i]) {
            /* found the new clone vertex */
            v = (igraph_integer_t)i;
            q = (igraph_real_t)VECTOR(quantcp)[i];
            break;
        }
    }
    IGRAPH_ASSERT(v >= 0);
    IGRAPH_ASSERT(q != 0.0);

    /* Now we know that v is a clone of some vertex. Determine the vertex that */
    /* v is a clone of. */
    for (i = 0; i < igraph_vector_size(&quant); i++) {
        if (VECTOR(quant)[i] == q) {
            /* found the vertex chosen for reproduction */
            u = (igraph_integer_t)i;
            break;
        }
    }
    IGRAPH_ASSERT(u >= 0);

    /* check that v is indeed a clone of u */
    if (VECTOR(quant)[u] != VECTOR(quantcp)[v]) {
        return IGRAPH_FAILURE;
    }
    if (VECTOR(strat)[u] != VECTOR(stratcp)[v]) {
        return IGRAPH_FAILURE;
    }

    igraph_destroy(&g);
    igraph_vector_destroy(&quant);
    igraph_vector_destroy(&quantcp);
    igraph_vector_destroy(&strat);
    igraph_vector_destroy(&stratcp);
    igraph_vector_destroy(&wgt);

    return IGRAPH_SUCCESS;
}

int main() {

    igraph_rng_seed(igraph_rng_default(), 1298);

    IGRAPH_ASSERT(error_tests() == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(moran_one_test() == IGRAPH_SUCCESS);

    VERIFY_FINALLY_STACK();

    return 0;
}
