/*
   igraph library.
   Copyright (C) 2021-2022  The igraph development team <igraph@igraph.org>

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

#include "test_utilities.h"

/* Helper functions for hypothesis testing for binomial and negative binomial distr. */

igraph_real_t log_binom_prob(igraph_int_t n, igraph_real_t p, igraph_int_t k) {
    return lgamma(1+n) - lgamma(1+k) - lgamma(1 - k + n) +
           k * log(p) + (n-k) * log(1-p);
}

igraph_real_t log_neg_binom_prob(igraph_int_t n, igraph_real_t p, igraph_int_t k) {
    return lgamma(k + n) - lgamma(1+k) - lgamma(n) +
           k * log(1-p) + n * log(p);
}

igraph_real_t binom_test(igraph_int_t n, igraph_real_t p, igraph_int_t k0) {
    igraph_real_t res = 0;
    igraph_real_t lP0 = log_binom_prob(n, p, k0);
    igraph_real_t last_lP = -IGRAPH_INFINITY;
    for (igraph_int_t k=0; k < n; k++) {
        igraph_real_t lP = log_binom_prob(n, p, k);
        if (lP <= lP0) {
            res += exp(lP);
            /* stop when relative change to 'res' is < 10^-12 */
            if (lP < last_lP && lP - log(res) < -12.0 * log(10)) break;
        }
        last_lP = lP;
    }
    return res;
}

igraph_real_t neg_binom_test(igraph_int_t n, igraph_real_t p, igraph_int_t k0) {
    igraph_real_t res = 0;
    igraph_real_t lP0 = log_neg_binom_prob(n, p, k0);
    igraph_real_t last_lP = -IGRAPH_INFINITY;
    for (igraph_int_t k=0; ; k++) {
        igraph_real_t lP = log_neg_binom_prob(n, p, k);
        if (lP <= lP0) {
            res += exp(lP);
            /* stop when relative change to 'res' is < 10^-12 */
            if (lP < last_lP && lP - log(res) < -12.0 * log(10)) break;
        }
        last_lP = lP;
    }
    return res;
}

void stress_test(void) {
    igraph_rng_seed(igraph_rng_default(), 137);

    for (igraph_int_t size=2; size < 5; size++) {
        for (igraph_int_t i=0; i < 100; i++) {
            igraph_t g;
            igraph_bool_t simple;

            igraph_erdos_renyi_game_gnp(&g, size, 0.5, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

            igraph_is_simple(&g, &simple, IGRAPH_DIRECTED);
            if (! simple) {
                printf("Erdos-Renyi GNP graph is not simple! size=%" IGRAPH_PRId ", i=%" IGRAPH_PRId ".\n",
                       size, i);
                print_graph(&g);
            }
            IGRAPH_ASSERT(simple);

            igraph_destroy(&g);
        }
    }

    for (igraph_int_t size=2; size < 5; size++) {
        for (igraph_int_t i=0; i < 100; i++) {
            igraph_t g;
            igraph_bool_t simple;

            igraph_erdos_renyi_game_gnp(&g, size, 0.5, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

            igraph_is_simple(&g, &simple, IGRAPH_DIRECTED);
            if (! simple) {
                printf("Erdos-Renyi GNP graph is not simple! size=%" IGRAPH_PRId ", i=%" IGRAPH_PRId ".\n",
                       size, i);
                print_graph(&g);
            }
            IGRAPH_ASSERT(simple);

            igraph_destroy(&g);
        }
    }

    VERIFY_FINALLY_STACK();
}

void check_gnp(
        igraph_int_t n, igraph_real_t p,
        igraph_bool_t directed,
        igraph_bool_t loops, igraph_bool_t multiple,
        igraph_bool_t edge_labeled) {

    igraph_t graph;
    igraph_bool_t has_loop, has_multi;
    igraph_edge_type_sw_t allowed_edge_types;

    allowed_edge_types = IGRAPH_SIMPLE_SW;
    if (loops) allowed_edge_types |= IGRAPH_LOOPS_SW;
    if (multiple) allowed_edge_types |= IGRAPH_MULTI_SW;

    igraph_erdos_renyi_game_gnp(&graph, n, p, directed, allowed_edge_types, edge_labeled);

    IGRAPH_ASSERT(igraph_is_directed(&graph) == directed);
    IGRAPH_ASSERT(igraph_vcount(&graph) == n);

    igraph_has_loop(&graph, &has_loop);
    igraph_has_multiple(&graph, &has_multi);

    if (!multiple) IGRAPH_ASSERT(!has_multi);
    if (!loops) IGRAPH_ASSERT(!has_loop);

    igraph_real_t complete_ecount;

    if (directed) {
        complete_ecount = loops ? (igraph_real_t) n * n : (igraph_real_t)  n * (n-1.0);
    } else {
        complete_ecount = loops ? (igraph_real_t)  n * (n+1.0) / 2.0 : (igraph_real_t)  n * (n-1) / 2.0;
    }

    if (p == 0 || complete_ecount == 0) {
        IGRAPH_ASSERT(igraph_ecount(&graph) == 0);
    } else if (!multiple && p == 1) {
        IGRAPH_ASSERT(igraph_ecount(&graph) == complete_ecount);
    } else {
        /* Edge count check using hypothesis testing.
         * With simple graphs, the edge count follows a binomial distribution,
         * with multigraphs, a negative binomial distribution. In both cases,
         * we perform an exact two-tailed test unless the graph is too large.
         */

        igraph_real_t pval = 1;
        igraph_real_t m = complete_ecount * p;
        igraph_real_t sd = multiple ? sqrt(m * (1+p)) : sqrt(m * (1-p));
        igraph_real_t dev = (m - igraph_ecount(&graph)) / sd;
        igraph_bool_t do_test = m < 10000 && ! edge_labeled;

        if (do_test) {
            if (multiple) {
                pval = neg_binom_test(complete_ecount, 1 / (1 + p), igraph_ecount(&graph));
            } else {
                pval = binom_test(complete_ecount, p, igraph_ecount(&graph));
            }
            printf("p-value=%.3f ", pval);
        } else {
            printf("p-value= ?    ");
        }

        /* Output stats, including the actual and expected edge counts,
         * as well as their difference in standard deviations. */
        printf("for n=%6" IGRAPH_PRId ", p=%5g, %s, %s, %s; "
               "ecount: %6" IGRAPH_PRId ", expected: %8.1f, % .2f SD\n",
               n, p,
               directed ? "  directed" : "undirected",
               loops    ?   "   loops" : "no-loops",
               multiple ?   "   multi" : "no-multi",
               igraph_ecount(&graph),
               m, dev);

        /* Failure is unlikely, but possible. If it happens, check manually. */
        if (do_test) {
            IGRAPH_ASSERT(pval > 0.01);
        } else if (!edge_labeled) {
            /* If the graph is too large for the exact test to run efficiently,
             * use a normal approximation, even though this is inaccurate with
             * small p. */
            IGRAPH_ASSERT(fabs(dev) < 2.576);
        }
    }

    igraph_destroy(&graph);
}

/* Check all parameter combinations */
void check_all_gnp(igraph_int_t n, igraph_real_t p) {

    if (p <= 1) {
        check_gnp(n, p, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE, IGRAPH_EDGE_UNLABELED);
        check_gnp(n, p, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE, IGRAPH_EDGE_UNLABELED);
        check_gnp(n, p, IGRAPH_UNDIRECTED, IGRAPH_LOOPS, IGRAPH_NO_MULTIPLE, IGRAPH_EDGE_UNLABELED);
        check_gnp(n, p, IGRAPH_DIRECTED, IGRAPH_LOOPS, IGRAPH_NO_MULTIPLE, IGRAPH_EDGE_UNLABELED);
    }

    check_gnp(n, p, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE, IGRAPH_EDGE_UNLABELED);
    check_gnp(n, p, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE, IGRAPH_EDGE_UNLABELED);
    check_gnp(n, p, IGRAPH_UNDIRECTED, IGRAPH_LOOPS, IGRAPH_MULTIPLE, IGRAPH_EDGE_UNLABELED);
    check_gnp(n, p, IGRAPH_DIRECTED, IGRAPH_LOOPS, IGRAPH_MULTIPLE, IGRAPH_EDGE_UNLABELED);

    check_gnp(n, p, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE, IGRAPH_EDGE_LABELED);
    check_gnp(n, p, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE, IGRAPH_EDGE_LABELED);
    check_gnp(n, p, IGRAPH_UNDIRECTED, IGRAPH_LOOPS, IGRAPH_MULTIPLE, IGRAPH_EDGE_LABELED);
    check_gnp(n, p, IGRAPH_DIRECTED, IGRAPH_LOOPS, IGRAPH_MULTIPLE, IGRAPH_EDGE_LABELED);
}

void test_examples(void) {

    /* Ensure that the test is deterministic */
    igraph_rng_seed(igraph_rng_default(), 42);

    check_all_gnp(0, 0.0);

    /* Empty graph */

    check_all_gnp(10, 0.0);

    /* Singleton, with loop if allowed */

    check_all_gnp(1, 1.0);
    check_all_gnp(1, 5.0);

    /* Complete graph */

    check_all_gnp(10, 1.0);

    /* Random graph */

    check_all_gnp(10, 0.5);
    check_all_gnp(10, 12.0);

    /* Create a couple of larger graphs too */

    check_all_gnp(100, 1.0 / 100);
    check_all_gnp(100, 2.0 / 100);
    check_all_gnp(100, 80.0 / 100);
    check_all_gnp(1000, 2.0 / 1000);
    check_all_gnp(100000, 2.0 / 100000);

    VERIFY_FINALLY_STACK();
}

int main(void) {

    test_examples();
    stress_test();

    return 0;
}
