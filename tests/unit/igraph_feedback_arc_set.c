/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

void check_fas(const igraph_t *graph, const igraph_vector_int_t *fas) {
    igraph_t g;
    igraph_bool_t is_acyclic = false;

    igraph_copy(&g, graph);
    igraph_delete_edges(&g, igraph_ess_vector(fas));

    igraph_is_acyclic(&g, &is_acyclic);
    IGRAPH_ASSERT(is_acyclic);

    if (! igraph_is_directed(graph)) {
        igraph_int_t comp_no;
        igraph_connected_components(graph, NULL, NULL, &comp_no, IGRAPH_WEAK);
        IGRAPH_ASSERT(igraph_ecount(graph) - igraph_vector_int_size(fas) == igraph_vcount(graph) - comp_no);
    }

    igraph_destroy(&g);
}

igraph_real_t weight(const igraph_vector_t *weights, const igraph_vector_int_t *fas) {
    const igraph_int_t size = igraph_vector_int_size(fas);
    if (!weights) {
        return (igraph_real_t) size;
    } else {
        igraph_real_t total = 0;
        for (igraph_int_t i=0; i < size; i++) {
            total += VECTOR(*weights)[ VECTOR(*fas)[i] ];
        }
        return total;
    }
}

void compare_methods(const igraph_t *graph, const igraph_vector_t *weights) {
    igraph_vector_int_t fas1, fas2;

    igraph_vector_int_init(&fas1, 0);
    igraph_vector_int_init(&fas2, 0);

    igraph_feedback_arc_set(graph, &fas1, weights, IGRAPH_FAS_EXACT_IP_TI);
    check_fas(graph, &fas1);

    igraph_feedback_arc_set(graph, &fas2, weights, IGRAPH_FAS_EXACT_IP_CG);
    check_fas(graph, &fas2);

    /* Check that the two IP methods return solutions of the same size/weight. */
    IGRAPH_ASSERT(igraph_almost_equals(weight(weights, &fas1), weight(weights, &fas2), 1e-10));

    igraph_feedback_arc_set(graph, &fas1, weights, IGRAPH_FAS_APPROX_EADES);
    check_fas(graph, &fas1);

    /* Check that the exact IP method returns a solution that is at most as large
     * as the solution from the heuristic. */
    IGRAPH_ASSERT(igraph_cmp_epsilon(weight(weights, &fas1), weight(weights, &fas2), 1e-10) >= 0);

    igraph_vector_int_destroy(&fas2);
    igraph_vector_int_destroy(&fas1);
}

void rand_weights(const igraph_t *graph, igraph_vector_t *weights) {
    const igraph_int_t ecount = igraph_ecount(graph);
    igraph_vector_resize(weights, ecount);
    for (igraph_int_t i=0; i < ecount; i++) {
        VECTOR(*weights)[i] = RNG_UNIF01();
    }
}

void test_undirected(void) {
    igraph_t graph;
    igraph_vector_int_t result;

    igraph_vector_int_init(&result, 0);

    /* Null graph */
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Singleton graph */
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Two isolated vertices */
    igraph_empty(&graph, 2, IGRAPH_UNDIRECTED);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Single vertex with loop */
    igraph_small(&graph, 1, IGRAPH_UNDIRECTED,
                 0,0,
                 -1);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 1);
    IGRAPH_ASSERT(VECTOR(result)[0] == 0);
    check_fas(&graph, &result);
    igraph_destroy(&graph);

    /* Random graph */
    igraph_erdos_renyi_game_gnm(&graph, 20, 40, IGRAPH_UNDIRECTED, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    check_fas(&graph, &result);
    igraph_destroy(&graph);

    igraph_vector_int_destroy(&result);

    VERIFY_FINALLY_STACK();
}

void test_directed(void) {
    igraph_t graph;
    igraph_vector_int_t result;
    igraph_vector_t weights;

    igraph_vector_int_init(&result, 0);
    igraph_vector_init(&weights, 0);

    /* Null graph */
    igraph_empty(&graph, 0, IGRAPH_DIRECTED);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_EXACT_IP_CG);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_EXACT_IP_TI);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Singleton graph */
    igraph_empty(&graph, 1, IGRAPH_DIRECTED);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_EXACT_IP_CG);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_EXACT_IP_TI);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Two isolated vertices */
    igraph_empty(&graph, 2, IGRAPH_DIRECTED);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_EXACT_IP_CG);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_EXACT_IP_TI);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Single vertex with loop */
    igraph_small(&graph, 1, IGRAPH_DIRECTED,
                 0,0,
                 -1);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 1);
    IGRAPH_ASSERT(VECTOR(result)[0] == 0);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_EXACT_IP_CG);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 1);
    IGRAPH_ASSERT(VECTOR(result)[0] == 0);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_EXACT_IP_TI);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 1);
    IGRAPH_ASSERT(VECTOR(result)[0] == 0);
    igraph_destroy(&graph);

    /* Kautz graph */
    igraph_kautz(&graph, 2, 2);
    compare_methods(&graph, NULL);
    rand_weights(&graph, &weights);
    compare_methods(&graph, &weights);
    igraph_destroy(&graph);

    /* De Bruijn graph */
    igraph_de_bruijn(&graph, 3, 2);
    compare_methods(&graph, NULL);
    rand_weights(&graph, &weights);
    compare_methods(&graph, &weights);
    igraph_destroy(&graph);

    /* Multigraph with loops and isolated vertices */
    igraph_small(&graph, 0, IGRAPH_DIRECTED,
                 1,2, 1,2, 2,3, 3,4, 3,4, 4,1, 5,4, 1,5, 5,6, 6,5, 6,7, 3,3, 3,3, 8,8,
                 -1);
    compare_methods(&graph, NULL);
    rand_weights(&graph, &weights);
    compare_methods(&graph, &weights);
    igraph_destroy(&graph);

    /* Large cycle graph */
    igraph_ring(&graph, 30, IGRAPH_DIRECTED, false, true);
    compare_methods(&graph, NULL);
    rand_weights(&graph, &weights);
    compare_methods(&graph, &weights);
    igraph_destroy(&graph);

    /* Acyclic graph */
    igraph_erdos_renyi_game_gnm(&graph, 20, 30, IGRAPH_DIRECTED, IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    igraph_to_directed(&graph, IGRAPH_TO_DIRECTED_ACYCLIC);
    compare_methods(&graph, NULL);
    rand_weights(&graph, &weights);
    compare_methods(&graph, &weights);
    igraph_destroy(&graph);

    /* Several random graphs */
    for (igraph_int_t m=10; m <= 60; m += 10) {
        igraph_erdos_renyi_game_gnm(&graph, 20, m, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
        compare_methods(&graph, NULL);
        rand_weights(&graph, &weights);
        compare_methods(&graph, &weights);
        igraph_destroy(&graph);
    }

    /* Stress test CG method.
     * This should run in a reasonable amount of time. */

    igraph_erdos_renyi_game_gnm(&graph, 200, 400, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_EXACT_IP_CG);
    igraph_destroy(&graph);

    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&result);

    VERIFY_FINALLY_STACK();
}

int main(void) {

    igraph_rng_seed(igraph_rng_default(), 137);

    test_undirected();
    test_directed();

    return 0;
}
