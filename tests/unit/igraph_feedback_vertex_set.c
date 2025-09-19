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

void check_fvs(const igraph_t *graph, const igraph_vector_int_t *fvs) {
    igraph_t g;
    igraph_bool_t is_acyclic = false;

    igraph_copy(&g, graph);
    igraph_delete_vertices(&g, igraph_vss_vector(fvs));

    igraph_is_acyclic(&g, &is_acyclic);
    IGRAPH_ASSERT(is_acyclic);

    igraph_destroy(&g);
}

igraph_real_t weight(const igraph_vector_t *weights, const igraph_vector_int_t *fvs) {
    const igraph_int_t size = igraph_vector_int_size(fvs);
    if (!weights) {
        return (igraph_real_t) size;
    } else {
        igraph_real_t total = 0;
        for (igraph_int_t i=0; i < size; i++) {
            total += VECTOR(*weights)[ VECTOR(*fvs)[i] ];
        }
        return total;
    }
}

void compare_methods(const igraph_t *graph, const igraph_vector_t *weights) {
    igraph_vector_int_t fvs;

    igraph_vector_int_init(&fvs, 0);

    igraph_feedback_vertex_set(graph, &fvs, weights, IGRAPH_FVS_EXACT_IP);
    check_fvs(graph, &fvs);

    igraph_vector_int_destroy(&fvs);
    igraph_vector_int_destroy(&fvs);
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
    igraph_vector_t weights;

    igraph_vector_int_init(&result, 0);
    igraph_vector_init(&weights, 0);

    /* Null graph */
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_feedback_vertex_set(&graph, &result, NULL, IGRAPH_FVS_EXACT_IP);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Singleton graph */
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_feedback_vertex_set(&graph, &result, NULL, IGRAPH_FVS_EXACT_IP);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Two isolated vertices */
    igraph_empty(&graph, 2, IGRAPH_UNDIRECTED);
    igraph_feedback_vertex_set(&graph, &result, NULL, IGRAPH_FVS_EXACT_IP);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Single vertex with loop */
    igraph_small(&graph, 1, IGRAPH_UNDIRECTED,
                 0,0,
                 -1);
    igraph_feedback_vertex_set(&graph, &result, NULL, IGRAPH_FVS_EXACT_IP);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 1);
    IGRAPH_ASSERT(VECTOR(result)[0] == 0);
    check_fvs(&graph, &result);
    igraph_destroy(&graph);

    /* Small graph */
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 0, 2, 2, 3, 1, 3, 0, 4, 2, 4, -1);
    igraph_feedback_vertex_set(&graph, &result, NULL, IGRAPH_FVS_EXACT_IP);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 1);
    IGRAPH_ASSERT(VECTOR(result)[0] == 2);
    check_fvs(&graph, &result);
    igraph_destroy(&graph);

    igraph_small(&graph, 5, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 0, 2, 3, 3, 1, 0, 4, 4, 2, 4, 3, -1);
    igraph_vector_resize(&weights, 5);
    VECTOR(weights)[0] = 3;
    VECTOR(weights)[1] = 3;
    VECTOR(weights)[2] = 2;
    VECTOR(weights)[3] = 2;
    VECTOR(weights)[4] = 1;
    igraph_feedback_vertex_set(&graph, &result, &weights, IGRAPH_FVS_EXACT_IP);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 2);
    igraph_vector_int_sort(&result);
    IGRAPH_ASSERT(VECTOR(result)[0] == 2);
    IGRAPH_ASSERT(VECTOR(result)[1] == 4);
    igraph_destroy(&graph);

    /* Random graph */
    igraph_erdos_renyi_game_gnm(&graph, 10, 20, IGRAPH_UNDIRECTED, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    igraph_feedback_vertex_set(&graph, &result, NULL, IGRAPH_FVS_EXACT_IP);
    check_fvs(&graph, &result);
    igraph_destroy(&graph);

    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&result);

    VERIFY_FINALLY_STACK();
}

void test_directed(void) {
    igraph_t graph;
    igraph_vector_int_t result;

    igraph_vector_int_init(&result, 0);

    /* Null graph */
    igraph_empty(&graph, 0, IGRAPH_DIRECTED);
    igraph_feedback_vertex_set(&graph, &result, NULL, IGRAPH_FVS_EXACT_IP);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Singleton graph */
    igraph_empty(&graph, 1, IGRAPH_DIRECTED);
    igraph_feedback_vertex_set(&graph, &result, NULL, IGRAPH_FVS_EXACT_IP);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Two isolated vertices */
    igraph_empty(&graph, 2, IGRAPH_DIRECTED);
    igraph_feedback_vertex_set(&graph, &result, NULL, IGRAPH_FVS_EXACT_IP);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Single vertex with loop */
    igraph_small(&graph, 1, IGRAPH_DIRECTED,
                 0,0,
                 -1);
    igraph_feedback_vertex_set(&graph, &result, NULL, IGRAPH_FVS_EXACT_IP);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 1);
    IGRAPH_ASSERT(VECTOR(result)[0] == 0);
    check_fvs(&graph, &result);
    igraph_destroy(&graph);

    /* Small graph */
    igraph_small(&graph, 6, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 0, 2, 3, 3, 1, 0, 4, 4, 2, 4, 3, -1);
    igraph_feedback_vertex_set(&graph, &result, NULL, IGRAPH_FVS_EXACT_IP);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 1);
    IGRAPH_ASSERT(VECTOR(result)[0] == 2);
    check_fvs(&graph, &result);
    igraph_destroy(&graph);

    /* Kautz graph, see https://doi.org/10.1016/j.disc.2006.09.010 Fig. 1 */
    igraph_kautz(&graph, 2, 2);
    igraph_feedback_vertex_set(&graph, &result, NULL, IGRAPH_FVS_EXACT_IP);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 5);
    check_fvs(&graph, &result);
    igraph_destroy(&graph);

    /* Random graph */
    igraph_erdos_renyi_game_gnm(&graph, 10, 20, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    igraph_feedback_vertex_set(&graph, &result, NULL, IGRAPH_FVS_EXACT_IP);
    check_fvs(&graph, &result);
    igraph_destroy(&graph);

    igraph_vector_int_destroy(&result);

    VERIFY_FINALLY_STACK();
}

int main(void) {

    igraph_rng_seed(igraph_rng_default(), 137);

    test_undirected();
    test_directed();

    return 0;
}
