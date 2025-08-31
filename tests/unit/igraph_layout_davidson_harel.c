/*
   igraph library.
   Copyright (C) 2014-2024  The igraph development team <igraph@igraph.org>

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

#include "layout/layout_internal.h"

#include "test_utilities.h"

int intersect(void) {

    igraph_real_t negative[][8] = {
        { 1, 2, 2, 2, 1, 1, 2, 1 }, /* 1 */
        { 1, 2, 1, 1, 2, 2, 2, 1 }, /* 2 */
        { 1, 0, 0, 1, 2, 0, 3, 1 }, /* 3 */
        { 1, 0, 1, 1, 0, 2, 2, 2 }, /* 4 */
        { 1, 0, 1, 2, 3, 1, 3, 3 }, /* 5 */
        { 0, 0, 0, 2, 1, 1, 1, 2 }, /* 6 */
        { 0, 1, 1, 1, 2, 0, 2, 3 }, /* 7 */
        { 0, 0, 5, 0, 2, 1, 4, 3 }, /* 8 */
        { 0, 0, 5, 5, 3, 2, 3, 2 }  /* 9 */
    };

    igraph_real_t positive[][8] = {
        { 0, 1, 2, 1, 1, 0, 1, 2 }, /* 10 */
        { 0, 2, 5, 2, 1, 1, 4, 3 }, /* 11 */
        { 0, 0, 0, 3, 0, 1, 5, 1 }, /* 12 */
        { 0, 4, 2, 6, 0, 4, 2, 2 }  /* 13 */
    };
    /* { 1,1,1,1, 1,1,0,0 }, /\* 14 *\/ */
    /* { 0,0,1,1, 1,1,1,1 }, /\* 15 *\/ */
    /* { 0,0,2,2, 1,1,1,1 }}; /\* 16 *\/ */

    int no_neg = sizeof(negative) / sizeof(igraph_real_t) / 8;
    int no_pos = sizeof(positive) / sizeof(igraph_real_t) / 8;

    for (int i = 0; i < no_neg; i++) {
        igraph_real_t *co = negative[i];
        if (igraph_i_layout_segments_intersect(
                co[0], co[1], co[2], co[3],
                co[4], co[5], co[6], co[7])) {
            return i + 1;
        }
    }

    for (int i = 0; i < no_pos; i++) {
        igraph_real_t *co = positive[i];
        if (!igraph_i_layout_segments_intersect(
                co[0], co[1], co[2], co[3],
                co[4], co[5], co[6], co[7])) {
            return no_neg + i + 1;
        }
    }

    return 0;
}

int distance(void) {

    igraph_real_t configs[][7] = {
        { 1, 1, 2, 0, 2, 3, 1.0 }, /* 1 */
        { 1, 1, 1, 0, 1, 3, 0.0 }, /* 2 */
        { 1, 1, 0, 1, 1, 0, 0.5 }, /* 3 */
        { 1, 2, 0, 0, 0, 1, 2.0 }, /* 4 */
        { 1, 0, 0, 1, 0, 2, 2.0 }, /* 5 */
        { 0, 0, 1, 1, 1, 2, 2.0 }, /* 6 */
        { 0, 3, 1, 1, 1, 2, 2.0 }  /* 7 */
    };

    int no = sizeof(configs) / sizeof(igraph_real_t) / 8;

    for (int i = 0; i < no; i++) {
        igraph_real_t *co = configs[i];
        igraph_real_t res = igraph_i_layout_point_segment_dist2(
            co[0], co[1], co[2], co[3], co[4], co[5]);
        if (fabs(res - co[6]) > 1e-12) {
            printf("%g\n", (double) res);
            return i + 1;
        }
    }

    return 0;
}

void check_layout_davidson_harel(void) {
    igraph_t g;
    igraph_matrix_t res;
    igraph_bool_t use_seed;
    igraph_int_t maxiter, fineiter;
    igraph_real_t cool_fact, weight_node_dist, weight_border;
    igraph_real_t weight_edge_lengths, weight_edge_crossings, weight_node_edge_dist;

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_matrix_init(&res, 0, 0);

    use_seed = 0;
    maxiter = 5;
    fineiter = 5;
    cool_fact = 0.75;
    weight_node_dist = 1.0;
    weight_border = 0.1;
    weight_edge_lengths = 0.5;
    weight_edge_crossings = 0.8;
    weight_node_edge_dist = 0.2;

    printf("Checking graph with no vertices\n");
    igraph_small(&g, 0, IGRAPH_DIRECTED, -1);
    igraph_layout_davidson_harel(&g, &res, use_seed, maxiter, fineiter, cool_fact,
            weight_node_dist, weight_border,
            weight_edge_lengths,
            weight_edge_crossings,
            weight_node_edge_dist);

    IGRAPH_ASSERT(igraph_matrix_nrow(&res) == 0);
    igraph_destroy(&g);

    printf("Checking graph with 10 vertices, no edges\n");
    igraph_small(&g, 10, IGRAPH_DIRECTED, -1);
    igraph_layout_davidson_harel(&g, &res, use_seed, maxiter, fineiter, cool_fact,
            weight_node_dist, weight_border,
            weight_edge_lengths,
            weight_edge_crossings,
            weight_node_edge_dist);
    IGRAPH_ASSERT(igraph_matrix_nrow(&res) == 10);
    IGRAPH_ASSERT(igraph_matrix_ncol(&res) == 2);
    IGRAPH_ASSERT(igraph_matrix_max(&res) < 20);
    IGRAPH_ASSERT(igraph_matrix_min(&res) > -20);
    igraph_destroy(&g);

    printf("Checking full graph with 10 vertices\n");
    igraph_full(&g, 10, IGRAPH_DIRECTED, 1);
    igraph_layout_davidson_harel(&g, &res, use_seed, maxiter, fineiter, cool_fact,
            weight_node_dist, weight_border,
            weight_edge_lengths,
            weight_edge_crossings,
            weight_node_edge_dist);
    IGRAPH_ASSERT(igraph_matrix_nrow(&res) == 10);
    IGRAPH_ASSERT(igraph_matrix_ncol(&res) == 2);
    IGRAPH_ASSERT(igraph_matrix_max(&res) < 20);
    IGRAPH_ASSERT(igraph_matrix_min(&res) > -20);
    igraph_destroy(&g);
    igraph_matrix_destroy(&res);
}

int main(void) {
    int res1, res2;

    res1 = intersect();
    if (res1 != 0) {
        printf("Unexpected result from intersect(), config %d.\n", res1);
        return res1;
    }
    res2 = distance() ;
    if (res2 != 0) {
        printf("Unexpected result from distance(), config %d.\n", res2);
        return res2;
    }

    check_layout_davidson_harel();

    VERIFY_FINALLY_STACK();

    return 0;
}
