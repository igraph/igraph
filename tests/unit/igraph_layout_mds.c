/*
   igraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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
#include <math.h>
#include <stdlib.h>

#include "test_utilities.h"

#define sqr(x) ((x)*(x))

static igraph_bool_t next_permutation(igraph_vector_int_t *perm) {
    igraph_int_t *a = VECTOR(*perm), n = igraph_vector_int_size(perm);
    igraph_int_t i = n - 2, j = n - 1, t;
    while (i >= 0 && a[i] >= a[i + 1]) i--;
    if (i < 0) return false;
    while (a[j] <= a[i]) j--;
    t = a[i]; a[i] = a[j]; a[j] = t;
    for (i++, j = n - 1; i < j; i++, j--) {
        t = a[i]; a[i] = a[j]; a[j] = t;
    }
    return true;
}

/* The doubly centered squared distance matrix of a path has rank 1, so for a 2D layout,
 * the second largest eigenvalue is a multiple eigenvalue 0. With BLIS as the BLAS,
 * DSYEVR fails to converge for some vertex orders of the path on 4 vertices. */
static void test_path_all_vertex_orders(void) {
    const igraph_int_t n = 4;
    igraph_t path, g;
    igraph_matrix_t coords, dist_mat;
    igraph_vector_int_t perm;

    igraph_path_graph(&path, n, IGRAPH_UNDIRECTED, false);
    igraph_matrix_init(&coords, 0, 0);
    igraph_matrix_init(&dist_mat, 0, 0);
    igraph_vector_int_init_range(&perm, 0, n);
    do {
        igraph_permute_vertices(&path, &g, &perm);
        igraph_layout_mds(&g, &coords, NULL, 2);
        /* A path embeds isometrically in a line, so the layout reproduces its distances. */
        igraph_distances(&g, NULL, &dist_mat, igraph_vss_all(), igraph_vss_all(), IGRAPH_ALL);
        for (igraph_int_t i = 0; i < n; i++) {
            for (igraph_int_t j = i + 1; j < n; j++) {
                double dist = sqrt(sqr(MATRIX(coords, i, 0) - MATRIX(coords, j, 0)) +
                                   sqr(MATRIX(coords, i, 1) - MATRIX(coords, j, 1)));
                IGRAPH_ASSERT(fabs(dist - MATRIX(dist_mat, i, j)) < 1e-8);
            }
        }
        igraph_destroy(&g);
    } while (next_permutation(&perm));

    igraph_vector_int_destroy(&perm);
    igraph_matrix_destroy(&dist_mat);
    igraph_matrix_destroy(&coords);
    igraph_destroy(&path);
}

int main(void) {
    igraph_t g;
    igraph_matrix_t coords, dist_mat;
    igraph_int_t i, j;

    igraph_rng_seed(igraph_rng_default(), 42); /* make tests deterministic */

    igraph_small(&g, 0, 0, -1);
    igraph_matrix_init(&coords, 0, 0);
    igraph_layout_mds(&g, &coords, 0, 2);
    print_matrix(&coords);
    igraph_matrix_destroy(&coords);
    igraph_destroy(&g);

    igraph_kary_tree(&g, 10, 2, IGRAPH_TREE_UNDIRECTED);
    igraph_matrix_init(&coords, 0, 0);
    igraph_layout_mds(&g, &coords, 0, 2);
    if (MATRIX(coords, 0, 0) > 0) {
        for (i = 0; i < igraph_matrix_nrow(&coords); i++) {
            MATRIX(coords, i, 0) *= -1;
        }
    }
    if (MATRIX(coords, 0, 1) < 0) {
        for (i = 0; i < igraph_matrix_nrow(&coords); i++) {
            MATRIX(coords, i, 1) *= -1;
        }
    }
    print_matrix(&coords);
    igraph_matrix_destroy(&coords);
    igraph_destroy(&g);

    igraph_full(&g, 8, IGRAPH_UNDIRECTED, 0);
    igraph_matrix_init(&coords, 8, 2);
    igraph_matrix_init(&dist_mat, 8, 8);
    for (i = 0; i < 8; i++)
        for (j = 0; j < 2; j++) {
            MATRIX(coords, i, j) = RNG_INTEGER(0, 1000);
        }
    for (i = 0; i < 8; i++)
        for (j = i + 1; j < 8; j++) {
            double dist_sq = 0.0;
            dist_sq += sqr(MATRIX(coords, i, 0) - MATRIX(coords, j, 0));
            dist_sq += sqr(MATRIX(coords, i, 1) - MATRIX(coords, j, 1));
            MATRIX(dist_mat, i, j) = sqrt(dist_sq);
            MATRIX(dist_mat, j, i) = sqrt(dist_sq);
        }
    igraph_layout_mds(&g, &coords, &dist_mat, 2);
    for (i = 0; i < 8; i++)
        for (j = i + 1; j < 8; j++) {
            double dist_sq = 0.0;
            dist_sq += sqr(MATRIX(coords, i, 0) - MATRIX(coords, j, 0));
            dist_sq += sqr(MATRIX(coords, i, 1) - MATRIX(coords, j, 1));
            if (fabs(sqrt(dist_sq) - MATRIX(dist_mat, i, j)) > 1e-2) {
                printf("dist(%" IGRAPH_PRId ", %" IGRAPH_PRId ") should be %.4f, but it is %.4f\n",
                       i, j, MATRIX(dist_mat, i, j), sqrt(dist_sq));
                return 1;
            }
        }
    igraph_matrix_destroy(&dist_mat);
    igraph_matrix_destroy(&coords);
    igraph_destroy(&g);

    test_path_all_vertex_orders();

    VERIFY_FINALLY_STACK();

    return 0;
}
