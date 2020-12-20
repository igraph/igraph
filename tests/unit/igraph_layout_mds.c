/* -*- mode: C -*-  */
/*
   IGraph library.
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

#include "test_utilities.inc"

#define sqr(x) ((x)*(x))

int main() {
    igraph_t g;
    igraph_matrix_t coords, dist_mat;
    int i, j;

    srand(42); /* make tests deterministic */

    igraph_tree(&g, 10, 2, IGRAPH_TREE_UNDIRECTED);
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
    igraph_matrix_print(&coords);
    igraph_matrix_destroy(&coords);
    igraph_destroy(&g);

    igraph_full(&g, 8, IGRAPH_UNDIRECTED, 0);
    igraph_matrix_init(&coords, 8, 2);
    igraph_matrix_init(&dist_mat, 8, 8);
    for (i = 0; i < 8; i++)
        for (j = 0; j < 2; j++) {
            MATRIX(coords, i, j) = rand() % 1000;
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
                printf("dist(%d,%d) should be %.4f, but it is %.4f\n",
                       i, j, MATRIX(dist_mat, i, j), sqrt(dist_sq));
                return 1;
            }
        }
    igraph_matrix_destroy(&dist_mat);
    igraph_matrix_destroy(&coords);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
