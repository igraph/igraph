/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA, 02138 USA

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#define SIZE (1000)

int main() {

    igraph_matrix_t M, M2;
    igraph_vector_t lambda;
    igraph_matrix_t V;
    igraph_vector_t groups;
    igraph_vector_t ivec;
    int i, j;
    int n;

    igraph_rng_seed(igraph_rng_default(), 42);

    /* Symmetric matrix, exponentially distributed elements */

    igraph_matrix_init(&M, SIZE, SIZE);
    n = igraph_matrix_nrow(&M);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            MATRIX(M, i, j) = igraph_rng_get_exp(igraph_rng_default(), 1);
        }
    }
    igraph_matrix_init(&M2, n, n);
    igraph_matrix_update(&M2, &M);
    igraph_matrix_transpose(&M2);
    igraph_matrix_add(&M, &M2);
    igraph_matrix_scale(&M, 0.5);
    igraph_matrix_destroy(&M2);

    /* Get first (most positive) two eigenvectors */

    igraph_vector_init(&lambda, 0);
    igraph_matrix_init(&V, 0, 0);
    igraph_lapack_dsyevr(&M, IGRAPH_LAPACK_DSYEV_SELECT, /*vl=*/ 0, /*vu=*/ 0,
                         /*vestimate=*/ 0, /*il=*/ n - 1, /*iu=*/ n,
                         /*abstol=*/ 0.0, /*values=*/ &lambda, /*vectors=*/ &V,
                         /*support=*/ 0);

    /* Grouping */

    igraph_vector_init(&groups, 0);
    igraph_vector_init(&ivec, 2);
    VECTOR(ivec)[0] = 2;
    VECTOR(ivec)[1] = 3;
    igraph_scg_grouping(&V, &groups, /*invervals=*/ 0,
                        /*intervals_vector=*/ &ivec, IGRAPH_SCG_SYMMETRIC,
                        IGRAPH_SCG_OPTIMUM, /*p=*/ 0, /*maxiter=*/ 100);

    igraph_vector_print(&groups);

    igraph_vector_destroy(&ivec);
    igraph_vector_destroy(&groups);
    igraph_vector_destroy(&lambda);
    igraph_matrix_destroy(&V);
    igraph_matrix_destroy(&M);

    return 0;
}
