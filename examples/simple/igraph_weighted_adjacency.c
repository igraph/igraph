/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2022  The igraph development team <igraph@igraph.org>

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
#include <stdarg.h>

int main(void) {
    igraph_matrix_t mat;
    igraph_t g;
    int m[4][4] = { { 0, 1, 2, 0 }, { 2, 0, 0, 1 }, { 0, 0, 1, 0 }, { 0, 1, 0, 0 } };
    igraph_vector_t weights;
    igraph_vector_int_t el;
    igraph_integer_t i, j, n;
    igraph_vector_int_init(&el, 0);
    igraph_vector_init(&weights, 0);

    igraph_matrix_init(&mat, 4, 4);
    for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
        MATRIX(mat, i, j) = m[i][j];
    }

    igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_DIRECTED, &weights, IGRAPH_LOOPS_ONCE);

    igraph_get_edgelist(&g, &el, 0);
    n = igraph_ecount(&g);

    for (i = 0, j = 0; i < n; i++, j += 2) {
        printf("%" IGRAPH_PRId " --> %" IGRAPH_PRId ": %g\n",
               VECTOR(el)[j], VECTOR(el)[j + 1], VECTOR(weights)[i]);
    }

    igraph_matrix_destroy(&mat);
    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&el);
    igraph_destroy(&g);

}
