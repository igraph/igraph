/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

void test_laplacian(igraph_t *g, const igraph_vector_t *w, igraph_bool_t dir, igraph_laplacian_normalization_t normalization) {
    igraph_matrix_t m;

    igraph_matrix_init(&m, 0, 0);

    igraph_get_laplacian(g, &m, IGRAPH_OUT, normalization, w);
    igraph_matrix_print(&m);

    igraph_matrix_destroy(&m);
}

int main() {
    igraph_t g_un, g_dir;
    igraph_vector_t weights;

    igraph_vector_init_int(&weights, 9, 1, 2, 3, 4, 5, 2, 2, 3, 3);
    igraph_small(&g_un, 6, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,4, 4,0, 1,1, 2,2, 1,2, 3,4, -1);
    igraph_small(&g_dir, 6, IGRAPH_DIRECTED, 0,1, 1,2, 2,3, 3,4, 4,0, 1,1, 2,2, 1,2, 3,4, -1);

    for (int i = 0; i < 8; i++) {
        igraph_bool_t is_normalized = i / 4;
        igraph_vector_t* v = ((i & 2) / 2 ? &weights : 0);
        igraph_bool_t dir = (i % 2 ? IGRAPH_DIRECTED : IGRAPH_UNDIRECTED);

        printf("=== %sormalized, %sweighted, %sdirected\n",
               (is_normalized ? "N" : "Unn"),
               (v != 0 ? "" : "un"),
               (dir == IGRAPH_DIRECTED ? "" : "un")
              );

        test_laplacian(dir ? &g_dir : &g_un, v, dir, is_normalized ? IGRAPH_LAPLACIAN_SYMMETRIC : IGRAPH_LAPLACIAN_UNNORMALIZED);
    }

    igraph_vector_destroy(&weights);
    igraph_destroy(&g_un);
    igraph_destroy(&g_dir);

    VERIFY_FINALLY_STACK();

    return 0;
}
