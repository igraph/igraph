/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

#include "test_utilities.inc"

int main() {
    igraph_t g, g_rev, g_test;
    igraph_bool_t same;
    igraph_matrix_t W;
    int i, j;

    /*    Directed, pentagram with ring, both clockwise    */
    igraph_matrix_init(&W, 1, 1);
    igraph_matrix_set(&W, 0, 0, 2);
    IGRAPH_ASSERT(igraph_extended_chordal_ring(&g, /* nodes */ 5, &W, 1 /*directed*/) == IGRAPH_SUCCESS);
    igraph_small(&g_test, 5, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 2, 1, 3, 2, 4, 3, 0, 4, 1, -1);
    IGRAPH_ASSERT(igraph_is_same_graph(&g, &g_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);

    /*     Use negative matrix value for same specification    */
    igraph_matrix_set(&W, 0, 0, -3);
    IGRAPH_ASSERT(igraph_extended_chordal_ring(&g_rev, /* nodes */ 5, &W, 1 /*directed*/) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_is_same_graph(&g_rev, &g_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&g);
    igraph_destroy(&g_rev);
    igraph_destroy(&g_test);
    igraph_matrix_destroy(&W);


    /*    From article, should give double edges for chords in igraph   */
    igraph_matrix_init(&W, 2, 2);
    int m[2][2] = {{4, 2}, 
                   {8, 10}};
    for (i=0; i < 2; i++) {
        for (j=0; j < 2; j++) {
            MATRIX(W, i, j) = m[i][j];
        }
    }
    IGRAPH_ASSERT(igraph_extended_chordal_ring(&g, /* nodes */ 12, &W, 0 /*undirected*/) == IGRAPH_SUCCESS);
    igraph_small(&g_test, 12, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 4,
                 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 0,
                 0, 4, 2, 6,  4, 8, 6, 10, 8, 0, 10, 2,
                 0, 4, 2, 6,  4, 8, 6, 10, 8, 0, 10, 2,
                 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 1,
                 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 1,
                 -1);

    IGRAPH_ASSERT(igraph_is_same_graph(&g, &g_test, &same) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(same);
    igraph_destroy(&g);
    igraph_destroy(&g_test);
    igraph_matrix_destroy(&W);

    VERIFY_FINALLY_STACK();
    return 0;
}
