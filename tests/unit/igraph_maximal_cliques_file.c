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
    igraph_t g_empty, g_lm;

    igraph_small(&g_empty, 0, 0, -1);
    igraph_small(&g_lm, 6, 1, 0,1, 0,2, 1,1, 1,2, 1,3, 2,0, 2,3, 3,4, 3,4, -1);

    printf("No vertices:\n");
    IGRAPH_ASSERT(igraph_maximal_cliques_file(&g_empty, stdout, /*min*/ 0, /*max*/ 0) == IGRAPH_SUCCESS);

    printf("\nGraph with loops and multiple edges:\n");
    IGRAPH_ASSERT(igraph_maximal_cliques_file(&g_lm, stdout, /*min*/ 0, /*max*/ 0) == IGRAPH_SUCCESS);

    printf("\nGraph with loops and multiple edges, minimum clique size 10:\n");
    IGRAPH_ASSERT(igraph_maximal_cliques_file(&g_lm, stdout, /*min*/ 10, /*max*/ 0) == IGRAPH_SUCCESS);

    printf("\nGraph with loops and multiple edges, maximum clique size 2:\n");
    IGRAPH_ASSERT(igraph_maximal_cliques_file(&g_lm, stdout, /*min*/ 0, /*max*/ 2) == IGRAPH_SUCCESS);

    igraph_destroy(&g_empty);
    igraph_destroy(&g_lm);

    VERIFY_FINALLY_STACK();
    return 0;
}
