/*
   igraph library.
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

void run_test(int n) {
    igraph_t graph;
    igraph_matrix_t layout;

    printf("\n%d-vertex graph:\n", n);
    igraph_empty(&graph, n, IGRAPH_UNDIRECTED);
    igraph_matrix_init(&layout, 0, 0);
    igraph_layout_sphere(&graph, &layout);
    IGRAPH_ASSERT(igraph_matrix_nrow(&layout) == igraph_vcount(&graph));
    IGRAPH_ASSERT(igraph_matrix_ncol(&layout) == 3);
    igraph_matrix_zapsmall(&layout, 0 /* automatic */);
    igraph_matrix_print(&layout);
    igraph_matrix_destroy(&layout);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
}

int main(void) {
    run_test(0);
    run_test(1);
    run_test(2);
    run_test(11);

    return 0;
}
