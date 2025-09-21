/*
   igraph library.
   Copyright (C) 2022 The igraph development team

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

#include "test_utilities.h"

int main(void) {
    igraph_t graph;

    igraph_small(&graph, 0, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,3, 3,1, 1,4,
                 -1);

    printf("Original graph:\n");
    print_graph(&graph);

    printf("Reverse one edge:\n");
    igraph_reverse_edges(&graph, igraph_ess_1(2));
    print_graph(&graph);

    printf("Reverse all edges:\n");
    igraph_reverse_edges(&graph, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    print_graph(&graph);

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
