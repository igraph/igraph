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

void print_and_destroy(igraph_t *g) {
    igraph_es_t eids;
    igraph_vector_int_t edges;

    igraph_es_all(&eids, IGRAPH_EDGEORDER_ID);
    igraph_vector_int_init(&edges, 0);

    printf("  row-wise order: ");
    igraph_edges(g, eids, &edges, 0);
    print_vector_int(&edges);

    printf("  column-wise order: ");
    igraph_edges(g, eids, &edges, 1);
    print_vector_int(&edges);

    igraph_destroy(g);
    igraph_vector_int_destroy(&edges);
}

int main(void) {
    igraph_t g;

    printf("No edges:\n");
    igraph_small(&g, 0, 1, -1);
    print_and_destroy(&g);

    printf("Directed graph with loops and multiple edges:\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED,
            0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 4,3, 4,3, -1);
    print_and_destroy(&g);

    printf("Undirected graph with loops and multiple edges:\n");
    igraph_small(&g, 6, IGRAPH_UNDIRECTED,
            0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 4,3, 4,3, -1);
    print_and_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
