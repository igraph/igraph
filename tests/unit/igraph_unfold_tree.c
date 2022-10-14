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

void print_and_destroy(igraph_t *g, igraph_vector_int_t *roots) {
    igraph_t tree;
    igraph_vector_int_t vertex_index;

    igraph_vector_int_init(&vertex_index, 0);

    igraph_unfold_tree(g, &tree, IGRAPH_OUT, roots, &vertex_index);

    printf("Tree:\n");
    print_graph_canon(&tree);
    printf("Vertex indices:\n");
    igraph_vector_int_print(&vertex_index);
    printf("\n");

    igraph_destroy(&tree);
    igraph_destroy(g);
    igraph_vector_int_destroy(roots);
    igraph_vector_int_destroy(&vertex_index);
}

int main(void) {
    igraph_t g, tree;
    igraph_vector_int_t roots;

    printf("Graph with no vertices\n");
    igraph_vector_int_init(&roots, 0);
    igraph_small(&g, 0, IGRAPH_UNDIRECTED, -1);
    print_and_destroy(&g, &roots);

    printf("Graph with loops, multiple edges and isolated vertex:\n");
    igraph_vector_int_init_int(&roots, 1, 0);
    igraph_small(&g, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
    print_and_destroy(&g, &roots);

    printf("Same graph, undirected:\n");
    igraph_vector_int_init_int(&roots, 1, 0);
    igraph_small(&g, 6, IGRAPH_UNDIRECTED, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
    print_and_destroy(&g, &roots);

    printf("Almost same graph, two roots:\n");
    igraph_vector_int_init_int(&roots, 2, 0, 5);
    igraph_small(&g, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, 5,5, 5,5, -1);
    print_and_destroy(&g, &roots);

    printf("Same graph, multiple roots in same tree:\n");
    igraph_vector_int_init_int(&roots, 5, 0, 0, 1, 2, 3);
    igraph_small(&g, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, 5,5, 5,5, -1);
    print_and_destroy(&g, &roots);

    VERIFY_FINALLY_STACK();

    printf("Check error for root out of bounds.\n");
    igraph_vector_int_init_int(&roots, 5, -1, 0, 1, 2, 3);
    igraph_small(&g, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, 5,5, 5,5, -1);
    CHECK_ERROR(igraph_unfold_tree(&g, &tree, IGRAPH_OUT, &roots, NULL), IGRAPH_EINVVID);
    igraph_destroy(&g);
    igraph_vector_int_destroy(&roots);

    VERIFY_FINALLY_STACK();
    return 0;
}
