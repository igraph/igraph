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


#define PRINT_DESTROY(name) \
    printf(name "\n"); \
    print_graph_canon(&graph); \
    igraph_destroy(&graph); \
    printf("\n");


int main(void) {
    igraph_t graph;

    igraph_kary_tree(&graph, 0, 1, IGRAPH_TREE_UNDIRECTED);
    PRINT_DESTROY("Null graph");

    igraph_kary_tree(&graph, 0, 1, IGRAPH_TREE_OUT);
    PRINT_DESTROY("Directed null graph");

    igraph_kary_tree(&graph, 1, 1, IGRAPH_TREE_UNDIRECTED);
    PRINT_DESTROY("Singleton graph");

    igraph_kary_tree(&graph, 3, 1, IGRAPH_TREE_OUT);
    PRINT_DESTROY("Path graph");

    igraph_kary_tree(&graph, 3, 2, IGRAPH_TREE_OUT);
    PRINT_DESTROY("Binary out-tree, n=3");

    igraph_kary_tree(&graph, 3, 2, IGRAPH_TREE_IN);
    PRINT_DESTROY("Binary in-tree, n=3");

    igraph_kary_tree(&graph, 14, 3, IGRAPH_TREE_OUT);
    PRINT_DESTROY("Ternary out-tree, n=14");

    VERIFY_FINALLY_STACK();

    return 0;
}
