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
#include "test_utilities.h"

#define PRINT_DESTROY(name) \
    printf(name "\n"); \
    print_graph_canon(&g); \
    igraph_destroy(&g); \
    igraph_vector_int_destroy(&v); \
    printf("\n");

int main(void) {

    igraph_t g;
    igraph_vector_int_t v;

    /**** Test with empty vector ****/
    // initialize variables for test
    igraph_vector_int_init_int(&v, 0);
    // test
    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_UNDIRECTED) == IGRAPH_SUCCESS);
    // root vertex always gets created, cannot create null graph with this function
    PRINT_DESTROY("Singleton graph");

    igraph_vector_int_init_int(&v, 0);
    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_OUT) == IGRAPH_SUCCESS);
    PRINT_DESTROY("Directed singleton graph");

    /**** Test with -1 value ****/
    // invalid number of children
    igraph_vector_int_init_int(&v, 1, -1);
    CHECK_ERROR(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_UNDIRECTED), IGRAPH_EINVAL);
    igraph_destroy(&g);
    igraph_vector_int_destroy(&v);

    /**** Test undirected symmetric graph with 1 child ****/
    // 1 edge, 2 vertices (root and child)
    igraph_vector_int_init_int(&v, 1, 1);
    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_UNDIRECTED) ==  IGRAPH_SUCCESS);
    PRINT_DESTROY("Undirected graph with 2 vertices and 1 edge");

    /**** Test directed symmetric graph with 1 child ****/
    // 1 edge, 2 vertices (root and child)
    igraph_vector_int_init_int(&v, 1, 1);
    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_IN) ==  IGRAPH_SUCCESS);
    PRINT_DESTROY("Directed graph with 2 vertices and 1 edge");

    /**** Test directed symmetric path graph with 1 child in each level ****/
    igraph_vector_int_init_int(&v, 3, 1, 1, 1);
    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_OUT) == IGRAPH_SUCCESS);
    PRINT_DESTROY("Directed path graph with 4 level and 1 child in each");

    /**** Test undirected symmetric path graph with 1 child in each level ****/
    igraph_vector_int_init_int(&v, 3, 1, 1, 1);
    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_UNDIRECTED) == IGRAPH_SUCCESS);
    PRINT_DESTROY("Undirected path graph with 4 level and 1 child in each");

    /**** Test directed symmetric graph as binary tree with 1 level ****/
    igraph_vector_int_init_int(&v, 1, 2);
    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_OUT) == IGRAPH_SUCCESS);
    PRINT_DESTROY("Binary out-tree with 3 vertices");

    /**** Test undirected symmetric graph as binary tree with 1 level ****/
    igraph_vector_int_init_int(&v, 1, 2);
    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_UNDIRECTED) == IGRAPH_SUCCESS);
    PRINT_DESTROY("Undirected binery tree with 3 vertices");


    /**** Test directed symmetric graph with 2 level, each 3 children  ****/
    igraph_vector_int_init_int(&v, 2, 3, 3);
    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_OUT) == IGRAPH_SUCCESS);
    PRINT_DESTROY("Symmetric out-tree with 2 level, 3 children in each");

    /**** Test undirected symmetric graph with 2 level, each 3 children  ****/
    igraph_vector_int_init_int(&v, 2, 3, 3);
    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_UNDIRECTED) == IGRAPH_SUCCESS);
    PRINT_DESTROY("Symmetric undirected tree with 2 level, 3 children in each");

    /**** Test directed symmetric graph with 2 level, first with 3 and second with 4 children  ****/
    igraph_vector_int_init_int(&v, 2, 3, 4);
    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_OUT) == IGRAPH_SUCCESS);
    PRINT_DESTROY("Symmetric out-tree with 2 level, 3 children in first, 4 in second");

    /**** Test undirected symmetric graph with 2 level, first with 4 and second with 3 children  ****/
    igraph_vector_int_init_int(&v, 2, 4, 3);
    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_UNDIRECTED) == IGRAPH_SUCCESS);
    PRINT_DESTROY("Undirected symmetric tree with 2 level, 4 children in first, 3 in second");

    /**** Test undirected symmetric graph with 3 level, first with 3, second with 4, third with 5 children  ****/
    igraph_vector_int_init_int(&v, 3, 3, 4, 5);
    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_UNDIRECTED) == IGRAPH_SUCCESS);
    PRINT_DESTROY("Undirected symmetric tree with 3 level, 3 children in first, 4 in second, 5 in third");

    VERIFY_FINALLY_STACK();
    return 0;
}
