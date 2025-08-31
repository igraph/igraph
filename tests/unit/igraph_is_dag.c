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

int main(void) {
    igraph_t graph;
    igraph_bool_t is_dag;

    // directed empty graphs
    for (igraph_int_t n=0; n < 3; n++ ){
        igraph_empty(&graph, n, IGRAPH_DIRECTED);
        igraph_is_dag(&graph, &is_dag);
        IGRAPH_ASSERT(is_dag);
        igraph_destroy(&graph);
    }

    // undirected graphs
    for (igraph_int_t n=0; n < 3; n++ ){
        igraph_empty(&graph, n, IGRAPH_UNDIRECTED);
        igraph_is_dag(&graph, &is_dag);
        IGRAPH_ASSERT(!is_dag);
        igraph_destroy(&graph);
    }

    // 1-path
    igraph_small(&graph, 2, IGRAPH_DIRECTED,
                 0,1,
                 -1);
    igraph_is_dag(&graph, &is_dag);
    IGRAPH_ASSERT(is_dag);
    igraph_destroy(&graph);

    // reciprocal edges -- not a DAG
    igraph_small(&graph, 2, IGRAPH_DIRECTED,
                 0,1, 1,0,
                 -1);
    igraph_is_dag(&graph, &is_dag);
    IGRAPH_ASSERT(! is_dag);
    igraph_destroy(&graph);

    // 4-cycle -- not a DAG
    igraph_small(&graph, 4, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,3, 3,0,
                 -1);
    igraph_is_dag(&graph, &is_dag);
    IGRAPH_ASSERT(!is_dag);
    igraph_destroy(&graph);

    // 4-cycle with outgoing edge -- not a DAG
    igraph_small(&graph, 5, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,3, 3,0, 0,4,
                 -1);
    igraph_is_dag(&graph, &is_dag);
    IGRAPH_ASSERT(!is_dag);
    igraph_destroy(&graph);

    // 4-cycle with incoming edge -- not a DAG
    igraph_small(&graph, 5, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,3, 3,0, 4,0,
                 -1);
    igraph_is_dag(&graph, &is_dag);
    IGRAPH_ASSERT(!is_dag);
    igraph_destroy(&graph);

    // X shape, DAG
    igraph_small(&graph, 5, IGRAPH_DIRECTED,
                 1,0, 2,0, 0,3, 0,4,
                 -1);
    igraph_is_dag(&graph, &is_dag);
    IGRAPH_ASSERT(is_dag);
    igraph_destroy(&graph);

    // self-loop present, not a DAG
    igraph_small(&graph, 5, IGRAPH_DIRECTED,
                 1,0, 2,0, 0,3, 0,4, 0,0,
                 -1);
    igraph_is_dag(&graph, &is_dag);
    IGRAPH_ASSERT(!is_dag);
    igraph_destroy(&graph);

    // singleton with self-loop
    igraph_small(&graph, 1, IGRAPH_DIRECTED,
                 0,0,
                 -1);
    igraph_is_dag(&graph, &is_dag);
    IGRAPH_ASSERT(!is_dag);
    igraph_destroy(&graph);

    // out-tree
    igraph_kary_tree(&graph, 6, 2, IGRAPH_TREE_OUT);
    igraph_is_dag(&graph, &is_dag);
    IGRAPH_ASSERT(is_dag);
    igraph_destroy(&graph);

    // in-tree
    igraph_kary_tree(&graph, 6, 2, IGRAPH_TREE_IN);
    igraph_is_dag(&graph, &is_dag);
    IGRAPH_ASSERT(is_dag);
    igraph_destroy(&graph);

    // undirected -- not a DAG
    igraph_kary_tree(&graph, 6, 2, IGRAPH_TREE_UNDIRECTED);
    igraph_is_dag(&graph, &is_dag);
    IGRAPH_ASSERT(!is_dag);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
    return 0;
}
