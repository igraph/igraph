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
    igraph_t graph;
    igraph_bool_t is_tree = 0, are_connected = 0;

    igraph_rng_seed(igraph_rng_default(), 74088);

    /* Undirected */

    IGRAPH_ASSERT(igraph_tree_game(&graph, 123, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_LERW) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_is_tree(&graph, &is_tree, NULL, IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(is_tree);
    igraph_destroy(&graph);

    IGRAPH_ASSERT(igraph_tree_game(&graph, 123, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_PRUFER) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_is_tree(&graph, &is_tree, NULL, IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(is_tree);
    igraph_destroy(&graph);

    /* Directed out-tree */

    IGRAPH_ASSERT(igraph_tree_game(&graph, 123, IGRAPH_DIRECTED, IGRAPH_RANDOM_TREE_LERW) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_is_tree(&graph, &is_tree, NULL, IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(is_tree);
    igraph_destroy(&graph);

    /* IGRAPH_RANDOM_TREE_PRUFER does not currently support directed graphs */

    /* Null graph */

    IGRAPH_ASSERT(igraph_tree_game(&graph, 0, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_LERW) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&graph) == 0);
    igraph_destroy(&graph);

    IGRAPH_ASSERT(igraph_tree_game(&graph, 0, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_PRUFER) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&graph) == 0);
    igraph_destroy(&graph);

    /* Singleton graph */

    IGRAPH_ASSERT(igraph_tree_game(&graph, 1, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_LERW) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&graph) == 1);
    igraph_destroy(&graph);

    IGRAPH_ASSERT(igraph_tree_game(&graph, 1, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_PRUFER) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&graph) == 1);
    igraph_destroy(&graph);

    /* P_2 */

    IGRAPH_ASSERT(igraph_tree_game(&graph, 2, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_LERW) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&graph) == 2);
    IGRAPH_ASSERT(igraph_ecount(&graph) == 1);
    IGRAPH_ASSERT(igraph_are_connected(&graph, 0, 1, &are_connected) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(are_connected);
    igraph_destroy(&graph);

    IGRAPH_ASSERT(igraph_tree_game(&graph, 2, IGRAPH_UNDIRECTED, IGRAPH_RANDOM_TREE_PRUFER) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&graph) == 2);
    IGRAPH_ASSERT(igraph_ecount(&graph) == 1);
    IGRAPH_ASSERT(igraph_are_connected(&graph, 0, 1, &are_connected) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(are_connected);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
    return 0;
}
