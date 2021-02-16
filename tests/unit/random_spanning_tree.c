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
    igraph_t graph, spanning_tree;
    igraph_vector_t tree_edges;
    igraph_bool_t is_tree;
    int err;

    igraph_rng_seed(igraph_rng_default(), 987);

    igraph_vector_init(&tree_edges, 0);

    /* This is guaranteed to create a connected graph. */
    igraph_barabasi_game(&graph, 100, 2, 2, NULL, 0, 1, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE, NULL);

    err = igraph_random_spanning_tree(&graph, &tree_edges, 0);
    IGRAPH_ASSERT(!err);

    IGRAPH_ASSERT(igraph_vector_size(&tree_edges) == igraph_vcount(&graph) - 1);

    err = igraph_subgraph_edges(&graph, &spanning_tree, igraph_ess_vector(&tree_edges), /* delete_vertices= */ 0);
    IGRAPH_ASSERT(!err);

    IGRAPH_ASSERT(igraph_vcount(&spanning_tree) == igraph_vcount(&graph));

    igraph_is_tree(&spanning_tree, &is_tree, NULL, IGRAPH_ALL);
    IGRAPH_ASSERT(is_tree);

    igraph_destroy(&spanning_tree);
    igraph_destroy(&graph);

    /* Non-connected forest graph. There is only one solution. */
    igraph_small(&graph, 4, IGRAPH_UNDIRECTED, 0,1, 2,3, -1);

    /* Find a spanning tree of the component containing vertex 0 */
    err = igraph_random_spanning_tree(&graph, &tree_edges, 0);
    IGRAPH_ASSERT(!err);

    IGRAPH_ASSERT(igraph_vector_size(&tree_edges) == 1);
    IGRAPH_ASSERT(VECTOR(tree_edges)[0] == 0);

    /* Find a spanning forest */
    err = igraph_random_spanning_tree(&graph, &tree_edges, -1);
    IGRAPH_ASSERT(!err);

    IGRAPH_ASSERT(igraph_vector_size(&tree_edges) == 2);
    IGRAPH_ASSERT(VECTOR(tree_edges)[0] == 0 && VECTOR(tree_edges)[1] == 1);

    igraph_destroy(&graph);

    igraph_vector_destroy(&tree_edges);

    VERIFY_FINALLY_STACK();
    return 0;
}
