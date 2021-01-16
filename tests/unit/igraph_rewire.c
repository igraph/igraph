/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2021  The igraph development team <igraph@igraph.org>

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

#include "operators/rewire_internal.h"

#include "test_utilities.inc"

static void check_rewiring(igraph_tree_mode_t tree_mode, igraph_bool_t use_adjlist, igraph_bool_t allow_loops, const char* description) {

    igraph_t g;
    igraph_vector_t indegree_before, outdegree_before, indegree_after, outdegree_after;

    igraph_tree(&g, 10, 3, tree_mode);

    igraph_vector_init(&indegree_before, 0);
    igraph_vector_init(&outdegree_before, 0);
    igraph_degree(&g, &indegree_before, igraph_vss_all(), IGRAPH_IN, 1);
    igraph_degree(&g, &outdegree_before, igraph_vss_all(), IGRAPH_OUT, 1);

    igraph_i_rewire(&g, 1000, allow_loops ? IGRAPH_REWIRING_SIMPLE_LOOPS : IGRAPH_REWIRING_SIMPLE, use_adjlist);

    igraph_vector_init(&indegree_after, 0);
    igraph_vector_init(&outdegree_after, 0);
    igraph_degree(&g, &indegree_after, igraph_vss_all(), IGRAPH_IN, 1);
    igraph_degree(&g, &outdegree_after, igraph_vss_all(), IGRAPH_OUT, 1);

    if ((!igraph_vector_all_e(&indegree_before, &indegree_after)) ||
        (!igraph_vector_all_e(&outdegree_before, &outdegree_after))) {

        printf("%s: graph degrees changed. Rewired graph is below.\n", description);
        print_graph(&g);

        abort();
    }

    igraph_destroy(&g);
    igraph_vector_destroy(&indegree_before);
    igraph_vector_destroy(&outdegree_before);
    igraph_vector_destroy(&indegree_after);
    igraph_vector_destroy(&outdegree_after);

}

int main() {
    igraph_rng_seed(igraph_rng_default(), 3925);

    /* Short test for the top-level igraph_rewire() functions (instead of igraph_i_rewire()). */
    {
        igraph_t graph;
        igraph_ring(&graph, 12, IGRAPH_UNDIRECTED, /* mutual= */ 0, /* circular= */ 1);
        igraph_rewire(&graph, 50, IGRAPH_REWIRING_SIMPLE);
        igraph_destroy(&graph);
    }

    check_rewiring(IGRAPH_TREE_OUT, 0, 0, "Directed, no loops, standard-method");
    check_rewiring(IGRAPH_TREE_OUT, 1, 0, "Directed, no loops, adjlist-method");
    check_rewiring(IGRAPH_TREE_OUT, 0, 1, "Directed, loops, standard-method");
    check_rewiring(IGRAPH_TREE_OUT, 1, 1, "Directed, loops, adjlist-method");
    check_rewiring(IGRAPH_TREE_UNDIRECTED, 0, 0, "Undirected, no loops, standard-method");
    check_rewiring(IGRAPH_TREE_UNDIRECTED, 1, 0, "Undirected, no loops, adjlist-method");
    check_rewiring(IGRAPH_TREE_UNDIRECTED, 0, 1, "Undirected, loops, standard-method");
    check_rewiring(IGRAPH_TREE_UNDIRECTED, 1, 1, "Undirected, loops, adjlist-method");

    VERIFY_FINALLY_STACK();
    return 0;
}
