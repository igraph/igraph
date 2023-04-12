/*
   IGraph library.
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

#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_vector.h"

/**
 * \function igraph_tree_from_parent_vector
 * \brief Constructs a tree or forest from a vector encoding the parent of each vertex.
 *
 * \experimental
 *
 * Rooted trees and forests are conveniently represented using a \p parents
 * vector where the ID of the parent of vertex \c v is stored in <code>parents[v]</code>.
 * This function serves to construct an igraph graph from a parent vector representation.
 * The result is guaranteed to be a forest or a tree. If the \p parents vector
 * is found to encode a cycle or a self-loop, an error is raised.
 *
 * </para><para>
 * Several igraph functions produce such vectors, such as graph traversal
 * functions (\ref igraph_bfs() and \ref igraph_dfs()), shortest path functions
 * that construct a shortest path tree, as well as some other specialized
 * functions like \ref igraph_dominator_tree() or \ref igraph_cohesive_blocks().
 * Vertices which do not have parents (i.e. roots) get a negative entry in the
 * \p parents vector.
 *
 * </para><para>
 * Use \ref igraph_bfs() or \ref igraph_dfs() to convert a forest into a parent
 * vector representation. For trees, i.e. forests with a single root, it is
 * more convenient to use \ref igraph_bfs_simple().
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param parents The parent vector. <code>parents[v]</code> is the ID of
 *    the parent vertex of \c v. <code>parents[v] < 0</code> indicates that
 *    \c v does not have a parent.
 * \param type Constant, gives whether to create a directed tree, and
 *        if this is the case, also its orientation. Possible values:
 *        \clist
 *        \cli IGRAPH_TREE_OUT
 *          directed tree, the edges point from the parents to their children.
 *        \cli IGRAPH_TREE_IN
 *          directed tree, the edges point from the children to their parents.
 *        \cli IGRAPH_TREE_UNDIRECTED undirected tree.
 *        \endclist
 * \return Error code.
 *
 * \sa \ref igraph_bfs(), \ref igraph_bfs_simple() for back-conversion;
 * \ref igraph_from_prufer() for creating trees from Prüfer sequences;
 * \ref igraph_is_tree() and \ref igraph_is_forest() to check if a graph
 * is a tree or forest.
 *
 * Time complexity: O(n) where n is the length of \p parents.
 */
igraph_error_t igraph_tree_from_parent_vector(
        igraph_t *graph,
        const igraph_vector_int_t *parents,
        igraph_tree_mode_t type) {

    const igraph_integer_t no_of_nodes = igraph_vector_int_size(parents);
    igraph_vector_int_t seen;
    igraph_vector_int_t edges;
    igraph_bool_t directed, intree;

    switch (type) {
    case IGRAPH_TREE_OUT:
        directed = true; intree = false; break;
    case IGRAPH_TREE_IN:
        directed = true; intree = true; break;
    case IGRAPH_TREE_UNDIRECTED:
        directed = false; intree = true; break;
    default:
        IGRAPH_ERROR("Invalid tree mode.", IGRAPH_EINVAL);
    }

    /* Catch null graph case */
    if (no_of_nodes == 0) {
        return igraph_empty(graph, 0, directed);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&seen, no_of_nodes);

    /* A tree has no_of_nodes - 1 edges but a forest has fewer. In order to support
     * the use case of extracting small sub-trees of large graphs, we only reserve
     * the full amount of memory needed for a tree when the graph is small.
     * This also eliminates the need to check for integer overflow. */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_nodes > 1024 ? 2048 : 2*(no_of_nodes-1));
    igraph_vector_int_clear(&edges);

    igraph_integer_t c=1;
    for (igraph_integer_t i=0; i < no_of_nodes; i++) {
        igraph_integer_t v = i;

        if (VECTOR(seen)[v]) continue;

        while (true) {
            igraph_integer_t u;

            VECTOR(seen)[v] = c; /* mark v as seen in the current round */
            u = VECTOR(*parents)[v];

            if (u < 0) {
                break; /* v is a root, stop traversal */
            }
            if (u >= no_of_nodes) {
                IGRAPH_ERROR("Invalid vertex ID in parent vector.", IGRAPH_EINVVID);
            }

            if (intree) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, v));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, u));
            } else {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, u));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, v));
            }

            if (VECTOR(seen)[u]) {
                if (VECTOR(seen)[u] == c) {
                    /* u was already seen in the current round, we found a cycle.
                     * We distinguish between self-loops, i.e. 1-cycles, and longer
                     * cycles in order to make the error message more useful. */
                    IGRAPH_ERROR(
                            u==v
                              ? "Found a self-loop while constructing tree from parent vector."
                              : "Found a cycle while constructing tree from parent vector.",
                            IGRAPH_EINVAL);
                }
                break; /* u was seen in a previous round, stop traversal */
            }

            v = u;
        }

        c++;
    }

    igraph_vector_int_destroy(&seen);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
