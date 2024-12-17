/*
  Copyright 2017 The Johns Hopkins University Applied Physics Laboratory LLC. All Rights Reserved.
  Copyright 2021 The igraph team.

  Truss algorithm for cohesive subgroups.

  Author: Alex Perrone
  Date: 2017-08-03
  Minor edits: The igraph team, 2021

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
  Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301 USA

*/

#include <vector>
#include <unordered_set>

#include "igraph_community.h"

#include "igraph_adjlist.h"
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_motifs.h"
#include "igraph_structural.h"

#include "core/exceptions.h"
#include "core/interruption.h"

using std::vector;
using std::unordered_set;


// Unpack the triangles as a vector of vertices to be a vector of edges.
// So, instead of the triangle specified as vertices [1, 2, 3], return the
// edges as [1, 2, 1, 3, 2, 3] so that the support can be computed.
static igraph_error_t igraph_truss_i_unpack(const igraph_vector_int_t *tri, igraph_vector_int_t *unpacked_tri) {
    igraph_integer_t num_triangles = igraph_vector_int_size(tri);

    IGRAPH_CHECK(igraph_vector_int_resize(unpacked_tri, 2 * num_triangles));

    for (igraph_integer_t i = 0, j = 0; i < num_triangles; i += 3, j += 6) {
        VECTOR(*unpacked_tri)[j]   = VECTOR(*unpacked_tri)[j+2] = VECTOR(*tri)[i];
        VECTOR(*unpacked_tri)[j+1] = VECTOR(*unpacked_tri)[j+4] = VECTOR(*tri)[i+1];
        VECTOR(*unpacked_tri)[j+3] = VECTOR(*unpacked_tri)[j+5] = VECTOR(*tri)[i+2];
    }

    return IGRAPH_SUCCESS;
}


// Compute the edge support, i.e. number of triangles each edge occurs in.
// Time complexity: O(m), where m is the number of edges listed in eid.
static void igraph_truss_i_compute_support(const igraph_vector_int_t *eid, igraph_vector_int_t *support) {
    igraph_integer_t m = igraph_vector_int_size(eid);
    for (igraph_integer_t i = 0; i < m; ++i) {
        VECTOR(*support)[VECTOR(*eid)[i]] += 1;
    }
}


/* internal function doing the computations once the support is defined */
static igraph_error_t igraph_i_trussness(const igraph_t *graph, igraph_vector_int_t *support,
                                         igraph_vector_int_t *trussness) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    igraph_adjlist_t adjlist;
    igraph_vector_int_t commonNeighbors;
    igraph_vector_bool_t completed;

    // C++ data structures
    vector< unordered_set<igraph_integer_t> > vec;

    // Allocate memory for result
    igraph_integer_t no_of_edges = igraph_vector_int_size(support);
    IGRAPH_CHECK(igraph_vector_int_resize(trussness, no_of_edges));
    if (no_of_edges == 0) {
        return IGRAPH_SUCCESS;
    }

    // Get max possible value = max entry in support.
    // This cannot be computed if there are no edges, hence the above check
    igraph_integer_t max = igraph_vector_int_max(support);

    // Initialize completed edges.
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&completed, no_of_edges);

    // The vector of levels. Each level of the vector is a set of edges initially
    // at that level of support, where support is # of triangles the edge is in.
    vec.resize(max + 1);

    // Add each edge to its appropriate level of support.
    for (igraph_integer_t i = 0; i < no_of_edges; ++i) {
        vec[VECTOR(*support)[i]].insert(i);  // insert edge i into its support level
    }

    // Record the trussness of edges at level 0. These edges are not part
    // of any triangles, so there's not much to do and we "complete" them
    for (auto edge : vec[0]) {
        VECTOR(*trussness)[edge] = 2;
        VECTOR(completed)[edge] = true;
    }

    // Initialize variables needed below.
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&commonNeighbors, 0);

    // Move through the levels, one level at a time, starting at first level.
    for (igraph_integer_t level = 1; level <= max; ++level) {

        /* Track down edges one at a time */
        while (!vec[level].empty()) {
            IGRAPH_ALLOW_INTERRUPTION();

            igraph_integer_t seed = *vec[level].begin();  // pull out the first edge
            vec[level].erase(seed);  // remove the first element

            /* Find the vertices of this edge */
            igraph_integer_t fromVertex = IGRAPH_FROM(graph, seed);
            igraph_integer_t toVertex = IGRAPH_TO(graph, seed);

            /* Find neighbors of both vertices. If they run into each other,
             * there is a triangle. We rely on the neighbor lists being sorted,
             * as guaranteed by igraph_adjlist_init(), when computing intersections. */
            igraph_vector_int_t *fromNeighbors = igraph_adjlist_get(&adjlist, fromVertex);
            igraph_vector_int_t *toNeighbors = igraph_adjlist_get(&adjlist, toVertex);
            igraph_vector_int_t *q1 = fromNeighbors;
            igraph_vector_int_t *q2 = toNeighbors;

            if (igraph_vector_int_size(q1) > igraph_vector_int_size(q2)) {
                // case: #fromNeighbors > #toNeigbors, so make q1 the smaller set.
                q1 = toNeighbors;
                q2 = fromNeighbors;
            }

            // Intersect the neighbors.
            IGRAPH_CHECK(igraph_vector_int_intersect_sorted(q1, q2, &commonNeighbors));

            /* Go over the overlapping neighbors and check each */
            igraph_integer_t ncommon = igraph_vector_int_size(&commonNeighbors);
            for (igraph_integer_t j = 0; j < ncommon; j++) {
                igraph_integer_t n = VECTOR(commonNeighbors)[j];  // the common neighbor
                igraph_integer_t e1, e2;
                IGRAPH_CHECK(igraph_get_eid(graph, &e1, fromVertex, n, IGRAPH_UNDIRECTED, /* error= */ true));
                IGRAPH_CHECK(igraph_get_eid(graph, &e2, toVertex, n, IGRAPH_UNDIRECTED, /* error= */ true));

                bool e1_complete = VECTOR(completed)[e1];
                bool e2_complete = VECTOR(completed)[e2];

                if (!e1_complete && !e2_complete) {
                    igraph_integer_t newLevel;

                    // Demote this edge, if higher than current level.
                    if (VECTOR(*support)[e1] > level) {
                        VECTOR(*support)[e1] -= 1;  // decrement the level
                        newLevel = VECTOR(*support)[e1];
                        vec[newLevel].insert(e1);
                        vec[newLevel + 1].erase(e1);  // the old level
                    }
                    // Demote this edge, if higher than current level.
                    if (VECTOR(*support)[e2] > level) {
                        VECTOR(*support)[e2] -= 1;  // decrement the level
                        newLevel = VECTOR(*support)[e2];
                        vec[newLevel].insert(e2);
                        vec[newLevel + 1].erase(e2);  // the old level
                    }
                }
            }
            // Record this edge; its level is its trussness.
            VECTOR(*trussness)[seed] = level + 2;
            VECTOR(completed)[seed] = true; // mark as complete
            igraph_vector_int_clear(&commonNeighbors);
        }  // end while
    }  // end for-loop over levels

    // Clean up.
    igraph_vector_int_destroy(&commonNeighbors);
    igraph_adjlist_destroy(&adjlist);
    igraph_vector_bool_destroy(&completed);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;

    IGRAPH_HANDLE_EXCEPTIONS_END;
}


/**
 * \function igraph_trussness
 * \brief Finding the "trussness" of the edges in a network.
 *
 * A k-truss is a subgraph in which every edge occurs in at least <code>k-2</code> triangles
 * in the subgraph. The trussness of an edge indicates the highest k-truss that
 * the edge occurs in.
 *
 * </para><para>
 * This function returns the highest \c k for each edge. If you are interested in
 * a particular k-truss subgraph, you can subset the graph to those edges
 * which are <code>&gt;= k</code> because each k-truss is a subgraph of a <code>(k–1)</code>-truss
 * Thus, to get all 4-trusses, take <code>k >= 4</code> because the 5-trusses, 6-trusses,
 * etc. need to be included.
 *
 * </para><para>
 * The current implementation of this function iteratively decrements support
 * of each edge using O(|E|) space and O(|E|^1.5) time. The implementation does
 * not support multigraphs; use \ref igraph_simplify() to collapse edges before
 * calling this function.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * See Algorithm 2 in:
 * Wang, Jia, and James Cheng. "Truss decomposition in massive networks."
 * Proceedings of the VLDB Endowment 5.9 (2012): 812-823.
 * https://doi.org/10.14778/2311906.2311909
 *
 * \param graph The input graph. Loop edges are allowed; multigraphs are not.
 * \param truss Pointer to initialized vector of truss values that will
 * indicate the highest k-truss each edge occurs in. It will be resized as
 * needed.
 * \return Error code.
 *
 * Time complexity: It should be O(|E|^1.5) according to the reference.
 */
igraph_error_t igraph_trussness(const igraph_t* graph, igraph_vector_int_t* trussness) {
    igraph_vector_int_t triangles, support, unpacked_triangles, eid;
    igraph_bool_t is_multigraph;

    /* Check whether the graph is a multigraph; trussness will not work for these */
    IGRAPH_CHECK(igraph_has_multiple(graph, &is_multigraph));
    if (! is_multigraph && igraph_is_directed(graph)) {
        /* Directed graphs with mutual edges are effectively multigraphs
         * when edge directions are ignored. */
        IGRAPH_CHECK(igraph_has_mutual(graph, &is_multigraph, /* loops */ false));
    }
    if (is_multigraph) {
        IGRAPH_ERROR("Trussness is not implemented for graphs with multi-edges.", IGRAPH_UNIMPLEMENTED);
    }

    /* Manage the stack to make it memory safe: do not change the order of
     * initialization of the following four vectors */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&support, igraph_ecount(graph));
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eid, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&unpacked_triangles, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&triangles, 0);

    // List the triangles as vertex triplets.
    IGRAPH_CHECK(igraph_list_triangles(graph, &triangles));

    // Unpack the triangles from vertex list to edge list.
    IGRAPH_CHECK(igraph_truss_i_unpack(&triangles, &unpacked_triangles));
    igraph_vector_int_destroy(&triangles);
    IGRAPH_FINALLY_CLEAN(1);

    // Get the edge IDs of the unpacked triangles. Note: a given eid can occur
    // multiple times in this list if it is in multiple triangles.
    IGRAPH_CHECK(igraph_get_eids(graph, &eid, &unpacked_triangles, /* directed = */ false, /* error = */ true));
    igraph_vector_int_destroy(&unpacked_triangles);
    IGRAPH_FINALLY_CLEAN(1);

    // Compute the support of the edges.
    igraph_truss_i_compute_support(&eid, &support);
    igraph_vector_int_destroy(&eid);
    IGRAPH_FINALLY_CLEAN(1);

    // Compute the trussness of the edges.
    IGRAPH_CHECK(igraph_i_trussness(graph, &support, trussness));
    igraph_vector_int_destroy(&support);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
