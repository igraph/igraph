/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

#include "igraph_path_cover.h"

#include "igraph_error.h"
#include "igraph_interface.h"

/**
 * \function igraph_minimum_chain_cover
 * \brief Minimum chain cover of a directed acyclic graph.
 *
 * A chain of a DAG is a sequence of vertices that is totally ordered by
 * reachability (each vertex must reach the next one, but not necessarily
 * via a direct edge, unlike a path). A chain cover is a set of chains that
 * together cover every vertex at least once (chains, like the paths of
 * \ref igraph_minimum_path_cover(), may share vertices). Since every path
 * is automatically a valid chain, a minimum path cover -- such as one
 * computed by \ref igraph_minimum_path_cover() -- is always also a minimum
 * chain cover of the same (minimum) cardinality.
 *
 * </para><para>
 * This baseline implementation simply reuses the supplied path cover as
 * the chain cover, and exists as a stable entry point for potential future
 * chain-cover-specific algorithms that could construct one more directly.
 *
 * \param graph The input graph. Currently unused by this baseline
 *        implementation, but kept for API stability and validation in
 *        future, non-trivial chain cover algorithms.
 * \param cover A path cover of \p graph, e.g. as computed by
 *        \ref igraph_minimum_path_cover().
 * \param chain_cover Pointer to an initialized list of integer vectors.
 *        The chains of the resulting chain cover are stored here on
 *        return.
 * \return Error code.
 *
 * Time complexity: O(|V|), where |V| is the total number of vertices
 * covered by \p cover.
 *
 * </para><para>
 * Ported from naive_chaincover_from_pathcover() in
 * https://github.com/algbio/PerformanceMPC/blob/main/src/mpc/cc.cpp. That a
 * path cover is automatically a chain cover of the same size (and hence
 * minimum) follows from the path/chain duality underlying Dilworth's
 * theorem: every path is totally ordered by reachability, i.e. is itself a
 * chain.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Dilworth RP: A Decomposition Theorem for Partially Ordered Sets.
 * Annals of Mathematics, Second Series, 51(1):161-166, 1950.
 *
 * \sa \ref igraph_minimum_path_cover().
 */
igraph_error_t igraph_minimum_chain_cover(
        const igraph_t *graph,
        const igraph_vector_int_list_t *cover,
        igraph_vector_int_list_t *chain_cover) {

    igraph_int_t no_of_paths = igraph_vector_int_list_size(cover);
    igraph_int_t i;

    IGRAPH_UNUSED(graph);

    igraph_vector_int_list_clear(chain_cover);
    IGRAPH_CHECK(igraph_vector_int_list_reserve(chain_cover, no_of_paths));
    for (i = 0; i < no_of_paths; i++) {
        IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(
                chain_cover, igraph_vector_int_list_get_ptr(cover, i)));
    }

    return IGRAPH_SUCCESS;
}
