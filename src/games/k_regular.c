/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2021 The igraph development team

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

#include "igraph_games.h"

#include "igraph_interface.h"

/**
 * \ingroup generators
 * \function igraph_k_regular_game
 * \brief Generates a random graph where each vertex has the same degree.
 *
 * This game generates a directed or undirected random graph where the
 * degrees of vertices are equal to a predefined constant k. For undirected
 * graphs, at least one of k and the number of vertices must be even.
 *
 * </para><para>
 * Currently, this game simply uses \ref igraph_degree_sequence_game with
 * the \c SIMPLE_NO_MULTIPLE method and appropriately constructed degree sequences.
 * Thefore, it does not sample uniformly: while it can generate all k-regular graphs
 * with the given number of vertices, it does not generate each one with the same
 * probability.
 *
 * \param graph        Pointer to an uninitialized graph object.
 * \param no_of_nodes  The number of nodes in the generated graph.
 * \param k            The degree of each vertex in an undirected graph, or
 *                     the out-degree and in-degree of each vertex in a
 *                     directed graph.
 * \param directed     Whether the generated graph will be directed.
 * \param multiple     Whether to allow multiple edges in the generated graph.
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid parameter; e.g., negative number of nodes,
 *                           or odd number of nodes and odd k for undirected
 *                           graphs.
 *         \c IGRAPH_ENOMEM: there is not enough memory for the operation.
 *
 * Time complexity: O(|V|+|E|) if \c multiple is true, otherwise not known.
 */
int igraph_k_regular_game(igraph_t *graph,
                          igraph_integer_t no_of_nodes, igraph_integer_t k,
                          igraph_bool_t directed, igraph_bool_t multiple) {
    igraph_vector_t degseq;
    igraph_degseq_t mode = multiple ? IGRAPH_DEGSEQ_SIMPLE : IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE;

    /* Note to self: we are not using IGRAPH_DEGSEQ_VL when multiple = false
     * because the VL method is not really good at generating k-regular graphs.
     * Actually, that's why we have added SIMPLE_NO_MULTIPLE. */

    if (no_of_nodes < 0) {
        IGRAPH_ERROR("number of nodes must be non-negative", IGRAPH_EINVAL);
    }
    if (k < 0) {
        IGRAPH_ERROR("degree must be non-negative", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&degseq, no_of_nodes);
    igraph_vector_fill(&degseq, k);
    IGRAPH_CHECK(igraph_degree_sequence_game(graph, &degseq, directed ? &degseq : 0, mode));

    igraph_vector_destroy(&degseq);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
