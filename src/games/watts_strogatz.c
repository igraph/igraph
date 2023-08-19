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

#include "igraph_constructors.h"
#include "igraph_interface.h"

/**
 * \function igraph_watts_strogatz_game
 * \brief The Watts-Strogatz small-world model.
 *
 * This function generates networks with the small-world property
 * based on a variant of the Watts-Strogatz model. The network is obtained
 * by first creating a periodic undirected lattice, then rewiring both
 * endpoints of each edge with probability \p p, while avoiding the
 * creation of multi-edges.
 *
 * </para><para>
 * This process differs from the original model of Watts and Strogatz
 * (see reference) in that it rewires \em both endpoints of edges. Thus in
 * the limit of <code>p=1</code>, we obtain a G(n,m) random graph with the
 * same number of vertices and edges as the original lattice. In comparison,
 * the original Watts-Strogatz model only rewires a single endpoint of each edge,
 * thus the network does not become fully random even for <code>p=1</code>.
 * For appropriate choices of \p p, both models exhibit the property of
 * simultaneously having short path lengths and high clustering.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Duncan J Watts and Steven H Strogatz:
 * Collective dynamics of <quote>small world</quote> networks, Nature
 * 393, 440-442, 1998.
 *
 * \param graph The graph to initialize.
 * \param dim The dimension of the lattice.
 * \param size The size of the lattice along each dimension.
 * \param nei The size of the neighborhood for each vertex. This is
 *    the same as the \p nei argument of \ref igraph_connect_neighborhood().
 * \param p The rewiring probability. A real number between zero and
 *   one (inclusive).
 * \param loops Logical, whether to generate loop edges.
 * \param multiple Logical, whether to allow multiple edges in the
 *   generated graph.
 * \return Error code.
 *
 * \sa \ref igraph_square_lattice(), \ref igraph_connect_neighborhood() and
 * \ref igraph_rewire_edges() can be used if more flexibility is
 * needed, e.g. a different type of lattice.
 *
 * Time complexity: O(|V|*d^o+|E|), |V| and |E| are the number of
 * vertices and edges, d is the average degree, o is the \p nei
 * argument.
 */
igraph_error_t igraph_watts_strogatz_game(igraph_t *graph, igraph_integer_t dim,
                               igraph_integer_t size, igraph_integer_t nei,
                               igraph_real_t p, igraph_bool_t loops,
                               igraph_bool_t multiple) {

    igraph_vector_int_t dimvector;
    igraph_vector_bool_t periodic;

    if (dim < 1) {
        IGRAPH_ERROR("WS game: dimension should be at least one", IGRAPH_EINVAL);
    }
    if (size < 1) {
        IGRAPH_ERROR("WS game: lattice size should be at least one",
                     IGRAPH_EINVAL);
    }
    if (p < 0 || p > 1) {
        IGRAPH_ERROR("WS game: rewiring probability should be between 0 and 1",
                     IGRAPH_EINVAL);
    }

    /* Create the lattice first */

    IGRAPH_VECTOR_INT_INIT_FINALLY(&dimvector, dim);
    igraph_vector_int_fill(&dimvector, size);

    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&periodic, dim);
    igraph_vector_bool_fill(&periodic, true);

    IGRAPH_CHECK(igraph_square_lattice(graph, &dimvector, nei, IGRAPH_UNDIRECTED,
                                /* mutual */ false, &periodic));

    igraph_vector_bool_destroy(&periodic);
    igraph_vector_int_destroy(&dimvector);
    IGRAPH_FINALLY_CLEAN(2);
    IGRAPH_FINALLY(igraph_destroy, graph);

    /* Rewire the edges then */

    IGRAPH_CHECK(igraph_rewire_edges(graph, p, loops, multiple));

    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
