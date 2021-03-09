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
 * This function generates a graph according to the Watts-Strogatz
 * model of small-world networks. The graph is obtained by creating a
 * circular undirected lattice and then rewire the edges randomly with
 * a constant probability.
 *
 * </para><para>See also: Duncan J Watts and Steven H Strogatz:
 * Collective dynamics of <quote>small world</quote> networks, Nature
 * 393, 440-442, 1998.
 *
 * \param graph The graph to initialize.
 * \param dim The dimension of the lattice.
 * \param size The size of the lattice along each dimension.
 * \param nei The size of the neighborhood for each vertex. This is
 *    the same as the \p nei argument of \ref
 *    igraph_connect_neighborhood().
 * \param p The rewiring probability. A real number between zero and
 *   one (inclusive).
 * \param loops Logical, whether to generate loop edges.
 * \param multiple Logical, whether to allow multiple edges in the
 *   generated graph.
 * \return Error code.
 *
 * \sa \ref igraph_lattice(), \ref igraph_connect_neighborhood() and
 * \ref igraph_rewire_edges() can be used if more flexibility is
 * needed, e.g. a different type of lattice.
 *
 * Time complexity: O(|V|*d^o+|E|), |V| and |E| are the number of
 * vertices and edges, d is the average degree, o is the \p nei
 * argument.
 */
int igraph_watts_strogatz_game(igraph_t *graph, igraph_integer_t dim,
                               igraph_integer_t size, igraph_integer_t nei,
                               igraph_real_t p, igraph_bool_t loops,
                               igraph_bool_t multiple) {

    igraph_vector_t dimvector;
    long int i;

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

    IGRAPH_VECTOR_INIT_FINALLY(&dimvector, dim);
    for (i = 0; i < dim; i++) {
        VECTOR(dimvector)[i] = size;
    }

    IGRAPH_CHECK(igraph_lattice(graph, &dimvector, nei, IGRAPH_UNDIRECTED,
                                0 /* mutual */, 1 /* circular */));
    igraph_vector_destroy(&dimvector);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_destroy, graph);

    /* Rewire the edges then */

    IGRAPH_CHECK(igraph_rewire_edges(graph, p, loops, multiple));

    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
