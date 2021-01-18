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
#include "igraph_nongraph.h"
#include "igraph_random.h"

/**
 * \function igraph_establishment_game
 * \brief Generates a graph with a simple growing model with vertex types.
 *
 * </para><para>
 * The simulation goes like this: a single vertex is added at each
 * time step. This new vertex tries to connect to \p k vertices in the
 * graph. The probability that such a connection is realized depends
 * on the types of the vertices involved.
 *
 * \param graph Pointer to an uninitialized graph.
 * \param nodes The number of vertices in the graph.
 * \param types The number of vertex types.
 * \param k The number of connections tried in each time step.
 * \param type_dist Vector giving the distribution of vertex types.
 * \param pref_matrix Matrix giving the connection probabilities for
 * different vertex types.
 * \param directed Logical, whether to generate a directed graph.
 * \return Error code.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|*k*log(|V|)), |V| is the number of vertices
 * and k is the \p k parameter.
 */
int igraph_establishment_game(igraph_t *graph, igraph_integer_t nodes,
                              igraph_integer_t types, igraph_integer_t k,
                              igraph_vector_t *type_dist,
                              igraph_matrix_t *pref_matrix,
                              igraph_bool_t directed) {

    long int i, j;
    igraph_vector_t edges;
    igraph_vector_t cumdist;
    igraph_vector_t potneis;
    igraph_real_t maxcum;
    igraph_vector_t nodetypes;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types + 1);
    IGRAPH_VECTOR_INIT_FINALLY(&potneis, k);
    IGRAPH_VECTOR_INIT_FINALLY(&nodetypes, nodes);

    VECTOR(cumdist)[0] = 0;
    for (i = 0; i < types; i++) {
        VECTOR(cumdist)[i + 1] = VECTOR(cumdist)[i] + VECTOR(*type_dist)[i];
    }
    maxcum = igraph_vector_tail(&cumdist);

    RNG_BEGIN();

    for (i = 0; i < nodes; i++) {
        igraph_real_t uni = RNG_UNIF(0, maxcum);
        long int type;
        igraph_vector_binsearch(&cumdist, uni, &type);
        VECTOR(nodetypes)[i] = type - 1;
    }

    for (i = k; i < nodes; i++) {
        long int type1 = (long int) VECTOR(nodetypes)[i];
        igraph_random_sample(&potneis, 0, i - 1, k);
        for (j = 0; j < k; j++) {
            long int type2 = (long int) VECTOR(nodetypes)[(long int)VECTOR(potneis)[j]];
            if (RNG_UNIF01() < MATRIX(*pref_matrix, type1, type2)) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, VECTOR(potneis)[j]));
            }
        }
    }

    RNG_END();

    igraph_vector_destroy(&nodetypes);
    igraph_vector_destroy(&potneis);
    igraph_vector_destroy(&cumdist);
    IGRAPH_FINALLY_CLEAN(3);
    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
