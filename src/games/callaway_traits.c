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
#include "igraph_memory.h"
#include "igraph_random.h"

/**
 * \function igraph_callaway_traits_game
 * \brief Simulates a growing network with vertex types.
 *
 * </para><para>
 * The different types of vertices prefer to connect other types of
 * vertices with a given probability.</para><para>
 *
 * </para><para>
 * The simulation goes like this: in each discrete time step a new
 * vertex is added to the graph. The type of this vertex is generated
 * based on \p type_dist. Then two vertices are selected uniformly
 * randomly from the graph. The probability that they will be
 * connected depends on the types of these vertices and is taken from
 * \p pref_matrix. Then another two vertices are selected and this is
 * repeated \p edges_per_step times in each time step.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * D. S. Callaway, J. E. Hopcroft, J. M. Kleinberg, M. E. J. Newman, and S. H. Strogatz,
 * Are randomly grown graphs really random?
 * Phys. Rev. E 64, 041902 (2001).
 * https://doi.org/10.1103/PhysRevE.64.041902
 *
 * \param graph Pointer to an uninitialized graph.
 * \param nodes The number of nodes in the graph.
 * \param types Number of node types.
 * \param edges_per_step The number of connections tried in each time step.
 * \param type_dist Vector giving the distribution of the vertex types.
 * If \c NULL, the distribution is assumed to be uniform.
 * \param pref_matrix Matrix giving the connection probabilities for
 * the vertex types.
 * \param directed Logical, whether to generate a directed graph.
 * \param node_type_vec An initialized vector or \c NULL.
 * If not \c NULL, the type of each node will be stored here.
 * \return Error code.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|*k*log(|V|)), |V| is the number of vertices,
 * k is \p edges_per_step.
 */
int igraph_callaway_traits_game(igraph_t *graph, igraph_integer_t nodes,
                                igraph_integer_t types, igraph_integer_t edges_per_step,
                                const igraph_vector_t *type_dist,
                                const igraph_matrix_t *pref_matrix,
                                igraph_bool_t directed,
                                igraph_vector_t *node_type_vec) {
    long int i, j;
    igraph_vector_t edges;
    igraph_vector_t cumdist;
    igraph_real_t maxcum;
    igraph_vector_t *nodetypes;

    /* Argument contracts */
    if(nodes < 0){
        IGRAPH_ERROR("The number of vertices must be non-negative.", IGRAPH_EINVAL);
    }

    if (types < 1) {
        IGRAPH_ERROR("The number of vertex types must be at least 1.", IGRAPH_EINVAL);
    }

    if (type_dist) {
        igraph_real_t lo;

        if (igraph_vector_size(type_dist) != types) {
            IGRAPH_ERROR("The vertex type distribution vector must agree in length with the number of types.",
                         IGRAPH_EINVAL);
        }

        lo = igraph_vector_min(type_dist);
        if (lo < 0) {
            IGRAPH_ERROR("The vertex type distribution vector must not contain negative values.", IGRAPH_EINVAL);
        }
        if (igraph_is_nan(lo)) {
            IGRAPH_ERROR("The vertex type distribution vector must not contain NaN.", IGRAPH_EINVAL);
        }
    }

    if (igraph_matrix_nrow(pref_matrix) != types || igraph_matrix_ncol(pref_matrix) != types) {
        IGRAPH_ERROR("The preference matrix must be square and agree in dimensions with the number of types.", IGRAPH_EINVAL);
    }

    {
        igraph_real_t lo, hi;
        igraph_matrix_minmax(pref_matrix, &lo, &hi);

        if (lo < 0 || hi > 1) {
            IGRAPH_ERROR("The preference matrix must contain probabilities in [0, 1].", IGRAPH_EINVAL);
        }
        if (igraph_is_nan(lo) || igraph_is_nan(hi)) {
            IGRAPH_ERROR("The preference matrix must not contain NaN.", IGRAPH_EINVAL);
        }
    }

    if (! directed && ! igraph_matrix_is_symmetric(pref_matrix)) {
        IGRAPH_ERROR("The preference matrix must be symmetric when generating undirected graphs.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types + 1);

    if (type_dist) {
        VECTOR(cumdist)[0] = 0;
        for (i = 0; i < types; ++i) {
            VECTOR(cumdist)[i + 1] = VECTOR(cumdist)[i] + VECTOR(*type_dist)[i];
        }
    } else {
        for (i = 0; i < types+1; ++i) {
            VECTOR(cumdist)[i] = i;
        }
    }
    maxcum = igraph_vector_tail(&cumdist);

    if (maxcum <= 0) {
        IGRAPH_ERROR("The vertex type distribution vector must contain at least one positive value.", IGRAPH_EINVAL);
    }

    if (node_type_vec) {
        nodetypes = node_type_vec;
        IGRAPH_CHECK(igraph_vector_resize(nodetypes, nodes));
    } else {
        nodetypes = IGRAPH_CALLOC(1, igraph_vector_t);
        if (! nodetypes) {
            IGRAPH_ERROR("Insufficient memory for callaway_traits_game.", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, nodetypes);
        IGRAPH_VECTOR_INIT_FINALLY(nodetypes, nodes);
    }

    RNG_BEGIN();

    for (i = 0; i < nodes; i++) {
        igraph_real_t uni = RNG_UNIF(0, maxcum);
        long int type;
        igraph_vector_binsearch(&cumdist, uni, &type);
        VECTOR(*nodetypes)[i] = type - 1;
    }

    for (i = 1; i < nodes; i++) {
        for (j = 0; j < edges_per_step; j++) {
            long int node1 = RNG_INTEGER(0, i);
            long int node2 = RNG_INTEGER(0, i);
            long int type1 = (long int) VECTOR(*nodetypes)[node1];
            long int type2 = (long int) VECTOR(*nodetypes)[node2];
            /*    printf("unif: %f, %f, types: %li, %li\n", uni1, uni2, type1, type2); */
            if (RNG_UNIF01() < MATRIX(*pref_matrix, type1, type2)) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, node1));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, node2));
            }
        }
    }

    RNG_END();

    if (! node_type_vec) {
        igraph_vector_destroy(nodetypes);
        IGRAPH_FREE(nodetypes);
        IGRAPH_FINALLY_CLEAN(2);
    }
    igraph_vector_destroy(&cumdist);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
