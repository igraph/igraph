/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2023  The igraph development team <igraph@igraph.org>

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

#include "igraph_community.h"

#include "igraph_memory.h"
#include "igraph_interface.h"
#include "igraph_iterators.h"

/**
 * \function igraph_coreness
 * \brief The coreness of the vertices in a graph.
 *
 * The k-core of a graph is a maximal subgraph in which each vertex
 * has at least degree k. (Degree here means the degree in the
 * subgraph of course.). The coreness of a vertex is the highest order
 * of a k-core containing the vertex.
 *
 * </para><para>
 * This function implements the algorithm presented in Vladimir
 * Batagelj, Matjaz Zaversnik: An O(m) Algorithm for Cores
 * Decomposition of Networks.
 * https://arxiv.org/abs/cs/0310049
 *
 * \param graph The input graph.
 * \param cores Pointer to an initialized vector, the result of the
 *        computation will be stored here. It will be resized as
 *        needed. For each vertex it contains the highest order of a
 *        core containing the vertex.
 * \param mode For directed graph it specifies whether to calculate
 *        in-cores, out-cores or the undirected version. It is ignored
 *        for undirected graphs. Possible values: \c IGRAPH_ALL
 *        undirected version, \c IGRAPH_IN in-cores, \c IGRAPH_OUT
 *        out-cores.
 * \return Error code.
 *
 * Time complexity: O(|E|), the number of edges.
 */

igraph_error_t igraph_coreness(const igraph_t *graph,
        igraph_vector_int_t *cores, igraph_neimode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t *bin, *vert, *pos;
    igraph_integer_t maxdeg;
    igraph_vector_int_t neis;
    igraph_neimode_t omode;

    if (mode != IGRAPH_ALL && mode != IGRAPH_OUT && mode != IGRAPH_IN) {
        IGRAPH_ERROR("Invalid mode in k-cores.", IGRAPH_EINVMODE);
    }
    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }
    omode = IGRAPH_REVERSE_MODE(mode);

    /* catch null graph */
    if (no_of_nodes == 0) {
        igraph_vector_int_clear(cores);
        return IGRAPH_SUCCESS;
    }

    vert = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(vert, "Insufficient memory for k-cores.");
    IGRAPH_FINALLY(igraph_free, vert);

    pos = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(pos, "Insufficient memory for k-cores.");
    IGRAPH_FINALLY(igraph_free, pos);

    /* maximum degree + degree of vertices */
    IGRAPH_CHECK(igraph_degree(graph, cores, igraph_vss_all(), mode, /* loops= */ true));

    /* null graph was already handled earlier, 'cores' is not empty */
    maxdeg = igraph_vector_int_max(cores);

    bin = IGRAPH_CALLOC(maxdeg + 1, igraph_integer_t);
    IGRAPH_CHECK_OOM(bin, "Insufficient memory for k-cores.");
    IGRAPH_FINALLY(igraph_free, bin);

    /* degree histogram */
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        bin[ VECTOR(*cores)[i] ] += 1;
    }

    /* start pointers */
    for (igraph_integer_t d = 0, start = 0; d <= maxdeg; d++) {
        igraph_integer_t k = bin[d];
        bin[d] = start;
        start += k;
    }

    /* sort in vert (and corrupt bin) */
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        pos[i] = bin[VECTOR(*cores)[i]];
        vert[pos[i]] = i;
        bin[VECTOR(*cores)[i]] += 1;
    }

    /* correct bin */
    for (igraph_integer_t d = maxdeg; d > 0; d--) {
        bin[d] = bin[d - 1];
    }
    bin[0] = 0;

    /* this is the main algorithm */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, maxdeg);
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_integer_t v = vert[i];
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, v, omode));
        igraph_integer_t nei_count = igraph_vector_int_size(&neis);
        for (igraph_integer_t j = 0; j < nei_count; j++) {
            igraph_integer_t u = VECTOR(neis)[j];
            if (VECTOR(*cores)[u] > VECTOR(*cores)[v]) {
                igraph_integer_t du = VECTOR(*cores)[u];
                igraph_integer_t pu = pos[u];
                igraph_integer_t pw = bin[du];
                igraph_integer_t w = vert[pw];
                if (u != w) {
                    pos[u] = pw;
                    pos[w] = pu;
                    vert[pu] = w;
                    vert[pw] = u;
                }
                bin[du] += 1;
                VECTOR(*cores)[u] -= 1;
            }
        }
    }

    igraph_vector_int_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_free(bin);
    igraph_free(pos);
    igraph_free(vert);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}
