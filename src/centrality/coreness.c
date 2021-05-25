/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "igraph_community.h"

#include "igraph_memory.h"
#include "igraph_interface.h"
#include "igraph_iterators.h"

/**
 * \function igraph_coreness
 * \brief Finding the coreness of the vertices in a network.
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

int igraph_coreness(const igraph_t *graph, igraph_vector_t *cores,
                    igraph_neimode_t mode) {

    long int no_of_nodes = igraph_vcount(graph);
    long int *bin, *vert, *pos;
    long int maxdeg;
    long int i, j = 0;
    igraph_vector_t neis;
    igraph_neimode_t omode;

    if (mode != IGRAPH_ALL && mode != IGRAPH_OUT && mode != IGRAPH_IN) {
        IGRAPH_ERROR("Invalid mode in k-cores", IGRAPH_EINVAL);
    }
    if (!igraph_is_directed(graph) || mode == IGRAPH_ALL) {
        mode = omode = IGRAPH_ALL;
    } else if (mode == IGRAPH_IN) {
        omode = IGRAPH_OUT;
    } else {
        omode = IGRAPH_IN;
    }

    if (no_of_nodes == 0) {
        igraph_vector_clear(cores);
        return IGRAPH_SUCCESS;
    }

    vert = IGRAPH_CALLOC(no_of_nodes, long int);
    if (vert == 0) {
        IGRAPH_ERROR("Cannot calculate k-cores", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, vert);
    pos = IGRAPH_CALLOC(no_of_nodes, long int);
    if (pos == 0) {
        IGRAPH_ERROR("Cannot calculate k-cores", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, pos);

    /* maximum degree + degree of vertices */
    IGRAPH_CHECK(igraph_degree(graph, cores, igraph_vss_all(), mode,
                               IGRAPH_LOOPS));
    maxdeg = (long int) igraph_vector_max(cores);

    bin = IGRAPH_CALLOC(maxdeg + 1, long int);
    if (bin == 0) {
        IGRAPH_ERROR("Cannot calculate k-cores", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, bin);

    /* degree histogram */
    for (i = 0; i < no_of_nodes; i++) {
        bin[ (long int)VECTOR(*cores)[i] ] += 1;
    }

    /* start pointers */
    j = 0;
    for (i = 0; i <= maxdeg; i++) {
        long int k = bin[i];
        bin[i] = j;
        j += k;
    }

    /* sort in vert (and corrupt bin) */
    for (i = 0; i < no_of_nodes; i++) {
        pos[i] = bin[(long int)VECTOR(*cores)[i]];
        vert[pos[i]] = i;
        bin[(long int)VECTOR(*cores)[i]] += 1;
    }

    /* correct bin */
    for (i = maxdeg; i > 0; i--) {
        bin[i] = bin[i - 1];
    }
    bin[0] = 0;

    /* this is the main algorithm */
    IGRAPH_VECTOR_INIT_FINALLY(&neis, maxdeg);
    for (i = 0; i < no_of_nodes; i++) {
        long int v = vert[i];
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) v, omode));
        for (j = 0; j < igraph_vector_size(&neis); j++) {
            long int u = (long int) VECTOR(neis)[j];
            if (VECTOR(*cores)[u] > VECTOR(*cores)[v]) {
                long int du = (long int) VECTOR(*cores)[u];
                long int pu = pos[u];
                long int pw = bin[du];
                long int w = vert[pw];
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

    igraph_vector_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_free(bin);
    igraph_free(pos);
    igraph_free(vert);
    IGRAPH_FINALLY_CLEAN(3);
    return 0;
}
