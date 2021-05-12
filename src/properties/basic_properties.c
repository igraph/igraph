/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_structural.h"

#include "igraph_interface.h"

/**
 * \section about_structural
 *
 * <para>These functions usually calculate some structural property
 * of a graph, like its diameter, the degree of the nodes, etc.</para>
 */

/**
 * \function igraph_density
 * Calculate the density of a graph.
 *
 * </para><para>The density of a graph is simply the ratio number of
 * edges and the number of possible edges. Note that density is
 * ill-defined for graphs with multiple and/or loop edges, so consider
 * calling \ref igraph_simplify() on the graph if you know that it
 * contains multiple or loop edges.
 * \param graph The input graph object.
 * \param res Pointer to a real number, the result will be stored
 *   here.
 * \param loops Logical constant, whether to include loops in the
 *   calculation. If this constant is TRUE then
 *   loop edges are thought to be possible in the graph (this does not
 *   necessarily mean that the graph really contains any loops). If
 *   this is FALSE then the result is only correct if the graph does not
 *   contain loops.
 * \return Error code.
 *
 * Time complexity: O(1).
 */
int igraph_density(const igraph_t *graph, igraph_real_t *res,
                   igraph_bool_t loops) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_real_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);

    if (no_of_nodes == 0) {
        *res = IGRAPH_NAN;
        return 0;
    }

    if (!loops) {
        if (no_of_nodes == 1) {
            *res = IGRAPH_NAN;
        } else if (directed) {
            *res = no_of_edges / no_of_nodes / (no_of_nodes - 1);
        } else {
            *res = no_of_edges / no_of_nodes * 2.0 / (no_of_nodes - 1);
        }
    } else {
        if (directed) {
            *res = no_of_edges / no_of_nodes / no_of_nodes;
        } else {
            *res = no_of_edges / no_of_nodes * 2.0 / (no_of_nodes + 1);
        }
    }

    return 0;
}

/**
 * \function igraph_diversity
 * Structural diversity index of the vertices
 *
 * This measure was defined in Nathan Eagle, Michael Macy and Rob
 * Claxton: Network Diversity and Economic Development, Science 328,
 * 1029--1031, 2010.
 *
 * </para><para>
 * It is simply the (normalized) Shannon entropy of the
 * incident edges' weights. D(i)=H(i)/log(k[i]), and
 * H(i) = -sum(p[i,j] log(p[i,j]), j=1..k[i]),
 * where p[i,j]=w[i,j]/sum(w[i,l], l=1..k[i]),  k[i] is the (total)
 * degree of vertex i, and w[i,j] is the weight of the edge(s) between
 * vertex i and j. The diversity of isolated vertices will be NaN
 * (not-a-number).
 *
 * </para><para>
 * The measure works only if the graph is undirected and has no multiple edges.
 * If the graph has multiple edges, simplify it first using \ref
 * igraph_simplify(). If the graph is directed, convert it into an undirected
 * graph with \ref igraph_to_undirected() .
 *
 * \param graph The undirected input graph.
 * \param weights The edge weights, in the order of the edge ids, must
 *    have appropriate length.
 * \param res An initialized vector, the results are stored here.
 * \param vids Vertex selector that specifies the vertices which to calculate
 *    the measure.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear.
 *
 */
int igraph_diversity(igraph_t *graph, const igraph_vector_t *weights,
                     igraph_vector_t *res, const igraph_vs_t vids) {

    int no_of_nodes = igraph_vcount(graph);
    int no_of_edges = igraph_ecount(graph);
    igraph_vector_t incident;
    igraph_vit_t vit;
    igraph_real_t s, ent, w;
    int i, j, k;
    igraph_bool_t has_multiple;

    if (igraph_is_directed(graph)) {
        IGRAPH_ERROR("Diversity measure works with undirected graphs only.", IGRAPH_EINVAL);
    }

    if (!weights) {
        IGRAPH_ERROR("Edge weights must be given.", IGRAPH_EINVAL);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid edge weight vector length.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_has_multiple(graph, &has_multiple));
    if (has_multiple) {
        IGRAPH_ERROR("Diversity measure works only if the graph has no multiple edges.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&incident, 10);

    if (igraph_vs_is_all(&vids)) {
        IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
        for (i = 0; i < no_of_nodes; i++) {
            s = ent = 0.0;
            IGRAPH_CHECK(igraph_incident(graph, &incident, i, /*mode=*/ IGRAPH_ALL));
            for (j = 0, k = (int) igraph_vector_size(&incident); j < k; j++) {
                w = VECTOR(*weights)[(long int)VECTOR(incident)[j]];
                s += w;
                ent += (w * log(w));
            }
            VECTOR(*res)[i] = (log(s) - ent / s) / log(k);
        }
    } else {
        IGRAPH_CHECK(igraph_vector_resize(res, 0));
        IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
        IGRAPH_FINALLY(igraph_vit_destroy, &vit);

        for (IGRAPH_VIT_RESET(vit), i = 0;
             !IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), i++) {
            long int v = IGRAPH_VIT_GET(vit);
            s = ent = 0.0;
            IGRAPH_CHECK(igraph_incident(graph, &incident, (igraph_integer_t) v,
                                         /*mode=*/ IGRAPH_ALL));
            for (j = 0, k = (int) igraph_vector_size(&incident); j < k; j++) {
                w = VECTOR(*weights)[(long int)VECTOR(incident)[j]];
                s += w;
                ent += (w * log(w));
            }
            IGRAPH_CHECK(igraph_vector_push_back(res, (log(s) - ent / s) / log(k)));
        }

        igraph_vit_destroy(&vit);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&incident);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \ingroup structural
 * \function igraph_reciprocity
 * \brief Calculates the reciprocity of a directed graph.
 *
 * </para><para>
 * The measure of reciprocity defines the proportion of mutual
 * connections, in a directed graph. It is most commonly defined as
 * the probability that the opposite counterpart of a directed edge is
 * also included in the graph. In adjacency matrix notation:
 * <code>sum(i, j, (A.*A')ij) / sum(i, j, Aij)</code>, where
 * <code>A.*A'</code> is the element-wise product of matrix
 * <code>A</code> and its transpose. This measure is
 * calculated if the \p mode argument is \c
 * IGRAPH_RECIPROCITY_DEFAULT.
 *
 * </para><para>
 * Prior to igraph version 0.6, another measure was implemented,
 * defined as the probability of mutual connection between a vertex
 * pair if we know that there is a (possibly non-mutual) connection
 * between them. In other words, (unordered) vertex pairs are
 * classified into three groups: (1) disconnected, (2)
 * non-reciprocally connected, (3) reciprocally connected.
 * The result is the size of group (3), divided by the sum of group
 * sizes (2)+(3). This measure is calculated if \p mode is \c
 * IGRAPH_RECIPROCITY_RATIO.
 *
 * \param graph The graph object.
 * \param res Pointer to an \c igraph_real_t which will contain the result.
 * \param ignore_loops Whether to ignore loop edges.
 * \param mode Type of reciprocity to calculate, possible values are
 *    \c IGRAPH_RECIPROCITY_DEFAULT and \c IGRAPH_RECIPROCITY_RATIO,
 *    please see their description above.
 * \return Error code:
 *         \c IGRAPH_EINVAL: graph has no edges
 *         \c IGRAPH_ENOMEM: not enough memory for
 *         temporary data.
 *
 * Time complexity: O(|V|+|E|), |V| is the number of vertices,
 * |E| is the number of edges.
 *
 * \example examples/simple/igraph_reciprocity.c
 */
int igraph_reciprocity(const igraph_t *graph, igraph_real_t *res,
                       igraph_bool_t ignore_loops,
                       igraph_reciprocity_t mode) {

    igraph_integer_t nonrec = 0, rec = 0, loops = 0;
    igraph_vector_t inneis, outneis;
    long int i;
    long int no_of_nodes = igraph_vcount(graph);

    if (mode != IGRAPH_RECIPROCITY_DEFAULT &&
        mode != IGRAPH_RECIPROCITY_RATIO) {
        IGRAPH_ERROR("Invalid reciprocity type", IGRAPH_EINVAL);
    }

    /* THIS IS AN EXIT HERE !!!!!!!!!!!!!! */
    if (!igraph_is_directed(graph)) {
        *res = 1.0;
        return 0;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&inneis, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&outneis, 0);

    for (i = 0; i < no_of_nodes; i++) {
        long int ip, op;
        igraph_neighbors(graph, &inneis, (igraph_integer_t) i, IGRAPH_IN);
        igraph_neighbors(graph, &outneis, (igraph_integer_t) i, IGRAPH_OUT);

        ip = op = 0;
        while (ip < igraph_vector_size(&inneis) &&
               op < igraph_vector_size(&outneis)) {
            if (VECTOR(inneis)[ip] < VECTOR(outneis)[op]) {
                nonrec += 1;
                ip++;
            } else if (VECTOR(inneis)[ip] > VECTOR(outneis)[op]) {
                nonrec += 1;
                op++;
            } else {

                /* loop edge? */
                if (VECTOR(inneis)[ip] == i) {
                    loops += 1;
                    if (!ignore_loops) {
                        rec += 1;
                    }
                } else {
                    rec += 1;
                }

                ip++;
                op++;
            }
        }
        nonrec += (igraph_vector_size(&inneis) - ip) +
                  (igraph_vector_size(&outneis) - op);
    }

    if (mode == IGRAPH_RECIPROCITY_DEFAULT) {
        if (ignore_loops) {
            *res = (igraph_real_t) rec / (igraph_ecount(graph) - loops);
        } else {
            *res = (igraph_real_t) rec / (igraph_ecount(graph));
        }
    } else if (mode == IGRAPH_RECIPROCITY_RATIO) {
        *res = (igraph_real_t) rec / (rec + nonrec);
    }

    igraph_vector_destroy(&inneis);
    igraph_vector_destroy(&outneis);
    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}
