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
 * \brief Calculate the density of a graph.
 *
 * The density of a graph is simply the ratio of the actual number of its
 * edges and the largest possible number of edges it could have.
 * The maximum number of edges depends on interpretation: are vertices
 * allowed to have a connection to themselves? This is controlled by the
 * \p loops parameter.
 *
 * </para><para>
 * Note that density is ill-defined for graphs which have multiple edges
 * between some pairs of vertices. Consider calling \ref igraph_simplify()
 * on such graphs. This function does not check whether the graph has
 * parallel edges. The result it returns for such graphs is not meaningful.
 *
 * \param graph The input graph object.
 * \param res Pointer to a real number, the result will be stored
 *   here. It must not have parallel edges.
 * \param loops Logical constant, whether to include self-loops in the
 *   calculation. If this constant is \c true then
 *   loop edges are thought to be possible in the graph (this does not
 *   necessarily mean that the graph really contains any loops). If
 *   this is \c false then the result is only correct if the graph does not
 *   contain loops.
 * \return Error code.
 *
 * Time complexity: O(1).
 */
igraph_error_t igraph_density(const igraph_t *graph, igraph_real_t *res,
                   igraph_bool_t loops) {

    igraph_real_t no_of_nodes = (igraph_real_t) igraph_vcount(graph);
    igraph_real_t no_of_edges = (igraph_real_t) igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);

    if (no_of_nodes == 0) {
        *res = IGRAPH_NAN;
        return IGRAPH_SUCCESS;
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

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_diversity
 * \brief Structural diversity index of the vertices.
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
 * (not-a-number), while that of vertices with a single connection
 * will be zero.
 *
 * </para><para>
 * The measure works only if the graph is undirected and has no multiple edges.
 * If the graph has multiple edges, simplify it first using \ref
 * igraph_simplify(). If the graph is directed, convert it into an undirected
 * graph with \ref igraph_to_undirected() .
 *
 * \param graph The undirected input graph.
 * \param weights The edge weights, in the order of the edge IDs, must
 *    have appropriate length. Weights must be non-negative.
 * \param res An initialized vector, the results are stored here.
 * \param vids Vertex selector that specifies the vertices which to calculate
 *    the measure.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear.
 *
 */
igraph_error_t igraph_diversity(const igraph_t *graph, const igraph_vector_t *weights,
                     igraph_vector_t *res, const igraph_vs_t vids) {

    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t k, i;
    igraph_vector_int_t incident;
    igraph_bool_t has_multiple;
    igraph_vit_t vit;

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

    if (no_of_edges > 0) {
        igraph_real_t minweight = igraph_vector_min(weights);
        if (minweight < 0) {
            IGRAPH_ERROR("Weight vector must be non-negative.", IGRAPH_EINVAL);
        } else if (isnan(minweight)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
        }
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&incident, 10);

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    igraph_vector_clear(res);
    IGRAPH_CHECK(igraph_vector_reserve(res, IGRAPH_VIT_SIZE(vit)));

    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        igraph_real_t d;
        igraph_integer_t v = IGRAPH_VIT_GET(vit);

        IGRAPH_CHECK(igraph_incident(graph, &incident, v, /*mode=*/ IGRAPH_ALL));
        k = igraph_vector_int_size(&incident); /* degree */

        /*
         * Non-normalized diversity is defined as
         * d = -sum_i w_i/s log (w_i/s)
         * where s = sum_i w_i. In order to avoid two passes through the w vector,
         * we use the equivalent formulation of
         * d = log s - (sum_i w_i log w_i) / s
         * However, this formulation may not give an exact 0.0 for some w when k=1,
         * due to roundoff errors (examples: w=3 or w=7). For this reason, we
         * special-case the computation for k=1 even for the unnormalized diversity
         * insted of just setting the normalization factor to 1 for this case.
         */
        if (k == 0) {
            d = IGRAPH_NAN;
        } else if (k == 1) {
            if (VECTOR(*weights)[0] > 0) d = 0.0; /* s > 0 */
            else d = IGRAPH_NAN; /* s == 0 */
        } else {
            igraph_real_t s = 0.0, ent = 0.0;
            for (i = 0; i < k; i++) {
                igraph_real_t w = VECTOR(*weights)[VECTOR(incident)[i]];
                if (w == 0) continue;
                s += w;
                ent += (w * log(w));
            }
            d = (log(s) - ent / s) / log(k);
        }

        igraph_vector_push_back(res, d); /* reserved */
    }

    igraph_vit_destroy(&vit);
    igraph_vector_int_destroy(&incident);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_reciprocity
 * \brief Calculates the reciprocity of a directed graph.
 *
 * The measure of reciprocity defines the proportion of mutual
 * connections, in a directed graph. It is most commonly defined as
 * the probability that the opposite counterpart of a randomly chosen
 * directed edge is also included in the graph. In adjacency matrix
 * notation: <code>1 - (sum_ij |A_ij - A_ji|) / (2 sum_ij A_ij)</code>.
 * In multigraphs, each parallel edges between two vertices must
 * have its own separate reciprocal edge, in accordance with the
 * above formula. This measure is calculated if the \p mode argument is
 * \c IGRAPH_RECIPROCITY_DEFAULT.
 *
 * </para><para>
 * For directed graphs with no edges, NaN is returned.
 * For undirected graphs, 1 is returned unconditionally.
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
 * \param ignore_loops Whether to ignore self-loops when counting edges.
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
igraph_error_t igraph_reciprocity(const igraph_t *graph, igraph_real_t *res,
                       igraph_bool_t ignore_loops,
                       igraph_reciprocity_t mode) {

    igraph_integer_t nonrec = 0, rec = 0, loops = 0;
    igraph_vector_int_t inneis, outneis;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    if (mode != IGRAPH_RECIPROCITY_DEFAULT &&
        mode != IGRAPH_RECIPROCITY_RATIO) {
        IGRAPH_ERROR("Invalid reciprocity type.", IGRAPH_EINVAL);
    }

    /* Undirected graphs has reciprocity 1.0 by definition. */
    if (!igraph_is_directed(graph)) {
        *res = 1.0;
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&inneis, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&outneis, 0);

    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_integer_t ip, op;
        IGRAPH_CHECK(igraph_neighbors(graph, &inneis, i, IGRAPH_IN));
        IGRAPH_CHECK(igraph_neighbors(graph, &outneis, i, IGRAPH_OUT));

        ip = op = 0;
        while (ip < igraph_vector_int_size(&inneis) &&
               op < igraph_vector_int_size(&outneis)) {
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
        nonrec += (igraph_vector_int_size(&inneis) - ip) +
                  (igraph_vector_int_size(&outneis) - op);
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

    igraph_vector_int_destroy(&inneis);
    igraph_vector_int_destroy(&outneis);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
