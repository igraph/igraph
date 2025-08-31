/*
   igraph library.
   Copyright (C) 2011-2025  The igraph development team <igraph@igraph.org>

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

#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_interrupt.h"

#include "core/exceptions.h"

#include "config.h"

#ifdef HAVE_INFOMAP
    #include "Infomap.h"
#endif

#include <climits>
#include <cmath>

static igraph_error_t infomap_get_membership(infomap::InfomapBase &infomap, igraph_vector_int_t *membership) {
    igraph_int_t n = infomap.numLeafNodes();

    IGRAPH_CHECK(igraph_vector_int_resize(membership, n));

    for (auto it(infomap.iterTreePhysical(1)); !it.isEnd(); ++it) {
        infomap::InfoNode &node = *it;
        if (node.isLeaf()) {
            // Note: We must use moduleIndex() and not moduleId(), as the latter
            // may be >= vcount, which causes igraph_reindex_membership() to fail.
            VECTOR(*membership)[node.physicalId] = it.moduleIndex();
        }
    }

    // Re-index membership
    IGRAPH_CHECK(igraph_reindex_membership(membership, NULL, NULL));

    return IGRAPH_SUCCESS;
}

static igraph_error_t convert_igraph_to_infomap(const igraph_t *graph,
                                      const igraph_vector_t *edge_weights,
                                      const igraph_vector_t *vertex_weights,
                                      infomap::Network &network) {

    igraph_int_t vcount = igraph_vcount(graph);
    igraph_int_t ecount = igraph_ecount(graph);

    if (vcount > UINT_MAX) {
        IGRAPH_ERROR("Graph has too many vertices for Infomap.", IGRAPH_EINVAL);
    }

    for (igraph_int_t v = 0; v < vcount; v++) {
        if (vertex_weights) {
            double weight = VECTOR(*vertex_weights)[v];
            if (weight < 0) {
                IGRAPH_ERRORF("Vertex weights must not be negative, got %g.",
                              IGRAPH_EINVAL, weight);
            }
            if (! std::isfinite(weight)) {
                IGRAPH_ERRORF("Vertex weights must not be infinite or NaN, got %g.",
                              IGRAPH_EINVAL, weight);
            }
            network.addNode(v, weight);
        } else {
            network.addNode(v);
        }
    }

    for (igraph_int_t e = 0; e < ecount; e++) {
        igraph_int_t v1 = IGRAPH_FROM(graph, e);
        igraph_int_t v2 = IGRAPH_TO(graph, e);

        if (edge_weights) {
            double weight = VECTOR(*edge_weights)[e];
            if (weight < 0) {
                IGRAPH_ERRORF("Edge weights must not be negative, got %g.",
                              IGRAPH_EINVAL, weight);
            }
            if (! std::isfinite(weight)) {
                IGRAPH_ERRORF("Edge weights must not be infinite or NaN, got %g.",
                              IGRAPH_EINVAL, weight);
            }
            network.addLink(v1, v2, weight);
        } else {
            network.addLink(v1, v2);
        }
    }

    return IGRAPH_SUCCESS;
}

// Needed in case C++'s bool is not compatible with igraph's igraph_bool_t
// which may happen in some configurations on some platforms, notable with R/igraph.
static bool infomap_allow_interruption() {
    return igraph_allow_interruption();
}

/**
 * \function igraph_community_infomap
 * \brief Community structure that minimizes the expected description length of a random walker trajectory.
 *
 * Implementation of the Infomap community detection algorithm of
 * Martin Rosvall and Carl T. Bergstrom. This algorithm takes edge directions
 * into account. For more details, see the visualization of the math and the
 * map generator at https://www.mapequation.org.
 *
 * </para><para>
 * Infomap is based on a random walker model similar to PageRank: the walker
 * either chooses out-edges to follow with probabilities proportional to edge
 * weights, or teleports to a random vertex with probability 0.15. Vertex weights
 * can be given to control the probability of choosing different vertices as
 * the target of the teleportation. In addition, Infomap can be regularized to
 * account for potential missing links.
 *
 * </para><para>
 * As of igraph 1.0, the Infomap library written by Daniel Edler, Anton Holmgren
 * and Martin Rosvall is used. See https://github.com/mapequation/infomap/.
 *
 * </para><para>
 * If you want to specify a random seed (as in the original
 * implementation) you can use \ref igraph_rng_seed().
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * M. Rosvall and C. T. Bergstrom:
 * Maps of information flow reveal community structure in complex networks,
 * PNAS 105, 1118 (2008).
 * https://dx.doi.org/10.1073/pnas.0706851105, https://arxiv.org/abs/0707.0609
 *
 * </para><para>
 * M. Rosvall, D. Axelsson, and C. T. Bergstrom:
 * The map equation,
 * Eur. Phys. J. Special Topics 178, 13 (2009).
 * https://dx.doi.org/10.1140/epjst/e2010-01179-1, https://arxiv.org/abs/0906.1405
 *
 * </para><para>
 * Smiljanić, Jelena, Daniel Edler, and Martin Rosvall: Mapping Flows on
 * Sparse Networks with Missing Links. Phys Rev E 102 (1–1): 012302 (2020).
 * https://doi.org/10.1103/PhysRevE.102.012302, https://arxiv.org/abs/2106.14798
 *
 * \param graph The input graph. Edge directions are taken into account.
 * \param edge_weights Numeric vector giving the weights of the edges.
 *     The random walker will favour edges with high weights over
 *     edges with low weights; the probability of picking a particular
 *     outbound edge from a node is directly proportional to its weight.
 *     If it is \c NULL then all edges will have equal
 *     weights. The weights are expected to be non-negative.
 * \param vertex_weights Numeric vector giving the weights of the vertices.
 *     Vertices with higher weights are favoured by the random walker
 *     when it teleports to a new vertex. The probability of picking a vertex
 *     when the random walker teleports is directly proportional to the weight
 *     of the vertex. If this argument is \c NULL then all vertices will have
 *     equal weights. Weights are expected to be positive.
 * \param nb_trials The number of attempts to partition the network
 *     (can be any integer value equal to or larger than 1).
 * \param is_regularized If true, adds a fully connected Bayesian prior network
 *     to avoid overfitting due to missing links.
 * \param regularization_strength Adjust relative strength of the Bayesian prior
 *     network used for regularization. This multiplies the default strength, a
 *     parameter of 1 hence uses the default regularization strength. Ignored
 *     when \p is_regularized is set to \c false.
 * \param membership Pointer to a vector. The membership vector is
 *     stored here. \c NULL means that the caller is not interested in the
 *     membership vector.
 * \param codelength Pointer to a real. If not \c NULL the code length of the
 *     partition is stored here.
 * \return Error code.
 *
 * \sa \ref igraph_community_spinglass(), \ref
 * igraph_community_edge_betweenness(), \ref igraph_community_walktrap().
 *
 * Time complexity: TODO.
 */

igraph_error_t igraph_community_infomap(
        const igraph_t *graph,
        const igraph_vector_t *edge_weights,
        const igraph_vector_t *vertex_weights,
        igraph_int_t nb_trials,
        igraph_bool_t is_regularized,
        igraph_real_t regularization_strength,
        igraph_vector_int_t *membership,
        igraph_real_t *codelength) {

#ifndef HAVE_INFOMAP
    IGRAPH_ERROR("Infomap is not available.", IGRAPH_UNIMPLEMENTED);
#else

    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_int_t ecount = igraph_ecount(graph);

    if (edge_weights) {
        if (igraph_vector_size(edge_weights) != ecount) {
            IGRAPH_ERROR("Length of edge weight vector does not match edge count.",
                         IGRAPH_EINVAL);
        }
    }

    if (vertex_weights) {
        if (igraph_vector_size(vertex_weights) != vcount) {
            IGRAPH_ERROR("Length of vertex weight vector does not match edge count.",
                         IGRAPH_EINVAL);
        }
    }

    if (nb_trials < 1) {
        IGRAPH_ERRORF("Number of trials must be at least 1, got %" IGRAPH_PRId ".",
                      IGRAPH_EINVAL,
                      nb_trials);
    }

    // Handle null graph
    if (vcount == 0) {
        if (membership) {
            IGRAPH_CHECK(igraph_vector_int_resize(membership, 0));
        }

        if (codelength) {
            *codelength = IGRAPH_NAN;
        }

        return IGRAPH_SUCCESS;
    }

    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    // Configure infomap
    infomap::Config conf;
    conf.twoLevel = true;
    conf.numTrials = nb_trials;
    conf.silent = true;
    conf.directed = igraph_is_directed(graph);
    conf.interruptionHandler = &infomap_allow_interruption;
    conf.regularized = is_regularized;
    conf.regularizationStrength = regularization_strength;

    infomap::InfomapBase infomap(conf);

    IGRAPH_CHECK(convert_igraph_to_infomap(graph, edge_weights, vertex_weights, infomap.network()));

    infomap.run();

    if (membership) {
        IGRAPH_CHECK(infomap_get_membership(infomap, membership));
    }

    if (codelength) {
        *codelength = infomap.codelength();
    }

    IGRAPH_HANDLE_EXCEPTIONS_END;

    return IGRAPH_SUCCESS;

#endif
}
