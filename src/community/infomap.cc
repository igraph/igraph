/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2025  The igraph development team
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

   ----
   The original version of this file was written by Martin Rosvall
   email: martin.rosvall@physics.umu.se
   homePage: http://www.tp.umu.se/~rosvall/

   It was integrated in igraph by Emmanuel Navarro
   email: navarro@irit.fr
   homePage: http://www.irit.fr/~Emmanuel.Navarro/
*/

#include "igraph_community.h"

#include "igraph_interface.h"
#include "igraph_random.h"

#include "core/exceptions.h"
#include "core/interruption.h"

#include "config.h"

#ifdef HAVE_INFOMAP
    #include "Infomap.h"
#endif

#include <cmath>
#include <vector>

// This is necessary for GCC 5 and earlier, where including <cmath>
// makes isnan() unusable without the std:: prefix, even if <math.h>
// was included as well.
using std::isnan;

static igraph_error_t infomap_get_membership(infomap::InfomapBase &infomap, igraph_vector_int_t *membership) {
    igraph_integer_t n = infomap.numLeafNodes();

    IGRAPH_CHECK(igraph_vector_int_resize(membership, n));

    for (auto it(infomap.iterTreePhysical(1)); !it.isEnd(); ++it) {
        infomap::InfoNode& node = *it;
        if (node.isLeaf()) {
            VECTOR(*membership)[node.physicalId] = it.moduleId();
        }
    }

    // Re-index membership
    IGRAPH_CHECK(igraph_reindex_membership(membership, 0, 0));

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_to_infomap(const igraph_t *graph,
                                       const igraph_vector_t *e_weights,
                                       const igraph_vector_t *v_weights,
                                       infomap::Network* network) {

    IGRAPH_HANDLE_EXCEPTIONS_BEGIN

    igraph_integer_t n = igraph_vcount(graph);
    igraph_integer_t m = igraph_ecount(graph);

    for (igraph_integer_t v = 0; v < n; v++)
    {
        network->addNode(v, v_weights != NULL ? VECTOR(*v_weights)[v] : 1);
    }

    for (igraph_integer_t e = 0; e < m; e++)
    {
        igraph_integer_t v1 = IGRAPH_FROM(graph, e);
        igraph_integer_t v2 = IGRAPH_TO(graph, e);
        network->addLink(v1, v2, e_weights != NULL ? VECTOR(*e_weights)[e] : 1);
    }

    IGRAPH_HANDLE_EXCEPTIONS_END;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_community_infomap
 * \brief Find community structure that minimizes the expected description length of a random walker trajectory.
 *
 * Implementation of the Infomap community detection algorithm of
 * Martin Rosvall and Carl T. Bergstrom. This algorithm takes edge directions
 * into account.
 *
 * </para><para>
 * For more details, see the visualization of the math and the map generator
 * at https://www.mapequation.org . The original paper describing the algorithm
 * is: M. Rosvall and C. T. Bergstrom, Maps of information flow reveal community
 * structure in complex networks, PNAS 105, 1118 (2008)
 * (https://dx.doi.org/10.1073/pnas.0706851105, https://arxiv.org/abs/0707.0609).
 * A more detailed paper about the algorithm is: M. Rosvall, D. Axelsson, and
 * C. T. Bergstrom, The map equation, Eur. Phys. J. Special Topics 178, 13 (2009).
 * (https://dx.doi.org/10.1140/epjst/e2010-01179-1, https://arxiv.org/abs/0906.1405)

 * </para><para>
 * The original C++ implementation of Martin Rosvall is used,
 * see http://www.tp.umu.se/~rosvall/downloads/infomap_undir.tgz .
 * Integration in igraph was done by Emmanuel Navarro (who is grateful to
 * Martin Rosvall and Carl T. Bergstrom for providing this source code).
 *
 * </para><para>
 * Note that the graph must not contain isolated vertices.
 *
 * </para><para>
 * If you want to specify a random seed (as in the original
 * implementation) you can use \ref igraph_rng_seed().
 *
 * \param graph The input graph. Edge directions are taken into account.
 * \param e_weights Numeric vector giving the weights of the edges.
 *     The random walker will favour edges with high weights over
 *     edges with low weights; the probability of picking a particular
 *     outbound edge from a node is directly proportional to its weight.
 *     If it is \c NULL then all edges will have equal
 *     weights. The weights are expected to be non-negative.
 * \param v_weights Numeric vector giving the weights of the vertices.
 *     Vertices with higher weights are favoured by the random walker
 *     when it needs to "teleport" to a new node after getting stuck in
 *     a sink node (i.e. a node with no outbound edges). The probability
 *     of picking a vertex when the random walker teleports is directly
 *     proportional to the weight of the vertex. If this argument is \c NULL
 *     then all vertices will have equal weights. Weights are expected
 *     to be positive.
 * \param nb_trials The number of attempts to partition the network
 *     (can be any integer value equal or larger than 1).
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

igraph_error_t igraph_community_infomap(const igraph_t * graph,
                             const igraph_vector_t *e_weights,
                             const igraph_vector_t *v_weights,
                             igraph_integer_t nb_trials,
                             igraph_vector_int_t *membership,
                             igraph_real_t *codelength) {

#ifndef HAVE_INFOMAP
    IGRAPH_ERROR("Infomap is not available.", IGRAPH_UNIMPLEMENTED);
#else
    if (igraph_vcount(graph) == 0) {
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
    conf.numTrials = 1;
    conf.silent = true;
    conf.directed = igraph_is_directed(graph);
    conf.interruptionHandler = &igraph_allow_interruption;

    infomap::InfomapBase infomap(conf);

    IGRAPH_CHECK(igraph_to_infomap(graph, e_weights, v_weights, &(infomap.network())));

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
