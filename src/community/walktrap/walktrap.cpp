/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

/* The original version of this file was written by Pascal Pons
   The original copyright notice follows here. The FSF address was
   fixed by Tamas Nepusz */

// File: walktrap.cpp
//-----------------------------------------------------------------------------
// Walktrap v0.2 -- Finds community structure of networks using random walks
// Copyright (C) 2004-2005 Pascal Pons
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//-----------------------------------------------------------------------------
// Author   : Pascal Pons
// Email    : pascal.pons@gmail.com
// Web page : http://www-rp.lip6.fr/~latapy/PP/walktrap.html
// Location : Paris, France
// Time     : June 2005
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include "walktrap_graph.h"
#include "walktrap_communities.h"

#include "igraph_community.h"
#include "igraph_components.h"
#include "igraph_interface.h"

#include "core/exceptions.h"
#include "core/interruption.h"

#include <climits>
#include <cmath>

// This is necessary for GCC 5 and earlier, where including <cmath>
// makes isnan() unusable without the std:: prefix, even if <math.h>
// was included as well.
using std::isnan;

using namespace igraph::walktrap;

/**
 * \function igraph_community_walktrap
 * \brief Community finding using a random walk based similarity measure.
 *
 * This function is the implementation of the Walktrap community
 * finding algorithm, see Pascal Pons, Matthieu Latapy: Computing
 * communities in large networks using random walks,
 * https://arxiv.org/abs/physics/0512106
 *
 * </para><para>
 * Currently the original C++ implementation is used in igraph,
 * see https://www-complexnetworks.lip6.fr/~latapy/PP/walktrap.html
 * We are grateful to Matthieu Latapy and Pascal Pons for providing this
 * source code.
 *
 * </para><para>
 * In contrast to the original implementation, isolated vertices are allowed
 * in the graph and they are assumed to have a single incident loop edge with
 * weight 1.
 *
 * \param graph The input graph, edge directions are ignored.
 * \param weights Numeric vector giving the weights of the edges.
 *     If it is a NULL pointer then all edges will have equal
 *     weights. The weights are expected to be positive.
 * \param steps Integer constant, the length of the random walks.
 *     Typically, good results are obtained with values between
 *     3-8 with 4-5 being a reasonable default.
 * \param merges Pointer to a matrix, the merges performed by the
 *     algorithm will be stored here (if not \c NULL). Each merge is a
 *     row in a two-column matrix and contains the IDs of the merged
 *     clusters. Clusters are numbered from zero and cluster numbers
 *     smaller than the number of nodes in the network belong to the
 *     individual vertices as singleton clusters. In each step a new
 *     cluster is created from two other clusters and its id will be
 *     one larger than the largest cluster id so far. This means that
 *     before the first merge we have \c n clusters (the number of
 *     vertices in the graph) numbered from zero to \c n-1. The first
 *     merge creates cluster \c n, the second cluster \c n+1, etc.
 * \param modularity Pointer to a vector. If not \c NULL then the
 *     modularity score of the current clustering is stored here after
 *     each merge operation.
 * \param membership Pointer to a vector. If not a \c NULL pointer, then
 *     the membership vector corresponding to the maximal modularity
 *     score is stored here.
 * \return Error code.
 *
 * \sa \ref igraph_community_spinglass(), \ref
 * igraph_community_edge_betweenness().
 *
 * Time complexity: O(|E||V|^2) in the worst case, O(|V|^2 log|V|) typically,
 * |V| is the number of vertices, |E| is the number of edges.
 *
 * \example examples/simple/walktrap.c
 */

igraph_error_t igraph_community_walktrap(const igraph_t *graph,
                              const igraph_vector_t *weights,
                              igraph_integer_t steps,
                              igraph_matrix_int_t *merges,
                              igraph_vector_t *modularity,
                              igraph_vector_int_t *membership) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t comp_count;
    igraph_matrix_int_t imerges, *pmerges = merges;
    igraph_vector_t imodularity, *pmodularity = modularity;

    if (steps <= 0) {
        IGRAPH_ERROR("Length of random walks must be positive for walktrap community detection.", IGRAPH_EINVAL);
    }

    if (steps > INT_MAX) {
        IGRAPH_ERROR("Length of random walks too large for walktrap community detection.", IGRAPH_EINVAL);
    }

    int length = steps;

    if (weights) {
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
        }

        if (no_of_edges > 0) {
            igraph_real_t minweight = igraph_vector_min(weights);
            if (minweight < 0) {
                IGRAPH_ERROR("Weight vector must be non-negative.", IGRAPH_EINVAL);
            } else if (isnan(minweight)) {
                IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
            }
        }
    }

    if (membership) {
        /* We need both 'modularity' and 'merges' to compute 'membership'.
         * If they were not provided by the called, we allocate these here. */

        if (! modularity) {
            IGRAPH_VECTOR_INIT_FINALLY(&imodularity, 0);
            pmodularity = &imodularity;
        }

        if (! merges) {
            IGRAPH_MATRIX_INT_INIT_FINALLY(&imerges, 0, 0);
            pmerges = &imerges;
        }
    }

    IGRAPH_HANDLE_EXCEPTIONS(
        Graph G;
        IGRAPH_CHECK(G.convert_from_igraph(graph, weights));

        if (pmerges || pmodularity) {
            IGRAPH_CHECK(igraph_connected_components(graph, /*membership=*/ NULL, /*csize=*/ NULL,
                                                     &comp_count, IGRAPH_WEAK));
        }
        if (pmerges) {
            IGRAPH_CHECK(igraph_matrix_int_resize(pmerges, no_of_nodes - comp_count, 2));
        }
        if (pmodularity) {
            IGRAPH_CHECK(igraph_vector_resize(pmodularity, no_of_nodes - comp_count + 1));
            igraph_vector_null(pmodularity);
        }
        Communities C(&G, length, pmerges, pmodularity);

        while (!C.H->is_empty()) {
            IGRAPH_ALLOW_INTERRUPTION();
            C.merge_nearest_communities();
        }
    );

    if (membership) {
        igraph_integer_t m;
        m = no_of_nodes > 0 ? igraph_vector_which_max(pmodularity) : 0;
        IGRAPH_CHECK(igraph_community_to_membership(pmerges, no_of_nodes,
                    /*steps=*/ m,
                    membership,
                    /*csize=*/ NULL));

        if (! merges) {
            igraph_matrix_int_destroy(&imerges);
            IGRAPH_FINALLY_CLEAN(1);
        }
        if (! modularity) {
            igraph_vector_destroy(&imodularity);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    /* The walktrap implementation cannot work with NaN values internally,
     * and produces 0 for the modularity of edgeless graphs. We correct
     * this to NaN in the last step for consistency. */
    if (modularity && no_of_edges == 0) {
        VECTOR(*modularity)[0] = IGRAPH_NAN;
    }

    return IGRAPH_SUCCESS;
}
