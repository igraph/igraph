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
#include "core/interruption.h"

using namespace igraph::walktrap;

/**
 * \function igraph_community_walktrap
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
 * \param merges Pointer to a matrix, the merges performed by the
 *     algorithm will be stored here (if not NULL). Each merge is a
 *     row in a two-column matrix and contains the ids of the merged
 *     clusters. Clusters are numbered from zero and cluster numbers
 *     smaller than the number of nodes in the network belong to the
 *     individual vertices as singleton clusters. In each step a new
 *     cluster is created from two other clusters and its id will be
 *     one larger than the largest cluster id so far. This means that
 *     before the first merge we have \c n clusters (the number of
 *     vertices in the graph) numbered from zero to \c n-1. The first
 *     merge creates cluster \c n, the second cluster \c n+1, etc.
 * \param modularity Pointer to a vector. If not NULL then the
 *     modularity score of the current clustering is stored here after
 *     each merge operation.
 * \param membership Pointer to a vector. If not a NULL pointer, then
 *     the membership vector corresponding to the maximal modularity
 *     score is stored here. If it is not a NULL pointer, then neither
 *     \p modularity nor \p merges may be NULL.
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

int igraph_community_walktrap(const igraph_t *graph,
                              const igraph_vector_t *weights,
                              int steps,
                              igraph_matrix_t *merges,
                              igraph_vector_t *modularity,
                              igraph_vector_t *membership) {

    long int no_of_nodes = (long int)igraph_vcount(graph);
    int length = steps;
    long max_memory = -1;

    if (membership && !(modularity && merges)) {
        IGRAPH_ERROR("Cannot calculate membership without modularity or merges",
                     IGRAPH_EINVAL);
    }

    Graph G;
    if (G.convert_from_igraph(graph, weights)) {
        IGRAPH_ERROR("Cannot convert igraph graph into walktrap format", IGRAPH_EINVAL);
    }

    if (merges) {
        igraph_integer_t no;
        IGRAPH_CHECK(igraph_clusters(graph, /*membership=*/ 0, /*csize=*/ 0,
                                     &no, IGRAPH_WEAK));
        IGRAPH_CHECK(igraph_matrix_resize(merges, no_of_nodes - no, 2));
    }
    if (modularity) {
        IGRAPH_CHECK(igraph_vector_resize(modularity, no_of_nodes));
        igraph_vector_null(modularity);
    }
    Communities C(&G, length, max_memory, merges, modularity);

    while (!C.H->is_empty()) {
        IGRAPH_ALLOW_INTERRUPTION();
        C.merge_nearest_communities();
    }

    if (membership) {
        long int m;
        m = no_of_nodes > 0 ? igraph_vector_which_max(modularity) : 0;
        IGRAPH_CHECK(igraph_community_to_membership(merges, no_of_nodes,
                     /*steps=*/ m,
                     membership,
                     /*csize=*/ NULL));
    }

    return IGRAPH_SUCCESS;
}
