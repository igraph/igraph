/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "core/exceptions.h"
#include "core/interruption.h"

#include "infomap_Node.h"
#include "infomap_FlowGraph.h"
#include "infomap_Greedy.h"

#include <cmath>
#include <vector>

// This is necessary for GCC 5 and earlier, where including <cmath>
// makes isnan() unusable without the std:: prefix, even if <math.h>
// was included as well.
using std::isnan;

/****************************************************************************/
static igraph_error_t infomap_partition(FlowGraph &fgraph, bool rcall) {

    // save the original graph
    FlowGraph cpy_fgraph(fgraph);

    igraph_integer_t Nnode = cpy_fgraph.Nnode;
    // "real" number of vertex, ie. number of vertex of the graph

    igraph_integer_t iteration = 0;
    double outer_oldCodeLength, newCodeLength;

    std::vector<igraph_integer_t> initial_move;
    bool initial_move_done = true;

    // re-use vector in loop for better performance
    std::vector<igraph_integer_t> subMoveTo;

    do { // Main loop
        outer_oldCodeLength = fgraph.codeLength;

        if (iteration > 0) {
            /**********************************************************************/
            //  FIRST PART: re-split the network (if need)
            // ===========================================

            // intial_move indicate current clustering
            initial_move.resize(Nnode);
            // new_cluster_id --> old_cluster_id (save curent clustering state)

            initial_move_done = false;

            subMoveTo.clear(); // enventual new partitionment of original graph

            if ((iteration % 2 == 0) && (fgraph.Nnode > 1)) {
                // 0/ Submodule movements : partition each module of the
                // current partition (rec. call)

                subMoveTo.resize(Nnode);
                // vid_cpy_fgraph  --> new_cluster_id (new partition)

                igraph_integer_t subModIndex = 0;

                for (igraph_integer_t i = 0 ; i < fgraph.Nnode ; i++) {
                    // partition each non trivial module
                    size_t sub_Nnode = fgraph.node[i].members.size();
                    if (sub_Nnode > 1) { // If the module is not trivial
                        const std::vector<igraph_integer_t> &sub_members = fgraph.node[i].members;

                        // extraction of the subgraph
                        FlowGraph sub_fgraph(cpy_fgraph, sub_members);
                        sub_fgraph.initiate();

                        // recursif call of partitionment on the subgraph
                        infomap_partition(sub_fgraph, true);

                        // Record membership changes
                        for (igraph_integer_t j = 0; j < sub_fgraph.Nnode; j++) {
                            for (const auto &v : sub_fgraph.node[j].members) {
                                subMoveTo[sub_members[v]] = subModIndex;
                            }
                            initial_move[subModIndex] = i;
                            subModIndex++;
                        }
                    } else {
                        subMoveTo[fgraph.node[i].members[0]] = subModIndex;
                        initial_move[subModIndex] = i;
                        subModIndex++;
                    }
                }
            } else {
                // 1/ Single-node movements : allows each node to move (again)
                // save current modules
                for (igraph_integer_t i = 0; i < fgraph.Nnode; i++) { // for each module
                    for (const auto &v : fgraph.node[i].members) { // for each vertex (of the module)
                        initial_move[v] = i;
                    }
                }
            }

            fgraph.back_to(cpy_fgraph);
            if (! subMoveTo.empty()) {
                Greedy cpy_greedy(&fgraph);

                cpy_greedy.setMove(subMoveTo);
                cpy_greedy.apply(false);
            }
        }
        /**********************************************************************/
        //  SECOND PART: greedy optimizing it self
        // ===========================================
        double oldCodeLength;

        do {
            // greedy optimizing object creation
            Greedy greedy(&fgraph);

            // Initial move to apply ?
            if (!initial_move_done && ! initial_move.empty()) {
                initial_move_done = true;
                greedy.setMove(initial_move);
            }

            oldCodeLength = greedy.codeLength;
            bool moved = true;
            double inner_oldCodeLength = 1000;

            while (moved) { // main greedy optimizing loop
                inner_oldCodeLength = greedy.codeLength;
                moved = greedy.optimize();

                if (fabs(greedy.codeLength - inner_oldCodeLength) < 1.0e-10)
                    // if the move does'n reduce the codelenght -> exit !
                {
                    moved = false;
                }
            }

            // transform the network to network of modules:
            greedy.apply(true);
            newCodeLength = greedy.codeLength;
        } while (oldCodeLength - newCodeLength >  1.0e-10);
        // while there is some improvement

        iteration++;
        if (!rcall) {
            IGRAPH_ALLOW_INTERRUPTION();
        }
    } while (outer_oldCodeLength - newCodeLength > 1.0e-10);

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
 *    stored here.
 * \param codelength Pointer to a real. If not NULL the code length of the
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

    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    if (e_weights) {
        const igraph_integer_t ecount = igraph_ecount(graph);
        if (igraph_vector_size(e_weights) != ecount) {
            IGRAPH_ERROR("Invalid edge weight vector length.", IGRAPH_EINVAL);
        }
        if (ecount > 0) {
            /* Allow both positive and zero weights.
             * The conversion to Infomap format will simply skip zero-weight edges/ */
            igraph_real_t minweight = igraph_vector_min(e_weights);
            if (minweight < 0) {
                IGRAPH_ERROR("Edge weights must not be negative.", IGRAPH_EINVAL);
            } else if (isnan(minweight)) {
                IGRAPH_ERROR("Edge weights must not be NaN values.", IGRAPH_EINVAL);
            }
        }
    }

    if (v_weights) {
        const igraph_integer_t vcount = igraph_vcount(graph);
        if (igraph_vector_size(v_weights) != vcount) {
            IGRAPH_ERROR("Invalid vertex weight vector length.", IGRAPH_EINVAL);
        }
        if (vcount > 0) {
            /* TODO: Currently we require strictly positive. Can this be
             * relaxed to non-negative values? */
            igraph_real_t minweight = igraph_vector_min(v_weights);
            if (minweight <= 0) {
                IGRAPH_ERROR("Vertex weights must be positive.", IGRAPH_EINVAL);
            } else if (isnan(minweight)) {
                IGRAPH_ERROR("Vertex weights must not be NaN values.", IGRAPH_EINVAL);
            }
        }
    }

    FlowGraph fgraph(graph, e_weights, v_weights);

    // compute stationary distribution
    fgraph.initiate();

    double shortestCodeLength = 1000.0;

    // create membership vector
    igraph_integer_t Nnode = fgraph.Nnode;
    IGRAPH_CHECK(igraph_vector_int_resize(membership, Nnode));

    for (igraph_integer_t trial = 0; trial < nb_trials; trial++) {
        FlowGraph cpy_fgraph(fgraph);

        //partition the network
        IGRAPH_CHECK(infomap_partition(cpy_fgraph, false));

        // if better than the better...
        if (cpy_fgraph.codeLength < shortestCodeLength) {
            shortestCodeLength = cpy_fgraph.codeLength;
            // ... store the partition
            for (igraph_integer_t i = 0 ; i < cpy_fgraph.Nnode ; i++) {
                size_t Nmembers = cpy_fgraph.node[i].members.size();
                for (size_t k = 0; k < Nmembers; k++) {
                    //cluster[ cpy_fgraph->node[i].members[k] ] = i;
                    VECTOR(*membership)[cpy_fgraph.node[i].members[k]] = i;
                }
            }
        }
    }

    *codelength = shortestCodeLength / log(2.0);

    IGRAPH_CHECK(igraph_reindex_membership(membership, NULL, NULL));

    return IGRAPH_SUCCESS;

    IGRAPH_HANDLE_EXCEPTIONS_END;
}
