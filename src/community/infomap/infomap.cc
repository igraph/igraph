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

#include <cmath>
#include "igraph_interface.h"
#include "igraph_community.h"
#include "core/interruption.h"


#include "infomap_Node.h"
#include "infomap_Greedy.h"

/****************************************************************************/
int infomap_partition(FlowGraph * fgraph, bool rcall) {
    Greedy * greedy;

    // save the original graph
    FlowGraph * cpy_fgraph = new FlowGraph(fgraph);
    IGRAPH_FINALLY(delete_FlowGraph, cpy_fgraph);

    int Nnode = cpy_fgraph->Nnode;
    // "real" number of vertex, ie. number of vertex of the graph

    int iteration = 0;
    double outer_oldCodeLength, newCodeLength;

    int *initial_move = NULL;
    bool initial_move_done = true;

    do { // Main loop
        outer_oldCodeLength = fgraph->codeLength;

        if (iteration > 0) {
            /**********************************************************************/
            //  FIRST PART: re-split the network (if need)
            // ===========================================

            // intial_move indicate current clustering
            initial_move = new int[Nnode];
            // new_cluster_id --> old_cluster_id (save curent clustering state)

            IGRAPH_FINALLY(operator delete [], initial_move);
            initial_move_done = false;

            int *subMoveTo = NULL; // enventual new partitionment of original graph

            if ((iteration % 2 == 0) && (fgraph->Nnode > 1)) {
                // 0/ Submodule movements : partition each module of the
                // current partition (rec. call)

                subMoveTo = new int[Nnode];
                // vid_cpy_fgraph  --> new_cluster_id (new partition)

                IGRAPH_FINALLY(operator delete [], subMoveTo);

                int subModIndex = 0;

                for (int i = 0 ; i < fgraph->Nnode ; i++) {
                    // partition each non trivial module
                    int sub_Nnode = fgraph->node[i]->members.size();
                    if (sub_Nnode > 1) { // If the module is not trivial
                        int *sub_members  = new int[sub_Nnode];      // id_sub --> id
                        IGRAPH_FINALLY(operator delete [], sub_members);

                        for (int j = 0 ; j < sub_Nnode ; j++) {
                            sub_members[j] = fgraph->node[i]->members[j];
                        }

                        // extraction of the subgraph
                        FlowGraph *sub_fgraph = new FlowGraph(cpy_fgraph, sub_Nnode,
                                                              sub_members);
                        IGRAPH_FINALLY(delete_FlowGraph, sub_fgraph);
                        sub_fgraph->initiate();

                        // recursif call of partitionment on the subgraph
                        infomap_partition(sub_fgraph, true);

                        // Record membership changes
                        for (int j = 0; j < sub_fgraph->Nnode; j++) {
                            int Nmembers = sub_fgraph->node[j]->members.size();
                            for (int k = 0; k < Nmembers; k++) {
                                subMoveTo[sub_members[sub_fgraph->node[j]->members[k]]] =
                                    subModIndex;
                            }
                            initial_move[subModIndex] = i;
                            subModIndex++;
                        }

                        delete sub_fgraph;
                        IGRAPH_FINALLY_CLEAN(1);
                        delete [] sub_members;
                        IGRAPH_FINALLY_CLEAN(1);
                    } else {
                        subMoveTo[fgraph->node[i]->members[0]] = subModIndex;
                        initial_move[subModIndex] = i;
                        subModIndex++;
                    }
                }
            } else {
                // 1/ Single-node movements : allows each node to move (again)
                // save current modules
                for (int i = 0; i < fgraph->Nnode; i++) { // for each module
                    int Nmembers = fgraph->node[i]->members.size(); // Module size
                    for (int j = 0; j < Nmembers; j++) { // for each vertex (of the module)
                        initial_move[fgraph->node[i]->members[j]] = i;
                    }
                }
            }

            fgraph->back_to(cpy_fgraph);
            if (subMoveTo) {
                Greedy *cpy_greedy = new Greedy(fgraph);
                IGRAPH_FINALLY(delete_Greedy, cpy_greedy);

                cpy_greedy->setMove(subMoveTo);
                cpy_greedy->apply(false);

                delete_Greedy(cpy_greedy);
                IGRAPH_FINALLY_CLEAN(1);
                delete [] subMoveTo;
                IGRAPH_FINALLY_CLEAN(1);
            }
        }
        /**********************************************************************/
        //  SECOND PART: greedy optimizing it self
        // ===========================================
        double oldCodeLength;

        do {
            // greedy optimizing object creation
            greedy = new Greedy(fgraph);
            IGRAPH_FINALLY(delete_Greedy, greedy);

            // Initial move to apply ?
            if (!initial_move_done && initial_move) {
                initial_move_done = true;
                greedy->setMove(initial_move);
            }

            oldCodeLength = greedy->codeLength;
            bool moved = true;
            int Nloops = 0;
            //int count = 0;
            double inner_oldCodeLength = 1000;

            while (moved) { // main greedy optimizing loop
                inner_oldCodeLength = greedy->codeLength;
                moved = greedy->optimize();

                Nloops++;
                //count++;

                if (fabs(greedy->codeLength - inner_oldCodeLength) < 1.0e-10)
                    // if the move does'n reduce the codelenght -> exit !
                {
                    moved = false;
                }

                //if (count == 10) {
                //  greedy->tune();
                //  count = 0;
                //}
            }

            // transform the network to network of modules:
            greedy->apply(true);
            newCodeLength = greedy->codeLength;

            // destroy greedy object
            delete greedy;
            IGRAPH_FINALLY_CLEAN(1);

        } while (oldCodeLength - newCodeLength >  1.0e-10);
        // while there is some improvement

        if (iteration > 0) {
            delete [] initial_move;
            IGRAPH_FINALLY_CLEAN(1);
        }

        iteration++;
        if (!rcall) {
            IGRAPH_ALLOW_INTERRUPTION();
        }
    } while (outer_oldCodeLength - newCodeLength > 1.0e-10);

    delete cpy_fgraph;
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_community_infomap
 * \brief Find community structure that minimizes the expected
 * description length of a random walker trajectory.
 *
 * Implementation of the InfoMap community detection algorithm.of
 * Martin Rosvall and Carl T. Bergstrom.
 *
 * See :
 * Visualization of the math and the map generator: www.mapequation.org
 * [2] The original paper: M. Rosvall and C. T. Bergstrom, Maps of
 * information flow reveal community structure in complex networks, PNAS
 * 105, 1118 (2008) [http://dx.doi.org/10.1073/pnas.0706851105 ,
 * http://arxiv.org/abs/0707.0609 ]
 * [3] A more detailed paper: M. Rosvall, D. Axelsson, and C. T. Bergstrom,
 * The map equation, Eur. Phys. J. Special Topics 178, 13 (2009).
 * [http://dx.doi.org/10.1140/epjst/e2010-01179-1 ,
 * http://arxiv.org/abs/0906.1405 ]

 * </para><para>
 * The original C++ implementation of Martin Rosvall is used,
 * see http://www.tp.umu.se/~rosvall/downloads/infomap_undir.tgz .
 * Intergation in igraph has be done by Emmanuel Navarro (who is grateful to
  * Martin Rosvall and Carl T. Bergstrom for providing this source code.)
 *
 * </para><para>
 * Note that the graph must not contain isolated vertices.
 *
 * </para><para>
 * If you want to specify a random seed (as in original
 * implementation) you can use \ref igraph_rng_seed().
 *
 * \param graph The input graph.
 * \param e_weights Numeric vector giving the weights of the edges.
 *     If it is a NULL pointer then all edges will have equal
 *     weights. The weights are expected to be positive.
 * \param v_weights Numeric vector giving the weights of the vertices.
 *     If it is a NULL pointer then all vertices will have equal
 *     weights. The weights are expected to be positive.
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
int igraph_community_infomap(const igraph_t * graph,
                             const igraph_vector_t *e_weights,
                             const igraph_vector_t *v_weights,
                             int nb_trials,
                             igraph_vector_t *membership,
                             igraph_real_t *codelength) {

    FlowGraph * fgraph = new FlowGraph(graph, e_weights, v_weights);
    IGRAPH_FINALLY(delete_FlowGraph, fgraph);

    // compute stationary distribution
    fgraph->initiate();

    FlowGraph * cpy_fgraph ;
    double shortestCodeLength = 1000.0;

    // create membership vector
    int Nnode = fgraph->Nnode;
    IGRAPH_CHECK(igraph_vector_resize(membership, Nnode));

    for (int trial = 0; trial < nb_trials; trial++) {
        cpy_fgraph = new FlowGraph(fgraph);
        IGRAPH_FINALLY(delete_FlowGraph, cpy_fgraph);

        //partition the network
        IGRAPH_CHECK(infomap_partition(cpy_fgraph, false));

        // if better than the better...
        if (cpy_fgraph->codeLength < shortestCodeLength) {
            shortestCodeLength = cpy_fgraph->codeLength;
            // ... store the partition
            for (int i = 0 ; i < cpy_fgraph->Nnode ; i++) {
                int Nmembers = cpy_fgraph->node[i]->members.size();
                for (int k = 0; k < Nmembers; k++) {
                    //cluster[ cpy_fgraph->node[i]->members[k] ] = i;
                    VECTOR(*membership)[cpy_fgraph->node[i]->members[k]] = i;
                }
            }
        }

        delete_FlowGraph(cpy_fgraph);
        IGRAPH_FINALLY_CLEAN(1);
    }

    *codelength = (igraph_real_t) shortestCodeLength / log(2.0);

    delete fgraph;
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_reindex_membership(membership, 0, 0));

    return IGRAPH_SUCCESS;
}
