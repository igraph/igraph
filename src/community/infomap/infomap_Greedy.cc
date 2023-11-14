/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
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

*/

#include "infomap_Greedy.h"

#include <algorithm>
#include <iterator>
#include <map>

using namespace std;

Greedy::Greedy(FlowGraph *fgraph) :
    graph(fgraph),
    Nnode(graph->Nnode),
    alpha(graph->alpha), // teleportation probability
    beta(1.0 - alpha),   // probability to take normal step

    node_index(Nnode),

    mod_empty(Nnode),
    mod_exit(Nnode),
    mod_size(Nnode),
    mod_danglingSize(Nnode),
    mod_teleportWeight(Nnode),
    mod_members(Nnode)
{
    nodeSize_log_nodeSize = graph->nodeSize_log_nodeSize;
    exit_log_exit         = graph->exit_log_exit;
    size_log_size         = graph->size_log_size;
    exitFlow              = graph->exitFlow;

    const std::vector<Node> &node = graph->node;
    for (igraph_integer_t i = 0; i < Nnode; i++) { // For each module
        node_index[i]         = i;
        mod_exit[i]           = node[i].exit;
        mod_size[i]           = node[i].size;

        mod_danglingSize[i]   = node[i].danglingSize;
        mod_teleportWeight[i] = node[i].teleportWeight;
        mod_members[i]        = node[i].members.size();
    }

    exit = plogp(exitFlow);

    codeLength = exit - 2.0 * exit_log_exit + size_log_size - nodeSize_log_nodeSize;
}


/** Greedy optimizing (as in Blodel and Al.) :
 * for each vertex (selected in a random order) compute the best possible move within neighborhood
 */
bool Greedy::optimize() {
    bool moved = false;
    const std::vector<Node> &node = graph->node;

    RNG_BEGIN();

    // Generate random enumeration of nodes
    vector<igraph_integer_t> randomOrder(Nnode);
    for (igraph_integer_t i = 0; i < Nnode; i++) {
        randomOrder[i] = i;
    }

    for (igraph_integer_t i = 0; i < Nnode - 1; i++) {
        igraph_integer_t randPos = RNG_INTEGER(i, Nnode - 1);
        // swap i & randPos
        igraph_integer_t tmp = randomOrder[i];
        randomOrder[i]       = randomOrder[randPos];
        randomOrder[randPos] = tmp;
    }

    igraph_uint_t offset = 1;
    vector<igraph_integer_t> redirect(Nnode, 0);
    vector<pair<igraph_integer_t, pair<double, double> > > flowNtoM(Nnode);

    for (igraph_integer_t k = 0; k < Nnode; k++) {

        // Pick nodes in random order
        igraph_integer_t flip = randomOrder[k];
        igraph_integer_t oldM = node_index[flip];

        // Reset offset when igraph_integer_t overflows
        if (offset > IGRAPH_INTEGER_MAX) {
            for (igraph_integer_t j = 0; j < Nnode; j++) {
                redirect[j] = 0;
            }
            offset = 1;
        }
        // Size of vector with module links
        igraph_integer_t NmodLinks = 0;
        // For all outLinks
        size_t NoutLinks = node[flip].outLinks.size();
        if (NoutLinks == 0) { //dangling node, add node to calculate flow below
            redirect[oldM] = offset + NmodLinks;
            flowNtoM[NmodLinks].first = oldM;
            flowNtoM[NmodLinks].second.first = 0.0;
            flowNtoM[NmodLinks].second.second = 0.0;
            NmodLinks++;
        } else {
            for (size_t j = 0; j < NoutLinks; j++) {
                igraph_integer_t nb_M = node_index[node[flip].outLinks[j].first];
                // index destination du lien
                double nb_flow = node[flip].outLinks[j].second;
                // wgt du lien
                if (redirect[nb_M] >= offset) {
                    flowNtoM[redirect[nb_M] - offset].second.first += nb_flow;
                } else {
                    redirect[nb_M] = offset + NmodLinks;
                    flowNtoM[NmodLinks].first = nb_M;
                    flowNtoM[NmodLinks].second.first = nb_flow;
                    flowNtoM[NmodLinks].second.second = 0.0;
                    NmodLinks++;
                }
            }
        }
        // For all inLinks
        size_t NinLinks = node[flip].inLinks.size();
        for (size_t j = 0; j < NinLinks; j++) {
            igraph_integer_t nb_M = node_index[node[flip].inLinks[j].first];
            double nb_flow = node[flip].inLinks[j].second;

            if (redirect[nb_M] >= offset) {
                flowNtoM[redirect[nb_M] - offset].second.second += nb_flow;
            } else {
                redirect[nb_M] = offset + NmodLinks;
                flowNtoM[NmodLinks].first = nb_M;
                flowNtoM[NmodLinks].second.first = 0.0;
                flowNtoM[NmodLinks].second.second = nb_flow;
                NmodLinks++;
            }
        }

        // For teleportation and dangling nodes
        for (igraph_integer_t j = 0; j < NmodLinks; j++) {
            igraph_integer_t newM = flowNtoM[j].first;
            if (newM == oldM) {
                flowNtoM[j].second.first  +=
                    (alpha * node[flip].size + beta * node[flip].danglingSize) *
                    (mod_teleportWeight[oldM] - node[flip].teleportWeight);
                flowNtoM[j].second.second +=
                    (alpha * (mod_size[oldM] - node[flip].size) +
                     beta * (mod_danglingSize[oldM] - node[flip].danglingSize)) *
                    node[flip].teleportWeight;
            } else {
                flowNtoM[j].second.first  +=
                    (alpha * node[flip].size + beta * node[flip].danglingSize) *
                    mod_teleportWeight[newM];
                flowNtoM[j].second.second +=
                    (alpha * mod_size[newM]   + beta * mod_danglingSize[newM]  ) *
                    node[flip].teleportWeight;
            }
        }

        // Calculate flow to/from own module (default value if no link to
        // own module)
        double outFlowOldM =
            (alpha * node[flip].size + beta * node[flip].danglingSize) *
            (mod_teleportWeight[oldM] - node[flip].teleportWeight) ;
        double inFlowOldM  =
            (alpha * (mod_size[oldM] - node[flip].size) +
             beta * (mod_danglingSize[oldM] - node[flip].danglingSize)) *
            node[flip].teleportWeight;
        if (redirect[oldM] >= offset) {
            outFlowOldM = flowNtoM[redirect[oldM] - offset].second.first;
            inFlowOldM  = flowNtoM[redirect[oldM] - offset].second.second;
        }

        // Option to move to empty module (if node not already alone)
        if (mod_members[oldM] > node[flip].members.size()) {
            if (Nempty > 0) {
                flowNtoM[NmodLinks].first = mod_empty[Nempty - 1];
                flowNtoM[NmodLinks].second.first = 0.0;
                flowNtoM[NmodLinks].second.second = 0.0;
                NmodLinks++;
            }
        }

        // Randomize link order for optimized search
        for (igraph_integer_t j = 0; j < NmodLinks - 1; j++) {
            igraph_integer_t randPos = RNG_INTEGER(j, NmodLinks - 1);
            igraph_integer_t tmp_M = flowNtoM[j].first;
            double tmp_outFlow = flowNtoM[j].second.first;
            double tmp_inFlow = flowNtoM[j].second.second;
            flowNtoM[j].first = flowNtoM[randPos].first;
            flowNtoM[j].second.first = flowNtoM[randPos].second.first;
            flowNtoM[j].second.second = flowNtoM[randPos].second.second;
            flowNtoM[randPos].first = tmp_M;
            flowNtoM[randPos].second.first = tmp_outFlow;
            flowNtoM[randPos].second.second = tmp_inFlow;
        }

        igraph_integer_t bestM = oldM;
        double best_outFlow = 0.0;
        double best_inFlow = 0.0;
        double best_delta = 0.0;

        // Find the move that minimizes the description length
        for (igraph_integer_t j = 0; j < NmodLinks; j++) {

            igraph_integer_t newM = flowNtoM[j].first;
            double outFlowNewM = flowNtoM[j].second.first;
            double inFlowNewM  = flowNtoM[j].second.second;

            if (newM != oldM) {

                double delta_exit = plogp(exitFlow + outFlowOldM + inFlowOldM -
                                          outFlowNewM - inFlowNewM) - exit;

                double delta_exit_log_exit = - plogp(mod_exit[oldM]) -
                                             plogp(mod_exit[newM]) +
                                             plogp(mod_exit[oldM] - node[flip].exit + outFlowOldM + inFlowOldM)
                                             + plogp(mod_exit[newM] + node[flip].exit - outFlowNewM -
                                                     inFlowNewM);

                double delta_size_log_size = - plogp(mod_exit[oldM] + mod_size[oldM])
                                             - plogp(mod_exit[newM] + mod_size[newM])
                                             + plogp(mod_exit[oldM] + mod_size[oldM] - node[flip].exit -
                                                     node[flip].size + outFlowOldM + inFlowOldM)
                                             + plogp(mod_exit[newM] + mod_size[newM] + node[flip].exit +
                                                     node[flip].size - outFlowNewM - inFlowNewM);

                double deltaL = delta_exit - 2.0 * delta_exit_log_exit +
                                delta_size_log_size;

                if (deltaL - best_delta < -1e-10) {
                    bestM = newM;
                    best_outFlow = outFlowNewM;
                    best_inFlow = inFlowNewM;
                    best_delta = deltaL;
                }
            }
        }

        // Make best possible move
        if (bestM != oldM) {
            //Update empty module vector
            if (mod_members[bestM] == 0) {
                Nempty--;
            }
            if (mod_members[oldM] == node[flip].members.size()) {
                mod_empty[Nempty] = oldM;
                Nempty++;
            }

            exitFlow -= mod_exit[oldM] + mod_exit[bestM];

            exit_log_exit -= plogp(mod_exit[oldM]) + plogp(mod_exit[bestM]);
            size_log_size -= plogp(mod_exit[oldM] + mod_size[oldM]) +
                             plogp(mod_exit[bestM] + mod_size[bestM]);

            mod_exit[oldM]            -= node[flip].exit - outFlowOldM -
                                         inFlowOldM;
            mod_size[oldM]            -= node[flip].size;
            mod_danglingSize[oldM]    -= node[flip].danglingSize;
            mod_teleportWeight[oldM]  -= node[flip].teleportWeight;
            mod_members[oldM]         -= node[flip].members.size();

            mod_exit[bestM]           += node[flip].exit - best_outFlow -
                                         best_inFlow;
            mod_size[bestM]           += node[flip].size;
            mod_danglingSize[bestM]   += node[flip].danglingSize;
            mod_teleportWeight[bestM] += node[flip].teleportWeight;
            mod_members[bestM]        += node[flip].members.size();

            exitFlow += mod_exit[oldM] + mod_exit[bestM];

            // Update terms in map equation

            exit_log_exit += plogp(mod_exit[oldM]) + plogp(mod_exit[bestM]);
            size_log_size += plogp(mod_exit[oldM] + mod_size[oldM]) +
                             plogp(mod_exit[bestM] + mod_size[bestM]);
            exit = plogp(exitFlow);

            // Update code length

            codeLength = exit - 2.0 * exit_log_exit + size_log_size -
                         nodeSize_log_nodeSize;

            node_index[flip] = bestM;
            moved = true;
        }
        offset += Nnode;
    }

    RNG_END();

    return moved;
}

/** Apply the move to the given network
 */
void Greedy::apply(bool sort) {

    //old fct prepare(sort)
    vector<igraph_integer_t> modSnode;  // will give IDs of no-empty modules (nodes)
    modSnode.reserve(Nnode);

    igraph_integer_t Nmod = 0;
    for (igraph_integer_t i = 0; i < Nnode; i++) {
        if (mod_members[i] > 0) {
            Nmod++;
            modSnode.push_back(i);
        }
    }

    if (sort) {
        // sort by mod_size
        std::sort(modSnode.begin(), modSnode.end(),
                  [&](size_t a, size_t b) { return mod_size[a] > mod_size[b]; } );
    }

    // Create the new graph
    FlowGraph tmp_fgraph(Nmod);
    vector<Node> &node_tmp = tmp_fgraph.node ;

    const vector<Node> &node = graph->node;

    vector<igraph_integer_t> nodeInMod(Nnode);

    // creation of new nodes
    for (igraph_integer_t i = 0; i < Nmod; i++) {
        node_tmp[i].members.clear(); // clear membership
        node_tmp[i].exit           =           mod_exit[modSnode[i]];
        node_tmp[i].size           =           mod_size[modSnode[i]];
        node_tmp[i].danglingSize   =   mod_danglingSize[modSnode[i]];
        node_tmp[i].teleportWeight = mod_teleportWeight[modSnode[i]];

        nodeInMod[modSnode[i]]      = i;
    }

    // Calculate outflow of links to different modules
    vector<map<igraph_integer_t, double> > outFlowNtoM(Nmod);

    for (igraph_integer_t i = 0; i < Nnode; i++) {
        igraph_integer_t i_M = nodeInMod[node_index[i]]; //final id of the module of the node i
        // add node members to the module
        copy( node[i].members.begin(), node[i].members.end(),
              back_inserter( node_tmp[i_M].members ) );

        for (const auto &link : node[i].outLinks) {
            igraph_integer_t nb         = link.first;
            igraph_integer_t nb_M       = nodeInMod[node_index[nb]];
            double nb_flow = link.second;
            if (nb != i) {
                // inserts key nb_M if it does not exist
                outFlowNtoM[i_M][nb_M] += nb_flow;
            }
        }
    }

    // Create outLinks at new level
    for (igraph_integer_t i = 0; i < Nmod; i++) {
        for (const auto &item : outFlowNtoM[i]) {
            if (item.first != i) {
                node_tmp[i].outLinks.emplace_back(item);
            }
        }
    }

    // Calculate inflow of links from different modules
    vector<map<igraph_integer_t, double> > inFlowNtoM(Nmod);

    for (igraph_integer_t i = 0; i < Nnode; i++) {
        igraph_integer_t i_M = nodeInMod[node_index[i]];
        for (const auto &inLink : node[i].inLinks) {
            igraph_integer_t nb         = inLink.first;
            igraph_integer_t nb_M       = nodeInMod[node_index[nb]];
            double nb_flow = inLink.second;
            if (nb != i) {
                // inserts key nb_M if it does not exist
                inFlowNtoM[i_M][nb_M] += nb_flow;
            }
        }
    }

    // Create inLinks at new level
    for (igraph_integer_t i = 0; i < Nmod; i++) {
        for (const auto &item : inFlowNtoM[i]) {
            if (item.first != i) {
                node_tmp[i].inLinks.emplace_back(item);
            }
        }
    }

    // Option to move to empty module
    mod_empty.clear();
    Nempty = 0;

    //swap node between tmp_graph and graph, then destroy tmp_fgraph
    graph->swap(tmp_fgraph);
    Nnode = Nmod;
}


/**
 * RAZ et recalcul :
 *  - mod_exit
 *  - mod_size
 *  - mod_danglingSize
 *  - mod_teleportWeight
 *  - mod_members
 *  and
 *  - exit_log_exit
 *  - size_log_size
 *  - exitFlow
 *  - exit
 *  - codeLength
 * according to **node / node[i]->index
 */


/* Compute the new CodeSize if modules are merged as indicated by moveTo
 */
void Greedy::setMove(const std::vector<igraph_integer_t> &moveTo) {
    const std::vector<Node> &node = graph->node;
    for (igraph_integer_t i = 0 ; i < Nnode ; i++) { // pour chaque module
        igraph_integer_t oldM = i;
        igraph_integer_t newM = moveTo[i];
        //printf("old -> new : %d -> %d \n", oldM, newM);
        if (newM != oldM) {

            // Si je comprend bien :
            // outFlow... : c'est le "flow" de i-> autre sommet du meme module
            // inFlow... : c'est le "flow" depuis un autre sommet du meme module --> i
            double outFlowOldM = (alpha * node[i].size + beta * node[i].danglingSize) *
                                 (mod_teleportWeight[oldM] - node[i].teleportWeight);
            double inFlowOldM  = (alpha * (mod_size[oldM] - node[i].size) +
                                  beta * (mod_danglingSize[oldM] -
                                          node[i].danglingSize)) *
                                 node[i].teleportWeight;
            double outFlowNewM = (alpha * node[i].size + beta * node[i].danglingSize)
                                 * mod_teleportWeight[newM];
            double inFlowNewM  = (alpha * mod_size[newM] +
                                  beta * mod_danglingSize[newM]) *
                                 node[i].teleportWeight;

            // For all outLinks
            for (const auto &outLink : node[i].outLinks) {
                igraph_integer_t nb_M = node_index[outLink.first];
                double nb_flow = outLink.second;
                if (nb_M == oldM) {
                    outFlowOldM += nb_flow;
                } else if (nb_M == newM) {
                    outFlowNewM += nb_flow;
                }
            }

            // For all inLinks
            for (const auto &inLink : node[i].inLinks) {
                igraph_integer_t nb_M = node_index[inLink.first];
                double nb_flow = inLink.second;
                if (nb_M == oldM) {
                    inFlowOldM += nb_flow;
                } else if (nb_M == newM) {
                    inFlowNewM += nb_flow;
                }
            }

            // Update empty module vector
            // RAZ de mod_empty et Nempty ds calibrate()
            if (mod_members[newM] == 0) {
                // si le nouveau etait vide, on a un vide de moins...
                Nempty--;
            }
            if (mod_members[oldM] == node[i].members.size()) {
                // si l'ancien avait la taille de celui qui bouge, un vide de plus
                mod_empty[Nempty] = oldM;
                Nempty++;
            }

            exitFlow -= mod_exit[oldM] + mod_exit[newM];
            exit_log_exit -= plogp(mod_exit[oldM]) + plogp(mod_exit[newM]);
            size_log_size -= plogp(mod_exit[oldM] + mod_size[oldM]) +
                             plogp(mod_exit[newM] + mod_size[newM]);

            mod_exit[oldM] -= node[i].exit - outFlowOldM - inFlowOldM;
            mod_size[oldM] -= node[i].size;
            mod_danglingSize[oldM] -= node[i].danglingSize;
            mod_teleportWeight[oldM] -= node[i].teleportWeight;
            mod_members[oldM] -= node[i].members.size();
            mod_exit[newM] += node[i].exit - outFlowNewM - inFlowNewM;
            mod_size[newM] += node[i].size;
            mod_danglingSize[newM] += node[i].danglingSize;
            mod_teleportWeight[newM] += node[i].teleportWeight;
            mod_members[newM] += node[i].members.size();

            exitFlow += mod_exit[oldM] + mod_exit[newM];
            exit_log_exit += plogp(mod_exit[oldM]) + plogp(mod_exit[newM]);
            size_log_size += plogp(mod_exit[oldM] + mod_size[oldM]) +
                             plogp(mod_exit[newM] + mod_size[newM]);
            exit = plogp(exitFlow);

            codeLength = exit - 2.0 * exit_log_exit + size_log_size -
                         nodeSize_log_nodeSize;

            node_index[i] = newM;

        }

    }
}
