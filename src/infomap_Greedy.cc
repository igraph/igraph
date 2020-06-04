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
#include <iterator>
#define plogp( x ) ( (x) > 0.0 ? (x)*log(x) : 0.0 )

using namespace std;

Greedy::Greedy(FlowGraph * fgraph) {
    graph = fgraph;
    Nnode = graph->Nnode;

    alpha = graph->alpha;// teleportation probability
    beta = 1.0 - alpha;  // probability to take normal step

    Nempty = 0;
    vector<int>(Nnode).swap(mod_empty);

    vector<int>(Nnode).swap(node_index);
    vector<double>(Nnode).swap(mod_exit);
    vector<double>(Nnode).swap(mod_size);
    vector<double>(Nnode).swap(mod_danglingSize);
    vector<double>(Nnode).swap(mod_teleportWeight);
    vector<int>(Nnode).swap(mod_members);

    nodeSize_log_nodeSize = graph->nodeSize_log_nodeSize;
    exit_log_exit         = graph->exit_log_exit;
    size_log_size         = graph->size_log_size;
    exitFlow              = graph->exitFlow;

    Node ** node = graph->node;
    for (int i = 0; i < Nnode; i++) { // For each module
        node_index[i]         = i;
        mod_exit[i]           = node[i]->exit;
        mod_size[i]           = node[i]->size;

        mod_danglingSize[i]   = node[i]->danglingSize;
        mod_teleportWeight[i] = node[i]->teleportWeight;
        mod_members[i]        = node[i]->members.size();
    }

    exit = plogp(exitFlow);

    codeLength = exit - 2.0 * exit_log_exit + size_log_size -
                 nodeSize_log_nodeSize;
}

Greedy::~Greedy() {
}

void delete_Greedy(Greedy *greedy) {
    delete greedy;
}


/** Greedy optimizing (as in Blodel and Al.) :
 * for each vertex (selected in a random order) compute the best possible move within neighborhood
 */
bool Greedy::optimize() {
    bool moved = false;
    Node ** node = graph->node;

    RNG_BEGIN();

    // Generate random enumeration of nodes
    vector<int> randomOrder(Nnode);
    for (int i = 0; i < Nnode; i++) {
        randomOrder[i] = i;
    }

    for (int i = 0; i < Nnode - 1; i++) {
        //int randPos = i ; //XXX
        int randPos = RNG_INTEGER(i, Nnode - 1);
        // swap i & randPos
        int tmp              = randomOrder[i];
        randomOrder[i]       = randomOrder[randPos];
        randomOrder[randPos] = tmp;
    }

    unsigned int offset = 1;
    vector<unsigned int> redirect(Nnode, 0);
    vector<pair<int, pair<double, double> > > flowNtoM(Nnode);

    for (int k = 0; k < Nnode; k++) {

        // Pick nodes in random order
        int flip = randomOrder[k];
        int oldM = node_index[flip];

        // Reset offset when int overflows
        if (offset > INT_MAX) {
            for (int j = 0; j < Nnode; j++) {
                redirect[j] = 0;
            }
            offset = 1;
        }
        // Size of vector with module links
        int NmodLinks = 0;
        // For all outLinks
        int NoutLinks = node[flip]->outLinks.size();
        if (NoutLinks == 0) { //dangling node, add node to calculate flow below
            redirect[oldM] = offset + NmodLinks;
            flowNtoM[NmodLinks].first = oldM;
            flowNtoM[NmodLinks].second.first = 0.0;
            flowNtoM[NmodLinks].second.second = 0.0;
            NmodLinks++;
        } else {
            for (int j = 0; j < NoutLinks; j++) {
                int nb_M       = node_index[node[flip]->outLinks[j].first];
                // index destination du lien
                double nb_flow = node[flip]->outLinks[j].second;
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
        int NinLinks = node[flip]->inLinks.size();
        for (int j = 0; j < NinLinks; j++) {
            int nb_M = node_index[node[flip]->inLinks[j].first];
            double nb_flow = node[flip]->inLinks[j].second;

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
        for (int j = 0; j < NmodLinks; j++) {
            int newM = flowNtoM[j].first;
            if (newM == oldM) {
                flowNtoM[j].second.first  +=
                    (alpha * node[flip]->size + beta * node[flip]->danglingSize) *
                    (mod_teleportWeight[oldM] - node[flip]->teleportWeight);
                flowNtoM[j].second.second +=
                    (alpha * (mod_size[oldM] - node[flip]->size) +
                     beta * (mod_danglingSize[oldM] - node[flip]->danglingSize)) *
                    node[flip]->teleportWeight;
            } else {
                flowNtoM[j].second.first  +=
                    (alpha * node[flip]->size + beta * node[flip]->danglingSize) *
                    mod_teleportWeight[newM];
                flowNtoM[j].second.second +=
                    (alpha * mod_size[newM]   + beta * mod_danglingSize[newM]  ) *
                    node[flip]->teleportWeight;
            }
        }

        // Calculate flow to/from own module (default value if no link to
        // own module)
        double outFlowOldM =
            (alpha * node[flip]->size + beta * node[flip]->danglingSize) *
            (mod_teleportWeight[oldM] - node[flip]->teleportWeight) ;
        double inFlowOldM  =
            (alpha * (mod_size[oldM] - node[flip]->size) +
             beta * (mod_danglingSize[oldM] - node[flip]->danglingSize)) *
            node[flip]->teleportWeight;
        if (redirect[oldM] >= offset) {
            outFlowOldM = flowNtoM[redirect[oldM] - offset].second.first;
            inFlowOldM  = flowNtoM[redirect[oldM] - offset].second.second;
        }

        // Option to move to empty module (if node not already alone)
        if (mod_members[oldM] > static_cast<int>(node[flip]->members.size())) {
            if (Nempty > 0) {
                flowNtoM[NmodLinks].first = mod_empty[Nempty - 1];
                flowNtoM[NmodLinks].second.first = 0.0;
                flowNtoM[NmodLinks].second.second = 0.0;
                NmodLinks++;
            }
        }

        // Randomize link order for optimized search
        for (int j = 0; j < NmodLinks - 1; j++) {
            //int randPos = j ; // XXX
            int randPos = RNG_INTEGER(j, NmodLinks - 1);
            int tmp_M = flowNtoM[j].first;
            double tmp_outFlow = flowNtoM[j].second.first;
            double tmp_inFlow = flowNtoM[j].second.second;
            flowNtoM[j].first = flowNtoM[randPos].first;
            flowNtoM[j].second.first = flowNtoM[randPos].second.first;
            flowNtoM[j].second.second = flowNtoM[randPos].second.second;
            flowNtoM[randPos].first = tmp_M;
            flowNtoM[randPos].second.first = tmp_outFlow;
            flowNtoM[randPos].second.second = tmp_inFlow;
        }

        int bestM = oldM;
        double best_outFlow = 0.0;
        double best_inFlow = 0.0;
        double best_delta = 0.0;

        // Find the move that minimizes the description length
        for (int j = 0; j < NmodLinks; j++) {

            int newM = flowNtoM[j].first;
            double outFlowNewM = flowNtoM[j].second.first;
            double inFlowNewM  = flowNtoM[j].second.second;

            if (newM != oldM) {

                double delta_exit = plogp(exitFlow + outFlowOldM + inFlowOldM -
                                          outFlowNewM - inFlowNewM) - exit;

                double delta_exit_log_exit = - plogp(mod_exit[oldM]) -
                                             plogp(mod_exit[newM]) +
                                             plogp(mod_exit[oldM] - node[flip]->exit + outFlowOldM + inFlowOldM)
                                             + plogp(mod_exit[newM] + node[flip]->exit - outFlowNewM -
                                                     inFlowNewM);

                double delta_size_log_size = - plogp(mod_exit[oldM] + mod_size[oldM])
                                             - plogp(mod_exit[newM] + mod_size[newM])
                                             + plogp(mod_exit[oldM] + mod_size[oldM] - node[flip]->exit -
                                                     node[flip]->size + outFlowOldM + inFlowOldM)
                                             + plogp(mod_exit[newM] + mod_size[newM] + node[flip]->exit +
                                                     node[flip]->size - outFlowNewM - inFlowNewM);

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
            if (mod_members[oldM] == static_cast<int>(node[flip]->members.size())) {
                mod_empty[Nempty] = oldM;
                Nempty++;
            }

            exitFlow -= mod_exit[oldM] + mod_exit[bestM];

            exit_log_exit -= plogp(mod_exit[oldM]) + plogp(mod_exit[bestM]);
            size_log_size -= plogp(mod_exit[oldM] + mod_size[oldM]) +
                             plogp(mod_exit[bestM] + mod_size[bestM]);

            mod_exit[oldM]            -= node[flip]->exit - outFlowOldM -
                                         inFlowOldM;
            mod_size[oldM]            -= node[flip]->size;
            mod_danglingSize[oldM]    -= node[flip]->danglingSize;
            mod_teleportWeight[oldM]  -= node[flip]->teleportWeight;
            mod_members[oldM]         -= node[flip]->members.size();

            mod_exit[bestM]           += node[flip]->exit - best_outFlow -
                                         best_inFlow;
            mod_size[bestM]           += node[flip]->size;
            mod_danglingSize[bestM]   += node[flip]->danglingSize;
            mod_teleportWeight[bestM] += node[flip]->teleportWeight;
            mod_members[bestM]        += node[flip]->members.size();

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
//void Greedy::level(Node ***node_tmp, bool sort) {

    //old fct prepare(sort)
    vector<int> modSnode;  // will give ids of no-empty modules (nodes)
    int Nmod = 0;
    if (sort) {
        multimap<double, int> Msize;
        for (int i = 0; i < Nnode; i++) {
            if (mod_members[i] > 0) {
                Nmod++;
                Msize.insert(pair<const double, int>(mod_size[i], i));
            }
        }
        for (multimap<double, int>::reverse_iterator it = Msize.rbegin();
             it != Msize.rend(); it++) {
            modSnode.push_back(it->second);
        }
    } else {
        for (int i = 0; i < Nnode; i++) {
            if (mod_members[i] > 0) {
                Nmod++;
                modSnode.push_back(i);
            }
        }
    }
    //modSnode[id_when_no_empty_node] = id_in_mod_tbl

    // Create the new graph
    FlowGraph * tmp_fgraph = new FlowGraph(Nmod);
    IGRAPH_FINALLY(delete_FlowGraph, tmp_fgraph);
    Node ** node_tmp = tmp_fgraph->node ;

    Node ** node = graph->node;

    vector<int> nodeInMod = vector<int>(Nnode);

    // creation of new nodes
    for (int i = 0; i < Nmod; i++) {
        //node_tmp[i] = new Node();
        vector<int>().swap(node_tmp[i]->members); // clear membership
        node_tmp[i]->exit           =           mod_exit[modSnode[i]];
        node_tmp[i]->size           =           mod_size[modSnode[i]];
        node_tmp[i]->danglingSize   =   mod_danglingSize[modSnode[i]];
        node_tmp[i]->teleportWeight = mod_teleportWeight[modSnode[i]];

        nodeInMod[modSnode[i]]      = i;
    }
    //nodeInMode[id_in_mod_tbl] = id_when_no_empty_node

    // Calculate outflow of links to different modules
    vector<map<int, double> > outFlowNtoM(Nmod);
    map<int, double>::iterator it_M;

    for (int i = 0; i < Nnode; i++) {
        int i_M = nodeInMod[node_index[i]]; //final id of the module of the node i
        // add node members to the module
        copy( node[i]->members.begin(), node[i]->members.end(),
              back_inserter( node_tmp[i_M]->members ) );

        int NoutLinks = node[i]->outLinks.size();
        for (int j = 0; j < NoutLinks; j++) {
            int nb         = node[i]->outLinks[j].first;
            int nb_M       = nodeInMod[node_index[nb]];
            double nb_flow = node[i]->outLinks[j].second;
            if (nb != i) {
                it_M = outFlowNtoM[i_M].find(nb_M);
                if (it_M != outFlowNtoM[i_M].end()) {
                    it_M->second += nb_flow;
                } else {
                    outFlowNtoM[i_M].insert(make_pair(nb_M, nb_flow));
                }
            }
        }
    }

    // Create outLinks at new level
    for (int i = 0; i < Nmod; i++) {
        for (it_M = outFlowNtoM[i].begin(); it_M != outFlowNtoM[i].end(); it_M++) {
            if (it_M->first != i) {
                node_tmp[i]->outLinks.push_back(make_pair(it_M->first, it_M->second));
            }
        }
    }

    // Calculate inflow of links from different modules
    vector<map<int, double> > inFlowNtoM(Nmod);

    for (int i = 0; i < Nnode; i++) {
        int i_M = nodeInMod[node_index[i]];
        int NinLinks = node[i]->inLinks.size();
        for (int j = 0; j < NinLinks; j++) {
            int nb         = node[i]->inLinks[j].first;
            int nb_M       = nodeInMod[node_index[nb]];
            double nb_flow = node[i]->inLinks[j].second;
            if (nb != i) {
                it_M = inFlowNtoM[i_M].find(nb_M);
                if (it_M != inFlowNtoM[i_M].end()) {
                    it_M->second += nb_flow;
                } else {
                    inFlowNtoM[i_M].insert(make_pair(nb_M, nb_flow));
                }
            }
        }
    }

    // Create inLinks at new level
    for (int i = 0; i < Nmod; i++) {
        for (it_M = inFlowNtoM[i].begin(); it_M != inFlowNtoM[i].end(); it_M++) {
            if (it_M->first != i) {
                node_tmp[i]->inLinks.push_back(make_pair(it_M->first, it_M->second));
            }
        }
    }

    // Option to move to empty module
    vector<int>().swap(mod_empty);
    Nempty = 0;

    //swap node between tmp_graph and graph, then destroy tmp_fgraph
    graph->swap(tmp_fgraph);
    Nnode = Nmod;

    delete tmp_fgraph;
    IGRAPH_FINALLY_CLEAN(1);
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
void Greedy::tune(void) {

    exit_log_exit = 0.0;
    size_log_size = 0.0;
    exitFlow = 0.0;

    for (int i = 0; i < Nnode; i++) {
        mod_exit[i] = 0.0;
        mod_size[i] = 0.0;
        mod_danglingSize[i] = 0.0;
        mod_teleportWeight[i] = 0.0;
        mod_members[i] = 0;
    }

    Node ** node = graph->node;
    // Update all values except contribution from teleportation
    for (int i = 0; i < Nnode; i++) {
        int i_M = node_index[i]; // module id of node i
        int Nlinks = node[i]->outLinks.size();

        mod_size[i_M]           += node[i]->size;
        mod_danglingSize[i_M]   += node[i]->danglingSize;
        mod_teleportWeight[i_M] += node[i]->teleportWeight;
        mod_members[i_M]++;

        for (int j = 0; j < Nlinks; j++) {
            int neighbor      = node[i]->outLinks[j].first;
            double neighbor_w = node[i]->outLinks[j].second;
            int neighbor_M    = node_index[neighbor];
            if (i_M != neighbor_M) { // neighbor in an other module
                mod_exit[i_M] += neighbor_w;
            }
        }
    }

    // Update contribution from teleportation
    for (int i = 0; i < Nnode; i++) {
        mod_exit[i] += (alpha * mod_size[i] + beta * mod_danglingSize[i]) *
                       (1.0 - mod_teleportWeight[i]);
    }

    for (int i = 0; i < Nnode; i++) {
        exit_log_exit += plogp(mod_exit[i]);
        size_log_size += plogp(mod_exit[i] + mod_size[i]);
        exitFlow += mod_exit[i];
    }
    exit = plogp(exitFlow);

    codeLength = exit - 2.0 * exit_log_exit + size_log_size -
                 nodeSize_log_nodeSize;
}


/* Compute the new CodeSize if modules are merged as indicated by moveTo
 */
void Greedy::setMove(int *moveTo) {
    //void Greedy::determMove(int *moveTo) {
    Node ** node = graph->node;
    //printf("setMove nNode:%d \n", Nnode);
    for (int i = 0 ; i < Nnode ; i++) { // pour chaque module
        int oldM = i;
        int newM = moveTo[i];
        //printf("old -> new : %d -> %d \n", oldM, newM);
        if (newM != oldM) {

            // Si je comprend bien :
            // outFlow... : c'est le "flow" de i-> autre sommet du meme module
            // inFlow... : c'est le "flow" depuis un autre sommet du meme module --> i
            double outFlowOldM = (alpha * node[i]->size + beta * node[i]->danglingSize) *
                                 (mod_teleportWeight[oldM] - node[i]->teleportWeight);
            double inFlowOldM  = (alpha * (mod_size[oldM] - node[i]->size) +
                                  beta * (mod_danglingSize[oldM] -
                                          node[i]->danglingSize)) *
                                 node[i]->teleportWeight;
            double outFlowNewM = (alpha * node[i]->size + beta * node[i]->danglingSize)
                                 * mod_teleportWeight[newM];
            double inFlowNewM  = (alpha * mod_size[newM] +
                                  beta * mod_danglingSize[newM]) *
                                 node[i]->teleportWeight;

            // For all outLinks
            int NoutLinks = node[i]->outLinks.size();
            for (int j = 0; j < NoutLinks; j++) {
                int nb_M = node_index[node[i]->outLinks[j].first];
                double nb_flow = node[i]->outLinks[j].second;
                if (nb_M == oldM) {
                    outFlowOldM += nb_flow;
                } else if (nb_M == newM) {
                    outFlowNewM += nb_flow;
                }
            }

            // For all inLinks
            int NinLinks = node[i]->inLinks.size();
            for (int j = 0; j < NinLinks; j++) {
                int nb_M = node_index[node[i]->inLinks[j].first];
                double nb_flow = node[i]->inLinks[j].second;
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
            if (mod_members[oldM] == static_cast<int>(node[i]->members.size())) {
                // si l'ancien avait la taille de celui qui bouge, un vide de plus
                mod_empty[Nempty] = oldM;
                Nempty++;
            }

            exitFlow -= mod_exit[oldM] + mod_exit[newM];
            exit_log_exit -= plogp(mod_exit[oldM]) + plogp(mod_exit[newM]);
            size_log_size -= plogp(mod_exit[oldM] + mod_size[oldM]) +
                             plogp(mod_exit[newM] + mod_size[newM]);

            mod_exit[oldM] -= node[i]->exit - outFlowOldM - inFlowOldM;
            mod_size[oldM] -= node[i]->size;
            mod_danglingSize[oldM] -= node[i]->danglingSize;
            mod_teleportWeight[oldM] -= node[i]->teleportWeight;
            mod_members[oldM] -= node[i]->members.size();
            mod_exit[newM] += node[i]->exit - outFlowNewM - inFlowNewM;
            mod_size[newM] += node[i]->size;
            mod_danglingSize[newM] += node[i]->danglingSize;
            mod_teleportWeight[newM] += node[i]->teleportWeight;
            mod_members[newM] += node[i]->members.size();

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


