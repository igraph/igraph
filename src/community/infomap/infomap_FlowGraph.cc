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

#include "infomap_FlowGraph.h"

using namespace std;

void FlowGraph::init(igraph_integer_t n, const igraph_vector_t *v_weights) {
    alpha = 0.15;
    beta  = 1.0 - alpha;
    Nnode = n;
    node.reserve(Nnode);
    if (v_weights) {
        for (igraph_integer_t i = 0; i < Nnode; i++) {
            node.emplace_back(i, VECTOR(*v_weights)[i]);
        }
    } else {
        for (igraph_integer_t i = 0; i < Nnode; i++) {
            node.emplace_back(i, 1.0);
        }
    }
}

FlowGraph::FlowGraph(igraph_integer_t n) {
    init(n, nullptr);
}

/* Build the graph from igraph_t object */
FlowGraph::FlowGraph(const igraph_t *graph,
                     const igraph_vector_t *e_weights,
                     const igraph_vector_t *v_weights) {

    igraph_integer_t n = igraph_vcount(graph);
    init(n, v_weights);

    bool directed = igraph_is_directed(graph);

    double linkWeight = 1.0;
    igraph_integer_t from, to;

    igraph_integer_t Nlinks = igraph_ecount(graph);
    if (!directed) {
        Nlinks = Nlinks * 2 ;
    }
    for (igraph_integer_t i = 0; i < Nlinks; i++) {
        if (!directed) { // not directed
            if (i % 2 == 0) {
                linkWeight = e_weights ? VECTOR(*e_weights)[i / 2] : 1.0;
                igraph_edge(graph, i / 2, &from, &to);
            } else {
                igraph_edge(graph, (i - 1) / 2, &to,   &from);
            }
        } else {         // directed
            linkWeight = e_weights ? VECTOR(*e_weights)[i] : 1.0;
            igraph_edge(graph, i, &from, &to);
        }

        // Populate node from igraph_graph
        // Negative edge weights were checked for already.
        // We skip adding zero-weight edges.
        if (linkWeight > 0.0) {
            if (from != to) {
                node[from].outLinks.emplace_back(to, linkWeight);
                node[to].inLinks.emplace_back(from, linkWeight);
            }
        }
    }
}

FlowGraph::FlowGraph(const FlowGraph &fgraph) {
    igraph_integer_t n = fgraph.Nnode;
    init(n, nullptr);
    for (igraph_integer_t i = 0; i < n; i++) {
        node[i] = fgraph.node[i];
    }

    //XXX: quid de danglings et Ndanglings?

    alpha = fgraph.alpha ;
    beta  = fgraph.beta ;

    exit = fgraph.exit;
    exitFlow = fgraph.exitFlow;
    exit_log_exit = fgraph.exit_log_exit;
    size_log_size = fgraph.size_log_size ;
    nodeSize_log_nodeSize = fgraph.nodeSize_log_nodeSize;

    codeLength = fgraph.codeLength;
}

/** construct a graph by extracting a subgraph from the given graph
 */
FlowGraph::FlowGraph(const FlowGraph &fgraph, const vector<igraph_integer_t> &sub_members) {
    igraph_integer_t sub_Nnode = sub_members.size();

    init(sub_Nnode, nullptr);

    //XXX: use set of integer to ensure that elements are sorted
    set<igraph_integer_t> sub_mem(sub_members.begin(), sub_members.end());

    auto it_mem = sub_mem.begin();

    vector<igraph_integer_t> sub_renumber(fgraph.Nnode, -1);
    // id --> sub_id

    for (igraph_integer_t j = 0; j < sub_Nnode; j++) {
        igraph_integer_t orig_nr = (*it_mem);

        node[j].teleportWeight = fgraph.node[orig_nr].teleportWeight;
        node[j].selfLink       = fgraph.node[orig_nr].selfLink;
        // Take care of self-link

        size_t orig_NoutLinks = fgraph.node[orig_nr].outLinks.size();
        size_t orig_NinLinks  = fgraph.node[orig_nr].inLinks.size();

        sub_renumber[orig_nr] = j;

        for (size_t k = 0; k < orig_NoutLinks; k++) {
            igraph_integer_t to = fgraph.node[orig_nr].outLinks[k].first;
            igraph_integer_t to_newnr = sub_renumber[to];
            double link_weight = fgraph.node[orig_nr].outLinks[k].second;

            if (to < orig_nr) {
                // we add links if the destination (to) has already be seen
                // (ie. smaller than current id) => orig

                if (sub_mem.find(to) != sub_mem.end()) {
                    // printf("%2d | %4d to %4d\n", j, orig_nr, to);
                    // printf("from %4d (%4d:%1.5f) to %4d (%4d)\n", j, orig_nr,
                    //        node[j].selfLink, to_newnr, to);
                    node[j].outLinks.emplace_back(to_newnr, link_weight);
                    node[to_newnr].inLinks.emplace_back(j, link_weight);
                }
            }
        }

        for (size_t k = 0; k < orig_NinLinks; k++) {
            igraph_integer_t to = fgraph.node[orig_nr].inLinks[k].first;
            igraph_integer_t to_newnr = sub_renumber[to];
            double link_weight = fgraph.node[orig_nr].inLinks[k].second;
            if (to < orig_nr) {
                if (sub_mem.find(to) != sub_mem.end()) {
                    node[j].inLinks.emplace_back(to_newnr, link_weight);
                    node[to_newnr].outLinks.emplace_back(j, link_weight);
                }
            }
        }
        it_mem++;
    }
}


/** Swap the graph with the one given
    the graph is "re" calibrate
    but NOT the given one.
 */
void FlowGraph::swap(FlowGraph &fgraph) noexcept {
    node.swap(fgraph.node);

    igraph_integer_t Nnode_tmp = fgraph.Nnode;
    fgraph.Nnode = Nnode;
    Nnode = Nnode_tmp;

    calibrate();
}

/** Initialisation of the graph, compute the flow inside the graph
 *   - count danglings nodes
 *   - normalized edge weights
 *   - Call eigenvector() to compute steady state distribution
 *   - call calibrate to compute codelenght
 */
void FlowGraph::initiate() {
    // Take care of dangling nodes, normalize outLinks, and calculate
    // total teleport weight
    Ndanglings = 0;
    double totTeleportWeight = 0.0;
    for (igraph_integer_t i = 0; i < Nnode; i++) {
        totTeleportWeight += node[i].teleportWeight;
    }

    for (igraph_integer_t i = 0; i < Nnode; i++) {
        node[i].teleportWeight /= totTeleportWeight;
        // normalize teleportation weight

        if (node[i].outLinks.empty() && (node[i].selfLink <= 0.0)) {
            danglings.push_back(i);
            Ndanglings++;
        } else { // Normalize the weights
            size_t NoutLinks = node[i].outLinks.size();
            double sum = node[i].selfLink; // Take care of self-links
            for (size_t j = 0; j < NoutLinks; j++) {
                sum += node[i].outLinks[j].second;
            }
            node[i].selfLink /= sum;
            for (size_t j = 0; j < NoutLinks; j++) {
                node[i].outLinks[j].second /= sum;
            }
        }
    }

    // Calculate steady state matrix
    eigenvector();

    // Update links to represent flow
    for (igraph_integer_t i = 0; i < Nnode; i++) {
        node[i].selfLink = beta * node[i].size * node[i].selfLink;
        //            (1 - \tau) *     \pi_i     *      P_{ii}

        if (!node[i].outLinks.empty()) {
            size_t NoutLinks = node[i].outLinks.size();
            for (size_t j = 0; j < NoutLinks; j++) {
                node[i].outLinks[j].second = beta * node[i].size *
                                              node[i].outLinks[j].second;
                //                      (1 - \tau) *     \pi_i     *          P_{ij}
            }

            // Update values for corresponding inlink
            for (size_t j = 0; j < NoutLinks; j++) {
                size_t NinLinks = node[node[i].outLinks[j].first].inLinks.size();
                for (size_t k = 0; k < NinLinks; k++) {
                    if (node[node[i].outLinks[j].first].inLinks[k].first == i) {
                        node[node[i].outLinks[j].first].inLinks[k].second =
                            node[i].outLinks[j].second;
                        k = NinLinks;
                    }
                }
            }
        }
    }

    // To be able to handle dangling nodes efficiently
    for (igraph_integer_t i = 0; i < Nnode; i++)
        if (node[i].outLinks.empty() && (node[i].selfLink <= 0.0)) {
            node[i].danglingSize = node[i].size;
        } else {
            node[i].danglingSize = 0.0;
        }

    nodeSize_log_nodeSize = 0.0 ;
    // The exit flow from each node at initiation
    for (igraph_integer_t i = 0; i < Nnode; i++) {
        node[i].exit = node[i].size // Proba to be on i
                        - (alpha * node[i].size + beta * node[i].danglingSize) *
                        node[i].teleportWeight // Proba teleport back to i
                        - node[i].selfLink;  // Proba stay on i

        // node[i].exit == q_{i\exit}
        nodeSize_log_nodeSize += plogp(node[i].size);
    }

    calibrate();
}


/* Compute steady state distribution (ie. PageRank) over the network
 * (for all i update node[i].size)
 */
void FlowGraph::eigenvector() {
    vector<double> size_tmp(Nnode, 1.0 / Nnode);

    int Niterations = 0;
    double danglingSize;

    double sqdiff = 1.0;
    double sqdiff_old;
    double sum;
    do {
        // Calculate dangling size
        danglingSize = 0.0;
        for (igraph_integer_t i = 0; i < Ndanglings; i++) {
            danglingSize += size_tmp[danglings[i]];
        }

        // Flow from teleportation
        for (igraph_integer_t i = 0; i < Nnode; i++) {
            node[i].size = (alpha + beta * danglingSize) * node[i].teleportWeight;
        }

        // Flow from network steps
        for (igraph_integer_t i = 0; i < Nnode; i++) {
            node[i].size += beta * node[i].selfLink * size_tmp[i];
            size_t Nlinks = node[i].outLinks.size();
            for (size_t j = 0; j < Nlinks; j++)
                node[node[i].outLinks[j].first].size += beta *
                        node[i].outLinks[j].second * size_tmp[i];
        }

        // Normalize
        sum = 0.0;
        for (igraph_integer_t i = 0; i < Nnode; i++) {
            sum += node[i].size;
        }
        sqdiff_old = sqdiff;
        sqdiff = 0.0;
        for (igraph_integer_t i = 0; i < Nnode; i++) {
            node[i].size /= sum;
            sqdiff += fabs(node[i].size - size_tmp[i]);
            size_tmp[i] = node[i].size;
        }
        Niterations++;

        if (sqdiff == sqdiff_old) {
            alpha += 1.0e-10;
            beta = 1.0 - alpha;
        }

    } while ((Niterations < 200) && (sqdiff > 1.0e-15 || Niterations < 50));

    danglingSize = 0.0;
    for (igraph_integer_t i = 0; i < Ndanglings; i++) {
        danglingSize += size_tmp[danglings[i]];
    }
    // cout << "done! (the error is " << sqdiff << " after " << Niterations
    //      << " iterations)" << endl;
}


/* Compute the codeLength of the given network
 * note: (in **node, one node == one module)
 */
void FlowGraph::calibrate() noexcept {
    exit_log_exit = 0.0;
    exitFlow = 0.0;
    size_log_size = 0.0;

    for (igraph_integer_t i = 0; i < Nnode; i++) { // For each module
        // own node/module codebook
        size_log_size += plogp(node[i].exit + node[i].size);

        // use of index codebook
        exitFlow      += node[i].exit;
        exit_log_exit += plogp(node[i].exit);
    }

    exit = plogp(exitFlow);

    codeLength = exit - 2.0 * exit_log_exit + size_log_size -
                 nodeSize_log_nodeSize;
}


/* Restore the data from the given FlowGraph object
 */
void FlowGraph::back_to(const FlowGraph &fgraph) {
    // delete current nodes and copy original ones
    Nnode = fgraph.Nnode;
    node = fgraph.node;

    // restore atributs
    alpha = fgraph.alpha ;
    beta  = fgraph.beta ;

    exit = fgraph.exit;
    exitFlow = fgraph.exitFlow;
    exit_log_exit = fgraph.exit_log_exit;
    size_log_size = fgraph.size_log_size ;
    nodeSize_log_nodeSize = fgraph.nodeSize_log_nodeSize;

    codeLength = fgraph.codeLength;
}
