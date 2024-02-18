/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

/* The original version of this file was written by Jörg Reichardt
   This file was modified by Vincent Traag
   The original copyright notice follows here */

/***************************************************************************
                          pottsmodel.cpp  -  description
                             -------------------
    begin                : Fri May 28 2004
    copyright            : (C) 2004 by
    email                :
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "pottsmodel_2.h"

#include "igraph_random.h"
#include "core/interruption.h"

#include <cstring>
#include <cmath>

using namespace std;

//#################################################################################################
PottsModel::PottsModel(network *n, igraph_integer_t qvalue, int m) :
        net(n), q(qvalue), operation_mode(m), Qmatrix(qvalue+1)
{
    DLList_Iter<NNode*> iter;
    const NNode *n_cur;
    igraph_integer_t *i_ptr;
    //needed in calculating modularity
    Qa     = new double[q + 1];
    //weights for each spin state needed in Monte Carlo process
    weights = new double[q + 1];
    //bookkeeping of occupation numbers of spin states or the number of links in community
    color_field = new double[q + 1];
    neighbours = new double[q + 1];

    num_of_nodes = net->node_list.Size();
    num_of_links = net->link_list.Size();

    n_cur = iter.First(&net->node_list);
    while (!iter.End()) {
        if (k_max < n_cur->Get_Degree()) {
            k_max = n_cur->Get_Degree();
        }
        i_ptr = new igraph_integer_t;
        *i_ptr = 0;
        new_spins.Push(i_ptr);
        i_ptr = new igraph_integer_t;
        *i_ptr = 0;
        previous_spins.Push(i_ptr);
        n_cur = iter.Next();
    }
}
//#######################################################
//Destructor of PottsModel
//########################################################
PottsModel::~PottsModel() {
    /* The DLItem destructor does not delete its item currently,
       because of some bad design. As a workaround, we delete them here
       by hand */
    new_spins.delete_items();
    previous_spins.delete_items();
    delete [] Qa;
    delete [] weights;
    delete [] color_field;
    delete [] neighbours;
}
//#####################################################
//Assing an initial random configuration of spins to nodes
//if called with negative argument or the spin used as argument
//when called with positve one.
//This may be handy, if you want to warm up the network.
//####################################################
igraph_integer_t PottsModel::assign_initial_conf(igraph_integer_t spin) {
    igraph_integer_t s;
    DLList_Iter<NNode*> iter;
    DLList_Iter<NLink*> l_iter;
    NNode *n_cur;
    const NLink *l_cur;
    double sum_weight;

    // initialize colorfield
    for (igraph_integer_t i = 0; i <= q; i++) {
        color_field[i] = 0.0;
    }
    //
    total_degree_sum = 0.0;
    n_cur = iter.First(&net->node_list);
    while (!iter.End()) {
        if (spin < 0) {
            s = RNG_INTEGER(1, q);
        } else {
            s = spin;
        }
        n_cur->Set_ClusterIndex(s);
        l_cur = l_iter.First(n_cur->Get_Links());
        sum_weight = 0;
        while (!l_iter.End()) {
            sum_weight += l_cur->Get_Weight(); //weight should be one, in case we are not using it.
            l_cur = l_iter.Next();
        }
        // we set the sum of the weights or the degree as the weight of the node, this way
        // we do not have to calculate it again.
        n_cur->Set_Weight(sum_weight);

        // in case we want all links to be contribute equally - parameter gamm=fixed
        if (operation_mode == 0) {
            color_field[s]++;
        } else {
            color_field[s] += sum_weight;
        }
        // or in case we want to use a weight of each link that is proportional to k_i\times k_j
        total_degree_sum += sum_weight;
        n_cur = iter.Next();
    }

    return net->node_list.Size();
}

//#####################################################################
// Q denotes the modularity of the network
// This function calculates it initially
// In the event of a spin changing its state, it only needs updating
// Note that Qmatrix and Qa are only counting! The normalization
// by num_of_links is done later
//####################################################################
double PottsModel::initialize_Qmatrix() {
    DLList_Iter<NLink*> l_iter;
    NLink *l_cur;
    igraph_integer_t i, j;
    //initialize with zeros
    num_of_links = net->link_list.Size();
    for (i = 0; i <= q; i++) {
        Qa[i] = 0.0;
        for (j = i; j <= q; j++) {
            Qmatrix[i][j] = 0.0;
            Qmatrix[j][i] = 0.0;
        }
    }
    //go over all links and make corresponding entries in Q matrix
    //An edge connecting state i wiht state j will get an entry in Qij and Qji
    l_cur = l_iter.First(&net->link_list);
    while (!l_iter.End()) {
        i = l_cur->Get_Start()->Get_ClusterIndex();
        j = l_cur->Get_End()->Get_ClusterIndex();
        Qmatrix[i][j] += l_cur->Get_Weight();
        Qmatrix[j][i] += l_cur->Get_Weight();

        l_cur = l_iter.Next();
    }
    //Finally, calculate sum over rows and keep in Qa
    for (i = 0; i <= q; i++) {
        for (j = 0; j <= q; j++) {
            Qa[i] += Qmatrix[i][j];
        }
    }
    return calculate_Q();
}
//####################################################################
// This function does the actual calculation of Q from the matrix
// The normalization by num_of_links is done here
//####################################################################
double PottsModel::calculate_Q() {
    double Q = 0.0;
    for (igraph_integer_t i = 0; i <= q; i++) {
        Q += Qmatrix[i][i] - Qa[i] * Qa[i] / double(2.0 * net->sum_weights);
    }
    Q /= double(2.0 * net->sum_weights);
    return Q;
}

//##########################################################################
// We would like to start from a temperature with at least 95 of all proposed
// spin changes accepted in 50 sweeps over the network
// The function returns the Temperature found
//#########################################################################
double PottsModel::FindStartTemp(double gamma, double prob, double ts) {
    double kT;
    kT = ts;
    //assing random initial condition
    assign_initial_conf(-1);
    //initialize Modularity matrix, from now on, it will be updated at every spin change
    initialize_Qmatrix();
    // the factor 1-1/q is important, since even, at infinite temperature,
    // only 1-1/q of all spins do change their state, since a randomly chooses new
    // state is with prob. 1/q the old state.
    while (acceptance < (1.0 - 1.0 / double(q)) * 0.95) { //want 95% acceptance
        kT = kT * 1.1;
        HeatBathParallelLookup(gamma, prob, kT, 50);
    }
    kT *= 1.1; // just to be sure...
    return kT;
}

//##############################################################
//This function does a parallel update at zero T
//Hence, it is really fast on easy problems
//max sweeps is the maximum number of sweeps it should perform,
//if it does not converge earlier
//##############################################################
igraph_integer_t PottsModel::HeatBathParallelLookupZeroTemp(double gamma, double prob, unsigned int max_sweeps) {
    DLList_Iter<NNode *> net_iter;
    DLList_Iter<NLink*> l_iter;
    DLList_Iter<igraph_integer_t*> i_iter, i_iter2;
    NNode *node, *n_cur;
    NLink *l_cur;
    unsigned int sweep;
    igraph_integer_t *SPIN, *P_SPIN, old_spin, new_spin, spin_opt;
    igraph_integer_t changes;
    double h, delta = 0, deltaE, deltaEmin, w, degree;
    bool cyclic = false;

    sweep = 0;
    changes = 1;
    while (sweep < max_sweeps && changes) {
        cyclic = true;
        sweep++;
        changes = 0;
        //Loop over all nodes
        node = net_iter.First(&net->node_list);
        SPIN = i_iter.First(&new_spins);
        while (!net_iter.End()) {
            // How many neighbors of each type?
            // set them all zero
            for (igraph_integer_t i = 0; i <= q; i++) {
                neighbours[i] = 0;
            }
            degree = node->Get_Weight();
            //Loop over all links (=neighbours)
            l_cur = l_iter.First(node->Get_Links());
            while (!l_iter.End()) {
                w = l_cur->Get_Weight();
                if (node == l_cur->Get_Start()) {
                    n_cur = l_cur->Get_End();
                } else {
                    n_cur = l_cur->Get_Start();
                }
                neighbours[n_cur->Get_ClusterIndex()] += w;
                l_cur = l_iter.Next();
            }
            //Search optimal Spin
            old_spin = node->Get_ClusterIndex();
            switch (operation_mode) {
            case 0: {
                delta = 1.0;
                break;
            }
            case 1: { //newman modularity
                prob = degree / total_degree_sum;
                delta = degree;
                break;
            }
            default:
                IGRAPH_FATAL("Must not reach here.");
            }


            spin_opt = old_spin;
            deltaEmin = 0.0;
            for (igraph_integer_t spin = 1; spin <= q; spin++) { // all possible spin states
                if (spin != old_spin) {
                    h = color_field[spin] + delta - color_field[old_spin];
                    deltaE = double(neighbours[old_spin] - neighbours[spin]) + gamma * prob * double(h);
                    if (deltaE < deltaEmin) {
                        spin_opt = spin;
                        deltaEmin = deltaE;
                    }
                }
            } // for spin

            //Put optimal spin on list for later update
            *SPIN = spin_opt;
            node = net_iter.Next();
            SPIN = i_iter.Next();
        } // while !net_iter.End()

        //-------------------------------
        //Now set all spins to new values
        node = net_iter.First(&net->node_list);
        SPIN = i_iter.First(&new_spins);
        P_SPIN = i_iter2.First(&previous_spins);
        while (!net_iter.End()) {
            old_spin = node->Get_ClusterIndex();
            new_spin = *SPIN;
            if (new_spin != old_spin) { // Do we really have a change??
                changes++;
                node->Set_ClusterIndex(new_spin);
                //this is important!!
                //In Parallel update, there occur cyclic attractors of size two
                //which then make the program run for ever
                if (new_spin != *P_SPIN) {
                    cyclic = false;
                }
                *P_SPIN = old_spin;
                color_field[old_spin]--;
                color_field[new_spin]++;

                //Qmatrix update
                //iteration over all neighbours
                l_cur = l_iter.First(node->Get_Links());
                while (!l_iter.End()) {
                    w = l_cur->Get_Weight();
                    if (node == l_cur->Get_Start()) {
                        n_cur = l_cur->Get_End();
                    } else {
                        n_cur = l_cur->Get_Start();
                    }
                    Qmatrix[old_spin][n_cur->Get_ClusterIndex()] -= w;
                    Qmatrix[new_spin][n_cur->Get_ClusterIndex()] += w;
                    Qmatrix[n_cur->Get_ClusterIndex()][old_spin] -= w;
                    Qmatrix[n_cur->Get_ClusterIndex()][new_spin] += w;
                    Qa[old_spin] -= w;
                    Qa[new_spin] += w;
                    l_cur = l_iter.Next();
                }  // while l_iter
            }
            node = net_iter.Next();
            SPIN = i_iter.Next();
            P_SPIN = i_iter2.Next();
        } // while (!net_iter.End())
    }  // while markov

    // In case of a cyclic attractor, we want to interrupt
    if (cyclic)  {
        acceptance = 0.0;
        return 0;
    } else {
        acceptance = double(changes) / double(num_of_nodes);
        return changes;
    }
}
//###################################################################################
//The same function as before, but rather than parallel update, it pics the nodes to update
//randomly
//###################################################################################
double PottsModel::HeatBathLookupZeroTemp(double gamma, double prob, unsigned int max_sweeps) {
    DLList_Iter<NLink*> l_iter;
    NNode *node, *n_cur;
    NLink *l_cur;
    igraph_integer_t new_spin, spin_opt, old_spin;
    unsigned int sweep;
    igraph_integer_t r;
    igraph_integer_t changes;
    double delta = 0, h, deltaE, deltaEmin, w, degree;

    sweep = 0;
    changes = 0;
    while (sweep < max_sweeps) {
        sweep++;
        //ueber alle Knoten im Netz
        for (igraph_integer_t n = 0; n < num_of_nodes; n++) {
            r = RNG_INTEGER(0, num_of_nodes - 1);
            node = net->node_list.Get(r);
            // Wir zaehlen, wieviele Nachbarn von jedem spin vorhanden sind
            // erst mal alles Null setzen
            for (igraph_integer_t i = 0; i <= q; i++) {
                neighbours[i] = 0;
            }
            degree = node->Get_Weight();
            //Loop over all links (=neighbours)
            l_cur = l_iter.First(node->Get_Links());
            while (!l_iter.End()) {
                w = l_cur->Get_Weight();
                if (node == l_cur->Get_Start()) {
                    n_cur = l_cur->Get_End();
                } else {
                    n_cur = l_cur->Get_Start();
                }
                neighbours[n_cur->Get_ClusterIndex()] += w;
                l_cur = l_iter.Next();
            }
            //Search optimal Spin
            old_spin = node->Get_ClusterIndex();
            switch (operation_mode) {
            case 0: {
                delta = 1.0;
                break;
            }
            case 1: { //newman modularity
                prob = degree / total_degree_sum;
                delta = degree;
                break;
            }
            default:
                IGRAPH_FATAL("Must not reach here.");
            }


            spin_opt = old_spin;
            deltaEmin = 0.0;
            for (igraph_integer_t spin = 1; spin <= q; spin++) { // alle moeglichen Spins
                if (spin != old_spin) {
                    h = color_field[spin] + delta - color_field[old_spin];
                    deltaE = double(neighbours[old_spin] - neighbours[spin]) + gamma * prob * double(h);
                    if (deltaE < deltaEmin) {
                        spin_opt = spin;
                        deltaEmin = deltaE;
                    }
                }
            } // for spin

            //-------------------------------
            //Now update the spins
            new_spin = spin_opt;
            if (new_spin != old_spin) { // Did we really change something??
                changes++;
                node->Set_ClusterIndex(new_spin);
                color_field[old_spin] -= delta;
                color_field[new_spin] += delta;

                //Qmatrix update
                //iteration over all neighbours
                l_cur = l_iter.First(node->Get_Links());
                while (!l_iter.End()) {
                    w = l_cur->Get_Weight();
                    if (node == l_cur->Get_Start()) {
                        n_cur = l_cur->Get_End();
                    } else {
                        n_cur = l_cur->Get_Start();
                    }
                    Qmatrix[old_spin][n_cur->Get_ClusterIndex()] -= w;
                    Qmatrix[new_spin][n_cur->Get_ClusterIndex()] += w;
                    Qmatrix[n_cur->Get_ClusterIndex()][old_spin] -= w;
                    Qmatrix[n_cur->Get_ClusterIndex()][new_spin] += w;
                    Qa[old_spin] -= w;
                    Qa[new_spin] += w;
                    l_cur = l_iter.Next();
                }  // while l_iter
            }
        } // for n
    }  // while markov

    acceptance = double(changes) / double(num_of_nodes) / double(sweep);
    return acceptance;
}
//#####################################################################################
//This function performs a parallel update at Terperature T
//#####################################################################################
igraph_integer_t PottsModel::HeatBathParallelLookup(double gamma, double prob, double kT, unsigned int max_sweeps) {
    DLList_Iter<NNode*> net_iter;
    DLList_Iter<NLink*> l_iter;
    DLList_Iter<igraph_integer_t*> i_iter, i_iter2;
    NNode *node, *n_cur;
    NLink *l_cur;
    igraph_integer_t new_spin, spin_opt, old_spin;
    igraph_integer_t *SPIN, *P_SPIN;
    unsigned int sweep;
    igraph_integer_t max_q;
    igraph_integer_t changes;
    double h, delta = 0, norm, r, beta, minweight, prefac = 0, w, degree;
    bool cyclic = false/*, found*/;
    igraph_integer_t number_of_nodes;

    sweep = 0;
    changes = 1;
    number_of_nodes = net->node_list.Size();
    while (sweep < max_sweeps && changes) {
        cyclic = true;
        sweep++;
        changes = 0;
        //Loop over all nodes
        node = net_iter.First(&net->node_list);
        SPIN = i_iter.First(&new_spins);
        while (!net_iter.End()) {
            // Initialize neighbours and weights
            for (igraph_integer_t i = 0; i <= q; i++) {
                neighbours[i] = 0;
                weights[i] = 0;
            }
            norm = 0.0;
            degree = node->Get_Weight();
            //Loop over all links (=neighbours)
            l_cur = l_iter.First(node->Get_Links());
            while (!l_iter.End()) {
                w = l_cur->Get_Weight();
                if (node == l_cur->Get_Start()) {
                    n_cur = l_cur->Get_End();
                } else {
                    n_cur = l_cur->Get_Start();
                }
                neighbours[n_cur->Get_ClusterIndex()] += w;
                l_cur = l_iter.Next();
            }
            //Search optimal Spin
            old_spin = node->Get_ClusterIndex();
            switch (operation_mode) {
            case 0: {
                prefac = 1.0;
                delta = 1.0;
                break;
            }
            case 1: { //newman modularity
                prefac = 1.0;
                prob = degree / total_degree_sum;
                delta = degree;
                break;
            }
            default:
                IGRAPH_FATAL("Must not reach here.");
            }
            spin_opt = old_spin;
            beta = 1.0 / kT * prefac;
            minweight = 0.0;
            weights[old_spin] = 0.0;
            for (igraph_integer_t spin = 1; spin <= q; spin++) { // loop over all possible new spins
                if (spin != old_spin) { // only if we have a different than old spin!
                    h = color_field[spin] + delta - color_field[old_spin];
                    weights[spin] = double(neighbours[old_spin] - neighbours[spin]) + gamma * prob * double(h);
                    if (weights[spin] < minweight) {
                        minweight = weights[spin];
                    }
                }
            }   // for spin
            for (igraph_integer_t spin = 1; spin <= q; spin++) { // loop over all possibe spins
                weights[spin] -= minweight;       // subtract minweight
                // to avoid numerical problems with large exponents
                weights[spin] = exp(-beta * weights[spin]);
                norm += weights[spin];
            }   // for spin

            //now choose a new spin
            r = RNG_UNIF(0, norm);
            new_spin = 1;
            while (new_spin <= q) {
                if (r <= weights[new_spin]) {
                    spin_opt = new_spin;
                    break;
                } else {
                    r -= weights[new_spin];
                }
                new_spin++;
            }
            //Put new spin on list
            *SPIN = spin_opt;

            node = net_iter.Next();
            SPIN = i_iter.Next();
        } // while !net_iter.End()

        //-------------------------------
        //now update all spins
        node = net_iter.First(&net->node_list);
        SPIN = i_iter.First(&new_spins);
        P_SPIN = i_iter2.First(&previous_spins);
        while (!net_iter.End()) {
            old_spin = node->Get_ClusterIndex();
            new_spin = *SPIN;
            if (new_spin != old_spin) { // Did we really change something??
                changes++;
                node->Set_ClusterIndex(new_spin);
                if (new_spin != *P_SPIN) {
                    cyclic = false;
                }
                *P_SPIN = old_spin;
                color_field[old_spin] -= delta;
                color_field[new_spin] += delta;

                //Qmatrix update
                //iteration over all neighbours
                l_cur = l_iter.First(node->Get_Links());
                while (!l_iter.End()) {
                    w = l_cur->Get_Weight();
                    if (node == l_cur->Get_Start()) {
                        n_cur = l_cur->Get_End();
                    } else {
                        n_cur = l_cur->Get_Start();
                    }
                    Qmatrix[old_spin][n_cur->Get_ClusterIndex()] -= w;
                    Qmatrix[new_spin][n_cur->Get_ClusterIndex()] += w;
                    Qmatrix[n_cur->Get_ClusterIndex()][old_spin] -= w;
                    Qmatrix[n_cur->Get_ClusterIndex()][new_spin] += w;
                    Qa[old_spin] -= w;
                    Qa[new_spin] += w;
                    l_cur = l_iter.Next();
                }  // while l_iter
            }
            node = net_iter.Next();
            SPIN = i_iter.Next();
            P_SPIN = i_iter2.Next();
        } // while (!net_iter.End())

    }  // while markov
    max_q = 0;
    for (igraph_integer_t i = 1; i <= q; i++) if (color_field[i] > max_q) {
            max_q = igraph_integer_t(color_field[i]);
        }

    //again, we would not like to end up in cyclic attractors
    if (cyclic && changes)  {
        acceptance = double(changes) / double(number_of_nodes);
        return 0;
    } else {
        acceptance = double(changes) / double(number_of_nodes);
        return changes;
    }
}
//##############################################################
// This is the function generally used for optimisation,
// as the parallel update has its flaws, due to the cyclic attractors
//##############################################################
double PottsModel::HeatBathLookup(double gamma, double prob, double kT, unsigned int max_sweeps) {
    DLList_Iter<NLink*> l_iter;
    NNode *node, *n_cur;
    NLink *l_cur;
    igraph_integer_t new_spin, spin_opt, old_spin;
    unsigned int sweep;
    igraph_integer_t max_q;
    igraph_integer_t rn;
    igraph_integer_t changes;
    double degree, w, delta = 0, h;
    double norm, r, beta, minweight, prefac = 0;
    igraph_integer_t number_of_nodes;
    sweep = 0;
    changes = 0;
    number_of_nodes = net->node_list.Size();
    while (sweep < max_sweeps) {
        sweep++;
        //loop over all nodes in network
        for (igraph_integer_t n = 0; n < number_of_nodes; n++) {
            rn = RNG_INTEGER(0, number_of_nodes - 1);

            node = net->node_list.Get(rn);
            // initialize the neighbours and the weights
            for (igraph_integer_t i = 0; i <= q; i++) {
                neighbours[i] = 0.0;
                weights[i] = 0.0;
            }
            norm = 0.0;
            degree = node->Get_Weight();
            //Loop over all links (=neighbours)
            l_cur = l_iter.First(node->Get_Links());
            while (!l_iter.End()) {
                w = l_cur->Get_Weight();
                if (node == l_cur->Get_Start()) {
                    n_cur = l_cur->Get_End();
                } else {
                    n_cur = l_cur->Get_Start();
                }
                neighbours[n_cur->Get_ClusterIndex()] += w;
                l_cur = l_iter.Next();
            }

            //Look for optimal spin

            old_spin = node->Get_ClusterIndex();
            switch (operation_mode) {
            case 0: {
                prefac = 1.0;
                delta = 1.0;
                break;
            }
            case 1:  {//newman modularity
                prefac = 1.0;
                prob = degree / total_degree_sum;
                delta = degree;
                break;
            }
            default:
                IGRAPH_FATAL("Must not reach here.");
            }
            spin_opt = old_spin;
            beta = 1.0 / kT * prefac;
            minweight = 0.0;
            weights[old_spin] = 0.0;
            for (igraph_integer_t spin = 1; spin <= q; spin++) { // all possible new spins
                if (spin != old_spin) { // except the old one!
                    h = color_field[spin] - (color_field[old_spin] - delta);
                    weights[spin] = neighbours[old_spin] - neighbours[spin] + gamma * prob * h;
                    if (weights[spin] < minweight) {
                        minweight = weights[spin];
                    }
                }
            }   // for spin
            for (igraph_integer_t spin = 1; spin <= q; spin++) { // all possible new spins
                weights[spin] -= minweight;       // subtract minweigt
                // for numerical stability
                weights[spin] = exp(-beta * weights[spin]);
                norm += weights[spin];
            }   // for spin


            //choose a new spin
            r = RNG_UNIF(0, norm);
            new_spin = 1;
            while (new_spin <= q) {
                if (r <= weights[new_spin]) {
                    spin_opt = new_spin;
                    break;
                } else {
                    r -= weights[new_spin];
                }
                new_spin++;
            }
            //-------------------------------
            //now set the new spin
            new_spin = spin_opt;
            if (new_spin != old_spin) { // Did we really change something??
                changes++;
                node->Set_ClusterIndex(new_spin);
                color_field[old_spin] -= delta;
                color_field[new_spin] += delta;

                //Qmatrix update
                //iteration over all neighbours
                l_cur = l_iter.First(node->Get_Links());
                while (!l_iter.End()) {
                    w = l_cur->Get_Weight();
                    if (node == l_cur->Get_Start()) {
                        n_cur = l_cur->Get_End();
                    } else {
                        n_cur = l_cur->Get_Start();
                    }
                    Qmatrix[old_spin][n_cur->Get_ClusterIndex()] -= w;
                    Qmatrix[new_spin][n_cur->Get_ClusterIndex()] += w;
                    Qmatrix[n_cur->Get_ClusterIndex()][old_spin] -= w;
                    Qmatrix[n_cur->Get_ClusterIndex()][new_spin] += w;
                    Qa[old_spin] -= w;
                    Qa[new_spin] += w;
                    l_cur = l_iter.Next();
                }  // while l_iter
            }
        } // for n
    }  // while markov
    max_q = 0;

    for (igraph_integer_t i = 1; i <= q; i++) if (color_field[i] > max_q) {
            max_q = igraph_integer_t(color_field[i] + 0.5);
        }

    acceptance = double(changes) / double(number_of_nodes) / double(sweep);
    return acceptance;
}

//###############################################################################################
//# Here we try to minimize the affinity to the rest of the network
//###############################################################################################
double PottsModel::FindCommunityFromStart(
        double gamma,
        const char *nodename,
        igraph_vector_int_t *result,
        igraph_real_t *cohesion,
        igraph_real_t *adhesion,
        igraph_integer_t *my_inner_links,
        igraph_integer_t *my_outer_links) const {
    DLList_Iter<NNode*> iter, iter2;
    DLList_Iter<NLink*> l_iter;
    DLList<NNode*> to_do;
    DLList<NNode*> community;
    NNode *start_node = nullptr, *n_cur, *neighbor, *max_aff_node, *node;
    NLink *l_cur;
    bool found = false, add = false, remove = false;
    double degree, delta_aff_add, delta_aff_rem, max_delta_aff, Ks = 0.0, Kr = 0, kis, kir, w;
    igraph_integer_t community_marker = 5;
    igraph_integer_t to_do_marker = 10;
    double inner_links = 0, outer_links = 0, aff_r, aff_s;

    // find the node in the network
    n_cur = iter.First(&net->node_list);
    while (!found && !iter.End()) {
        if (0 == strcmp(n_cur->Get_Name(), nodename)) {
            start_node = n_cur;
            found = true;
            community.Push(start_node);
            start_node->Set_Marker(community_marker);
            Ks = start_node->Get_Weight();
            Kr = total_degree_sum - start_node->Get_Weight();
        }
        n_cur = iter.Next();
    }
    if (!found) {
        return -1;
    }
    //#############################
    // initialize the to_do list and community with the neighbours of start node
    //#############################
    neighbor = iter.First(start_node->Get_Neighbours());
    while (!iter.End()) {
        community.Push(neighbor);
        neighbor->Set_Marker(community_marker);
        Ks += neighbor->Get_Weight();
        Kr -= neighbor->Get_Weight();
        neighbor = iter.Next();
    }
    node = iter.First(&community);
    while (!iter.End()) {
        //now add at the second neighbors to the to_do list
        neighbor = iter2.First(node->Get_Neighbours());
        while (!iter2.End()) {
            if (neighbor->Get_Marker() != community_marker && neighbor->Get_Marker() != to_do_marker) {
                to_do.Push(neighbor);
                neighbor->Set_Marker(to_do_marker);
            }
            neighbor = iter2.Next();
        }
        node = iter.Next();
    }

    //#############
    //repeat, as long as we are still adding nodes to the communtiy
    //#############
    add = true;
    remove = true;
    while (add || remove) {
        //#############################
        //calculate the affinity changes of all nodes for adding every node in the to_do list to the community
        //##############################

        IGRAPH_ALLOW_INTERRUPTION(); /* This is not clean.... */

        max_delta_aff = 0.0;
        max_aff_node = nullptr;
        add = false;
        node = iter.First(&to_do);
        while (!iter.End()) {
            //printf("Checking Links of %s\n",node->Get_Name());
            degree = node->Get_Weight();
            kis = 0.0;
            kir = 0.0;
            // For every of the neighbors, check, count the links to the community
            l_cur = l_iter.First(node->Get_Links());
            while (!l_iter.End()) {
                w = l_cur->Get_Weight();
                if (node == l_cur->Get_Start()) {
                    n_cur = l_cur->Get_End();
                } else {
                    n_cur = l_cur->Get_Start();
                }
                if (n_cur->Get_Marker() == community_marker) {
                    kis += w; //the weight/number of links to the community
                } else {
                    kir += w; //the weight/number of links to the rest of the network
                }
                l_cur = l_iter.Next();
            }
            aff_r = kir - gamma / total_degree_sum * (Kr - degree) * degree;
            aff_s = kis - gamma / total_degree_sum * Ks * degree;
            delta_aff_add = aff_r - aff_s;
            if (delta_aff_add <= max_delta_aff) {
                max_delta_aff = delta_aff_add;
                max_aff_node = node;
                add = true;
            }
            node = iter.Next();
        }
        //################
        //calculate the affinity changes for removing every single node from the community
        //################
        inner_links = 0;
        outer_links = 0;
        remove = false;
        node = iter.First(&community);
        while (!iter.End()) {
            //printf("Checking Links of %s\n",node->Get_Name());
            degree = node->Get_Weight();
            kis = 0.0;
            kir = 0.0;
            // For every of the neighbors, check, count the links to the community
            l_cur = l_iter.First(node->Get_Links());
            while (!l_iter.End()) {
                w = l_cur->Get_Weight();
                if (node == l_cur->Get_Start()) {
                    n_cur = l_cur->Get_End();
                } else {
                    n_cur = l_cur->Get_Start();
                }
                if (n_cur->Get_Marker() == community_marker) {
                    kis += w;
                    inner_links += w; //summing all w gives twice the number of inner links(weights)
                } else {
                    kir += w;
                    outer_links += w;
                }
                l_cur = l_iter.Next();
            }
            aff_r = kir - gamma / total_degree_sum * Kr * degree;
            aff_s = kis - gamma / total_degree_sum * (Ks - degree) * degree;
            delta_aff_rem = aff_s - aff_r;
            // we should not remove the nodes, we have just added
            if (delta_aff_rem < max_delta_aff) {
                max_delta_aff = delta_aff_rem ;
                max_aff_node = node;
                remove = true;
                add = false;
            }
            node = iter.Next();
        }
        inner_links = inner_links * 0.5;
        //################
        // Now check, whether we want to remove or add a node
        //################
        if (add) {
            //################
            //add the node of maximum affinity to the community
            //###############
            community.Push(max_aff_node);
            max_aff_node->Set_Marker(community_marker);
            //delete node from to_do
            to_do.fDelete(max_aff_node);
            //update the sum of degrees in the community
            Ks += max_aff_node->Get_Weight();
            Kr -= max_aff_node->Get_Weight();
            //now add all neighbors of this node, that are not already
            //in the to_do list or in the community
            neighbor = iter.First(max_aff_node->Get_Neighbours());
            while (!iter.End()) {
                if (neighbor->Get_Marker() != community_marker && neighbor->Get_Marker() != to_do_marker) {
                    to_do.Push(neighbor);
                    neighbor->Set_Marker(to_do_marker);
                    //printf("Adding node %s to to_do list.\n",neighbor->Get_Name());
                }
                neighbor = iter.Next();
            }
        }
        if (remove) {
            //################
            //remove those with negative affinities
            //################
            community.fDelete(max_aff_node);
            max_aff_node->Set_Marker(to_do_marker);
            //update the sum of degrees in the community
            Ks -= max_aff_node->Get_Weight();
            Kr += max_aff_node->Get_Weight();
            //add the node to to_do again
            to_do.Push(max_aff_node);
        }
        IGRAPH_ALLOW_INTERRUPTION(); /* This is not clean.... */
    }
    //###################
    //write the node in the community to a file
    //###################
    if (cohesion) {
        *cohesion = inner_links - gamma / total_degree_sum * Ks * Ks * 0.5;
    }
    if (adhesion) {
        *adhesion = outer_links - gamma / total_degree_sum * Ks * Kr;
    }
    if (my_inner_links) {
        *my_inner_links = inner_links;
    }
    if (my_outer_links) {
        *my_outer_links = outer_links;
    }
    if (result) {
        node = iter.First(&community);
        igraph_vector_int_clear(result);
        while (!iter.End()) {
            IGRAPH_CHECK(igraph_vector_int_push_back(result, node->Get_Index()));
            node = iter.Next();
        }
    }
    igraph_integer_t size = community.Size();
    return size;
}

//################################################################################################
// this Function writes the clusters to disk
//################################################################################################
igraph_integer_t PottsModel::WriteClusters(igraph_real_t *modularity,
                               igraph_real_t *temperature,
                               igraph_vector_int_t *csize,
                               igraph_vector_int_t *membership,
                               double kT, double gamma) const {
    const NNode *n_cur, *n_cur2;
    DLList_Iter<NNode*> iter, iter2;
    HugeArray<int> inner_links;
    HugeArray<int> outer_links;
    HugeArray<int> nodes;

    if (temperature) {
        *temperature = kT;
    }

    if (csize || membership || modularity) {
        // TODO: count the number of clusters
        for (igraph_integer_t spin = 1; spin <= q; spin++) {
            inner_links[spin] = 0;
            outer_links[spin] = 0;
            nodes[spin] = 0;
            n_cur = iter.First(&net->node_list);
            while (!iter.End()) {
                if (n_cur->Get_ClusterIndex() == spin) {
                    nodes[spin]++;
                    n_cur2 = iter2.First(n_cur->Get_Neighbours());
                    while (!iter2.End()) {
                        if (n_cur2->Get_ClusterIndex() == spin) {
                            inner_links[spin]++;
                        } else {
                            outer_links[spin]++;
                        }
                        n_cur2 = iter2.Next();
                    }
                }
                n_cur = iter.Next();
            }
        }
    }
    if (modularity) {
        *modularity = 0.0;
        for (igraph_integer_t spin = 1; spin <= q; spin++) {
            if (nodes[spin] > 0) {
                double t1 = inner_links[spin] / net->sum_weights / 2.0;
                double t2 = (inner_links[spin] + outer_links[spin]) /
                            net->sum_weights / 2.0;
                *modularity += t1;
                *modularity -= gamma * t2 * t2;
            }
        }
    }
    if (csize) {
        igraph_vector_int_clear(csize);
        for (igraph_integer_t spin = 1; spin <= q; spin++) {
            if (nodes[spin] > 0) {
                inner_links[spin] /= 2;
                IGRAPH_CHECK(igraph_vector_int_push_back(csize, nodes[spin]));
            }
        }
    }

    //die Elemente der Cluster
    if (membership) {
        igraph_integer_t no = -1;
        IGRAPH_CHECK(igraph_vector_int_resize(membership, num_of_nodes));
        for (igraph_integer_t spin = 1; spin <= q; spin++) {
            if (nodes[spin] > 0) {
                no++;
            }
            n_cur = iter.First(&net->node_list);
            while (!iter.End()) {
                if (n_cur->Get_ClusterIndex() == spin) {
                    VECTOR(*membership)[ n_cur->Get_Index() ] = no;
                }
                n_cur = iter.Next();
            }
        }
    }

    return num_of_nodes;
}

//#################################################################################################
PottsModelN::PottsModelN(network *n, igraph_integer_t num_communities, bool directed) :
    net(n), q(num_communities), num_nodes(net->node_list.Size()), is_directed(directed)
{ }
//#######################################################
//Destructor of PottsModel
//########################################################
PottsModelN::~PottsModelN() {
    delete [] degree_pos_in;
    delete [] degree_neg_in;
    delete [] degree_pos_out;
    delete [] degree_neg_out;

    delete [] degree_community_pos_in;
    delete [] degree_community_neg_in;
    delete [] degree_community_pos_out;
    delete [] degree_community_neg_out;

    delete [] weights;
    delete [] neighbours;
    delete [] csize;

    delete [] spin;
}

void PottsModelN::assign_initial_conf(bool init_spins) {
    igraph_integer_t s;
    DLList_Iter<NLink*> l_iter;
    const NNode *n_cur;
    const NLink *l_cur;

    if (init_spins) {
        // Free the arrays before (re-)allocating them
        // These arrays are initialized to NULL, so it is safe to delete even before allocation
        delete [] degree_pos_in;
        delete [] degree_neg_in;
        delete [] degree_pos_out;
        delete [] degree_neg_out;

        delete [] spin;

        //Bookkeeping of the various degrees (positive/negative) and (in/out)
        degree_pos_in   = new double[num_nodes]; //Postive indegree of the nodes (or sum of weights)
        degree_neg_in   = new double[num_nodes]; //Negative indegree of the nodes (or sum of weights)
        degree_pos_out  = new double[num_nodes]; //Postive outdegree of the nodes (or sum of weights)
        degree_neg_out  = new double[num_nodes]; //Negative outdegree of the nodes (or sum of weights)

        spin            = new igraph_integer_t[num_nodes]; //The spin state of each node
    }

    if (is_init) {
        delete [] degree_community_pos_in;
        delete [] degree_community_neg_in;
        delete [] degree_community_pos_out;
        delete [] degree_community_neg_out;

        delete [] weights;
        delete [] neighbours;
        delete [] csize;
    }

    is_init = true;

    //Bookkeep of occupation numbers of spin states or the number of links in community...
    degree_community_pos_in     = new double[q + 1]; //Positive sum of indegree for communities
    degree_community_neg_in     = new double[q + 1]; //Negative sum of indegree for communities
    degree_community_pos_out    = new double[q + 1]; //Positive sum of outegree for communities
    degree_community_neg_out    = new double[q + 1]; //Negative sum of outdegree for communities

    //...and of weights and neighbours for in the HeathBathLookup
    weights                     = new double[q + 1]; //The weights for changing to another spin state
    neighbours                  = new double[q + 1]; //The number of neighbours (or weights) in different spin states
    csize                       = new igraph_integer_t[q + 1]; //The number of nodes in each community


    //Initialize communities
    for (igraph_integer_t i = 0; i <= q; i++) {
        degree_community_pos_in[i]  = 0.0;
        degree_community_neg_in[i]  = 0.0;
        degree_community_pos_out[i] = 0.0;
        degree_community_neg_out[i] = 0.0;

        csize[i]                    = 0;
    }

    //Initialize vectors
    if (init_spins) {
        for (igraph_integer_t i = 0; i < num_nodes; i++) {
            degree_pos_in[i]    = 0.0;
            degree_neg_in[i]    = 0.0;
            degree_pos_out[i]   = 0.0;
            degree_neg_out[i]   = 0.0;

#ifdef SPINGLASS_DEBUG
            printf("Initializing spin %d", i);
#endif
            spin[i] = 0;
        }
    }
    m_p = 0.0;
    m_n = 0.0;
    //Set community for each node, and
    //correctly store it in the bookkeeping

    double sum_weight_pos_in, sum_weight_pos_out, sum_weight_neg_in, sum_weight_neg_out;

    for (igraph_integer_t v = 0; v < num_nodes; v++) {
        if (init_spins) {
            s = RNG_INTEGER(1, q);  //The new spin s
            spin[v] = s;
        } else {
            s = spin[v];
        }

#ifdef SPINGLASS_DEBUG
        printf("Spin %d assigned to node %d.\n", s, v);
#endif

        n_cur               =  net->node_list.Get(v);

        l_cur               = l_iter.First(n_cur->Get_Links());

        sum_weight_pos_in   = 0.0;
        sum_weight_pos_out  = 0.0;
        sum_weight_neg_in   = 0.0;
        sum_weight_neg_out  = 0.0;

        while (!l_iter.End()) {
            double w = l_cur->Get_Weight();
            if (l_cur->Get_Start() == n_cur) //From this to other, so outgoing link
                if (w > 0) {
                    sum_weight_pos_out += w;    //Increase positive outgoing weight
                } else {
                    sum_weight_neg_out -= w;    //Increase negative outgoing weight
                } else if (w > 0) {
                sum_weight_pos_in += w;    //Increase positive incoming weight
            } else {
                sum_weight_neg_in -= w;    //Increase negative incoming weight
            }

            l_cur = l_iter.Next();
        }

        if (!is_directed) {
            double sum_weight_pos       = sum_weight_pos_out + sum_weight_pos_in;
            sum_weight_pos_out   = sum_weight_pos;
            sum_weight_pos_in    = sum_weight_pos;
            double sum_weight_neg = sum_weight_neg_out + sum_weight_neg_in;
            sum_weight_neg_out   = sum_weight_neg;
            sum_weight_neg_in    = sum_weight_neg;
        }

        if (init_spins) {
            //Set the degrees correctly
            degree_pos_in[v]    = sum_weight_pos_in;
            degree_neg_in[v]    = sum_weight_neg_in;
            degree_pos_out[v]   = sum_weight_pos_out;
            degree_neg_out[v]   = sum_weight_neg_out;
        }

        //Correct the community bookkeeping
        degree_community_pos_in[s]  += sum_weight_pos_in;
        degree_community_neg_in[s]  += sum_weight_neg_in;
        degree_community_pos_out[s] += sum_weight_pos_out;
        degree_community_neg_out[s] += sum_weight_neg_out;

        //Community just increased
        csize[s]++;

        //Sum the weights (notice that sum of indegrees equals sum of outdegrees)
        m_p += sum_weight_pos_in;
        m_n += sum_weight_neg_in;
    }

#ifdef SPINGLASS_DEBUG
    printf("Done assigning.\n");
#endif
}
//##############################################################
// This is the function generally used for optimisation,
// as the parallel update has its flaws, due to the cyclic attractors
//##############################################################
double PottsModelN::HeatBathLookup(double gamma, double lambda, double t, unsigned int max_sweeps) {
#ifdef SPINGLASS_DEBUG
    printf("Starting sweep at temperature %f.\n", t);
#endif
    DLList_Iter<NLink*> l_iter;
    const NNode *node, *n_cur;
    const NLink *l_cur;
    /* The new_spin contains the spin to which we will update,
     * the spin_opt is the optional spin we will consider and
     * the old_spin is the spin of the node we are currently
     * changing.
     */
    igraph_integer_t new_spin, spin_opt, old_spin;
    unsigned int sweep; //current sweep
    igraph_integer_t changes/*, problemcount*/; //Number of changes and number of problems encountered

    double exp_old_spin; //The expectation value for the old spin
    double exp_spin; //The expectation value for the other spin(s)
    igraph_integer_t v; //The node we will be investigating

    //The variables required for the calculations
    double delta_pos_out, delta_pos_in, delta_neg_out, delta_neg_in;
    double k_v_pos_out, k_v_pos_in, k_v_neg_out, k_v_neg_in;

    //weight of edge
    double w;

    double beta = 1.0 / t; //Weight for probabilities
    double r = 0.0; //random number used for assigning new spin

    double maxweight = 0.0;
    double sum_weights = 0.0; //sum_weights for normalizing the probabilities

    sweep = 0;
    changes = 0;
    double m_pt = m_p;
    double m_nt = m_n;

    if (m_pt < 0.001) {
        m_pt = 1;
    }

    if (m_nt < 0.001) {
        m_nt = 1;
    }

    while (sweep < max_sweeps) {
        sweep++;
        //loop over all nodes in network
        for (igraph_integer_t n = 0; n < num_nodes; n++) {
            //Look for a random node
            v = RNG_INTEGER(0, num_nodes - 1);
            //We will be investigating node v

            node = net->node_list.Get(v);

            /*******************************************/
            // initialize the neighbours and the weights
            // problemcount = 0;
            for (igraph_integer_t i = 0; i <= q; i++) {
                neighbours[i] = 0.0;
                weights[i] = 0.0;
            }

            //Loop over all links (=neighbours)
            l_cur = l_iter.First(node->Get_Links());
            while (!l_iter.End()) {
                w = l_cur->Get_Weight();
                if (node == l_cur->Get_Start()) {
                    n_cur = l_cur->Get_End();
                } else {
                    n_cur = l_cur->Get_Start();
                }
                //Add the link to the correct cluster
                neighbours[spin[n_cur->Get_Index()]] += w;
                l_cur = l_iter.Next();
            }
            //We now have the weight of the (in and out) neighbours
            //in each cluster available to us.
            /*******************************************/
            old_spin = spin[v];

            //Look for optimal spin

            //Set the appropriate variable
            delta_pos_out   = degree_pos_out[v];
            delta_pos_in    = degree_pos_in[v];
            delta_neg_out   = degree_neg_out[v];
            delta_neg_in    = degree_neg_in[v];

            k_v_pos_out     = gamma * delta_pos_out / m_pt;
            k_v_pos_in      = gamma * delta_pos_in / m_pt;
            k_v_neg_out     = lambda * delta_neg_out / m_nt;
            k_v_neg_in      = lambda * delta_neg_in / m_nt;

            //The expectation value for the old spin
            if (is_directed)
                exp_old_spin = (k_v_pos_out * (degree_community_pos_in[old_spin] - delta_pos_in) -
                                k_v_neg_out * (degree_community_neg_in[old_spin] - delta_neg_in)) +
                               (k_v_pos_in * (degree_community_pos_out[old_spin] - delta_pos_out) -
                                k_v_neg_in * (degree_community_neg_out[old_spin] - delta_neg_out));
            else
                exp_old_spin = (k_v_pos_out * (degree_community_pos_in[old_spin] - delta_pos_in) -
                                k_v_neg_out * (degree_community_neg_in[old_spin] - delta_neg_in));

            /*******************************************/
            //Calculating probabilities for each transition to another
            //community.

            maxweight = 0.0;
            weights[old_spin] = 0.0;

            for (spin_opt = 1; spin_opt <= q; spin_opt++) { // all possible new spins
                if (spin_opt != old_spin) { // except the old one!
                    if (is_directed)
                        exp_spin = (k_v_pos_out * degree_community_pos_in[spin_opt] - k_v_neg_out * degree_community_neg_in[spin_opt]) +
                                   (k_v_pos_in * degree_community_pos_out[spin_opt] - k_v_neg_in * degree_community_neg_out[spin_opt]);
                    else {
                        exp_spin = (k_v_pos_out * degree_community_pos_in[spin_opt] - k_v_neg_out * degree_community_neg_in[spin_opt]);
                    }

                    weights[spin_opt] = (neighbours[spin_opt] - exp_spin) - (neighbours[old_spin] - exp_old_spin);

                    if (weights[spin_opt] > maxweight) {
                        maxweight = weights[spin_opt];
                    }
                }
            }   // for spin

            //Calculate exp. prob. an
            sum_weights = 0.0;
            for (spin_opt = 1; spin_opt <= q; spin_opt++) { // all possible new spins
                weights[spin_opt] -= maxweight;  //subtract maxweight for numerical stability (otherwise overflow).
                weights[spin_opt]  = exp(beta * weights[spin_opt]);
                sum_weights   += weights[spin_opt];
            }   // for spin
            /*******************************************/


            /*******************************************/
            //Choose a new spin dependent on the calculated probabilities
            r = RNG_UNIF(0, sum_weights);
            new_spin = 1;

            while (new_spin <= q) {
                if (r <= weights[new_spin]) {
                    spin_opt = new_spin; //We have found are new spin
                    break;
                } else {
                    r -= weights[new_spin];    //Perhaps the next spin is the one we want
                }

                new_spin++;
            }

            new_spin = spin_opt;
            //If there wasn't a problem we should have found
            //our new spin.
            /*******************************************/


            /*******************************************/
            //The new spin is available to us, so change
            //all the appropriate counters.
            if (new_spin != old_spin) { // Did we really change something??
                changes++;
                spin[v] = new_spin;

                //The new spin increase by one, and the old spin decreases by one
                csize[new_spin]++; csize[old_spin]--;

                //Change the sums of degree for the old spin...
                degree_community_pos_in[old_spin]   -= delta_pos_in;
                degree_community_neg_in[old_spin]   -= delta_neg_in;
                degree_community_pos_out[old_spin]  -= delta_pos_out;
                degree_community_neg_out[old_spin]  -= delta_neg_out;

                //...and for the new spin
                degree_community_pos_in[new_spin]   += delta_pos_in;
                degree_community_neg_in[new_spin]   += delta_neg_in;
                degree_community_pos_out[new_spin]  += delta_pos_out;
                degree_community_neg_out[new_spin]  += delta_neg_out;
            }

            //We have no change a node from old_spin to new_spin
            /*******************************************/

        } // for n
    }  // while sweep
#ifdef SPINGLASS_DEBUG
    printf("Done %d sweeps.\n", max_sweeps);
    printf("%ld changes made for %d nodes.\n", changes, num_nodes);
    printf("Last node is %d and last random number is %f with sum of weights %f with spin %d.\n", v, r, sum_weights, old_spin);
#endif

    return (double(changes) / double(num_nodes) / double(sweep));
}

//We need to begin at a suitable temperature. That is, a temperature at which
//enough nodes may change their initially assigned communties
double PottsModelN::FindStartTemp(double gamma, double lambda, double ts) {
    double kT;
    kT = ts;
    //assing random initial condition
    assign_initial_conf(true);
    // the factor 1-1/q is important, since even, at infinite temperature,
    // only 1-1/q of all spins do change their state, since a randomly chooses new
    // state is with prob. 1/q the old state.
    double acceptance = 0.0;
    while (acceptance < (1.0 - 1.0 / double(q)) * 0.95) { //want 95% acceptance
        kT = kT * 1.1;
        acceptance = HeatBathLookup(gamma, lambda, kT, 50);
    }
    kT *= 1.1; // just to be sure...
    return kT;
}

igraph_integer_t PottsModelN::WriteClusters(igraph_real_t *modularity,
                                igraph_real_t *temperature,
                                igraph_vector_int_t *community_size,
                                igraph_vector_int_t *membership,
                                igraph_matrix_t *adhesion,
                                igraph_matrix_t *normalised_adhesion,
                                igraph_real_t *polarization,
                                double t,
                                double d_p,
                                double d_n) {

#ifdef SPINGLASS_DEBUG
    printf("Start writing clusters.\n");
#endif
    //Reassign each community so that we retrieve a community assignment 1 through num_communities
    auto *cluster_assign = new igraph_integer_t[q + 1];
    for (igraph_integer_t i = 0; i <= q; i++) {
        cluster_assign[i] = 0;
    }

    igraph_integer_t num_clusters = 0;

    //Find out what the new communities will be
    for (igraph_integer_t i = 0; i < num_nodes; i++) {
        igraph_integer_t s = spin[i];
        if (cluster_assign[s] == 0) {
            num_clusters++;
            cluster_assign[s] = num_clusters;
#ifdef SPINGLASS_DEBUG
            printf("Setting cluster %d to %d.\n", s, num_clusters);
#endif
        }
    }

    //And now assign each node to its new community
    q = num_clusters;
    for (igraph_integer_t i = 0; i < num_nodes; i++) {
#ifdef SPINGLASS_DEBUG
        printf("Setting node %d to %d.\n", i, cluster_assign[spin[i]]);
#endif
        igraph_integer_t s = cluster_assign[spin[i]];
        spin[i] = s;
#ifdef SPINGLASS_DEBUG
        printf("Have set node %d to %d.\n", i, s);
#endif
    }
    assign_initial_conf(false);

    delete [] cluster_assign;

    if (temperature) {
        *temperature = t;
    }

    if (community_size) {
        //Initialize the vector
        IGRAPH_CHECK(igraph_vector_int_resize(community_size, q));
        for (igraph_integer_t spin_opt = 1; spin_opt <= q; spin_opt++) {
            //Set the community size
            VECTOR(*community_size)[spin_opt - 1] = csize[spin_opt];
        }
    }

    //Set the membership
    if (membership) {
        IGRAPH_CHECK(igraph_vector_int_resize(membership, num_nodes));
        for (igraph_integer_t i = 0; i < num_nodes; i++) {
            VECTOR(*membership)[ i ] = spin[i] - 1;
        }
    }

    double Q = 0.0; //Modularity
    if (adhesion) {
        IGRAPH_CHECK(igraph_matrix_resize(adhesion, q, q));
        IGRAPH_CHECK(igraph_matrix_resize(normalised_adhesion, q, q));

        double **num_links_pos = nullptr;
        double **num_links_neg = nullptr;
        //memory allocated for elements of rows.
        num_links_pos = new double *[q + 1] ;
        num_links_neg = new double *[q + 1] ;

        //memory allocated for  elements of each column.
        for ( igraph_integer_t i = 0 ; i < q + 1 ; i++) {
            num_links_pos[i] = new double[q + 1];
            num_links_neg[i] = new double[q + 1];
        }



        //Init num_links
        for (igraph_integer_t i = 0; i <= q; i++) {
            for (igraph_integer_t j = 0; j <= q; j++) {
                num_links_pos[i][j] = 0.0;
                num_links_neg[i][j] = 0.0;
            }
        }

        DLList_Iter<NLink*> iter_l;
        const NLink *l_cur = iter_l.First(&net->link_list);

        double w = 0.0;

        while (!iter_l.End()) {
            w = l_cur->Get_Weight();
            igraph_integer_t a = spin[l_cur->Get_Start()->Get_Index()];
            igraph_integer_t b = spin[l_cur->Get_End()->Get_Index()];
            if (w > 0) {
                num_links_pos[a][b] += w;
                if (!is_directed && a != b) { //Only one edge is defined in case it is undirected
                    num_links_pos[b][a] += w;
                }
            } else {
                num_links_neg[a][b] -= w;
                if (!is_directed && a != b) { //Only one edge is defined in case it is undirected
                    num_links_neg[b][a] -= w;
                }
            }

            l_cur = iter_l.Next();
        } //while links

#ifdef SPINGLASS_DEBUG
        printf("d_p: %f\n", d_p);
        printf("d_n: %f\n", d_n);
#endif

        double expected = 0.0;
        double a = 0.0;
        double normal_a = 0.0;

        double delta, u_p, u_n;
        double max_expected, max_a;

        //We don't take into account the lambda or gamma for
        //computing the modularity and adhesion, since they
        //are then incomparable to other definitions.
        for (igraph_integer_t i = 1; i <= q; i++) {
            for (igraph_integer_t j = 1; j <= q; j++) {
                if (!is_directed && i == j)
                    expected    = degree_community_pos_out[i] * degree_community_pos_in[j] / (m_p == 0 ? 1 : 2 * m_p)
                                  - degree_community_neg_out[i] * degree_community_neg_in[j] / (m_n == 0 ? 1 : 2 * m_n);
                else
                    expected    = degree_community_pos_out[i] * degree_community_pos_in[j] / (m_p == 0 ? 1 : m_p)
                                  - degree_community_neg_out[i] * degree_community_neg_in[j] / (m_n == 0 ? 1 : m_n);

                a           = (num_links_pos[i][j] - num_links_neg[i][j]) - expected;

                if (i == j) { //cohesion
                    if (is_directed) {
                        delta = d_p * csize[i] * (csize[i] - 1);    //Maximum amount
                    } else {
                        delta = d_p * csize[i] * (csize[i] - 1) / 2;    //Maximum amount
                    }

                    u_p     = delta - num_links_pos[i][i]; //Add as many positive links we can
                    u_n     = -num_links_neg[i][i]; //Delete as many negative links we can
                    Q      += a;
                } else { //adhesion
                    if (is_directed) {
                        delta = d_n * csize[i] * csize[j] * 2;    //Maximum amount
                    } else {
                        delta = d_n * csize[i] * csize[j];    //Maximum amount
                    }

                    u_p     = -num_links_pos[i][j]; //Delete as many positive links we can
                    u_n     = delta - num_links_neg[i][j]; //Add as many negative links we can
                }

                if (!is_directed && i == j)
                    max_expected    = (degree_community_pos_out[i] + u_p) * (degree_community_pos_in[j] + u_p) / ((m_p + u_p) == 0 ? 1 : 2 * (m_p + u_p))
                                      - (degree_community_neg_out[i] - u_n) * (degree_community_neg_in[j] + u_n) / ((m_n + u_n) == 0 ? 1 : 2 * (m_n + u_n));
                else
                    max_expected    = (degree_community_pos_out[i] + u_p) * (degree_community_pos_in[j] + u_p) / ((m_p + u_p) == 0 ? 1 : m_p + u_p)
                                      - (degree_community_neg_out[i] - u_n) * (degree_community_neg_in[j] + u_n) / ((m_n + u_n) == 0 ? 1 : m_n + u_n);
                max_a           = ((num_links_pos[i][j] + u_p) - (num_links_neg[i][j] + u_n)) - max_expected;


                //In cases where we haven't actually found a ground state
                //the adhesion/cohesion *might* not be negative/positive,
                //hence the maximum adhesion and cohesion might behave quite
                //strangely. In order to prevent that, we limit them to 1 in
                //absolute value, and prevent from dividing by zero (even if
                //chuck norris would).
                if (i == j) {
                    normal_a = a / (max_a == 0 ? a : max_a);
                } else {
                    normal_a = -a / (max_a == 0 ? a : max_a);
                }

                if (normal_a > 1) {
                    normal_a = 1;
                } else if (normal_a < -1) {
                    normal_a = -1;
                }

                MATRIX(*adhesion, i - 1, j - 1) = a;
                MATRIX(*normalised_adhesion, i - 1, j - 1) = normal_a;
            } //for j
            //printf("\n");
        } //for i

        //free the allocated memory
        for ( igraph_integer_t i = 0 ; i < q + 1 ; i++ ) {
            delete [] num_links_pos[i] ;
            delete [] num_links_neg[i];
        }
        delete [] num_links_pos ;
        delete [] num_links_neg ;

    } //adhesion

    if (modularity) {
        if (is_directed) {
            *modularity = Q / (m_p + m_n);
        } else {
            *modularity = 2 * Q / (m_p + m_n);    //Correction for the way m_p and m_n are counted. Modularity is 1/m, not 1/2m
        }
    }

    if (polarization) {
        double sum_ad = 0.0;
        for (igraph_integer_t i = 0; i < q; i++) {
            for (igraph_integer_t j = 0; j < q; j++) {
                if (i != j) {
                    sum_ad -= MATRIX(*normalised_adhesion, i, j);
                }
            }
        }
        *polarization = sum_ad / (q * q - q);
    }
#ifdef SPINGLASS_DEBUG
    printf("Finished writing cluster.\n");
#endif
    return num_nodes;
}
