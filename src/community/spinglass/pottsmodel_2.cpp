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
#include "NetRoutines.h"

#include "igraph_random.h"
#include "core/interruption.h"

#include <cstring>
#include <cmath>

using namespace std;

//#################################################################################################
PottsModel::PottsModel(network *n, unsigned int qvalue, int m) : Qmatrix(qvalue+1), acceptance(0)
{
    DLList_Iter<NNode*> iter;
    NNode *n_cur;
    unsigned int *i_ptr;
    net = n;
    q = qvalue;
    operation_mode = m;
    k_max = 0;
    //needed in calculating modularity
    Qa     = new double[q + 1];
    //weights for each spin state needed in Monte Carlo process
    weights = new double[q + 1];
    //bookkeeping of occupation numbers of spin states or the number of links in community
    color_field = new double[q + 1];
    neighbours = new double[q + 1];

    num_of_nodes = net->node_list->Size();
    num_of_links = net->link_list->Size();

    n_cur = iter.First(net->node_list);
    //these lists are needed to keep track of spin states for parallel update mode
    new_spins = new DL_Indexed_List<unsigned int*>();
    previous_spins = new DL_Indexed_List<unsigned int*>();
    while (!iter.End()) {
        if (k_max < n_cur->Get_Degree()) {
            k_max = n_cur->Get_Degree();
        }
        i_ptr = new unsigned int;
        *i_ptr = 0;
        new_spins->Push(i_ptr);
        i_ptr = new unsigned int;
        *i_ptr = 0;
        previous_spins->Push(i_ptr);
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
    new_spins->delete_items();
    previous_spins->delete_items();
    delete new_spins;
    delete previous_spins;
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
unsigned long PottsModel::assign_initial_conf(int spin) {
    int s;
    DLList_Iter<NNode*> iter;
    DLList_Iter<NLink*> l_iter;
    NNode *n_cur;
    NLink *l_cur;
    double sum_weight;
    double av_k_squared = 0.0;
    double av_k = 0.0;
//   printf("Assigning initial configuration...\n");
    // initialize colorfield
    for (unsigned int i = 0; i <= q; i++) {
        color_field[i] = 0.0;
    }
    //
    total_degree_sum = 0.0;
    n_cur = iter.First(net->node_list);
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
        av_k_squared += sum_weight * sum_weight;
        av_k += sum_weight;

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
    av_k_squared /= double(net->node_list->Size());
    av_k /= double(net->node_list->Size());
    // total_degree_sum-=av_k_squared/av_k;
//   printf("Total Degree Sum=2M=%f\n",total_degree_sum);
    return net->node_list->Size();
}
//#####################################################################
//If I ever manage to write a decent LookUp function, it will be here
//#####################################################################
unsigned long PottsModel::initialize_lookup(double kT, double gamma) {
    IGRAPH_UNUSED(kT);
    IGRAPH_UNUSED(gamma);
    /*
    double beta;
    // the look-up table contains all entries of exp(-beta(-neighbours+gamma*h))
    // as needed in the HeatBath algorithm
    beta=1.0/kT;
    for (long w=0; w<=k_max+num_of_nodes; w++)
    {
       neg_lookup[w]=exp(-beta*-w
    }
    delta_ij[0]=1.0;
    for (long w=-num_of_nodes-k_max; w<=k_max+num_of_nodes; w++)
    {

    }

    // wenn wir spaeter exp(-1/kT*gamma*(nk+1-nj) fuer eine spin-flip von j nach k benoetigen schauen wir nur noch hier nach
    for (unsigned long n=1; n<=num_of_nodes; n++)
    {
      gamma_term[n]=exp(-double(n)/kT*gamma);
    }
    gamma_term[0]=1.0;
    */
    return 1;
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
    unsigned int i, j;
    //initialize with zeros
    num_of_links = net->link_list->Size();
    for (i = 0; i <= q; i++) {
        Qa[i] = 0.0;
        for (j = i; j <= q; j++) {
            Qmatrix[i][j] = 0.0;
            Qmatrix[j][i] = 0.0;
        }
    }
    //go over all links and make corresponding entries in Q matrix
    //An edge connecting state i wiht state j will get an entry in Qij and Qji
    l_cur = l_iter.First(net->link_list);
    while (!l_iter.End()) {
        i = l_cur->Get_Start()->Get_ClusterIndex();
        j = l_cur->Get_End()->Get_ClusterIndex();
        //printf("%d %d\n",i,j);
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
    for (unsigned int i = 0; i <= q; i++) {
        Q += Qmatrix[i][i] - Qa[i] * Qa[i] / double(2.0 * net->sum_weights);
        if ((Qa[i] < 0.0) || Qmatrix[i][i] < 0.0) {
//         printf("Negatives Qa oder Qii\n\n\n");
            //printf("Press any key to continue\n\n");
            //cin >> Q;
        }
    }
    Q /= double(2.0 * net->sum_weights);
    return Q;
}
double PottsModel::calculate_genQ(double gamma) {
    double Q = 0.0;
    for (unsigned int i = 0; i <= q; i++) {
        Q += Qmatrix[i][i] - gamma * Qa[i] * Qa[i] / double(2.0 * net->sum_weights);
        if ((Qa[i] < 0.0) || Qmatrix[i][i] < 0.0) {
//         printf("Negatives Qa oder Qii\n\n\n");
            //printf("Press any key to continue\n\n");
            //cin >> Q;
        }
    }
    Q /= double(2.0 * net->sum_weights);
    return Q;
}
//#######################################################################
// This function calculates the Energy for the standard Hamiltonian
// given a particular value of gamma and the current spin states
// #####################################################################
double PottsModel::calculate_energy(double gamma) {
    double e = 0.0;
    DLList_Iter<NLink*> l_iter;
    NLink *l_cur;
    l_cur = l_iter.First(net->link_list);
    //every in-cluster edge contributes -1
    while (!l_iter.End()) {
        if (l_cur->Get_Start()->Get_ClusterIndex() == l_cur->Get_End()->Get_ClusterIndex()) {
            e--;
        }
        l_cur = l_iter.Next();
    }
    //and the penalty term contributes according to cluster sizes
    for (unsigned int i = 1; i <= q; i++) {
        e += gamma * 0.5 * double(color_field[i]) * double((color_field[i] - 1));
    }
    energy = e;
    return e;
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
        // if I ever have a lookup table, it will need initialization for every kT
        //initialize_lookup(kT,k_max,net->node_list->Size());
        HeatBathParallelLookup(gamma, prob, kT, 50);
//        printf("kT=%f acceptance=%f\n", kT, acceptance);
    }
    kT *= 1.1; // just to be sure...
//   printf("Starting with acceptance ratio: %1.6f bei kT=%2.4f\n",acceptance,kT);
    return kT;
}

//##############################################################
//This function does a parallel update at zero T
//Hence, it is really fast on easy problems
//max sweeps is the maximum number of sweeps it should perform,
//if it does not converge earlier
//##############################################################
long PottsModel::HeatBathParallelLookupZeroTemp(double gamma, double prob, unsigned int max_sweeps) {
    DLList_Iter<NNode*> iter, net_iter;
    DLList_Iter<NLink*> l_iter;
    DLList_Iter<unsigned int*> i_iter, i_iter2;
    NNode *node, *n_cur;
    NLink *l_cur;
    unsigned int *SPIN, *P_SPIN, new_spin, spin_opt, old_spin, spin, sweep;
    // long h; // degree;
    unsigned long changes;
    double h, delta = 0, deltaE, deltaEmin, w, degree;
    //HugeArray<double> neighbours;
    bool cyclic = false;

    sweep = 0;
    changes = 1;
    while (sweep < max_sweeps && changes) {
        cyclic = true;
        sweep++;
        changes = 0;
        //Loop over all nodes
        node = net_iter.First(net->node_list);
        SPIN = i_iter.First(new_spins);
        while (!net_iter.End()) {
            // How many neigbors of each type?
            // set them all zero
            for (unsigned int i = 0; i <= q; i++) {
                neighbours[i] = 0;
            }
            degree = node->Get_Weight();
            //Loop over all links (=neighbours)
            l_cur = l_iter.First(node->Get_Links());
            while (!l_iter.End()) {
                //printf("%s %s\n",node->Get_Name(),n_cur->Get_Name());
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
            //degree=node->Get_Degree();
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
            }


            spin_opt = old_spin;
            deltaEmin = 0.0;
            for (spin = 1; spin <= q; spin++) { // all possible spin states
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
        node = net_iter.First(net->node_list);
        SPIN = i_iter.First(new_spins);
        P_SPIN = i_iter2.First(previous_spins);
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
//       printf("Cyclic attractor!\n");
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
    DLList_Iter<NNode*> iter;
    DLList_Iter<NLink*> l_iter;
    DLList_Iter<unsigned int*> i_iter, i_iter2;
    NNode *node, *n_cur;
    NLink *l_cur;
    unsigned int new_spin, spin_opt, old_spin, spin, sweep;
    long r;// degree;
    unsigned long changes;
    double delta = 0, h, deltaE, deltaEmin, w, degree;
    //HugeArray<int> neighbours;

    sweep = 0;
    changes = 0;
    while (sweep < max_sweeps) {
        sweep++;
        //ueber alle Knoten im Netz
        for (unsigned long n = 0; n < num_of_nodes; n++) {
            r = -1;
            while ((r < 0) || (r > (long)num_of_nodes - 1)) {
                r = RNG_INTEGER(0, num_of_nodes - 1);
            }
            /* r=long(double(num_of_nodes*double(rand())/double(RAND_MAX+1.0)));*/
            node = net->node_list->Get(r);
            // Wir zaehlen, wieviele Nachbarn von jedem spin vorhanden sind
            // erst mal alles Null setzen
            for (unsigned int i = 0; i <= q; i++) {
                neighbours[i] = 0;
            }
            degree = node->Get_Weight();
            //Loop over all links (=neighbours)
            l_cur = l_iter.First(node->Get_Links());
            while (!l_iter.End()) {
                //printf("%s %s\n",node->Get_Name(),n_cur->Get_Name());
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
            //degree=node->Get_Degree();
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
            }


            spin_opt = old_spin;
            deltaEmin = 0.0;
            for (spin = 1; spin <= q; spin++) { // alle moeglichen Spins
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
long PottsModel::HeatBathParallelLookup(double gamma, double prob, double kT, unsigned int max_sweeps) {
    DLList_Iter<NNode*> iter, net_iter;
    DLList_Iter<NLink*> l_iter;
    DLList_Iter<unsigned int*> i_iter, i_iter2;
    NNode *node, *n_cur;
    NLink *l_cur;
    unsigned int new_spin, spin_opt, old_spin;
    unsigned int *SPIN, *P_SPIN;
    unsigned int sweep;
    long max_q;
    unsigned long changes, /*degree,*/ problemcount;
    //HugeArray<int> neighbours;
    double h, delta = 0, norm, r, beta, minweight, prefac = 0, w, degree;
    bool cyclic = false, found;
    unsigned long number_of_nodes;

    sweep = 0;
    changes = 1;
    number_of_nodes = net->node_list->Size();
    while (sweep < max_sweeps && changes) {
        cyclic = true;
        sweep++;
        changes = 0;
        //Loop over all nodes
        node = net_iter.First(net->node_list);
        SPIN = i_iter.First(new_spins);
        while (!net_iter.End()) {
            // Initialize neighbours and weights
            problemcount = 0;
            for (unsigned int i = 0; i <= q; i++) {
                neighbours[i] = 0;
                weights[i] = 0;
            }
            norm = 0.0;
            degree = node->Get_Weight();
            //Loop over all links (=neighbours)
            l_cur = l_iter.First(node->Get_Links());
            while (!l_iter.End()) {
                //printf("%s %s\n",node->Get_Name(),n_cur->Get_Name());
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
            //degree=node->Get_Degree();
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
            }
            spin_opt = old_spin;
            beta = 1.0 / kT * prefac;
            minweight = 0.0;
            weights[old_spin] = 0.0;
            for (unsigned spin = 1; spin <= q; spin++) { // loop over all possible new spins
                if (spin != old_spin) { // only if we have a different than old spin!
                    h = color_field[spin] + delta - color_field[old_spin];
                    weights[spin] = double(neighbours[old_spin] - neighbours[spin]) + gamma * prob * double(h);
                    if (weights[spin] < minweight) {
                        minweight = weights[spin];
                    }
                }
            }   // for spin
            for (unsigned spin = 1; spin <= q; spin++) { // loop over all possibe spins
                weights[spin] -= minweight;       // subtract minweight
                // to avoid numerical problems with large exponents
                weights[spin] = exp(-beta * weights[spin]);
                norm += weights[spin];
            }   // for spin

            //now choose a new spin
            r = RNG_UNIF(0, norm);
            /* norm*double(rand())/double(RAND_MAX + 1.0); */
            new_spin = 1;
            found = false;
            while (!found && new_spin <= q) {
                if (r <= weights[new_spin]) {
                    spin_opt = new_spin;
                    found = true;
                    break;
                } else {
                    r -= weights[new_spin];
                }
                new_spin++;
            }
            if (!found) {
//         printf(".");
                problemcount++;
            }
            //Put new spin on list
            *SPIN = spin_opt;

            node = net_iter.Next();
            SPIN = i_iter.Next();
        } // while !net_iter.End()

        //-------------------------------
        //now update all spins
        node = net_iter.First(net->node_list);
        SPIN = i_iter.First(new_spins);
        P_SPIN = i_iter2.First(previous_spins);
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
    for (unsigned int i = 1; i <= q; i++) if (color_field[i] > max_q) {
            max_q = long(color_field[i]);
        }

    //again, we would not like to end up in cyclic attractors
    if (cyclic && changes)  {
//       printf("Cyclic attractor!\n");
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
    DLList_Iter<NNode*> iter;
    DLList_Iter<NLink*> l_iter;
    DLList_Iter<unsigned int*> i_iter, i_iter2;
    NNode *node, *n_cur;
    NLink *l_cur;
    unsigned int new_spin, spin_opt, old_spin;
    unsigned int sweep;
    long max_q, rn;
    unsigned long changes, /*degree,*/ problemcount;
    double degree, w, delta = 0, h;
    //HugeArray<int> neighbours;
    double norm, r, beta, minweight, prefac = 0;
    bool found;
    long int number_of_nodes;
    sweep = 0;
    changes = 0;
    number_of_nodes = net->node_list->Size();
    while (sweep < max_sweeps) {
        sweep++;
        //loop over all nodes in network
        for (int n = 0; n < number_of_nodes; n++) {
            rn = -1;
            while ((rn < 0) || (rn > number_of_nodes - 1)) {
                rn = RNG_INTEGER(0, number_of_nodes - 1);
            }
            /* rn=long(double(number_of_nodes*double(rand())/double(RAND_MAX+1.0))); */

            node = net->node_list->Get(rn);
            // initialize the neighbours and the weights
            problemcount = 0;
            for (unsigned int i = 0; i <= q; i++) {
                neighbours[i] = 0.0;
                weights[i] = 0.0;
            }
            norm = 0.0;
            degree = node->Get_Weight();
            //Loop over all links (=neighbours)
            l_cur = l_iter.First(node->Get_Links());
            while (!l_iter.End()) {
                //printf("%s %s\n",node->Get_Name(),n_cur->Get_Name());
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
            //degree=node->Get_Degree();
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
            }
            spin_opt = old_spin;
            beta = 1.0 / kT * prefac;
            minweight = 0.0;
            weights[old_spin] = 0.0;
            for (unsigned spin = 1; spin <= q; spin++) { // all possible new spins
                if (spin != old_spin) { // except the old one!
                    h = color_field[spin] - (color_field[old_spin] - delta);
                    weights[spin] = neighbours[old_spin] - neighbours[spin] + gamma * prob * h;
                    if (weights[spin] < minweight) {
                        minweight = weights[spin];
                    }
                }
            }   // for spin
            for (unsigned spin = 1; spin <= q; spin++) { // all possible new spins
                weights[spin] -= minweight;       // subtract minweigt
                // for numerical stability
                weights[spin] = exp(-beta * weights[spin]);
                norm += weights[spin];
            }   // for spin


            //choose a new spin
            /*      r = norm*double(rand())/double(RAND_MAX + 1.0); */
            r = RNG_UNIF(0, norm);
            new_spin = 1;
            found = false;
            while (!found && new_spin <= q) {
                if (r <= weights[new_spin]) {
                    spin_opt = new_spin;
                    found = true;
                    break;
                } else {
                    r -= weights[new_spin];
                }
                new_spin++;
            }
            if (!found) {
//         printf(".");
                problemcount++;
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

    for (unsigned int i = 1; i <= q; i++) if (color_field[i] > max_q) {
            max_q = long(color_field[i] + 0.5);
        }

    acceptance = double(changes) / double(number_of_nodes) / double(sweep);
    return acceptance;
}

//###############################################################################################
//# Here we try to minimize the affinity to the rest of the network
//###############################################################################################
double PottsModel::FindCommunityFromStart(double gamma, double prob,
        char *nodename,
        igraph_vector_t *result,
        igraph_real_t *cohesion,
        igraph_real_t *adhesion,
        igraph_integer_t *my_inner_links,
        igraph_integer_t *my_outer_links) {
    DLList_Iter<NNode*> iter, iter2;
    DLList_Iter<NLink*> l_iter;
    DLList<NNode*>* to_do;
    DLList<NNode*>* community;
    NNode *start_node = NULL, *n_cur, *neighbor, *max_aff_node, *node;
    NLink *l_cur;
    bool found = false, add = false, remove = false;
    double degree, delta_aff_add, delta_aff_rem, max_delta_aff, Ks = 0.0, Kr = 0, kis, kir, w;
    long community_marker = 5;
    long to_do_marker = 10;
    double inner_links = 0, outer_links = 0, aff_r, aff_s;

    IGRAPH_UNUSED(prob);

    to_do = new DLList<NNode*>;
    community = new DLList<NNode*>;

    // find the node in the network
    n_cur = iter.First(net->node_list);
    while (!found && !iter.End()) {
        if (0 == strcmp(n_cur->Get_Name(), nodename)) {
            start_node = n_cur;
            found = true;
            start_node->Set_Affinity(0.0);
            community->Push(start_node);
            start_node->Set_Marker(community_marker);
            Ks = start_node->Get_Weight();
            Kr = total_degree_sum - start_node->Get_Weight();
        }
        n_cur = iter.Next();
    }
    if (!found) {
//      printf("%s not found found. Aborting.\n",nodename);
//      fprintf(file,"%s not found found. Aborting.\n",nodename);
        delete to_do;
        delete community;
        return -1;
    }
    //#############################
    // initialize the to_do list and community with the neighbours of start node
    //#############################
    neighbor = iter.First(start_node->Get_Neighbours());
    while (!iter.End()) {
//     printf("Adding node %s to comunity.\n",neighbor->Get_Name());
        community->Push(neighbor);
        neighbor->Set_Marker(community_marker);
        Ks += neighbor->Get_Weight();
        Kr -= neighbor->Get_Weight();
        neighbor = iter.Next();
    }
    node = iter.First(community);
    while (!iter.End()) {
        //now add at the second neighbors to the to_do list
        neighbor = iter2.First(node->Get_Neighbours());
        while (!iter2.End()) {
            if ((long)neighbor->Get_Marker() != community_marker && (long)neighbor->Get_Marker() != to_do_marker) {
                to_do->Push(neighbor);
                neighbor->Set_Marker(to_do_marker);
//  printf("Adding node %s to to_do list.\n",neighbor->Get_Name());
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
        max_aff_node = NULL;
        add = false;
        node = iter.First(to_do);
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
                if ((long)n_cur->Get_Marker() == community_marker) {
                    kis += w; //the weight/number of links to the community
                } else {
                    kir += w; //the weight/number of links to the rest of the network
                }
                l_cur = l_iter.Next();
            }
            aff_r = kir - gamma / total_degree_sum * (Kr - degree) * degree;
            aff_s = kis - gamma / total_degree_sum * Ks * degree;
            delta_aff_add = aff_r - aff_s;
            //  if (aff_s>=aff_r && delta_aff_add<=max_delta_aff) {
            if (delta_aff_add <= max_delta_aff) {
                node->Set_Affinity(aff_s);
                max_delta_aff = delta_aff_add;
                max_aff_node = node;
                add = true;
            }
            //printf("%s in to_do list with affinity %f\n",node->Get_Name(),node->Get_Affinity());
            node = iter.Next();
        }
        //################
        //calculate the affinity changes for removing every single node from the community
        //################
        inner_links = 0;
        outer_links = 0;
        remove = false;
        node = iter.First(community);
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
                if ((long)n_cur->Get_Marker() == community_marker) {
                    kis += w;
                    inner_links += w; //summing all w gives twice the number of inner links(weights)
                } else {
                    kir += w;
                    outer_links += w;
                }
                l_cur = l_iter.Next();
            }
//  if (kir+kis!=degree) {  printf("error kir=%f\tkis=%f\tk=%f\n",kir,kis,degree); }
            aff_r = kir - gamma / total_degree_sum * Kr * degree;
            aff_s = kis - gamma / total_degree_sum * (Ks - degree) * degree;
            delta_aff_rem = aff_s - aff_r;
            node->Set_Affinity(aff_s);
            // we should not remove the nodes, we have just added
            if (delta_aff_rem < max_delta_aff) {
                max_delta_aff = delta_aff_rem ;
                max_aff_node = node;
                remove = true;
                add = false;
            }
            //printf("%s in to_do list with affinity %f\n",node->Get_Name(),node->Get_Affinity());
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
            community->Push(max_aff_node);
            max_aff_node->Set_Marker(community_marker);
            //delete node from to_do
            to_do->fDelete(max_aff_node);
            //update the sum of degrees in the community
            Ks += max_aff_node->Get_Weight();
            Kr -= max_aff_node->Get_Weight();
//  printf("Adding node %s to community with affinity of %f delta_aff: %f.\n",max_aff_node->Get_Name(), max_aff_node->Get_Affinity(),max_delta_aff);
            //now add all neighbors of this node, that are not already
            //in the to_do list or in the community
            neighbor = iter.First(max_aff_node->Get_Neighbours());
            while (!iter.End()) {
                if ((long)neighbor->Get_Marker() != community_marker && (long)neighbor->Get_Marker() != to_do_marker) {
                    to_do->Push(neighbor);
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
            community->fDelete(max_aff_node);
            max_aff_node->Set_Marker(to_do_marker);
            //update the sum of degrees in the community
            Ks -= max_aff_node->Get_Weight();
            Kr += max_aff_node->Get_Weight();
            //add the node to to_do again
            to_do->Push(max_aff_node);
//  printf("Removing node %s from community with affinity of %f delta_aff: %f.\n",max_aff_node->Get_Name(), max_aff_node->Get_Affinity(),max_delta_aff);
        }
        IGRAPH_ALLOW_INTERRUPTION(); /* This is not clean.... */
    }
    //###################
    //write the node in the community to a file
    //###################
    // TODO return this instead of writing it
//   fprintf(file,"Number_of_nodes:\t%d\n",community->Size());
//   fprintf(file,"Inner_Links:\t%f\n",inner_links);
//   fprintf(file,"Outer_Links:\t%f\n",Ks-2*inner_links);
//   fprintf(file,"Cohesion:\t%f\n",inner_links-gamma/total_degree_sum*Ks*Ks*0.5);
//   fprintf(file,"Adhesion:\t%f\n",outer_links-gamma/total_degree_sum*Ks*Kr);
//   fprintf(file,"\n");
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
        node = iter.First(community);
        igraph_vector_resize(result, 0);
        while (!iter.End()) {
            // printf("%s in community.\n",node->Get_Name());
            // fprintf(file,"%s\t%f\n",node->Get_Name(),node->Get_Affinity());
            IGRAPH_CHECK(igraph_vector_push_back(result, node->Get_Index()));
            node = iter.Next();
        }
    }
//   printf("%d nodes in community around %s\n",community->Size(),start_node->Get_Name());
//   fclose(file);
    unsigned int size = community->Size();
    delete to_do;
    delete community;
    return size;
}

//################################################################################################
// this Function writes the clusters to disk
//################################################################################################
long PottsModel::WriteClusters(igraph_real_t *modularity,
                               igraph_real_t *temperature,
                               igraph_vector_t *csize,
                               igraph_vector_t *membership,
                               double kT, double gamma) {
    NNode *n_cur, *n_cur2;
    /*
    double a1,a2,a3,p,p1,p2;
    long n,N,lin,lout;
    */
    DLList_Iter<NNode*> iter, iter2;
    HugeArray<int> inner_links;
    HugeArray<int> outer_links;
    HugeArray<int> nodes;

    //den Header schreiben
//   p=2.0*double(num_of_links)/double(num_of_nodes)/double(num_of_nodes-1);
//   fprintf(file,"      Nodes=\t%lu\n",num_of_nodes);
//   fprintf(file,"      Links=\t%lu\n",num_of_links);
//   fprintf(file,"          q=\t%d\n",q);
//   fprintf(file,"          p=\t%f\n",p);
//   fprintf(file," Modularity=\t%f\n",calculate_Q());
//   fprintf(file,"Temperature=\t%f\n", kT);
//   fprintf(file,"Cluster\tNodes\tInnerLinks\tOuterLinks\tp_in\tp_out\t<Ln(#comm.)>\n");

    if (temperature) {
        *temperature = kT;
    }

    if (csize || membership || modularity) {
        // TODO: count the number of clusters
        for (unsigned int spin = 1; spin <= q; spin++) {
            inner_links[spin] = 0;
            outer_links[spin] = 0;
            nodes[spin] = 0;
            n_cur = iter.First(net->node_list);
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
        for (unsigned int spin = 1; spin <= q; spin++) {
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
        igraph_vector_resize(csize, 0);
        for (unsigned int spin = 1; spin <= q; spin++) {
            if (nodes[spin] > 0) {
                inner_links[spin] /= 2;
                //    fprintf(file,"Cluster\tNodes\tInnerLinks\tOuterLinks\tp_in\tp_out\n");
                /*
                N=num_of_nodes;
                n=nodes[spin];
                lin=inner_links[spin];
                lout=outer_links[spin];
                a1=N*log((double)N)-n*log((double)n)*(N-n)*log((double)N-n);
                if ((lin==long(n*(n-1)*0.5+0.5)) || (n==1)) a2=0.0;
                else a2=(n*(n-1)*0.5    )*log((double)n*(n-1)*0.5    )-(n*(n-1)*0.5    )-
                   (n*(n-1)*0.5-lin)*log((double)n*(n-1)*0.5-lin)+(n*(n-1)*0.5-lin)-
                   lin*log((double)lin            )+lin;
                */

                /*
                if ((lout==n*(N-n)) || n==N) a3=0.0;
                else a3=(n*(N-n)     )*log((double)n*(N-n)     )-(n*(N-n))-
                   (n*(N-n)-lout)*log((double)n*(N-n)-lout)+(n*(N-n)-lout)-
                   lout*log((double)lout        )+lout;
                */

                /*
                p1=(lin+lout)*log((double)p);
                p2=(0.5*n*(n-1)-lin + n*(N-n)-lout)*log((double)1.0-p);
                */
                //       fprintf(file,"%d\t%d\t%d\t%d\t%f\t%f\t%f\n",spin,nodes[spin], inner_links[spin], outer_links[spin], p_in, p_out,log_num_exp);
                IGRAPH_CHECK(igraph_vector_push_back(csize, nodes[spin]));
            }
        }
        //   fprintf(file,"\n");
    }

    //die Elemente der Cluster
    if (membership) {
        long int no = -1;
        IGRAPH_CHECK(igraph_vector_resize(membership, num_of_nodes));
        for (unsigned int spin = 1; spin <= q; spin++) {
            if (nodes[spin] > 0) {
                no++;
            }
            n_cur = iter.First(net->node_list);
            while (!iter.End()) {
                if (n_cur->Get_ClusterIndex() == spin) {
                    //         fprintf(file,"%d\t%s\n",spin,n_cur->Get_Name());
                    VECTOR(*membership)[ n_cur->Get_Index() ] = no;
                }
                n_cur = iter.Next();
            }
        }
    }

    return num_of_nodes;
}
//################################################################################################
//This function writes the soft clusters after a gamma sweep
//that is, it groups every node together that was found in
// more than threshold percent together with the other node
// in the same cluster
//################################################################################################
// Does not work at the moment !!!
//################################################################################################
// long PottsModel::WriteSoftClusters(char *filename, double threshold)
// {
//   FILE *file;
//   NNode *n_cur, *n_cur2;
//   DLList_Iter<NNode*> iter, iter2;
//   DL_Indexed_List<ClusterList<NNode*>*> *cl_list, *old_clusterlist;
//   ClusterList<NNode*> *cl_cur;

//   double max;

//   file=fopen(filename,"w");
//   if (!file) {
//     printf("Could not open %s for writing.\n",filename);
//     return -1;
//   }

//   max=correlation[0]->Get(0);
//   //printf("max=%f\n",max);
//   cl_list=new DL_Indexed_List<ClusterList<NNode*>*>();

//   n_cur=iter.First(net->node_list);
//   while (!iter.End())
//   {
//     cl_cur=new ClusterList<NNode*>();
//     cl_list->Push(cl_cur);
//     n_cur2=iter2.First(net->node_list);
//     while (!iter2.End())
//     {
//       if (double(correlation[n_cur->Get_Index()]->Get(n_cur2->Get_Index()))/max>threshold)
//         cl_cur->Push(n_cur2);
//       n_cur2=iter2.Next();
//     }
//     n_cur=iter.Next();
//   }
//   old_clusterlist=net->cluster_list;
//   net->cluster_list=cl_list;
//   clear_all_markers(net);
//   //printf("Es gibt %d Cluster\n",cl_list->Size());
//   reduce_cliques2(net, false, 15);
//   //printf("Davon bleiben %d Cluster uebrig\n",cl_list->Size());
//   clear_all_markers(net);
//   while (net->cluster_list->Size()){
//     cl_cur=net->cluster_list->Pop();
//     while (cl_cur->Size())
//     {
//       n_cur=cl_cur->Pop();
//       fprintf(file,"%s\n",n_cur->Get_Name());
//       //printf("%s\n",n_cur->Get_Name());
//     }
//     fprintf(file,"\n");
//   }
//   net->cluster_list=old_clusterlist;
//   fclose(file);

//   return 1;
// }
//#############################################################################
// Performs a gamma sweep
//#############################################################################
double PottsModel::GammaSweep(double gamma_start, double gamma_stop, double prob, unsigned int steps, bool non_parallel, int repetitions) {
    double stepsize;
    double kT = 0.5, kT_start;
    long changes;
    double gamma, acc;
    NNode *n_cur, *n_cur2;
    DLList_Iter<NNode*> iter, iter2;

    stepsize = (gamma_stop - gamma_start) / double(steps);

    n_cur = iter.First(net->node_list);
    while (!iter.End()) {
        correlation[n_cur->Get_Index()] = new HugeArray<double>();
        n_cur2 = iter2.First(net->node_list);
        while (!iter2.End()) {
            correlation[n_cur->Get_Index()]->Set(n_cur->Get_Index()) = 0.0;
            n_cur2 = iter2.Next();
        }
        n_cur = iter.Next();
    }

    for (unsigned int n = 0; n <= steps; n++) {
        assign_initial_conf(-1);
        initialize_Qmatrix();
        gamma = gamma_start + stepsize * n;
        kT = 0.5;
        acceptance = 0.5;
        while (acceptance < (1.0 - 1.0 / double(q)) * 0.95) { //wollen 95% Acceptance
            kT *= 1.1;
            //initialize_lookup(kT,kmax,net->node_list->Size());
            if (!non_parallel) {
                HeatBathParallelLookup(gamma, prob, kT, 25);
            } else {
                HeatBathLookup(gamma, prob, kT, 25);
            }
            // printf("kT=%f acceptance=%f\n", kT, acceptance);
        }
        // printf("Starting with gamma=%f\n", gamma);
        kT_start = kT;

        for (int i = 0; i < repetitions; i++) {
            changes = 1;
            kT = kT_start;
            assign_initial_conf(-1);
            initialize_Qmatrix();
            while ((changes > 0) && (kT > 0.01)) {
                kT = kT * 0.99;
                //initialize_lookup(kT,kmax,net->node_list->Size());
                if (!non_parallel) {
                    changes = HeatBathParallelLookup(gamma, prob, kT, 50);
                    // printf("kT: %f   \t Changes %li\n",kT, changes);
                } else {
                    acc = HeatBathLookup(gamma, prob, kT, 50);
                    if (acc > (1.0 - 1.0 / double(q)) * 0.01) {
                        changes = 1;
                    } else {
                        changes = 0;
                    }
                    // printf("kT: %f   Acceptance: %f\n",kT, acc);
                }
            }
            // printf("Finisched with acceptance: %1.6f bei kT=%2.4f und gamma=%2.4f\n",acceptance,kT, gamma);
//      fprintf(file,"%f\t%f\n",gamma_,acceptance);
//      fprintf(file2,"%f\t%f\n",gamma_,kT);
            //   fprintf(file3,"%f\t%d\n",gamma_,count_clusters(5));

            //Die Correlation berechnen
            n_cur = iter.First(net->node_list);
            while (!iter.End()) {
                n_cur2 = iter2.First(net->node_list);
                while (!iter2.End()) {
                    if (n_cur->Get_ClusterIndex() == n_cur2->Get_ClusterIndex()) {
                        correlation[n_cur->Get_Index()]->Set(n_cur2->Get_Index()) += 0.5;
                    }
                    n_cur2 = iter2.Next();
                }
                n_cur = iter.Next();
            }
        } // for i
    } //for n
    return kT;
}
//#############################################################################
//Performs a Gamma sweep at zero T
//#############################################################################
double PottsModel::GammaSweepZeroTemp(double gamma_start, double gamma_stop, double prob, unsigned int steps, bool non_parallel, int repetitions) {
    double stepsize;
    long changes;
    double gamma = gamma_start, acc;
    long runs;
    NNode *n_cur, *n_cur2;
    DLList_Iter<NNode*> iter, iter2;

    stepsize = (gamma_stop - gamma_start) / double(steps);

    n_cur = iter.First(net->node_list);
    while (!iter.End()) {
        correlation[n_cur->Get_Index()] = new HugeArray<double>();
        n_cur2 = iter2.First(net->node_list);
        while (!iter2.End()) {
            correlation[n_cur->Get_Index()]->Set(n_cur->Get_Index()) = 0.0;
            n_cur2 = iter2.Next();
        }
        n_cur = iter.Next();
    }

    for (unsigned int n = 0; n <= steps; n++) {
        assign_initial_conf(-1);
        initialize_Qmatrix();
        gamma = gamma_start + stepsize * n;
        // printf("Starting with gamma=%f\n", gamma);
        for (int i = 0; i < repetitions; i++) {
            changes = 1;
            assign_initial_conf(-1);
            initialize_Qmatrix();
            runs = 0;
            while (changes > 0 && runs < 250) {
                //initialize_lookup(kT,kmax,net->node_list->Size());
                if (!non_parallel) {
                    changes = HeatBathParallelLookupZeroTemp(gamma, prob, 1);
                    // printf("Changes %li\n", changes);
                } else {
                    acc = HeatBathLookupZeroTemp(gamma, prob, 1);
                    if (acc > (1.0 - 1.0 / double(q)) * 0.01) {
                        changes = 1;
                    } else {
                        changes = 0;
                    }
                    // printf("Acceptance: %f\n", acc);
                }
                runs++;
            }
            // printf("Finisched with Modularity: %1.6f bei Gamma=%1.6f\n",calculate_Q(), gamma);
//      fprintf(file,"%f\t%f\n",gamma_,acceptance);
//      fprintf(file2,"%f\t%f\n",gamma_,kT);
            //   fprintf(file3,"%f\t%d\n",gamma_,count_clusters(5));

            //Die Correlation berechnen
            n_cur = iter.First(net->node_list);
            while (!iter.End()) {
                n_cur2 = iter2.First(net->node_list);
                while (!iter2.End()) {
                    if (n_cur->Get_ClusterIndex() == n_cur2->Get_ClusterIndex()) {
                        correlation[n_cur->Get_Index()]->Set(n_cur2->Get_Index()) += 0.5;
                        correlation[n_cur2->Get_Index()]->Set(n_cur->Get_Index()) += 0.5;
                    }
                    n_cur2 = iter2.Next();
                }
                n_cur = iter.Next();
            }
        } // for i
    } //for n
    return gamma;
}
//#######################################################################
//-----------------------------------------------------------------------
//#######################################################################
// This function writes the Correlation Matrix that results from a
// Gamma-Sweep, this matrix is used to make ps files of it.
// ######################################################################
// long PottsModel::WriteCorrelationMatrix(char *filename)
// {
//   FILE *file, *file2;
//   char filename2[255];
//   NNode *n_cur, *n_cur2;
//   DLList_Iter<NNode*> iter, iter2;

//   sprintf(filename2,"%s.mat",filename);
//   file=fopen(filename,"w");
//   if (!file) {
//     printf("Could not open %s for writing.\n",filename);
//     return -1;
//   }
//   file2=fopen(filename2,"w");
//   if (!file2) {
//     printf("Could not open %s for writing.\n",filename2);
//     return -1;
//   }
//   //write the header in one line
//   n_cur=iter.First(net->node_list);
//   while (!iter.End())
//   {
//       fprintf(file, "\t%s",n_cur->Get_Name());
//       n_cur=iter.Next();
//   }
//   fprintf(file, "\n");

//   //fprintf(file, "%d\t%d\n",net->node_list->Size(),net->node_list->Size());

//   long r=0,c=0;
//   n_cur=iter.First(net->node_list);
//   while (!iter.End())
//   {
//     fprintf(file, "%s",n_cur->Get_Name());
//     r++;
//     n_cur2=iter2.First(net->node_list);
//     while (!iter2.End())
//     {
//       c++;
//       fprintf(file,"\t%f",correlation[n_cur->Get_Index()]->Get(n_cur2->Get_Index()));
//       fprintf(file2,"%li\t%li\t%f\n",r,c,correlation[n_cur->Get_Index()]->Get(n_cur2->Get_Index()));
//       n_cur2=iter2.Next();
//     }
//     fprintf(file,"\n");
//     n_cur=iter.Next();
//   }
//   fclose(file);
//   fclose(file2);
//   return 1;
// }
//##############################################################################

//#################################################################################################
PottsModelN::PottsModelN(network *n, unsigned int num_communities, bool directed) :
    degree_pos_in(NULL), degree_neg_in(NULL),
    degree_pos_out(NULL), degree_neg_out(NULL),
    degree_community_pos_in(NULL), degree_community_neg_in(NULL),
    degree_community_pos_out(NULL), degree_community_neg_out(NULL),
    csize(NULL), spin(NULL), neighbours(NULL), weights(NULL)
{
    //Set internal variable
    net = n;
    q   = num_communities;

    is_directed = directed;

    is_init = false;

    num_nodes   = net->node_list->Size();
}
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
#ifdef SPINGLASS_DEBUG
    printf("Start assigning.\n");
#endif
    unsigned int s;
    DLList_Iter<NNode*> iter;
    DLList_Iter<NLink*> l_iter;
    NNode *n_cur;
    NLink *l_cur;


    if (init_spins) {
#ifdef SPINGLASS_DEBUG
        printf("Initializing spin.\n");
#endif
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

        spin            = new unsigned int[num_nodes]; //The spin state of each node
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
    csize                       = new unsigned int[q + 1]; //The number of nodes in each community


    //Initialize communities
    for (unsigned int i = 0; i <= q; i++) {
        degree_community_pos_in[i]  = 0.0;
        degree_community_neg_in[i]  = 0.0;
        degree_community_pos_out[i] = 0.0;
        degree_community_neg_out[i] = 0.0;

        csize[i]                    = 0;
    }

    //Initialize vectors
    if (init_spins) {
        for (unsigned int i = 0; i < num_nodes; i++) {
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
    //double av_w = 0.0, av_k=0.0;
    //int l = 0;
#ifdef SPINGLASS_DEBUG
    printf("Visiting each node.\n");
#endif
    for (unsigned int v = 0; v < num_nodes; v++) {
        if (init_spins) {
            s = RNG_INTEGER(1, q);  //The new spin s
            spin[v] = (unsigned int)s;
        } else {
            s = spin[v];
        }

#ifdef SPINGLASS_DEBUG
        printf("Spin %d assigned to node %d.\n", s, v);
#endif

        n_cur               =  net->node_list->Get(v);

        l_cur               = l_iter.First(n_cur->Get_Links());

        sum_weight_pos_in   = 0.0;
        sum_weight_pos_out  = 0.0;
        sum_weight_neg_in   = 0.0;
        sum_weight_neg_out  = 0.0;

        while (!l_iter.End()) {
            double w = l_cur->Get_Weight();
            //av_w = (av_w*l + w)/(l+1); //Average weight
            //l++;
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

        //av_k = (av_k*l + sum_weight_pos_in)/(l+1); //Average k

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
    DLList_Iter<NNode*> iter;
    DLList_Iter<NLink*> l_iter;
    DLList_Iter<unsigned int*> i_iter, i_iter2;
    NNode *node, *n_cur;
    NLink *l_cur;
    /* The new_spin contains the spin to which we will update,
     * the spin_opt is the optional spin we will consider and
     * the old_spin is the spin of the node we are currently
     * changing.
     */
    unsigned int new_spin, spin_opt, old_spin;
    unsigned int sweep; //current sweep
    unsigned long changes, problemcount; //Number of changes and number of problems encountered

    double exp_old_spin; //The expectation value for the old spin
    double exp_spin; //The expectation value for the other spin(s)
    int v; //The node we will be investigating

    //The variables required for the calculations
    double delta_pos_out, delta_pos_in, delta_neg_out, delta_neg_in;
    double k_v_pos_out, k_v_pos_in, k_v_neg_out, k_v_neg_in;

    //weight of edge
    double w;

    double beta = 1 / t; //Weight for probabilities
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
        for (unsigned int n = 0; n < num_nodes; n++) {
            //Look for a random node
            v = RNG_INTEGER(0, num_nodes - 1);
            //We will be investigating node v

            node = net->node_list->Get(v);

            /*******************************************/
            // initialize the neighbours and the weights
            problemcount = 0;
            for (unsigned int i = 0; i <= q; i++) {
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
                weights[spin_opt]  = exp((double)(beta * weights[spin_opt]));
                sum_weights   += weights[spin_opt];
            }   // for spin
            /*******************************************/


            /*******************************************/
            //Choose a new spin dependent on the calculated probabilities
            r = RNG_UNIF(0, sum_weights);
            new_spin = 1;

            bool found = false;
            while (!found && new_spin <= q) {
                if (r <= weights[new_spin]) {
                    spin_opt = new_spin; //We have found are new spin
                    found = true;
                    break;
                } else {
                    r -= weights[new_spin];    //Perhaps the next spin is the one we want
                }

                new_spin++;
            }

            //Some weird thing happened. We haven't found a new spin
            //while that shouldn't be the case. Numerical problems?
            if (!found) {
                problemcount++;
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

long PottsModelN::WriteClusters(igraph_real_t *modularity,
                                igraph_real_t *temperature,
                                igraph_vector_t *community_size,
                                igraph_vector_t *membership,
                                igraph_matrix_t *adhesion,
                                igraph_matrix_t *normalised_adhesion,
                                igraph_real_t *polarization,
                                double t,
                                double d_p,
                                double d_n,
                                double gamma,
                                double lambda) {
    IGRAPH_UNUSED(gamma);
    IGRAPH_UNUSED(lambda);
#ifdef SPINGLASS_DEBUG
    printf("Start writing clusters.\n");
#endif
    //Reassign each community so that we retrieve a community assignment 1 through num_communities
    unsigned int *cluster_assign = new unsigned int[q + 1];
    for (unsigned int i = 0; i <= q; i++) {
        cluster_assign[i] = 0;
    }

    int num_clusters = 0;

    //Find out what the new communities will be
    for (unsigned int i = 0; i < num_nodes; i++) {
        unsigned int s = spin[i];
        if (cluster_assign[s] == 0) {
            num_clusters++;
            cluster_assign[s] = num_clusters;
#ifdef SPINGLASS_DEBUG
            printf("Setting cluster %d to %d.\n", s, num_clusters);
#endif
        }
    }


    /*
    DLList_Iter<NNode*> iter;
    NNode *n_cur=iter.First(net->node_list);
    n_cur = iter.First(net->node_list);
    */

    //And now assign each node to its new community
    q = num_clusters;
    for (unsigned int i = 0; i < num_nodes; i++) {
#ifdef SPINGLASS_DEBUG
        printf("Setting node %d to %d.\n", i, cluster_assign[spin[i]]);
#endif
        unsigned int s = cluster_assign[spin[i]];
        spin[i] = s;
#ifdef SPINGLASS_DEBUG
        printf("Have set node %d to %d.\n", i, s);
#endif
    }
    assign_initial_conf(false);

    delete[] cluster_assign;

    if (temperature) {
        *temperature = t;
    }

    if (community_size) {
        //Initialize the vector
        IGRAPH_CHECK(igraph_vector_resize(community_size, q));
        for (unsigned int spin_opt = 1; spin_opt <= q; spin_opt++) {
            //Set the community size
            VECTOR(*community_size)[spin_opt - 1] = csize[spin_opt];
        }
    }

    //Set the membership
    if (membership) {
        IGRAPH_CHECK(igraph_vector_resize(membership, num_nodes));
        for (unsigned int i = 0; i < num_nodes; i++) {
            VECTOR(*membership)[ i ] = spin[i] - 1;
        }
    }

    double Q = 0.0; //Modularity
    if (adhesion) {
        IGRAPH_CHECK(igraph_matrix_resize(adhesion, q, q));
        IGRAPH_CHECK(igraph_matrix_resize(normalised_adhesion, q, q));

        double **num_links_pos = NULL;
        double **num_links_neg = NULL;
        //memory allocated for elements of rows.
        num_links_pos = new double *[q + 1] ;
        num_links_neg = new double *[q + 1] ;

        //memory allocated for  elements of each column.
        for ( unsigned int i = 0 ; i < q + 1 ; i++) {
            num_links_pos[i] = new double[q + 1];
            num_links_neg[i] = new double[q + 1];
        }



        //Init num_links
        for (unsigned int i = 0; i <= q; i++) {
            for (unsigned int j = 0; j <= q; j++) {
                num_links_pos[i][j] = 0.0;
                num_links_neg[i][j] = 0.0;
            }
        }

        DLList_Iter<NLink*> iter_l;
        NLink *l_cur = iter_l.First(net->link_list);

        double w = 0.0;

        while (!iter_l.End()) {
            w = l_cur->Get_Weight();
            unsigned int a = spin[l_cur->Get_Start()->Get_Index()];
            unsigned int b =  spin[l_cur->Get_End()->Get_Index()];
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
        for (unsigned int i = 1; i <= q; i++) {
            for (unsigned int j = 1; j <= q; j++) {
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
                //printf("%f/%f %d/%d\t", num_links_pos[i][j], num_links_neg[i][j], csize[i], csize[j]);
                //printf("%f/%f - %f(%f)\t", u_p, u_n, expected, max_expected);
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
        for ( unsigned int i = 0 ; i < q + 1 ; i++ ) {
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
        for (unsigned int i = 0; i < q; i++) {
            for (unsigned int j = 0; j < q; j++) {
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
