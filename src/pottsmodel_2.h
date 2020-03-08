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
                          pottsmodel.h  -  description
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

#ifndef POTTSMODEL_H
#define POTTSMODEL_H

#include "NetDataTypes.h"

#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"

#define qmax 500

class PottsModel {
private:
    //  HugeArray<double> neg_gammalookup;
    //  HugeArray<double> pos_gammalookup;
    DL_Indexed_List<unsigned int*> *new_spins;
    DL_Indexed_List<unsigned int*> *previous_spins;
    HugeArray<HugeArray<double>*> correlation;
    network *net;
    unsigned int q;
    unsigned int operation_mode;
    FILE *Qfile, *Magfile;
    double Qmatrix[qmax + 1][qmax + 1];
    double* Qa;
    double* weights;
    double total_degree_sum;
    unsigned long num_of_nodes;
    unsigned long num_of_links;
    unsigned long k_max;
    double energy;
    double acceptance;
    double *neighbours;
public:
    PottsModel(network *net, unsigned int q, int norm_by_degree);
    ~PottsModel();
    double* color_field;
    unsigned long assign_initial_conf(int spin);
    unsigned long initialize_lookup(double kT, double gamma);
    double initialize_Qmatrix(void);
    double calculate_Q(void);
    double calculate_genQ(double gamma);
    double FindStartTemp(double gamma, double prob,  double ts);
    long   HeatBathParallelLookupZeroTemp(double gamma, double prob, unsigned int max_sweeps);
    double HeatBathLookupZeroTemp(double gamma, double prob, unsigned int max_sweeps);
    long   HeatBathParallelLookup(double gamma, double prob, double kT, unsigned int max_sweeps);
    double HeatBathLookup(double gamma, double prob, double kT, unsigned int max_sweeps);
    double GammaSweep(double gamma_start, double gamma_stop, double prob, unsigned int steps, bool non_parallel = true, int repetitions = 1);
    double GammaSweepZeroTemp(double gamma_start, double gamma_stop, double prob, unsigned int steps, bool non_parallel = true, int repetitions = 1);
    long   WriteCorrelationMatrix(char *filename);
    double calculate_energy(double gamma);
    long   WriteClusters(igraph_real_t *modularity,
                         igraph_real_t *temperature,
                         igraph_vector_t *csize, igraph_vector_t *membership,
                         double kT, double gamma);
    long   WriteSoftClusters(char *filename, double threshold);
    double Get_Energy(void) {
        return energy;
    }
    double FindCommunityFromStart(double gamma, double prob, char *nodename,
                                  igraph_vector_t *result,
                                  igraph_real_t *cohesion,
                                  igraph_real_t *adhesion,
                                  igraph_integer_t *inner_links,
                                  igraph_integer_t *outer_links);
};


class PottsModelN {
private:
    //  HugeArray<double> neg_gammalookup;
    //  HugeArray<double> pos_gammalookup;
    DL_Indexed_List<unsigned int*> *new_spins;
    DL_Indexed_List<unsigned int*> *previous_spins;
    HugeArray<HugeArray<double>*> correlation;
    network *net;

    unsigned int q; //number of communities
    double m_p; //number of positive ties (or sum of degrees), this equals the number of edges only if it is undirected and each edge has a weight of 1
    double m_n; //number of negative ties (or sum of degrees)
    unsigned int num_nodes; //number of nodes
    bool is_directed;

    bool is_init;

    double *degree_pos_in; //Postive indegree of the nodes (or sum of weights)
    double *degree_neg_in; //Negative indegree of the nodes (or sum of weights)
    double *degree_pos_out; //Postive outdegree of the nodes (or sum of weights)
    double *degree_neg_out; //Negative outdegree of the nodes (or sum of weights)

    double *degree_community_pos_in; //Positive sum of indegree for communities
    double *degree_community_neg_in; //Negative sum of indegree for communities
    double *degree_community_pos_out; //Positive sum of outegree for communities
    double *degree_community_neg_out; //Negative sum of outdegree for communities

    unsigned int *csize; //The number of nodes in each community
    unsigned int *spin; //The membership of each node

    double *neighbours; //Array of neighbours of a vertex in each community
    double *weights; //Weights of all possible transitions to another community

public:
    PottsModelN(network *n, unsigned int num_communities, bool directed);
    ~PottsModelN();
    void assign_initial_conf(bool init_spins);
    double FindStartTemp(double gamma, double lambda, double ts);
    double HeatBathLookup(double gamma, double lambda, double t, unsigned int max_sweeps);
    double HeatBathJoin(double gamma, double lambda);
    double HeatBathLookupZeroTemp(double gamma, double lambda, unsigned int max_sweeps);
    long WriteClusters(igraph_real_t *modularity,
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
                       double lambda);
};

#endif
