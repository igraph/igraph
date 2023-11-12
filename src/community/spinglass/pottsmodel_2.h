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

// Simple matrix class with heap allocation, allowing mat[i][j] indexing.
class SimpleMatrix {
    double *data;
    const size_t n;

public:
    explicit SimpleMatrix(size_t n_) : n(n_) { data = new double[n*n]; }
    SimpleMatrix(const SimpleMatrix &) = delete;
    ~SimpleMatrix() { delete [] data; }

    // Return a pointer to the i'th column, which can be indexed into using a second [] operator.
    // We assume column-major storage.
    double *operator [] (size_t i) { return &(data[n*i]); }
};

class PottsModel {
private:
    //these lists are needed to keep track of spin states for parallel update mode
    DL_Indexed_List<igraph_integer_t*> new_spins;
    DL_Indexed_List<igraph_integer_t*> previous_spins;

    HugeArray<HugeArray<double>*> correlation;
    network *net;
    igraph_integer_t q;
    unsigned int operation_mode;
    SimpleMatrix Qmatrix;
    double* Qa;
    double* weights;
    double total_degree_sum;
    igraph_integer_t num_of_nodes;
    igraph_integer_t num_of_links;
    igraph_integer_t k_max = 0;
    double acceptance = 0;
    double* neighbours;
    double* color_field;
public:
    PottsModel(network *net, igraph_integer_t q, int norm_by_degree);
    ~PottsModel();

    igraph_integer_t assign_initial_conf(igraph_integer_t spin);

    double initialize_Qmatrix();
    double calculate_Q();

    double FindStartTemp(double gamma, double prob,  double ts);
    igraph_integer_t HeatBathParallelLookupZeroTemp(double gamma, double prob, unsigned int max_sweeps);
    double HeatBathLookupZeroTemp(double gamma, double prob, unsigned int max_sweeps);
    igraph_integer_t HeatBathParallelLookup(double gamma, double prob, double kT, unsigned int max_sweeps);
    double HeatBathLookup(double gamma, double prob, double kT, unsigned int max_sweeps);

    igraph_integer_t WriteClusters(igraph_real_t *modularity,
                                   igraph_real_t *temperature,
                                   igraph_vector_int_t *csize, igraph_vector_int_t *membership,
                                   double kT, double gamma) const;

    double FindCommunityFromStart(double gamma, const char *nodename,
                                  igraph_vector_int_t *result,
                                  igraph_real_t *cohesion,
                                  igraph_real_t *adhesion,
                                  igraph_integer_t *inner_links,
                                  igraph_integer_t *outer_links) const;
};


class PottsModelN {
private:
    HugeArray<HugeArray<double>*> correlation;
    network *net;

    igraph_integer_t q; //number of communities
    double m_p; //number of positive ties (or sum of degrees), this equals the number of edges only if it is undirected and each edge has a weight of 1
    double m_n; //number of negative ties (or sum of degrees)
    igraph_integer_t num_nodes; //number of nodes
    bool is_directed;

    bool is_init = false;

    double *degree_pos_in = nullptr; //Postive indegree of the nodes (or sum of weights)
    double *degree_neg_in = nullptr; //Negative indegree of the nodes (or sum of weights)
    double *degree_pos_out = nullptr; //Postive outdegree of the nodes (or sum of weights)
    double *degree_neg_out = nullptr; //Negative outdegree of the nodes (or sum of weights)

    double *degree_community_pos_in = nullptr; //Positive sum of indegree for communities
    double *degree_community_neg_in = nullptr; //Negative sum of indegree for communities
    double *degree_community_pos_out = nullptr; //Positive sum of outegree for communities
    double *degree_community_neg_out = nullptr; //Negative sum of outdegree for communities

    igraph_integer_t *csize = nullptr; //The number of nodes in each community
    igraph_integer_t *spin = nullptr; //The membership of each node

    double *neighbours = nullptr; //Array of neighbours of a vertex in each community
    double *weights = nullptr; //Weights of all possible transitions to another community

public:
    PottsModelN(network *n, igraph_integer_t num_communities, bool directed);
    ~PottsModelN();
    void assign_initial_conf(bool init_spins);
    double FindStartTemp(double gamma, double lambda, double ts);
    double HeatBathLookup(double gamma, double lambda, double t, unsigned int max_sweeps);
    igraph_integer_t WriteClusters(igraph_real_t *modularity,
                       igraph_real_t *temperature,
                       igraph_vector_int_t *community_size,
                       igraph_vector_int_t *membership,
                       igraph_matrix_t *adhesion,
                       igraph_matrix_t *normalised_adhesion,
                       igraph_real_t *polarization,
                       double t,
                       double d_p,
                       double d_n);
};

#endif
