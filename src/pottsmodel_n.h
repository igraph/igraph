/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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

/* The original version of this file was written by Jˆrg Reichardt 
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

#ifndef POTTSMODEL_HN
#define POTTSMODEL_HN

#include "NetDataTypes.h"

#include "igraph.h"

#define qmax 500

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
