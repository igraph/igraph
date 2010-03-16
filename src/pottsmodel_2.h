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

/* The original version of this file was written by Jörg Reichardt 
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

#include "igraph.h"

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
    double Qmatrix[qmax+1][qmax+1];
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
    double GammaSweep(double gamma_start, double gamma_stop, double prob, unsigned int steps, bool non_parallel=true, int repetitions=1);
    double GammaSweepZeroTemp(double gamma_start, double gamma_stop, double prob, unsigned int steps, bool non_parallel=true, int repetitions=1);
    long   WriteCorrelationMatrix(char *filename);
    double calculate_energy(double gamma);
    long   WriteClusters(igraph_real_t *modularity,
			 igraph_real_t *temperature, 
			 igraph_vector_t *csize, igraph_vector_t *membership,
			 double kT, double gamma);
    long   WriteSoftClusters(char *filename, double threshold);
    double Get_Energy(void) { return energy;}
    double FindCommunityFromStart(double gamma, double prob, char *nodename,
				  igraph_vector_t *result, 
				  igraph_real_t *cohesion,
				  igraph_real_t *adhesion,
				  igraph_integer_t *inner_links,
				  igraph_integer_t *outer_links);
};

#endif
