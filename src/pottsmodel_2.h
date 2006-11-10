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
    double FindStartTemp(double gamma, double prob,  double ts);
    long   HeatBathParallelLookupZeroTemp(double gamma, double prob, unsigned int max_sweeps);
    double HeatBathLookupZeroTemp(double gamma, double prob, unsigned int max_sweeps);
    long   HeatBathParallelLookup(double gamma, double prob, double kT, unsigned int max_sweeps);
    double HeatBathLookup(double gamma, double prob, double kT, unsigned int max_sweeps);
    double GammaSweep(double gamma_start, double gamma_stop, double prob, unsigned int steps, bool non_parallel=true, int repetitions=1);
    double GammaSweepZeroTemp(double gamma_start, double gamma_stop, double prob, unsigned int steps, bool non_parallel=true, int repetitions=1);
    long   WriteCorrelationMatrix(char *filename);
    double calculate_energy(double gamma);
    long   WriteClusters(igraph_real_t *q, igraph_real_t *modularity,
			 igraph_real_t *temperature, 
			 igraph_vector_t *csize, igraph_vector_t *membership,
			 double kT);
    long   WriteSoftClusters(char *filename, double threshold);
    double Get_Energy(void) { return energy;}
    double FindCommunityFromStart(double gamma, double prob, char *nodename);
};

#endif
