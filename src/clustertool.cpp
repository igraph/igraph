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
                          main.cpp  -  description
                             -------------------
    begin                : Tue Jul 13 11:26:47 CEST 2004
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "NetDataTypes.h"
#include "NetRoutines.h"
#include "pottsmodel_2.h"

#include "igraph.h"
#include "error.h"
#include "random.h"

int igraph_spinglass_community(const igraph_t *graph,
			       const igraph_vector_t *weights,
			       igraph_real_t *modularity,
			       igraph_real_t *temperature,
			       igraph_vector_t *membership, 
			       igraph_vector_t *csize, 
			       igraph_integer_t spins,
			       igraph_bool_t parupdate,
			       igraph_real_t starttemp,
			       igraph_real_t stoptemp,
			       igraph_real_t coolfact,
			       igraph_spincomm_update_t update_rule,
			       igraph_real_t gamma) {

  unsigned long changes, runs;
  igraph_bool_t use_weights=0;
  bool zeroT;
  double Q, kT, acc, prob;
  ClusterList<NNode*> *cl_cur;
  network *net;
  PottsModel *pm;

  /* Check arguments */

  if (spins < 2 || spins > 500) {
    IGRAPH_ERROR("Invalid number of spins", IGRAPH_EINVAL);
  }
  if (update_rule != IGRAPH_SPINCOMM_UPDATE_SIMPLE &&
      update_rule != IGRAPH_SPINCOMM_UPDATE_CONFIG) {
    IGRAPH_ERROR("Invalid update rule", IGRAPH_EINVAL);
  }
  if (weights) {
    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
      IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }
    use_weights=1;
  }
  if (coolfact < 0 || coolfact>=1.0) {
    IGRAPH_ERROR("Invalid cooling factor", IGRAPH_EINVAL);
  }
  if (gamma < 0.0) {
    IGRAPH_ERROR("Invalid gamme value", IGRAPH_EINVAL);
  }
  if (starttemp/stoptemp<1.0) {
    IGRAPH_ERROR("starttemp should be larger in absolute value than stoptemp",
		 IGRAPH_EINVAL);
  }
  
  /* Check whether we have a single component */
  igraph_bool_t conn;
  IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_WEAK));
  if (!conn) {
    IGRAPH_ERROR("Cannot work with unconnected graph", IGRAPH_EINVAL);
  }

  net = new network;
  net->node_list   =new DL_Indexed_List<NNode*>();
  net->link_list   =new DL_Indexed_List<NLink*>();
  net->cluster_list=new DL_Indexed_List<ClusterList<NNode*>*>();

  /* Transform the igraph_t */
  IGRAPH_CHECK(igraph_i_read_network(graph, weights,
				     net, 0.0, use_weights, 0));

  prob=2.0*net->sum_weights/double(net->node_list->Size())
    /double(net->node_list->Size()-1);

  pm=new PottsModel(net,(unsigned int)spins,update_rule);

  /* initialize the random number generator */
  RNG_BEGIN();
  
  if ((stoptemp==0.0) && (starttemp==0.0)) zeroT=true; else zeroT=false;
  if (!zeroT) kT=pm->FindStartTemp(gamma, prob, starttemp); else kT=stoptemp;
  /* assign random initial configuration */
  pm->assign_initial_conf(-1);
  Q=pm->initialize_Qmatrix();
  runs=0;
  changes=1;

  while (changes>0 && (kT/stoptemp>1.0 || zeroT && runs<150)) {
    runs++;
    if (!zeroT) {
      kT*=coolfact;
      if (parupdate) { 
	changes=pm->HeatBathParallelLookup(gamma, prob, kT, 50);
      } else {
	acc=pm->HeatBathLookup(gamma, prob, kT, 50);
	if (acc<(1.0-1.0/double(spins))*0.01) {
	  changes=0; 
	} else { 
	  changes=1;
	}
      }
    } else {
      if (parupdate) { 
	changes=pm->HeatBathParallelLookupZeroTemp(gamma, prob, 50);
      } else {
	acc=pm->HeatBathLookupZeroTemp(gamma, prob, 50);
	/* less than 1 percent acceptance ratio */
	if (acc<(1.0-1.0/double(spins))*0.01) {
	  changes=0; 
	} else { 
	  changes=1;
	}
      }
    }
  } /* while loop */

  pm->WriteClusters(modularity, temperature, csize, membership, kT);

  while (net->link_list->Size()) delete net->link_list->Pop();
  while (net->node_list->Size()) delete net->node_list->Pop();
  while (net->cluster_list->Size())
    {
      cl_cur=net->cluster_list->Pop();
      while (cl_cur->Size()) cl_cur->Pop();
      delete cl_cur;
    }
  
  RNG_END();

  return 0;
}

int igraph_spinglass_my_community(const igraph_t *graph,
				  const igraph_vector_t *weights,
				  igraph_integer_t vertex,
				  igraph_vector_t *community,
				  igraph_integer_t spins,
				  igraph_bool_t parupdate,
				  igraph_real_t starttemp,
				  igraph_real_t stoptemp,
				  igraph_real_t coolfact,
				  igraph_spincomm_update_t update_rule,
				  igraph_real_t gamma) {

  igraph_bool_t use_weights=0;
  double prob;
  ClusterList<NNode*> *cl_cur;
  network *net;
  PottsModel *pm;
  char startnode[255];

  /* Check arguments */

  if (spins < 2 || spins > 500) {
    IGRAPH_ERROR("Invalid number of spins", IGRAPH_EINVAL);
  }
  if (update_rule != IGRAPH_SPINCOMM_UPDATE_SIMPLE &&
      update_rule != IGRAPH_SPINCOMM_UPDATE_CONFIG) {
    IGRAPH_ERROR("Invalid update rule", IGRAPH_EINVAL);
  }
  if (weights) {
    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
      IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }
    use_weights=1;
  }
  if (coolfact < 0 || coolfact>=1.0) {
    IGRAPH_ERROR("Invalid cooling factor", IGRAPH_EINVAL);
  }
  if (gamma < 0.0) {
    IGRAPH_ERROR("Invalid gamme value", IGRAPH_EINVAL);
  }
  if (starttemp/stoptemp<1.0) {
    IGRAPH_ERROR("starttemp should be larger in absolute value than stoptemp",
		 IGRAPH_EINVAL);
  }
  
  /* Check whether we have a single component */
  igraph_bool_t conn;
  IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_WEAK));
  if (!conn) {
    IGRAPH_ERROR("Cannot work with unconnected graph", IGRAPH_EINVAL);
  }

  net = new network;
  net->node_list   =new DL_Indexed_List<NNode*>();
  net->link_list   =new DL_Indexed_List<NLink*>();
  net->cluster_list=new DL_Indexed_List<ClusterList<NNode*>*>();

  /* Transform the igraph_t */
  IGRAPH_CHECK(igraph_i_read_network(graph, weights,
				     net, 0.0, use_weights, 0));

  prob=2.0*net->sum_weights/double(net->node_list->Size())
    /double(net->node_list->Size()-1);

  pm=new PottsModel(net,(unsigned int)spins,update_rule);

  /* initialize the random number generator */
  RNG_BEGIN();

  /* to be exected, if we want to find the community around a particular node*/
  /* the initial conf is needed, because otherwise, 
     the degree of the nodes is not in the weight property, stupid!!! */
  pm->assign_initial_conf(-1);
  snprintf(startnode, 255, "%li", (long int)vertex+1);
  pm->FindCommunityFromStart(gamma, prob, startnode, community);
  
  while (net->link_list->Size()) delete net->link_list->Pop();
  while (net->node_list->Size()) delete net->node_list->Pop();
  while (net->cluster_list->Size())
    {
      cl_cur=net->cluster_list->Pop();
      while (cl_cur->Size()) cl_cur->Pop();
      delete cl_cur;
    }
  
  RNG_END();

  return 0;
}
