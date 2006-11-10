/* -*- mode: C -*-  *//* 
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
			       igraph_real_t *myq,
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

// int main(int argc, char *argv[])
// {
//   char infile[255], filename[255], *name;
  char *name;
  bool error=false;
  int q, r, s, cl_fileindex=1;;
  unsigned long dummy, net_nodes, net_links, net_components, changes, runs;
  bool custom_q,custom_c,custom_ts,custom_tm,custom_g, custom_gs, custom_ge, use_weights;
  bool custom_r, custom_s, custom_th, custom_l, custom_np, sweep, cluster, zeroT, findComm;
  bool write_always;
  int operation_mode;
  float c, ts=1.0, tm=0.01, g, gs=0.0, ge=1.0, th=0.5, l=0.0;
  double Q, kT, acc, prob;
  time_t start, end;
  ClusterList<NNode*> *cl_cur;
  network *net;
  char startnode[255];
  PottsModel *pm;

  net = new network;
  net->node_list   =new DL_Indexed_List<NNode*>();
  net->link_list   =new DL_Indexed_List<NLink*>();
  net->cluster_list=new DL_Indexed_List<ClusterList<NNode*>*>();

  custom_q=custom_c=custom_ts=custom_tm=custom_g=use_weights=false;
  custom_gs=custom_ge=custom_r=custom_s=custom_th=custom_l=sweep=cluster=false;
  write_always=false;
  findComm=false;
  operation_mode=0;
  custom_np=true;

  // Transfer the igraph parameters into original clustertool parameters

  if (spins < 2 || spins > 500) {
    IGRAPH_ERROR("Invalid number of spins", IGRAPH_EINVAL);
  }
  q=(int)spins; custom_q=true;
  if (update_rule != IGRAPH_SPINCOMM_UPDATE_SIMPLE &&
      update_rule != IGRAPH_SPINCOMM_UPDATE_CONFIG) {
    IGRAPH_ERROR("Invalid update rule", IGRAPH_EINVAL);
  }
  operation_mode=update_rule;
  if (parupdate) {
    custom_np=false;
  }
  if (weights) {
    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
      IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }
    use_weights=true;
  }
  write_always=false;
  if (coolfact < 0 || coolfact>=1.0) {
    IGRAPH_ERROR("Invalid cooling factor", IGRAPH_EINVAL);
  }
  c=coolfact; custom_c=true; cluster=true;
  findComm=false;		// This will be a separate function (?)
  custom_l=false;		// No weight limit, remove edges before instead
  sweep=false;			// No sweeps, do it "by hand" instead
  ts=starttemp; custom_ts=true;
  tm=stoptemp; custom_tm=true;
  if (gamma < 0.0) {
    IGRAPH_ERROR("Invalid gamme value", IGRAPH_EINVAL);
  }
  g=gamma; custom_g=true;
  if (ts/tm<1.0) {
    IGRAPH_ERROR("starttemp should be larger in absolute value than stoptemp",
		 IGRAPH_EINVAL);
  }
  
  if (!sweep && !cluster && !findComm) cluster=true;  

  if (!custom_l)  l=0.0;

  /* Check whether we have a single component */
  igraph_bool_t conn;
  IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_WEAK));
  if (!conn) {
    IGRAPH_ERROR("Cannot work with unconnected graph", IGRAPH_EINVAL);
  }

  // Transform the igraph_t  
  // TODO: what is the _bi version??? bipartite graph?
  IGRAPH_CHECK(igraph_i_read_network(graph, weights,
				     net, 0.0, use_weights, 0));
/*   read_network_mtx("/tmp/edgelist.mtx",net,l,use_weights); */
  
  /* TODO: what is this? */
  /* write_remaining_clusters(net, -1, filename); */
/*   name=new char[25]; */
/*   net_components=find_components(net, net_nodes, net_links, &name[0]);  */
/*   clear_all_markers(net); */
/*   reduce_to_component(net,name); */
/*   clear_all_markers(net); */

  prob=2.0*net->sum_weights/double(net->node_list->Size())/double(net->node_list->Size()-1);

  if (!custom_q)  q=25;
  if (!custom_c)  c=0.99;
  if (!custom_ts) ts=1.0;
  if (!custom_tm) tm=0.01;
  if (!custom_g)  g=1.0;
  
  if (!custom_gs) gs=g; else gs*=g; //we assume repetion at gamma=p
  if (!custom_ge) ge=g; else ge*=g;
  if (!custom_th) th=0.5;
  if (!custom_r)  r=20;
  if (!custom_s)  s=10;
  
  pm=new PottsModel(net,q,operation_mode);
  //initialize the random number generator
  RNG_BEGIN();

  // TOOD findComm will be a separate function
    // to be exected, if we want to find the community around a particular node
/*     if (findComm) { */
/*       sprintf(filename,"%s_%s_%s",infile,startnode,"com.txt"); */

      //the initial conf is needed, because otherwise, the degree of the nodes is not in the weight property, stupid!!!
/*       pm->assign_initial_conf(-1); */
/*       pm->FindCommunityFromStart(g,prob,startnode); */
/*     } */

/*   printf("Starting with q=%d, l=%f, c=%f, ts=%f, tm=%f, gamma=%f\n",q,l,c,ts,tm,g); */

  if ((tm==0.0) && (ts==0.0)) zeroT=true; else zeroT=false;
/*   time(&start); */
  if (!zeroT) kT=pm->FindStartTemp(g,prob, ts); else kT=tm;
  //assign random initial configuration
  pm->assign_initial_conf(-1);
  Q=pm->initialize_Qmatrix();
  changes=net->node_list->Size();
  runs=0;
  changes=1;
  while ((changes>0) && ((kT/tm>1.0) || (zeroT && runs<150)))
    {
      runs++;
      if (!zeroT) {
	kT*=c;
	if (!custom_np) changes=pm->HeatBathParallelLookup(g, prob, kT, 50);
	else {
	  acc=pm->HeatBathLookup(g, prob, kT, 50);
	  if (acc<(1.0-1.0/double(q))*0.01) changes=0; else changes=1;
	}
      } else if (!custom_np) changes=pm->HeatBathParallelLookupZeroTemp(g,prob, 50);
      else {
	acc=pm->HeatBathLookupZeroTemp(g,prob, 50);
	//less than 1 percent acceptance ratio
	if (acc<(1.0-1.0/double(q))*0.01) changes=0; else changes=1;
      }
      
/*       for (int spin=1; spin<=q; spin++) */
/*         { */
/* 	  printf("%lu ",long(pm->color_field[spin])); */
/*         } */
/*       printf("\n"); */
/*       Q=pm->calculate_Q(); */
/*       if (!custom_np) printf("kT=%f\tQ=%f\t Changes: %d\n",kT, Q, changes); */
/*       else printf("kT=%f\tQ=%f\t acceptance: %2.3f\n",kT, Q, acc); */
/*       if (write_always) { */
/* 	printf("Writing Clusters....\n"); */
/* 	cl_fileindex=(cl_fileindex +1) % 2; */
/* 	  sprintf(filename,"%s_CLUSTERES_%d.txt",infile,cl_fileindex); */
/*       pm->WriteClusters(myq, modularity, temperature, csize, membership, kT); */
/*       } */
    }
/*   time(&end); */
/*   printf("Time needed %f\n", difftime(end,start)); */
/*       if (!write_always) { */
/* 	printf("Writing Clusters....\n"); */
/* 	sprintf(filename,"%s_%s",infile,"CLUSTERS.txt"); */
/* 	pm->WriteClusters(filename, kT); */
/* 	printf("Done.\n"); */
/*       } */
/*     } */
    // to be executed in case of sweep
/*     if (sweep) { */
/*         time(&start); */
/*         printf("Sweeping at zero T with q=%d, l=%f, th=%f from g=%f to g=%f in %d steps with %d repetitions\n",q,l,th,gs,ge,s,r); */
/*         pm->GammaSweepZeroTemp(gs, ge, prob, s, custom_np, r); */
/*         sprintf(filename,"%s_%s",infile,"GammaSweep.dat"); */
/*         pm->WriteCorrelationMatrix(filename); */
/*         sprintf(filename,"%s_%s",infile,"SoftClusters.txt"); */
/*         //pm->WriteSoftClusters(filename,th); */
/*         time(&end); */
/*         printf("Time needed %f\n", difftime(end,start)); */
/*     } */

  pm->WriteClusters(myq, modularity, temperature, csize, membership, kT);

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
