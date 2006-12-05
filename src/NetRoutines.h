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
                          NetRoutines.h  -  description
                             -------------------
    begin                : Tue Oct 28 2003
    copyright            : (C) 2003 by Joerg Reichardt
    email                : reichardt@mitte
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef NETROUTINES_H
#define NETROUTINES_H

#include "NetDataTypes.h"
#include "igraph.h"

int igraph_i_read_network(const igraph_t *graph,
			  const igraph_vector_t *weights,
			  network *net, float limit,
			  igraph_bool_t use_weights,
			  unsigned int states);

unsigned long read_marker_list(char*, DLList<char*>*);
int mark_special_nodes(DLList<NNode*>*, DLList<char*>*, RGBcolor);
unsigned long write_tulip_file(char*, long, network* net, char* rootname=NULL);
unsigned long  write_pajek_file(char*, unsigned long,network*) ;
void net_stats(network *net, double &average_k, double &average_c, char *pvk_file=NULL, char *cvk_file=NULL, bool do_mail=false);
unsigned long group_clusters(DLList<ClusterList<NNode*>*>*, DLList<ClusterList<NNode*>*>* , FILE*, unsigned int, long);
void reduce_cliques(DLList<ClusterList<NNode*>*>*, FILE *file);
void reduce_cliques2(network*, bool,  long );
unsigned long iterate_tree_hierarchy(NNode *parent, unsigned long depth,FILE *file, bool as_list);
unsigned long mark_tree_nodes(network *net, char* rootname=NULL);
unsigned long write_subtree_list(DLList<NNode*> *node_list,char *filename, char *root_name);
unsigned long write_tree_hierarchy(network *net,char *filename, char *rootname=NULL);
unsigned long       write_tree_nsf(network *net,char *filename, char *rootname=NULL);
unsigned long write_tree_clustering(network *net,char *filename);
unsigned long write_tree_similarity(network *net,char *filename);
unsigned long write_arity_file(network *net, char *filename);
unsigned long write_remaining_clusters(network *net, long cluster_except, char *filename);
unsigned long write_subnetwork(char *filename, DL_Indexed_List<NLink*> *link_list, bool with_neighbours, RGBcolor marker_c, float limit);
unsigned long write_cluster(char *filename, network *net, unsigned long cl_index);
unsigned long find_component(char* startname, network *net, unsigned long &links, unsigned int marker);
long get_distance(NNode *start, NNode *end, unsigned long marker);
double calculate_average_shortest_path_length(network *net, long &diameter);
double assortativity(network *net);
unsigned long min(unsigned long a, unsigned long b);
unsigned long find_components(network *net, unsigned long &nodes, unsigned long &links, char *max_name);
void clear_all_markers(network *net);
unsigned long reduce_to_component(network *net, char* startname);
void initialize_Pajek_colors(char Pajek_Color[97][20]);

#endif

