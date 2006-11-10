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
                          NetRoutines.cpp  -  description
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "NetRoutines.h"

#include "igraph.h"

//#############################################################################
unsigned long re_mark(NNode* parent, unsigned long index)
{
  NNode* next_node;
  DLList_Iter<NNode*> iter;
  next_node=iter.First(parent->Get_Neighbours());
  while (!(iter.End()))
  {
   if (next_node->Get_Marker()>parent->Get_Marker())
    {
      index++;
  //    printf("Re_index %u auf %u\n", next_node->Get_Index(),index);
      next_node->Set_Marker(index);
      index=re_mark(next_node,index);
     }
  next_node=iter.Next();
  }
  return(index);
}
//################################################################################
unsigned long mark_tree_nodes(network *net, char* rootname)
{
  unsigned long num_of_nodes, max_degree=0, running_index, index;
  NNode *root, *n_cur;
  DLList_Iter<NNode*> iter;
  bool found;
  
  //Re-Indizieren
  //aber ueber den Marker!!

  num_of_nodes=net->node_list->Size();
  //erst mal alle marker hochsetzen
  //und die Wurzel finden
  found=false;
  n_cur=iter.First(net->node_list);
  if (!rootname) printf("No root node given!\n");
  while (!(iter.End())) {
       n_cur->Set_Marker(num_of_nodes+1);
        if (!rootname) {
          if (max_degree<n_cur->Get_Degree()) {
            root=n_cur;
            max_degree=n_cur->Get_Degree();
          }
        } else {
          if (strcmp(n_cur->Get_Name(),rootname)==0) {
            found=true;
            root=n_cur;
          }
        }
       n_cur=iter.Next();
   }
   if ((rootname) && (!found)) printf("Rootnode %s not found!\n",rootname);
   printf("Re-indexing...\n");
   running_index=0;
   if (!root) {
      root=iter.First(net->node_list);
   }
   printf("Setting Root in cluster %lu to: %s\n", root->Get_ClusterIndex(), root->Get_Name());
   root->Set_Marker(running_index);
   index=re_mark(root,root->Get_Marker());
   printf("Done\n");
   return index;
}
//##############################################################################
unsigned long make_spanning_tree(char *filename, network *net, DLList<char*> *ex_list, char* rootname, bool tree, float limit, bool use_weights,unsigned int states)
{
  FILE *file1;
  float Links;
  char S1[255],S2[255];//,Sh[15];
  int result=0;
  unsigned long index=0, max_k=0, min_k=99999999;
  double sum_weight=0.0, av_k=0.0, min_weight=1e60, max_weight=-1e60;
  unsigned long cluster_index=0;
  unsigned long num_of_clusters=0;
  unsigned long c1,c2,ch;
  ClusterList<NNode*> *cluster_item_list1,*cluster_item_list2;
  NNode *node1,*node2;
  NNode *n_cur;
  DLList_Iter<NNode*> iter;
  DLList_Iter<ClusterList<NNode*>*> c_iter;
  DLList_Iter<char*> ex_iter;
  char* ex_item;
  bool excluded;
  bool A_in_tree, B_in_tree;

  file1=fopen(filename,"r");
  if (!file1) {
    printf("Error Opening File %s for reading!\n", filename);
    return 0;
  }
  printf("Reading file %s....\n",filename);


  do
  {
    excluded=false;
    result=fscanf(file1,"%s %s %f\n",S1,S2,&Links);
   if (Links>=limit)
   {
     // In case we do not use the weights, set them to one, this will make everything else easier...
     if (!use_weights) Links=1.0;
   

     //Befindet sich einer der beiden in der Exclusion_List, gibt es die ueberhaupt??
    if (result>0 && ex_list)
    {
      ex_item=ex_iter.First(ex_list);
      while (!excluded && !ex_iter.End())
      {
       if (strcmp(ex_item,S1)==0 || strcmp(ex_item,S2)==0) excluded=true;
       ex_item=ex_iter.Next();
      }
    }     
    if (result>0 && !excluded)
      {
	if (Links>max_weight) max_weight=Links;
	if (Links<min_weight) min_weight=Links;
	sum_weight+=Links;            
	A_in_tree=B_in_tree=false;
        //printf("%s\t%s\t%d\n",S1,S2,Links);

        n_cur=iter.First(net->node_list);
        while (!(iter.End()) && !(A_in_tree && B_in_tree)) {
            if ((!strcmp(n_cur->Get_Name(), S1)) && (strlen(n_cur->Get_Name())==strlen(S1)))
              {
               A_in_tree=true;
               node1=n_cur;
              }
           if ((!strcmp(n_cur->Get_Name(), S2)) && (strlen(n_cur->Get_Name())==strlen(S2)))
              {
               B_in_tree=true;
               node2=n_cur;
              }
            n_cur=iter.Next();
        }

        if (A_in_tree)
        {// A ist im Baum zu finden
            if (B_in_tree) { // B finden wir auch
              c1=node1->Get_ClusterIndex();
              c2=node2->Get_ClusterIndex();
              if (c1!=c2)
              { //wir verbinden zwei alleinstehende Cluster...
                if (c1>c2) { ch=c1; c1=c2; c2=ch;}// strcpy(Sh,S1); strcpy(S1,S2); strcpy(S2,Sh);}
         //       printf("\tVerbinde %s in Cluster %u und %s in Cluster %u.\n", S1,c1,S2,c2);
                node1->Connect_To(node2,Links);
                //wir gehen durch die Cluster_item_list von node2
                //faerben um und schieben den Knoten in Cluster c1
                cluster_item_list1=net->cluster_list->Get(c1);
                cluster_item_list2=net->cluster_list->Get(c2);
                 //auf C1 umfaerben...
                while (cluster_item_list2->Size()) {
                       n_cur=cluster_item_list2->Pop();
                       n_cur->Set_ClusterIndex(c1);
                       cluster_item_list1->Push(n_cur);
                }
                //net->cluster_list->fDelete(cluster_item_list2);
                //delete c2;
                num_of_clusters--;
              } else {
                    //printf("%s und %s hatten wir schon im gleichen Cluster!\n", S1, S2);
		if (!tree) if (!node1->Connect_To(node2,Links)) {
		  printf("Found a double link between %s and %s\n",node1->Get_Name(), node2->Get_Name());
		  sum_weight-=Links;
		}   //if (!tree) node1->Connect_To(node2,Links);  // wenn wir ein ganzes Netz wollen, muessen wir nur diese Links hinzufuegen
              }
              
            } else { // A, aber nicht B ist im Baum zu finden
            c1=node1->Get_ClusterIndex();
 //           printf("Fuege %s zum Baum an %s in Cluster %d hinzu!\n",S2, S1, c1);
            node2 = new NNode(index,c1,net->link_list,S2,states);
            //der erste Index soll 0 sein, weil das auch die erste Feldpostion ist
            //das ist wichtig, weil sonst nicht der richtige Index an der richtigen Stelle steht!!
            //das gleiche gilt auch fuer die ClusterIndices!!
            net->node_list->Push(node2);  
            index++;
            cluster_item_list1=net->cluster_list->Get(c1);
            cluster_item_list1->Push(node2);

            node1->Connect_To(node2,Links);
	    }

        } else { // A ist nicht im Baum zu finden
            if (B_in_tree) { // nur B ist im Baum zu finden
            c2=node2->Get_ClusterIndex();
//            printf("Fuege %s zum Baum an %s in Cluster %u hinzu!\n",S1, S2, c2);
            node1 = new NNode(index,c2,net->link_list,S1,states);
            net->node_list->Push(node1);
            index++;
            cluster_item_list2=net->cluster_list->Get(c2);
            cluster_item_list2->Push(node1);

            node1->Connect_To(node2,Links);

            } else {  //keiner von beiden im Baum
    
            num_of_clusters++;
            cluster_item_list1 = net->cluster_list->Push(new ClusterList<NNode*>());
            node1 = net->node_list->Push(new NNode(index,cluster_index,net->link_list,S1,states));
            index++;
            node2 = net->node_list->Push(new NNode(index,cluster_index,net->link_list,S2,states));
            index++;
            cluster_item_list1->Push(node1);
            cluster_item_list1->Push(node2);
            cluster_index++;

            node1->Connect_To(node2,Links);

 //           printf("Neuer Cluster %u mit %s und %s!\n",cluster_index,S1,S2);

            }
        }
      }
   }  // if links>limit
  } while (result>0);
  fclose(file1);

  node1=iter.First(net->node_list);
  while (!iter.End())
  {
    if (node1->Get_Degree()>max_k) max_k=node1->Get_Degree();
    if (node1->Get_Degree()<min_k) min_k=node1->Get_Degree();
    av_k+=node1->Get_Degree();
    node1=iter.Next();
  }
  net->av_k=av_k/double(net->node_list->Size());
  net->sum_weights=sum_weight;
  net->av_weight=sum_weight/double(net->link_list->Size());
  net->min_k=min_k;
  net->max_k=max_k;
  net->min_weight=min_weight;
  net->max_weight=max_weight;
  //  printf("Sum of Link weights: %f\n",sum_weight);
  printf("%d nodes and %d links read in %d connected components.\n",net->node_list->Size(),net->link_list->Size(), num_of_clusters);
  printf("Min - Average - Maximum Connectivity: %d - %f - %d\n",net->min_k,net->av_k, net->max_k);
  printf("Min - Average - Maximum link weight : %f - %f - %f\n",net->min_weight,net->av_weight,net->max_weight);     
  return (net->node_list->Size());
}

//#################################################################################
unsigned long read_network(char *filename, network *net, float limit, DLList<char*> *exlist, bool use_weights,unsigned int states)
{
  return make_spanning_tree(filename, net, exlist, NULL, false, limit, use_weights,states);
}
//#################################################################################

int igraph_i_read_network(const igraph_t *graph, 
			  const igraph_vector_t *weights,
			  network *net, float limit,
			  igraph_bool_t use_weights, 
			  unsigned int states) {
  
  float Links;
  double av_k=0.0, sum_weight=0.0, min_weight=1e60, max_weight=-1e60;
  unsigned long min_k=999999999, max_k=0;
  long max_index=0;
  char name[255];
  int result=0;
  NNode *node1,*node2;
  DLList_Iter<NNode*> iter;
  igraph_vector_t edgelist;
  long int no_of_nodes=(long int)igraph_vcount(graph);
  long int no_of_edges=(long int)igraph_ecount(graph);
  long int ii;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edgelist, no_of_edges*2);
  IGRAPH_CHECK(igraph_get_edgelist(graph, &edgelist, 0 /* rowwise */));

  for (ii=0; ii<no_of_edges; ii++) {
    long int i1=(long int)VECTOR(edgelist)[2*ii]+1;
    long int i2=(long int)VECTOR(edgelist)[2*ii+1]+1;
    igraph_real_t Links;
    if (use_weights) {
      Links=VECTOR(*weights)[ii];
    } else {
      Links=1.0;
    }
    // From the original source
    if (max_index<i1) {
      for (int i=max_index; i<i1; i++)
	net->node_list->Push(new NNode(i,0,net->link_list,"",states));
      max_index=i1;
    }
    if (max_index<i2) {
      for (int i=max_index; i<i2; i++)
	net->node_list->Push(new NNode(i,0,net->link_list,"",states));
      max_index=i2;
    }
    
    node1=net->node_list->Get(i1-1);
    sprintf(name,"%d",i1);
    node1->Set_Name(name);
    
    node2=net->node_list->Get(i2-1);
    sprintf(name,"%d",i2);
    node2->Set_Name(name);
    
    node1->Connect_To(node2,Links);
    
    if (Links<min_weight) min_weight=Links;
    if (Links>max_weight) max_weight=Links;
    sum_weight+=Links;
  }

  IGRAPH_FINALLY_CLEAN(1);
  igraph_vector_destroy(&edgelist);
  
  node1=iter.First(net->node_list);
  while (!iter.End())
  {
    if (node1->Get_Degree()>max_k) max_k=node1->Get_Degree();
    if (node1->Get_Degree()<min_k) min_k=node1->Get_Degree();
    av_k+=node1->Get_Degree();
    node1=iter.Next();
  }
  net->av_k=av_k/double(net->node_list->Size());
  net->sum_weights=sum_weight;
  net->av_weight=sum_weight/double(net->link_list->Size());
  net->min_k=min_k;
  net->max_k=max_k;
  net->min_weight=min_weight;
  net->max_weight=max_weight;
  net->sum_bids=0;
  net->min_bids=0;
  net->max_bids=0;
  
  return 0;
}

unsigned long read_network_mtx(char *filename, network *net, float limit, bool use_weights,unsigned int states)
{

  FILE *file1;
  float Links;
  double av_k=0.0, sum_weight=0.0, min_weight=1e60, max_weight=-1e60;
  unsigned long min_k=999999999, max_k=0;
  long i1, i2, max_index=0;
  char name[255];
  int result=0;
  NNode *node1,*node2;
  DLList_Iter<NNode*> iter;

  file1=fopen(filename,"r");
  if (!file1) {
    printf("Error Opening File %s!\n",filename);
    return(0);
  }
  printf("Reading file %s.....\n",filename);
  while (0<(result=fscanf(file1,"%ld %ld %f\n",&i1,&i2,&Links)))
  {
             //haben wir schon so viele Knoten im Netz?
             //wenn nicht, dann Speicher belegen
             //der Index ist eins kleiner als der Name des Knotens, weil diese bei 1 und nicht
             //bei Null anfangen
       if (Links>=limit)
       {      
         if (!use_weights) Links=1.0;    
	 if (max_index<i1)
             {
                for (int i=max_index; i<i1; i++)
		  net->node_list->Push(new NNode(i,0,net->link_list,"",states));
                max_index=i1;
             }
             if (max_index<i2)
             {
                for (int i=max_index; i<i2; i++)
		  net->node_list->Push(new NNode(i,0,net->link_list,"",states));
                max_index=i2;
             }
                 
             node1=net->node_list->Get(i1-1);
             sprintf(name,"%d",i1);
             node1->Set_Name(name);
            
             node2=net->node_list->Get(i2-1);
             sprintf(name,"%d",i2);
             node2->Set_Name(name);
             
             node1->Connect_To(node2,Links);

	     if (Links<min_weight) min_weight=Links;
	     if (Links>max_weight) max_weight=Links;
	     sum_weight+=Links;
       }
    //  printf("Zeile gelesen\n");
  };
  fclose(file1);
  node1=iter.First(net->node_list);
  while (!iter.End())
  {
    if (node1->Get_Degree()>max_k) max_k=node1->Get_Degree();
    if (node1->Get_Degree()<min_k) min_k=node1->Get_Degree();
    av_k+=node1->Get_Degree();
    node1=iter.Next();
  }
  net->av_k=av_k/double(net->node_list->Size());
  net->sum_weights=sum_weight;
  net->av_weight=sum_weight/double(net->link_list->Size());
  net->min_k=min_k;
  net->max_k=max_k;
  net->min_weight=min_weight;
  net->max_weight=max_weight;
  net->sum_bids=0;
  net->min_bids=0;
  net->max_bids=0;

  printf("%d nodes and %d links read.\n",net->node_list->Size(),net->link_list->Size());
  printf("Min - Average - Maximum Connectivity: %d - %f - %d\n",net->min_k,net->av_k, net->max_k);
  printf("Min - Average - Maximum link weight : %f - %f - %f\n",net->min_weight,net->av_weight,net->max_weight);  
  return (net->node_list->Size());
}

//#################################################################################
unsigned long read_network_mtx_bi(char *filename, network *net, float limit, bool use_weights,unsigned int states)
{

  FILE *file1;
  float Links;
  double av_k=0.0, sum_weight=0.0, min_weight=1e60, max_weight=-1e60;
  unsigned long min_k=999999999, max_k=0, min_bids=99999999, max_bids=0, sum_bids=0;
  long i1, i2, bids1,bids2,max_index=0;
  char name[255];
  int result=0;
  NNode *node1,*node2;
  DLList_Iter<NNode*> iter;

  file1=fopen(filename,"r");
  if (!file1) {
    printf("Error Opening File %s!\n",filename);
    return(0);
  }
  printf("Reading file %s.....\n",filename);
  while (0<(result=fscanf(file1,"%ld %ld %f %d %d\n",&i1,&i2,&Links,&bids1,&bids2)))
  {
             //haben wir schon so viele Knoten im Netz?
             //wenn nicht, dann Speicher belegen
             //der Index ist eins kleiner als der Name des Knotens, weil diese bei 1 und nicht
             //bei Null anfangen
       if (Links>=limit)
       {      
         if (!use_weights) Links=1.0;    
	 if (max_index<i1)
             {
                for (int i=max_index; i<i1; i++)
		  net->node_list->Push(new NNode(i,0,net->link_list,"",states));
                max_index=i1;
             }
             if (max_index<i2)
             {
                for (int i=max_index; i<i2; i++)
		  net->node_list->Push(new NNode(i,0,net->link_list,"",states));
                max_index=i2;
             }
                 
             node1=net->node_list->Get(i1-1);
             sprintf(name,"%d",i1);
             node1->Set_Name(name);
	     if (node1->Get_Affiliations()==0) {
	       node1->Set_Affiliations(bids1);
	       sum_bids+=bids1;
	     } else {
	       if (node1->Get_Affiliations()!=bids1) {
		 printf("Error: node %s has different number of bids (%ld and %ld)!\n",node1->Get_Name(),node1->Get_Affiliations(),bids1);
	       }
	     }
            
             node2=net->node_list->Get(i2-1);
             sprintf(name,"%d",i2);
             node2->Set_Name(name);
	     if (node2->Get_Affiliations()==0) {
	       node2->Set_Affiliations(bids2);
	       sum_bids+=bids2;
	     } else {
	       if (node2->Get_Affiliations()!=bids2) {
		 printf("Error: node %s has different number of bids (%ld and %ld)!\n",node2->Get_Name(),node2->Get_Affiliations(),bids2);
	       }
	     }
             
             node1->Connect_To(node2,Links);

	     if (Links<min_weight) min_weight=Links;
	     if (Links>max_weight) max_weight=Links;

	     if (bids1<min_bids) min_bids=bids1;
	     if (bids2>max_bids) max_bids=bids1;

	     if (bids2<min_bids) min_bids=bids2;
	     if (bids2>max_bids) max_bids=bids2;


	     sum_weight+=Links;
	     
       }
    //  printf("Zeile gelesen\n");
  };
  fclose(file1);
  node1=iter.First(net->node_list);
  while (!iter.End())
  {
    if (node1->Get_Degree()>max_k) max_k=node1->Get_Degree();
    if (node1->Get_Degree()<min_k) min_k=node1->Get_Degree();
    av_k+=node1->Get_Degree();
    node1=iter.Next();
  }
  net->av_k=av_k/double(net->node_list->Size());
  net->sum_weights=sum_weight;
  net->sum_bids=sum_bids;
  net->av_weight=sum_weight/double(net->link_list->Size());
  net->min_k=min_k;
  net->max_k=max_k;
  net->min_weight=min_weight;
  net->max_weight=max_weight;
  net->min_bids=min_bids;
  net->max_bids=max_bids;
  net->av_bids=sum_bids/double(net->node_list->Size());
  printf("%d nodes and %d links read.\n",net->node_list->Size(),net->link_list->Size());
  printf("Min - Average - Maximum Connectivity: %d - %f - %d\n",net->min_k,net->av_k, net->max_k);
  printf("Min - Average - Maximum link weight : %f - %f - %f\n",net->min_weight,net->av_weight,net->max_weight);  
  printf("Min - Average - Maximum number of bids : %d - %f - %d\n",net->min_bids,net->av_bids,net->max_bids);  
  

  return (net->node_list->Size());
}

//#################################################################################
unsigned long read_marker_list(char *filename, DLList<char*> *mark_list)
{
  FILE *file;
  char *marker;
  unsigned long num_of_markers=0;
  
  file=fopen(filename,"r");
  if (!file) {
    printf("Error opening file %s for reading!\n",filename);
    return 0;
  }
  printf("Reading file %s...",filename);
  marker=new char[255];
  while (0<fscanf(file,"%s\n",marker))
  {
    num_of_markers++;
    mark_list->Push(marker);
    marker=new char[255];
  }
  printf("Done.\n    %lu marker read.\n",num_of_markers);
  return num_of_markers;
  delete marker;
}

//###################################################################################
int mark_special_nodes(DLList<NNode*> *node_list, DLList<char*> *mark_list, RGBcolor c)
{
  NNode *n_cur;
  char *marker;
  unsigned int num_of_markers=0;
  DLList_Iter<NNode*> *iter;
  DLList_Iter<char*> *c_iter;
  iter = new DLList_Iter<NNode*>();
  c_iter = new DLList_Iter<char*>();
  RGBcolor mix_color;
  
  n_cur=iter->First(node_list);
  while (!(iter->End())) {
      marker=c_iter->First(mark_list);
      while (!(c_iter->End())) {
        if (0==strcmp(n_cur->Get_Name(),(char*)marker)) {
          num_of_markers++;
          mix_color=n_cur->Get_Color();
          if (strcmp(mix_color.pajek_c,"Green")!=0) {
            mix_color.red=c.red;
            mix_color.blue=c.blue;
            mix_color.green=c.green;
            strcpy(mix_color.pajek_c,"Black");
            n_cur->Set_Color(mix_color);            
          } else {
            n_cur->Set_Color(c);
          }
        }
        marker=c_iter->Next();
      }
      n_cur=iter->Next();
  }
  delete iter;
  delete c_iter;
  printf("    %u nodes marked.\n", num_of_markers);
  return num_of_markers;
}
//#############################################################################
// siehe Newman, Girvan, Mixing patterns and community structure in networks
// oder Newman: Assortative Mixing in Networks
//##############################################################################
double assortativity(network *net)
{
  HugeArray<double> qvonk;
  HugeArray<HugeArray<double>*> excess_degree;
  double average_q=0.0, average_q_sq=0.0, variance=0.0, r, ejk;
  unsigned long max_degree=0, sum_excess_degree=0;
  unsigned long i,j,k;
  DLList_Iter<NNode*> iter;
  DLList_Iter<NLink*> l_iter;
  NLink *l_cur;
  NNode *n_cur;
  //ermitteln des maximalen Degree
  n_cur=iter.First(net->node_list);
  while (!iter.End())
  {
    if (n_cur->Get_Degree()>max_degree) max_degree=n_cur->Get_Degree();
    n_cur=iter.Next();
  }
  //aufstellen und initialisieren der excess_degree matrix
  for (i=0; i<=max_degree; i++) {
      excess_degree[i]=new HugeArray<double>();
      qvonk[i]=0.0;
      for (j=0; j<max_degree; j++) excess_degree[i]->Set(j)=0.0;
    }
  //fuellen der Matrix  
  l_cur=l_iter.First(net->link_list);
  while (!l_iter.End())
  {
    i=l_cur->Get_Start()->Get_Degree()-1;
    j=l_cur->Get_End()->Get_Degree()-1;    
    excess_degree[i]->Set(j)++;
    excess_degree[j]->Set(i)++;
    sum_excess_degree+=2;
    l_cur=l_iter.Next();
  }
  //normieren, damit die Summe ueber die Matrix 1 ergibt
  for (j=0; j<max_degree; j++) 
      for (k=0; k<max_degree; k++)
      {
        excess_degree[j]->Set(k)/=double(sum_excess_degree);
        qvonk[j]+=excess_degree[j]->Get(k);
      }
        
  //ausrechnen der Varianz von q von k
  for (k=0; k<=max_degree; k++)
  {
    average_q_sq+=double(k*k)*qvonk[k];
    average_q+=double(k)*qvonk[k];
  }
  variance=average_q_sq-average_q*average_q;
  
  //Berechnen des Assortativity-Koeffizienten
  r=0.0;
  for (j=0; j<max_degree; j++) {
      for (k=0; k<max_degree; k++)
      {
        ejk=excess_degree[j]->Get(k);
        r+=j*k*(ejk-(qvonk[j]*qvonk[k]));
      }
   }
   if (variance>1e-8) r/=variance; else printf("Varianz der Degree-Korrelation ist Null!\n");
//   printf("Assortativity-Coefficient ist: %f\n",r);
   for (i=0; i<=max_degree; i++) delete excess_degree[i];
   return r;
}
//#############################################################################
void net_stats(network *net, double &average_k, double &average_c, char *pvk_file, char *cvk_file, bool do_mail)
{  // Jetzt schauen wir uns die Statistik an
  unsigned long num_of_nodes=0;
  NNode *n_cur;
  unsigned long degree=0, max_degree=0;
  double clustering=0;
  HugeArray<int> Pvonk;
  HugeArray<double> Cvonk;
  DLList_Iter<NNode*> iter;


  average_c=0;
  average_k=0;

  //wir basteln eine List fuer P von k und C von k
  for (unsigned long i=0; i<=net->node_list->Size(); i++)
  {
    Pvonk[i]=0;
    Cvonk[i]=0.0;
  }

 // printf("%s\n","Calculating statistics and writing output files.....");
  num_of_nodes=net->node_list->Size();
 // printf("Total number of nodes: %lu\n",num_of_nodes);

  n_cur=iter.First(net->node_list);
  while (!(iter.End())) {
      degree=n_cur->Get_Degree();
      average_k+=degree;
      if (degree>max_degree) max_degree=degree;
      Pvonk[degree]++;       // die k-te stelle der Liste wird um 1 erhoeht
      if (cvk_file)
      {
          clustering=n_cur->Get_Clustering();
          average_c+=clustering;
          Cvonk[degree]+=clustering; // wir addieren den gefundenen ClusteringCoeff zur k-ten Stelle der Liste
      }                              // gemittelt wird spaeter
     // printf("%u\t%u\t%f\n",n_cur->Get_Index(),n_cur->Get_Degree(),n_cur->Get_Clustering());
      n_cur=iter.Next();
  }
  //###################################################################################
  // jetzt schreiben wir die Files uns lassen alles mit Gnuplot ausgeben.

  FILE *file1, *file2;
  if (pvk_file) {
    file1=fopen(pvk_file,"w");
    if (!file1) {
      printf("Error opening %s for writing.\n",pvk_file);
      return;
    }
  }
  if (cvk_file) {
    file2=fopen(cvk_file,"w");
    if (!file2) {
      printf("Error opening %s for writing.\n",cvk_file);
      return;
    }
  }
  
  double Pvk;
  double Cvk;
  unsigned long nvk;
  unsigned long rem_nodes;
 
  rem_nodes=num_of_nodes;          // noch uebrigen Knoten
  for (unsigned long i=0; i<=max_degree; i++)
  {
    if (Pvonk[i])              // wenn es Elemente mit Degree k gibt, dann
      {
        nvk=Pvonk[i];
        Cvk=Cvonk[i]/double(nvk);    //mittelung ueber alle Knoten mit diesem k

        Pvk=double(rem_nodes)/double(num_of_nodes);  //Kumulative Berechnung
        rem_nodes-=nvk;
        if (pvk_file) fprintf(file1,"%lu\t%e\n",i,Pvk);
        if (cvk_file) fprintf(file2,"%lu\t%e\n",i,Cvk);
      }
  }
  if (pvk_file) fclose(file1);
  if (cvk_file) fclose(file2);

//  printf("%s\n","Done!");

  average_c/=double(num_of_nodes);
  average_k/=double(num_of_nodes);
//  printf("C_average=%1.6f\n",average_c);
//  printf("k_average=%4.2f\n",average_k);
  if (do_mail) {system("mail -s\"Simulation Fertig\" reichardt@izbi.uni-leipzig.de < /home/reichardt/FertigMail.txt");}
 // if (do_pvk) {system("gnuplot test.gnu");}
 // if (do_cvk) {system("gnuplot test2.gnu");}
//  delete Pvonk;
//  delete Cvonk;
}
//#####################################################################################

unsigned long write_tulip_file(char *filename, long clusterindex, network *net, char* rootname)
{
  NNode *n_cur;
  NLink *l_cur;
  unsigned long written_nodes=0;

  DLList_Iter<NNode*> iter;
  DLList_Iter<NLink*> l_iter;
  
  
  FILE *file;
  file=fopen(filename,"w");
  if (!file) {
      printf("Error opening file %s for writing.\n",filename);
      return 0;
  }

  mark_tree_nodes(net, rootname);
  
  // Schreiben der Knoten
  fprintf(file,"(nodes");
  n_cur=iter.First(net->node_list);
  while (!(iter.End())) {
      if ((n_cur->Get_ClusterIndex()==clusterindex) || (clusterindex<0))
      {
        written_nodes++;
        fprintf(file," %lu",n_cur->Get_Marker());
      }
      n_cur=iter.Next();
  }
  fprintf(file," )\n");

  //Schreiben der Links
  l_cur=l_iter.First(net->link_list);
  unsigned long link_index=0;
  while (!(l_iter.End())) {
    if ((clusterindex<0) || ((l_cur->Get_Start()->Get_ClusterIndex()==clusterindex) && (l_cur->Get_End()->Get_ClusterIndex()==clusterindex)))
    {
       link_index++;
       if (l_cur->Get_Start()->Get_Marker()>l_cur->Get_End()->Get_Marker())
       {
         fprintf(file,"(edge %lu %lu %lu )\n",link_index,l_cur->Get_End()->Get_Marker(), l_cur->Get_Start()->Get_Marker());
       } else {
         fprintf(file,"(edge %lu %lu %lu )\n",link_index,l_cur->Get_Start()->Get_Marker(), l_cur->Get_End()->Get_Marker());
       }
    } 
    l_cur=l_iter.Next();
  }
  
 //und jetzt die Groesse der Nodes
 fprintf(file,"(property 0 size \"viewSize\"\n");
 fprintf(file,"(default \"(0.100000,0.100000,0.100000)\" \"(1.000000,1.000000,1.000000)\" )\n");
 n_cur=iter.First(net->node_list);
  while (!(iter.End())) {
     if ((n_cur->Get_ClusterIndex()==clusterindex) || (clusterindex<0))
      {
        if (n_cur->Get_Color().red==255) {
        fprintf(file,"(node  %lu \"(%1.6f,%1.6f,%1.6f)\")\n",n_cur->Get_Marker(), 0.3,0.3,0.3);
        }
      }
      n_cur=iter.Next();
  }
 fprintf(file,")\n");
 
 //und jetzt die Farben der Nodes
 fprintf(file,"(property 0 color \"viewColors\"\n");
 fprintf(file,"(default \"(0,0,0,255)\" \"(0,0,0,255)\" )\n");
 n_cur=iter.First(net->node_list);
  while (!(iter.End())) {
     if ((n_cur->Get_ClusterIndex()==clusterindex) || (clusterindex<0))
      {
        fprintf(file,"(node  %lu \"(%u,%u,%u,%u)\")\n",n_cur->Get_Marker(), n_cur->Get_Color().red,n_cur->Get_Color().green,n_cur->Get_Color().blue,255);
      }
      n_cur=iter.Next();
  }
 fprintf(file,")\n");

 //und jetzt die Namen der Nodes
 fprintf(file,"(property 0 string \"viewLabel\"\n");
 fprintf(file,"(default \"\" \"\")\n");
  n_cur=iter.First(net->node_list);
  while (!(iter.End())) {
     if (((n_cur->Get_ClusterIndex()==clusterindex) || (clusterindex<0)) && (strcmp(n_cur->Get_Name(),"")!=0))
      {
        fprintf(file,"(node  %lu \"%s\")\n",n_cur->Get_Marker(),n_cur->Get_Name());
      }
      n_cur=iter.Next();
  }
  fprintf(file,")\n");
  fprintf(file,"(displaying \n");
  fprintf(file,"(bool \"_viewLabel\" true)\n");
  fprintf(file,"(uint \"_FontsType\" 1)\n");
  fprintf(file,")\n");

  //und jetzt die Staerke der Links

  fclose(file);
  clear_all_markers(net);
  printf("Tulip-Graph with %lu nodes and %lu edges written!\n",written_nodes,link_index);
  return written_nodes;
}
//############################################################################################
unsigned long write_pajek_file(char *filename, unsigned long clusterindex, network *net)
{
  NNode *n_cur;
  NLink *l_cur;
  unsigned long written_nodes=0, new_index=0;
  DLList_Iter<NNode*> iter;
  DLList_Iter<NLink*> l_iter;

  FILE *file;
  file=fopen(filename,"w");
  if (!file) {
      printf("Error opening file %s for writing.\n",filename);
      return 0;
  }

  // Schreiben der Knoten
  if (clusterindex==0)
    fprintf(file,"*Vertices %lu %c",net->node_list->Size(),13);
  else {
     written_nodes=0;
     n_cur=iter.First(net->node_list);
     while (!iter.End()) {
       if (n_cur->Get_ClusterIndex()==clusterindex) written_nodes++;
       n_cur=iter.Next();
     }
     fprintf(file,"*Vertices %lu %c",written_nodes,13);
  }
  written_nodes=0;
  n_cur=iter.First(net->node_list);
  while (!iter.End()) {
       if ((n_cur->Get_ClusterIndex()==clusterindex) || (clusterindex==0))
      {
        written_nodes++;
        new_index++;
        n_cur->Set_Index(new_index);
        fprintf(file,"%lu \"%s\" ic %s %c",n_cur->Get_Index(),n_cur->Get_Name(),n_cur->Get_Color().pajek_c,13);
      }
      n_cur=iter.Next();
  }

  //Schreiben der Links
  if (clusterindex) fprintf(file,"*Arcs %c",13); else fprintf(file,"*Edges %c",13);
  l_cur=l_iter.First(net->link_list);
  unsigned long link_index=0;
  while (!l_iter.End()) {
    if ((clusterindex==0) || ((l_cur->Get_Start()->Get_ClusterIndex()==clusterindex) && (l_cur->Get_End()->Get_ClusterIndex()==clusterindex)))
    {
       link_index++;
       if (l_cur->Get_Start_Index()>l_cur->Get_End_Index())
       {
         fprintf(file,"%lu %lu %f %c",l_cur->Get_End_Index(), l_cur->Get_Start_Index(),l_cur->Get_Weight(),13);
       } else {
         fprintf(file,"%lu %lu %f %c",l_cur->Get_Start_Index(), l_cur->Get_End_Index(),l_cur->Get_Weight(),13);
       }
    }  
      l_cur=l_iter.Next();
  }

  fclose(file);
  printf("Pajek-Graph %s with %lu nodes and %lu edges written!\n",filename,written_nodes,link_index);
  return written_nodes;
}
//#####################################################################################
unsigned long group_clusters(DLList<ClusterList<NNode*>*> *cluster_list, DLList<ClusterList<NNode*>*> *global_cluster_list, FILE *file, unsigned int depth, long marker)
{
  unsigned long size, pos;
  char format_string[256];
  bool in_subset;
  ClusterList<NNode*> *c_cur, *largest_c;
  DLList<ClusterList<NNode*>*> *subsets;
  DLList_Iter<ClusterList<NNode*>*> c_iter, sub_iter;
  DLList_Iter<NNode*> iter;
  NNode *n_cur;

  subsets=new DLList<ClusterList<NNode*>*>();

  if (!(cluster_list->Size())) return 0;
  //wir suchen die groesste Teilmenge in der Liste der Cluster
  size=0;
  c_cur=c_iter.First(cluster_list);
  while (!(c_iter.End()))
  {
    if ((c_cur->Size()>size) && (c_cur->Get_Marker()!=marker))
    {
      size=c_cur->Size();
      largest_c=c_cur;
    }
    c_cur=c_iter.Next();
  }
  // printf("Groesster Cluster hat %u Elemente.\n",largest_c->Size());
  // Wir markieren den groessten Cluster der Clusterliste
  // Gibt es ueberhaupt noch einen grossen Cluster in der List??
  if (size) {
    largest_c->Set_Marker(marker);
    //und eroeffnen dessen Subgruppen

    c_cur=c_iter.First(global_cluster_list);
    while (!(c_iter.End()))
    {
      if ((c_cur->Get_Marker()!=marker) && (*c_cur<*largest_c))    //alle echten, unmarkierten Teilcluster von largest_c
      {
        subsets->Push(c_cur);
//      printf("Teilmenge gefunden mit %u Elementen!\n", c_cur->Size());
      }
      c_cur=c_iter.Next();
    }
    strcpy(format_string,"");
    for (unsigned int i=1; i<=depth; i++)
    {
      strcat(format_string,"\t");
    }
    //wenn wir echte Teilmengen gefunden haben, dann verfahren wir rekursiv
    if (subsets->Size())
    {
      if (depth==0) fprintf(file,"\n");
      fprintf(file,"%s[ ",format_string);
      //bevor es weitergeht, geben wir erst einmal alle elemente von largest_c largest,
      //die keinen eigenen Cluster bilden und somit nicht in subsets auftauchen
      pos=0;
      n_cur=iter.First(largest_c);
      while (!(iter.End()))
      {
        in_subset=false;
        c_cur=sub_iter.First(subsets);
        while (!(sub_iter.End()) && !in_subset)
        {
          if (c_cur->Is_In_List(n_cur)) {
            in_subset=true;
            pos++;
         }
          c_cur=sub_iter.Next();
        }
        if (!in_subset) fprintf(file,"%s ",n_cur->Get_Name());
        n_cur=iter.Next();
      }
      if (pos) fprintf(file,"\n");
       // Dann wenden wir den Algorithmus auf die verbliebenen Untergruppen an.
      while (0<group_clusters(subsets, subsets, file, depth+1,marker)) {};
      fprintf(file,"%s]\n",format_string);
    } else {   //ansonsten zeigen wir die Teilmenge
      //    printf("largest_c Size %u\n", largest_c->Size());
      if (depth==0) fprintf(file,"\n");
      fprintf(file, "%s(",format_string);
      n_cur=iter.First(largest_c);
      while (!(iter.End()))
      {
         fprintf(file,"%s",n_cur->Get_Name());
         n_cur=iter.Next();
         if (n_cur) fprintf(file," ");
      }
//    printf("Schreibe %u SubClusterElemente.\n",pos);
     fprintf(file,")\n");
    }
    // dem unmarkierten Rest der Clusterliste geht es weiter
    
    group_clusters(cluster_list,global_cluster_list, file, depth,marker);
 }
 while (subsets->Size()) subsets->Pop(); // die Subsetliste muss geloescht werden
 delete subsets;
 return size;
}
//###############################################################################################################
void reduce_cliques(DLList<ClusterList<NNode*>*> *global_cluster_list, FILE *file)
{
  unsigned long size;
  ClusterList<NNode*> *c_cur, *largest_c;
  DLList<ClusterList<NNode*>*> *subsets;
  DLList_Iter<ClusterList<NNode*>*> c_iter, sub_iter;
  DLList_Iter<NNode*> iter;
  NNode *n_cur;

  if (!(global_cluster_list->Size())) return;
  //wir suchen den groessten Cluster

  c_cur=c_iter.First(global_cluster_list);
  size=0;
  while (!(c_iter.End()))
  {
    if (c_cur->Size()>size)
      {
        size=c_cur->Size();
        largest_c=c_cur;
      }
    c_cur=c_iter.Next();
  }
// printf("Groesster Cluster hat %u Elemente.\n",largest_c->Size());

  //Schauen, ob es Teilmengen gibt, die ebenfalls gefunden wurden
  subsets=new DLList<ClusterList<NNode*>*>();
  c_cur=c_iter.First(global_cluster_list);
  while (!(c_iter.End()))
  {
    if ((*c_cur<*largest_c || *c_cur==*largest_c) && c_cur!=largest_c)   //alle echten Teilcluster von largest_c und die doppelten
    {
      subsets->Push(c_cur);
    }
    c_cur=c_iter.Next();
  }
  // die gefundenen Subsets werden aus der cluster_liste geloescht
  while (subsets->Size())
  {
    global_cluster_list->fDelete(subsets->Pop());
  }
  delete subsets;
  // Dann schreiben wir den groessten Cluster in das File
  fprintf(file,"Energie: %1.12f   Nodes:%3lu    -   ",largest_c->Get_Energy(),largest_c->Size());
  
  n_cur=iter.First(largest_c);
  while (!(iter.End()))
  {
     fprintf(file,"%s",n_cur->Get_Name());
     n_cur=iter.Next();
     if (n_cur) fprintf(file,", ");
   }
   fprintf(file,"\n");
   
   
  //Schliesslich schmeissen wir noch den eben gefundenen groessten Cluster raus
 global_cluster_list->fDelete(largest_c);
  //und dann geht es von vorn mit der Reduzierten ClusterListe los
  reduce_cliques(global_cluster_list, file);

}
//##################################################################################
void reduce_cliques2(network *net, bool only_double, long marker)
{
  unsigned long size;
  ClusterList<NNode*> *c_cur, *largest_c;
  DLList_Iter<ClusterList<NNode*>*> c_iter;
  do
  {
    //wir suchen den groessten, nicht markierten Cluster
    size=0;
    c_cur=c_iter.First(net->cluster_list);
    while (!(c_iter.End()))
    {
      if ((c_cur->Size()>size) && (c_cur->Get_Marker()!=marker))
        {
          size=c_cur->Size();
          largest_c=c_cur;
        }
      c_cur=c_iter.Next();
    }
    // printf("Groesster Cluster hat %u Elemente.\n",largest_c->Size());
    //Schauen, ob es Teilmengen gibt, die ebenfalls gefunden wurden
    c_cur=c_iter.First(net->cluster_list);
    while (!(c_iter.End()))
    {
      if (((!only_double && (*c_cur<*largest_c)) || (*c_cur==*largest_c)) && (c_cur!=largest_c))   //alle echten Teilcluster von largest_c und die doppelten
      {
        net->cluster_list->fDelete(c_cur);
        while (c_cur->Get_Candidates()->Size()) c_cur->Get_Candidates()->Pop();
        while (c_cur->Size()) c_cur->Pop();    // die knoten aber nicht loeschen!!
        delete c_cur;    // nicht vergessen, die global geloeschte Clusterliste zu loeschen
      }
      c_cur=c_iter.Next();      
    }
    //Schliesslich markieren wir noch den eben gefundenen groessten Cluster
    largest_c->Set_Marker(marker);
  } while (size);
}

//##################################################################################################
unsigned long iterate_tree_hierarchy(NNode *parent, unsigned long depth,FILE *file, bool as_list)
{
    NNode* next_node;
    char format[255];
    unsigned long newdepth, maxdepth;
    DLList_Iter<NNode*> *iter;
    maxdepth=newdepth=depth;
    iter=new DLList_Iter<NNode*>;
    next_node=iter->First(parent->Get_Neighbours());
    while (!(iter->End()))
    {
      if (next_node->Get_Marker()>parent->Get_Marker())  // wir gehen nach unten
      {
        if (!as_list) fprintf(file,"\n");
        strcpy(format,"");
        if (!as_list) for (unsigned long i=1;i<=depth;i++) strcat(format,"\t");
        if (!as_list) strcat(format,"%s"); else strcpy(format,"%s\n");
        fprintf(file,format, next_node->Get_Name());
        newdepth=iterate_tree_hierarchy(next_node,depth+1, file, as_list);
        if (maxdepth<newdepth) maxdepth=newdepth;
     }
    next_node=iter->Next();
    }
    delete iter;
    return maxdepth;
}
//##################################################################################################
unsigned long iterate_nsf_hierarchy(NNode *parent, unsigned long depth,FILE *file)
{
    NNode* next_node;
    unsigned long newdepth, maxdepth;
    bool first=true;
    DLList_Iter<NNode*> *iter;
    maxdepth=newdepth=depth;
    iter=new DLList_Iter<NNode*>;
    next_node=iter->First(parent->Get_Neighbours());
    while (!(iter->End()))
    {
      if (next_node->Get_Marker()>parent->Get_Marker())  // wir gehen nach unten
      {
        if (first) fprintf(file,",(");                 // eine Neue Klammer auf
        if (first) fprintf(file,"%s",next_node->Get_Name());  // nur vor dem ersten kein Komma
          else fprintf(file,",%s",next_node->Get_Name());     // sonst immer mit Komma
        first=false;
        newdepth=iterate_nsf_hierarchy(next_node,depth+1, file);
        if (maxdepth<newdepth) maxdepth=newdepth;
     }
    next_node=iter->Next();
    }
    if (!first) fprintf(file,")");                     //hat es ueberhaupt einen gegeben?
                                                       //dann klamer zu!
    delete iter;
    return maxdepth;
}
//----------------------------------------------------------------------------------------------
unsigned long write_tree_hierarchy(network *net,char *filename, char* rootname)
{
  FILE *file;
  NNode *root;
  unsigned long maxdepth=0;
  DLList_Iter<NNode*> iter;

  file=fopen(filename,"w");
  if (!file) {
       printf("Error opening %s for writing!\n", filename);
      return 0;
  }

  mark_tree_nodes(net, rootname);
  
  root=iter.First(net->node_list);
  while ((root->Get_Marker()!=0) && !(iter.End())) root=iter.Next();

  fprintf(file,"%s",root->Get_Name());
  maxdepth=iterate_tree_hierarchy(root,1,file, false);  // eine hierarchische List soll generiert werden
  fprintf(file,"\n");

  fclose(file);
  clear_all_markers(net);
  printf("Hierarchietiefe des Baumes: %lu Ebenen.\n",maxdepth);
  return maxdepth;
}
//----------------------------------------------------------------------------------------------
unsigned long write_tree_nsf(network *net,char *filename, char* rootname)
{
  FILE *file;
  NNode *root;
  unsigned long maxdepth=0;
  DLList_Iter<NNode*> iter;

  file=fopen(filename,"w");
  if (!file) {
       printf("Error opening %s for writing!\n", filename);
      return 0;
  }
  
  mark_tree_nodes(net, rootname);

  root=iter.First(net->node_list);
  while ((root->Get_Marker()!=0) && !(iter.End())) root=iter.Next();
  
  fprintf(file,"(%s",root->Get_Name());
  maxdepth=iterate_nsf_hierarchy(root,1,file);  // eine hierarchische List soll generiert werden
  fprintf(file,");\n");

  fclose(file);
  clear_all_markers(net);
  printf("Hierarchietiefe des Baumes: %lu Ebenen.\n",maxdepth);
  return maxdepth;
}


//---------------------------------------------------------------------------------------------------------------------
unsigned long write_subtree_list(DLList<NNode*> *node_list,char *filename, char *root_name)
{
  FILE *file;
  NNode *root;
  unsigned long maxdepth=0;
  DLList_Iter<NNode*> iter;

  file=fopen(filename,"w");
  if (!file) {
       printf("Error opening %s for writing!\n", filename);
      return 0;
  }
  root=iter.First(node_list);
  while ((strcmp(root_name,root->Get_Name())!=0) && !(iter.End())) root=iter.Next();
  if (iter.End()) return 0;

  fprintf(file,"%s\n",root->Get_Name());
  maxdepth=iterate_tree_hierarchy(root,1,file, true);  // eine richtige Liste soll erstellt werden
  fprintf(file,"\n");

  fclose(file);

  printf("Hierarchietiefe des Baumes: %lu Ebenen.\n",maxdepth);
  return maxdepth;

}


//##################################################################################################

unsigned long write_arity_file(network *net, char *filename)
{
  FILE *file;
  NNode *n_cur, *max_node;
  unsigned long max_arity;
  unsigned long sum_arity=0;
  char max_name[15];
  DLList_Iter<NNode*> iter;

  file=fopen(filename,"w");
  if (!file) {
    printf("Error opening %s for writing.\n",filename);

    return 0;
  }
  clear_all_markers(net);
  do {
  n_cur=iter.First(net->node_list);
  max_arity=0;
  strcpy(max_name,"");
  max_node=NULL;
  while (!iter.End()) {
    if (n_cur->Get_Marker()!=123) {
      if (n_cur->Get_Degree()>max_arity)  //nach Degree
      {
        max_arity=n_cur->Get_Degree();
        max_node=n_cur;
      }  else if (n_cur->Get_Degree()==max_arity)
        {
        if (0>strcmp(n_cur->Get_Name(),max_name))    //nach alfabet ordnen
        {
        max_node=n_cur;
        }
      }
    }  
    n_cur=iter.Next();
  }
  if (max_node)
  {
    max_node->Set_Marker(123);
    fprintf(file,"%s\t%lu\n",max_node->Get_Name(),max_node->Get_Degree());
    sum_arity+=max_arity;
  }
  } while (max_arity);
  fclose(file);
  clear_all_markers(net);
  printf("Haben insgesamt %lu Verbindungen\n",sum_arity);
  return sum_arity;
}
//###########################################################################
unsigned long write_remaining_clusters(network *net, long cluster_except, char *filename)
{
  FILE *file;
  ClusterList<NNode*> *c_cur;
  NNode *n_cur;
  unsigned long rem_c=0;
  DLList_Iter<ClusterList<NNode*>*> c_iter;
  DLList_Iter<NNode*> iter;
  file=fopen(filename,"w");
  if (!file) {
    printf("Error opening file %s for writing!\n", filename);
    return 0;
  }
  
  c_cur=c_iter.First(net->cluster_list);
  while (!(c_iter.End()))
  {
    if (c_cur->Size())
    {
      n_cur=iter.First(c_cur);
      if (n_cur->Get_ClusterIndex()!=cluster_except)
      {
        //fprintf(file,"[ ");
        while(!iter.End())
        {
          fprintf(file,"%s\n",n_cur->Get_Name());
          n_cur=iter.Next();
        }
        fprintf(file,"\n");
      }
      rem_c++;
    }
    c_cur=c_iter.Next();
  }
  fclose(file);
  return rem_c;
}
//###########################################################################################################
unsigned long write_subnetwork(char *filename, DL_Indexed_List<NLink*> *link_list, bool with_neighbours, RGBcolor marker_c, float limit)
{
  FILE *file;
  DLList_Iter<NLink*> *iter;
  NLink *l_cur;
  NNode *n1, *n2;
  unsigned long count=0;
  float weight;
  bool m1, m2;

  file=fopen(filename,"w");
  if (!file) {
    printf("Error opening File %s for writing.\n",filename);
    return 0;
  }
  
  iter=new DLList_Iter<NLink*>();
  l_cur=iter->First(link_list);
  while (!(iter->End()))
  {
      n1=l_cur->Get_Start();
      n2=l_cur->Get_End();
      weight=l_cur->Get_Weight();
      // entweder beide markiert oder einer und  nachbar    
      m1=((n1->Get_Color().red==marker_c.red) && (n1->Get_Color().blue==marker_c.blue) && (n1->Get_Color().green==marker_c.green));
      m2=((n2->Get_Color().red==marker_c.red) && (n2->Get_Color().blue==marker_c.blue) && (n2->Get_Color().green==marker_c.green));
      if ((weight>limit) && ((m1 && m2) || (with_neighbours && (m2 || m1))))
      {
          fprintf(file, "%s\t%s\t%f\n",n1->Get_Name(),n2->Get_Name(), weight);
          count++;
      }
      l_cur=iter->Next();
  }
  printf("Subnet %s with %lu links written.\n",filename, count);
  return count;
}
//################################################################################################################3
unsigned long write_cluster(char *filename, network *net, unsigned long cl_index)
{
  FILE *file;
  DLList_Iter<NLink*> l_iter;
  NLink *l_cur;
  NNode *n1, *n2;
  unsigned long count=0;
  float weight;

  file=fopen(filename,"w");
  if (!file) {
    printf("Error opening File %s for writing.\n",filename);
    return 0;
  }

  l_cur=l_iter.First(net->link_list);
  while (!(l_iter.End()))
  {
      n1=l_cur->Get_Start();
      n2=l_cur->Get_End();
      weight=l_cur->Get_Weight();
      if ((n1->Get_ClusterIndex()==cl_index) || (n2->Get_ClusterIndex()==cl_index))
      {
          fprintf(file, "%s %s %f %c",n1->Get_Name(),n2->Get_Name(), weight, 13);
          count++;
      }
      l_cur=l_iter.Next();
  }
  printf("Cluster mit %s with %lu links written.\n",filename, count);
  return count;
}
//################################################################################################################3
unsigned long write_tree_clustering(network *net,char *filename)
{
  FILE *file;
  NLink *l_cur, *mw_l;
  NNode *n1, *n2;
  unsigned long count=0;
  double weight, mw;
  unsigned long max_weight;
  double c1,c2;
  bool finished;
  DLList_Iter<NLink*> iter;

  file=fopen(filename,"w");
  if (!file) {
    printf("Error opening File %s for writing.\n",filename);
    return 0;
  }
  printf("Calculating the link weights for %lu links ...\n", net->link_list->Size());
  
  //max_weight ist das Gewicht, das Knoten mit einem Clustering von Null bekommen
    max_weight=net->node_list->Size()*(net->node_list->Size()-1);
  //iteriere ueber alle Links
  l_cur=iter.First(net->link_list);
  while (!(iter.End()))
  {
      n1=l_cur->Get_Start();
      n2=l_cur->Get_End();
      c1=n1->Get_Clustering();  if (c1<0.0) printf("Fehler - negatives Clustering!!\n");
      c2=n2->Get_Clustering();  if (c2<0.0) printf("Fehler - negatives Clustering!!\n");
      //entweder max_weight oder max von 1/c
      if ((c1<1e-8) || (c2<1e-8)) weight=max_weight;
      else if (c1<c2) weight=1.0/c1; else weight=1.0/c2;

      l_cur->Set_Weight((unsigned long int)weight);
      l_cur=iter.Next();
  }
  printf("Writing File ...\n");
  do
  {
    finished=true;
    mw=0.0;
    l_cur=iter.First(net->link_list);
    while (!(iter.End()))
    {
        if ((l_cur->Get_Marker()==0) && (l_cur->Get_Weight()>=mw))
        {
          mw_l=l_cur;
          mw=mw_l->Get_Weight();
          finished=false;
        }
        l_cur=iter.Next();
    }
    if (!finished) {
      fprintf(file,"%s %s %lu\n",mw_l->Get_Start()->Get_Name(),mw_l->Get_End()->Get_Name(),(unsigned long int) mw);
      mw_l->Set_Marker(10);
      count++;
    }
  } while (!finished);
  fclose(file);
  printf("Cluster_Tree_Net %s with %lu links written.\n",filename, count);
  return count;
}
//##############################################################################################################
unsigned long min(unsigned long a, unsigned long b)
{
  if (a<b) return a;
  else return b;
}

//##############################################################################################################
unsigned long write_tree_similarity(network *net,char *filename)
{
  FILE *file;
  NLink *l_cur, *mw_l;
  NNode *n1, *n2, *neighbour;
  unsigned long count=0, common_neighbours;
  
  double weight, mw;
  unsigned long max_weight;
  bool finished;
  DLList_Iter<NLink*> iter;
  DLList_Iter<NNode*> n_iter;

  file=fopen(filename,"w");
  if (!file) {
    printf("Error opening File %s for writing.\n",filename);
    return 0;
  }
  printf("Calculating the link weights for %lu links ...\n", net->link_list->Size());

  //max_weight ist das Gewicht, das Knoten mit einem Clustering von Null bekommen
    max_weight=net->node_list->Size()*(net->node_list->Size()-1);
  //iteriere ueber alle Links
  l_cur=iter.First(net->link_list);
  while (!(iter.End()))
  {
      n1=l_cur->Get_Start();
      n2=l_cur->Get_End();
      common_neighbours=1;
      neighbour=n_iter.First(n2->Get_Neighbours());    //iteriere ueber alle nachbarn von n2
      while (!n_iter.End())
      {
        if (n1->Get_Neighbours()->Is_In_List(neighbour)) common_neighbours++;
        neighbour=n_iter.Next();
      }
      
      weight=double(common_neighbours)/double(min(n1->Get_Neighbours()->Size(),n2->Get_Neighbours()->Size()))*max_weight;
      l_cur->Set_Weight((unsigned long)weight);
      l_cur=iter.Next();
  }
  printf("Writing File ...\n");
  do
  {
    finished=true;
    mw=0.0;
    l_cur=iter.First(net->link_list);
    while (!(iter.End()))
    {
        if ((l_cur->Get_Marker()==0) && (l_cur->Get_Weight()>=mw))
        {
          mw_l=l_cur;
          mw=mw_l->Get_Weight();
          finished=false;
        }
        l_cur=iter.Next();
    }
    if (!finished) {
      fprintf(file,"%s %s %lu\n",mw_l->Get_Start()->Get_Name(),mw_l->Get_End()->Get_Name(),(unsigned long int) mw);
      mw_l->Set_Marker(10);
      count++;
    }
  } while (!finished);
  fclose(file);
  printf("Similarity_Tree_Net %s with %lu links written.\n",filename, count);
  return count;
}

unsigned long find_component(char* startname, network *net, unsigned long &links, unsigned int marker)
{
  unsigned long nodes;
  DLList<NNode*> *list;
  DLList_Iter<NNode*> iter;
  DLList_Iter<NLink*> l_iter;
  NLink *l_cur;
  NNode *n_cur, *start;
  bool found=false;
  n_cur=iter.First(net->node_list);
  while (!iter.End() && !found)
  {
    if (strcmp(startname,n_cur->Get_Name())==0)
    {
      found=true;
      start=n_cur;
    }
    n_cur=iter.Next();
  }
  
  if (!found) return 0;

  list=new DLList<NNode*>();
  nodes=0;
  list->Enqueue(start);
  while (list->Size())
  {
    start=list->Dequeue();
    start->Set_Marker(marker);
    nodes++;
    n_cur=iter.First(start->Get_Neighbours());
    while (!iter.End())
    {
       if (n_cur->Get_Marker()!=marker) {
         list->Enqueue(n_cur);
         n_cur->Set_Marker(marker);
       }
       n_cur=iter.Next();
    }
  }
  
  l_cur=l_iter.First(net->link_list);
  links=0;
  while (!l_iter.End())
  {
    if (l_cur->Get_Start()->Get_Marker()==marker && l_cur->Get_End()->Get_Marker()==marker) links++;
    l_cur=l_iter.Next();
  }
//  printf("Es gibt %lu Knoten und %lu links in der Komponente.\n",nodes, links);
  delete list;
  return nodes;
}
//#################################################################################################
//Berechnet die Kuerzeste Distanz zwischen zwei Knoten per Breitensuche
//liefert -1 falls die beiden knoten nicht in der gleichen Komponente liegen
//wichtig ist die Marker wieder zurueckzusetzen und fuer jede Suche einen neuen marker zu verwenden!!
//#################################################################################################
long get_distance(NNode *start, NNode *end, unsigned long marker)
{
  DLList<NNode*> *to_go;
  DLList_Iter<NNode*> iter;
  NNode *n_cur, *n_cur2;
  bool found=false;
  long level_members;
  unsigned long dist=0;

  if (end==start) {
    return 0;
  }

  to_go=new DLList<NNode*>();

  start->Set_Marker(marker);  //diesen Knoten haben wir besucht
  to_go->Enqueue(start);
  level_members=to_go->Size();
  while (to_go->Size() && !found)  // solange noch was zu suchen ist und wir nichts gefunden haben
  {
    n_cur=to_go->Dequeue();        //nimm n_cur aus Schlange
    level_members--;               //reduziere Anzahl der Knoten auf dem Level von n_cur
    n_cur2=iter.First(n_cur->Get_Neighbours());  //packe alle nachbarn von n_cur auf Schlange,
    while (!iter.End() && !found)                //aber nur, falls noch nicht besucht
    {
      if (n_cur2->Get_Marker()!=marker)          //noch nicht besucht
        if (n_cur2!=end) {                       //und nicht gesuchter knoten
          to_go->Enqueue(n_cur2);                //dann auf Schlange
          n_cur2->Set_Marker(marker);            //und als besucht markieren
        } else {
          found=true;                       //ansonsten muss es der gesuchte sein!
          if (level_members) dist++;        //falls noch andere auf diesem Niveau uebrig sind
                                            //wird die distanz unten nicht erhoeht
                                            //deshalb diese korrektur
        }
      n_cur2=iter.Next();
    }
    // falls n_cur der letzte auf seinem level war, dann sind alle noch verbliebenen
    // Knoten in der Schlagen auf genau einem Level hoeher
    if (!level_members) {
       dist++;            //jedes Mal, wenn wir eine Stufe tiefer gestiegen sind
       level_members=to_go->Size();
    }
  }
  while (to_go->Size()) to_go->Pop();
  delete to_go;
  if (found) return dist; else return -1;
}
//###########################################################
double calculate_average_shortest_path_length(network *net, long &diameter)
{
  double av_dist=0;
  long dist, pairs=0;
  diameter=0;
  for (unsigned long i=0; i<net->node_list->Size()-1; i++)
    for (unsigned long j=i+1; j<net->node_list->Size(); j++)
    {
      dist=get_distance(net->node_list->Get(i), net->node_list->Get(j), pairs);
      if (dist>0) av_dist+=dist;
      if (dist>diameter) diameter=dist;
      pairs++;
    }
  if (pairs!=0.5*net->node_list->Size()*(net->node_list->Size()-1))
    printf("Error in Calculation of shortest average path length\n");
  av_dist/=double(pairs);
  return av_dist;
}
//#######################################################################################
unsigned long find_components(network *net, unsigned long &nodes, unsigned long &links, char *max_name)
{
  NNode *n_cur;
  DLList_Iter<NNode*> iter;
  unsigned int marker=0;
  unsigned long components=0;
  unsigned long max_nodes=0, max_links=0;

  printf("Checking for disconnected components...\n");
  n_cur=iter.First(net->node_list);
  while (!iter.End())
  {
    if (n_cur->Get_Marker()==0) {
        marker++;
        components++;
        nodes=find_component(n_cur->Get_Name(),net,links,marker);
        if (nodes>max_nodes) {
            max_nodes=nodes;
            max_links=links;
            strcpy(max_name,n_cur->Get_Name());
        }
    }
    n_cur=iter.Next();
  }
//  printf("Die groesste der %lu Komponenten hat: %lu Knoten und %lu Links\n",components, max_nodes, max_links);
  nodes=max_nodes;
  links=max_links;
  printf("Done.\n");
  return components;
}
//################################################################
void clear_all_markers(network *net)
{
  DLList_Iter<NNode*> iter;
  NNode *n_cur;
  n_cur=iter.First(net->node_list);
  while (!iter.End())
  {
    n_cur->Set_Marker(0);
    n_cur=iter.Next();
  }
}
//################################################################################
unsigned long reduce_to_component(network *net, char* startname)
{
  unsigned long new_index=0;
  double av_k=0, av_weight=0.0, min_weight=1e60, max_weight=-1e60, w;
  unsigned long marker=5, max_k=0,  min_k=999999999;
  DLList<NNode*> *list;
  DLList_Iter<NNode*> iter;
  DLList_Iter<NLink*> l_iter;
  DL_Indexed_List<NNode*> *component_nodes;
  DL_Indexed_List<NLink*> *component_links;
  NLink *l_cur;
  NNode *n_cur, *start;
  bool found=false;
  
  printf("Reducing nework to component around %s ...\n",startname);
  n_cur=iter.First(net->node_list);
  while (!iter.End() && !found)
  {
    if (strcmp(startname,n_cur->Get_Name())==0)
    {
      found=true;
      start=n_cur;
    }
    n_cur=iter.Next();
  }

  if (!found) return 0;

  component_nodes=new DL_Indexed_List<NNode*>();
  component_links=new DL_Indexed_List<NLink*>();

  list=new DLList<NNode*>();
  list->Enqueue(start);
  while (list->Size())
  {
    start=list->Dequeue();
    start->Set_Marker(marker);
    component_nodes->Push(start);
    if (start->Get_Degree()>max_k) max_k=start->Get_Degree();
    if (start->Get_Degree()<min_k) min_k=start->Get_Degree();

    av_k+=start->Get_Degree();

    start->Set_Index(new_index);
    new_index++;
    n_cur=iter.First(start->Get_Neighbours());
    while (!iter.End())
    {
       if (n_cur->Get_Marker()!=marker) {
         list->Enqueue(n_cur);
         n_cur->Set_Marker(marker);
       }
       n_cur=iter.Next();
    }
  }

  l_cur=l_iter.First(net->link_list);
  while (!l_iter.End())
  {
    if (l_cur->Get_Start()->Get_Marker()==marker && l_cur->Get_End()->Get_Marker()==marker) {
      component_links->Push(l_cur);
      w=l_cur->Get_Weight();
      if (w<min_weight) min_weight=w;
      if (w>max_weight) max_weight=w;
      av_weight+=w;
     }
    l_cur=l_iter.Next();
  }
  delete list;
  while (net->node_list->Size()) net->node_list->Pop();
  while (net->link_list->Size()) net->link_list->Pop();
  net->node_list=component_nodes;
  net->link_list=component_links;
  net->av_weight=av_weight/double(net->link_list->Size());
  net->av_k     =     av_k/double(net->node_list->Size());
  net->min_weight=min_weight;
  net->max_weight=max_weight;
  net->max_k=max_k;
  net->min_k=min_k;
  printf("Component around %s has %lu nodes and %lu links.\n",startname, net->node_list->Size(),net->link_list->Size());
  printf("Min - Average - Maximum Connectivity: %d - %f - %d\n",net->min_k,net->av_k, net->max_k);
  printf("Min - Average - Maximum link weight : %f - %f - %f\n",net->min_weight,net->av_weight,net->max_weight);
  return net->node_list->Size();
}
//#####################################################################################################
void initialize_Pajek_colors(char Pajek_Color[97][20])
{
  strcpy(Pajek_Color[1],"GreenYellow");
  strcpy(Pajek_Color[2],"Yellow");
  strcpy(Pajek_Color[3],"GoldenRod");
  strcpy(Pajek_Color[4],"Dendelion");
  strcpy(Pajek_Color[5],"Apricot");
  strcpy(Pajek_Color[6],"Peach");
  strcpy(Pajek_Color[7],"Melon");
  strcpy(Pajek_Color[8],"YellowOrange");
  strcpy(Pajek_Color[9],"Orange");
  strcpy(Pajek_Color[10],"BurntOrange");
  strcpy(Pajek_Color[11],"BitterSweet");
  strcpy(Pajek_Color[12],"RedOrange");
  strcpy(Pajek_Color[13],"Mahagony");
  strcpy(Pajek_Color[14],"Maroon");
  strcpy(Pajek_Color[15],"BrickRed");
  strcpy(Pajek_Color[16],"Red");
  strcpy(Pajek_Color[17],"OrangeRed");
  strcpy(Pajek_Color[18],"RubineRed");
  strcpy(Pajek_Color[19],"WildStrawberry");
  strcpy(Pajek_Color[20],"Salmon");
  strcpy(Pajek_Color[21],"CarnationPink");
  strcpy(Pajek_Color[22],"Magenta");
  strcpy(Pajek_Color[23],"VioletRed");
  strcpy(Pajek_Color[24],"Rhodamine");
  strcpy(Pajek_Color[25],"Mulberry");
  strcpy(Pajek_Color[26],"RedViolet");
  strcpy(Pajek_Color[27],"Gray05");
  strcpy(Pajek_Color[28],"Gray20");
  strcpy(Pajek_Color[29],"Gray35");
  strcpy(Pajek_Color[30],"Gray55");
  strcpy(Pajek_Color[31],"Gray70");
  strcpy(Pajek_Color[32],"Gray85");
  strcpy(Pajek_Color[33],"Fuchsia");
  strcpy(Pajek_Color[34],"Lavender");
  strcpy(Pajek_Color[35],"Thistle");
  strcpy(Pajek_Color[36],"Orchid");
  strcpy(Pajek_Color[37],"DarkOrchid");
  strcpy(Pajek_Color[38],"Purple");
  strcpy(Pajek_Color[39],"Plum");
  strcpy(Pajek_Color[40],"Violet");
  strcpy(Pajek_Color[41],"RoyalPurple");
  strcpy(Pajek_Color[42],"BlueViolet");
  strcpy(Pajek_Color[43],"Periwinkle");
  strcpy(Pajek_Color[44],"CadetBlue");
  strcpy(Pajek_Color[45],"CornflowerBlue");
  strcpy(Pajek_Color[46],"MidnightBlue");
  strcpy(Pajek_Color[47],"NavyBlue");
  strcpy(Pajek_Color[48],"RoyalBlue");
  strcpy(Pajek_Color[49],"Blue");
  strcpy(Pajek_Color[50],"Cerulean");
  strcpy(Pajek_Color[51],"Cyan");
  strcpy(Pajek_Color[52],"ProcessBlue");
  strcpy(Pajek_Color[53],"SkyBlue");
  strcpy(Pajek_Color[54],"Turquoise");
  strcpy(Pajek_Color[55],"TealBlue");
  strcpy(Pajek_Color[56],"Aquamarine");
  strcpy(Pajek_Color[57],"BlueGreen");
  strcpy(Pajek_Color[58],"Emerald");
  strcpy(Pajek_Color[59],"Gray10");
  strcpy(Pajek_Color[60],"Gray25");
  strcpy(Pajek_Color[61],"Gray40");
  strcpy(Pajek_Color[62],"Gray60");
  strcpy(Pajek_Color[63],"Gray75");
  strcpy(Pajek_Color[64],"Gray90");
  strcpy(Pajek_Color[65],"JungleGreen");
  strcpy(Pajek_Color[66],"SeaGreen");
  strcpy(Pajek_Color[67],"Green");
  strcpy(Pajek_Color[68],"ForestGreen");
  strcpy(Pajek_Color[69],"PineGreen");
  strcpy(Pajek_Color[70],"LimeGreen");
  strcpy(Pajek_Color[71],"YellowGreen");
  strcpy(Pajek_Color[72],"SpringGreen");
  strcpy(Pajek_Color[73],"OliveGreen");
  strcpy(Pajek_Color[74],"RawSienna");
  strcpy(Pajek_Color[75],"Sepia");
  strcpy(Pajek_Color[76],"Brown");
  strcpy(Pajek_Color[77],"Tan");
  strcpy(Pajek_Color[78],"Gray");
  strcpy(Pajek_Color[79],"Black");
  strcpy(Pajek_Color[80],"White");
  strcpy(Pajek_Color[81],"LightYellow");
  strcpy(Pajek_Color[82],"LightCyan");
  strcpy(Pajek_Color[83],"LightMagenta");
  strcpy(Pajek_Color[84],"LightPurple");
  strcpy(Pajek_Color[85],"LightGreen");
  strcpy(Pajek_Color[86],"LightOrange");
  strcpy(Pajek_Color[87],"Canary");
  strcpy(Pajek_Color[88],"LFadedGreen");
  strcpy(Pajek_Color[89],"Pink");
  strcpy(Pajek_Color[90],"LSkyBlue");
  strcpy(Pajek_Color[91],"Gray15");
  strcpy(Pajek_Color[92],"Gray30");
  strcpy(Pajek_Color[93],"Gray45");
  strcpy(Pajek_Color[94],"Gray65");
  strcpy(Pajek_Color[95],"Gray80");
  strcpy(Pajek_Color[96],"Gray95");
 }
