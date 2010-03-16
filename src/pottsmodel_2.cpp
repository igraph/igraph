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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "pottsmodel_2.h"
#include "NetRoutines.h"

using namespace std;

#include "random.h"
#include "config.h"

//#################################################################################################
PottsModel::PottsModel(network *n, unsigned int qvalue, int m) : acceptance(0)
{
  DLList_Iter<NNode*> iter;
  NNode *n_cur;
  unsigned int *i_ptr;
  net=n;
  q=qvalue;
  operation_mode=m;
  k_max=0;
  //needed in calculating modularity
  Qa     =new double[q+1];
  //weights for each spin state needed in Monte Carlo process
  weights=new double[q+1];
  //bookkeeping of occupation numbers of spin states or the number of links in community
  color_field=new double[q+1];
  neighbours=new double[q+1];

  num_of_nodes=net->node_list->Size();
  num_of_links=net->link_list->Size();
  
  n_cur=iter.First(net->node_list);
  //these lists are needed to keep track of spin states for parallel update mode
  new_spins=new DL_Indexed_List<unsigned int*>();
  previous_spins=new DL_Indexed_List<unsigned int*>();
  while (!iter.End())
  {
    if (k_max<n_cur->Get_Degree()) k_max=n_cur->Get_Degree();
    i_ptr=new unsigned int;
    *i_ptr=0;
    new_spins->Push(i_ptr);
    i_ptr=new unsigned int;
    *i_ptr=0;
    previous_spins->Push(i_ptr);
    n_cur=iter.Next();
  }
  return;
}
//#######################################################
//Destructor of PottsModel
//########################################################
PottsModel::~PottsModel()
{
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
  return;
}
//#####################################################
//Assing an initial random configuration of spins to nodes
//if called with negative argument or the spin used as argument
//when called with positve one.
//This may be handy, if you want to warm up the network.
//####################################################
unsigned long PottsModel::assign_initial_conf(int spin)
{
  int s;
  DLList_Iter<NNode*> iter;
  DLList_Iter<NLink*> l_iter;
  NNode *n_cur;
  NLink *l_cur;
  double sum_weight;
  double av_k_squared=0.0;
  double av_k=0.0;
//   printf("Assigning initial configuration...\n");
  // initialize colorfield
  for (unsigned int i=0; i<=q; i++) color_field[i]=0.0;
  //
  total_degree_sum=0.0;
  n_cur=iter.First(net->node_list);
  while (!iter.End())
  {
    if (spin<0) s=RNG_INTEGER(1,q); else s=spin;
    n_cur->Set_ClusterIndex(s);
      l_cur=l_iter.First(n_cur->Get_Links());
      sum_weight=0;
      while (!l_iter.End())
      {
        sum_weight+=l_cur->Get_Weight();   //weight should be one, in case we are not using it.
        l_cur=l_iter.Next();
      }
      // we set the sum of the weights or the degree as the weight of the node, this way
      // we do not have to calculate it again.
      n_cur->Set_Weight(sum_weight);
      av_k_squared+=sum_weight*sum_weight;
      av_k+=sum_weight;

    // in case we want all links to be contribute equally - parameter gamm=fixed
    if (operation_mode==0) {
      color_field[s]++;
    } else {
      color_field[s]+=sum_weight;
    }
    // or in case we want to use a weight of each link that is proportional to k_i\times k_j
      total_degree_sum+=sum_weight;
    n_cur=iter.Next();
  }
  av_k_squared/=double(net->node_list->Size());
          av_k/=double(net->node_list->Size());
  // total_degree_sum-=av_k_squared/av_k;
//   printf("Total Degree Sum=2M=%f\n",total_degree_sum);
  return net->node_list->Size();
}
//#####################################################################
//If I ever manage to write a decent LookUp function, it will be here
//#####################################################################
unsigned long PottsModel::initialize_lookup(double kT, double gamma)
{ /*
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
// Q denotes the modulary of the network
// This function calculates it initially 
// In the event of a spin changing its state, it only needs updating
// Note that Qmatrix and Qa are only counting! The normalization 
// by num_of_links is done later
//####################################################################
double PottsModel::initialize_Qmatrix(void)
{
  DLList_Iter<NLink*> l_iter;
  NLink *l_cur;
  unsigned int i,j;
  //initialize with zeros
  num_of_links=net->link_list->Size();
  for (i=0; i<=q; i++)
  {
      Qa[i]=0.0;
      for (j=i; j<=q; j++) {
          Qmatrix[i][j]=0.0;
          Qmatrix[j][i]=0.0;
      }
  }
  //go over all links and make corresponding entries in Q matrix
  //An edge connecting state i wiht state j will get an entry in Qij and Qji
  l_cur=l_iter.First(net->link_list);
  while (!l_iter.End())
  {
    i=l_cur->Get_Start()->Get_ClusterIndex();
    j=l_cur->Get_End()->Get_ClusterIndex();
    //printf("%d %d\n",i,j);
    Qmatrix[i][j]+=l_cur->Get_Weight();
    Qmatrix[j][i]+=l_cur->Get_Weight();

    l_cur=l_iter.Next();
  }
  //Finally, calculate sum over rows and keep in Qa
  for (i=0; i<=q; i++)
  {
     for (j=0; j<=q; j++) Qa[i]+=Qmatrix[i][j];
  }
  return calculate_Q();
}
//####################################################################
// This function does the actual calculation of Q from the matrix 
// The normalization by num_of_links is done here
//####################################################################
double PottsModel::calculate_Q()
{
  double Q=0.0;
  for (unsigned int i=0; i<=q; i++)
  {
    Q+=Qmatrix[i][i]-Qa[i]*Qa[i]/double(2.0*net->sum_weights);
    if ((Qa[i]<0.0) || Qmatrix[i][i]<0.0) {
//         printf("Negatives Qa oder Qii\n\n\n");
        //printf("Press any key to continue\n\n");
        //cin >> Q;
    }
  }
  Q/=double(2.0*net->sum_weights);
  return Q;
}
double PottsModel::calculate_genQ(double gamma)
{
  double Q=0.0;
  for (unsigned int i=0; i<=q; i++)
  {
    Q+=Qmatrix[i][i]-gamma*Qa[i]*Qa[i]/double(2.0*net->sum_weights);
    if ((Qa[i]<0.0) || Qmatrix[i][i]<0.0) {
//         printf("Negatives Qa oder Qii\n\n\n");
        //printf("Press any key to continue\n\n");
        //cin >> Q;
    }
  }
  Q/=double(2.0*net->sum_weights);
  return Q;
}
//#######################################################################
// This function calculates the Energy for the standard Hamiltonian
// given a particular value of gamma and the current spin states
// #####################################################################
double PottsModel::calculate_energy(double gamma)
{
  double e=0.0;
  DLList_Iter<NLink*> l_iter;
  NLink *l_cur;
  l_cur=l_iter.First(net->link_list);
  //every in-cluster edge contributes -1
  while (!l_iter.End())
  {
    if (l_cur->Get_Start()->Get_ClusterIndex()==l_cur->Get_End()->Get_ClusterIndex()) e--;;
   l_cur=l_iter.Next();
  }
  //and the penalty term contributes according to cluster sizes
  for (unsigned int i=1; i<=q; i++)
  {
      e+=gamma*0.5*double(color_field[i])*double((color_field[i]-1));
  }
  energy=e;
  return e;
}
//##########################################################################
// We would like to start from a temperature with at least 95 of all proposed 
// spin changes accepted in 50 sweeps over the network
// The function returns the Temperature found
//#########################################################################
double PottsModel::FindStartTemp(double gamma, double prob, double ts)
{
  double kT;
  kT=ts;
  //assing random initial condition
  assign_initial_conf(-1);
  //initialize Modularity matrix, from now on, it will be updated at every spin change
  initialize_Qmatrix();
  // the factor 1-1/q is important, since even, at infinite temperature,
  // only 1-1/q of all spins do change their state, since a randomly chooses new
  // state is with prob. 1/q the old state.
  while (acceptance<(1.0-1.0/double(q))*0.95)      //want 95% acceptance
  {
       kT=kT*1.1;
       // if I ever have a lookup table, it will need initialization for every kT
       //initialize_lookup(kT,k_max,net->node_list->Size());
       HeatBathParallelLookup(gamma,prob, kT,50);
//        printf("kT=%f acceptance=%f\n", kT, acceptance);
  }
  kT*=1.1; // just to be sure...
//   printf("Starting with acceptance ratio: %1.6f bei kT=%2.4f\n",acceptance,kT);
  return kT;
}

//##############################################################
//This function does a parallel update at zero T
//Hence, it is really fast on easy problems
//max sweeps is the maximum number of sweeps it should perform,
//if it does not converge earlier
//##############################################################
long PottsModel::HeatBathParallelLookupZeroTemp(double gamma, double prob, unsigned int max_sweeps)
{
  DLList_Iter<NNode*> iter, net_iter;
  DLList_Iter<NLink*> l_iter;
  DLList_Iter<unsigned int*> i_iter, i_iter2;
  NNode *node, *n_cur;
  NLink *l_cur;
  unsigned int *SPIN, *P_SPIN, new_spin, spin_opt, old_spin, spin, sweep;
  // long h; // degree;
  unsigned long changes;
  double h, delta=0, deltaE, deltaEmin, w, degree;
  //HugeArray<double> neighbours;
  bool cyclic=0;
  
  sweep=0;
  changes=1;  
  while (sweep<max_sweeps && changes)
  {
    cyclic=true;
    sweep++;
    changes=0;
    //Loop over all nodes
    node=net_iter.First(net->node_list);
    SPIN=i_iter.First(new_spins);
    while (!net_iter.End())
    {
      // How many neigbors of each type?
      // set them all zero
      for (unsigned int i=0; i<=q; i++) neighbours[i]=0;
      degree=node->Get_Weight();
      //Loop over all links (=neighbours)
      l_cur=l_iter.First(node->Get_Links());
      while (!l_iter.End())
      {
        //printf("%s %s\n",node->Get_Name(),n_cur->Get_Name());
	w=l_cur->Get_Weight();
	if (node==l_cur->Get_Start()) {
	  n_cur=l_cur->Get_End();
	} else { 
	  n_cur=l_cur->Get_Start(); 
	}
        neighbours[n_cur->Get_ClusterIndex()]+=w;
        l_cur=l_iter.Next();
      }
      //Search optimal Spin      
      old_spin=node->Get_ClusterIndex();
      //degree=node->Get_Degree();
      switch (operation_mode) {
      case 0: { 
	delta=1.0; 
	break;
      } 
      case 1: { //newman modularity 
	prob=degree/total_degree_sum; 
	delta=degree; 
	break;
       }
      }


      spin_opt=old_spin;
      deltaEmin=0.0;
      for (spin=1; spin<=q; spin++)  // all possible spin states
      {
          if (spin!=old_spin)
          {
            h=color_field[spin]+delta-color_field[old_spin];
            deltaE=double(neighbours[old_spin]-neighbours[spin])+gamma*prob*double(h);
            if (deltaE<deltaEmin) {
              spin_opt=spin;
              deltaEmin=deltaE;
            }
          }
      } // for spin

     //Put optimal spin on list for later update 
     *SPIN=spin_opt;     
     node=net_iter.Next();
     SPIN=i_iter.Next();
    } // while !net_iter.End()

    //-------------------------------
    //Now set all spins to new values
    node=net_iter.First(net->node_list);
    SPIN=i_iter.First(new_spins);
    P_SPIN=i_iter2.First(previous_spins);    
    while (!net_iter.End())
    {
      old_spin=node->Get_ClusterIndex();
      new_spin=*SPIN;
      if (new_spin!=old_spin) // Do we really have a change??
      {
        changes++;
        node->Set_ClusterIndex(new_spin);
	//this is important!!
	//In Parallel update, there occur cyclic attractors of size two
	//which then make the program run for ever
        if (new_spin!=*P_SPIN) cyclic=false;
        *P_SPIN=old_spin;
        color_field[old_spin]--;
        color_field[new_spin]++;

        //Qmatrix update
        //iteration over all neighbours
        l_cur=l_iter.First(node->Get_Links());
        while (!l_iter.End())
        {
	   w=l_cur->Get_Weight();
	   if (node==l_cur->Get_Start()) {
	     n_cur=l_cur->Get_End();
	   } else { 
	     n_cur=l_cur->Get_Start(); 
	   }
          Qmatrix[old_spin][n_cur->Get_ClusterIndex()]-=w;
          Qmatrix[new_spin][n_cur->Get_ClusterIndex()]+=w;
          Qmatrix[n_cur->Get_ClusterIndex()][old_spin]-=w;
          Qmatrix[n_cur->Get_ClusterIndex()][new_spin]+=w;
          Qa[old_spin]-=w;
          Qa[new_spin]+=w;
          l_cur=l_iter.Next();
        }  // while l_iter
      }
      node=net_iter.Next();
      SPIN=i_iter.Next();
      P_SPIN=i_iter2.Next();
    } // while (!net_iter.End())
  }  // while markov

  // In case of a cyclic attractor, we want to interrupt
  if (cyclic)  {
//       printf("Cyclic attractor!\n");
      acceptance=0.0;
      return 0;
  } else {
    acceptance=double(changes)/double(num_of_nodes);
    return changes;
  }
}
//###################################################################################
//The same function as before, but rather than parallel update, it pics the nodes to update 
//randomly
//###################################################################################
double PottsModel::HeatBathLookupZeroTemp(double gamma, double prob, unsigned int max_sweeps)
{
  DLList_Iter<NNode*> iter;
  DLList_Iter<NLink*> l_iter;
  DLList_Iter<unsigned int*> i_iter, i_iter2;
  NNode *node, *n_cur;
  NLink *l_cur;
  unsigned int new_spin, spin_opt, old_spin, spin, sweep;
  long r;// degree;
  unsigned long changes;
  double delta=0, h, deltaE, deltaEmin,w,degree;
  //HugeArray<int> neighbours;
  bool cyclic;

  sweep=0;
  changes=0;
  while (sweep<max_sweeps)
  {
    cyclic=true;
    sweep++;
    //ueber alle Knoten im Netz
    for (unsigned long n=0; n<num_of_nodes; n++)
    {
      r=-1;
      while ((r<0) || (r>(long)num_of_nodes-1))
	r=RNG_INTEGER(0,num_of_nodes-1);
      /* r=long(double(num_of_nodes*double(rand())/double(RAND_MAX+1.0)));*/
      node=net->node_list->Get(r);
      // Wir zaehlen, wieviele Nachbarn von jedem spin vorhanden sind
      // erst mal alles Null setzen
      for (unsigned int i=0; i<=q; i++) neighbours[i]=0;
      degree=node->Get_Weight();
      //Loop over all links (=neighbours)
      l_cur=l_iter.First(node->Get_Links());
      while (!l_iter.End())
      {
        //printf("%s %s\n",node->Get_Name(),n_cur->Get_Name());
	w=l_cur->Get_Weight();
	if (node==l_cur->Get_Start()) {
	  n_cur=l_cur->Get_End();
	} else { 
	  n_cur=l_cur->Get_Start(); 
	}
        neighbours[n_cur->Get_ClusterIndex()]+=w;
        l_cur=l_iter.Next();
      }
      //Search optimal Spin      
      old_spin=node->Get_ClusterIndex();
      //degree=node->Get_Degree();
      switch (operation_mode) {
      case 0: { 
	delta=1.0; 
	break;
      } 
      case 1: { //newman modularity 
	prob=degree/total_degree_sum; 
	delta=degree; 
	break;
       }
      }


      spin_opt=old_spin;
      deltaEmin=0.0;
      for (spin=1; spin<=q; spin++)  // alle moeglichen Spins
      {
          if (spin!=old_spin)
          {
            h=color_field[spin]+delta-color_field[old_spin];
            deltaE=double(neighbours[old_spin]-neighbours[spin])+gamma*prob*double(h);
            if (deltaE<deltaEmin) {
              spin_opt=spin;
              deltaEmin=deltaE;
            }
          }
      } // for spin

      //-------------------------------
      //Now update the spins
      new_spin=spin_opt;
      if (new_spin!=old_spin) // Did we really change something??
      {
        changes++;
        node->Set_ClusterIndex(new_spin);
        color_field[old_spin]-=delta;
        color_field[new_spin]+=delta;

        //Qmatrix update
        //iteration over all neighbours
        l_cur=l_iter.First(node->Get_Links());
        while (!l_iter.End())
        {
	   w=l_cur->Get_Weight();
	   if (node==l_cur->Get_Start()) {
	     n_cur=l_cur->Get_End();
	   } else { 
	     n_cur=l_cur->Get_Start(); 
	   }
          Qmatrix[old_spin][n_cur->Get_ClusterIndex()]-=w;
          Qmatrix[new_spin][n_cur->Get_ClusterIndex()]+=w;
          Qmatrix[n_cur->Get_ClusterIndex()][old_spin]-=w;
          Qmatrix[n_cur->Get_ClusterIndex()][new_spin]+=w;
          Qa[old_spin]-=w;
          Qa[new_spin]+=w;
          l_cur=l_iter.Next();
        }  // while l_iter
       }
    } // for n
  }  // while markov

  acceptance=double(changes)/double(num_of_nodes)/double(sweep);
  return acceptance;
}
//#####################################################################################
//This function performs a parallel update at Terperature T
//#####################################################################################
long PottsModel::HeatBathParallelLookup(double gamma, double prob, double kT, unsigned int max_sweeps)
{
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
  double h, delta=0, norm, r, beta,minweight, mag, prefac=0,w, degree;
  bool cyclic=0, found;
  unsigned long num_of_nodes;

  sweep=0;
  changes=1;
  num_of_nodes=net->node_list->Size();
  while (sweep<max_sweeps && changes)
  {
    cyclic=true;
    sweep++;
    changes=0;
    //Loop over all nodes
    node=net_iter.First(net->node_list);
    SPIN=i_iter.First(new_spins);
    while (!net_iter.End())
    {
      // Initialize neighbours and weights
      problemcount=0;
      for (unsigned int i=0; i<=q; i++) {
        neighbours[i]=0;
        weights[i]=0;
      }
      norm=0.0;
      degree=node->Get_Weight();
      //Loop over all links (=neighbours)
      l_cur=l_iter.First(node->Get_Links());
      while (!l_iter.End())
      {
        //printf("%s %s\n",node->Get_Name(),n_cur->Get_Name());
	w=l_cur->Get_Weight();
	if (node==l_cur->Get_Start()) {
	  n_cur=l_cur->Get_End();
	} else { 
	  n_cur=l_cur->Get_Start(); 
	}
        neighbours[n_cur->Get_ClusterIndex()]+=w;
        l_cur=l_iter.Next();
      }
      //Search optimal Spin      
      old_spin=node->Get_ClusterIndex();
      //degree=node->Get_Degree();
      switch (operation_mode) {
      case 0: { 
	prefac=1.0; 
	delta=1.0; 
	break;
      } 
      case 1: { //newman modularity 
	prefac=1.0; 
	prob=degree/total_degree_sum; 
	delta=degree; 
	break;
       }
      }
      spin_opt=old_spin;
      beta=1.0/kT*prefac;
      minweight=0.0;
      weights[old_spin]=0.0;
      for (unsigned spin=1; spin<=q; spin++)  // loop over all possible new spins
      {
          if (spin!=old_spin) // only if we have a different than old spin!
          {
            h=color_field[spin]+delta-color_field[old_spin];
            weights[spin]=double(neighbours[old_spin]-neighbours[spin])+gamma*prob*double(h);
            if (weights[spin]<minweight) minweight=weights[spin];
          }
      }   // for spin      
      for (unsigned spin=1; spin<=q; spin++)  // loop over all possibe spins
      {
            weights[spin]-=minweight;         // subtract minweight
                                              // to avoid numerical problems with large exponents
            weights[spin]=exp(-beta*weights[spin]);
            norm+=weights[spin];
      }   // for spin

     //now choose a new spin
     r = RNG_UNIF(0, norm);
     /* norm*double(rand())/double(RAND_MAX + 1.0); */
     new_spin=1;
     found=false;
     while (!found && new_spin<=q) {
       if (r<=weights[new_spin]) {
         spin_opt=new_spin;
         found=true;
         break;
       } else r-=weights[new_spin];
       new_spin++;
     }
     if (!found) {
//         printf(".");
        problemcount++;
     }
     //Put new spin on list
     *SPIN=spin_opt;

     node=net_iter.Next();
     SPIN=i_iter.Next();
    } // while !net_iter.End()

    //-------------------------------
    //now update all spins
    node=net_iter.First(net->node_list);
    SPIN=i_iter.First(new_spins);
    P_SPIN=i_iter2.First(previous_spins);
    while (!net_iter.End())
    {
      old_spin=node->Get_ClusterIndex();
      new_spin=*SPIN;
      if (new_spin!=old_spin) // Did we really change something??
      {
        changes++;
        node->Set_ClusterIndex(new_spin);
        if (new_spin!=*P_SPIN) cyclic=false;
        *P_SPIN=old_spin;
        color_field[old_spin]-=delta;
        color_field[new_spin]+=delta;

        //Qmatrix update
        //iteration over all neighbours
        l_cur=l_iter.First(node->Get_Links());
        while (!l_iter.End())
        {
	   w=l_cur->Get_Weight();
	   if (node==l_cur->Get_Start()) {
	     n_cur=l_cur->Get_End();
	   } else { 
	     n_cur=l_cur->Get_Start(); 
	   }
          Qmatrix[old_spin][n_cur->Get_ClusterIndex()]-=w;
          Qmatrix[new_spin][n_cur->Get_ClusterIndex()]+=w;
          Qmatrix[n_cur->Get_ClusterIndex()][old_spin]-=w;
          Qmatrix[n_cur->Get_ClusterIndex()][new_spin]+=w;
          Qa[old_spin]-=w;
          Qa[new_spin]+=w;
          l_cur=l_iter.Next();
        }  // while l_iter
      }
      node=net_iter.Next();
      SPIN=i_iter.Next();
      P_SPIN=i_iter2.Next();
    } // while (!net_iter.End())

  }  // while markov
  max_q=0;
  for (unsigned int i=1; i<=q; i++) if (color_field[i]>max_q) max_q=long(color_field[i]);
  mag=(double(max_q*q)/double(num_of_nodes)-1.0)/double(q-1);

  //again, we would not like to end up in cyclic attractors
  if (cyclic && changes)  {
//       printf("Cyclic attractor!\n");
      acceptance=double(changes)/double(num_of_nodes);
      return 0;
  } else {
    acceptance=double(changes)/double(num_of_nodes);
    return changes;
  }
}
//##############################################################
// This is the function generally used for optimisation, 
// as the parallel update has its flaws, due to the cyclic attractors
//##############################################################
double PottsModel::HeatBathLookup(double gamma, double prob, double kT, unsigned int max_sweeps)
{
  DLList_Iter<NNode*> iter;
  DLList_Iter<NLink*> l_iter;
  DLList_Iter<unsigned int*> i_iter, i_iter2;
  NNode *node, *n_cur;
  NLink *l_cur;
  unsigned int new_spin, spin_opt, old_spin;
  unsigned int sweep;
  long max_q, rn;
  unsigned long changes, /*degree,*/ problemcount;
  double degree,w, delta=0, h;
  //HugeArray<int> neighbours;
  double norm, r, beta,minweight, mag, prefac=0;
  bool cyclic, found;
  long int num_of_nodes;
  sweep=0;
  changes=0;
  num_of_nodes=net->node_list->Size();
  while (sweep<max_sweeps)
  {
    cyclic=true;
    sweep++;
    //loop over all nodes in network
    for (int n=0; n<num_of_nodes; n++)
    {
      rn=-1;
      while ((rn<0) || (rn>num_of_nodes-1))
	rn=RNG_INTEGER(0, num_of_nodes-1);
      /* rn=long(double(num_of_nodes*double(rand())/double(RAND_MAX+1.0))); */
        
      node=net->node_list->Get(rn);
      // initialize the neighbours and the weights
      problemcount=0;
      for (unsigned int i=0; i<=q; i++) {
        neighbours[i]=0.0;
        weights[i]=0.0;
      }
      norm=0.0;
      degree=node->Get_Weight();
      //Loop over all links (=neighbours)
      l_cur=l_iter.First(node->Get_Links());
      while (!l_iter.End())
      {
        //printf("%s %s\n",node->Get_Name(),n_cur->Get_Name());
	w=l_cur->Get_Weight();
	if (node==l_cur->Get_Start()) {
	  n_cur=l_cur->Get_End();
	} else { 
	  n_cur=l_cur->Get_Start(); 
	}
        neighbours[n_cur->Get_ClusterIndex()]+=w;
        l_cur=l_iter.Next();
      }
      
      //Look for optimal spin
      
      old_spin=node->Get_ClusterIndex();
      //degree=node->Get_Degree();
      switch (operation_mode) {
      case 0: { 
	prefac=1.0; 
	delta=1.0; 
	break;
      } 
      case 1:  {//newman modularity 
	prefac=1.0; 
	prob=degree/total_degree_sum; 
	delta=degree; 
	break;
       }
      }
      spin_opt=old_spin;
      beta=1.0/kT*prefac;
      minweight=0.0;
      weights[old_spin]=0.0;
      for (unsigned spin=1; spin<=q; spin++)  // all possible new spins
      {
          if (spin!=old_spin) // except the old one!
          {
            h=color_field[spin]-(color_field[old_spin]-delta);
            weights[spin]=neighbours[old_spin]-neighbours[spin]+gamma*prob*h;
            if (weights[spin]<minweight) minweight=weights[spin];
          }
      }   // for spin      
      for (unsigned spin=1; spin<=q; spin++)  // all possible new spins
      {
            weights[spin]-=minweight;         // subtract minweigt
                                              // for numerical stability
            weights[spin]=exp(-beta*weights[spin]);
            norm+=weights[spin];
      }   // for spin
      

     //choose a new spin
/*      r = norm*double(rand())/double(RAND_MAX + 1.0); */
     r=RNG_UNIF(0, norm);
     new_spin=1;
     found=false;
     while (!found && new_spin<=q) {
       if (r<=weights[new_spin]) {
         spin_opt=new_spin;
         found=true;
         break;
       } else r-=weights[new_spin];
       new_spin++;
     }
     if (!found) {
//         printf(".");
        problemcount++;
     }
    //-------------------------------
    //now set the new spin
    new_spin=spin_opt;
    if (new_spin!=old_spin) // Did we really change something??
    {
        changes++;
        node->Set_ClusterIndex(new_spin);
        color_field[old_spin]-=delta;
        color_field[new_spin]+=delta;

        //Qmatrix update
        //iteration over all neighbours
        l_cur=l_iter.First(node->Get_Links());
        while (!l_iter.End())
        {
	   w=l_cur->Get_Weight();
	   if (node==l_cur->Get_Start()) {
	     n_cur=l_cur->Get_End();
	   } else { 
	     n_cur=l_cur->Get_Start(); 
	   }
          Qmatrix[old_spin][n_cur->Get_ClusterIndex()]-=w;
          Qmatrix[new_spin][n_cur->Get_ClusterIndex()]+=w;
          Qmatrix[n_cur->Get_ClusterIndex()][old_spin]-=w;
          Qmatrix[n_cur->Get_ClusterIndex()][new_spin]+=w;
          Qa[old_spin]-=w;
          Qa[new_spin]+=w;
          l_cur=l_iter.Next();
        }  // while l_iter
      }
    } // for n
  }  // while markov
  max_q=0;

  for (unsigned int i=1; i<=q; i++) if (color_field[i]>max_q) max_q=long(color_field[i]+0.5);
  mag=(double(max_q*q)/double(num_of_nodes)-1.0)/double(q-1);

    acceptance=double(changes)/double(num_of_nodes)/double(sweep);
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
  NNode *start_node=0, *n_cur, *neighbor, *max_aff_node, *node;
  NLink *l_cur;
  bool found=false, add=false, remove=false;
  double degree, delta_aff_add, delta_aff_rem, max_delta_aff, Ks=0.0, Kr=0, kis, kir, w;
  long community_marker=5;
  long to_do_marker=10;
  double inner_links=0, outer_links=0, aff_r, aff_s, css;

  to_do=new DLList<NNode*>;
  community=new DLList<NNode*>;
  
  // find the node in the network
  n_cur=iter.First(net->node_list);
  while (!found && !iter.End())
  {
    if (0==strcmp(n_cur->Get_Name(),nodename)) {
      start_node=n_cur;
      found=true;
      start_node->Set_Affinity(0.0);
      community->Push(start_node);
      start_node->Set_Marker(community_marker);
      Ks=start_node->Get_Weight();
      Kr=total_degree_sum-start_node->Get_Weight();
    }
    n_cur=iter.Next();
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
  neighbor=iter.First(start_node->Get_Neighbours());
  while (!iter.End()) {
//     printf("Adding node %s to comunity.\n",neighbor->Get_Name());
    community->Push(neighbor);
    neighbor->Set_Marker(community_marker);
    Ks+=neighbor->Get_Weight();
    Kr-=neighbor->Get_Weight();
    neighbor=iter.Next();
  }
  node=iter.First(community);
  while (!iter.End()) {
    //now add at the second neighbors to the to_do list
    neighbor=iter2.First(node->Get_Neighbours());
    while (!iter2.End()) {
      if ((long)neighbor->Get_Marker()!=community_marker && (long)neighbor->Get_Marker()!=to_do_marker) {
	to_do->Push(neighbor);
	neighbor->Set_Marker(to_do_marker);
// 	printf("Adding node %s to to_do list.\n",neighbor->Get_Name());
      }
      neighbor=iter2.Next();
    }
    node=iter.Next();
  }
  
  //#############
  //repeat, as long as we are still adding nodes to the communtiy
  //#############
  add=true;
  remove=true;
  while (add || remove) {
      //#############################
      //calculate the affinity changes of all nodes for adding every node in the to_do list to the community
      //##############################

      IGRAPH_ALLOW_INTERRUPTION(); /* This is not clean.... */

      max_delta_aff=0.0;
      max_aff_node=NULL;
      add=false;
      node=iter.First(to_do);
      while (!iter.End()) {
	//printf("Checking Links of %s\n",node->Get_Name());
	degree=node->Get_Weight();
	kis=0.0;
	kir=0.0;
	// For every of the neighbors, check, count the links to the community
	l_cur=l_iter.First(node->Get_Links());
	while (!l_iter.End())
	{
	  w=l_cur->Get_Weight();
	  if (node==l_cur->Get_Start()) {
	    n_cur=l_cur->Get_End();
	  } else { 
	    n_cur=l_cur->Get_Start(); 
	  }
	  if ((long)n_cur->Get_Marker()==community_marker) {
	    kis+=w;  //the weight/number of links to the community
	  } else {
	    kir+=w;  //the weight/number of links to the rest of the network
	  }
	  l_cur=l_iter.Next();
	}
	aff_r=kir-gamma/total_degree_sum*(Kr-degree)*degree;
	aff_s=kis-gamma/total_degree_sum*Ks*degree;
	css=-gamma/total_degree_sum*degree*degree*0.5;
	//printf("Node: %s\taff_r:%f\taff_s:%f\tcss:%f\tsum:%f\n",node->Get_Name(),aff_r,aff_s,css,aff_r+aff_s+2.0*css);
	delta_aff_add=aff_r-aff_s;
	// 	if (aff_s>=aff_r && delta_aff_add<=max_delta_aff) {
	if (delta_aff_add<=max_delta_aff) {
	    node->Set_Affinity(aff_s);
	    max_delta_aff=delta_aff_add;
	    max_aff_node=node;
	    add=true;
	}
	//printf("%s in to_do list with affinity %f\n",node->Get_Name(),node->Get_Affinity());
	node=iter.Next();
      }
      //################
      //calculate the affinity changes for removing every single node from the community
      //################
      inner_links=0;
      outer_links=0;
      remove=false;
      node=iter.First(community);
      while (!iter.End()) {
	//printf("Checking Links of %s\n",node->Get_Name());
	degree=node->Get_Weight();
	kis=0.0;
	kir=0.0;
	// For every of the neighbors, check, count the links to the community
	l_cur=l_iter.First(node->Get_Links());
	while (!l_iter.End())
	{
	  w=l_cur->Get_Weight();
	  if (node==l_cur->Get_Start()) {
	    n_cur=l_cur->Get_End();
	  } else { 
	    n_cur=l_cur->Get_Start(); 
	  }
	  if ((long)n_cur->Get_Marker()==community_marker) {
	    kis+=w;
	    inner_links+=w;  //summing all w gives twice the number of inner links(weights)
	  } else {
	    kir+=w;
	    outer_links+=w;
	  }
	  l_cur=l_iter.Next();
	}
// 	if (kir+kis!=degree) {  printf("error kir=%f\tkis=%f\tk=%f\n",kir,kis,degree); }
	aff_r=kir-gamma/total_degree_sum*Kr*degree;
	aff_s=kis-gamma/total_degree_sum*(Ks-degree)*degree;
	//css=-gamma/total_degree_sum*degree*degree*0.5;
	//printf("Node: %s\taff_r:%f\taff_s:%f\tcss:%f\tsum:%f\n",node->Get_Name(),aff_r,aff_s,css,aff_r+aff_s+2.0*css);
	delta_aff_rem=aff_s-aff_r;
	node->Set_Affinity(aff_s);
	// we should not remove the nodes, we have just added
	if (delta_aff_rem<max_delta_aff) {
	  max_delta_aff=delta_aff_rem ;
	  max_aff_node=node;
	  remove=true;
	  add=false;
	}
	//printf("%s in to_do list with affinity %f\n",node->Get_Name(),node->Get_Affinity());
	node=iter.Next();
      }
      inner_links=inner_links*0.5;
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
	Ks+=max_aff_node->Get_Weight();
	Kr-=max_aff_node->Get_Weight();
// 	printf("Adding node %s to community with affinity of %f delta_aff: %f.\n",max_aff_node->Get_Name(), max_aff_node->Get_Affinity(),max_delta_aff);
	//now add all neighbors of this node, that are not already 
	//in the to_do list or in the community
	neighbor=iter.First(max_aff_node->Get_Neighbours());
	while (!iter.End()) {
	  if ((long)neighbor->Get_Marker()!=community_marker && (long)neighbor->Get_Marker()!=to_do_marker) {
	    to_do->Push(neighbor);
	    neighbor->Set_Marker(to_do_marker);
	    //printf("Adding node %s to to_do list.\n",neighbor->Get_Name());
	  }
	  neighbor=iter.Next();
	}
      }
      if (remove) {
	//################
	//remove those with negative affinities
	//################
	community->fDelete(max_aff_node);
	max_aff_node->Set_Marker(to_do_marker);
	//update the sum of degrees in the community
	Ks-=max_aff_node->Get_Weight();
	Kr+=max_aff_node->Get_Weight();
	//add the node to to_do again
	to_do->Push(max_aff_node);
// 	printf("Removing node %s from community with affinity of %f delta_aff: %f.\n",max_aff_node->Get_Name(), max_aff_node->Get_Affinity(),max_delta_aff);
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
    *cohesion=inner_links-gamma/total_degree_sum*Ks*Ks*0.5;
  } 
  if (adhesion) {
    *adhesion=outer_links-gamma/total_degree_sum*Ks*Kr;
  }
  if (my_inner_links) {
    *my_inner_links=inner_links;
  }
  if (my_outer_links) {
    *my_outer_links=outer_links;
  }
  if (result) {
    node=iter.First(community);
    igraph_vector_resize(result, 0);
    while (!iter.End()) {
      // printf("%s in community.\n",node->Get_Name());
      // fprintf(file,"%s\t%f\n",node->Get_Name(),node->Get_Affinity());
      IGRAPH_CHECK(igraph_vector_push_back(result, node->Get_Index()));
      node=iter.Next();
    }
  }
//   printf("%d nodes in community around %s\n",community->Size(),start_node->Get_Name());
//   fclose(file);
  unsigned int size=community->Size();
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
			       double kT, double gamma)
{
  NNode *n_cur, *n_cur2;
  bool found;
  double p_in, p_out, log_num_exp, a1,a2,a3,p,p1,p2;
  long n,N,lin,lout;
  DLList_Iter<NNode*> iter, iter2;
  HugeArray<int> inner_links;
  HugeArray<int> outer_links;
  HugeArray<int> nodes;
  
  //den Header schreiben
  p=2.0*double(num_of_links)/double(num_of_nodes)/double(num_of_nodes-1);
//   fprintf(file,"      Nodes=\t%lu\n",num_of_nodes);
//   fprintf(file,"      Links=\t%lu\n",num_of_links);
//   fprintf(file,"          q=\t%d\n",q);
//   fprintf(file,"          p=\t%f\n",p);
//   fprintf(file," Modularity=\t%f\n",calculate_Q());
//   fprintf(file,"Temperature=\t%f\n", kT);
//   fprintf(file,"Cluster\tNodes\tInnerLinks\tOuterLinks\tp_in\tp_out\t<Ln(#comm.)>\n");
  
  if (modularity)  { *modularity=calculate_genQ(gamma); }
  if (temperature) { *temperature=kT; }

  if (csize || membership) {
    // TODO: count the number of clusters
    for (unsigned int spin=1; spin<=q; spin++)
      {
	inner_links[spin]=0;
	outer_links[spin]=0;
	nodes[spin]=0;
	n_cur=iter.First(net->node_list);
	while (!iter.End())
	  {
	    if (n_cur->Get_ClusterIndex()==spin)
	      {
		nodes[spin]++;
		n_cur2=iter2.First(n_cur->Get_Neighbours());
		while (!iter2.End())
		  {
		    if (n_cur2->Get_ClusterIndex()==spin) inner_links[spin]++;
		    else outer_links[spin]++;
		    n_cur2=iter2.Next();
		  }
	      }
	    n_cur=iter.Next();
	  }
      }
  }
  if (csize) {
    igraph_vector_resize(csize, 0);
    for (unsigned int spin=1; spin<=q; spin++)
      {
	if (nodes[spin]>0)
	  {
	    inner_links[spin]/=2;
	    //    fprintf(file,"Cluster\tNodes\tInnerLinks\tOuterLinks\tp_in\tp_out\n");
	    if (nodes[spin]>1)
	      p_in=2.0*double(inner_links[spin])/double(nodes[spin])/double(nodes[spin]-1);
	    else p_in=0.0;
	    p_out=double(outer_links[spin])/double(nodes[spin])/double(num_of_nodes-nodes[spin]);
	    N=num_of_nodes;
	    n=nodes[spin];
	    lin=inner_links[spin];
	    lout=outer_links[spin];
	    a1=N*log((double)N)-n*log((double)n)*(N-n)*log((double)N-n);
	    if ((lin==long(n*(n-1)*0.5+0.5)) || (n==1)) a2=0.0;
	    else a2=(n*(n-1)*0.5    )*log((double)n*(n-1)*0.5    )-(n*(n-1)*0.5    )-
		   (n*(n-1)*0.5-lin)*log((double)n*(n-1)*0.5-lin)+(n*(n-1)*0.5-lin)-
		   lin*log((double)lin            )+lin;
	    
	    if ((lout==n*(N-n)) || n==N) a3=0.0;
	    else a3=(n*(N-n)     )*log((double)n*(N-n)     )-(n*(N-n))-
		   (n*(N-n)-lout)*log((double)n*(N-n)-lout)+(n*(N-n)-lout)-
		   lout*log((double)lout        )+lout;
	    p1=(lin+lout)*log((double)p);
	    p2=(0.5*n*(n-1)-lin + n*(N-n)-lout)*log((double)1.0-p);
	    log_num_exp=a1+a2+a3+p1+p2;
	    //       fprintf(file,"%d\t%d\t%d\t%d\t%f\t%f\t%f\n",spin,nodes[spin], inner_links[spin], outer_links[spin], p_in, p_out,log_num_exp);
	    IGRAPH_CHECK(igraph_vector_push_back(csize, nodes[spin]));
	  }
      }
    //   fprintf(file,"\n");
  }
  
  //die Elemente der Cluster
  if (membership) {
    long int no=-1;
    IGRAPH_CHECK(igraph_vector_resize(membership, num_of_nodes));
    for (unsigned int spin=1; spin<=q; spin++)
      {
	if (nodes[spin]>0) {
	  no++;
	}
	found=false;
	n_cur=iter.First(net->node_list);
	while (!iter.End())
	  {
	    if (n_cur->Get_ClusterIndex()==spin)
	      {
		//         fprintf(file,"%d\t%s\n",spin,n_cur->Get_Name());
		VECTOR(*membership)[ n_cur->Get_Index() ]=no;
		found=true;
	      }
	    n_cur=iter.Next();
	  }
	//     if (found) fprintf(file,"\n");
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
long PottsModel::WriteSoftClusters(char *filename, double threshold)
{
  FILE *file;
  NNode *n_cur, *n_cur2;
  DLList_Iter<NNode*> iter, iter2;
  DL_Indexed_List<ClusterList<NNode*>*> *cl_list, *old_clusterlist;
  ClusterList<NNode*> *cl_cur;
  
  double max;
  
  file=fopen(filename,"w");
  if (!file) {
    printf("Could not open %s for writing.\n",filename);
    return -1;
  }

  max=correlation[0]->Get(0);
  //printf("max=%f\n",max);
  cl_list=new DL_Indexed_List<ClusterList<NNode*>*>();
  
  n_cur=iter.First(net->node_list);
  while (!iter.End())
  {
    cl_cur=new ClusterList<NNode*>();
    cl_list->Push(cl_cur);
    n_cur2=iter2.First(net->node_list);
    while (!iter2.End())
    {
      if (double(correlation[n_cur->Get_Index()]->Get(n_cur2->Get_Index()))/max>threshold)
        cl_cur->Push(n_cur2);
      n_cur2=iter2.Next();
    }
    n_cur=iter.Next();
  }
  old_clusterlist=net->cluster_list;
  net->cluster_list=cl_list;
  clear_all_markers(net);
  //printf("Es gibt %d Cluster\n",cl_list->Size());
  reduce_cliques2(net, false, 15);
  //printf("Davon bleiben %d Cluster uebrig\n",cl_list->Size());
  clear_all_markers(net);
  while (net->cluster_list->Size()){
    cl_cur=net->cluster_list->Pop();
    while (cl_cur->Size())
    {
      n_cur=cl_cur->Pop();
      fprintf(file,"%s\n",n_cur->Get_Name());
      //printf("%s\n",n_cur->Get_Name());
    }
    fprintf(file,"\n");
  }
  net->cluster_list=old_clusterlist;
  fclose(file);

  return 1;
}
//#############################################################################
// Performs a gamma sweep
//#############################################################################
double PottsModel::GammaSweep(double gamma_start, double gamma_stop, double prob, unsigned int steps, bool non_parallel, int repetitions)
{
  double stepsize;
  double kT, kT_start;
  long changes;
  double gamma, acc;
  NNode *n_cur, *n_cur2;
  DLList_Iter<NNode*> iter, iter2;
  
  stepsize=(gamma_stop-gamma_start)/double(steps);
  
  n_cur=iter.First(net->node_list);
  while (!iter.End())
  {
    correlation[n_cur->Get_Index()]=new HugeArray<double>();
    n_cur2=iter2.First(net->node_list);
    while (!iter2.End())
    {
      correlation[n_cur->Get_Index()]->Set(n_cur->Get_Index())=0.0;
      n_cur2=iter2.Next();
    }
    n_cur=iter.Next();
  }

  for (unsigned int n=0; n<=steps; n++)
  {
    assign_initial_conf(-1);
    initialize_Qmatrix();
    gamma=gamma_start+stepsize*n;
    kT=0.5;
    acceptance=0.5;
    while (acceptance<(1.0-1.0/double(q))*0.95)     //wollen 95% Acceptance
    {
          kT*=1.1;
          //initialize_lookup(kT,kmax,net->node_list->Size());
          if (!non_parallel) HeatBathParallelLookup(gamma,prob, kT,25);
          else HeatBathLookup(gamma,prob, kT,25);
          printf("kT=%f acceptance=%f\n", kT, acceptance);
    }
    printf("Starting with gamma=%f\n", gamma);
    kT_start=kT;
    
    for (int i=0; i<repetitions; i++)
    {
      changes=1;
      kT=kT_start;
      assign_initial_conf(-1);
      initialize_Qmatrix();
      while ((changes>0) && (kT>0.01)) 
      {
          kT=kT*0.99;
          //initialize_lookup(kT,kmax,net->node_list->Size());
          if (!non_parallel) {
	    changes=HeatBathParallelLookup(gamma, prob, kT, 50);
              printf("kT: %f   \t Changes %li\n",kT, changes);
          } else {
	    acc=HeatBathLookup(gamma, prob, kT, 50);
             if (acc>(1.0-1.0/double(q))*0.01) changes=1; else changes=0;
             printf("kT: %f   Acceptance: %f\n",kT, acc);
          }
      }
      printf("Finisched with acceptance: %1.6f bei kT=%2.4f und gamma=%2.4f\n",acceptance,kT, gamma);
//      fprintf(file,"%f\t%f\n",gamma_,acceptance);
//      fprintf(file2,"%f\t%f\n",gamma_,kT);
   //   fprintf(file3,"%f\t%d\n",gamma_,count_clusters(5));

      //Die Correlation berechnen
      n_cur=iter.First(net->node_list);
      while (!iter.End())
      {
        n_cur2=iter2.First(net->node_list);
        while (!iter2.End())
        {
          if (n_cur->Get_ClusterIndex()==n_cur2->Get_ClusterIndex())
          {
            correlation[n_cur->Get_Index()]->Set(n_cur2->Get_Index())+=0.5;
          }
          n_cur2=iter2.Next();
        }
        n_cur=iter.Next();
      }
    } // for i
} //for n
  return kT;
}
//#############################################################################
//Performs a Gamma sweep at zero T
//#############################################################################
double PottsModel::GammaSweepZeroTemp(double gamma_start, double gamma_stop, double prob, unsigned int steps, bool non_parallel, int repetitions)
{
  double stepsize;
  long changes;
  double gamma, acc;
  long runs;
  NNode *n_cur, *n_cur2;
  DLList_Iter<NNode*> iter, iter2;

  stepsize=(gamma_stop-gamma_start)/double(steps);

  n_cur=iter.First(net->node_list);
  while (!iter.End())
  {
    correlation[n_cur->Get_Index()]=new HugeArray<double>();
    n_cur2=iter2.First(net->node_list);
    while (!iter2.End())
    {
      correlation[n_cur->Get_Index()]->Set(n_cur->Get_Index())=0.0;
      n_cur2=iter2.Next();
    }
    n_cur=iter.Next();
  }

  for (unsigned int n=0; n<=steps; n++)
  {
    assign_initial_conf(-1);
    initialize_Qmatrix();
    gamma=gamma_start+stepsize*n;
    printf("Starting with gamma=%f\n", gamma);
    for (int i=0; i<repetitions; i++)
    {
      changes=1;
      assign_initial_conf(-1);
      initialize_Qmatrix();
      runs=0;
      while (changes>0 && runs<250)
      {
          //initialize_lookup(kT,kmax,net->node_list->Size());
          if (!non_parallel) {
	    changes=HeatBathParallelLookupZeroTemp(gamma, prob, 1);
              printf("Changes %li\n", changes);
          } else {
            acc=HeatBathLookupZeroTemp(gamma, prob, 1);
            if (acc>(1.0-1.0/double(q))*0.01) changes=1; else changes=0;
            printf("Acceptance: %f\n", acc);
          }
          runs++;
      }
      printf("Finisched with Modularity: %1.6f bei Gamma=%1.6f\n",calculate_Q(), gamma);
//      fprintf(file,"%f\t%f\n",gamma_,acceptance);
//      fprintf(file2,"%f\t%f\n",gamma_,kT);
   //   fprintf(file3,"%f\t%d\n",gamma_,count_clusters(5));

      //Die Correlation berechnen
      n_cur=iter.First(net->node_list);
      while (!iter.End())
      {
        n_cur2=iter2.First(net->node_list);
        while (!iter2.End())
        {
          if (n_cur->Get_ClusterIndex()==n_cur2->Get_ClusterIndex())
          {
            correlation[n_cur->Get_Index()]->Set(n_cur2->Get_Index())+=0.5;
            correlation[n_cur2->Get_Index()]->Set(n_cur->Get_Index())+=0.5;
          }
          n_cur2=iter2.Next();
        }
        n_cur=iter.Next();
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
long PottsModel::WriteCorrelationMatrix(char *filename)
{
  FILE *file, *file2;
  char filename2[255];
  NNode *n_cur, *n_cur2;
  DLList_Iter<NNode*> iter, iter2;

  sprintf(filename2,"%s.mat",filename);
  file=fopen(filename,"w");
  if (!file) {
    printf("Could not open %s for writing.\n",filename);
    return -1;
  }
  file2=fopen(filename2,"w");
  if (!file2) {
    printf("Could not open %s for writing.\n",filename2);
    return -1;
  }
  //write the header in one line
  n_cur=iter.First(net->node_list);
  while (!iter.End())
  {  
      fprintf(file, "\t%s",n_cur->Get_Name());
      n_cur=iter.Next();
  }    
  fprintf(file, "\n");

  //fprintf(file, "%d\t%d\n",net->node_list->Size(),net->node_list->Size());

  long r=0,c=0;
  n_cur=iter.First(net->node_list);
  while (!iter.End())
  {
    fprintf(file, "%s",n_cur->Get_Name());
    r++;
    n_cur2=iter2.First(net->node_list);
    while (!iter2.End())
    {
      c++;
      fprintf(file,"\t%f",correlation[n_cur->Get_Index()]->Get(n_cur2->Get_Index()));
      fprintf(file2,"%li\t%li\t%f\n",r,c,correlation[n_cur->Get_Index()]->Get(n_cur2->Get_Index()));
      n_cur2=iter2.Next();
    }
    fprintf(file,"\n");
    n_cur=iter.Next();
  }
  fclose(file);
  fclose(file2);
  return 1;
}
//##############################################################################
