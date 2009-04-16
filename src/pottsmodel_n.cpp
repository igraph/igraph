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
#undef DEBUG

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pottsmodel_2.h"
#include "NetRoutines.h"

#include "random.h"


//#################################################################################################
PottsModelN::PottsModelN(network *n, unsigned int num_communities, bool directed)
{
	//Set internal variable
	net	= n;
	q	= num_communities;
	
	is_directed = directed;
	
	is_init = false;
	
	num_nodes	= net->node_list->Size();
}
//#######################################################
//Destructor of PottsModel
//########################################################
PottsModelN::~PottsModelN()
{
	delete degree_pos_in;
	delete degree_neg_in;
	delete degree_pos_out;
	delete degree_neg_out;
	
	delete degree_community_pos_in;
	delete degree_community_neg_in;	
	delete degree_community_pos_out;
	delete degree_community_neg_out;
	
	delete weights;
	delete neighbours;
	delete csize;
	
	delete spin;
	
	return;
}

void PottsModelN::assign_initial_conf(bool init_spins)
{
	#ifdef DEBUG
	printf("Start assigning.\n");
	#endif
	int s;
	DLList_Iter<NNode*> iter;
	DLList_Iter<NLink*> l_iter;
	NNode *n_cur;
	NLink *l_cur;
	

	if(init_spins)
	{
		#ifdef DEBUG
		printf("Initializing spin.\n");
		#endif
		//Bookkeeping of the various degrees (positive/negative) and (in/out)
		degree_pos_in	= new double[num_nodes]; //Postive indegree of the nodes (or sum of weights)
		degree_neg_in	= new double[num_nodes]; //Negative indegree of the nodes (or sum of weights)
		degree_pos_out	= new double[num_nodes]; //Postive outdegree of the nodes (or sum of weights)
		degree_neg_out	= new double[num_nodes]; //Negative outdegree of the nodes (or sum of weights)	
				
		spin			= new unsigned int[num_nodes]; //The spin state of each node
	}
	
	if (is_init)
	{
		delete degree_community_pos_in;
		delete degree_community_neg_in;
		delete degree_community_pos_out;
		delete degree_community_neg_out;
		
		delete weights;
		delete neighbours;
		delete csize;
	}
	
	is_init = true;
	
	//Bookkeep of occupation numbers of spin states or the number of links in community...
	degree_community_pos_in		= new double[q+1]; //Positive sum of indegree for communities
	degree_community_neg_in		= new double[q+1]; //Negative sum of indegree for communities
	degree_community_pos_out	= new double[q+1];//Positive sum of outegree for communities
	degree_community_neg_out	= new double[q+1]; //Negative sum of outdegree for communities	

	//...and of weights and neighbours for in the HeathBathLookup
	weights						= new double[q+1]; //The weights for changing to another spin state
	neighbours					= new double[q+1]; //The number of neighbours (or weights) in different spin states
	csize						= new unsigned int[q+1]; //The number of nodes in each community	
	

	//Initialize communities
	for (unsigned int i=0; i<=q; i++) 
	{
		degree_community_pos_in[i]	= 0.0;
		degree_community_neg_in[i]	= 0.0;	
		degree_community_pos_out[i]	= 0.0;
		degree_community_neg_out[i]	= 0.0;
		
		csize[i]					= 0.0;
	}
	
	//Initialize vectors
	if (init_spins)
	{	
		for (unsigned int i = 0; i < num_nodes; i++)
		{
			degree_pos_in[i]	= 0.0;
			degree_neg_in[i]	= 0.0;
			degree_pos_out[i]	= 0.0;
			degree_neg_out[i]	= 0.0;
			
			#ifdef DEBUG
			printf("Initializing spin %d", i);
			#endif
			spin[i]	= 0;
		}
	}
	m_p=0.0;
	m_n=0.0;
	//Set community for each node, and 
	//correctly store it in the bookkeeping
	
	double sum_weight_pos_in, sum_weight_pos_out, sum_weight_neg_in, sum_weight_neg_out;
	//double av_w = 0.0, av_k=0.0;
	//int l = 0;
	#ifdef DEBUG
	printf("Visiting each node.\n");
	#endif
	for (unsigned int v = 0; v < num_nodes; v++)
	{
		if (init_spins)
		{
			s = RNG_INTEGER(1, q);  //The new spin s
			spin[v] = (unsigned int)s;
		}
		else
			s = spin[v];
		
		#ifdef DEBUG		
		printf("Spin %d assigned to node %d.\n", s, v);
		#endif
			
		n_cur				=  net->node_list->Get(v);
		
		l_cur				= l_iter.First(n_cur->Get_Links());
		
		sum_weight_pos_in	= 0.0;
		sum_weight_pos_out	= 0.0;
		sum_weight_neg_in	= 0.0;
		sum_weight_neg_out	= 0.0;
		
		while (!l_iter.End())
		{
			double w = l_cur->Get_Weight();
			//av_w = (av_w*l + w)/(l+1); //Average weight
			//l++;
			if (l_cur->Get_Start() == n_cur) //From this to other, so outgoing link
				if (w > 0)
					sum_weight_pos_out += w;   //Increase positive outgoing weight
				else
					sum_weight_neg_out -= w;	//Increase negative outgoing weight
			else
				if (w > 0)
					sum_weight_pos_in += w;   //Increase positive incoming weight
				else
					sum_weight_neg_in -= w;	//Increase negative incoming weight			
			
			l_cur=l_iter.Next();
		}
		
		if (!is_directed)
		{
			double sum_weight_pos		= sum_weight_pos_out + sum_weight_pos_in; 
				   sum_weight_pos_out	= sum_weight_pos; 
				   sum_weight_pos_in	= sum_weight_pos;
			double sum_weight_neg = sum_weight_neg_out + sum_weight_neg_in;
				   sum_weight_neg_out	= sum_weight_neg; 
				   sum_weight_neg_in	= sum_weight_neg;			
		}
		
		//av_k = (av_k*l + sum_weight_pos_in)/(l+1); //Average k
		
		if (init_spins)
		{
			//Set the degrees correctly
			degree_pos_in[v]	= sum_weight_pos_in;
			degree_neg_in[v]	= sum_weight_neg_in;
			degree_pos_out[v]	= sum_weight_pos_out;
			degree_neg_out[v]	= sum_weight_neg_out;
		}
		
		//Correct the community bookkeeping
		degree_community_pos_in[s]	+= sum_weight_pos_in;
		degree_community_neg_in[s]	+= sum_weight_neg_in;
		degree_community_pos_out[s]	+= sum_weight_pos_out;
		degree_community_neg_out[s]	+= sum_weight_neg_out;
		
		//Community just increased
		csize[s]++;
	
		//Sum the weights (notice that sum of indegrees equals sum of outdegrees)
		m_p += sum_weight_pos_in;
		m_n += sum_weight_neg_in;
	}
	
	#ifdef DEBUG
	printf("Done assigning.\n");
	#endif

	return;
}
//##############################################################
// This is the function generally used for optimisation, 
// as the parallel update has its flaws, due to the cyclic attractors
//##############################################################
double PottsModelN::HeatBathLookup(double gamma, double lambda, double t, unsigned int max_sweeps)
{
	#ifdef DEBUG
	printf("Starting sweep at temperature %f.\n", t);
	#endif
	DLList_Iter<NNode*> iter;
	DLList_Iter<NLink*> l_iter;
	DLList_Iter<unsigned int*> i_iter, i_iter2;
	NNode *node, *n_cur;
	NLink *l_cur;
	/* The new_spin contains the spin to which we will update,
	 * the spin_opt is the optional spin we will consider and
	 * the old_spin is the spin of the node we are currently
	 * changing.
	 */  
	unsigned int new_spin, spin_opt, old_spin;
	unsigned int sweep; //current sweep
	unsigned long changes, problemcount; //Number of changes and number of problems encountered

	double exp_old_spin; //The expectation value for the old spin
	double exp_spin; //The expectation value for the other spin(s)
	int v; //The node we will be investigating
	
	//The variables required for the calculations
	double delta_pos_out, delta_pos_in, delta_neg_out, delta_neg_in;
	double k_v_pos_out, k_v_pos_in, k_v_neg_out, k_v_neg_in;
	
	//weight of edge
	double w;
	
	double beta = 1/t; //Weight for probabilities
	double r = 0.0; //random number used for assigning new spin
	
	double maxweight = 0.0;
	double sum_weights = 0.0; //sum_weights for normalizing the probabilities
	
	sweep=0;
	changes=0;
	double m_pt = m_p;
	double m_nt = m_n;
	
	if (m_pt < 0.001)
		m_pt = 1;
		
	if (m_nt < 0.001)
		m_nt = 1;
		
	while (sweep<max_sweeps)
	{
		sweep++;
		//loop over all nodes in network
		for (unsigned int n = 0; n < num_nodes; n++)
		{
			//Look for a random node
			v = RNG_INTEGER(0, num_nodes-1);
			//We will be investigating node v
			
			node=net->node_list->Get(v);
			
			/*******************************************/
			// initialize the neighbours and the weights
			problemcount=0;
			for (unsigned int i=0; i<=q; i++) {
				neighbours[i]=0.0;
				weights[i]=0.0;
			}

			//Loop over all links (=neighbours)
			l_cur=l_iter.First(node->Get_Links());
			while (!l_iter.End())
			{
				w=l_cur->Get_Weight();
				if (node==l_cur->Get_Start()) {
					n_cur=l_cur->Get_End();
				} else { 
					n_cur=l_cur->Get_Start(); 
				}
				//Add the link to the correct cluster
				neighbours[spin[n_cur->Get_Index()]]+=w;
				l_cur=l_iter.Next();
			}
			//We now have the weight of the (in and out) neighbours 
			//in each cluster available to us.
			/*******************************************/
			old_spin=spin[v];
						
			//Look for optimal spin
						
			//Set the appropriate variable
			delta_pos_out	= degree_pos_out[v];
			delta_pos_in	= degree_pos_in[v];
			delta_neg_out	= degree_neg_out[v];
			delta_neg_in	= degree_neg_in[v];
			
			k_v_pos_out		= gamma*delta_pos_out/m_pt;
			k_v_pos_in		= gamma*delta_pos_in/m_pt;
			k_v_neg_out		= lambda*delta_neg_out/m_nt;
			k_v_neg_in		= lambda*delta_neg_in/m_nt;			
			
			//The expectation value for the old spin
			if (is_directed)
				exp_old_spin = (k_v_pos_out * (degree_community_pos_in[old_spin] - delta_pos_in) - 
								k_v_neg_out * (degree_community_neg_in[old_spin] - delta_neg_in)) + 
							   (k_v_pos_in * (degree_community_pos_out[old_spin] - delta_pos_out) - 
								k_v_neg_in * (degree_community_neg_out[old_spin] - delta_neg_out));
			else
				exp_old_spin = (k_v_pos_out * (degree_community_pos_in[old_spin] - delta_pos_in) - 
								k_v_neg_out * (degree_community_neg_in[old_spin] - delta_neg_in));			

			/*******************************************/
			//Calculating probabilities for each transition to another 
			//community.			
			
			maxweight=0.0;
			weights[old_spin]=0.0;
														
			for (spin_opt=1; spin_opt<=q; spin_opt++)  // all possible new spins
			{
				if (spin_opt!=old_spin) // except the old one!
				{
					if (is_directed)
						exp_spin = (k_v_pos_out * degree_community_pos_in[spin_opt] - k_v_neg_out * degree_community_neg_in[spin_opt]) + 
								   (k_v_pos_in * degree_community_pos_out[spin_opt] - k_v_neg_in * degree_community_neg_out[spin_opt]);
					else
						exp_spin = (k_v_pos_out * degree_community_pos_in[spin_opt] - k_v_neg_out * degree_community_neg_in[spin_opt]);
						
					weights[spin_opt] = (neighbours[spin_opt] - exp_spin) - (neighbours[old_spin] - exp_old_spin);
					
					if (weights[spin_opt] > maxweight)
						maxweight = weights[spin_opt];
				}
			}   // for spin      
			
			//Calculate exp. prob. an
			sum_weights = 0.0;
			for (spin_opt=1; spin_opt<=q; spin_opt++)  // all possible new spins
			{
				weights[spin_opt] -= maxweight;  //subtract maxweight for numerical stability (otherwise overflow).
				weights[spin_opt]  = exp(beta*weights[spin_opt]);
				sum_weights   += weights[spin_opt];
			}   // for spin
			/*******************************************/
						
			
			/*******************************************/
			//Choose a new spin dependent on the calculated probabilities
			r = RNG_UNIF(0, sum_weights);
			new_spin = 1;
			
			bool found = false;
			while (!found && new_spin <= q) 
			{
				if (r <= weights[new_spin]) 
				{
					spin_opt = new_spin; //We have found are new spin
					found = true;
					break;
				} 
				else 
					r -= weights[new_spin]; //Perhaps the next spin is the one we want

				new_spin++;
			}
			
			//Some weird thing happened. We haven't found a new spin
			//while that shouldn't be the case. Numerical problems?
			if (!found) 
				problemcount++;

			new_spin=spin_opt;
			//If there wasn't a problem we should have found
			//our new spin.
			/*******************************************/
			
			
			/*******************************************/
			//The new spin is available to us, so change
			//all the appropriate counters.
			if (new_spin!=old_spin) // Did we really change something??
			{
				changes++;
				spin[v] = new_spin;
				
				//The new spin increase by one, and the old spin decreases by one
				csize[new_spin]++; csize[old_spin]--;
				
				//Change the sums of degree for the old spin...		
				degree_community_pos_in[old_spin]	-= delta_pos_in;
				degree_community_neg_in[old_spin]	-= delta_neg_in;
				degree_community_pos_out[old_spin]	-= delta_pos_out;
				degree_community_neg_out[old_spin]	-= delta_neg_out;
				
				//...and for the new spin
				degree_community_pos_in[new_spin]	+= delta_pos_in;
				degree_community_neg_in[new_spin]	+= delta_neg_in;
				degree_community_pos_out[new_spin]	+= delta_pos_out;
				degree_community_neg_out[new_spin]	+= delta_neg_out;				
			}
			
			//We have no change a node from old_spin to new_spin
			/*******************************************/
		
		} // for n
	}  // while sweep
	#ifdef DEBUG
	printf("Done %d sweeps.\n", max_sweeps);
	printf("%d changes made for %d nodes.\n", changes, num_nodes);
	printf("Last node is %d and last random number is %f with sum of weights %f with spin %d.\n", v, r, sum_weights, old_spin);
	#endif
	
    return (double(changes)/double(num_nodes)/double(sweep));
}

//We need to begin at a suitable temperature. That is, a temperature at which
//enough nodes may change their initially assigned communties
double PottsModelN::FindStartTemp(double gamma, double lambda, double ts)
{
  double kT;
  kT=ts;
  //assing random initial condition
  assign_initial_conf(true);
  // the factor 1-1/q is important, since even, at infinite temperature,
  // only 1-1/q of all spins do change their state, since a randomly chooses new
  // state is with prob. 1/q the old state.
  double acceptance = 0.0;
  while (acceptance<(1.0-1.0/double(q))*0.95)      //want 95% acceptance
  {
       kT=kT*1.1;
       acceptance=HeatBathLookup(gamma,lambda, kT,50);
  }
  kT*=1.1; // just to be sure...
  return kT;
}

long PottsModelN::WriteClusters(igraph_real_t *modularity,
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
							   double lambda)
{
	#ifdef DEBUG
	printf("Start writing clusters.\n");
	#endif
	//Reassign each community so that we retrieve a community assignment 1 through num_communities	
	unsigned int *cluster_assign = new unsigned int[q+1];
	for (unsigned int i = 0; i <= q; i++)
	{
		cluster_assign[i] = 0;
	}
	
	int num_clusters = 0;
	
	//Find out what the new communities will be
	for (unsigned int i = 0; i < num_nodes; i++)
	{
		int s = spin[i];
		if (cluster_assign[s] == 0)
		{
			num_clusters++;		
			cluster_assign[s] = num_clusters;
			#ifdef DEBUG
			printf("Setting cluster %d to %d.\n", s, num_clusters);
			#endif
		}
	}
	

	DLList_Iter<NNode*> iter;
	NNode *n_cur=iter.First(net->node_list);		
	n_cur = iter.First(net->node_list);
	//And now assign each node to its new community
	q = num_clusters;	
	for (unsigned int i = 0; i < num_nodes; i++)
	{
		#ifdef DEBUG
		printf("Setting node %d to %d.\n", i, cluster_assign[spin[i]]);
		#endif
		unsigned int s = cluster_assign[spin[i]];
		spin[i] = s;
		#ifdef DEBUG
		printf("Have set node %d to %d.\n", i, s);
		#endif
	}
	assign_initial_conf(false);
	
	delete cluster_assign;
		
	if (temperature) { *temperature=t; }
	
	if (community_size)
	{
		//Initialize the vector
		IGRAPH_CHECK(igraph_vector_resize(community_size, q));
		for (unsigned int spin_opt = 1; spin_opt <= q; spin_opt++)
		{
			//Set the community size
			VECTOR(*community_size)[spin_opt-1]=csize[spin_opt];
		}
	}
	
	//Set the membership	
	if (membership) 
	{
		IGRAPH_CHECK(igraph_vector_resize(membership, num_nodes));		
		for (unsigned int i = 0; i < num_nodes; i++)
		{
			VECTOR(*membership)[ i ]= spin[i];
		}
	}
	
	double Q = 0.0; //Modularity
	if (adhesion)
	{
		IGRAPH_CHECK(igraph_matrix_resize(adhesion, q, q));
		IGRAPH_CHECK(igraph_matrix_resize(normalised_adhesion, q, q));
		
		double **num_links_pos = 0;
		double **num_links_neg = 0;
		//memory allocated for elements of rows.
		num_links_pos = new double *[q+1] ;
		num_links_neg = new double *[q+1] ;

		//memory allocated for  elements of each column.
		for( int i = 0 ; i < q+1 ; i++)
		{
			num_links_pos[i] = new double[q+1];
			num_links_neg[i] = new double[q+1];
		}
		

		
		//Init num_links
		for (unsigned int i = 0; i <= q; i++)
		{
			for (unsigned int j = 0; j <= q; j++)
			{
				num_links_pos[i][j] = 0.0;
				num_links_neg[i][j] = 0.0;
			}
		}
		
		DLList_Iter<NLink*> iter_l;
		NLink *l_cur = iter_l.First(net->link_list);
		
		double w = 0.0;
		
		while (!iter_l.End())
		{
			w = l_cur->Get_Weight();
			unsigned int a = spin[l_cur->Get_Start()->Get_Index()];
			unsigned int b =  spin[l_cur->Get_End()->Get_Index()];
			if (w > 0)
			{
				num_links_pos[a][b] += w;
				if (!is_directed && a != b) //Only one edge is defined in case it is undirected
					num_links_pos[b][a] += w;
			}
			else
			{
				num_links_neg[a][b] -= w;
				if (!is_directed && a != b) //Only one edge is defined in case it is undirected
					num_links_neg[b][a] -= w;
			}
			
			l_cur = iter_l.Next();
		} //while links
		
		#ifdef DEBUG
		printf("d_p: %f\n", d_p);
		printf("d_n: %f\n", d_n);
		#endif
	
		double expected = 0.0;
		double a = 0.0;
		double normal_a = 0.0;
		
		double delta, u_p, u_n;
		double max_expected, max_a;
		
		//We don't take into account the lambda or gamma for
		//computing the modularity and adhesion, since they
		//are then incomparable to other definitions.
		for (unsigned int i = 1; i <= q; i++)
		{
			for (unsigned int j = 1; j <= q; j++)
			{
				if (!is_directed && i == j)
					expected	= degree_community_pos_out[i] * degree_community_pos_in[j]/(m_p == 0 ? 1 : 2*m_p)
								- degree_community_neg_out[i] * degree_community_neg_in[j]/(m_n == 0 ? 1 : 2*m_n);
				else
					expected	= degree_community_pos_out[i] * degree_community_pos_in[j]/(m_p == 0 ? 1: m_p)
								- degree_community_neg_out[i] * degree_community_neg_in[j]/(m_n == 0 ? 1 : m_n);
				
				a			= (num_links_pos[i][j] - num_links_neg[i][j]) - expected;
				
				if (i == j) //cohesion
				{
					if (is_directed)
						delta	= d_p * csize[i] * (csize[i] - 1); //Maximum amount
					else
						delta	= d_p * csize[i] * (csize[i] - 1)/2; //Maximum amount
					
					u_p		= delta - num_links_pos[i][i]; //Add as many positive links we can
					u_n		= -num_links_neg[i][i]; //Delete as many negative links we can				
					Q	   += a;
				}
				else //adhesion
				{
					if (is_directed)
						delta	= d_n * csize[i] * csize[j]*2; //Maximum amount
					else
						delta	= d_n * csize[i] * csize[j]; //Maximum amount

					u_p		= -num_links_pos[i][j]; //Delete as many positive links we can
					u_n		= delta - num_links_neg[i][j]; //Add as many negative links we can
				}
				
				if (!is_directed && i == j)
					max_expected	= (degree_community_pos_out[i] + u_p) * (degree_community_pos_in[j] + u_p)/((m_p + u_p) == 0 ? 1 : 2*(m_p + u_p)) 
									- (degree_community_neg_out[i] - u_n) * (degree_community_neg_in[j] + u_n)/((m_n + u_n) == 0 ? 1 : 2*(m_n + u_n));
				else
					max_expected	= (degree_community_pos_out[i] + u_p) * (degree_community_pos_in[j] + u_p)/((m_p + u_p) == 0 ? 1 : m_p + u_p) 
									- (degree_community_neg_out[i] - u_n) * (degree_community_neg_in[j] + u_n)/((m_n + u_n) == 0 ? 1 : m_n + u_n);
				//printf("%f/%f %d/%d\t", num_links_pos[i][j], num_links_neg[i][j], csize[i], csize[j]);
				//printf("%f/%f - %f(%f)\t", u_p, u_n, expected, max_expected);
				max_a			= ((num_links_pos[i][j] + u_p) - (num_links_neg[i][j] + u_n)) - max_expected;				
				
				
				//In cases where we haven't actually found a ground state
				//the adhesion/cohesion *might* not be negative/positive,
				//hence the maximum adhesion and cohesion might behave quite
				//strangely. In order to prevent that, we limit them to 1 in
				//absolute value, and prevent from dividing by zero (even if
				//chuck norris would).
				if (i == j)
					normal_a = a/(max_a == 0 ? a : max_a);
				else
					normal_a = -a/(max_a == 0 ? a : max_a);
					
				if (normal_a > 1)
					normal_a = 1;
				else if (normal_a < -1)
					normal_a = -1;
					
				MATRIX(*adhesion, i - 1, j - 1) = a;
				MATRIX(*normalised_adhesion, i - 1, j - 1) = normal_a;
			} //for j
			//printf("\n");
		} //for i
		
		//free the allocated memory
		for( int i = 0 ; i < q+1 ; i++ )
		{
			delete [] num_links_pos[i] ;
			delete [] num_links_neg[i];
		}
		delete [] num_links_pos ;	
		delete [] num_links_neg ;
		
	} //adhesion
	
	if (modularity)  
	{ 
		if (is_directed)
			*modularity=Q/(m_p + m_n); 
		else
			*modularity=2*Q/(m_p + m_n); //Correction for the way m_p and m_n are counted. Modularity is 1/m, not 1/2m
	}
	
	if (polarization) 
	{ 
		double sum_ad = 0.0;
		for (unsigned int i = 0; i < q; i++)
		{
			for (unsigned int j = 0; j < q; j++)
			{
				if (i != j)
				{
					sum_ad -= MATRIX(*normalised_adhesion, i, j);
				}
			}
		}
		*polarization= sum_ad/(q*q - q);
	}
	#ifdef DEBUG	
	printf("Finished writing cluster.\n");
	#endif
	return num_nodes;
}
