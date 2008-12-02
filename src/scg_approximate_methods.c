/*
 *  SCGlib : A C library for the spectral coarse graining of matrices
 *	as described in the paper: Shrinking Matrices while preserving their
 *	eigenpairs with Application to the Spectral Coarse Graining of Graphs.
 *	Preprint available at <http://people.epfl.ch/david.morton>
 *  
 *	Copyright (C) 2008 David Morton de Lachapelle <david.morton@a3.epfl.ch>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301 USA
 *
 *  DESCRIPTION
 *	-----------
 *    The intervals_method and intervals_plus_kmeans implements the 
 *	  methods of sec. 5.3.2 and sec. 5.3.3 of the above reference.
 *	  They take an eigenvector 'v' as parameter and a vector 'breaks'
 *	  of length 'nb', which provide the intervals used to cut 'v'.
 *	  Then all components of 'v' that fall into the same interval are
 *	  assigned the same group label in 'gr'. The group labels are
 *	  positive consecutive integers starting from 0.
 *	  The intervals_method function is adapted from bincode of the R
 *	  base package.
 *	  The intervals_plus_kmeans is initialized with regularly-spaced
 *	  breaks, which rougly corresponds to the intervals_method. Then
 *	  kmeans minimizes iteratively the objective function until it gets
 *	  stuck in a (usually) local minimum, or until 'itermax' is reached.
 *	  So far, the breaks_computation function allows computation of 
 *	  constant bins, as used in intervals_method, and of equidistant
 *	  centers as used in intervals_plus_kmeans.
 */

#include "scg_headers.h"

int igraph_i_scg_intervals_plus_kmeans(const igraph_vector_t *v, 
				       igraph_vector_long_t *gr,
				       const unsigned int n,
				       const unsigned int n_interv,
				       const unsigned int maxiter)
{
	unsigned int i;
	int converge;
	igraph_vector_t centers;
	igraph_vector_init(&centers, n_interv);
	igraph_i_scg_breaks_computation(v,n,&centers,n_interv,2);
	converge = igraph_i_scg_kmeans_Lloyd(v, n, 1, &centers, n_interv, gr, maxiter);
	
	/*renumber the groups*/
	for(i=0; i<n; i++) VECTOR(*gr)[i] = VECTOR(*gr)[i]-1 + FIRST_GROUP_NB;

	igraph_vector_destroy(&centers);

	return converge;
}
										
void igraph_i_scg_intervals_method(const igraph_vector_t *v, 
				   igraph_vector_long_t *gr, 
				   const unsigned int n, 
				   const unsigned int n_interv)
{
	unsigned int i, lo, hi, new;
	const unsigned int lft = 1;
	const unsigned int include_border = 1;
			
	igraph_vector_t breaks;
	igraph_vector_init(&breaks, n_interv+1);
	
	igraph_i_scg_breaks_computation(v, n, &breaks, n_interv+1, 1);

	for(i = 0; i < n; i++) {
	    lo = 0;
	    hi = n_interv;
	    if(VECTOR(*v)[i] <  VECTOR(breaks)[lo] || VECTOR(breaks)[hi] < VECTOR(*v)[i] ||
	       (VECTOR(*v)[i] == VECTOR(breaks)[lft ? hi : lo] && !include_border)) ;
	    else {
			while(hi - lo >= 2) {
			  new = (hi + lo)/2;
			  if(VECTOR(*v)[i] > VECTOR(breaks)[new] || (lft && VECTOR(*v)[i] == VECTOR(breaks)[new]))
				lo = new;
		    	else
				hi = new;
			}
			VECTOR(*gr)[i] = lo + FIRST_GROUP_NB;
	    }
	}
    igraph_vector_destroy(&breaks);
}

int igraph_i_scg_breaks_computation(const igraph_vector_t *v,
				    const unsigned int n, 
				    igraph_vector_t *breaks,
				    const unsigned int nb,
				    const unsigned int method)
{
	unsigned int i;
	igraph_real_t eps, vmin,vmax;
	vmin = igraph_vector_min(v);
	vmax = igraph_vector_max(v);
	
	if(vmax==vmin)
	  IGRAPH_ERROR("There is only one (repeated) value in argument 'v' "
		       "of bin_size_computation()", IGRAPH_EINVAL);
	if(nb<2)
	  IGRAPH_ERROR("'nb' in bin_size_computation() must be >= 2", 
		       IGRAPH_EINVAL);
			
	switch(method)
	{	//constant bins for fixed-size intervals method
		case 1:
			eps = (vmax-vmin)/(igraph_real_t)(nb-1);
			VECTOR(*breaks)[0] = vmin;
			for(i=1; i<nb-1; i++) VECTOR(*breaks)[i]=VECTOR(*breaks)[i-1]+eps;
			VECTOR(*breaks)[nb-1] = vmax;
		break;
		//equidistant centers for kmeans
		case 2:
			eps = (vmax-vmin)/(igraph_real_t)nb;
			VECTOR(*breaks)[0] = vmin + eps/2.;
			for(i=1; i<nb; i++) VECTOR(*breaks)[i] = VECTOR(*breaks)[i-1]+eps;
		break;
		//TODO: implement logarithmic binning for power-law-like distributions
		
		default:
		  IGRAPH_ERROR("Choose a method to compute the breaks in breaks_computation(): "
			       "1-constant bins (intervals method),2-equidistant centers", 
			       IGRAPH_EINVAL);
	}
	return 0;
}
