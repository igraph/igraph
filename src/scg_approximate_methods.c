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

INT intervals_plus_kmeans(const REAL *v, UINT *gr, const UINT n,
							const UINT n_interv, const UINT maxiter)
{
	UINT i;
	INT converge;
	REAL *centers = real_vector(n_interv);
	breaks_computation(v,n,centers,n_interv,2);
	converge = kmeans_Lloyd(v, n, 1, centers, n_interv, (INT*) gr,maxiter);
	
	/*renumber the groups*/
	for(i=0; i<n; i++) gr[i] = gr[i]-1 + FIRST_GROUP_NB;

	free_real_vector(centers);

	return converge;
}
										
void intervals_method(const REAL *v, UINT *gr, const UINT n, const UINT n_interv)
{
	UINT i, lo, hi, new;
	const UINT lft = 1;
	const UINT include_border = 1;
			
	REAL *breaks = real_vector(n_interv+1);
	
	breaks_computation(v, n, breaks, n_interv+1, 1);

    for(i = 0; i < n; i++) {
	    lo = 0;
	    hi = n_interv;
	    if(v[i] <  breaks[lo] || breaks[hi] < v[i] ||
	       (v[i] == breaks[lft ? hi : lo] && !include_border)) ;
	    else {
			while(hi - lo >= 2) {
		    	new = (hi + lo)/2;
		    	if(v[i] > breaks[new] || (lft && v[i] == breaks[new]))
				lo = new;
		    	else
				hi = new;
			}
		gr[i] = lo + FIRST_GROUP_NB;
	    }
	}
	free_real_vector(breaks);
}

void breaks_computation(const REAL *v,const UINT n, REAL *breaks,
						const UINT nb,const UINT method)
{
	UINT i;
	REAL eps, vmin,vmax;
	vmin = min_real_vector(v,n);
	vmax = max_real_vector(v,n);
	
	if(vmax==vmin)
		error("There is only one (repeated) value in argument 'v'\
				of bin_size_computation()");
	if(nb<2)
		error("'nb' in bin_size_computation() must be >= 2");
			
	switch(method)
	{	//constant bins for fixed-size intervals method
		case 1:
			eps = (vmax-vmin)/(REAL)(nb-1);
			breaks[0] = vmin;
			for(i=1; i<nb-1; i++) breaks[i]=breaks[i-1]+eps;
			breaks[nb-1] = vmax;
		break;
		//equidistant centers for kmeans
		case 2:
			eps = (vmax-vmin)/(REAL)nb;
			breaks[0] = vmin + eps/2.;
			for(i=1; i<nb; i++) breaks[i] = breaks[i-1]+eps;
		break;
		//TODO: implement logarithmic binning for power-law-like distributions
		
		default:
			error("Choose a method to compute the breaks in breaks_computation():\
					1-constant bins (intervals method),2-equidistant centers");
	}
}
