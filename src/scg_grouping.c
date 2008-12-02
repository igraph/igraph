/*
 *  Scglib : A C library for the spectral coarse graining of matrices
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
 *    The grouping function takes as argument 'nev' eigenvectors and 
 *	  and tries to minimize the eigenpair shifts induced by the coarse
 *	  graining (Section 5 of the above reference). The eigenvectors are
 *	  stored in a 'nev'x'n' matrix 'v'.
 *	  The 'algo' parameter can take the following values
 *		1  ->  Optimal method (sec. 5.3.1)
 *		2  ->  Intervals+k-means (sec. 5.3.3)
 *		3  ->  Intervals (sec. 5.3.2)
 *		4  ->  Exact SCG (sec. 5.4.1--last paragraph)
 *	  'nt' is a vector of length 'nev' giving either the size of the 
 *	  partitions (if algo = 1) or the number of intervals to cut the
 *	  eigenvectors if algo = 2 or algo = 3. When algo = 4 this parameter
 *	  is ignored. 'maxiter' fixes the maximum number of iterations of
 *	  the k-means algorithm, and is only considered when algo = 2.
 *	  All the algorithms try to find a minimizing partition of
 *	  ||v_i-Pv_i|| where P is a problem-specific projector and v_i denotes
 *	  the eigenvectors stored in v. The final partition is worked out
 *	  as decribed in Method 1 of Section 5.4.2.
 *	  'matrix' provides the type of SCG (i.e. the form of P). So far,
 *	  the options are those described in section 6, that is:
 *		1  ->  Symmetric (sec. 6.1)
 *		2  ->  Laplacian (sec. 6.2)
 *		3  ->  Stochastic (sec. 6.3)
 *	  In the stochastic case, a valid distribution probability 'p' must be
 *	  provided. In all other cases, 'p' is ignored and can be set to NULL.
 *	  The group labels in the final partition are given in 'gr' as positive
 *	  consecutive integers starting from 0.
 */

#include "scg_headers.h"
#include "error.h"

int igraph_scg_grouping(igraph_real_t **v, unsigned int *gr, 
			const unsigned int n, const unsigned int *nt, 
			const unsigned int nev,
			const unsigned int matrix, const igraph_real_t *p, 
			const unsigned int algo, const unsigned int maxiter)
{
	unsigned int i,j;
	unsigned int **gr_mat = igraph_uint_matrix(nev,n);
	unsigned int **gr_mat_t = igraph_uint_matrix(n,nev);
	
	switch (algo)
	{
		case 1:
			for(i=0; i<nev; i++)
				igraph_i_scg_optimal_partition(v[i], gr_mat[i], n, nt[i], matrix, p);              
			break;
			
		case 2:
			for(i=0; i<nev; i++){
				if(!igraph_i_scg_intervals_plus_kmeans(v[i], gr_mat[i], n, nt[i], maxiter))
				  IGRAPH_WARNING("kmeans did not converge");
			}
			break;
			
		case 3:
			for(i=0; i<nev; i++)
				igraph_i_scg_intervals_method(v[i], gr_mat[i], n, nt[i]);
			break;
			
		case 4:
			for(i=0; i<nev; i++)
				igraph_i_scg_exact_coarse_graining(v[i],gr_mat[i],n);
			break;	
	 
		default:
			igraph_free_uint_matrix(gr_mat, nev);
			igraph_free_uint_matrix(gr_mat_t, n);
			IGRAPH_ERROR("Choose a grouping method: 1-Optimal, 2-Fixed_size intervals+kmeans, "
				     "3-Fixed_size intervals, 4-Exact coarse graining", 
				     IGRAPH_EINVAL);
	}
	
	//If only one vector copy the groups and jump out
 	if(nev==1){
		for(i=0; i<n; i++)
			gr[i] = gr_mat[0][i];
			
		igraph_free_uint_matrix(gr_mat, nev);
		igraph_free_uint_matrix(gr_mat_t, n);
		
		return 0;
	}
	
	//Otherwise works out the final groups as decribed in section 5.4.2
		//First, works with the tranpose of gr_mat
	for(i=0; i<n; i++)
		for(j=0; j<nev; j++)
			gr_mat_t[i][j] = gr_mat[j][i];
	igraph_free_uint_matrix(gr_mat, nev);
	
		//Then computes the final groups. Use qsort for speed
	igraph_i_scg_groups_t *g = (igraph_i_scg_groups_t*)igraph_Calloc(n, igraph_i_scg_groups_t);
	for(i=0; i<n; i++){
		g[i].ind = i;
		g[i].n = nev;
		g[i].gr = gr_mat_t[i];
	}
		
	qsort(g, n, sizeof(igraph_i_scg_groups_t), igraph_i_scg_compare_groups);
	unsigned int gr_nb = FIRST_GROUP_NB;
	gr[g[0].ind] = gr_nb;
	for(i=1; i<n; i++){
		if(igraph_i_scg_compare_groups(&g[i], &g[i-1]) != 0) gr_nb++;
		gr[g[i].ind] = gr_nb;
	}
	igraph_Free(g);
	igraph_free_uint_matrix(gr_mat_t,n);
	
	return 0;
}



