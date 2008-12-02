/*
 *  SCGlib : A C library for the spectral coarse graining of matrices
 *	as described in the paper: Shrinking Matrices while preserving their
 *	eigenpairs with Application to the Spectral Coarse Graining of Graphs.
 *	Preprint available at <http://people.epfl.ch/david.morton>.
 *  
 *	Copyright (C) 2008 David Morton de Lachapelle
 *	<david.mortondelachapelle@swissquote.ch>
 *  <david.morton@epfl.ch>
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
 *    This is to show the SCGlib in action through a simple example.
 *	  Three vectors to be preserved are drawn at random and placed
 *	  row-wise in a matrix 'v'. The partition obtained by each grouping 
 *	  method is displayed.
 *	  
 *	  Use dedicated libraries (arpack, linpack, lapack, etc) to compute
 *	  actual eigenvectors/eigenvalues of matrices
 */

#include "igraph.h"

#include <stdlib.h>

#define RDN (rand() % 10000 + 1)/10000.

int main (unsigned int argc, const char * argv[]) {
	
	unsigned int n = 30;
	unsigned int nev = 3;
	igraph_vector_t nt; 
	igraph_real_t ntv[] = {3,3,2};
	//Type of SCG according to section 6 of the above reference
	//1:Symmetric (sec. 6.1); 2:Laplacian (sec. 6.2); 3:Stochastic (sec. 6.3)
	unsigned int matrix = IGRAPH_SCG_MATRIX_SYMMETRIC;
	
	igraph_matrix_t v;
	igraph_vector_t gr;

	unsigned int i,j;

	unsigned int algo;
	unsigned int maxiter = 100; //ignored when algo not equal to 2

	igraph_vector_view(&nt, ntv, sizeof(ntv)/sizeof(igraph_real_t));
	igraph_vector_init(&gr, n);
	igraph_matrix_init(&v, n, nev);
	
	srand(10);
	for(i=0;i<nev;i++)
		for(j=0;j<n;j++)
		       MATRIX(v,j,i) = RDN;
	
	//Algorithm used in the coarse graining:
	//1:Optimal method (sec. 5.3.1); 2:Intervals+k-means (sec. 5.3.3);
	//3:Intervals (sec. 5.3.2); 4:Exact SCG (sec. 5.4.1--last paragraph)
	for(algo=1; algo<=4; algo++){
	        igraph_scg_grouping(&v, &gr, &nt, matrix, /*p=*/ 0, algo, maxiter);
		printf("\nAlgo %i, %i groups:", algo, (int)igraph_vector_max(&gr)+1);
		for(i=0; i<n; i++)
		  printf(" %li", (long int)VECTOR(gr)[i]);
		printf("\n");
	}
	printf("\n");
				
	igraph_vector_destroy(&gr);
	igraph_matrix_destroy(&v);
	//free_real_vector(p);
	return 0;
}

