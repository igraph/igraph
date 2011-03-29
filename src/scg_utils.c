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
 *    This files contains the data structures and error handing
 *	  functions used throughout the SCGlib.  
 */

#include "igraph_error.h" 
 
#include "scg_headers.h"

/*to be used with qsort and struct ind_val arrays */  
int igraph_i_compare_ind_val(const void *a, const void *b)
{
	INDVAL *arg1 = (INDVAL *) a;
	INDVAL *arg2 = (INDVAL *) b;
	
	if( arg1->val < arg2->val ) return -1;
	else if( arg1->val == arg2->val ) return 0;
	else return 1;
}

/*to be used with qsort and struct groups*/  
int igraph_i_compare_groups(const void *a,const void *b)
{
		GROUPS *arg1 = (GROUPS *) a;
		GROUPS *arg2 = (GROUPS *) b;
		UINT i;
		for(i=0; i<arg1->n; i++){
			if(arg1->gr[i]>arg2->gr[i]) return 1;
			else if(arg1->gr[i]<arg2->gr[i]) return -1;
		}
		return 0;
}

/*to be used with qsort and real_vectors */  
int igraph_i_compare_real(const void *a, const void *b)
{
	REAL arg1 = * (REAL *) a;
	REAL arg2 = * (REAL *) b;
	
	if(arg1 < arg2) return -1;
	else if(arg1 == arg2) return 0;
	else return 1;
}

/*to be used with qsort and integer vectors */ 
int igraph_i_compare_int(const void *a, const void *b)
{
	INT arg1 = * (INT *) a;
	INT arg2 = * (INT *) b;
	return (arg1 -arg2);
}

/* allocate a REAL symmetrix matrix with dimension size x size in vector format*/ 
REAL *igraph_i_real_sym_matrix(const UINT size) 
{ 
	REAL *S;
	S = (REAL *) CALLOC(size*(size+1)/2,sizeof(REAL)); 
	if (!S) igraph_error("\nallocation failure in real_sym_matrix()\n",
			     __FILE__, __LINE__, IGRAPH_ENOMEM); 
	return S;
}

/* allocate a REAL matrix with dimension nrow x ncol*/ 
REAL **igraph_i_real_matrix(const UINT nrow, const UINT ncol) 
{ 
	UINT i;
	REAL **M; 
	//allocate pointers to rows
	M = (REAL **) CALLOC(nrow,sizeof(REAL*) );
	if (!M) igraph_error("row allocation failure in real_matrix()", 
			     __FILE__, __LINE__, IGRAPH_ENOMEM); 
	for(i=0;i<nrow;i++){
		M[i] = (REAL *) CALLOC(ncol,sizeof(REAL) ); 
		if(!M[i]) igraph_error("column allocation failure in real_matrix()", __FILE__, __LINE__, IGRAPH_ENOMEM);
	}
	for(i=0;i<nrow;i++)
		M[i] = (REAL *) CALLOC(ncol,sizeof(REAL) );
	return M; 
}
void igraph_i_free_real_matrix(REAL **M,const UINT nrow)
{
	int i;
	for(i=0; i<nrow; i++)
		FREE(M[i]);
	FREE(M);
}

/* allocate an UINT matrix with dimension nrow x ncol*/
UINT **igraph_i_uint_matrix(const UINT nrow, const UINT ncol) 
{ 
	UINT i;
	UINT **M; 
	M = (UINT **) CALLOC(nrow, sizeof(UINT*));
	if (!M) igraph_error("row allocation failure in uint_matrix()", 
			     __FILE__, __LINE__, IGRAPH_ENOMEM); 
	for(i=0;i<nrow;i++){
		M[i] = (UINT *) CALLOC(ncol, sizeof(UINT)); 
		if(!M[i]) igraph_error("column allocation failure in uint_matrix()", __FILE__, __LINE__, IGRAPH_ENOMEM);
	}
	for(i=0;i<nrow;i++)
		M[i] = (UINT *) CALLOC(ncol, sizeof(UINT));
	return M; 
}
void igraph_i_free_uint_matrix(UINT **M, const UINT nrow)
{
	UINT i;
	for(i=0; i<nrow; i++)
		FREE(M[i]);
	FREE(M);
}

/* allocate a REAL vector of size n*/ 
REAL *igraph_i_real_vector(const UINT n)
{ 
	REAL *v;
	v = (REAL *) CALLOC(n, sizeof(REAL)); 
	if (!v) igraph_error("allocation failure in real_vector()", 
			     __FILE__, __LINE__, IGRAPH_ENOMEM);
	return v;
}
REAL igraph_i_min_real_vector(const REAL *v,const UINT n)
{
	UINT i;
	REAL temp = v[0];
	for(i=1; i<n; i++)
		if( v[i] < temp) temp=v[i];
	return temp;
}
REAL igraph_i_max_real_vector(const REAL *v,const UINT n)
{
		UINT i;
		REAL temp = v[0];
		for(i=1; i<n; i++)
			if( v[i] > temp) temp=v[i];
		return temp;
}

/* allocate a UINT vector of size n*/
UINT *igraph_i_uint_vector(const UINT n) 
{ 
	UINT *v;
	v = (UINT *) CALLOC(n, sizeof(UINT));
	if (!v) igraph_error("allocation failure in uint_vector()", 
			     __FILE__, __LINE__, IGRAPH_ENOMEM); 
	return v;
}
UINT igraph_i_min_uint_vector(const UINT *v, const UINT n)
{
	UINT i;
	UINT temp = v[0];
	for(i=1; i<n; i++)
		if( v[i] < temp) temp=v[i];
	return temp;
}
UINT igraph_i_max_uint_vector(const UINT *v, const UINT n)
{
	UINT i;
	UINT temp = v[0];
	for(i=1; i<n; i++)
		if( v[i] > temp) temp=v[i];
	return temp;
}

