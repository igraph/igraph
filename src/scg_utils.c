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
 
 #include "scg_headers.h"

/*to be used with qsort and struct ind_val arrays */  
int igraph_i_scg_compare_ind_val(const void *a, const void *b)
{
	igraph_i_scg_indval_t *arg1 = (igraph_i_scg_indval_t *) a;
	igraph_i_scg_indval_t *arg2 = (igraph_i_scg_indval_t *) b;
	
	if( arg1->val < arg2->val ) return -1;
	else if( arg1->val == arg2->val ) return 0;
	else return 1;
}

/*to be used with qsort and struct groups*/  
int igraph_i_scg_compare_groups(const void *a,const void *b)
{
		igraph_i_scg_groups_t *arg1 = (igraph_i_scg_groups_t *) a;
		igraph_i_scg_groups_t *arg2 = (igraph_i_scg_groups_t *) b;
		unsigned int i;
		for(i=0; i<arg1->n; i++){
			if(arg1->gr[i]>arg2->gr[i]) return 1;
			else if(arg1->gr[i]<arg2->gr[i]) return -1;
		}
		return 0;
}

/*to be used with qsort and real_vectors */  
int igraph_i_scg_compare_real(const void *a, const void *b)
{
	igraph_real_t arg1 = * (igraph_real_t *) a;
	igraph_real_t arg2 = * (igraph_real_t *) b;
	
	if(arg1 < arg2) return -1;
	else if(arg1 == arg2) return 0;
	else return 1;
}

/*to be used with qsort and integer vectors */ 
int igraph_i_scg_compare_int(const void *a, const void *b)
{
	int arg1 = * (int *) a;
	int arg2 = * (int *) b;
	return (arg1 -arg2);
}

/* allocate a igraph_real_t symmetrix matrix with dimension size x size in vector format*/ 
igraph_real_t *igraph_real_sym_matrix(const unsigned int size) 
{ 
	igraph_real_t *S;
	S = (igraph_real_t *) igraph_Calloc(size*(size+1)/2,igraph_real_t); 
	if (!S) IGRAPH_ERROR("allocation failure in real_sym_matrix()", 
			    IGRAPH_EINVAL); 
	return S;
}

/* allocate a igraph_real_t matrix with dimension nrow x ncol*/ 
igraph_real_t **igraph_real_matrix(const unsigned int nrow, const unsigned int ncol) 
{ 
	unsigned int i;
	igraph_real_t **M; 
	//allocate pointers to rows
	M = (igraph_real_t **) igraph_Calloc(nrow,igraph_real_t*);
	if (!M) 
	  IGRAPH_ERROR("row allocation failure in real_matrix()", IGRAPH_EINVAL); 
	for(i=0;i<nrow;i++){
	  M[i] = (igraph_real_t *) igraph_Calloc(ncol,igraph_real_t); 
	  if(!M[i])
	    IGRAPH_ERROR("column allocation failure in real_matrix()", 
			 IGRAPH_EINVAL);
	}
	for(i=0;i<nrow;i++)
	  M[i] = (igraph_real_t *) igraph_Calloc(ncol,igraph_real_t);
	return M; 
}
void igraph_free_real_matrix(igraph_real_t **M,const unsigned int nrow)
{
	int i;
	for(i=0; i<nrow; i++)
		igraph_Free(M[i]);
	igraph_Free(M);
}

/* allocate an unsigned int matrix with dimension nrow x ncol*/
unsigned int **igraph_uint_matrix(const unsigned int nrow, const unsigned int ncol) 
{ 
	unsigned int i;
	unsigned int **M; 
	M = (unsigned int **) igraph_Calloc(nrow, unsigned int*);
	if (!M) 
	  IGRAPH_ERROR("row allocation failure in uint_matrix()",
		       IGRAPH_EINVAL); 
	for(i=0;i<nrow;i++){
	  M[i] = (unsigned int *) igraph_Calloc(ncol, unsigned int); 
	  if(!M[i]) IGRAPH_ERROR("column allocation failure in uint_matrix()", 
				 IGRAPH_EINVAL);
	}
	for(i=0;i<nrow;i++)
		M[i] = (unsigned int *) igraph_Calloc(ncol, unsigned int);
	return M; 
}
void igraph_free_uint_matrix(unsigned int **M, const unsigned int nrow)
{
	unsigned int i;
	for(i=0; i<nrow; i++)
		igraph_Free(M[i]);
	igraph_Free(M);
}

/* allocate a igraph_real_t vector of size n*/ 
igraph_real_t *igraph_real_vector(const unsigned int n)
{ 
	igraph_real_t *v;
	v = (igraph_real_t *) igraph_Calloc(n, igraph_real_t); 
	if (!v) 
	  IGRAPH_ERROR("allocation failure in real_vector()", 
		       IGRAPH_EINVAL);
	return v;
}
igraph_real_t igraph_min_real_vector(const igraph_real_t *v,const unsigned int n)
{
	unsigned int i;
	igraph_real_t temp = v[0];
	for(i=1; i<n; i++)
		if( v[i] < temp) temp=v[i];
	return temp;
}
igraph_real_t igraph_max_real_vector(const igraph_real_t *v,const unsigned int n)
{
		unsigned int i;
		igraph_real_t temp = v[0];
		for(i=1; i<n; i++)
			if( v[i] > temp) temp=v[i];
		return temp;
}

/* allocate a unsigned int vector of size n*/
unsigned int *igraph_uint_vector(const unsigned int n) 
{ 
	unsigned int *v;
	v = (unsigned int *) igraph_Calloc(n, unsigned int);
	if (!v) IGRAPH_ERROR("allocation failure in uint_vector()", 
			     IGRAPH_EINVAL); 
	return v;
}
unsigned int igraph_min_uint_vector(const unsigned int *v, const unsigned int n)
{
	unsigned int i;
	unsigned int temp = v[0];
	for(i=1; i<n; i++)
		if( v[i] < temp) temp=v[i];
	return temp;
}
unsigned int igraph_max_uint_vector(const unsigned int *v, const unsigned int n)
{
	unsigned int i;
	unsigned int temp = v[0];
	for(i=1; i<n; i++)
		if( v[i] > temp) temp=v[i];
	return temp;
}

/* #ifdef R_COMPIL */
/* 	void scg_r_wrapper(double *v, int *gr, int *n, int *nt, */
/* 					int *nev, int *nmatrix, int *nalgo, double *p, int *maxiter) */
/* 	{ */
/* 		unsigned int i; */
/* 		igraph_real_t **M = (igraph_real_t**) igraph_Calloc((unsigned int)*nev, igraph_real_t* );  */
/* 		if(!M)  */
/* 		  IGRAPH_ERROR("row allocation failure in real_matrix()",  */
/* 			       IGRAPH_EINVAL);  */
/* 		for(i=0;i<*nev;i++) M[i] = &v[(*n) * i]; */

/* 		grouping(M,(unsigned int*)gr,(unsigned int)*n,(unsigned int*)nt,(unsigned int)*nev, */
/* 				(unsigned int)*nmatrix, (igraph_real_t*)p,(unsigned int)*nalgo,(unsigned int)*maxiter); */
/* 		igraph_Free(M); */
/* 	} */
/* #endif */
