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
int compare_ind_val(const void *a, const void *b)
{
	INDVAL *arg1 = (INDVAL *) a;
	INDVAL *arg2 = (INDVAL *) b;
	
	if( arg1->val < arg2->val ) return -1;
	else if( arg1->val == arg2->val ) return 0;
	else return 1;
}

/*to be used with qsort and struct groups*/  
int compare_groups(const void *a,const void *b)
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
int compare_real(const void *a, const void *b)
{
	REAL arg1 = * (REAL *) a;
	REAL arg2 = * (REAL *) b;
	
	if(arg1 < arg2) return -1;
	else if(arg1 == arg2) return 0;
	else return 1;
}

/*to be used with qsort and integer vectors */ 
int compare_int(const void *a, const void *b)
{
	INT arg1 = * (INT *) a;
	INT arg2 = * (INT *) b;
	return (arg1 -arg2);
}

/* allocate a REAL symmetrix matrix with dimension size x size in vector format*/ 
REAL *real_sym_matrix(const UINT size) 
{ 
	REAL *S;
	S = (REAL *) CALLOC(size*(size+1)/2,sizeof(REAL)); 
	if (!S) error("\nallocation failure in real_sym_matrix()\n"); 
	return S;
}

/* allocate a REAL matrix with dimension nrow x ncol*/ 
REAL **real_matrix(const UINT nrow, const UINT ncol) 
{ 
	UINT i;
	REAL **M; 
	//allocate pointers to rows
	M = (REAL **) CALLOC(nrow,sizeof(REAL*) );
	if (!M) error("row allocation failure in real_matrix()"); 
	for(i=0;i<nrow;i++){
		M[i] = (REAL *) CALLOC(ncol,sizeof(REAL) ); 
		if(!M[i]) error("column allocation failure in real_matrix()");
	}
	for(i=0;i<nrow;i++)
		M[i] = (REAL *) CALLOC(ncol,sizeof(REAL) );
	return M; 
}
void free_real_matrix(REAL **M,const UINT nrow)
{
	int i;
	for(i=0; i<nrow; i++)
		FREE(M[i]);
	FREE(M);
}

/* allocate an UINT matrix with dimension nrow x ncol*/
UINT **uint_matrix(const UINT nrow, const UINT ncol) 
{ 
	UINT i;
	UINT **M; 
	M = (UINT **) CALLOC(nrow, sizeof(UINT*));
	if (!M) error("row allocation failure in uint_matrix()"); 
	for(i=0;i<nrow;i++){
		M[i] = (UINT *) CALLOC(ncol, sizeof(UINT)); 
		if(!M[i]) error("column allocation failure in uint_matrix()");
	}
	for(i=0;i<nrow;i++)
		M[i] = (UINT *) CALLOC(ncol, sizeof(UINT));
	return M; 
}
void free_uint_matrix(UINT **M, const UINT nrow)
{
	UINT i;
	for(i=0; i<nrow; i++)
		FREE(M[i]);
	FREE(M);
}

/* allocate a REAL vector of size n*/ 
REAL *real_vector(const UINT n)
{ 
	REAL *v;
	v = (REAL *) CALLOC(n, sizeof(REAL)); 
	if (!v) error("allocation failure in real_vector()");
	return v;
}
REAL min_real_vector(const REAL *v,const UINT n)
{
	UINT i;
	REAL temp = v[0];
	for(i=1; i<n; i++)
		if( v[i] < temp) temp=v[i];
	return temp;
}
REAL max_real_vector(const REAL *v,const UINT n)
{
		UINT i;
		REAL temp = v[0];
		for(i=1; i<n; i++)
			if( v[i] > temp) temp=v[i];
		return temp;
}

/* allocate a UINT vector of size n*/
UINT *uint_vector(const UINT n) 
{ 
	UINT *v;
	v = (UINT *) CALLOC(n, sizeof(UINT));
	if (!v) error("allocation failure in uint_vector()"); 
	return v;
}
UINT min_uint_vector(const UINT *v, const UINT n)
{
	UINT i;
	UINT temp = v[0];
	for(i=1; i<n; i++)
		if( v[i] < temp) temp=v[i];
	return temp;
}
UINT max_uint_vector(const UINT *v, const UINT n)
{
	UINT i;
	UINT temp = v[0];
	for(i=1; i<n; i++)
		if( v[i] > temp) temp=v[i];
	return temp;
}

#ifndef R_COMPIL
//if defined R_COMPIL error() and warning() defined in <R.h>
	void error(const char error_text[])
	{
		fprintf(stderr,"\nRun-time error...\n"); 
		fprintf(stderr,"%s\n",error_text); 
		fprintf(stderr,"...now exiting to system...\n"); 
		exit(1); 
	}
	void warning(const char warning_text[]) 
	{
		printf("\nWarning: %s \n",warning_text); 
	}
#endif

#ifdef R_COMPIL
	void scg_r_wrapper(double *v, int *gr, int *n, int *nt,
					int *nev, int *nmatrix, int *nalgo, double *p, int *maxiter)
	{
		UINT i;
		REAL **M = (REAL**) CALLOC( (UINT)*nev, sizeof(REAL*) ); 
		if(!M) error("row allocation failure in real_matrix()"); 
		for(i=0;i<*nev;i++) M[i] = &v[(*n) * i];

		grouping(M,(UINT*)gr,(UINT)*n,(UINT*)nt,(UINT)*nev,
				(UINT)*nmatrix, (REAL*)p,(UINT)*nalgo,(UINT)*maxiter);
		FREE(M);
	}
#endif
