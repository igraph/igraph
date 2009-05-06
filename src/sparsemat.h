/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2009  Gabor Csardi <csardi@rmki.kfki.hu>
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

#ifndef IGRAPH_SPARSEMAT_H
#define IGRAPH_SPARSEMAT_H

#include "cs/cs.h"
#include "types.h"

typedef struct igraph_sparsemat_t {
  cs *cs;
} igraph_sparsemat_t;

typedef enum { IGRAPH_SPARSEMAT_TRIPLET, 
	       IGRAPH_SPARSEMAT_CC        } igraph_sparsemat_type_t;

int igraph_sparsemat_init(igraph_sparsemat_t *A, int rows, int cols, int nzmax);
void igraph_sparsemat_destroy(igraph_sparsemat_t *A);
int igraph_sparsemat_realloc(igraph_sparsemat_t *A, int nzmax);

long int igraph_sparsemat_nrow(const igraph_sparsemat_t *A);
long int igraph_sparsemat_ncol(const igraph_sparsemat_t *B);
igraph_sparsemat_type_t igraph_sparsemat_type(const igraph_sparsemat_t *A);
igraph_bool_t igraph_sparsemat_is_triplet(const igraph_sparsemat_t *A);
igraph_bool_t igraph_sparsemat_is_cc(const igraph_sparsemat_t *A);

int igraph_sparsemat_entry(igraph_sparsemat_t *A, int row, int col, 
			   igraph_real_t elem);
int igraph_sparsemat_compress(const igraph_sparsemat_t *A, 
			      igraph_sparsemat_t *res);
int igraph_sparsemat_transpose(const igraph_sparsemat_t *A, 
			       igraph_sparsemat_t *res, int values);
int igraph_sparsemat_dupl(igraph_sparsemat_t *A);
int igraph_sparsemat_fkeep(igraph_sparsemat_t *A, 
			   int (*fkeep)(int, int, igraph_real_t, void*),
			   void *other);
int igraph_sparsemat_dropzeros(igraph_sparsemat_t *A);
int igraph_sparsemat_multiply(const igraph_sparsemat_t *A,
			      const igraph_sparsemat_t *B,
			      igraph_sparsemat_t *res);
int igraph_sparsemat_add(const igraph_sparsemat_t *A, 
			 const igraph_sparsemat_t *B,
			 igraph_real_t alpha,
			 igraph_real_t beta,
			 igraph_sparsemat_t *res);
int igraph_sparsemat_gaxpy(const igraph_sparsemat_t *A,
			   const igraph_vector_t *x,
			   igraph_vector_t *res);

int igraph_sparsemat_lsolve(const igraph_sparsemat_t *A,
			    const igraph_vector_t *b,
			    igraph_vector_t *res);
int igraph_sparsemat_ltsolve(const igraph_sparsemat_t *A,
			     const igraph_vector_t *b,
			     igraph_vector_t *res);
int igraph_sparsemat_usolve(const igraph_sparsemat_t *A,
			    const igraph_vector_t *b,
			    igraph_vector_t *res);
int igraph_sparsemat_utsolve(const igraph_sparsemat_t *A,
			     const igraph_vector_t *b,
			     igraph_vector_t *res);

int igraph_sparsemat_cholsol(const igraph_sparsemat_t *A,
			     const igraph_vector_t *b,
			     igraph_vector_t *res, 
			     int order);

int igraph_sprasemat_lusol(const igraph_sparsemat_t *A,
			   const igraph_vector_t *b,
			   igraph_vector_t *res,
			   int order,
			   igraph_real_t tol);
#endif
