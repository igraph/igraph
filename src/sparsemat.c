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

#include "sparsemat.h"
#include "error.h"

int igraph_sparsemat_init(igraph_sparsemat_t *A, int rows, int cols, int nzmax) {

  if (rows < 0) { 
    IGRAPH_ERROR("Negative number of rows", IGRAPH_EINVAL);
  }
  if (cols < 0) {
    IGRAPH_ERROR("Negative number of columns", IGRAPH_EINVAL);
  }
  
  A->cs=cs_spalloc( rows, cols, nzmax, /*values=*/ 1, 
		  /*triplet=*/ 1);
  if (!A->cs) {
    IGRAPH_ERROR("Cannot allocate memory for sparse matrix", IGRAPH_ENOMEM);
  }

  return 0;
}

void igraph_sparsemat_destroy(igraph_sparsemat_t *A) {
  cs_spfree(A->cs);
}

int igraph_sparsemat_realloc(igraph_sparsemat_t *A, int nzmax) {
  return !cs_sprealloc(A->cs, nzmax);
}

long int igraph_sparsemat_nrow(const igraph_sparsemat_t *A) {
  return A->cs->m;
}

long int igraph_sparsemat_ncol(const igraph_sparsemat_t *A) {
  return A->cs->n;
}

igraph_sparsemat_type_t igraph_sparsemat_type(const igraph_sparsemat_t *A) {
  return A->cs->nz < 0 ? IGRAPH_SPARSEMAT_CC : IGRAPH_SPARSEMAT_TRIPLET;
}

igraph_bool_t igraph_sparsemat_is_triplet(const igraph_sparsemat_t *A) {
  return A->cs->nz >= 0;
}

igraph_bool_t igraph_sparsemat_is_cc(const igraph_sparsemat_t *A) {
  return A->cs->nz < 0;
}

/* Limitations:
   - p and q must be of length nrow and ncol, respectively
   - no duplicated indices
*/

int igraph_sparsemat_index(const igraph_sparsemat_t *A,
			   const igraph_vector_int_t *p, 
			   const igraph_vector_int_t *q,
			   igraph_sparsemat_t *res) {

  long int nrow=A->cs->m, ncol=A->cs->n;
  long int plen=igraph_vector_int_size(p);
  igraph_vector_int_t pinv;
  long int i;
  
  /* We invert the permutation by hand */
  IGRAPH_CHECK(igraph_vector_int_init(&pinv, nrow));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &pinv);
  for (i=0; i<nrow; i++) {
    VECTOR(pinv)[ VECTOR(*p)[i] ] = i;
  }
  
  /* And call the permutation routine */
  if (! (res->cs = cs_permute(A->cs, VECTOR(pinv), VECTOR(*q), /*values=*/ 1))) {
    IGRAPH_ERROR("Cannot index sparse matrix", IGRAPH_FAILURE);
  }
  
  igraph_vector_int_destroy(&pinv);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

int igraph_sparsemat_entry(igraph_sparsemat_t *A, int row, int col, 
			   igraph_real_t elem) {
  
  if (!cs_entry(A->cs, row, col, elem)) {
    IGRAPH_ERROR("Cannot add entry to sparse matrix", 
		 IGRAPH_FAILURE);
  }

  return 0;
}

int igraph_sparsemat_compress(const igraph_sparsemat_t *A, 
			      igraph_sparsemat_t *res) {

  if (! (res->cs=cs_compress(A->cs)) ) {
    IGRAPH_ERROR("Cannot compress sparse matrix", IGRAPH_FAILURE);
  }

  return 0;
}


int igraph_sparsemat_transpose(const igraph_sparsemat_t *A, 
			       igraph_sparsemat_t *res, 
			       int values) {

  if (! (res->cs=cs_transpose(A->cs, values)) ) {
    IGRAPH_ERROR("Cannot transpose sparse matrix", IGRAPH_FAILURE);
  }

  return 0;
}


int igraph_sparsemat_dupl(igraph_sparsemat_t *A) {

  if (!cs_dupl(A->cs)) {
    IGRAPH_ERROR("Cannot transpose sparse matrix", IGRAPH_FAILURE);
  }
  
  return 0;
}


int igraph_sparsemat_fkeep(igraph_sparsemat_t *A, 
			   int (*fkeep)(int, int, igraph_real_t, void*),
			   void *other) {
  
  if (!cs_fkeep(A->cs, fkeep, other)) {
    IGRAPH_ERROR("Cannot filter sparse matrix", IGRAPH_FAILURE);
  }

  return 0;
}


int igraph_sparsemat_dropzeros(igraph_sparsemat_t *A) {

  if (!cs_dropzeros(A->cs)) {
    IGRAPH_ERROR("Cannot drop zeros from sparse matrix", IGRAPH_FAILURE);
  }

  return 0;
}


int igraph_sparsemat_multiply(const igraph_sparsemat_t *A,
			      const igraph_sparsemat_t *B,
			      igraph_sparsemat_t *res) {

  if (! (res->cs=cs_multiply(A->cs, B->cs))) {
    IGRAPH_ERROR("Cannot multiply matrices", IGRAPH_FAILURE);
  }

  return 0;
}


int igraph_sparsemat_add(const igraph_sparsemat_t *A, 
			 const igraph_sparsemat_t *B,
			 igraph_real_t alpha,
			 igraph_real_t beta,
			 igraph_sparsemat_t *res) {

  if (! (res->cs=cs_add(A->cs, B->cs, alpha, beta))) {
    IGRAPH_ERROR("Cannot add matrices", IGRAPH_FAILURE);
  }
  
  return 0;
}

int igraph_sparsemat_gaxpy(const igraph_sparsemat_t *A,
			   const igraph_vector_t *x,
			   igraph_vector_t *res) {

  if (A->cs->n != igraph_vector_size(x) || 
      A->cs->m != igraph_vector_size(res)) {
    IGRAPH_ERROR("Invalid matrix/vector size for multiplication",
		 IGRAPH_EINVAL);
  }
  
  if (! (cs_gaxpy(A->cs, VECTOR(*x), VECTOR(*res)))) {
    IGRAPH_ERROR("Cannot perform sparse matrix vector multiplication",
		 IGRAPH_FAILURE);
  }
  
  return 0;
}

int igraph_sparsemat_lsolve(const igraph_sparsemat_t *A,
			    const igraph_vector_t *b,
			    igraph_vector_t *res) {

  if (A->cs->m != A->cs->n) {
    IGRAPH_ERROR("Cannot perform lower triangular solve", IGRAPH_NONSQUARE);
  }

  if (res != b) {
    IGRAPH_CHECK(igraph_vector_update(res, b));
  }

  if (! cs_lsolve(A->cs, VECTOR(*res))) {
    IGRAPH_ERROR("Cannot perform lower triangular solve", IGRAPH_FAILURE);
  }
  
  return 0;
}

int igraph_sparsemat_ltsolve(const igraph_sparsemat_t *A,
			     const igraph_vector_t *b,
			     igraph_vector_t *res) {
  
  if (A->cs->m != A->cs->n) {
    IGRAPH_ERROR("Cannot perform transposed lower triangular solve",
		 IGRAPH_NONSQUARE);
  }
  
  if (res != b) {
    IGRAPH_CHECK(igraph_vector_update(res,b));
  }

  if (!cs_ltsolve(A->cs, VECTOR(*res))) {
    IGRAPH_ERROR("Cannot perform lower triangular solve", IGRAPH_FAILURE);
  }
  
  return 0;
}

int igraph_sparsemat_usolve(const igraph_sparsemat_t *A,
			    const igraph_vector_t *b,
			    igraph_vector_t *res) {

  if (A->cs->m != A->cs->n) {
    IGRAPH_ERROR("Cannot perform upper triangular solve", IGRAPH_NONSQUARE);
  }

  if (res != b) {
    IGRAPH_CHECK(igraph_vector_update(res, b));
  }

  if (! cs_usolve(A->cs, VECTOR(*res))) {
    IGRAPH_ERROR("Cannot perform upper triangular solve", IGRAPH_FAILURE);
  }
  
  return 0;
}

int igraph_sparsemat_utsolve(const igraph_sparsemat_t *A,
			     const igraph_vector_t *b,
			     igraph_vector_t *res) {
  
  if (A->cs->m != A->cs->n) {
    IGRAPH_ERROR("Cannot perform transposed upper triangular solve",
		 IGRAPH_NONSQUARE);
  }

  if (res != b) { 
    IGRAPH_CHECK(igraph_vector_update(res,b));
  }

  if (!cs_utsolve(A->cs, VECTOR(*res))) {
    IGRAPH_ERROR("Cannot perform transposed upper triangular solve", 
		 IGRAPH_FAILURE);
  }
  
  return 0;
}

int igraph_sparsemat_cholsol(const igraph_sparsemat_t *A,
			     const igraph_vector_t *b,
			     igraph_vector_t *res, 
			     int order) {
  
  if (A->cs->m != A->cs->n) {
    IGRAPH_ERROR("Cannot perform sparse symmetric solve",
		 IGRAPH_NONSQUARE);
  }

  if (res != b) { 
    IGRAPH_CHECK(igraph_vector_update(res,b));
  }

  if (! cs_cholsol(order, A->cs, VECTOR(*res))) {
    IGRAPH_ERROR("Cannot perform sparse symmetric solve", IGRAPH_FAILURE);
  }
  
  return 0;
}

int igraph_sparsemat_lusol(const igraph_sparsemat_t *A,
			   const igraph_vector_t *b,
			   igraph_vector_t *res,
			   int order,
			   igraph_real_t tol) {
  
  if (A->cs->m != A->cs->n) {
    IGRAPH_ERROR("Cannot perform LU solve",
		 IGRAPH_NONSQUARE);
  }

  if (res != b) { 
    IGRAPH_CHECK(igraph_vector_update(res,b));
  }

  if (! cs_lusol(order, A->cs, VECTOR(*res), tol)) {
    IGRAPH_ERROR("Cannot perform LU solve", IGRAPH_FAILURE);
  }
  
  return 0;
}

int igraph_i_sparsemat_cc(igraph_t *graph, const igraph_sparsemat_t *A,
			  igraph_bool_t directed) {

  igraph_vector_t edges;
  long int no_of_nodes=A->cs->m;
  long int no_of_edges=A->cs->nzmax;
  int *p=A->cs->p;
  int *i=A->cs->i;
  long int from=0;
  long int to=0;
  long int e=0;
  
  if (no_of_nodes != A->cs->n) {
    IGRAPH_ERROR("Cannot create graph object", IGRAPH_NONSQUARE);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges*2);
  
  while (*p < no_of_edges) {
    while (to < *(p+1)) {
      VECTOR(edges)[e++] = from;
      VECTOR(edges)[e++] = (*i);
      to++;
      i++;      
    }
    from++;
    p++;
  }

  IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_i_sparsemat_triplet(igraph_t *graph, const igraph_sparsemat_t *A,
			       igraph_bool_t directed) {

  igraph_vector_t edges;
  long int no_of_nodes=A->cs->m;
  long int no_of_edges=A->cs->nz;
  int *i=A->cs->p;
  int *j=A->cs->i;
  long int e;
  
  if (no_of_nodes != A->cs->n) {
    IGRAPH_ERROR("Cannot create graph object", IGRAPH_NONSQUARE);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges*2);
  
  for (e=0; e<2*no_of_edges; i++, j++) {
    VECTOR(edges)[e++] = (*i);
    VECTOR(edges)[e++] = (*j);
  }
  
  IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_sparsemat(igraph_t *graph, const igraph_sparsemat_t *A,
		     igraph_bool_t directed) {
  
  if (A->cs->nz < 0) {
    return(igraph_i_sparsemat_cc(graph, A, directed));
  } else {
    return(igraph_i_sparsemat_triplet(graph, A, directed));
  }
}

int igraph_get_sparsemat(const igraph_t *graph, igraph_sparsemat_t *res) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  long int i;
  
  IGRAPH_CHECK(igraph_sparsemat_init(res, no_of_nodes, no_of_nodes, 
				     no_of_edges));
  for (i=0; i<no_of_edges; i++) {
    IGRAPH_CHECK(igraph_sparsemat_entry(res, 
					IGRAPH_FROM(graph, i), 
					IGRAPH_TO(graph, i),
					1.0));
  }
  
  return 0;
}
