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

int igraph_sparsemat_permute(const igraph_sparsemat_t *A,
			     const igraph_vector_int_t *p, 
			     const igraph_vector_int_t *q,
			     igraph_sparsemat_t *res) {

  long int nrow=A->cs->m, ncol=A->cs->n;
  long int plen=igraph_vector_int_size(p);
  igraph_vector_int_t pinv;
  long int i;

  if (nrow != igraph_vector_int_size(p)) {
    IGRAPH_ERROR("Invalid row permutation length", IGRAPH_FAILURE);
  }
  if (ncol != igraph_vector_int_size(q)) {
    IGRAPH_ERROR("Invalud column permutation length", IGRAPH_FAILURE);
  }

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

int igraph_sparsemat_index(const igraph_sparsemat_t *A,
			   const igraph_vector_int_t *p,
			   const igraph_vector_int_t *q,
			   igraph_sparsemat_t *res,
			   igraph_real_t *constres) {

  igraph_sparsemat_t II, JJ, II2, JJ2, tmp;
  long int nrow=A->cs->m;
  long int ncol=A->cs->n;
  long int idx_rows=igraph_vector_int_size(p);
  long int idx_cols=igraph_vector_int_size(q);
  long int k;

  igraph_sparsemat_t *myres=res, mres;

  if (!res && (idx_rows != 1 || idx_cols != 1)) {
    IGRAPH_ERROR("Sparse matrix indexing: must give `res' if not a "
		 "single element is selected", IGRAPH_EINVAL);
  }

  if (!res) {
    myres=&mres;
  }

  /* Create first index matrix */
  IGRAPH_CHECK(igraph_sparsemat_init(&II2, idx_rows, nrow, idx_rows));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &II2);
  for (k=0; k<idx_rows; k++) {
    igraph_sparsemat_entry(&II2, k, VECTOR(*p)[k], 1.0);
  }
  IGRAPH_CHECK(igraph_sparsemat_compress(&II2, &II));
  igraph_sparsemat_destroy(&II2);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &II);

  /* Create second index matrix */
  IGRAPH_CHECK(igraph_sparsemat_init(&JJ2, ncol, idx_cols, idx_cols));
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &JJ2);
  for (k=0; k<idx_cols; k++) {
    igraph_sparsemat_entry(&JJ2, VECTOR(*q)[k], k, 1.0);
  }
  IGRAPH_CHECK(igraph_sparsemat_compress(&JJ2, &JJ));
  igraph_sparsemat_destroy(&JJ2);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &JJ);
  
  /* Multiply */
  IGRAPH_CHECK(igraph_sparsemat_multiply(&II, A, &tmp));
  igraph_sparsemat_destroy(&II);
  IGRAPH_FINALLY_CLEAN(1);
  IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmp);
  IGRAPH_CHECK(igraph_sparsemat_multiply(&tmp, &JJ, myres));
  igraph_sparsemat_destroy(&tmp);
  igraph_sparsemat_destroy(&JJ);
  IGRAPH_FINALLY_CLEAN(2);

  if (constres) {
    if (myres->cs->p [1] != 0) {
      *constres = myres->cs->x [0];
    } else {
      *constres = 0.0;
    }
  }

  if (!res) {
    igraph_sparsemat_destroy(myres);
  }

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
  igraph_bool_t directed=igraph_is_directed(graph);
  long int nzmax= directed ? no_of_edges : no_of_edges*2;
  long int i;  
  
  IGRAPH_CHECK(igraph_sparsemat_init(res, no_of_nodes, no_of_nodes, 
				     nzmax));
  
  for (i=0; i<no_of_edges; i++) {
    long int from=IGRAPH_FROM(graph, i);
    long int to=IGRAPH_TO(graph, i);
    IGRAPH_CHECK(igraph_sparsemat_entry(res, from, to, 1.0));
    if (!directed && from != to) {
      IGRAPH_CHECK(igraph_sparsemat_entry(res, to, from, 1.0));
    }
  }    
  
  return 0;
}

#define CHECK(x) if ((x)<0) { IGRAPH_ERROR("Cannot write to file", IGRAPH_EFILE); }

int igraph_sparsemat_print(const igraph_sparsemat_t *A,
			   FILE *outstream) {
  
  if (A->cs->nz < 0) {
    /* CC */
    int j, p;
    for (j=0; j<A->cs->n; j++) {
      CHECK(fprintf(outstream, "col %i: locations %i to %i\n", 
		    j, A->cs->p[j], A->cs->p[j+1]-1));
      for (p=A->cs->p[j]; p < A->cs->p[j+1]; p++) {
	CHECK(fprintf(outstream, "%i : %g\n", A->cs->i[p], A->cs->x[p]));
      }
    }
  } else {
    /* Triplet */
    int p;
    for (p=0; p<A->cs->nz; p++) {
      CHECK(fprintf(outstream, "%i %i : %g\n", 
		    A->cs->i[p], A->cs->p[p], A->cs->x[p]));
    }
  }
  
  return 0;
}

#undef CHECK

int igraph_i_sparsemat_eye_triplet(igraph_sparsemat_t *A, int n, int nzmax, 
				   igraph_real_t value) {
  long int i;

  IGRAPH_CHECK(igraph_sparsemat_init(A, n, n, nzmax));
  
  for (i=0; i<n; i++) {
    igraph_sparsemat_entry(A, i, i, value);
  }
  
  return 0;
}

int igraph_i_sparsemat_eye_cc(igraph_sparsemat_t *A, int n,
			      igraph_real_t value) {
  long int i;

  if (! (A->cs = cs_spalloc(n, n, n, /*values=*/ 1, /*triplet=*/ 0)) ) {
    IGRAPH_ERROR("Cannot create eye sparse matrix", IGRAPH_FAILURE);
  }
  
  for (i=0; i<n; i++) {
    A->cs->p [i] = i;
    A->cs->i [i] = i;
    A->cs->x [i] = value;
  }
  A->cs->p [n] = n;
  
  return 0;
}

int igraph_sparsemat_eye(igraph_sparsemat_t *A, int n, int nzmax,
			 igraph_real_t value,
			 igraph_bool_t compress) {
  if (compress) {
    return(igraph_i_sparsemat_eye_cc(A, n, value));
  } else {
    return(igraph_i_sparsemat_eye_triplet(A, n, nzmax, value));
  }
}

int igraph_i_sparsemat_diag_triplet(igraph_sparsemat_t *A, int nzmax,
				    const igraph_vector_t *values) {

  long int i, n=igraph_vector_size(values);

  IGRAPH_CHECK(igraph_sparsemat_init(A, n, n, nzmax));
  
  for (i=0; i<n; i++) {
    igraph_sparsemat_entry(A, i, i, VECTOR(*values)[i]);
  }
  
  return 0;
  
}

int igraph_i_sparsemat_diag_cc(igraph_sparsemat_t *A,
			       const igraph_vector_t *values) {

  long int i, n=igraph_vector_size(values);

  if (! (A->cs = cs_spalloc(n, n, n, /*values=*/ 1, /*triplet=*/ 0)) ) {
    IGRAPH_ERROR("Cannot create eye sparse matrix", IGRAPH_FAILURE);
  }
  
  for (i=0; i<n; i++) {
    A->cs->p [i] = i;
    A->cs->i [i] = i;
    A->cs->x [i] = VECTOR(*values)[i];
  }
  A->cs->p [n] = n;
  
  return 0;

}

int igraph_sparsemat_diag(igraph_sparsemat_t *A, int nzmax,
			  const igraph_vector_t *values,
			  igraph_bool_t compress) {

  if (compress) {
    return(igraph_i_sparsemat_diag_cc(A, values));
  } else {
    return(igraph_i_sparsemat_diag_triplet(A, nzmax, values));
  }
}
