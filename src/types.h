/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#ifndef REST_TYPES_H
#define REST_TYPES_H

#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif

#include "error.h"

typedef double integer_t;
typedef double real_t;
typedef int    bool_t;

/* -------------------------------------------------- */
/* double ended queue, very useful                    */
/* -------------------------------------------------- */

/**
 * \defgroup dqueue Double ended queue data type.
 * \ingroup internal
 */

typedef struct s_dqueue {
  real_t *begin;
  real_t *end;
  real_t *stor_begin;
  real_t *stor_end;
} dqueue_t;

#define DQUEUE_NULL { 0,0,0,0 }
#define DQUEUE_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(dqueue_init(v, size)); \
  IGRAPH_FINALLY(dqueue_destroy, v); } while (0)

int dqueue_init    (dqueue_t* q, long int size);
void dqueue_destroy (dqueue_t* q);
bool_t dqueue_empty   (dqueue_t* q);
void dqueue_clear   (dqueue_t* q);
bool_t dqueue_full    (dqueue_t* q);
long int dqueue_size    (dqueue_t* q);
real_t dqueue_pop     (dqueue_t* q);
real_t dqueue_pop_back(dqueue_t* q);
real_t dqueue_head    (dqueue_t* q);
real_t dqueue_back    (dqueue_t* q);
int dqueue_push    (dqueue_t* q, real_t elem);

/* -------------------------------------------------- */
/* Flexible vector                                    */
/* -------------------------------------------------- */

/** \defgroup vector Vector, dealing with arrays efficiently.
 * \ingroup types
 * 
 * The <code>vector_t</code> data type is a simple and efficient
 * interface to arrays containing real numbers. It is something
 * similar as (but much simpler than) the <code>vector</code> template
 * in the C++ standard library.
 *
 * Vectors are used extensively in \a igraph, all functions which
 * expects or returns a list of numbers use <code>vector_t</code> to
 * achive this.
 * 
 * The <code>vector_t</code> type susally uses <code>O(n)</code> space
 * to store <code>n</code> elements. Sometimes it uses more, this is
 * because vectors can shrink, but even if they shrink, the current
 * implementation does not free a single bit of memory.
 */

/**
 * vector_t:
 * 
 * Flecible arrays
 */
typedef struct s_vector {
   /*< private >*/
  real_t* stor_begin;
  real_t* stor_end;
  real_t* end;
} vector_t;

#define VECTOR_NULL { 0,0,0 }
#define VECTOR_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(vector_init(v, size)); \
  IGRAPH_FINALLY(vector_destroy, v); } while (0)

/**
 * \ingroup vector
 * \brief Accessing an element of a vector
 * 
 * Usage: 
 * \verbatim VECTOR(v)[0] \endverbatim 
 * to access the first element of the vector, you can also use this in
 * assignments, like: 
 * \verbatim VECTOR(v)[10]=5; \endverbatim
 *
 * Note that there are no range checks right now.
 * This functionality might be redefined later as a real function
 * instead of a <code>#define</code>. 
 * @param v The vector object.
 * 
 * Time complexity: <code>O(1)</code>.
 */
#define VECTOR(v) ((v).stor_begin) /* DIRTY */
int vector_init      (vector_t* v, long int size);
int vector_init_copy (vector_t* v, real_t* data, long int length);
int vector_init_seq(vector_t *v, real_t from, real_t to);
int vector_init_real(vector_t *v, int no, ...);
int vector_init_int(vector_t *v, int no, ...);
int vector_init_real_end(vector_t *v, real_t endmark, ...);
int vector_init_int_end(vector_t *v, int endmark, ...);
const vector_t *vector_view (const vector_t *v, const real_t *data, 
			     long int length);
void vector_destroy   (vector_t* v);
int vector_reserve   (vector_t* v, long int size);
bool_t vector_empty     (const vector_t* v);
long int vector_size      (const vector_t* v);
void vector_clear     (vector_t* v);
void vector_null      (vector_t* v);
int vector_push_back (vector_t* v, real_t e);
real_t vector_e         (const vector_t* v, long int pos);
real_t*vector_e_ptr  (const vector_t* v, long int pos);
void vector_set       (vector_t* v, long int pos, real_t value);
real_t vector_tail(const vector_t *v);
real_t vector_pop_back(vector_t* v);
int vector_order(const vector_t* v, vector_t* res, integer_t maxval);
void vector_sort(vector_t *v);
int vector_resize(vector_t* v, long int newsize);
real_t vector_max(const vector_t* v);
void vector_copy_to(const vector_t *v, real_t* to);
int vector_copy(vector_t *to, const vector_t *from);
real_t vector_sum(const vector_t *v);
real_t vector_prod(const vector_t *v);
void vector_remove_section(vector_t *v, long int from, long int to);
int vector_move_interval(vector_t *v, long int begin, long int end, 
			 long int to);
void vector_remove(vector_t *v, long int elem);
void vector_permdelete(vector_t *v, long int *index, long int nremove);
void vector_remove_negidx(vector_t *v, const vector_t *neg, long int nremove);
bool_t vector_isininterval(const vector_t *v, real_t low, real_t high);
bool_t vector_any_smaller(const vector_t *v, real_t limit);
bool_t vector_is_equal(const vector_t *lhs, const vector_t *rhs);
bool_t vector_binsearch(const vector_t *v, real_t what, long int *pos);

/* -------------------------------------------------- */
/* Flexible vector, storing pointers                  */
/* -------------------------------------------------- */

/** \defgroup vectorptr Vector, storing pointers efficiently
 * \ingroup internal
 * 
 */
typedef struct s_vector_ptr {
  void** stor_begin;
  void** stor_end;
  void** end;
} vector_ptr_t;

#define VECTOR_PTR_NULL { 0,0,0 }
#define VECTOR_PTR_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(vector_ptr_init(v, size)); \
  IGRAPH_FINALLY(vector_ptr_destroy, v); } while (0)

int vector_ptr_init      (vector_ptr_t* v, long int size);
int vector_ptr_init_copy (vector_ptr_t* v, void** data, long int length);
const vector_ptr_t *vector_ptr_view (const vector_ptr_t *v, 
				     void *const *data, long int length);
void vector_ptr_destroy   (vector_ptr_t* v);
void vector_ptr_free_all   (vector_ptr_t* v);
void vector_ptr_destroy_all   (vector_ptr_t* v);
int vector_ptr_reserve   (vector_ptr_t* v, long int size);
bool_t vector_ptr_empty     (const vector_ptr_t* v);
long int vector_ptr_size      (const vector_ptr_t* v);
void vector_ptr_clear     (vector_ptr_t* v);
void vector_ptr_null      (vector_ptr_t* v);
int vector_ptr_push_back (vector_ptr_t* v, void* e);
void* vector_ptr_e         (const vector_ptr_t* v, long int pos);
void vector_ptr_set       (vector_ptr_t* v, long int pos, void* value);
int vector_ptr_resize(vector_ptr_t* v, long int newsize);
void vector_ptr_copy_to(const vector_ptr_t *v, void** to);
int vector_ptr_copy(vector_ptr_t *to, const vector_ptr_t *from);
void vector_ptr_remove(vector_ptr_t *v, long int pos);

/* -------------------------------------------------- */
/* Matrix, very similar to vector                     */
/* -------------------------------------------------- */

/** \defgroup matrix Matrix type, for storing real matrices efficienly. 
 * \ingroup types
 * 
 * This type is just an interface to vector.
 *
 * The <code>matrix_t</code> type ususally stores <code>n</code>
 * elements in <code>O(n)</code> space, but not always, see the
 * documentation of the vector type.
 */
typedef struct s_matrix {
  vector_t data;
  long int nrow, ncol;
} matrix_t;

#define MATRIX_NULL { VECTOR_NULL, 0, 0 }
#define MATRIX_INIT_FINALLY(m, nr, nc) \
  do { IGRAPH_CHECK(matrix_init(m, nr, nc)); \
  IGRAPH_FINALLY(matrix_destroy, m); } while (0)

/**
 * \ingroup matrix
 * \brief Accessing an element of a matrix.
 *
 * Note that there are no range checks right now. 
 * This functionality might be redefines as a proper function later. 
 * @param m The matrix object.
 * @param i The index of the row, starting with zero.
 * @param j The index of the column, starting with zero.
 *
 * Time complexity: <code>O(1)</code>.
 */
#define MATRIX(m,i,j) ((m).data.stor_begin[(m).nrow*(j)+(i)])
int matrix_init(matrix_t *m, long int nrow, long int ncol);
void matrix_destroy(matrix_t *m);
int matrix_resize(matrix_t *m, long int nrow, long int ncol);
long int matrix_size(const matrix_t *m);
long int matrix_nrow(const matrix_t *m);
long int matrix_ncol(const matrix_t *m);
int matrix_copy_to(const matrix_t *m, real_t *to);
int matrix_null(matrix_t *m);
int matrix_add_cols(matrix_t *m, long int n);
int matrix_add_rows(matrix_t *m, long int n);
int matrix_remove_col(matrix_t *m, long int col);
int matrix_permdelete_rows(matrix_t *m, long int *index, long int nremove);
int matrix_delete_rows_neg(matrix_t *m, vector_t *neg, long int nremove);
int matrix_copy(matrix_t *to, const matrix_t *from);

/* -------------------------------------------------- */
/* Plain stack                                        */
/* -------------------------------------------------- */

/**
 * \defgroup stack Stack data type.
 * \ingroup internal
 */

typedef struct s_stack {
  real_t* stor_begin;
  real_t* stor_end;
  real_t* end;
} igraph_stack_t;

#define IGRAPH_STACK_NULL { 0,0,0 }

int igraph_stack_init       (igraph_stack_t* s, long int size);
void igraph_stack_destroy    (igraph_stack_t* s);
int igraph_stack_reserve    (igraph_stack_t* s, long int size);
bool_t igraph_stack_empty      (igraph_stack_t* s);
long int igraph_stack_size       (igraph_stack_t* s);
void igraph_stack_clear      (igraph_stack_t* s);
int igraph_stack_push       (igraph_stack_t* s, real_t elem);
real_t igraph_stack_pop        (igraph_stack_t* s);

/* -------------------------------------------------- */
/* Heap                                               */
/* -------------------------------------------------- */

/**
 * \defgroup heap Heap data type.
 * \ingroup internal
 */

typedef struct s_heap {
  real_t* stor_begin;
  real_t* stor_end;
  real_t* end;
  int destroy;
} heap_t;

#define HEAP_NULL { 0,0,0 }

int heap_init           (heap_t* h, long int size);
int heap_init_array     (heap_t *t, real_t* data, long int len);
void heap_destroy        (heap_t* h);
bool_t heap_empty          (heap_t* h);
int heap_push           (heap_t* h, real_t elem);
real_t heap_max       (heap_t* h);
real_t heap_delete_max(heap_t* h);
long int heap_size      (heap_t* h);
int heap_reserve        (heap_t* h, long int size);

void heap_i_build(real_t* arr, long int size, long int head);
void heap_i_shift_up(real_t* arr, long int size, long int elem);
void heap_i_sink(real_t* arr, long int size, long int head);
void heap_i_switch(real_t* arr, long int e1, long int e2);

/* -------------------------------------------------- */
/* Indexed heap                                       */
/* -------------------------------------------------- */

/**
 * \defgroup indheap Indexed heap data type.
 * \ingroup internal
 */

typedef struct s_indheap {
  real_t* stor_begin;
  real_t* stor_end;
  real_t* end;
  int destroy;
  long int* index_begin;
} indheap_t;

#define INDHEAP_NULL { 0,0,0,0,0 }

int indheap_init           (indheap_t* h, long int size);
int indheap_init_array     (indheap_t *t, real_t* data, long int len);
void indheap_destroy        (indheap_t* h);
bool_t indheap_empty          (indheap_t* h);
int indheap_push           (indheap_t* h, real_t elem);
real_t indheap_max       (indheap_t* h);
real_t indheap_delete_max(indheap_t* h);
long int indheap_size      (indheap_t* h);
int indheap_reserve        (indheap_t* h, long int size);
long int indheap_max_index(indheap_t *h);

void indheap_i_build(indheap_t* h, long int head);
void indheap_i_shift_up(indheap_t* h, long int elem);
void indheap_i_sink(indheap_t* h, long int head);
void indheap_i_switch(indheap_t* h, long int e1, long int e2);

/* -------------------------------------------------- */
/* Doubly indexed heap                                */
/* -------------------------------------------------- */

/* This is a heap containing double elements and 
   two indices, its intended usage is the storage of
   weighted edges.
*/

/**
 * \defgroup doubleindheap Doubly indexed heap data type.
 * \ingroup internal
 */

typedef struct s_indheap_d {
  real_t* stor_begin;
  real_t* stor_end;
  real_t* end;
  int destroy;
  long int* index_begin;
  long int* index2_begin;
} d_indheap_t;


#define D_INDHEAP_NULL { 0,0,0,0,0,0 }

int d_indheap_init           (d_indheap_t* h, long int size);
void d_indheap_destroy        (d_indheap_t* h);
bool_t d_indheap_empty          (d_indheap_t* h);
int d_indheap_push           (d_indheap_t* h, real_t elem, 
			      long int idx, long int idx2);
real_t d_indheap_max       (d_indheap_t* h);
real_t d_indheap_delete_max(d_indheap_t* h);
long int d_indheap_size      (d_indheap_t* h);
int d_indheap_reserve        (d_indheap_t* h, long int size);
void d_indheap_max_index(d_indheap_t *h, long int *idx, long int *idx2);

void d_indheap_i_build(d_indheap_t* h, long int head);
void d_indheap_i_shift_up(d_indheap_t* h, long int elem);
void d_indheap_i_sink(d_indheap_t* h, long int head);
void d_indheap_i_switch(d_indheap_t* h, long int e1, long int e2);

/**
 * \defgroup strvector Vector of strings
 * \ingroup internal
 */

typedef struct s_igraph_strvector {
  char **data;
  long int len;
} igraph_strvector_t;

#define IGRAPH_STRVECTOR_NULL { 0,0 }
#define IGRAPH_STRVECTOR_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(igraph_strvector_init(v, size)); \
  IGRAPH_FINALLY( (igraph_finally_func_t*) igraph_strvector_destroy, v); } while (0)

int igraph_strvector_init(igraph_strvector_t *sv, long int len);
void igraph_strvector_destroy(igraph_strvector_t *sv);
long int igraph_strvector_size(const igraph_strvector_t *sv);
void igraph_strvector_get(const igraph_strvector_t *sv, 
			  long int idx, char **value);
int igraph_strvector_set(igraph_strvector_t *sv, long int idx, 
			 const char *value);
void igraph_strvector_remove_section(igraph_strvector_t *v, long int from, 
				     long int to);
void igraph_strvector_remove(igraph_strvector_t *v, long int elem);
void igraph_strvector_move_interval(igraph_strvector_t *v, long int begin, 
				   long int end, long int to);
int igraph_strvector_copy(igraph_strvector_t *to, 
			  const igraph_strvector_t *from);
int igraph_strvector_resize(igraph_strvector_t* v, long int newsize);
int igraph_strvector_add(igraph_strvector_t *v, const char *value);
void igraph_strvector_permdelete(igraph_strvector_t *v, long int *index, 
				 long int nremove);
void igraph_strvector_remove_negidx(igraph_strvector_t *v, const vector_t *neg,
				    long int nremove);
  
/**
 * \defgroup igraphtrie Trie data type
 * \ingroup internal
 */

typedef struct s_igraph_trie_node {
  igraph_strvector_t strs;
  vector_ptr_t children;
  vector_t values;
} igraph_trie_node_t;

typedef struct s_igraph_trie {
  igraph_strvector_t strs;
  vector_ptr_t children;
  vector_t values;
  long int maxvalue;
  bool_t storekeys;
  igraph_strvector_t keys;
} igraph_trie_t;

#define IGRAPH_TRIE_NULL { IGRAPH_STRVECTOR_NULL, VECTOR_PTR_NULL, \
                           VECTOR_NULL, 0, 0, IGRAPH_STRVECTOR_NULL }
#define IGRAPH_TRIE_INIT_FINALLY(tr, sk) \
  do { IGRAPH_CHECK(igraph_trie_init(tr, sk)); \
  IGRAPH_FINALLY(igraph_trie_destroy, tr); } while (0)

int igraph_trie_init(igraph_trie_t *t, bool_t storekeys);
void igraph_trie_destroy(igraph_trie_t *t);
int igraph_trie_get(igraph_trie_t *t, const char *key, long int *id);
int igraph_trie_get2(igraph_trie_t *t, const char *key, long int length, 
		     long int *id);
void igraph_trie_idx(igraph_trie_t *t, long int idx, char **str);
long int igraph_trie_size(igraph_trie_t *t);

#endif
