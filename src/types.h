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
 * Double ended queue data type.
 * \ingroup internal
 */

typedef struct s_dqueue {
  real_t *begin;
  real_t *end;
  real_t *stor_begin;
  real_t *stor_end;
} igraph_dqueue_t;

#define IGRAPH_DQUEUE_NULL { 0,0,0,0 }
#define IGRAPH_DQUEUE_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(igraph_dqueue_init(v, size)); \
  IGRAPH_FINALLY(igraph_dqueue_destroy, v); } while (0)

int igraph_dqueue_init    (igraph_dqueue_t* q, long int size);
void igraph_dqueue_destroy (igraph_dqueue_t* q);
bool_t igraph_dqueue_empty   (igraph_dqueue_t* q);
void igraph_dqueue_clear   (igraph_dqueue_t* q);
bool_t igraph_dqueue_full    (igraph_dqueue_t* q);
long int igraph_dqueue_size    (igraph_dqueue_t* q);
real_t igraph_dqueue_pop     (igraph_dqueue_t* q);
real_t igraph_dqueue_pop_back(igraph_dqueue_t* q);
real_t igraph_dqueue_head    (igraph_dqueue_t* q);
real_t igraph_dqueue_back    (igraph_dqueue_t* q);
int igraph_dqueue_push    (igraph_dqueue_t* q, real_t elem);

/* -------------------------------------------------- */
/* Flexible vector                                    */
/* -------------------------------------------------- */

/** 
 * Vector, dealing with arrays efficiently.
 * \ingroup types
 */

typedef struct s_vector {
  real_t* stor_begin;
  real_t* stor_end;
  real_t* end;
} igraph_vector_t;

#define IGRAPH_VECTOR_NULL { 0,0,0 }
#define IGRAPH_VECTOR_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(igraph_vector_init(v, size)); \
  IGRAPH_FINALLY(igraph_vector_destroy, v); } while (0)

/**
 * \ingroup vector
 * \define VECTOR
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
 * \param v The vector object.
 * 
 * Time complexity: O(1).
 */
#define VECTOR(v) ((v).stor_begin) /* DIRTY */
int igraph_vector_init      (igraph_vector_t* v, long int size);
int igraph_vector_init_copy (igraph_vector_t* v, real_t* data, long int length);
int igraph_vector_init_seq(igraph_vector_t *v, real_t from, real_t to);
int igraph_vector_init_real(igraph_vector_t *v, int no, ...);
int igraph_vector_init_int(igraph_vector_t *v, int no, ...);
int igraph_vector_init_real_end(igraph_vector_t *v, real_t endmark, ...);
int igraph_vector_init_int_end(igraph_vector_t *v, int endmark, ...);
const igraph_vector_t *igraph_vector_view (const igraph_vector_t *v, const real_t *data, 
			     long int length);
void igraph_vector_destroy   (igraph_vector_t* v);
int igraph_vector_reserve   (igraph_vector_t* v, long int size);
bool_t igraph_vector_empty     (const igraph_vector_t* v);
long int igraph_vector_size      (const igraph_vector_t* v);
void igraph_vector_clear     (igraph_vector_t* v);
void igraph_vector_null      (igraph_vector_t* v);
int igraph_vector_push_back (igraph_vector_t* v, real_t e);
real_t igraph_vector_e         (const igraph_vector_t* v, long int pos);
real_t*igraph_vector_e_ptr  (const igraph_vector_t* v, long int pos);
void igraph_vector_set       (igraph_vector_t* v, long int pos, real_t value);
real_t igraph_vector_tail(const igraph_vector_t *v);
real_t igraph_vector_pop_back(igraph_vector_t* v);
int igraph_vector_order(const igraph_vector_t* v, igraph_vector_t* res, integer_t maxval);
void igraph_vector_sort(igraph_vector_t *v);
int igraph_vector_resize(igraph_vector_t* v, long int newsize);
real_t igraph_vector_max(const igraph_vector_t* v);
long int igraph_vector_which_max(const igraph_vector_t* v);
void igraph_vector_copy_to(const igraph_vector_t *v, real_t* to);
int igraph_vector_copy(igraph_vector_t *to, const igraph_vector_t *from);
real_t igraph_vector_sum(const igraph_vector_t *v);
real_t igraph_vector_prod(const igraph_vector_t *v);
void igraph_vector_remove_section(igraph_vector_t *v, long int from, long int to);
int igraph_vector_move_interval(igraph_vector_t *v, long int begin, long int end, 
			 long int to);
void igraph_vector_remove(igraph_vector_t *v, long int elem);
void igraph_vector_permdelete(igraph_vector_t *v, long int *index, long int nremove);
void igraph_vector_remove_negidx(igraph_vector_t *v, const igraph_vector_t *neg, long int nremove);
bool_t igraph_vector_isininterval(const igraph_vector_t *v, real_t low, real_t high);
bool_t igraph_vector_any_smaller(const igraph_vector_t *v, real_t limit);
bool_t igraph_vector_is_equal(const igraph_vector_t *lhs, const igraph_vector_t *rhs);
bool_t igraph_vector_binsearch(const igraph_vector_t *v, real_t what, long int *pos);

/* -------------------------------------------------- */
/* Flexible vector, storing pointers                  */
/* -------------------------------------------------- */

/** 
 * Vector, storing pointers efficiently
 * \ingroup internal
 * 
 */
typedef struct s_vector_ptr {
  void** stor_begin;
  void** stor_end;
  void** end;
} igraph_vector_ptr_t;

#define IGRAPH_VECTOR_PTR_NULL { 0,0,0 }
#define IGRAPH_VECTOR_PTR_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(igraph_vector_ptr_init(v, size)); \
  IGRAPH_FINALLY(igraph_vector_ptr_destroy, v); } while (0)

int igraph_vector_ptr_init      (igraph_vector_ptr_t* v, long int size);
int igraph_vector_ptr_init_copy (igraph_vector_ptr_t* v, void** data, long int length);
const igraph_vector_ptr_t *igraph_vector_ptr_view (const igraph_vector_ptr_t *v, 
				     void *const *data, long int length);
void igraph_vector_ptr_destroy   (igraph_vector_ptr_t* v);
void igraph_vector_ptr_free_all   (igraph_vector_ptr_t* v);
void igraph_vector_ptr_destroy_all   (igraph_vector_ptr_t* v);
int igraph_vector_ptr_reserve   (igraph_vector_ptr_t* v, long int size);
bool_t igraph_vector_ptr_empty     (const igraph_vector_ptr_t* v);
long int igraph_vector_ptr_size      (const igraph_vector_ptr_t* v);
void igraph_vector_ptr_clear     (igraph_vector_ptr_t* v);
void igraph_vector_ptr_null      (igraph_vector_ptr_t* v);
int igraph_vector_ptr_push_back (igraph_vector_ptr_t* v, void* e);
void* igraph_vector_ptr_e         (const igraph_vector_ptr_t* v, long int pos);
void igraph_vector_ptr_set       (igraph_vector_ptr_t* v, long int pos, void* value);
int igraph_vector_ptr_resize(igraph_vector_ptr_t* v, long int newsize);
void igraph_vector_ptr_copy_to(const igraph_vector_ptr_t *v, void** to);
int igraph_vector_ptr_copy(igraph_vector_ptr_t *to, const igraph_vector_ptr_t *from);
void igraph_vector_ptr_remove(igraph_vector_ptr_t *v, long int pos);

/* -------------------------------------------------- */
/* Matrix, very similar to vector                     */
/* -------------------------------------------------- */

/** 
 * \section about_igraph_matrix_t_objects About \type igraph_matrix_t objects
 * 
 * <para>This type is just an interface to vector.</para>
 *
 * <para>The \type igraph_matrix_t type ususally stores n
 * elements in O(n) space, but not always, see the documentation of
 * the vector type.</para>
 */
typedef struct s_matrix {
  igraph_vector_t data;
  long int nrow, ncol;
} igraph_matrix_t;

#define IGRAPH_MATRIX_NULL { IGRAPH_VECTOR_NULL, 0, 0 }
#define IGRAPH_MATRIX_INIT_FINALLY(m, nr, nc) \
  do { IGRAPH_CHECK(igraph_matrix_init(m, nr, nc)); \
  IGRAPH_FINALLY(igraph_matrix_destroy, m); } while (0)

/**
 * \ingroup matrix
 * \define MATRIX
 * \brief Accessing an element of a matrix.
 *
 * Note that there are no range checks right now. 
 * This functionality might be redefines as a proper function later. 
 * \param m The matrix object.
 * \param i The index of the row, starting with zero.
 * \param j The index of the column, starting with zero.
 *
 * Time complexity: O(1).
 */
#define MATRIX(m,i,j) ((m).data.stor_begin[(m).nrow*(j)+(i)])
int igraph_matrix_init(igraph_matrix_t *m, long int nrow, long int ncol);
void igraph_matrix_destroy(igraph_matrix_t *m);
int igraph_matrix_resize(igraph_matrix_t *m, long int nrow, long int ncol);
long int igraph_matrix_size(const igraph_matrix_t *m);
long int igraph_matrix_nrow(const igraph_matrix_t *m);
long int igraph_matrix_ncol(const igraph_matrix_t *m);
int igraph_matrix_copy_to(const igraph_matrix_t *m, real_t *to);
int igraph_matrix_null(igraph_matrix_t *m);
int igraph_matrix_add_cols(igraph_matrix_t *m, long int n);
int igraph_matrix_add_rows(igraph_matrix_t *m, long int n);
int igraph_matrix_remove_col(igraph_matrix_t *m, long int col);
int igraph_matrix_permdelete_rows(igraph_matrix_t *m, long int *index, long int nremove);
int igraph_matrix_delete_rows_neg(igraph_matrix_t *m, igraph_vector_t *neg, long int nremove);
int igraph_matrix_copy(igraph_matrix_t *to, const igraph_matrix_t *from);

/* -------------------------------------------------- */
/* Plain stack                                        */
/* -------------------------------------------------- */

/**
 * Stack data type.
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
 * Heap data type.
 * \ingroup internal
 */

typedef struct s_heap {
  real_t* stor_begin;
  real_t* stor_end;
  real_t* end;
  int destroy;
} igraph_heap_t;

#define IGRAPH_HEAP_NULL { 0,0,0 }

int igraph_heap_init           (igraph_heap_t* h, long int size);
int igraph_heap_init_array     (igraph_heap_t *t, real_t* data, long int len);
void igraph_heap_destroy        (igraph_heap_t* h);
bool_t igraph_heap_empty          (igraph_heap_t* h);
int igraph_heap_push           (igraph_heap_t* h, real_t elem);
real_t igraph_heap_max       (igraph_heap_t* h);
real_t igraph_heap_delete_max(igraph_heap_t* h);
long int igraph_heap_size      (igraph_heap_t* h);
int igraph_heap_reserve        (igraph_heap_t* h, long int size);

void igraph_heap_i_build(real_t* arr, long int size, long int head);
void igraph_heap_i_shift_up(real_t* arr, long int size, long int elem);
void igraph_heap_i_sink(real_t* arr, long int size, long int head);
void igraph_heap_i_switch(real_t* arr, long int e1, long int e2);

/* -------------------------------------------------- */
/* Indexed heap                                       */
/* -------------------------------------------------- */

/**
 * Indexed heap data type.
 * \ingroup internal
 */

typedef struct s_indheap {
  real_t* stor_begin;
  real_t* stor_end;
  real_t* end;
  int destroy;
  long int* index_begin;
} igraph_indheap_t;

#define IGRAPH_INDHEAP_NULL { 0,0,0,0,0 }

int igraph_indheap_init           (igraph_indheap_t* h, long int size);
int igraph_indheap_init_array     (igraph_indheap_t *t, real_t* data, long int len);
void igraph_indheap_destroy        (igraph_indheap_t* h);
bool_t igraph_indheap_empty          (igraph_indheap_t* h);
int igraph_indheap_push           (igraph_indheap_t* h, real_t elem);
real_t igraph_indheap_max       (igraph_indheap_t* h);
real_t igraph_indheap_delete_max(igraph_indheap_t* h);
long int igraph_indheap_size      (igraph_indheap_t* h);
int igraph_indheap_reserve        (igraph_indheap_t* h, long int size);
long int igraph_indheap_max_index(igraph_indheap_t *h);

void igraph_indheap_i_build(igraph_indheap_t* h, long int head);
void igraph_indheap_i_shift_up(igraph_indheap_t* h, long int elem);
void igraph_indheap_i_sink(igraph_indheap_t* h, long int head);
void igraph_indheap_i_switch(igraph_indheap_t* h, long int e1, long int e2);

/* -------------------------------------------------- */
/* Doubly indexed heap                                */
/* -------------------------------------------------- */

/* This is a heap containing double elements and 
   two indices, its intended usage is the storage of
   weighted edges.
*/

/**
 * Doubly indexed heap data type.
 * \ingroup internal
 */

typedef struct s_indheap_d {
  real_t* stor_begin;
  real_t* stor_end;
  real_t* end;
  int destroy;
  long int* index_begin;
  long int* index2_begin;
} igraph_d_indheap_t;


#define IGRAPH_D_INDHEAP_NULL { 0,0,0,0,0,0 }

int igraph_d_indheap_init           (igraph_d_indheap_t* h, long int size);
void igraph_d_indheap_destroy        (igraph_d_indheap_t* h);
bool_t igraph_d_indheap_empty          (igraph_d_indheap_t* h);
int igraph_d_indheap_push           (igraph_d_indheap_t* h, real_t elem, 
			      long int idx, long int idx2);
real_t igraph_d_indheap_max       (igraph_d_indheap_t* h);
real_t igraph_d_indheap_delete_max(igraph_d_indheap_t* h);
long int igraph_d_indheap_size      (igraph_d_indheap_t* h);
int igraph_d_indheap_reserve        (igraph_d_indheap_t* h, long int size);
void igraph_d_indheap_max_index(igraph_d_indheap_t *h, long int *idx, long int *idx2);

void igraph_d_indheap_i_build(igraph_d_indheap_t* h, long int head);
void igraph_d_indheap_i_shift_up(igraph_d_indheap_t* h, long int elem);
void igraph_d_indheap_i_sink(igraph_d_indheap_t* h, long int head);
void igraph_d_indheap_i_switch(igraph_d_indheap_t* h, long int e1, long int e2);

/**
 * Vector of strings
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
void igraph_strvector_remove_negidx(igraph_strvector_t *v, const igraph_vector_t *neg,
				    long int nremove);
  
/**
 * Trie data type
 * \ingroup internal
 */

typedef struct s_igraph_trie_node {
  igraph_strvector_t strs;
  igraph_vector_ptr_t children;
  igraph_vector_t values;
} igraph_trie_node_t;

typedef struct s_igraph_trie {
  igraph_strvector_t strs;
  igraph_vector_ptr_t children;
  igraph_vector_t values;
  long int maxvalue;
  bool_t storekeys;
  igraph_strvector_t keys;
} igraph_trie_t;

#define IGRAPH_TRIE_NULL { IGRAPH_STRVECTOR_NULL, IGRAPH_VECTOR_PTR_NULL, \
                           IGRAPH_VECTOR_NULL, 0, 0, IGRAPH_STRVECTOR_NULL }
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
