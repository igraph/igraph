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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#ifndef REST_TYPES_H
#define REST_TYPES_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif

#include "error.h"

typedef double igraph_integer_t;
typedef double igraph_real_t;
typedef int    igraph_bool_t;

/* -------------------------------------------------- */
/* double ended queue, very useful                    */
/* -------------------------------------------------- */

/**
 * Double ended queue data type.
 * \ingroup internal
 */

typedef struct s_dqueue {
  igraph_real_t *begin;
  igraph_real_t *end;
  igraph_real_t *stor_begin;
  igraph_real_t *stor_end;
} igraph_dqueue_t;

#define IGRAPH_DQUEUE_NULL { 0,0,0,0 }
#define IGRAPH_DQUEUE_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(igraph_dqueue_init(v, size)); \
  IGRAPH_FINALLY(igraph_dqueue_destroy, v); } while (0)

int igraph_dqueue_init    (igraph_dqueue_t* q, long int size);
void igraph_dqueue_destroy (igraph_dqueue_t* q);
igraph_bool_t igraph_dqueue_empty   (igraph_dqueue_t* q);
void igraph_dqueue_clear   (igraph_dqueue_t* q);
igraph_bool_t igraph_dqueue_full    (igraph_dqueue_t* q);
long int igraph_dqueue_size    (igraph_dqueue_t* q);
igraph_real_t igraph_dqueue_pop     (igraph_dqueue_t* q);
igraph_real_t igraph_dqueue_pop_back(igraph_dqueue_t* q);
igraph_real_t igraph_dqueue_head    (igraph_dqueue_t* q);
igraph_real_t igraph_dqueue_back    (igraph_dqueue_t* q);
int igraph_dqueue_push    (igraph_dqueue_t* q, igraph_real_t elem);

/* -------------------------------------------------- */
/* Flexible vector                                    */
/* -------------------------------------------------- */

/** 
 * Vector, dealing with arrays efficiently.
 * \ingroup types
 */

typedef struct s_vector {
  igraph_real_t* stor_begin;
  igraph_real_t* stor_end;
  igraph_real_t* end;
} igraph_vector_t;

#define IGRAPH_VECTOR_NULL { 0,0,0 }
#define IGRAPH_VECTOR_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(igraph_vector_init(v, size)); \
  IGRAPH_FINALLY(igraph_vector_destroy, v); } while (0)

/**
 * \ingroup vector
 * \define VECTOR
 * \brief Accessing an element of a vector.
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
#define VECTOR(v) ((v).stor_begin) 
int igraph_vector_init      (igraph_vector_t* v, long int size);
int igraph_vector_init_copy (igraph_vector_t* v, igraph_real_t* data, long int length);
int igraph_vector_init_seq(igraph_vector_t *v, igraph_real_t from, igraph_real_t to);
int igraph_vector_init_real(igraph_vector_t *v, int no, ...);
int igraph_vector_init_int(igraph_vector_t *v, int no, ...);
int igraph_vector_init_real_end(igraph_vector_t *v, igraph_real_t endmark, ...);
int igraph_vector_init_int_end(igraph_vector_t *v, int endmark, ...);
const igraph_vector_t *igraph_vector_view (const igraph_vector_t *v, const igraph_real_t *data, 
			     long int length);
void igraph_vector_destroy   (igraph_vector_t* v);
int igraph_vector_reserve   (igraph_vector_t* v, long int size);
igraph_bool_t igraph_vector_empty     (const igraph_vector_t* v);
long int igraph_vector_size      (const igraph_vector_t* v);
void igraph_vector_clear     (igraph_vector_t* v);
void igraph_vector_null      (igraph_vector_t* v);
int igraph_vector_push_back (igraph_vector_t* v, igraph_real_t e);
igraph_real_t igraph_vector_e         (const igraph_vector_t* v, long int pos);
igraph_real_t*igraph_vector_e_ptr  (const igraph_vector_t* v, long int pos);
void igraph_vector_set       (igraph_vector_t* v, long int pos, igraph_real_t value);
igraph_real_t igraph_vector_tail(const igraph_vector_t *v);
igraph_real_t igraph_vector_pop_back(igraph_vector_t* v);
int igraph_vector_order(const igraph_vector_t* v, const igraph_vector_t *v2,
			igraph_vector_t* res, igraph_integer_t maxval);
int igraph_vector_order1(const igraph_vector_t* v, 
			igraph_vector_t* res, igraph_integer_t maxval);
int igraph_vector_order2(igraph_vector_t *v);
void igraph_vector_sort(igraph_vector_t *v);
int igraph_vector_resize(igraph_vector_t* v, long int newsize);
igraph_real_t igraph_vector_max(const igraph_vector_t* v);
long int igraph_vector_which_max(const igraph_vector_t* v);
igraph_real_t igraph_vector_min(const igraph_vector_t* v);
long int igraph_vector_which_min(const igraph_vector_t* v);
void igraph_vector_copy_to(const igraph_vector_t *v, igraph_real_t* to);
int igraph_vector_copy(igraph_vector_t *to, const igraph_vector_t *from);
igraph_real_t igraph_vector_sum(const igraph_vector_t *v);
igraph_real_t igraph_vector_prod(const igraph_vector_t *v);
void igraph_vector_remove_section(igraph_vector_t *v, long int from, long int to);
int igraph_vector_move_interval(igraph_vector_t *v, long int begin, long int end, 
			 long int to);
void igraph_vector_remove(igraph_vector_t *v, long int elem);
void igraph_vector_permdelete(igraph_vector_t *v, long int *index, long int nremove);
void igraph_vector_remove_negidx(igraph_vector_t *v, const igraph_vector_t *neg, long int nremove);
igraph_bool_t igraph_vector_isininterval(const igraph_vector_t *v, igraph_real_t low, igraph_real_t high);
igraph_bool_t igraph_vector_any_smaller(const igraph_vector_t *v, igraph_real_t limit);
igraph_bool_t igraph_vector_is_equal(const igraph_vector_t *lhs, const igraph_vector_t *rhs);
igraph_bool_t igraph_vector_binsearch(const igraph_vector_t *v, igraph_real_t what, long int *pos);
igraph_bool_t igraph_vector_binsearch2(const igraph_vector_t *v, igraph_real_t what);
void igraph_vector_multiply(igraph_vector_t *v, igraph_real_t by);
igraph_bool_t igraph_vector_contains(const igraph_vector_t *v, igraph_real_t e);
igraph_bool_t igraph_vector_search(const igraph_vector_t *v, long int from, igraph_real_t what, 
			    long int *pos);
int igraph_vector_filter_smaller(igraph_vector_t *v, igraph_real_t elem);
int igraph_vector_append(igraph_vector_t *to, const igraph_vector_t *from);
int igraph_vector_get_interval(const igraph_vector_t *v, igraph_vector_t *res,
			       long int from, long int to);
int igraph_vector_rank(const igraph_vector_t *v, igraph_vector_t *res, 
		       long int nodes);

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
igraph_bool_t igraph_vector_ptr_empty     (const igraph_vector_ptr_t* v);
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
 * <para>This type is just an interface to \type igraph_vector_t.</para>
 *
 * <para>The \type igraph_matrix_t type usually stores n
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
int igraph_matrix_copy_to(const igraph_matrix_t *m, igraph_real_t *to);
int igraph_matrix_null(igraph_matrix_t *m);
int igraph_matrix_add_cols(igraph_matrix_t *m, long int n);
int igraph_matrix_add_rows(igraph_matrix_t *m, long int n);
int igraph_matrix_remove_col(igraph_matrix_t *m, long int col);
int igraph_matrix_permdelete_rows(igraph_matrix_t *m, long int *index, long int nremove);
int igraph_matrix_delete_rows_neg(igraph_matrix_t *m, igraph_vector_t *neg, long int nremove);
int igraph_matrix_copy(igraph_matrix_t *to, const igraph_matrix_t *from);
igraph_real_t igraph_matrix_max(const igraph_matrix_t *m);
void igraph_matrix_multiply(igraph_matrix_t *m, igraph_real_t by);
int igraph_matrix_select_rows(const igraph_matrix_t *m, igraph_matrix_t *res, 
			      const igraph_vector_t *rows);
int igraph_matrix_get_col(const igraph_matrix_t *m, igraph_vector_t *res,
			  long int index);
igraph_real_t igraph_matrix_sum(const igraph_matrix_t *m);

/* -------------------------------------------------- */
/* 3D array                                           */
/* -------------------------------------------------- */

typedef struct s_array3 {
  igraph_vector_t data;
  long int n1, n2, n3, n1n2;
} igraph_array3_t;

#define ARRAY3(m,i,j,k) ((m).data.stor_begin[(m).n1n2*(k)+(m).n1*(j)+(i)])
int igraph_array3_init(igraph_array3_t *a, long int n1, long int n2, 
		       long int n3);
void igraph_array3_destroy(igraph_array3_t *a);
long int igraph_array3_size(const igraph_array3_t *a);
long int igraph_array3_n(const igraph_array3_t *a, long int idx);
int igraph_array3_resize(igraph_array3_t *a, long int n1, long int n2, 
			 long int n3);
void igraph_array3_null(igraph_array3_t *a);

/* -------------------------------------------------- */
/* Plain stack                                        */
/* -------------------------------------------------- */

/**
 * Stack data type.
 * \ingroup internal
 */

typedef struct s_stack {
  igraph_real_t* stor_begin;
  igraph_real_t* stor_end;
  igraph_real_t* end;
} igraph_stack_t;

#define IGRAPH_STACK_NULL { 0,0,0 }

int igraph_stack_init       (igraph_stack_t* s, long int size);
void igraph_stack_destroy    (igraph_stack_t* s);
int igraph_stack_reserve    (igraph_stack_t* s, long int size);
igraph_bool_t igraph_stack_empty      (igraph_stack_t* s);
long int igraph_stack_size       (igraph_stack_t* s);
void igraph_stack_clear      (igraph_stack_t* s);
int igraph_stack_push       (igraph_stack_t* s, igraph_real_t elem);
igraph_real_t igraph_stack_pop        (igraph_stack_t* s);
igraph_real_t igraph_stack_top        (const igraph_stack_t* s);

/* -------------------------------------------------- */
/* Heap                                               */
/* -------------------------------------------------- */

/**
 * Heap data type.
 * \ingroup internal
 */

typedef struct s_heap {
  igraph_real_t* stor_begin;
  igraph_real_t* stor_end;
  igraph_real_t* end;
  int destroy;
} igraph_heap_t;

#define IGRAPH_HEAP_NULL { 0,0,0 }

int igraph_heap_init           (igraph_heap_t* h, long int size);
int igraph_heap_init_array     (igraph_heap_t *t, igraph_real_t* data, long int len);
void igraph_heap_destroy        (igraph_heap_t* h);
igraph_bool_t igraph_heap_empty          (igraph_heap_t* h);
int igraph_heap_push           (igraph_heap_t* h, igraph_real_t elem);
igraph_real_t igraph_heap_max       (igraph_heap_t* h);
igraph_real_t igraph_heap_delete_max(igraph_heap_t* h);
long int igraph_heap_size      (igraph_heap_t* h);
int igraph_heap_reserve        (igraph_heap_t* h, long int size);

void igraph_heap_i_build(igraph_real_t* arr, long int size, long int head);
void igraph_heap_i_shift_up(igraph_real_t* arr, long int size, long int elem);
void igraph_heap_i_sink(igraph_real_t* arr, long int size, long int head);
void igraph_heap_i_switch(igraph_real_t* arr, long int e1, long int e2);

/* -------------------------------------------------- */
/* Indexed heap                                       */
/* -------------------------------------------------- */

/**
 * Indexed heap data type.
 * \ingroup internal
 */

typedef struct s_indheap {
  igraph_real_t* stor_begin;
  igraph_real_t* stor_end;
  igraph_real_t* end;
  int destroy;
  long int* index_begin;
} igraph_indheap_t;

#define IGRAPH_INDHEAP_NULL { 0,0,0,0,0 }

int igraph_indheap_init           (igraph_indheap_t* h, long int size);
int igraph_indheap_init_array     (igraph_indheap_t *t, igraph_real_t* data, long int len);
void igraph_indheap_destroy        (igraph_indheap_t* h);
igraph_bool_t igraph_indheap_empty          (igraph_indheap_t* h);
int igraph_indheap_push           (igraph_indheap_t* h, igraph_real_t elem);
igraph_real_t igraph_indheap_max       (igraph_indheap_t* h);
igraph_real_t igraph_indheap_delete_max(igraph_indheap_t* h);
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
  igraph_real_t* stor_begin;
  igraph_real_t* stor_end;
  igraph_real_t* end;
  int destroy;
  long int* index_begin;
  long int* index2_begin;
} igraph_d_indheap_t;


#define IGRAPH_D_INDHEAP_NULL { 0,0,0,0,0,0 }

int igraph_d_indheap_init           (igraph_d_indheap_t* h, long int size);
void igraph_d_indheap_destroy        (igraph_d_indheap_t* h);
igraph_bool_t igraph_d_indheap_empty          (igraph_d_indheap_t* h);
int igraph_d_indheap_push           (igraph_d_indheap_t* h, igraph_real_t elem, 
			      long int idx, long int idx2);
igraph_real_t igraph_d_indheap_max       (igraph_d_indheap_t* h);
igraph_real_t igraph_d_indheap_delete_max(igraph_d_indheap_t* h);
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
  igraph_bool_t storekeys;
  igraph_strvector_t keys;
} igraph_trie_t;

#define IGRAPH_TRIE_NULL { IGRAPH_STRVECTOR_NULL, IGRAPH_VECTOR_PTR_NULL, \
                           IGRAPH_VECTOR_NULL, 0, 0, IGRAPH_STRVECTOR_NULL }
#define IGRAPH_TRIE_INIT_FINALLY(tr, sk) \
  do { IGRAPH_CHECK(igraph_trie_init(tr, sk)); \
  IGRAPH_FINALLY(igraph_trie_destroy, tr); } while (0)

int igraph_trie_init(igraph_trie_t *t, igraph_bool_t storekeys);
void igraph_trie_destroy(igraph_trie_t *t);
int igraph_trie_get(igraph_trie_t *t, const char *key, long int *id);
int igraph_trie_check(igraph_trie_t *t, const char *key, long int *id);
int igraph_trie_get2(igraph_trie_t *t, const char *key, long int length, 
		     long int *id);
void igraph_trie_idx(igraph_trie_t *t, long int idx, char **str);
int igraph_trie_getkeys(igraph_trie_t *t, const igraph_strvector_t **strv);
long int igraph_trie_size(igraph_trie_t *t);

typedef struct {
  igraph_vector_t v;
  long int size;
  long int offset;
} igraph_psumtree_t;
int igraph_psumtree_init(igraph_psumtree_t *t, long int size);
void igraph_psumtree_destroy(igraph_psumtree_t *t);
igraph_real_t igraph_psumtree_get(const igraph_psumtree_t *t, long int idx);
long int igraph_psumtree_size(const igraph_psumtree_t *t);
int igraph_psumtree_search(const igraph_psumtree_t *t, long int *idx,
			   igraph_real_t elem);
int igraph_psumtree_update(igraph_psumtree_t *t, long int idx, 
			   igraph_real_t new_value);
igraph_real_t igraph_psumtree_sum(const igraph_psumtree_t *t);

/**
 * 2d grid containing points
 */

typedef struct igraph_2dgrid_t {
  igraph_matrix_t *coords;
  igraph_real_t minx, maxx, deltax;
  igraph_real_t miny, maxy, deltay;
  long int stepsx, stepsy;
  igraph_matrix_t startidx;
  igraph_vector_t next;
  igraph_vector_t prev;
  igraph_real_t massx, massy;		/* The sum of the coordinates */
  long int vertices;		/* Number of active vertices  */
} igraph_2dgrid_t;

int igraph_2dgrid_init(igraph_2dgrid_t *grid, igraph_matrix_t *coords, 
		       igraph_real_t minx, igraph_real_t maxx, igraph_real_t deltax,
		       igraph_real_t miny, igraph_real_t maxy, igraph_real_t deltay);
void igraph_2dgrid_destroy(igraph_2dgrid_t *grid);
void igraph_2dgrid_add(igraph_2dgrid_t *grid, long int elem, 
		       igraph_real_t xc, igraph_real_t yc);
void igraph_2dgrid_add2(igraph_2dgrid_t *grid, long int elem);
void igraph_2dgrid_move(igraph_2dgrid_t *grid, long int elem, 
			igraph_real_t xc, igraph_real_t yc);
void igraph_2dgrid_getcenter(const igraph_2dgrid_t *grid, 
			     igraph_real_t *massx, igraph_real_t *massy);
igraph_bool_t igraph_2dgrid_in(const igraph_2dgrid_t *grid, long int elem);
igraph_real_t igraph_2dgrid_dist(const igraph_2dgrid_t *grid, 
			  long int e1, long int e2);
int igraph_2dgrid_neighbors(igraph_2dgrid_t *grid, igraph_vector_t *eids, 
			    igraph_integer_t vid, igraph_real_t r);

typedef struct igraph_2dgrid_iterator_t {
  long int vid, x, y;
  long int nei;
  long int nx[4], ny[4], ncells;
} igraph_2dgrid_iterator_t;

void igraph_2dgrid_reset(igraph_2dgrid_t *grid, igraph_2dgrid_iterator_t *it);
igraph_integer_t igraph_2dgrid_next(igraph_2dgrid_t *grid, 
			      igraph_2dgrid_iterator_t *it);
igraph_integer_t igraph_2dgrid_next_nei(igraph_2dgrid_t *grid,
				 igraph_2dgrid_iterator_t *it);

/* Another type of grid, each cell is owned by exactly one graph */

typedef struct igraph_i_layout_mergegrid_t {
  long int *data;
  long int stepsx, stepsy;
  igraph_real_t minx, maxx, deltax;
  igraph_real_t miny, maxy, deltay;
} igraph_i_layout_mergegrid_t;

int igraph_i_layout_mergegrid_init(igraph_i_layout_mergegrid_t *grid,
				   igraph_real_t minx, igraph_real_t maxx, long int stepsx,
				   igraph_real_t miny, igraph_real_t maxy, long int stepsy);
void igraph_i_layout_mergegrid_destroy(igraph_i_layout_mergegrid_t *grid);

int igraph_i_layout_merge_place_sphere(igraph_i_layout_mergegrid_t *grid,
				       igraph_real_t x, igraph_real_t y, igraph_real_t r,
				       long int id);

long int igraph_i_layout_mergegrid_get(igraph_i_layout_mergegrid_t *grid,
				       igraph_real_t x, igraph_real_t y);

long int igraph_i_layout_mergegrid_get_sphere(igraph_i_layout_mergegrid_t *g,
					      igraph_real_t x, igraph_real_t y, igraph_real_t r);

/* string -> string hash table */

typedef struct igraph_hashtable_t {
  igraph_trie_t keys;
  igraph_strvector_t elements;
  igraph_strvector_t defaults;
} igraph_hashtable_t;

int igraph_hashtable_init(igraph_hashtable_t *ht);
void igraph_hashtable_destroy(igraph_hashtable_t *ht);
int igraph_hashtable_addset(igraph_hashtable_t *ht,
			    const char *key, const char *def, 
			    const char *elem);
int igraph_hashtable_addset2(igraph_hashtable_t *ht,
			     const char *key, const char *def,
			     const char *elem, int elemlen);
int igraph_hashtable_get(igraph_hashtable_t *ht,
			 const char *key, char **elem);
int igraph_hashtable_getkeys(igraph_hashtable_t *ht, 
			     const igraph_strvector_t **sv);
int igraph_hashtable_reset(igraph_hashtable_t *ht);

/* Buckets, needed for the maximum flow algorithm */

typedef struct igraph_buckets_t {
  igraph_vector_t bptr;
  igraph_vector_t buckets;
  igraph_integer_t max, no;
} igraph_buckets_t;

int igraph_buckets_init(igraph_buckets_t *b, long int bsize, long int size);
void igraph_buckets_destroy(igraph_buckets_t *b);
long int igraph_buckets_popmax(igraph_buckets_t *b);
igraph_bool_t igraph_buckets_empty(const igraph_buckets_t *b);
void igraph_buckets_add(igraph_buckets_t *b, long int bucket,
		       igraph_real_t elem);

/* Special maximum heap, needed for the minimum cut algorithm */

typedef struct igraph_i_cutheap_t {
  igraph_vector_t heap;
  igraph_vector_t index;
  igraph_vector_t hptr;
  long int dnodes;
} igraph_i_cutheap_t;

int igraph_i_cutheap_init(igraph_i_cutheap_t *ch, igraph_integer_t nodes);
void igraph_i_cutheap_destroy(igraph_i_cutheap_t *ch);
igraph_bool_t igraph_i_cutheap_empty(igraph_i_cutheap_t *ch);
igraph_integer_t igraph_i_cutheap_active_size(igraph_i_cutheap_t *ch);
igraph_integer_t igraph_i_cutheap_size(igraph_i_cutheap_t *ch);
igraph_real_t igraph_i_cutheap_maxvalue(igraph_i_cutheap_t *ch);
igraph_integer_t igraph_i_cutheap_popmax(igraph_i_cutheap_t *ch);
int igraph_i_cutheap_update(igraph_i_cutheap_t *ch, igraph_integer_t index,
			    igraph_real_t add);
int igraph_i_cutheap_reset_undefine(igraph_i_cutheap_t *ch, long int vertex);

/* -------------------------------------------------- */
/* Flexible set                                       */
/* -------------------------------------------------- */

/** 
 * Set containing integer numbers regardless of the order
 * \ingroup types
 */

typedef struct s_set {
  igraph_real_t* stor_begin;
  igraph_real_t* stor_end;
  igraph_real_t* end;
} igraph_set_t;

#define IGRAPH_SET_NULL { 0,0,0 }
#define IGRAPH_SET_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(igraph_set_init(v, size)); \
  IGRAPH_FINALLY(igraph_set_destroy, v); } while (0)

int igraph_set_init      (igraph_set_t* set, long int size);
void igraph_set_destroy   (igraph_set_t* set);
igraph_bool_t igraph_set_inited   (igraph_set_t* set);
int igraph_set_reserve   (igraph_set_t* set, long int size);
igraph_bool_t igraph_set_empty     (const igraph_set_t* set);
void igraph_set_clear      (igraph_set_t* set);
long int igraph_set_size      (const igraph_set_t* set);
int igraph_set_add (igraph_set_t* v, igraph_integer_t e);
igraph_bool_t igraph_set_contains (igraph_set_t* set, igraph_integer_t e);
igraph_bool_t igraph_set_iterate (igraph_set_t* set, long int* state,
				  igraph_integer_t* element);

__END_DECLS

#endif

