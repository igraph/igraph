/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#ifndef IGRAPH_TYPES_INTERNAL_H
#define IGRAPH_TYPES_INTERNAL_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
    #define __BEGIN_DECLS extern "C" {
    #define __END_DECLS }
#else
    #define __BEGIN_DECLS /* empty */
    #define __END_DECLS /* empty */
#endif

#include "igraph_types.h"
#include "igraph_matrix.h"
#include "igraph_stack.h"
#include "igraph_strvector.h"
#include "igraph_vector.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

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
int igraph_indheap_clear(igraph_indheap_t *h);
igraph_bool_t igraph_indheap_empty          (igraph_indheap_t* h);
int igraph_indheap_push           (igraph_indheap_t* h, igraph_real_t elem);
int igraph_indheap_push_with_index(igraph_indheap_t* h, long int idx, igraph_real_t elem);
int igraph_indheap_modify(igraph_indheap_t* h, long int idx, igraph_real_t elem);
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

/* -------------------------------------------------- */
/* Two-way indexed heap                               */
/* -------------------------------------------------- */

/* This is a smart indexed heap. In addition to the "normal" indexed heap
   it allows to access every element through its index in O(1) time.
   In other words, for this heap the _modify operation is O(1), the
   normal heap does this in O(n) time.... */

typedef struct igraph_2wheap_t {
    long int size;
    igraph_vector_t data;
    igraph_vector_long_t index;
    igraph_vector_long_t index2;
} igraph_2wheap_t;

int igraph_2wheap_init(igraph_2wheap_t *h, long int size);
void igraph_2wheap_destroy(igraph_2wheap_t *h);
int igraph_2wheap_clear(igraph_2wheap_t *h);
int igraph_2wheap_push_with_index(igraph_2wheap_t *h,
                                  long int idx, igraph_real_t elem);
igraph_bool_t igraph_2wheap_empty(const igraph_2wheap_t *h);
long int igraph_2wheap_size(const igraph_2wheap_t *h);
long int igraph_2wheap_max_size(const igraph_2wheap_t *h);
igraph_real_t igraph_2wheap_max(const igraph_2wheap_t *h);
long int igraph_2wheap_max_index(const igraph_2wheap_t *h);
igraph_real_t igraph_2wheap_deactivate_max(igraph_2wheap_t *h);
igraph_bool_t igraph_2wheap_has_elem(const igraph_2wheap_t *h, long int idx);
igraph_bool_t igraph_2wheap_has_active(const igraph_2wheap_t *h, long int idx);
igraph_real_t igraph_2wheap_get(const igraph_2wheap_t *h, long int idx);
igraph_real_t igraph_2wheap_delete_max(igraph_2wheap_t *h);
igraph_real_t igraph_2wheap_delete_max_index(igraph_2wheap_t *h, long int *idx);
int igraph_2wheap_modify(igraph_2wheap_t *h, long int idx, igraph_real_t elem);
int igraph_2wheap_check(igraph_2wheap_t *h);

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
    igraph_real_t massx, massy;       /* The sum of the coordinates */
    long int vertices;        /* Number of active vertices  */
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
    igraph_vector_long_t bptr;
    igraph_vector_long_t buckets;
    igraph_integer_t max, no;
} igraph_buckets_t;

int igraph_buckets_init(igraph_buckets_t *b, long int bsize, long int size);
void igraph_buckets_destroy(igraph_buckets_t *b);
void igraph_buckets_clear(igraph_buckets_t *b);
long int igraph_buckets_popmax(igraph_buckets_t *b);
long int igraph_buckets_pop(igraph_buckets_t *b, long int bucket);
igraph_bool_t igraph_buckets_empty(const igraph_buckets_t *b);
igraph_bool_t igraph_buckets_empty_bucket(const igraph_buckets_t *b,
        long int bucket);
void igraph_buckets_add(igraph_buckets_t *b, long int bucket,
                        long int elem);

typedef struct igraph_dbuckets_t {
    igraph_vector_long_t bptr;
    igraph_vector_long_t next, prev;
    igraph_integer_t max, no;
} igraph_dbuckets_t;

int igraph_dbuckets_init(igraph_dbuckets_t *b, long int bsize, long int size);
void igraph_dbuckets_destroy(igraph_dbuckets_t *b);
void igraph_dbuckets_clear(igraph_dbuckets_t *b);
long int igraph_dbuckets_popmax(igraph_dbuckets_t *b);
long int igraph_dbuckets_pop(igraph_dbuckets_t *b, long int bucket);
igraph_bool_t igraph_dbuckets_empty(const igraph_dbuckets_t *b);
igraph_bool_t igraph_dbuckets_empty_bucket(const igraph_dbuckets_t *b,
        long int bucket);
void igraph_dbuckets_add(igraph_dbuckets_t *b, long int bucket,
                         long int elem);
void igraph_dbuckets_delete(igraph_dbuckets_t *b, long int bucket,
                            long int elem);

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
    igraph_integer_t* stor_begin;
    igraph_integer_t* stor_end;
    igraph_integer_t* end;
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

/* -------------------------------------------------- */
/* Vectorlist, fixed length                           */
/* -------------------------------------------------- */

typedef struct igraph_fixed_vectorlist_t {
    igraph_vector_t *vecs;
    igraph_vector_ptr_t v;
    long int length;
} igraph_fixed_vectorlist_t;

void igraph_fixed_vectorlist_destroy(igraph_fixed_vectorlist_t *l);
int igraph_fixed_vectorlist_convert(igraph_fixed_vectorlist_t *l,
                                    const igraph_vector_t *from,
                                    long int size);

__END_DECLS

#endif
