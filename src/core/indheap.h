/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2020  The igraph development team

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

#ifndef IGRAPH_CORE_INDHEAP_H
#define IGRAPH_CORE_INDHEAP_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_vector.h"

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

IGRAPH_PRIVATE_EXPORT int igraph_d_indheap_init(igraph_d_indheap_t *h, long int size);
IGRAPH_PRIVATE_EXPORT void igraph_d_indheap_destroy(igraph_d_indheap_t *h);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_d_indheap_empty(igraph_d_indheap_t *h);
IGRAPH_PRIVATE_EXPORT int igraph_d_indheap_push(igraph_d_indheap_t *h, igraph_real_t elem,
                                                long int idx, long int idx2);
IGRAPH_PRIVATE_EXPORT igraph_real_t igraph_d_indheap_max(igraph_d_indheap_t *h);
IGRAPH_PRIVATE_EXPORT igraph_real_t igraph_d_indheap_delete_max(igraph_d_indheap_t *h);
IGRAPH_PRIVATE_EXPORT long int igraph_d_indheap_size(igraph_d_indheap_t *h);
IGRAPH_PRIVATE_EXPORT int igraph_d_indheap_reserve(igraph_d_indheap_t *h, long int size);
IGRAPH_PRIVATE_EXPORT void igraph_d_indheap_max_index(igraph_d_indheap_t *h, long int *idx, long int *idx2);

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

IGRAPH_PRIVATE_EXPORT int igraph_2wheap_init(igraph_2wheap_t *h, long int size);
IGRAPH_PRIVATE_EXPORT void igraph_2wheap_destroy(igraph_2wheap_t *h);
IGRAPH_PRIVATE_EXPORT int igraph_2wheap_clear(igraph_2wheap_t *h);
IGRAPH_PRIVATE_EXPORT int igraph_2wheap_push_with_index(igraph_2wheap_t *h,
                                                        long int idx, igraph_real_t elem);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_2wheap_empty(const igraph_2wheap_t *h);
IGRAPH_PRIVATE_EXPORT long int igraph_2wheap_size(const igraph_2wheap_t *h);
IGRAPH_PRIVATE_EXPORT long int igraph_2wheap_max_size(const igraph_2wheap_t *h);
IGRAPH_PRIVATE_EXPORT igraph_real_t igraph_2wheap_max(const igraph_2wheap_t *h);
IGRAPH_PRIVATE_EXPORT long int igraph_2wheap_max_index(const igraph_2wheap_t *h);
IGRAPH_PRIVATE_EXPORT igraph_real_t igraph_2wheap_deactivate_max(igraph_2wheap_t *h);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_2wheap_has_elem(const igraph_2wheap_t *h, long int idx);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_2wheap_has_active(const igraph_2wheap_t *h, long int idx);
IGRAPH_PRIVATE_EXPORT igraph_real_t igraph_2wheap_get(const igraph_2wheap_t *h, long int idx);
IGRAPH_PRIVATE_EXPORT igraph_real_t igraph_2wheap_delete_max(igraph_2wheap_t *h);
IGRAPH_PRIVATE_EXPORT igraph_real_t igraph_2wheap_delete_max_index(igraph_2wheap_t *h, long int *idx);
IGRAPH_PRIVATE_EXPORT int igraph_2wheap_modify(igraph_2wheap_t *h, long int idx, igraph_real_t elem);
IGRAPH_PRIVATE_EXPORT int igraph_2wheap_check(igraph_2wheap_t *h);

__END_DECLS

#endif
