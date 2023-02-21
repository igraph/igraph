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

#ifndef IGRAPH_CORE_BUCKETS_H
#define IGRAPH_CORE_BUCKETS_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/* Buckets, needed for the maximum flow algorithm */

typedef struct igraph_buckets_t {
    igraph_vector_int_t bptr;
    igraph_vector_int_t buckets;
    igraph_integer_t max, no;
} igraph_buckets_t;

igraph_error_t igraph_buckets_init(igraph_buckets_t *b, igraph_integer_t bsize, igraph_integer_t size);
void igraph_buckets_destroy(igraph_buckets_t *b);
void igraph_buckets_clear(igraph_buckets_t *b);
igraph_integer_t igraph_buckets_popmax(igraph_buckets_t *b);
igraph_integer_t igraph_buckets_pop(igraph_buckets_t *b, igraph_integer_t bucket);
igraph_bool_t igraph_buckets_empty(const igraph_buckets_t *b);
igraph_bool_t igraph_buckets_empty_bucket(const igraph_buckets_t *b,
        igraph_integer_t bucket);
void igraph_buckets_add(igraph_buckets_t *b, igraph_integer_t bucket,
                        igraph_integer_t elem);

typedef struct igraph_dbuckets_t {
    igraph_vector_int_t bptr;
    igraph_vector_int_t next, prev;
    igraph_integer_t max, no;
} igraph_dbuckets_t;

igraph_error_t igraph_dbuckets_init(igraph_dbuckets_t *b, igraph_integer_t bsize, igraph_integer_t size);
void igraph_dbuckets_destroy(igraph_dbuckets_t *b);
void igraph_dbuckets_clear(igraph_dbuckets_t *b);
igraph_integer_t igraph_dbuckets_popmax(igraph_dbuckets_t *b);
igraph_integer_t igraph_dbuckets_pop(igraph_dbuckets_t *b, igraph_integer_t bucket);
igraph_bool_t igraph_dbuckets_empty(const igraph_dbuckets_t *b);
igraph_bool_t igraph_dbuckets_empty_bucket(const igraph_dbuckets_t *b,
        igraph_integer_t bucket);
void igraph_dbuckets_add(igraph_dbuckets_t *b, igraph_integer_t bucket,
                         igraph_integer_t elem);
void igraph_dbuckets_delete(igraph_dbuckets_t *b, igraph_integer_t bucket,
                            igraph_integer_t elem);

__END_DECLS

#endif
