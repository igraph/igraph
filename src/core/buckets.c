/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "core/buckets.h"

/* The igraph_buckets_t data structure can store at most 'size'
 * unique integers in 'bsize' buckets. It has the following simple
 * operations (in addition to _init() and _destroy():
 * - _add() adding an element to the given bucket.
 * - _popmax() removing an element from the bucket with the highest
 *   id.
 *   Currently buckets work as stacks, last-in-first-out mode.
 * - _empty() queries whether the buckets is empty.
 *
 * Internal representation: we use a vector to create single linked
 * lists, and another vector that points to the starting element of
 * each bucket. Zero means the end of the chain. So bucket i contains
 * elements bptr[i], buckets[bptr[i]], buckets[buckets[bptr[i]]],
 * etc., until a zero is found.
 *
 * We also keep the total number of elements in the buckets and the
 * id of the non-empty bucket with the highest id, to facilitate the
 * _empty() and _popmax() operations.
 */

igraph_error_t igraph_buckets_init(igraph_buckets_t *b, igraph_integer_t bsize, igraph_integer_t size) {
    IGRAPH_VECTOR_INT_INIT_FINALLY(&b->bptr, bsize);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&b->buckets, size);
    b->max = -1; b->no = 0;
    IGRAPH_FINALLY_CLEAN(2);
    return IGRAPH_SUCCESS;
}

void igraph_buckets_destroy(igraph_buckets_t *b) {
    igraph_vector_int_destroy(&b->bptr);
    igraph_vector_int_destroy(&b->buckets);
}

igraph_integer_t igraph_buckets_popmax(igraph_buckets_t *b) {
    /* Precondition: there is at least a non-empty bucket */
    /* Search for the highest bucket first */
    igraph_integer_t max;
    while ( (max = VECTOR(b->bptr)[b->max]) == 0) {
        b->max --;
    }
    VECTOR(b->bptr)[b->max] = VECTOR(b->buckets)[max - 1];
    b->no--;

    return max - 1;
}

igraph_integer_t igraph_buckets_pop(igraph_buckets_t *b, igraph_integer_t bucket) {
    igraph_integer_t ret = VECTOR(b->bptr)[bucket] - 1;
    VECTOR(b->bptr)[bucket] = VECTOR(b->buckets)[ret];
    b->no--;
    return ret;
}

igraph_bool_t igraph_buckets_empty(const igraph_buckets_t *b) {
    return (b->no == 0);
}

igraph_bool_t igraph_buckets_empty_bucket(const igraph_buckets_t *b,
        igraph_integer_t bucket) {
    return VECTOR(b->bptr)[bucket] == 0;
}

void igraph_buckets_add(igraph_buckets_t *b, igraph_integer_t bucket,
                        igraph_integer_t elem) {

    VECTOR(b->buckets)[elem] = VECTOR(b->bptr)[bucket];
    VECTOR(b->bptr)[bucket] = elem + 1;
    if (bucket > b->max) {
        b->max = bucket;
    }
    b->no++;
}

void igraph_buckets_clear(igraph_buckets_t *b) {
    igraph_vector_int_null(&b->bptr);
    igraph_vector_int_null(&b->buckets);
    b->max = -1;
    b->no = 0;
}

igraph_error_t igraph_dbuckets_init(igraph_dbuckets_t *b, igraph_integer_t bsize, igraph_integer_t size) {
    IGRAPH_VECTOR_INT_INIT_FINALLY(&b->bptr, bsize);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&b->next, size);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&b->prev, size);
    b->max = -1; b->no = 0;
    IGRAPH_FINALLY_CLEAN(3);
    return IGRAPH_SUCCESS;
}

void igraph_dbuckets_destroy(igraph_dbuckets_t *b) {
    igraph_vector_int_destroy(&b->bptr);
    igraph_vector_int_destroy(&b->next);
    igraph_vector_int_destroy(&b->prev);
}

void igraph_dbuckets_clear(igraph_dbuckets_t *b) {
    igraph_vector_int_null(&b->bptr);
    igraph_vector_int_null(&b->next);
    igraph_vector_int_null(&b->prev);
    b->max = -1;
    b->no = 0;
}

igraph_integer_t igraph_dbuckets_popmax(igraph_dbuckets_t *b) {
    while ( VECTOR(b->bptr)[b->max] == 0) {
        b->max--;
    }
    return igraph_dbuckets_pop(b, b->max);
}

igraph_integer_t igraph_dbuckets_pop(igraph_dbuckets_t *b, igraph_integer_t bucket) {
    igraph_integer_t ret = VECTOR(b->bptr)[bucket] - 1;
    igraph_integer_t next = VECTOR(b->next)[ret];
    VECTOR(b->bptr)[bucket] = next;
    if (next != 0) {
        VECTOR(b->prev)[next - 1] = 0;
    }

    b->no--;
    return ret;
}

igraph_bool_t igraph_dbuckets_empty(const igraph_dbuckets_t *b) {
    return (b->no == 0);
}

igraph_bool_t igraph_dbuckets_empty_bucket(const igraph_dbuckets_t *b,
        igraph_integer_t bucket) {
    return VECTOR(b->bptr)[bucket] == 0;
}

void igraph_dbuckets_add(igraph_dbuckets_t *b, igraph_integer_t bucket,
                         igraph_integer_t elem) {
    igraph_integer_t oldfirst = VECTOR(b->bptr)[bucket];
    VECTOR(b->bptr)[bucket] = elem + 1;
    VECTOR(b->next)[elem] = oldfirst;
    if (oldfirst != 0) {
        VECTOR(b->prev)[oldfirst - 1] = elem + 1;
    }
    if (bucket > b->max) {
        b->max = bucket;
    }
    b->no++;
}

/* Remove an arbitrary element */

void igraph_dbuckets_delete(igraph_dbuckets_t *b, igraph_integer_t bucket,
                            igraph_integer_t elem) {
    if (VECTOR(b->bptr)[bucket] == elem + 1) {
        /* First element in bucket */
        igraph_integer_t next = VECTOR(b->next)[elem];
        if (next != 0) {
            VECTOR(b->prev)[next - 1] = 0;
        }
        VECTOR(b->bptr)[bucket] = next;
    } else {
        igraph_integer_t next = VECTOR(b->next)[elem];
        igraph_integer_t prev = VECTOR(b->prev)[elem];
        if (next != 0) {
            VECTOR(b->prev)[next - 1] = prev;
        }
        if (prev != 0) {
            VECTOR(b->next)[prev - 1] = next;
        }
    }
    b->no--;
}
