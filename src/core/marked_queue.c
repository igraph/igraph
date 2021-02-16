/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "core/marked_queue.h"

#define BATCH_MARKER -1

int igraph_marked_queue_init(igraph_marked_queue_t *q,
                             long int size) {
    IGRAPH_CHECK(igraph_dqueue_init(&q->Q, 0));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &q->Q);
    IGRAPH_CHECK(igraph_vector_long_init(&q->set, size));
    q->mark = 1;
    q->size = 0;
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

void igraph_marked_queue_destroy(igraph_marked_queue_t *q) {
    igraph_vector_long_destroy(&q->set);
    igraph_dqueue_destroy(&q->Q);
}

void igraph_marked_queue_reset(igraph_marked_queue_t *q) {
    igraph_dqueue_clear(&q->Q);
    q->size = 0;
    q->mark += 1;
    if (q->mark == 0) {
        igraph_vector_long_null(&q->set);
        q->mark += 1;
    }
}

igraph_bool_t igraph_marked_queue_empty(const igraph_marked_queue_t *q) {
    return q->size == 0;
}

long int igraph_marked_queue_size(const igraph_marked_queue_t *q) {
    return q->size;
}

igraph_bool_t igraph_marked_queue_iselement(const igraph_marked_queue_t *q,
        long int elem) {
    return (VECTOR(q->set)[elem] == q->mark);
}

int igraph_marked_queue_push(igraph_marked_queue_t *q, long int elem) {
    if (VECTOR(q->set)[elem] != q->mark) {
        IGRAPH_CHECK(igraph_dqueue_push(&q->Q, elem));
        VECTOR(q->set)[elem] = q->mark;
        q->size += 1;
    }
    return 0;
}

int igraph_marked_queue_start_batch(igraph_marked_queue_t *q) {
    IGRAPH_CHECK(igraph_dqueue_push(&q->Q, BATCH_MARKER));
    return 0;
}

void igraph_marked_queue_pop_back_batch(igraph_marked_queue_t *q) {
    long int size = igraph_dqueue_size(&q->Q);
    long int elem;
    while (size > 0 &&
           (elem = (long int) igraph_dqueue_pop_back(&q->Q)) != BATCH_MARKER) {
        VECTOR(q->set)[elem] = 0;
        size--;
        q->size--;
    }
}

#ifndef USING_R
int igraph_marked_queue_print(const igraph_marked_queue_t *q) {
    IGRAPH_CHECK(igraph_dqueue_print(&q->Q));
    return 0;
}
#endif

int igraph_marked_queue_fprint(const igraph_marked_queue_t *q, FILE *file) {
    IGRAPH_CHECK(igraph_dqueue_fprint(&q->Q, file));
    return 0;
}

int igraph_marked_queue_as_vector(const igraph_marked_queue_t *q,
                                  igraph_vector_t *vec) {
    long int i, p, n = igraph_dqueue_size(&q->Q);
    IGRAPH_CHECK(igraph_vector_resize(vec, q->size));
    for (i = 0, p = 0; i < n; i++) {
        igraph_real_t e = igraph_dqueue_e(&q->Q, i);
        if (e != BATCH_MARKER) {
            VECTOR(*vec)[p++] = e;
        }
    }
    return 0;
}
