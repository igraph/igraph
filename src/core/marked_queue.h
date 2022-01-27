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

#ifndef IGRAPH_MARKED_QUEUE_H
#define IGRAPH_MARKED_QUEUE_H

#include "igraph_decls.h"
#include "igraph_vector.h"
#include "igraph_dqueue.h"

#include <stdio.h>

__BEGIN_DECLS

/* This is essentially a double ended queue, with some extra features:
   (1) The is-element? operation is fast, O(1). This requires that we
       know a limit for the number of elements in the queue.
   (2) We can insert elements in batches, and the whole batch can be
      removed at once.

  Currently only the top-end operations are implemented, so the queue
  is essentially a stack.
*/

typedef struct igraph_marked_queue_int_t {
    igraph_dqueue_int_t Q;
    igraph_vector_int_t set;
    igraph_integer_t mark;
    igraph_integer_t size;
} igraph_marked_queue_int_t;

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_marked_queue_int_init(igraph_marked_queue_int_t *q,
                                                              igraph_integer_t size);
IGRAPH_PRIVATE_EXPORT void igraph_marked_queue_int_destroy(igraph_marked_queue_int_t *q);
IGRAPH_PRIVATE_EXPORT void igraph_marked_queue_int_reset(igraph_marked_queue_int_t *q);

IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_marked_queue_int_empty(const igraph_marked_queue_int_t *q);
IGRAPH_PRIVATE_EXPORT igraph_integer_t igraph_marked_queue_int_size(const igraph_marked_queue_int_t *q);

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_marked_queue_int_print(const igraph_marked_queue_int_t *q);
IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_marked_queue_int_fprint(const igraph_marked_queue_int_t *q, FILE *file);

IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_marked_queue_int_iselement(const igraph_marked_queue_int_t *q,
                                                                  igraph_integer_t elem);
IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_marked_queue_int_push(igraph_marked_queue_int_t *q, igraph_integer_t elem);

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_marked_queue_int_start_batch(igraph_marked_queue_int_t *q);
IGRAPH_PRIVATE_EXPORT void igraph_marked_queue_int_pop_back_batch(igraph_marked_queue_int_t *q);

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_marked_queue_int_as_vector(const igraph_marked_queue_int_t *q,
                                                                   igraph_vector_int_t *vec);

__END_DECLS

#endif
