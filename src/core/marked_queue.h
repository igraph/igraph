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

#include "igraph_vector.h"
#include "igraph_dqueue.h"

#include <stdio.h>

/* This is essentially a double ended queue, with some extra features:
   (1) The is-element? operation is fast, O(1). This requires that we
       know a limit for the number of elements in the queue.
   (2) We can insert elements in batches, and the whole batch can be
      removed at once.

  Currently only the top-end operations are implemented, so the queue
  is essentially a stack.
*/

typedef struct igraph_marked_queue_t {
    igraph_dqueue_t Q;
    igraph_vector_long_t set;
    long int mark;
    long int size;
} igraph_marked_queue_t;

int igraph_marked_queue_init(igraph_marked_queue_t *q,
                             long int size);
void igraph_marked_queue_destroy(igraph_marked_queue_t *q);
void igraph_marked_queue_reset(igraph_marked_queue_t *q);

igraph_bool_t igraph_marked_queue_empty(const igraph_marked_queue_t *q);
long int igraph_marked_queue_size(const igraph_marked_queue_t *q);
int igraph_marked_queue_print(const igraph_marked_queue_t *q);
int igraph_marked_queue_fprint(const igraph_marked_queue_t *q, FILE *file);

igraph_bool_t igraph_marked_queue_iselement(const igraph_marked_queue_t *q,
        long int elem);

int igraph_marked_queue_push(igraph_marked_queue_t *q, long int elem);

int igraph_marked_queue_start_batch(igraph_marked_queue_t *q);
void igraph_marked_queue_pop_back_batch(igraph_marked_queue_t *q);

int igraph_marked_queue_as_vector(const igraph_marked_queue_t *q,
                                  igraph_vector_t *vec);

#endif
