/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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

#include <igraph.h>
#include "core/marked_queue.h"

#include "test_utilities.inc"

int main() {
    igraph_marked_queue_t Q;
    long int i;

    igraph_marked_queue_init(&Q, 100);
    for (i = 0; i < 50; i++) {
        igraph_marked_queue_push(&Q, i);
        if (!igraph_marked_queue_iselement(&Q, i)) {
            return 4;
        }
        if (! ((i + 1) % 5)) {
            igraph_marked_queue_start_batch(&Q);
        }
    }

    for (i = 1; i < 50; i++) {
        if (!igraph_marked_queue_iselement(&Q, i)) {
            printf("Problem with %li.\n", i);
            return 3;
        }
    }

    for (i = 0; i <= 50 / 5; i++) {
        if (igraph_marked_queue_empty(&Q)) {
            return 1;
        }
        igraph_marked_queue_pop_back_batch(&Q);
    }
    if (!igraph_marked_queue_empty(&Q)) {
        return 2;
    }

    igraph_marked_queue_destroy(&Q);

    VERIFY_FINALLY_STACK();

    return 0;
}
