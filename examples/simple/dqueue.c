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

#include <igraph.h>

int main() {

    igraph_dqueue_t q;
    int i;

    /* igraph_dqueue_init, igraph_dqueue_destroy, igraph_dqueue_empty */
    igraph_dqueue_init(&q, 5);
    if (!igraph_dqueue_empty(&q)) {
        return 1;
    }
    igraph_dqueue_destroy(&q);

    /* igraph_dqueue_push, igraph_dqueue_pop */
    igraph_dqueue_init(&q, 4);
    igraph_dqueue_push(&q, 1);
    igraph_dqueue_push(&q, 2);
    igraph_dqueue_push(&q, 3);
    igraph_dqueue_push(&q, 4);
    if (igraph_dqueue_pop(&q) != 1) {
        return 2;
    }
    if (igraph_dqueue_pop(&q) != 2) {
        return 3;
    }
    if (igraph_dqueue_pop(&q) != 3) {
        return 4;
    }
    if (igraph_dqueue_pop(&q) != 4) {
        return 5;
    }
    igraph_dqueue_destroy(&q);

    /* igraph_dqueue_clear, igraph_dqueue_size */
    igraph_dqueue_init(&q, 0);
    if (igraph_dqueue_size(&q) != 0) {
        return 6;
    }
    igraph_dqueue_clear(&q);
    if (igraph_dqueue_size(&q) != 0) {
        return 7;
    }
    for (i = 0; i < 10; i++) {
        igraph_dqueue_push(&q, i);
    }
    igraph_dqueue_clear(&q);
    if (igraph_dqueue_size(&q) != 0) {
        return 8;
    }
    igraph_dqueue_destroy(&q);

    /* TODO: igraph_dqueue_full */

    /* igraph_dqueue_head, igraph_dqueue_back, igraph_dqueue_pop_back */
    igraph_dqueue_init(&q, 0);
    for (i = 0; i < 10; i++) {
        igraph_dqueue_push(&q, i);
    }
    for (i = 0; i < 10; i++) {
        if (igraph_dqueue_head(&q) != 0) {
            return 9;
        }
        if (igraph_dqueue_back(&q) != 9 - i) {
            return 10;
        }
        if (igraph_dqueue_pop_back(&q) != 9 - i) {
            return 11;
        }
    }
    igraph_dqueue_destroy(&q);

    /* print */
    igraph_dqueue_init(&q, 4);
    igraph_dqueue_push(&q, 1);
    igraph_dqueue_push(&q, 2);
    igraph_dqueue_push(&q, 3);
    igraph_dqueue_push(&q, 4);
    igraph_dqueue_pop(&q);
    igraph_dqueue_pop(&q);
    igraph_dqueue_push(&q, 5);
    igraph_dqueue_push(&q, 6);
    igraph_dqueue_print(&q);

    igraph_dqueue_clear(&q);
    igraph_dqueue_print(&q);

    igraph_dqueue_destroy(&q);

    if (IGRAPH_FINALLY_STACK_SIZE() != 0) {
        return 12;
    }

    return 0;
}
