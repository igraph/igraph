/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.inc"

int main() {

    igraph_stack_t st;
    int i;

    /* igraph_stack_init, igraph_stack_destroy */
    igraph_stack_init(&st, 0);
    igraph_stack_destroy(&st);
    igraph_stack_init(&st, 10);
    igraph_stack_destroy(&st);

    /* igraph_stack_reserve */
    igraph_stack_init(&st, 0);
    igraph_stack_reserve(&st, 10);
    igraph_stack_reserve(&st, 5);

    /* igraph_stack_empty */
    if (!igraph_stack_empty(&st)) {
        return 1;
    }
    igraph_stack_push(&st, 1);
    if (igraph_stack_empty(&st)) {
        return 2;
    }

    /* igraph_stack_size */
    if (igraph_stack_size(&st) != 1) {
        return 3;
    }
    for (i = 0; i < 10; i++) {
        igraph_stack_push(&st, i);
    }
    if (igraph_stack_size(&st) != 11) {
        return 4;
    }

    /* igraph_stack_clear */
    igraph_stack_clear(&st);
    if (!igraph_stack_empty(&st)) {
        return 5;
    }
    igraph_stack_push(&st, 100);
    if (igraph_stack_pop(&st) != 100) {
        return 6;
    }
    igraph_stack_clear(&st);
    igraph_stack_clear(&st);

    /* igraph_stack_push, igraph_stack_pop */
    for (i = 0; i < 100; i++) {
        igraph_stack_push(&st, 100 - i);
    }
    for (i = 0; i < 100; i++) {
        if (igraph_stack_pop(&st) != i + 1) {
            return 7;
        }
    }
    if (!igraph_stack_empty(&st)) {
        return 8;
    }

    igraph_stack_destroy(&st);

    VERIFY_FINALLY_STACK();

    return 0;
}
