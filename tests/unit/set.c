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
#include <stdlib.h>

#include "core/set.h"

#include "test_utilities.inc"

void print_set(igraph_set_t *set, FILE *f) {
    long int state = 0;
    igraph_integer_t element;
    while (igraph_set_iterate(set, &state, &element)) {
        fprintf(f, " %li", (long int) element);
    }
    fprintf(f, "\n");
}

int main() {

    igraph_set_t set;
    int i;

    /* simple init */
    igraph_set_init(&set, 0);
    igraph_set_destroy(&set);

    /* addition, igraph_set_size */
    igraph_set_init(&set, 10);
    i = 10;
    while (igraph_set_size(&set) < 10) {
        igraph_set_add(&set, 2 * i);
        i--;
    }
    while (igraph_set_size(&set) < 21) {
        igraph_set_add(&set, 2 * i + 1);
        i++;
    }
    print_set(&set, stdout);

    /* adding existing element */
    igraph_set_add(&set, 8);
    if (igraph_set_size(&set) != 21) {
        return 4;
    }

    /* igraph_set_contains */
    if (igraph_set_contains(&set, 42) || !igraph_set_contains(&set, 7)) {
        return 3;
    }

    /* igraph_set_empty, igraph_set_clear */
    if (igraph_set_empty(&set)) {
        return 1;
    }
    igraph_set_clear(&set);
    if (!igraph_set_empty(&set)) {
        return 2;
    }
    igraph_set_destroy(&set);

    VERIFY_FINALLY_STACK();

    return 0;
}

