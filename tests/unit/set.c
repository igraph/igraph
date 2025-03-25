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

#include "test_utilities.h"

void print_set(igraph_set_t *set, FILE *f) {
    igraph_integer_t state = 0;
    igraph_integer_t element;
    while (igraph_set_iterate(set, &state, &element)) {
        fprintf(f, " %" IGRAPH_PRId , element);
    }
    fprintf(f, "\n");
}

int main(void) {

    igraph_set_t set;
    igraph_integer_t i;

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
    
    /* igraph_set_remove */
    // remove element in the middle
    igraph_set_remove(&set, 10);
    if (igraph_set_size(&set) != 20 || igraph_set_contains(&set, 10)) {
        return 5;
    }
    // remove element not in the set - above, below, middle
    igraph_set_remove(&set, 25);
    igraph_set_remove(&set, 10);
    igraph_set_remove(&set, -25);
    if (igraph_set_size(&set) != 20) {
        return 6;
    }
    // remove last element
    igraph_set_remove(&set, 21);
    if (igraph_set_size(&set) != 19 || igraph_set_contains(&set, 21)) {
        return 7;
    }
    // remove first element
    igraph_set_remove(&set, 1);
    if (igraph_set_size(&set) != 18 || igraph_set_contains(&set, 1)) {
        return 8;
    }
    // remove elements around gap
    igraph_set_remove(&set, 9);
    igraph_set_remove(&set, 11);
    if (igraph_set_size(&set) != 16 || igraph_set_contains(&set, 9) || igraph_set_contains(&set, 11)) {
        return 9;
    }

    /* igraph_set_empty, igraph_set_clear */
    if (igraph_set_empty(&set)) {
        return 1;
    }
    igraph_set_clear(&set);
    if (!igraph_set_empty(&set)) {
        return 2;
    }

    /* igraph_set_empty for small sets */
    igraph_set_remove(&set, 10);
    if (igraph_set_size(&set) != 0) { return 12; }
    igraph_set_add(&set, 10);
    igraph_set_remove(&set, 10);
    if (igraph_set_size(&set) != 0) { return 13; }
    igraph_set_destroy(&set);

    VERIFY_FINALLY_STACK();

    return 0;
}
