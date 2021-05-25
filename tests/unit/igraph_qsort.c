/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA 02139, USA

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

#include "test_utilities.inc"

int comp(const void *a, const void *b) {
    igraph_real_t *aa = (igraph_real_t *) a;
    igraph_real_t *bb = (igraph_real_t *) b;

    if (*aa < *bb) {
        return -1;
    } else if (*aa > *bb) {
        return 1;
    }

    return 0;
}

int main() {
    const int len = 100;
    igraph_vector_t v;
    int i;

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_vector_init(&v, len);
    for (i = 0; i < len; i++) {
        VECTOR(v)[i] = i;
    }
    igraph_vector_shuffle(&v);

    igraph_qsort(VECTOR(v), igraph_vector_size(&v), sizeof(VECTOR(v)[0]), comp);

    igraph_vector_print(&v);

    igraph_vector_destroy(&v);

    VERIFY_FINALLY_STACK();

    return 0;
}
