/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, USA 02139

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
    igraph_vector_t v;

    igraph_vector_init_seq(&v, 1, 1000);
    IGRAPH_ASSERT(igraph_vector_capacity(&v) == 1000);

    igraph_vector_push_back(&v, 1001);
    IGRAPH_ASSERT(igraph_vector_capacity(&v) == 2000);

    igraph_vector_resize_min(&v);
    IGRAPH_ASSERT(igraph_vector_capacity(&v) == igraph_vector_size(&v));

    igraph_vector_destroy(&v);

    /* regression test for #1479 -- calling resize_min() on an empty vector */
    igraph_vector_init_seq(&v, 1, 1000);
    igraph_vector_clear(&v);
    igraph_vector_resize_min(&v);
    IGRAPH_ASSERT(igraph_vector_capacity(&v) == 0);
    IGRAPH_ASSERT(igraph_vector_size(&v) == 0);
    igraph_vector_destroy(&v);

    VERIFY_FINALLY_STACK();

    return 0;
}
