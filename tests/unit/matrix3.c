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
    igraph_matrix_t m;

    igraph_matrix_init(&m, 10, 10);
    if (igraph_matrix_capacity(&m) != 100) {
        return 1;
    }

    igraph_matrix_add_cols(&m, 5);
    igraph_matrix_resize(&m, 5, 5);
    igraph_matrix_resize_min(&m);
    if (igraph_matrix_capacity(&m) != igraph_matrix_size(&m)) {
        return 2;
    }

    igraph_matrix_destroy(&m);

    VERIFY_FINALLY_STACK();

    return 0;
}
