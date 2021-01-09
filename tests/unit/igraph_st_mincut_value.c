/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

int main() {

    igraph_t g;
    igraph_vector_t capacity;
    igraph_real_t value;

    igraph_small(&g, 6, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 1, 2, 1, 3, 2, 4, 3, 4, 3, 5, 4, 5, -1);
    igraph_vector_init_int_end(&capacity, -1, 5, 2, 2, 3, 4, 1, 2, 5, -1);

    igraph_st_mincut_value(&g, &value, 0, 5, &capacity);

    igraph_vector_destroy(&capacity);
    igraph_destroy(&g);

    if (value == 7) {
        return 0;
    } else {
        return 1;
    }
}
