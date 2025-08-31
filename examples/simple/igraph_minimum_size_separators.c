/*
   igraph library.
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

int main(void) {
    igraph_setup();
    igraph_t g;
    igraph_vector_int_list_t sep;

    igraph_small(&g, 7, IGRAPH_UNDIRECTED,
                 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0,
                 -1);
    igraph_vector_int_list_init(&sep, 0);
    igraph_minimum_size_separators(&g, &sep);

    for (igraph_int_t i = 0; i < igraph_vector_int_list_size(&sep); i++) {
        igraph_vector_int_t* v = igraph_vector_int_list_get_ptr(&sep, i);
        igraph_vector_int_print(v);
    }

    igraph_vector_int_list_destroy(&sep);
    igraph_destroy(&g);

    return 0;
}
