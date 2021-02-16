/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.inc"

int main() {
    igraph_t g;
    igraph_vector_int_t res, res_all;
    long int i;

    igraph_small(&g, 6, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 5,
                 0, 3, 3, 4, 4, 5,
                 3, 2, 3, 5,
                 -1);

    igraph_vector_int_init(&res, 0);    

    for (i = 0; i <= 5; i++) {
        igraph_get_all_simple_paths(&g, &res, 0, igraph_vss_1(5), i, IGRAPH_ALL);

        printf("Paths for cutoff %li:\n", i);
        igraph_vector_int_print(&res);
    }

    igraph_vector_int_init(&res_all, 0);

    igraph_get_all_simple_paths(&g, &res_all, 0, igraph_vss_1(5), -1, IGRAPH_ALL);

    IGRAPH_ASSERT(igraph_vector_int_all_e(&res, &res_all) && "Paths of all lengths does not equal result for maximum cutoff.");

    igraph_vector_int_destroy(&res_all);
    igraph_vector_int_destroy(&res);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
