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
    igraph_t g;
    igraph_vector_int_list_t cliques;
    igraph_int_t no;

    /* Initialize the library. */
    igraph_setup();

    igraph_small(&g, 9, IGRAPH_UNDIRECTED,
                 0,1, 0,2, 1,2, 2,3, 3,4, 3, 5, 4,5, 5,6, 6,7, -1);
    igraph_vector_int_list_init(&cliques, 0);
    igraph_maximal_cliques(&g, &cliques, /*min_size=*/ IGRAPH_UNLIMITED,
                           /*max_size=*/ IGRAPH_UNLIMITED,
                           /*max_results=*/ IGRAPH_UNLIMITED);
    igraph_maximal_cliques_count(&g, &no, /*min_size=*/ IGRAPH_UNLIMITED,
                                 /*max_size=*/ IGRAPH_UNLIMITED);

    if (no != igraph_vector_int_list_size(&cliques)) {
        return 1;
    }

    for (igraph_int_t i = 0; i < igraph_vector_int_list_size(&cliques); i++) {
        igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(&cliques, i);
        igraph_vector_int_print(v);
    }

    igraph_vector_int_list_destroy(&cliques);
    igraph_destroy(&g);

    return 0;
}
