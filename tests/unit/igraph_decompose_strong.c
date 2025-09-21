/*
   igraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdlib.h>

#include "test_utilities.h"

int main(void) {

    igraph_t ring, g;
    igraph_graph_list_t complist;
    igraph_int_t i;

    igraph_graph_list_init(&complist, 0);

    /* A directed ring, a single strongly connected component */
    igraph_ring(&ring, 10, IGRAPH_DIRECTED, 0, 1);
    igraph_decompose(&ring, &complist, IGRAPH_STRONG, -1, 0);
    print_graph(igraph_graph_list_get_ptr(&complist, 0));
    igraph_destroy(&ring);

    /* a toy graph, three components maximum, with at least 2 vertices each */
    /* 0 >->  1 >-> 3 >-> 4
       ^      v
       \< 2 < /           */
    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 0,
                 1, 3,
                 3, 4,
                 -1);
    igraph_decompose(&g, &complist, IGRAPH_STRONG, 3, 2);
    for (i = 0; i < igraph_graph_list_size(&complist); i++) {
        print_graph(igraph_graph_list_get_ptr(&complist, i));
    }
    igraph_destroy(&g);

    igraph_graph_list_destroy(&complist);

    VERIFY_FINALLY_STACK();

    return 0;
}
