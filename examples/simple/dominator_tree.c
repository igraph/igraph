/*
   igraph library.
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
#include <stdio.h>

int main(void) {
    igraph_t g, domtree;
    igraph_vector_int_t dom;
    igraph_vector_int_t leftout;

    /* Initialize the library. */
    igraph_setup();

    igraph_vector_int_init(&dom, 0);
    igraph_vector_int_init(&leftout, 0);

    igraph_small(&g, 10, IGRAPH_DIRECTED,
                 0, 9,
                 1, 0, 1, 2,
                 2, 3, 2, 7,
                 3, 1,
                 4, 1, 4, 3,
                 5, 2, 5, 3, 5, 4, 5, 8,
                 6, 5, 6, 9,
                 8, 7,
                 -1);

    igraph_dominator_tree(&g, /*root=*/ 9, &dom, &domtree,
                          &leftout, /*mode=*/ IGRAPH_IN);
    igraph_vector_int_print(&dom);
    igraph_vector_int_print(&leftout);
    igraph_write_graph_edgelist(&domtree, stdout);

    igraph_vector_int_destroy(&dom);
    igraph_vector_int_destroy(&leftout);
    igraph_destroy(&domtree);
    igraph_destroy(&g);

    return 0;
}
