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
#include "test_utilities.inc"

int main() {
    igraph_t g;
    long int n_vertices = 10;

    printf("Undirected, no loops\n");
    igraph_full(&g, n_vertices, 0 /*undirected*/, 0/*no loops*/);
    print_graph(&g, stdout);
    igraph_destroy(&g);

    printf("Directed, no loops\n");
    igraph_full(&g, n_vertices, 1 /*directed*/, 0/*no loops*/);
    print_graph(&g, stdout);
    igraph_destroy(&g);

    printf("Undirected, with loops\n");
    igraph_full(&g, n_vertices, 0 /*undirected*/, 1/*loops*/);
    print_graph(&g, stdout);
    igraph_destroy(&g);

    printf("Directed, with loops\n");
    igraph_full(&g, n_vertices, 1 /*directed*/, 1/*loops*/);
    print_graph(&g, stdout);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;

}
