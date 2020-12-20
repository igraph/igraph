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

void igraph_warnings_ignore(const char *reason, const char *file,
                            int line, int igraph_errno) {
    /* Do nothing */
}

int main() {
    igraph_t g;
    igraph_vector_t vpath, epath;
    igraph_vector_t w;

    /* Unweighted */

    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4, 0, 3,
                 -1);
    igraph_vector_init(&vpath, 0);
    igraph_vector_init(&epath, 0);
    igraph_get_shortest_path(&g, &vpath, &epath, 0, 4, IGRAPH_OUT);
    igraph_vector_print(&vpath);
    igraph_vector_print(&epath);

    igraph_get_shortest_path(&g, &vpath, &epath, 0, 0, IGRAPH_OUT);
    igraph_vector_print(&vpath);
    igraph_vector_print(&epath);

    igraph_set_warning_handler(igraph_warnings_ignore);
    igraph_get_shortest_path(&g, &vpath, &epath, 4, 0, IGRAPH_OUT);
    igraph_vector_print(&vpath);
    igraph_vector_print(&epath);
    igraph_set_warning_handler(igraph_warning_handler_print);

    igraph_get_shortest_path(&g, &vpath, &epath, 4, 0, IGRAPH_ALL);
    igraph_vector_print(&vpath);
    igraph_vector_print(&epath);

    /* Weighted */

    igraph_vector_init(&w, 5);
    VECTOR(w)[0] = 1;
    VECTOR(w)[1] = 1;
    VECTOR(w)[2] = 1;
    VECTOR(w)[3] = 1;
    VECTOR(w)[4] = 3.1;

    igraph_get_shortest_path_dijkstra(&g, &vpath, &epath, 0, 4, &w, IGRAPH_OUT);
    igraph_vector_print(&vpath);
    igraph_vector_print(&epath);

    igraph_vector_destroy(&w);
    igraph_vector_destroy(&epath);
    igraph_vector_destroy(&vpath);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}

