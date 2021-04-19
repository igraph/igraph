/* -*- mode: C -*-  */
/*
   IGraph library.
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

int main() {

    igraph_t g;
    int ret;

    igraph_atlas(&g, 45);
    igraph_write_graph_edgelist(&g, stdout);
    printf("\n");
    igraph_destroy(&g);

    igraph_atlas(&g, 0);
    igraph_write_graph_edgelist(&g, stdout);
    printf("\n");
    igraph_destroy(&g);

    igraph_atlas(&g, 1252);
    igraph_write_graph_edgelist(&g, stdout);
    printf("\n");
    igraph_destroy(&g);

    igraph_set_error_handler(igraph_error_handler_ignore);
    ret = igraph_atlas(&g, -1);
    if (ret != IGRAPH_EINVAL) {
        return 1;
    }

    ret = igraph_atlas(&g, 1253);
    if (ret != IGRAPH_EINVAL) {
        return 2;
    }

    return 0;
}
