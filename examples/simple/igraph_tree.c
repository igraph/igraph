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

int main() {
    igraph_t graph;
    igraph_bool_t res;

    /* Create a directed binary tree on 15 nodes,
       with edges pointing towards the root. */
    igraph_tree(&graph, 15, 2, IGRAPH_TREE_IN);

    igraph_is_tree(&graph, &res, NULL, IGRAPH_IN);
    printf("Is it an in-tree? %s\n", res ? "Yes" : "No");

    igraph_is_tree(&graph, &res, NULL, IGRAPH_OUT);
    printf("Is it an out-tree? %s\n", res ? "Yes" : "No");

    igraph_destroy(&graph);

    return 0;
}
