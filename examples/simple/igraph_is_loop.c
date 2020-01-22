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

void print_vector(igraph_vector_bool_t *v, FILE *f) {
    long int i;
    for (i = 0; i < igraph_vector_bool_size(v); i++) {
        fprintf(f, " %i", (int) VECTOR(*v)[i]);
    }
    fprintf(f, "\n");
}

int main() {

    igraph_t graph;
    igraph_vector_bool_t v;

    igraph_vector_bool_init(&v, 0);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 1, 0, 1, 1, 0, 3, 4, 11, 10, -1);
    igraph_is_loop(&graph, &v, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    print_vector(&v, stdout);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 0, 0, 1, 1, 2, 2, 2, 3, 2, 4, 2, 5, 2, 6, 2, 2, 0, 0, -1);
    igraph_is_loop(&graph, &v, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    print_vector(&v, stdout);
    igraph_destroy(&graph);

    igraph_vector_bool_destroy(&v);

    return 0;
}
