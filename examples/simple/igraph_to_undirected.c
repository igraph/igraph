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

    igraph_vector_t v;
    igraph_t g;

    igraph_vector_init_int(&v, 2, 5, 5);
    igraph_lattice(&g, &v, 1, IGRAPH_DIRECTED, 1 /*mutual*/, 0 /*circular*/);
    igraph_to_undirected(&g, IGRAPH_TO_UNDIRECTED_COLLAPSE,
                         /*edge_comb=*/ 0);
    igraph_write_graph_edgelist(&g, stdout);

    igraph_destroy(&g);
    igraph_vector_destroy(&v);

    printf("---\n");

    igraph_small(&g, 10, IGRAPH_DIRECTED,
                 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 3,
                 5, 6, 6, 5, 6, 7, 6, 7, 7, 6, 7, 8, 7, 8, 8, 7, 8, 7, 8, 8, 9, 9, 9, 9,
                 -1);
    igraph_to_undirected(&g, IGRAPH_TO_UNDIRECTED_MUTUAL,
                         /*edge_comb=*/ 0);
    igraph_write_graph_edgelist(&g, stdout);
    igraph_destroy(&g);

    return 0;
}
