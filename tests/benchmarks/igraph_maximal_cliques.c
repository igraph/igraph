/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "bench.h"

void free_result(igraph_vector_ptr_t *res) {
    long int i, n;

    n = igraph_vector_ptr_size(res);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(*res)[i];
        igraph_vector_destroy(v);
        igraph_free(v);
    }

    igraph_vector_ptr_resize(res, 0);
}

int main() {

    igraph_t g;
    igraph_real_t toremovev[] = {  2609,  2098, 14517,  7540, 19560,  8855,
                                   5939, 14947,   441, 16976, 19642,  4188,
                                   15447, 11837,  2333,  7309, 18539, 14099,
                                   14264,  9240
                                };
    igraph_vector_t toremove;
    igraph_vector_ptr_t res;

    BENCH_INIT();

    igraph_vector_view(&toremove, toremovev,
                       sizeof(toremovev) / sizeof(igraph_real_t));
    igraph_full(&g, 200, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_delete_edges(&g, igraph_ess_vector(&toremove));

    igraph_vector_ptr_init(&res, 0);

    BENCH(" 1 Maximal cliques of almost complete graph",
          igraph_maximal_cliques(&g, &res, /* min_size= */ 0,
                                 /* max_size= */ 0);
         );

    igraph_destroy(&g);

    free_result(&res);
    igraph_vector_ptr_destroy(&res);

    return 0;
}
