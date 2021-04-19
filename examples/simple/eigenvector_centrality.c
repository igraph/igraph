/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph.h"

#include <math.h>

int main() {

    igraph_t g;
    igraph_vector_t v, weights;
    long int i;
    igraph_real_t value;
    igraph_arpack_options_t options;

    igraph_star(&g, 100, IGRAPH_STAR_UNDIRECTED, 0);

    igraph_arpack_options_init(&options);
    igraph_vector_init(&v, 0);
    igraph_eigenvector_centrality(&g, &v, &value, /*directed=*/ 0,
                                  /*scale=*/1, /*weights=*/0,
                                  &options);

    if (options.info != 0) {
        return 1;
    }

    for (i = 0; i < igraph_vector_size(&v); i++) {
        printf(" %.4f", fabs(VECTOR(v)[i]));
    }
    printf("\n");

    igraph_destroy(&g);

    /* Special cases: check for empty graph */
    igraph_empty(&g, 10, 0);
    igraph_eigenvector_centrality(&g, &v, &value, 0, 0, 0, &options);
    if (value != 0.0) {
        return 1;
    }
    for (i = 0; i < igraph_vector_size(&v); i++) {
        printf(" %.2f", fabs(VECTOR(v)[i]));
    }
    printf("\n");
    igraph_destroy(&g);

    /* Special cases: check for full graph, zero weights */
    igraph_full(&g, 10, 0, 0);
    igraph_vector_init(&weights, 45);
    igraph_vector_fill(&weights, 0);
    igraph_eigenvector_centrality(&g, &v, &value, 0, 0, &weights, &options);
    igraph_vector_destroy(&weights);
    if (value != 0.0) {
        return 2;
    }
    for (i = 0; i < igraph_vector_size(&v); i++) {
        printf(" %.2f", fabs(VECTOR(v)[i]));
    }
    printf("\n");
    igraph_destroy(&g);

    igraph_vector_destroy(&v);

    return 0;
}
