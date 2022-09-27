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

int main(void) {

    igraph_t g;
    igraph_vector_t v, weights;
    igraph_real_t value;
    igraph_arpack_options_t options;

    igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);
    igraph_vector_init(&weights, 9);
    igraph_vector_fill(&weights, 1);

    igraph_arpack_options_init(&options);
    igraph_vector_init(&v, 0);
    igraph_eigenvector_centrality(&g, &v, &value, /*directed=*/ 0,
                                  /*scale=*/1, &weights,
                                  &options);

    if (options.info != 0) {
        return 1;
    }

    printf("eigenvalue: %g\n", value);
    printf("eigenvector:\n");
    igraph_vector_print(&v);

    igraph_destroy(&g);
    igraph_vector_destroy(&v);
    igraph_vector_destroy(&weights);

    return 0;
}
