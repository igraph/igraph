/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
#include <math.h>

int main(void) {

    igraph_t g;
    igraph_real_t cent;
    igraph_arpack_options_t arpack_options;

    /****************************/
    /* undirected star */
    igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, /*center=*/ 0);

    igraph_centralization_degree(&g, /*res=*/ NULL,
                                 /*mode=*/ IGRAPH_ALL, IGRAPH_NO_LOOPS,
                                 &cent, /*theoretical_max=*/ NULL,
                                 /*normalized=*/ true);
    if (cent != 1.0) {
        fprintf(stderr, "undirected star, degree: %g\n", cent);
        return 21;
    }

    igraph_centralization_betweenness(&g, /*res=*/ NULL,
                                      IGRAPH_UNDIRECTED, &cent,
                                      /*theoretical_max=*/ NULL,
                                      /*normalized=*/ true);
    if (cent != 1.0) {
        fprintf(stderr, "undirected star, betweenness: %g\n", cent);
        return 22;
    }

    igraph_centralization_closeness(&g, /*res=*/ NULL,
                                    IGRAPH_ALL, &cent,
                                    /*theoretical_max=*/ NULL,
                                    /*normalized=*/ true);

    if (!igraph_almost_equals(cent, 1.0, 1e-8)) {
        fprintf(stderr, "undirected star, closeness: %g\n", cent);
        return 23;
    }

    igraph_destroy(&g);

    /****************************/
    /* single dyad */

    igraph_small(&g, /*n=*/ 10, /*directed=*/ 0,
                 0, 1, -1);

    igraph_arpack_options_init(&arpack_options);
    igraph_centralization_eigenvector_centrality(
                &g,
                /*vector=*/ NULL,
                /*value=*/ NULL,
                IGRAPH_DIRECTED,
                /*scale=*/ true,
                &arpack_options, &cent,
                /*theoretical_max=*/ NULL,
                /*normalized=*/ true);

    if (!igraph_almost_equals(cent, 1.0, 1e-8)) {
        fprintf(stderr, "dyad, eigenvector centrality: %g\n", cent);
        return 24;
    }

    igraph_destroy(&g);

    return 0;
}
