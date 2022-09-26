/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>
#include "test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_real_t cent;
    igraph_arpack_options_t arpack_options;

    igraph_star(&g, 10, IGRAPH_STAR_IN, /*center=*/ 0);
    igraph_centralization_degree(&g,
                                 /*res=*/ NULL,
                                 /*mode=*/ IGRAPH_IN, IGRAPH_NO_LOOPS,
                                 &cent, /*theoretical_max=*/ NULL,
                                 /*normalized=*/ true);
    printf("in-star, degree: %g\n", cent);

    igraph_centralization_betweenness(&g, /*res=*/ NULL,
                                      IGRAPH_UNDIRECTED, &cent,
                                      /*theoretical_max=*/ NULL,
                                      /*normalized=*/ true);
    printf("in-star, betweenness: %g\n", cent);
    igraph_destroy(&g);


    igraph_star(&g, 10, IGRAPH_STAR_OUT, /*center=*/ 0);
    igraph_centralization_degree(&g, /*res=*/ NULL,
                                 /*mode=*/ IGRAPH_OUT, IGRAPH_NO_LOOPS,
                                 &cent, /*theoretical_max=*/ NULL,
                                 /*normalized=*/ true);
    printf("out-star, degree: %g\n", cent);

    igraph_centralization_betweenness(&g, /*res=*/ NULL,
                                      IGRAPH_UNDIRECTED, &cent,
                                      /*theoretical_max=*/ NULL,
                                      /*normalized=*/ true);
    printf("out-star, betweenness: %g\n", cent);
    igraph_destroy(&g);

    igraph_small(&g, /*n=*/ 10, /*directed=*/ 0, 0, 1, -1);
    igraph_arpack_options_init(&arpack_options);
    igraph_centralization_eigenvector_centrality(
                &g,
                /*vector=*/ NULL,
                /*value=*/ NULL,
                IGRAPH_DIRECTED,
                /*scale=*/ false,
                &arpack_options, &cent,
                /*theoretical_max=*/ NULL,
                /*normalized=*/ true);

    printf("dyad, eigenvector centrality, not scaled: %g\n", cent);

    igraph_destroy(&g);
    VERIFY_FINALLY_STACK();
    return 0;
}
