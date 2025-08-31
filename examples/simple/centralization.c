/*
   igraph library.
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

int main(void) {

    igraph_t graph;
    igraph_real_t cent;

    /* Initialize the library. */
    igraph_setup();

    /* Create an undirected star graph, which is the most centralized graph
     * with several common centrality scores. */
    printf("undirected star graph:\n");
    igraph_star(&graph, 10, IGRAPH_STAR_UNDIRECTED, /*center=*/ 0);

    igraph_centralization_degree(&graph, /*res=*/ NULL,
                                 /*mode=*/ IGRAPH_ALL, IGRAPH_NO_LOOPS,
                                 &cent, /*theoretical_max=*/ NULL,
                                 /*normalized=*/ true);
    printf("degree centralization: %g\n", cent);

    igraph_centralization_betweenness(&graph, /*res=*/ NULL,
                                      IGRAPH_UNDIRECTED, &cent,
                                      /*theoretical_max=*/ NULL,
                                      /*normalized=*/ true);
    printf("betweenness centralization: %g\n", cent);


    igraph_centralization_closeness(&graph, /*res=*/ NULL,
                                    IGRAPH_ALL, &cent,
                                    /*theoretical_max=*/ NULL,
                                    /*normalized=*/ true);
    printf("closeness centralization: %g\n", cent);

    igraph_destroy(&graph);

    /* With eigenvector centrality, the most centralized structure is
     * a graph containing a single edge. */
    printf("\ngraph with single edge:\n");
    igraph_small(&graph, /*n=*/ 10, IGRAPH_UNDIRECTED,
                 0,1, -1);

    igraph_centralization_eigenvector_centrality(
                &graph,
                /*vector=*/ NULL,
                /*value=*/ NULL,
                /*mode=*/ IGRAPH_ALL,
                /*options=*/ NULL,
                &cent,
                /*theoretical_max=*/ NULL,
                /*normalized=*/ true);
    printf("eigenvector centralization: %g\n", cent);

    igraph_destroy(&graph);

    return 0;
}
