/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_t graph;
    igraph_vector_int_t component_sizes;

    /* Initialize the library. */
    igraph_setup();

    igraph_rng_seed(igraph_rng_default(), 42); /* make program deterministic */

    /* Sample a graph from the Erdős-Rényi G(n,p) model */

    igraph_erdos_renyi_game_gnp(
            &graph, /* n= */ 100, /* p= */ 0.01,
            IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    /* Compute the fraction of vertices contained within the largest connected component */

    igraph_vector_int_init(&component_sizes, 0);
    igraph_connected_components(&graph, NULL, &component_sizes, NULL, IGRAPH_STRONG);

    printf(
        "Fraction of vertices in giant component: %g\n",
        ((double) igraph_vector_int_max(&component_sizes)) / igraph_vcount(&graph)
    );

    /* Clean up data structures when no longer needed */

    igraph_vector_int_destroy(&component_sizes);
    igraph_destroy(&graph);

    return 0;
}
