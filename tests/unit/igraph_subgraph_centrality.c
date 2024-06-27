/*
   IGraph library.
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

#include "test_utilities.h"

int main(void) {
    igraph_t graph;
    igraph_vector_t sc, sc2;

    igraph_vector_init(&sc, 0);
    igraph_vector_init(&sc2, 0);

    igraph_famous(&graph, "frucht");
    igraph_subgraph_centrality(&graph, &sc, true);
    print_vector(&sc);

    igraph_destroy(&graph);

    /* Test that edge directions are ignored, and that self-loops are handled correctly. */

    igraph_de_bruijn(&graph, 3, 2);
    igraph_subgraph_centrality(&graph, &sc, true);
    print_vector(&sc);

    igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_EACH, NULL);
    igraph_subgraph_centrality(&graph, &sc2, true);

    IGRAPH_ASSERT(igraph_vector_all_e(&sc, &sc2));

    igraph_destroy(&graph);

    igraph_vector_destroy(&sc2);
    igraph_vector_destroy(&sc);

    return 0;
}
