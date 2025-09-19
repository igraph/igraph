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
#include <math.h>

#include "test_utilities.h"

int main(void) {
    igraph_t graph;
    igraph_real_t k;

    igraph_empty(&graph, 0, IGRAPH_DIRECTED);
    igraph_mean_degree(&graph, &k, true);
    IGRAPH_ASSERT(isnan(k));
    igraph_destroy(&graph);

    igraph_empty(&graph, 10, IGRAPH_UNDIRECTED);
    igraph_mean_degree(&graph, &k, true);
    IGRAPH_ASSERT(k == 0);
    igraph_destroy(&graph);

    igraph_ring(&graph, 5, IGRAPH_DIRECTED, false, true);
    igraph_mean_degree(&graph, &k, true);
    IGRAPH_ASSERT(k == 1);
    igraph_destroy(&graph);

    igraph_ring(&graph, 5, IGRAPH_UNDIRECTED, false, true);
    igraph_mean_degree(&graph, &k, true);
    IGRAPH_ASSERT(k == 2);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 0,1, 1,1, -1);
    igraph_mean_degree(&graph, &k, true);
    IGRAPH_ASSERT(k == 2);
    igraph_mean_degree(&graph, &k, false);
    IGRAPH_ASSERT(k == 1);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
