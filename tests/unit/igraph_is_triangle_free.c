/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

void check_triangle_free(const char *name, const igraph_t *graph) {
    igraph_real_t count;
    igraph_bool_t triangle_free;

    igraph_is_triangle_free(graph, &triangle_free);
    igraph_count_triangles(graph, &count);

    printf("%s: %s\n", name, triangle_free ? "yes": " no");

    IGRAPH_ASSERT((!!count) == (!triangle_free));
}

int main(void) {
    igraph_t graph;
    igraph_bool_t simple;

    igraph_setup();

    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    check_triangle_free("Null graph", &graph);
    igraph_destroy(&graph);

    igraph_full(&graph, 1, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    check_triangle_free("Singleton with loop", &graph);
    igraph_destroy(&graph);

    igraph_cycle_graph(&graph, 3, IGRAPH_UNDIRECTED, false);
    check_triangle_free("C_3", &graph);
    igraph_destroy(&graph);

    igraph_cycle_graph(&graph, 4, IGRAPH_UNDIRECTED, false);
    check_triangle_free("C_4", &graph);
    igraph_destroy(&graph);

    igraph_full(&graph, 4, IGRAPH_UNDIRECTED, false);
    check_triangle_free("K4", &graph);
    igraph_is_simple(&graph, &simple, IGRAPH_UNDIRECTED); /* populate cache */
    check_triangle_free("K4", &graph);
    igraph_destroy(&graph);

    igraph_full(&graph, 4, IGRAPH_DIRECTED, false);
    check_triangle_free("Directed K4", &graph);
    igraph_is_simple(&graph, &simple, IGRAPH_UNDIRECTED); /* populate cache */
    check_triangle_free("Directed K4", &graph);
    igraph_destroy(&graph);

    igraph_small(&graph, 4, IGRAPH_UNDIRECTED,
        0,1, 1,2, 2,3, 3,0, 0,2,
        -1);
    check_triangle_free("Extremal + 1 on 4 vertices", &graph);
    igraph_is_simple(&graph, &simple, IGRAPH_UNDIRECTED); /* populate cache */
    check_triangle_free("Extremal + 1 on 4 vertices", &graph);
    igraph_destroy(&graph);

    igraph_small(&graph, 4, IGRAPH_DIRECTED,
        0,1, 1,2, 2,3, 3,0, 0,3,
        -1);
    check_triangle_free("Directed with mutual, no triangles, high edge count", &graph);
    igraph_is_simple(&graph, &simple, IGRAPH_UNDIRECTED); /* populate cache */
    check_triangle_free("Directed with mutual, no triangles, high edge count", &graph);
    igraph_destroy(&graph);

    igraph_small(&graph, 4, IGRAPH_UNDIRECTED,
        0,1, 1,2, 2,3, 3,0, 0,1, 2,2,
        -1);
    check_triangle_free("Non-simple, no triangle, high edge count", &graph);
    igraph_is_simple(&graph, &simple, IGRAPH_UNDIRECTED); /* populate cache */
    check_triangle_free("Non-simple, no triangle, high edge count", &graph);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
