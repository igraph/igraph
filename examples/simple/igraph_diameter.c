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
    igraph_real_t result;
    igraph_int_t from, to;
    igraph_vector_int_t path, path_edge;

    /* Initialize the library. */
    igraph_setup();

    igraph_barabasi_game(&graph, 30, /*power=*/ 1, 30, 0, 0, /*A=*/ 1,
                         IGRAPH_DIRECTED, IGRAPH_BARABASI_BAG,
                         /*start_from=*/ NULL);
    igraph_diameter(&graph, NULL, &result, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, 1);

    /* printf("diameter: %g\n", result); */

    igraph_destroy(&graph);

    igraph_ring(&graph, 10, IGRAPH_DIRECTED, 0, 0);
    igraph_vector_int_init(&path, 0);
    igraph_vector_int_init(&path_edge, 0);
    igraph_diameter(&graph, NULL, &result, &from, &to, &path, &path_edge, IGRAPH_DIRECTED, 1);
    printf(
        "diameter: %g, from %" IGRAPH_PRId " to %" IGRAPH_PRId "\n",
        result, from, to
    );
    igraph_vector_int_print(&path);
    igraph_vector_int_print(&path_edge);

    igraph_vector_int_destroy(&path);
    igraph_vector_int_destroy(&path_edge);
    igraph_destroy(&graph);

    return 0;
}
