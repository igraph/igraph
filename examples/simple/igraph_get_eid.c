/*
   igraph library.
   Copyright (C) 2006-2024  The igraph development team <igraph@igraph.org>

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
    igraph_int_t eid;
    igraph_vector_int_t hist;

    /* Initialize the library. */
    igraph_setup();

    /* DIRECTED */

    igraph_star(&graph, 10, IGRAPH_STAR_OUT, 0);

    igraph_vector_int_init(&hist, 9);

    for (igraph_int_t i = 1; i < 10; i++) {
        igraph_get_eid(&graph, &eid, 0, i, IGRAPH_DIRECTED, /*error=*/ true);
        VECTOR(hist)[eid] = 1;
    }

    igraph_vector_int_print(&hist);

    igraph_vector_int_destroy(&hist);
    igraph_destroy(&graph);

    /* UNDIRECTED */

    igraph_star(&graph, 10, IGRAPH_STAR_UNDIRECTED, 0);

    igraph_vector_int_init(&hist, 9);

    for (igraph_int_t i = 1; i < 10; i++) {
        igraph_get_eid(&graph, &eid, 0, i, IGRAPH_UNDIRECTED, /*error=*/ true);
        VECTOR(hist)[eid] += 1;
        igraph_get_eid(&graph, &eid, i, 0, IGRAPH_DIRECTED, /*error=*/ true);
        VECTOR(hist)[eid] += 1;
    }
    igraph_vector_int_print(&hist);

    igraph_vector_int_destroy(&hist);
    igraph_destroy(&graph);

    return 0;
}
