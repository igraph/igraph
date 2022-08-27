/* -*- mode: C -*-  */
/* IGraph library.
   Copyright (C) 2007-2021  The igraph development team <igraph@igraph.org>

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
#include <string.h>
#include <stdlib.h>

int main() {
    igraph_t graph;
    igraph_vector_t y;

    /* Turn on attribute handling. */
    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_small(&graph, 3, IGRAPH_DIRECTED, 0,1, 1,2, -1);

    /* Set graph attributes. */
    /* numeric */
    SETGAN(&graph, "id", 10);
    /* string */
    SETGAS(&graph, "name", "toy");
    /* boolean */
    SETGAB(&graph, "is_regular", false);

    /* Set edge string attribute. */
    SETEAS(&graph, "color", 1, "RED");

    /* Set vertex attributes as vector. */
    igraph_vector_init(&y, igraph_vcount(&graph));
    igraph_vector_fill(&y, 1.23);
    SETVANV(&graph, "y", &y);
    igraph_vector_destroy(&y);

    /* Set single vertex numeric attribute. */
    SETVAN(&graph, "y", 0, -1);

    /* Delete graph attribute. */
    DELGA(&graph, "is_regular");

    /* Print the final result. */
    igraph_print_attributes(&graph);

    /* Delete all remaining attributes. */
    DELALL(&graph);

    /* Destroy the graph. */
    igraph_destroy(&graph);

    return 0;
}
