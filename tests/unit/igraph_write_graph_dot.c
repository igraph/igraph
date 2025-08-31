/*
   igraph library.
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

    igraph_set_warning_handler(igraph_warning_handler_ignore);
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_small(&g, 6, IGRAPH_UNDIRECTED,
            0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 4,3, 4,3, -1);
    /* Set vertex attributes */
    SETVAN(&g, "VAN", 0, -1);
    SETVAN(&g, "VAN", 1, 2.1);
    SETVAN(&g, "VAN", 2, 1.23e-6);
    SETVAN(&g, "VAN", 3, 1e7);

    SETVAS(&g, "VAS", 0, "foo");
    SETVAS(&g, "VAS", 1, "bar");

    SETVAB(&g, "VAB", 0, 1);
    SETVAB(&g, "VAB", 1, 0);

    /* Set edge attributes */
    SETEAN(&g, "EAN", 0, -100.1);
    SETEAN(&g, "EAN", 2, 100.0);

    SETEAS(&g, "EAS", 0, "Blue");
    SETEAS(&g, "EAS", 2, "RED");

    SETEAB(&g, "EAB", 0, 1);
    SETEAB(&g, "EAB", 2, 0);

    igraph_write_graph_dot(&g, stdout);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
