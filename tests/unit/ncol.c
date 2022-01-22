/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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
#include <unistd.h>
#include "test_utilities.inc"

int main() {
    igraph_t g;
    char filename[] = "ncol.tmpXXXXXX";
    FILE *file;

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_small(&g, 4, IGRAPH_DIRECTED,
            0, 1, 0, 1, 0, 1,
            1, 2, 2, 3, 3, 0,
            -1);

    SETEAN(&g, "edge_attr", 0, 10);
    SETVAS(&g, "vertex_attr", 0, "vertex_name0");
    SETVAS(&g, "vertex_attr", 1, "vertex_name1");
    SETVAS(&g, "vertex_attr", 2, "vertex_name2");
    SETVAS(&g, "vertex_attr", 3, "vertex_name3");

    mktemp(filename);;
    file = fopen(filename, "w");

    igraph_write_graph_ncol(&g, file, "vertex_attr", "edge_attr");
    fclose(file);

    igraph_destroy(&g);

    file = fopen(filename, "r");
    igraph_read_graph_ncol(&g, file, NULL, 1, IGRAPH_ADD_WEIGHTS_YES,
            IGRAPH_DIRECTED);

    igraph_write_graph_graphml(&g, stdout, 0);

    fclose(file);
    unlink(filename);
    igraph_destroy(&g);
    VERIFY_FINALLY_STACK();
    return 0;
}
