/*
   igraph library.
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
#include "test_utilities.h"

int main(void) {
    igraph_t g_in;
    igraph_t g_out;
    FILE *file;
    igraph_bool_t same;

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_small(&g_in, 4, IGRAPH_DIRECTED,
            0, 1, 0, 1, 0, 1,
            1, 2, 2, 3, 3, 0,
            -1);

    SETEAN(&g_in, "edge_attr", 0, 12.3); /* edge 0-1 */
    SETEAN(&g_in, "edge_attr", 4, -IGRAPH_INFINITY); /* edge 2-3 */
    SETVAS(&g_in, "vertex_attr", 0, "vertex_name0");
    SETVAS(&g_in, "vertex_attr", 1, "vertex_name1");
    SETVAS(&g_in, "vertex_attr", 2, "vertex_name2");
    SETVAS(&g_in, "vertex_attr", 3, "vertex_name3");

    char filename[] = "ncol.tmp";
    file = fopen(filename, "w");
    IGRAPH_ASSERT(file != NULL); /* make sure that the file was created successfully */

    igraph_write_graph_ncol(&g_in, file, "vertex_attr", "edge_attr");
    fclose(file);


    file = fopen(filename, "r");
    IGRAPH_ASSERT(file);
    igraph_read_graph_ncol(&g_out, file, NULL, 1, IGRAPH_ADD_WEIGHTS_YES,
            IGRAPH_DIRECTED);

    /* is_same_graph() checks that vertex and edge lists are the same,
     * but does not verify attributes. */
    igraph_is_same_graph(&g_in, &g_out, &same);
    if (!same) {
        IGRAPH_FATAL("Written and read graph are not the same.");
    }

    print_attributes(&g_out);

    fclose(file);
    unlink(filename);
    igraph_destroy(&g_in);
    igraph_destroy(&g_out);

    VERIFY_FINALLY_STACK();

    return 0;
}
