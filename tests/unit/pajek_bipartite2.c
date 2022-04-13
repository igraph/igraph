/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.h"

int main() {
    igraph_t graph;
    FILE *input;

    /* turn on attribute handling */
    igraph_set_attribute_table(&igraph_cattribute_table);

    /* first file, without marginals */

    input = fopen("pajek_bip.net", "r");
    if (input == 0) {
        return 1;
    }

    igraph_read_graph_pajek(&graph, input);
    fclose(input);

    print_attributes(&graph);

    igraph_destroy(&graph);

    /* second file, with marginals */

    printf("---\n");

    input = fopen("pajek_bip2.net", "r");
    if (input == 0) {
        return 1;
    }

    igraph_read_graph_pajek(&graph, input);
    fclose(input);

    print_attributes(&graph);

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
