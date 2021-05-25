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

#include "test_utilities.inc"

int main() {
    igraph_t graph;
    igraph_vector_bool_t type;
    igraph_bool_t typev[] = { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1 };

    /* turn on attribute handling */
    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_ring(&graph, 10, IGRAPH_UNDIRECTED, /*mutual=*/ 0, /*circular=*/ 1);
    igraph_vector_bool_view(&type, typev, sizeof(typev) / sizeof(igraph_bool_t));
    SETVABV(&graph, "type", &type);

    igraph_write_graph_pajek(&graph, stdout);

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
