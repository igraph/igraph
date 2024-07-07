/*
   IGraph library.
   Copyright (C) 2021-2022  The igraph development team

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
#include <cstdlib>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    igraph_t graph;
    igraph_vector_int_t edges;

    igraph_set_warning_handler(igraph_warning_handler_ignore);

    if (Size % 2 == 1 || Size > 65280) {
        return 0;
    }

    igraph_vector_int_init(&edges, Size);
    for (size_t i=0; i < Size; ++i) {
        VECTOR(edges)[i] = Data[i];
    }

    if (! igraph_create(&graph, &edges, 0, IGRAPH_DIRECTED)) {
        igraph_integer_t conn;

        /* Enable connectivity checks in order to try to force the fuzzer
         * to find connected graphs. Disconnected graphs result in low coverage. */
        igraph_edge_connectivity(&graph, &conn, 1);
        igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_COLLAPSE, nullptr);
        igraph_edge_connectivity(&graph, &conn, 1);

        igraph_destroy(&graph);
    }

    igraph_vector_int_destroy(&edges);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;  // Non-zero return values are reserved for future use.
}
