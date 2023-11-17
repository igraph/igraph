/*
   IGraph library.
   Copyright (C) 2021-2023  The igraph development team

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

#include "igraph.h"
#include <cstdio>

extern "C"
int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {

    igraph_set_error_handler(igraph_error_handler_ignore);
    igraph_set_warning_handler(igraph_warning_handler_ignore);

    // Turn on attribute handling
    igraph_set_attribute_table(&igraph_cattribute_table);

    // Read input file
    FILE *ifile = fmemopen((void*) data, size, "r");
    if (!ifile) {
        return 0;
    }

    igraph_strvector_t problem;
    IGRAPH_ASSERT(igraph_strvector_init(&problem, 0) == IGRAPH_SUCCESS);

    igraph_vector_int_t label;
    IGRAPH_ASSERT(igraph_vector_int_init(&label, 0) == IGRAPH_SUCCESS);

    igraph_vector_t capacity;
    IGRAPH_ASSERT(igraph_vector_init(&capacity, 0) == IGRAPH_SUCCESS);

    igraph_integer_t source, target;

    // Do the fuzzing
    igraph_t g;
    if (igraph_read_graph_dimacs_flow(&g, ifile, &problem, &label, &source, &target, &capacity, IGRAPH_DIRECTED) == IGRAPH_SUCCESS) {
        // Clean up
        igraph_destroy(&g);
    }

    // no need to call igraph_destroy() if igraph_read_graph_dl() returns an
    // error code as we don't have a valid graph object in that case

    igraph_vector_destroy(&capacity);
    igraph_vector_int_destroy(&label);
    igraph_strvector_destroy(&problem);

    fclose(ifile);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;
}
