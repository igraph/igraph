/*
   IGraph library.
   Copyright (C) 2024  The igraph development team

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
#include <cstdio>

// This fuzzer checks two things:
//  - Test writers for formats that igraph can only write but not read
//  - Test that files that igraph writes can be read back

#define CHECK_ERR(expr) \
    do { \
        igraph_error_handler_t *handler = igraph_set_error_handler(igraph_error_handler_abort); \
        expr; \
        igraph_set_error_handler(handler); \
    } while (0)

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

    // Do the fuzzing
    igraph_t g;
    if (igraph_read_graph_gml(&g, ifile) == IGRAPH_SUCCESS) {
        FILE *file;
        igraph_t g2;

        file = tmpfile();
        IGRAPH_ASSERT(file != NULL);
        // Do not check error because:
        // "Vertex attribute values cannot contain newline characters."
        igraph_write_graph_leda(&g, file, "label", "weight");
        fclose(file);

        file = tmpfile();
        IGRAPH_ASSERT(file != NULL);
        CHECK_ERR(igraph_write_graph_dot(&g, file));
        fclose(file);

        file = tmpfile();
        IGRAPH_ASSERT(file != NULL);
        CHECK_ERR(igraph_write_graph_gml(&g, file, IGRAPH_WRITE_GML_DEFAULT_SW, NULL, "no one"));
        rewind(file);
        CHECK_ERR(igraph_read_graph_gml(&g2, file));
        igraph_destroy(&g2);
        fclose(file);

        // Reading Pajek files back is disabled because of
        // https://github.com/igraph/igraph/issues/2560
        file = tmpfile();
        IGRAPH_ASSERT(file != NULL);
        CHECK_ERR(igraph_write_graph_pajek(&g, file));
        /*
        rewind(file);
        CHECK_ERR(igraph_read_graph_pajek(&g2, file));
        igraph_destroy(&g2);
        */
        fclose(file);

        file = tmpfile();
        IGRAPH_ASSERT(file != NULL);
        // Do not check error when writing GraphML because string attributes
        // may contain forbidden control characters.
        if (igraph_write_graph_graphml(&g, file, false) == IGRAPH_SUCCESS) {
            rewind(file);
            // Do not check error when reading because strings may not be
            // in a valid encoding, which confuses libxml2.
            if (igraph_read_graph_graphml(&g2, file, 0) == IGRAPH_SUCCESS) {
                igraph_destroy(&g2);
            }
        }
        fclose(file);

        file = tmpfile();
        IGRAPH_ASSERT(file != NULL);
        // Do not check error when writing LGL because not all possible
        // vertex names are supported by this format.
        if (igraph_write_graph_lgl(&g, file, "label", "weight", false) == IGRAPH_SUCCESS) {
            rewind(file);
            CHECK_ERR(igraph_read_graph_lgl(&g2, file, true, IGRAPH_ADD_WEIGHTS_IF_PRESENT, true));
            igraph_destroy(&g2);
        }
        fclose(file);

        file = tmpfile();
        IGRAPH_ASSERT(file != NULL);
        // Do not check error when writing NCOL because not all possible
        // vertex names are supported by this format.
        if (igraph_write_graph_ncol(&g, file, "label", "weight") == IGRAPH_SUCCESS) {
            rewind(file);
            CHECK_ERR(igraph_read_graph_ncol(&g2, file, NULL, true, IGRAPH_ADD_WEIGHTS_IF_PRESENT, true));
            igraph_destroy(&g2);
        }
        fclose(file);

        // Clean up
        igraph_destroy(&g);
    }

    // no need to call igraph_destroy() if igraph_read_graph_gml() returns an
    // error code as we don't have a valid graph object in that case

    fclose(ifile);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;
}
