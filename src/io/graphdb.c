/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2020  The igraph development team

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

#include "igraph_foreign.h"

#include "igraph_constructors.h"
#include "core/interruption.h"

/* Read a little-endian encoded 16-bit unsigned word.
 * Returns negative value on failure. */
static int igraph_i_read_graph_graphdb_getword(FILE *instream) {
    int b1, b2;
    unsigned char c1, c2;
    b1 = fgetc(instream);
    b2 = fgetc(instream);
    if (b1 != EOF && b2 != EOF) {
        c1 = (unsigned char) b1; c2 = (unsigned char) b2;
        return c1 | (c2 << 8);
    } else {
        return -1;
    }
}

/* Determine whether the read failed due to an input error or end-of-file condition.
 * Must only be called after a read failure, always returns a non-success error code. */
static igraph_error_t handle_input_error(FILE *instream) {
    if (feof(instream)) {
        IGRAPH_ERROR("Unexpected end of file, truncated graphdb file.", IGRAPH_PARSEERROR);
    } else {
        IGRAPH_ERROR("Cannot read from file.", IGRAPH_EFILE);
    }
}

/**
 * \function igraph_read_graph_graphdb
 * \brief Read a graph in the binary graph database format.
 *
 * This is a binary format, used in the ARG Graph Database
 * for isomorphism testing. For more information, see
 * https://mivia.unisa.it/datasets/graph-database/arg-database/
 *
 * </para><para>
 * From the graph database homepage:
 * </para>
 *
 * \blockquote <para>
 * The graphs are stored in a compact binary format, one graph per
 * file. The file is composed of 16 bit words, which are represented
 * using the so-called little-endian convention, i.e. the least
 * significant byte of the word is stored first.</para>
 *
 * <para>
 * Then, for each node, the file contains the list of edges coming
 * out of the node itself. The list is represented by a word encoding
 * its length, followed by a word for each edge, representing the
 * destination node of the edge. Node numeration is 0-based, so the
 * first node of the graph has index 0.</para> \endblockquote
 *
 * <para>
 * As of igraph 0.10, only unlabelled graphs are implemented.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param instream The stream to read from. It should be opened
 *    in binary mode.
 * \param directed Logical scalar, whether to create a directed graph.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the
 * number of edges.
 *
 * \example examples/simple/igraph_read_graph_graphdb.c
 */

igraph_error_t igraph_read_graph_graphdb(igraph_t *graph, FILE *instream,
                              igraph_bool_t directed) {

    const igraph_integer_t nodes = igraph_i_read_graph_graphdb_getword(instream);
    if (nodes < 0) {
        IGRAPH_CHECK(handle_input_error(instream));
    }

    igraph_vector_int_t edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 100);
    igraph_vector_int_clear(&edges);

    for (igraph_integer_t i = 0; i < nodes; i++) {
        igraph_integer_t len = igraph_i_read_graph_graphdb_getword(instream);
        if (len < 0) {
            IGRAPH_CHECK(handle_input_error(instream));
        }
        for (igraph_integer_t j = 0; j < len; j++) {
            igraph_integer_t to = igraph_i_read_graph_graphdb_getword(instream);
            if (to < 0) {
                IGRAPH_CHECK(handle_input_error(instream));
            }
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
            IGRAPH_ALLOW_INTERRUPTION();
        }
    }
    if (fgetc(instream) != EOF) {
        IGRAPH_ERROR("Extra bytes at end of graphdb file.", IGRAPH_PARSEERROR);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
