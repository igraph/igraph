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

#include "core/trie.h"

static int igraph_i_read_graph_graphdb_getword(FILE *instream) {
    int b1, b2;
    unsigned char c1, c2;
    b1 = fgetc(instream);
    b2 = fgetc(instream);
    if (b1 != EOF) {
        c1 = (unsigned char) b1; c2 = (unsigned char) b2;
        return c1 | (c2 << 8);
    } else {
        return -1;
    }
}

/**
 * \function igraph_read_graph_graphdb
 * \brief Read a graph in the binary graph database format.
 *
 * This is a binary format, used in the graph database
 * for isomorphism testing. From the (now defunct) graph database
 * homepage:
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
 * Only unlabelled graphs are implemented.
 * \param graph Pointer to an uninitialized graph object.
 * \param instream The stream to read from.
 * \param directed Logical scalar, whether to create a directed graph.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the
 * number of edges.
 *
 * \example examples/simple/igraph_read_graph_graphdb.c
 */

int igraph_read_graph_graphdb(igraph_t *graph, FILE *instream,
                              igraph_bool_t directed) {

    igraph_vector_t edges;
    long int nodes;
    long int i, j;
    igraph_bool_t end = 0;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    nodes = igraph_i_read_graph_graphdb_getword(instream);
    if (nodes < 0) {
        IGRAPH_ERROR("Can't read from file", IGRAPH_EFILE);
    }
    for (i = 0; !end && i < nodes; i++) {
        long int len = igraph_i_read_graph_graphdb_getword(instream);
        if (len < 0) {
            end = 1;
            break;
        }
        for (j = 0; ! end && j < len; j++) {
            long int to = igraph_i_read_graph_graphdb_getword(instream);
            if (to < 0) {
                end = 1;
                break;
            }
            IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
            IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
        }
    }

    if (end) {
        IGRAPH_ERROR("Truncated graphdb file", IGRAPH_EFILE);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) nodes,
                               directed));
    igraph_vector_destroy(&edges);

    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
