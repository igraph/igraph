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

#include "igraph_attributes.h"
#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_vector_ptr.h"
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

/**
 * Checks if a value read succeeded (failure being indicated by negative \p value) and classifies
 * the failure as end-of-file or read error.
 *
 * \param value Expected to be a value returned by igraph_i_read_graph_graphdb_getword()
 * \param instream The input stream used to read \p value from.
 */
static igraph_error_t handle_input_error(igraph_integer_t value, FILE *instream) {
    if (value < 0) {
        if (feof(instream)) {
            IGRAPH_ERROR("Unexpected end of file, truncated graphdb file.", IGRAPH_PARSEERROR);
        } else {
            IGRAPH_ERROR("Cannot read from file.", IGRAPH_EFILE);
        }
    }
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_read_graph_graphdb
 * \brief Read a graph in the unlabeled binary graph database format.
 *
 * This is a binary format for storing unlabeled graphs,
 * used in the ARG Graph Database
 * for isomorphism testing. For more information, see
 * https://mivia.unisa.it/datasets/graph-database/arg-database/
 * </para>
 *
 * <para>
 * This function reads the unlabeled version of the format.
 * See \ref igraph_read_graph_graphdb() for reading the labeled
 * version.</para>
 *
 * <para>
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

    /* Read vertex count */
    const igraph_integer_t vcount = igraph_i_read_graph_graphdb_getword(instream);
    IGRAPH_CHECK(handle_input_error(vcount, instream));

    /* Read edges */

    igraph_vector_int_t edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 100);
    igraph_vector_int_clear(&edges);

    for (igraph_integer_t i = 0; i < vcount; i++) {
        igraph_integer_t len = igraph_i_read_graph_graphdb_getword(instream);
        IGRAPH_CHECK(handle_input_error(len, instream));
        for (igraph_integer_t j = 0; j < len; j++) {
            igraph_integer_t to = igraph_i_read_graph_graphdb_getword(instream);
            IGRAPH_CHECK(handle_input_error(to, instream));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
            IGRAPH_ALLOW_INTERRUPTION();
        }
    }

    /* Check that the complete file has been read */
    if (fgetc(instream) != EOF) {
        IGRAPH_ERROR("Extra bytes at end of graphdb file.", IGRAPH_PARSEERROR);
    }

    /* Create graph */
    IGRAPH_CHECK(igraph_create(graph, &edges, vcount, directed));

    /* Clean up */
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_read_graph_graphdb_labeled
 * \brief Read a graph in the labeled binary graph database format.
 *
 * This is a binary format for storing vertex and edge labeled graphs,
 * used in the ARG Graph Database
 * for isomorphism testing. For more information, see
 * https://mivia.unisa.it/datasets/graph-database/arg-database/
 * </para>
 *
 * <para>
 * This function reads the labeled version of the format. Vertex and
 * edge labels will be stored in numeric vertex and edge attributes
 * named "label". To read the unlabeled version of the format, see
 * \ref igraph_read_graph_graphdb().
 * </para>
 *
 * <para>
 * From the graph database homepage:
 * </para>
 *
 * \blockquote <para>
 * The file is composed of 16 bit words, which are represented using
 * the so-called little-endian convention, i.e. the least significant
 * byte of the word is stored first.</para>
 *
 * <para>
 * The first word contains the number of nodes in the graph; this means
 * that this format can deal with graph of up to 65535 nodes (however,
 * actual graphs in the database are quite smaller, up to 100 nodes).
 * Then follow N words, that are the attribute of corresponding node.
 * Then, for each node i, the file contains the number of edges coming
 * out of the node itself (Ei) followed by 2*Ei words that code the edges.
 * In particular, for each edge, the first word represents the destination
 * node of the edge and the second word represents the edge attribute.
 * Node numeration is 0-based, so the first node of the graph has index 0.
 * </para> \endblockquote
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param instream The stream to read from. It should be opened
 *    in binary mode.
 * \param directed Logical scalar, whether to create a directed graph.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the
 * number of edges.
 */

igraph_error_t igraph_read_graph_graphdb_labeled(igraph_t *graph, FILE *instream,
                              igraph_bool_t directed) {

    /* Read vertex count */
    const igraph_integer_t vcount = igraph_i_read_graph_graphdb_getword(instream);
    IGRAPH_CHECK(handle_input_error(vcount, instream));

    /* The largest possible value of 'vcount' is 65535 so there is no
     * risk of an extreme size memory allocation due to a corrupted file. */
    igraph_vector_t vertex_labels;
    IGRAPH_VECTOR_INIT_FINALLY(&vertex_labels, vcount);

    /* Read vertex labels */
    for (igraph_integer_t i=0; i < vcount; i++) {
        igraph_integer_t label = igraph_i_read_graph_graphdb_getword(instream);
        IGRAPH_CHECK(handle_input_error(label, instream));
        VECTOR(vertex_labels)[i] = label;
    }

    /* Read edges and edge labels */

    igraph_vector_t edge_labels;
    IGRAPH_VECTOR_INIT_FINALLY(&edge_labels, 50);
    igraph_vector_clear(&edge_labels);

    igraph_vector_int_t edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 100);
    igraph_vector_int_clear(&edges);

    for (igraph_integer_t i=0; i < vcount; i++) {
        igraph_integer_t len = igraph_i_read_graph_graphdb_getword(instream);
        IGRAPH_CHECK(handle_input_error(len, instream));
        for (igraph_integer_t j = 0; j < len; j++) {
            igraph_integer_t to = igraph_i_read_graph_graphdb_getword(instream);
            IGRAPH_CHECK(handle_input_error(to, instream));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));

            igraph_integer_t label = igraph_i_read_graph_graphdb_getword(instream);
            IGRAPH_CHECK(handle_input_error(label, instream));
            IGRAPH_CHECK(igraph_vector_push_back(&edge_labels, label));

            IGRAPH_ALLOW_INTERRUPTION();
        }
    }

    /* Check that the complete file has been read */
    if (fgetc(instream) != EOF) {
        IGRAPH_ERROR("Extra bytes at end of graphdb file.", IGRAPH_PARSEERROR);
    }

    /* Set up attribute records */

    igraph_vector_ptr_t vattrs, eattrs;
    igraph_attribute_record_t vrec, erec;

    IGRAPH_VECTOR_PTR_INIT_FINALLY(&vattrs, 1);
    vrec.name = "label";
    vrec.type = IGRAPH_ATTRIBUTE_NUMERIC;
    vrec.value = &vertex_labels;
    VECTOR(vattrs)[0] = &vrec;

    IGRAPH_VECTOR_PTR_INIT_FINALLY(&eattrs, 1);
    erec.name = "label";
    erec.type = IGRAPH_ATTRIBUTE_NUMERIC;
    erec.value = &edge_labels;
    VECTOR(eattrs)[0] = &erec;

    /* Create graph */
    IGRAPH_CHECK(igraph_empty(graph, 0, directed));
    IGRAPH_FINALLY(igraph_destroy, graph);
    IGRAPH_CHECK(igraph_add_vertices(graph, vcount, &vattrs));
    IGRAPH_CHECK(igraph_add_edges(graph, &edges, &eattrs));

    /* Clean up */
    igraph_vector_ptr_destroy(&eattrs);
    igraph_vector_ptr_destroy(&vattrs);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&edge_labels);
    igraph_vector_destroy(&vertex_labels);
    IGRAPH_FINALLY_CLEAN(6); /* +1 for the graph */

    return IGRAPH_SUCCESS;
}
