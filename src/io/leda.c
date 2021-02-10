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
#include "igraph_interface.h"
#include "igraph_iterators.h"

#include "graph/attributes.h"

#include <string.h>

#define CHECK(cmd) do { ret=cmd; if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE); } while (0)

/**
 * \function igraph_write_graph_leda
 * \brief Write a graph in LEDA native graph format.
 *
 * This function writes a graph to an output stream in LEDA format.
 * See http://www.algorithmic-solutions.info/leda_guide/graphs/leda_native_graph_fileformat.html
 *
 * </para><para>
 * The support for the LEDA format is very basic at the moment; igraph
 * writes only the LEDA graph section which supports one selected vertex
 * and edge attribute and no layout information or visual attributes.
 *
 * \param graph The graph to write to the stream.
 * \param outstream The stream.
 * \param vertex_attr_name The name of the vertex attribute whose values
 *                         are to be stored in the output or \c NULL if no
 *                         vertex attribute has to be stored.
 * \param edge_attr_name   The name of the edge attribute whose values
 *                         are to be stored in the output or \c NULL if no
 *                         edge attribute has to be stored.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the number of vertices and edges in the
 * graph.
 *
 * \example examples/simple/igraph_write_graph_leda.c
 */

int igraph_write_graph_leda(const igraph_t *graph, FILE *outstream,
                            const char* vertex_attr_name,
                            const char* edge_attr_name) {
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_eit_t it;
    long int i = 0;
    int ret;
    igraph_attribute_type_t vertex_attr_type = IGRAPH_ATTRIBUTE_DEFAULT;
    igraph_attribute_type_t edge_attr_type = IGRAPH_ATTRIBUTE_DEFAULT;
    igraph_integer_t from, to, rev;

    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_FROM),
                                   &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);

    /* Check if we have the vertex attribute */
    if (vertex_attr_name &&
        !igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX, vertex_attr_name)) {
        vertex_attr_name = 0;
        IGRAPH_WARNING("specified vertex attribute does not exist");
    }
    if (vertex_attr_name) {
        IGRAPH_CHECK(igraph_i_attribute_gettype(graph, &vertex_attr_type,
                                                IGRAPH_ATTRIBUTE_VERTEX, vertex_attr_name));
        if (vertex_attr_type != IGRAPH_ATTRIBUTE_NUMERIC &&
            vertex_attr_type != IGRAPH_ATTRIBUTE_STRING) {
            vertex_attr_name = 0; vertex_attr_type = IGRAPH_ATTRIBUTE_DEFAULT;
            IGRAPH_WARNING("specified vertex attribute must be numeric or string");
        }
    }

    /* Check if we have the edge attribute */
    if (edge_attr_name &&
        !igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE, edge_attr_name)) {
        edge_attr_name = 0;
        IGRAPH_WARNING("specified edge attribute does not exist");
    }
    if (edge_attr_name) {
        IGRAPH_CHECK(igraph_i_attribute_gettype(graph, &edge_attr_type,
                                                IGRAPH_ATTRIBUTE_EDGE, edge_attr_name));
        if (edge_attr_type != IGRAPH_ATTRIBUTE_NUMERIC &&
            edge_attr_type != IGRAPH_ATTRIBUTE_STRING) {
            edge_attr_name = 0; edge_attr_type = IGRAPH_ATTRIBUTE_DEFAULT;
            IGRAPH_WARNING("specified edge attribute must be numeric or string");
        }
    }

    /* Start writing header */
    CHECK(fprintf(outstream, "LEDA.GRAPH\n"));

    switch (vertex_attr_type) {
    case IGRAPH_ATTRIBUTE_NUMERIC:
        CHECK(fprintf(outstream, "float\n"));
        break;
    case IGRAPH_ATTRIBUTE_STRING:
        CHECK(fprintf(outstream, "string\n"));
        break;
    default:
        CHECK(fprintf(outstream, "void\n"));
    }

    switch (edge_attr_type) {
    case IGRAPH_ATTRIBUTE_NUMERIC:
        CHECK(fprintf(outstream, "float\n"));
        break;
    case IGRAPH_ATTRIBUTE_STRING:
        CHECK(fprintf(outstream, "string\n"));
        break;
    default:
        CHECK(fprintf(outstream, "void\n"));
    }

    CHECK(fprintf(outstream, "%d\n", (igraph_is_directed(graph) ? -1 : -2)));

    /* Start writing vertices */
    CHECK(fprintf(outstream, "# Vertices\n"));
    CHECK(fprintf(outstream, "%ld\n", no_of_nodes));

    if (vertex_attr_type == IGRAPH_ATTRIBUTE_NUMERIC) {
        /* Vertices with numeric attributes */
        igraph_vector_t values;

        IGRAPH_VECTOR_INIT_FINALLY(&values, no_of_nodes);
        IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(
                         graph, vertex_attr_name, igraph_vss_all(), &values));

        for (i = 0; i < no_of_nodes; i++) {
            CHECK(fprintf(outstream, "|{"));
            CHECK(igraph_real_fprintf_precise(outstream, VECTOR(values)[i]));
            CHECK(fprintf(outstream, "}|\n"));
        }

        igraph_vector_destroy(&values);
        IGRAPH_FINALLY_CLEAN(1);
    } else if (vertex_attr_type == IGRAPH_ATTRIBUTE_STRING) {
        /* Vertices with string attributes */
        igraph_strvector_t values;

        IGRAPH_CHECK(igraph_strvector_init(&values, no_of_nodes));
        IGRAPH_FINALLY(igraph_strvector_destroy, &values);

        IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(
                         graph, vertex_attr_name, igraph_vss_all(), &values));

        for (i = 0; i < no_of_nodes; i++) {
            const char* str = STR(values, i);
            if (strchr(str, '\n') != 0) {
                IGRAPH_ERROR("edge attribute values cannot contain newline characters",
                             IGRAPH_EINVAL);
            }
            CHECK(fprintf(outstream, "|{%s}|\n", str));
        }

        igraph_strvector_destroy(&values);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        /* Vertices with no attributes */
        for (i = 0; i < no_of_nodes; i++) {
            CHECK(fprintf(outstream, "|{}|\n"));
        }
    }

    CHECK(fprintf(outstream, "# Edges\n"));
    CHECK(fprintf(outstream, "%ld\n", no_of_edges));

    if (edge_attr_type == IGRAPH_ATTRIBUTE_NUMERIC) {
        /* Edges with numeric attributes */
        igraph_vector_t values;
        IGRAPH_VECTOR_INIT_FINALLY(&values, no_of_nodes);
        IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(
                         graph, edge_attr_name, igraph_ess_all(IGRAPH_EDGEORDER_ID), &values));
        while (!IGRAPH_EIT_END(it)) {
            long int eid = IGRAPH_EIT_GET(it);
            igraph_edge(graph, (igraph_integer_t) eid, &from, &to);
            igraph_get_eid(graph, &rev, to, from, 1, 0);
            if (rev == IGRAPH_EIT_GET(it)) {
                rev = -1;
            }
            CHECK(fprintf(outstream, "%ld %ld %ld |{",
                          (long int) from + 1, (long int) to + 1,
                          (long int) rev + 1));
            CHECK(igraph_real_fprintf_precise(outstream, VECTOR(values)[eid]));
            CHECK(fprintf(outstream, "}|\n"));
            IGRAPH_EIT_NEXT(it);
        }
        igraph_vector_destroy(&values);
        IGRAPH_FINALLY_CLEAN(1);
    } else if (edge_attr_type == IGRAPH_ATTRIBUTE_STRING) {
        /* Edges with string attributes */
        igraph_strvector_t values;
        IGRAPH_CHECK(igraph_strvector_init(&values, no_of_nodes));
        IGRAPH_FINALLY(igraph_strvector_destroy, &values);
        IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(
                         graph, edge_attr_name, igraph_ess_all(IGRAPH_EDGEORDER_ID), &values));
        while (!IGRAPH_EIT_END(it)) {
            long int eid = IGRAPH_EIT_GET(it);
            const char* str = STR(values, eid);
            igraph_edge(graph, (igraph_integer_t) eid, &from, &to);
            igraph_get_eid(graph, &rev, to, from, 1, 0);
            if (rev == IGRAPH_EIT_GET(it)) {
                rev = -1;
            }
            if (strchr(str, '\n') != 0) {
                IGRAPH_ERROR("edge attribute values cannot contain newline characters",
                             IGRAPH_EINVAL);
            }
            CHECK(fprintf(outstream, "%ld %ld %ld |{%s}|\n",
                          (long int) from + 1, (long int) to + 1,
                          (long int) rev + 1, str));
            IGRAPH_EIT_NEXT(it);
        }
        igraph_strvector_destroy(&values);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        /* Edges with no attributes */
        while (!IGRAPH_EIT_END(it)) {
            igraph_edge(graph, IGRAPH_EIT_GET(it), &from, &to);
            igraph_get_eid(graph, &rev, to, from, 1, 0);
            if (rev == IGRAPH_EIT_GET(it)) {
                rev = -1;
            }
            CHECK(fprintf(outstream, "%ld %ld %ld |{}|\n",
                          (long int) from + 1, (long int) to + 1,
                          (long int) rev + 1));
            IGRAPH_EIT_NEXT(it);
        }
    }

    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

#undef CHECK
