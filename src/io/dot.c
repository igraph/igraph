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
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_version.h"

#include "graph/attributes.h"
#include "internal/hacks.h" /* strcasecmp & strdup */
#include "math/safe_intop.h" /* IGRAPH_MAX_EXACT_REAL */

#include <ctype.h>
#include <string.h>

#define CHECK(cmd) do { int ret=cmd; if (ret<0) IGRAPH_ERROR("Writing DOT format failed.", IGRAPH_EFILE); } while (0)

static igraph_error_t dot_escape(const char *orig, igraph_vector_char_t* result) {
    /* do we have to escape the string at all? */
    igraph_integer_t i, j, len = strlen(orig), newlen = 0;
    igraph_bool_t need_quote = false, is_number = true;

    /* first, check whether the string is equal to some reserved word, or empty */
    if (!strcasecmp(orig, "graph") || !strcasecmp(orig, "digraph") ||
        !strcasecmp(orig, "node") || !strcasecmp(orig, "edge") ||
        !strcasecmp(orig, "strict") || !strcasecmp(orig, "subgraph") || len == 0) {
        need_quote = true;
        is_number = false;
    }

    /* next, check whether we need to escape the string for any other reason.
     * Also update is_number and newlen */
    for (i = 0; i < len; i++) {
        if (isdigit(orig[i])) {
            newlen++;
        } else if (orig[i] == '-' && i == 0) {
            newlen++;
        } else if (orig[i] == '.') {
            if (is_number) {
                newlen++;
            } else {
                need_quote = true;
                newlen++;
            }
        } else if (orig[i] == '_') {
            is_number = false; newlen++;
        } else if (orig[i] == '\\' || orig[i] == '"' || orig[i] == '\n') {
            need_quote = true; is_number = false; newlen += 2; /* will be escaped */
        } else if (isalpha(orig[i])) {
            is_number = false; newlen++;
        } else {
            is_number = false; need_quote = true; newlen++;
        }
    }
    if (is_number && len > 0 && orig[len - 1] == '.') {
        is_number = false;
    }
    if (!is_number && isdigit(orig[0])) {
        need_quote = true;
    }

    if (is_number || !need_quote) {
        IGRAPH_CHECK(igraph_vector_char_resize(result, newlen + 1));
        memcpy(VECTOR(*result), orig, newlen);
    } else {
        newlen += 2;
        IGRAPH_CHECK(igraph_vector_char_resize(result, newlen + 1));
        VECTOR(*result)[0] = '"';
        VECTOR(*result)[newlen - 1] = '"';

        /* Escape quotes, backslashes and newlines.
         * Even though the format spec at https://graphviz.org/doc/info/lang.html
         * claims that only quotes need escaping, escaping backslashes appears to
         * be necessary as well for GraphViz to render labels correctly.
         * Tested with GraphViz 2.50. */
        for (i = 0, j = 1; i < len; i++) {
            if (orig[i] == '\n') {
                VECTOR(*result)[j++] = '\\';
                VECTOR(*result)[j++] = 'n';
                continue;
            }
            if (orig[i] == '\\' || orig[i] == '"') {
                VECTOR(*result)[j++] = '\\';
            }
            VECTOR(*result)[j++] = orig[i];
        }
    }
    VECTOR(*result)[newlen] = 0;

    return IGRAPH_SUCCESS;
}

/* Writes exactly representable integral values in standard integer notation, without decimal points or e-notation.
 * Floating point values that are written with e-notation are quoted, otherwise the Graphviz parser cannot handle them.
 */
static igraph_error_t fprint_integral_or_precise(FILE *file, igraph_real_t x, igraph_vector_char_t *buf) {
    if (fabs(x) <= IGRAPH_MAX_EXACT_REAL && floor(x) == x) {
        /* write exactly representable integral values in standard integer notation;
         * the above conditional skips +-Inf and NaN */
        CHECK(fprintf(file, "%.f", x));
    } else {
        /* write as precise float and quote if necessary */
        char str[50]; /* large enough to hold any precisely printed real */
        CHECK(igraph_real_snprintf_precise(str, sizeof(str) / sizeof(str[0]), x));
        IGRAPH_CHECK(dot_escape(str, buf));
        CHECK(fputs(VECTOR(*buf), file));
    }
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_write_graph_dot
 * \brief Write the graph to a stream in DOT format.
 *
 * </para><para>
 * DOT is the format used by the widely known GraphViz software, see
 * http://www.graphviz.org for details. The grammar of the DOT format
 * can be found here: http://www.graphviz.org/doc/info/lang.html
 *
 * </para><para>
 * This is only a preliminary implementation, no visualization
 * information is written.
 *
 * </para><para>
 * This format is meant solely for interoperability with Graphviz.
 * It is not recommended for data exchange or archival.
 *
 * \param graph The graph to write to the stream.
 * \param outstream The stream to write the file to.
 *
 * Time complexity: should be proportional to the number of characters written
 * to the file.
 *
 * \sa \ref igraph_write_graph_graphml() for a more modern format.
 *
 * \example examples/simple/dot.c
 */
igraph_error_t igraph_write_graph_dot(const igraph_t *graph, FILE* outstream) {
    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    const igraph_integer_t no_of_edges = igraph_ecount(graph);
    const igraph_bool_t directed = igraph_is_directed(graph);
    const char *edgeop = directed ? "->" : "--";
    igraph_strvector_t gnames, vnames, enames;
    igraph_vector_int_t gtypes, vtypes, etypes;
    igraph_vector_t numv;
    igraph_strvector_t strv;
    igraph_vector_bool_t boolv;
    igraph_vector_char_t buf, buf2;

    IGRAPH_STRVECTOR_INIT_FINALLY(&gnames, 0);
    IGRAPH_STRVECTOR_INIT_FINALLY(&vnames, 0);
    IGRAPH_STRVECTOR_INIT_FINALLY(&enames, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&gtypes, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vtypes, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&etypes, 0);
    IGRAPH_CHECK(igraph_i_attribute_get_info(graph,
                 &gnames, &gtypes,
                 &vnames, &vtypes,
                 &enames, &etypes));

    IGRAPH_VECTOR_INIT_FINALLY(&numv, 0);
    IGRAPH_STRVECTOR_INIT_FINALLY(&strv, 0);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&boolv, 0);

    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&buf, 0);
    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&buf2, 0);

    CHECK(fprintf(outstream, "/* Created by igraph %s */\n", IGRAPH_VERSION));

    if (directed) {
        CHECK(fprintf(outstream, "digraph {\n"));
    } else {
        CHECK(fprintf(outstream, "graph {\n"));
    }

    /* Write the graph attributes */
    if (igraph_vector_int_size(&gtypes) > 0) {
        CHECK(fprintf(outstream, "  graph [\n"));
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(&gtypes); i++) {
            const char *name;
            name = igraph_strvector_get(&gnames, i);
            IGRAPH_CHECK(dot_escape(name, &buf));
            if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_graph_attr(graph, name, &numv));
                CHECK(fprintf(outstream, "    %s=", VECTOR(buf)));
                IGRAPH_CHECK(fprint_integral_or_precise(outstream, VECTOR(numv)[0], &buf));
                CHECK(fputc('\n', outstream));
            } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
                const char *s;
                IGRAPH_CHECK(igraph_i_attribute_get_string_graph_attr(graph, name, &strv));
                s = igraph_strvector_get(&strv, 0);
                IGRAPH_CHECK(dot_escape(s, &buf2));
                CHECK(fprintf(outstream, "    %s=%s\n", VECTOR(buf), VECTOR(buf2)));
            } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
                IGRAPH_CHECK(igraph_i_attribute_get_bool_graph_attr(graph, name, &boolv));
                CHECK(fprintf(outstream, "    %s=%d\n", VECTOR(buf), VECTOR(boolv)[0] ? 1 : 0));
                IGRAPH_WARNING("Boolean graph attribute was converted to numeric.");
            } else {
                IGRAPH_WARNING("A non-numeric, non-string, non-boolean graph attribute was ignored.");
            }
        }
        CHECK(fprintf(outstream, "  ];\n"));
    }

    /* Write the vertices */
    if (igraph_vector_int_size(&vtypes) > 0) {
        for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
            CHECK(fprintf(outstream, "  %" IGRAPH_PRId " [\n", i));
            for (igraph_integer_t j = 0; j < igraph_vector_int_size(&vtypes); j++) {
                const char *name;
                name = igraph_strvector_get(&vnames, j);
                IGRAPH_CHECK(dot_escape(name, &buf));
                if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
                    IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(graph, name, igraph_vss_1(i), &numv));
                    CHECK(fprintf(outstream, "    %s=", VECTOR(buf)));
                    IGRAPH_CHECK(fprint_integral_or_precise(outstream, VECTOR(numv)[0], &buf));
                    CHECK(fputc('\n', outstream));
                } else if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_STRING) {
                    const char *s;
                    IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, name, igraph_vss_1(i), &strv));
                    s = igraph_strvector_get(&strv, 0);
                    IGRAPH_CHECK(dot_escape(s, &buf2));
                    CHECK(fprintf(outstream, "    %s=%s\n", VECTOR(buf), VECTOR(buf2)));
                } else if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_BOOLEAN) {
                    IGRAPH_CHECK(igraph_i_attribute_get_bool_vertex_attr(graph, name, igraph_vss_1(i), &boolv));
                    CHECK(fprintf(outstream, "    %s=%d\n", VECTOR(buf), VECTOR(boolv)[0] ? 1 : 0));
                    IGRAPH_WARNING("A boolean vertex attribute was converted to numeric.");
                } else {
                    IGRAPH_WARNING("A non-numeric, non-string, non-boolean vertex attribute was ignored.");
                }
            }
            CHECK(fprintf(outstream, "  ];\n"));
        }
    } else {
        for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
            CHECK(fprintf(outstream, "  %" IGRAPH_PRId ";\n", i));
        }
    }
    CHECK(fprintf(outstream, "\n"));

    /* Write the edges */
    if (igraph_vector_int_size(&etypes) > 0) {
        for (igraph_integer_t i = 0; i < no_of_edges; i++) {
            igraph_integer_t from = IGRAPH_FROM(graph, i);
            igraph_integer_t to = IGRAPH_TO(graph, i);
            CHECK(fprintf(outstream, "  %" IGRAPH_PRId " %s %" IGRAPH_PRId " [\n", from, edgeop, to));
            for (igraph_integer_t j = 0; j < igraph_vector_int_size(&etypes); j++) {
                const char *name;
                name = igraph_strvector_get(&enames, j);
                IGRAPH_CHECK(dot_escape(name, &buf));
                if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
                    IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph,
                                 name, igraph_ess_1(i), &numv));
                    CHECK(fprintf(outstream, "    %s=", VECTOR(buf)));
                    IGRAPH_CHECK(fprint_integral_or_precise(outstream, VECTOR(numv)[0], &buf));
                    CHECK(fputc('\n', outstream));
                } else if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_STRING) {
                    const char *s;
                    IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph,
                                 name, igraph_ess_1(i), &strv));
                    s = igraph_strvector_get(&strv, 0);
                    IGRAPH_CHECK(dot_escape(s, &buf2));
                    CHECK(fprintf(outstream, "    %s=%s\n", VECTOR(buf), VECTOR(buf2)));
                } else if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_BOOLEAN) {
                    IGRAPH_CHECK(igraph_i_attribute_get_bool_edge_attr(graph,
                                 name, igraph_ess_1(i), &boolv));
                    CHECK(fprintf(outstream, "    %s=%d\n", VECTOR(buf), VECTOR(boolv)[0] ? 1 : 0));
                    IGRAPH_WARNING("A boolean edge attribute was converted to numeric.");
                } else {
                    IGRAPH_WARNING("A non-numeric, non-string graph attribute ignored.");
                }
            }
            CHECK(fprintf(outstream, "  ];\n"));
        }
    } else {
        for (igraph_integer_t i = 0; i < no_of_edges; i++) {
            igraph_integer_t from = IGRAPH_FROM(graph, i);
            igraph_integer_t to = IGRAPH_TO(graph, i);
            CHECK(fprintf(outstream, "  %" IGRAPH_PRId " %s %" IGRAPH_PRId ";\n", from, edgeop, to));
        }
    }
    CHECK(fprintf(outstream, "}\n"));

    igraph_vector_char_destroy(&buf2);
    igraph_vector_char_destroy(&buf);
    igraph_vector_bool_destroy(&boolv);
    igraph_strvector_destroy(&strv);
    igraph_vector_destroy(&numv);
    igraph_vector_int_destroy(&etypes);
    igraph_vector_int_destroy(&vtypes);
    igraph_vector_int_destroy(&gtypes);
    igraph_strvector_destroy(&enames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&gnames);
    IGRAPH_FINALLY_CLEAN(11);

    return IGRAPH_SUCCESS;
}

#undef CHECK
