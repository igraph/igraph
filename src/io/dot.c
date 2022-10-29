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

#include <ctype.h>
#include <string.h>

#define CHECK(cmd) do { int ret=cmd; if (ret<0) IGRAPH_ERROR("Writing DOT format failed.", IGRAPH_EFILE); } while (0)

static igraph_error_t igraph_i_dot_escape(const char *orig, char **result) {
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
                need_quote = 1;
                newlen++;
            }
        } else if (orig[i] == '_') {
            is_number = 0; newlen++;
        } else if (orig[i] == '\\' || orig[i] == '"' || orig[i] == '\n') {
            need_quote = 1; is_number = 0; newlen += 2; /* will be escaped */
        } else if (isalpha(orig[i])) {
            is_number = 0; newlen++;
        } else {
            is_number = 0; need_quote = 1; newlen++;
        }
    }
    if (is_number && len > 0 && orig[len - 1] == '.') {
        is_number = 0;
    }
    if (!is_number && isdigit(orig[0])) {
        need_quote = 1;
    }

    if (is_number || !need_quote) {
        *result = strdup(orig);
        if (!*result) {
            IGRAPH_ERROR("Writing DOT format failed.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
    } else {
        *result = IGRAPH_CALLOC(newlen + 3, char);
        if (!*result) {
            IGRAPH_ERROR("Writing DOT format failed.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        (*result)[0] = '"';
        (*result)[newlen + 1] = '"';
        (*result)[newlen + 2] = '\0';
        /* Escape quotes, backslashes and newlines.
         * Even though the format spec at https://graphviz.org/doc/info/lang.html
         * claims that only quotes need escaping, escaping backslashes appears to
         * be necessary as well for GraphViz to render labels correctly.
         * Tested with GraphViz 2.50. */
        for (i = 0, j = 1; i < len; i++) {
            if (orig[i] == '\n') {
                (*result)[j++] = '\\';
                (*result)[j++] = 'n';
                continue;
            }
            if (orig[i] == '\\' || orig[i] == '"') {
                (*result)[j++] = '\\';
            }
            (*result)[j++] = orig[i];
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_write_graph_dot
 * \brief Write the graph to a stream in DOT format.
 *
 * DOT is the format used by the widely known GraphViz software, see
 * http://www.graphviz.org for details. The grammar of the DOT format
 * can be found here: http://www.graphviz.org/doc/info/lang.html
 *
 * </para><para>This is only a preliminary implementation, no visualization
 * information is written.
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
    igraph_integer_t i, j;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    char edgeop[3];
    igraph_strvector_t gnames, vnames, enames;
    igraph_vector_int_t gtypes, vtypes, etypes;
    igraph_vector_t numv;
    igraph_strvector_t strv;
    igraph_vector_bool_t boolv;

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

    IGRAPH_VECTOR_INIT_FINALLY(&numv, 1);
    IGRAPH_STRVECTOR_INIT_FINALLY(&strv, 1);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&boolv, 1);

    CHECK(fprintf(outstream, "/* Created by igraph %s */\n",
                  IGRAPH_VERSION));

    if (igraph_is_directed(graph)) {
        CHECK(fprintf(outstream, "digraph {\n"));
        strcpy(edgeop, "->");
    } else {
        CHECK(fprintf(outstream, "graph {\n"));
        strcpy(edgeop, "--");
    }

    /* Write the graph attributes */
    if (igraph_vector_int_size(&gtypes) > 0) {
        CHECK(fprintf(outstream, "  graph [\n"));
        for (i = 0; i < igraph_vector_int_size(&gtypes); i++) {
            const char *name;
            char *newname;
            name = igraph_strvector_get(&gnames, i);
            IGRAPH_CHECK(igraph_i_dot_escape(name, &newname));
            IGRAPH_FINALLY(igraph_free, newname);
            if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
                IGRAPH_CHECK(igraph_i_attribute_get_numeric_graph_attr(graph, name, &numv));
                if (VECTOR(numv)[0] == (igraph_integer_t)VECTOR(numv)[0]) {
                    CHECK(fprintf(outstream, "    %s=%" IGRAPH_PRId "\n", newname, (igraph_integer_t)VECTOR(numv)[0]));
                } else {
                    CHECK(fprintf(outstream, "    %s=", newname));
                    CHECK(igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]));
                    CHECK(fputc('\n', outstream));
                }
            } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
                const char *s;
                char *news;
                IGRAPH_CHECK(igraph_i_attribute_get_string_graph_attr(graph, name, &strv));
                s = igraph_strvector_get(&strv, 0);
                IGRAPH_CHECK(igraph_i_dot_escape(s, &news));
                CHECK(fprintf(outstream, "    %s=%s\n", newname, news));
                IGRAPH_FREE(news);
            } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
                IGRAPH_CHECK(igraph_i_attribute_get_bool_graph_attr(graph, name, &boolv));
                CHECK(fprintf(outstream, "    %s=%d\n", newname, VECTOR(boolv)[0] ? 1 : 0));
                IGRAPH_WARNING("A boolean graph attribute was converted to numeric");
            } else {
                IGRAPH_WARNING("A non-numeric, non-string, non-boolean graph attribute ignored");
            }
            IGRAPH_FREE(newname);
            IGRAPH_FINALLY_CLEAN(1);
        }
        CHECK(fprintf(outstream, "  ];\n"));
    }

    /* Write the vertices */
    if (igraph_vector_int_size(&vtypes) > 0) {
        for (i = 0; i < no_of_nodes; i++) {
            CHECK(fprintf(outstream, "  %" IGRAPH_PRId " [\n", i));
            for (j = 0; j < igraph_vector_int_size(&vtypes); j++) {
                const char *name;
                char *newname;
                name = igraph_strvector_get(&vnames, j);
                IGRAPH_CHECK(igraph_i_dot_escape(name, &newname));
                IGRAPH_FINALLY(igraph_free, newname);
                if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
                    IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(graph, name, igraph_vss_1(i), &numv));
                    if (VECTOR(numv)[0] == floor(VECTOR(numv)[0])) {
                        CHECK(fprintf(outstream, "    %s=%g\n", newname, VECTOR(numv)[0]));
                    } else {
                        CHECK(fprintf(outstream, "    %s=", newname));
                        CHECK(igraph_real_fprintf_precise(outstream,
                                                          VECTOR(numv)[0]));
                        CHECK(fputc('\n', outstream));
                    }
                } else if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_STRING) {
                    const char *s;
                    char *news;
                    IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, name, igraph_vss_1(i), &strv));
                    s = igraph_strvector_get(&strv, 0);
                    IGRAPH_CHECK(igraph_i_dot_escape(s, &news));
                    CHECK(fprintf(outstream, "    %s=%s\n", newname, news));
                    IGRAPH_FREE(news);
                } else if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_BOOLEAN) {
                    IGRAPH_CHECK(igraph_i_attribute_get_bool_vertex_attr(graph, name, igraph_vss_1(i), &boolv));
                    CHECK(fprintf(outstream, "    %s=%d\n", newname, VECTOR(boolv)[0] ? 1 : 0));
                    IGRAPH_WARNING("A boolean vertex attribute was converted to numeric");
                } else {
                    IGRAPH_WARNING("A non-numeric, non-string, non-boolean vertex attribute was ignored");
                }
                IGRAPH_FREE(newname);
                IGRAPH_FINALLY_CLEAN(1);
            }
            CHECK(fprintf(outstream, "  ];\n"));
        }
    } else {
        for (i = 0; i < no_of_nodes; i++) {
            CHECK(fprintf(outstream, "  %" IGRAPH_PRId ";\n", i));
        }
    }
    CHECK(fprintf(outstream, "\n"));

    /* Write the edges */
    if (igraph_vector_int_size(&etypes) > 0) {
        for (i = 0; i < no_of_edges; i++) {
            igraph_integer_t from = IGRAPH_FROM(graph, i);
            igraph_integer_t to = IGRAPH_TO(graph, i);
            CHECK(fprintf(outstream, "  %" IGRAPH_PRId " %s %" IGRAPH_PRId " [\n", from, edgeop, to));
            for (j = 0; j < igraph_vector_int_size(&etypes); j++) {
                const char *name;
                char *newname;
                name = igraph_strvector_get(&enames, j);
                IGRAPH_CHECK(igraph_i_dot_escape(name, &newname));
                IGRAPH_FINALLY(igraph_free, newname);
                if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
                    IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph,
                                 name, igraph_ess_1(i), &numv));
                    if (VECTOR(numv)[0] == floor(VECTOR(numv)[0])) {
                        CHECK(fprintf(outstream, "    %s=%g\n", newname, VECTOR(numv)[0]));
                    } else {
                        CHECK(fprintf(outstream, "    %s=", newname));
                        CHECK(igraph_real_fprintf_precise(outstream, VECTOR(numv)[0]));
                        CHECK(fputc('\n', outstream));
                    }
                } else if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_STRING) {
                    const char *s;
                    char *news;
                    IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph,
                                 name, igraph_ess_1(i), &strv));
                    s = igraph_strvector_get(&strv, 0);
                    IGRAPH_CHECK(igraph_i_dot_escape(s, &news));
                    CHECK(fprintf(outstream, "    %s=%s\n", newname, news));
                    IGRAPH_FREE(news);
                } else if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_BOOLEAN) {
                    IGRAPH_CHECK(igraph_i_attribute_get_bool_edge_attr(graph,
                                 name, igraph_ess_1(i), &boolv));
                    CHECK(fprintf(outstream, "    %s=%d\n", newname, VECTOR(boolv)[0] ? 1 : 0));
                    IGRAPH_WARNING("A boolean edge attribute was converted to numeric");
                } else {
                    IGRAPH_WARNING("A non-numeric, non-string graph attribute ignored");
                }
                IGRAPH_FREE(newname);
                IGRAPH_FINALLY_CLEAN(1);
            }
            CHECK(fprintf(outstream, "  ];\n"));
        }
    } else {
        for (i = 0; i < no_of_edges; i++) {
            igraph_integer_t from = IGRAPH_FROM(graph, i);
            igraph_integer_t to = IGRAPH_TO(graph, i);
            CHECK(fprintf(outstream, "  %" IGRAPH_PRId " %s %" IGRAPH_PRId ";\n", from, edgeop, to));
        }
    }
    CHECK(fprintf(outstream, "}\n"));

    igraph_vector_bool_destroy(&boolv);
    igraph_strvector_destroy(&strv);
    igraph_vector_destroy(&numv);
    igraph_vector_int_destroy(&etypes);
    igraph_vector_int_destroy(&vtypes);
    igraph_vector_int_destroy(&gtypes);
    igraph_strvector_destroy(&enames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&gnames);
    IGRAPH_FINALLY_CLEAN(9);

    return IGRAPH_SUCCESS;
}

#undef CHECK
