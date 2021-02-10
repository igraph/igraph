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

#include "graph/attributes.h"

#include "lgl-header.h"

int igraph_lgl_yylex_init_extra (igraph_i_lgl_parsedata_t* user_defined,
                                 void* scanner);
void igraph_lgl_yylex_destroy (void *scanner );
int igraph_lgl_yyparse (igraph_i_lgl_parsedata_t* context);
void igraph_lgl_yyset_in  (FILE * in_str, void* yyscanner );

/**
 * \ingroup loadsave
 * \function igraph_read_graph_lgl
 * \brief Reads a graph from an <code>.lgl</code> file
 *
 * </para><para>
 * The <code>.lgl</code> format is used by the Large Graph
 * Layout visualization software
 * (http://lgl.sourceforge.net), it can
 * describe undirected optionally weighted graphs. From the LGL
 * manual:
 *
 * \blockquote <para>The second format is the LGL file format
 * (<code>.lgl</code> file
 * suffix). This is yet another graph file format that tries to be as
 * stingy as possible with space, yet keeping the edge file in a human
 * readable (not binary) format. The format itself is like the
 * following:
 * \verbatim # vertex1name
vertex2name [optionalWeight]
vertex3name [optionalWeight] \endverbatim
 * Here, the first vertex of an edge is preceded with a pound sign
 * '#'.  Then each vertex that shares an edge with that vertex is
 * listed one per line on subsequent lines.</para> \endblockquote
 *
 * </para><para>
 * LGL cannot handle loop and multiple edges or directed graphs, but
 * in \a igraph it is not an error to have multiple and loop edges.
 * \param graph Pointer to an uninitialized graph object.
 * \param instream A stream, it should be readable.
 * \param names Logical value, if TRUE the symbolic names of the
 *        vertices will be added to the graph as a vertex attribute
 *        called \quote name\endquote.
 * \param weights Whether to add the weights of the edges to the
 *        graph as an edge attribute called \quote weight\endquote.
 *        \c IGRAPH_ADD_WEIGHTS_YES adds the weights (even if they
 *        are not present in the file, in this case they are assumed
 *        to be zero). \c IGRAPH_ADD_WEIGHTS_NO does not add any
 *        edge attribute. \c IGRAPH_ADD_WEIGHTS_IF_PRESENT adds the
 *        attribute if and only if there is at least one explicit
 *        edge weight in the input file.
 * \param directed Whether to create a directed graph. As this format
 *        was originally used only for undirected graphs there is no
 *        information in the file about the directedness of the graph.
 *        Set this parameter to \c IGRAPH_DIRECTED or \c
 *        IGRAPH_UNDIRECTED to create a directed or undirected graph.
 * \return Error code:
 *         \c IGRAPH_PARSEERROR: if there is a
 *         problem reading the file, or the file is syntactically
 *         incorrect.
 *
 * Time complexity:
 * O(|V|+|E|log(|V|)) if we neglect
 * the time required by the parsing. As usual
 * |V| is the number of vertices,
 * while |E| is the number of edges.
 *
 * \sa \ref igraph_read_graph_ncol(), \ref igraph_write_graph_lgl()
 *
 * \example examples/simple/igraph_read_graph_lgl.c
 */
int igraph_read_graph_lgl(igraph_t *graph, FILE *instream,
                          igraph_bool_t names,
                          igraph_add_weights_t weights,
                          igraph_bool_t directed) {

    igraph_vector_t edges = IGRAPH_VECTOR_NULL, ws = IGRAPH_VECTOR_NULL;
    igraph_trie_t trie = IGRAPH_TRIE_NULL;
    igraph_vector_ptr_t name, weight;
    igraph_vector_ptr_t *pname = 0, *pweight = 0;
    igraph_attribute_record_t namerec, weightrec;
    const char *namestr = "name", *weightstr = "weight";
    igraph_i_lgl_parsedata_t context;

    IGRAPH_VECTOR_INIT_FINALLY(&ws, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_TRIE_INIT_FINALLY(&trie, names);

    context.has_weights = 0;
    context.vector = &edges;
    context.weights = &ws;
    context.trie = &trie;
    context.eof = 0;

    igraph_lgl_yylex_init_extra(&context, &context.scanner);
    IGRAPH_FINALLY(igraph_lgl_yylex_destroy, context.scanner);

    igraph_lgl_yyset_in(instream, context.scanner);

    if (igraph_lgl_yyparse(&context)) {
        if (context.errmsg[0] != 0) {
            IGRAPH_ERROR(context.errmsg, IGRAPH_PARSEERROR);
        } else {
            IGRAPH_ERROR("Cannot read LGL file", IGRAPH_PARSEERROR);
        }
    }

    IGRAPH_CHECK(igraph_empty(graph, 0, directed));
    IGRAPH_FINALLY(igraph_destroy, graph);

    if (names) {
        const igraph_strvector_t *namevec;
        IGRAPH_CHECK(igraph_vector_ptr_init(&name, 1));
        IGRAPH_FINALLY(igraph_vector_ptr_destroy, &name);
        pname = &name;
        igraph_trie_getkeys(&trie, &namevec); /* dirty */
        namerec.name = namestr;
        namerec.type = IGRAPH_ATTRIBUTE_STRING;
        namerec.value = namevec;
        VECTOR(name)[0] = &namerec;
    }

    if (weights == IGRAPH_ADD_WEIGHTS_YES ||
        (weights == IGRAPH_ADD_WEIGHTS_IF_PRESENT && context.has_weights)) {
        IGRAPH_CHECK(igraph_vector_ptr_init(&weight, 1));
        IGRAPH_FINALLY(igraph_vector_ptr_destroy, &weight);
        pweight = &weight;
        weightrec.name = weightstr;
        weightrec.type = IGRAPH_ATTRIBUTE_NUMERIC;
        weightrec.value = &ws;
        VECTOR(weight)[0] = &weightrec;
    }

    IGRAPH_CHECK(igraph_add_vertices(graph, (igraph_integer_t)
                                     igraph_trie_size(&trie), pname));
    IGRAPH_CHECK(igraph_add_edges(graph, &edges, pweight));

    if (pweight) {
        igraph_vector_ptr_destroy(pweight);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (pname) {
        igraph_vector_ptr_destroy(pname);
        IGRAPH_FINALLY_CLEAN(1);
    }
    igraph_trie_destroy(&trie);
    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&ws);
    igraph_lgl_yylex_destroy(context.scanner);
    IGRAPH_FINALLY_CLEAN(5);

    return 0;
}

/**
 * \ingroup loadsave
 * \function igraph_write_graph_lgl
 * \brief Writes the graph to a file in <code>.lgl</code> format
 *
 * </para><para>
 * <code>.lgl</code> is a format used by LGL, see \ref
 * igraph_read_graph_lgl() for details.
 *
 * </para><para>
 * Note that having multiple or loop edges in an
 * <code>.lgl</code> file breaks the  LGL software but \a igraph
 * does not check for this condition.
 * \param graph The graph to write.
 * \param outstream The stream object to write to, it should be
 *        writable.
 * \param names The name of the vertex attribute, if symbolic names
 *        are written to the file. If not supply 0 here.
 * \param weights The name of the edge attribute, if they are also
 *        written to the file. If you don't want weights supply 0
 *        here.
 * \param isolates Logical, if TRUE isolated vertices are also written
 *        to the file. If FALSE they will be omitted.
 * \return Error code:
 *         \c IGRAPH_EFILE if there is an error
 *         writing the file.
 *
 * Time complexity: O(|E|), the
 * number of edges if \p isolates is
 * FALSE, O(|V|+|E|) otherwise. All
 * file operations are expected to have time complexity
 * O(1).
 *
 * \sa \ref igraph_read_graph_lgl(), \ref igraph_write_graph_ncol()
 *
 * \example examples/simple/igraph_write_graph_lgl.c
 */
int igraph_write_graph_lgl(const igraph_t *graph, FILE *outstream,
                           const char *names, const char *weights,
                           igraph_bool_t isolates) {
    igraph_eit_t it;
    long int actvertex = -1;
    igraph_attribute_type_t nametype, weighttype;

    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_FROM),
                                   &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);

    /* Check if we have the names attribute */
    if (names && !igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX,
            names)) {
        names = 0;
        IGRAPH_WARNING("names attribute does not exists");
    }
    if (names) {
        IGRAPH_CHECK(igraph_i_attribute_gettype(graph, &nametype,
                                                IGRAPH_ATTRIBUTE_VERTEX, names));
        if (nametype != IGRAPH_ATTRIBUTE_NUMERIC && nametype != IGRAPH_ATTRIBUTE_STRING) {
            IGRAPH_WARNING("ignoring names attribute, unknown attribute type");
            names = 0;
        }
    }

    /* Check the weights as well */
    if (weights && !igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE,
            weights)) {
        weights = 0;
        IGRAPH_WARNING("weights attribute does not exists");
    }
    if (weights) {
        IGRAPH_CHECK(igraph_i_attribute_gettype(graph, &weighttype,
                                                IGRAPH_ATTRIBUTE_EDGE, weights));
        if (weighttype != IGRAPH_ATTRIBUTE_NUMERIC && weighttype != IGRAPH_ATTRIBUTE_STRING) {
            IGRAPH_WARNING("ignoring weights attribute, unknown attribute type");
            weights = 0;
        }
    }

    if (names == 0 && weights == 0) {
        /* No names, no weights */
        while (!IGRAPH_EIT_END(it)) {
            igraph_integer_t from, to;
            int ret;
            igraph_edge(graph, IGRAPH_EIT_GET(it), &from, &to);
            if (from == actvertex) {
                ret = fprintf(outstream, "%li\n", (long int)to);
            } else {
                actvertex = from;
                ret = fprintf(outstream, "# %li\n%li\n", (long int)from, (long int)to);
            }
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
            IGRAPH_EIT_NEXT(it);
        }
    } else if (weights == 0) {
        /* No weights but use names */
        igraph_strvector_t nvec;
        IGRAPH_CHECK(igraph_strvector_init(&nvec, igraph_vcount(graph)));
        IGRAPH_FINALLY(igraph_strvector_destroy, &nvec);
        IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, names,
                     igraph_vss_all(),
                     &nvec));
        while (!IGRAPH_EIT_END(it)) {
            igraph_integer_t edge = IGRAPH_EIT_GET(it);
            igraph_integer_t from, to;
            int ret = 0;
            char *str1, *str2;
            igraph_edge(graph, edge, &from, &to);
            igraph_strvector_get(&nvec, to, &str2);

            if (from == actvertex) {
                ret = fprintf(outstream, "%s\n", str2);
            } else {
                actvertex = from;
                igraph_strvector_get(&nvec, from, &str1);
                ret = fprintf(outstream, "# %s\n%s\n", str1, str2);
            }
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
            IGRAPH_EIT_NEXT(it);
        }
        IGRAPH_FINALLY_CLEAN(1);
    } else if (names == 0) {
        igraph_strvector_t wvec;
        IGRAPH_CHECK(igraph_strvector_init(&wvec, igraph_ecount(graph)));
        IGRAPH_FINALLY(igraph_strvector_destroy, &wvec);
        IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph, weights,
                     igraph_ess_all(IGRAPH_EDGEORDER_ID),
                     &wvec));
        /* No names but weights */
        while (!IGRAPH_EIT_END(it)) {
            igraph_integer_t edge = IGRAPH_EIT_GET(it);
            igraph_integer_t from, to;
            int ret = 0;
            char *str1;
            igraph_edge(graph, edge, &from, &to);
            igraph_strvector_get(&wvec, edge, &str1);
            if (from == actvertex) {
                ret = fprintf(outstream, "%li %s\n", (long)to, str1);
            } else {
                actvertex = from;
                ret = fprintf(outstream, "# %li\n%li %s\n", (long)from, (long)to, str1);
            }
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
            IGRAPH_EIT_NEXT(it);
        }
        igraph_strvector_destroy(&wvec);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        /* Both names and weights */
        igraph_strvector_t nvec, wvec;
        IGRAPH_CHECK(igraph_strvector_init(&wvec, igraph_ecount(graph)));
        IGRAPH_FINALLY(igraph_strvector_destroy, &wvec);
        IGRAPH_CHECK(igraph_strvector_init(&nvec, igraph_vcount(graph)));
        IGRAPH_FINALLY(igraph_strvector_destroy, &nvec);
        IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph, weights,
                     igraph_ess_all(IGRAPH_EDGEORDER_ID),
                     &wvec));
        IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, names,
                     igraph_vss_all(),
                     &nvec));
        while (!IGRAPH_EIT_END(it)) {
            igraph_integer_t edge = IGRAPH_EIT_GET(it);
            igraph_integer_t from, to;
            int ret = 0;
            char *str1, *str2, *str3;
            igraph_edge(graph, edge, &from, &to);
            igraph_strvector_get(&nvec, to, &str2);
            igraph_strvector_get(&wvec, edge, &str3);
            if (from == actvertex) {
                ret = fprintf(outstream, "%s ", str2);
            } else {
                actvertex = from;
                igraph_strvector_get(&nvec, from, &str1);
                ret = fprintf(outstream, "# %s\n%s ", str1, str2);
            }
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
            ret = fprintf(outstream, "%s\n", str3);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
            IGRAPH_EIT_NEXT(it);
        }
        igraph_strvector_destroy(&nvec);
        igraph_strvector_destroy(&wvec);
        IGRAPH_FINALLY_CLEAN(2);
    }

    if (isolates) {
        long int nov = igraph_vcount(graph);
        long int i;
        int ret = 0;
        igraph_vector_t deg;
        igraph_strvector_t nvec;
        char *str;

        IGRAPH_VECTOR_INIT_FINALLY(&deg, 1);
        IGRAPH_CHECK(igraph_strvector_init(&nvec, 1));
        IGRAPH_FINALLY(igraph_strvector_destroy, &nvec);
        for (i = 0; i < nov; i++) {
            igraph_degree(graph, &deg, igraph_vss_1((igraph_integer_t) i),
                          IGRAPH_ALL, IGRAPH_LOOPS);
            if (VECTOR(deg)[0] == 0) {
                if (names == 0) {
                    ret = fprintf(outstream, "# %li\n", i);
                } else {
                    IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, names,
                                 igraph_vss_1((igraph_integer_t) i), &nvec));
                    igraph_strvector_get(&nvec, 0, &str);
                    ret = fprintf(outstream, "# %s\n", str);
                }
            }
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
        }
        igraph_strvector_destroy(&nvec);
        igraph_vector_destroy(&deg);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
