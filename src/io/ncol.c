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

#include "ncol-header.h"

int igraph_ncol_yylex_init_extra (igraph_i_ncol_parsedata_t* user_defined,
                                  void* scanner);
void igraph_ncol_yylex_destroy (void *scanner );
int igraph_ncol_yyparse (igraph_i_ncol_parsedata_t* context);
void igraph_ncol_yyset_in  (FILE * in_str, void* yyscanner );

/**
 * \ingroup loadsave
 * \function igraph_read_graph_ncol
 * \brief Reads a <code>.ncol</code> file used by LGL.
 *
 * Also useful for creating graphs from \quote named\endquote (and
 * optionally weighted) edge lists.
 *
 * </para><para>
 * This format is used by the Large Graph Layout program
 * (http://lgl.sourceforge.net), and it is simply a
 * symbolic weighted edge list. It is a simple text file with one edge
 * per line. An edge is defined by two symbolic vertex names separated
 * by whitespace. (The symbolic vertex names themselves cannot contain
 * whitespace. They might follow by an optional number, this will be
 * the weight of the edge; the number can be negative and can be in
 * scientific notation. If there is no weight specified to an edge it
 * is assumed to be zero.
 *
 * </para><para>
 * The resulting graph is always undirected.
 * LGL cannot deal with files which contain multiple or loop edges,
 * this is however not checked here, as \a igraph is happy with
 * these.
 * \param graph Pointer to an uninitialized graph object.
 * \param instream Pointer to a stream, it should be readable.
 * \param predefnames Pointer to the symbolic names of the vertices in
 *        the file. If \c NULL is given here then vertex ids will be
 *        assigned to vertex names in the order of their appearance in
 *        the \c .ncol file. If it is not \c NULL and some unknown
 *        vertex names are found in the \c .ncol file then new vertex
 *        ids will be assigned to them.
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
 *          problem reading
 *         the file, or the file is syntactically incorrect.
 *
 * Time complexity:
 * O(|V|+|E|log(|V|)) if we neglect
 * the time required by the parsing. As usual
 * |V| is the number of vertices,
 * while |E| is the number of edges.
 *
 * \sa \ref igraph_read_graph_lgl(), \ref igraph_write_graph_ncol()
 */
int igraph_read_graph_ncol(igraph_t *graph, FILE *instream,
                           igraph_strvector_t *predefnames,
                           igraph_bool_t names,
                           igraph_add_weights_t weights,
                           igraph_bool_t directed) {

    igraph_vector_t edges, ws;
    igraph_trie_t trie = IGRAPH_TRIE_NULL;
    igraph_integer_t no_of_nodes;
    long int no_predefined = 0;
    igraph_vector_ptr_t name, weight;
    igraph_vector_ptr_t *pname = 0, *pweight = 0;
    igraph_attribute_record_t namerec, weightrec;
    const char *namestr = "name", *weightstr = "weight";
    igraph_i_ncol_parsedata_t context;

    IGRAPH_CHECK(igraph_empty(graph, 0, directed));
    IGRAPH_FINALLY(igraph_destroy, graph);
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    IGRAPH_TRIE_INIT_FINALLY(&trie, names);
    IGRAPH_VECTOR_INIT_FINALLY(&ws, 0);

    /* Add the predefined names, if any */
    if (predefnames != 0) {
        long int i, id, n;
        char *key;
        n = no_predefined = igraph_strvector_size(predefnames);
        for (i = 0; i < n; i++) {
            igraph_strvector_get(predefnames, i, &key);
            igraph_trie_get(&trie, key, &id);
            if (id != i) {
                IGRAPH_WARNING("reading NCOL file, duplicate entry in predefnames");
                no_predefined--;
            }
        }
    }

    context.has_weights = 0;
    context.vector = &edges;
    context.weights = &ws;
    context.trie = &trie;
    context.eof = 0;

    igraph_ncol_yylex_init_extra(&context, &context.scanner);
    IGRAPH_FINALLY(igraph_ncol_yylex_destroy, context.scanner);

    igraph_ncol_yyset_in(instream, context.scanner);

    if (igraph_ncol_yyparse(&context)) {
        if (context.errmsg[0] != 0) {
            IGRAPH_ERROR(context.errmsg, IGRAPH_PARSEERROR);
        } else {
            IGRAPH_ERROR("Cannot read NCOL file", IGRAPH_PARSEERROR);
        }
    }

    if (predefnames != 0 &&
        igraph_trie_size(&trie) != no_predefined) {
        IGRAPH_WARNING("unknown vertex/vertices found, predefnames extended");
    }

    if (names) {
        const igraph_strvector_t *namevec;
        IGRAPH_CHECK(igraph_vector_ptr_init(&name, 1));
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
        pweight = &weight;
        weightrec.name = weightstr;
        weightrec.type = IGRAPH_ATTRIBUTE_NUMERIC;
        weightrec.value = &ws;
        VECTOR(weight)[0] = &weightrec;
    }

    if (igraph_vector_empty(&edges)) {
        no_of_nodes = 0;
    } else {
        no_of_nodes = igraph_vector_max(&edges) + 1;
    }

    IGRAPH_CHECK(igraph_add_vertices(graph, no_of_nodes, pname));
    IGRAPH_CHECK(igraph_add_edges(graph, &edges, pweight));

    if (pname) {
        igraph_vector_ptr_destroy(pname);
    }
    if (pweight) {
        igraph_vector_ptr_destroy(pweight);
    }
    igraph_vector_destroy(&ws);
    igraph_trie_destroy(&trie);
    igraph_vector_destroy(&edges);
    igraph_ncol_yylex_destroy(context.scanner);
    IGRAPH_FINALLY_CLEAN(5);

    return 0;
}

/**
 * \ingroup loadsave
 * \function igraph_write_graph_ncol
 * \brief Writes the graph to a file in <code>.ncol</code> format
 *
 * </para><para>
 * <code>.ncol</code> is a format used by LGL, see \ref
 * igraph_read_graph_ncol() for details.
 *
 * </para><para>
 * Note that having multiple or loop edges in an
 * <code>.ncol</code> file breaks the  LGL software but
 * \a igraph does not check for this condition.
 * \param graph The graph to write.
 * \param outstream The stream object to write to, it should be
 *        writable.
 * \param names The name of the vertex attribute, if symbolic names
 *        are written to the file. If not, supply 0 here.
 * \param weights The name of the edge attribute, if they are also
 *        written to the file. If you don't want weights, supply 0
 *        here.
 * \return Error code:
 *         \c IGRAPH_EFILE if there is an error writing the
 *         file.
 *
 * Time complexity: O(|E|), the
 * number of edges. All file operations are expected to have time
 * complexity O(1).
 *
 * \sa \ref igraph_read_graph_ncol(), \ref igraph_write_graph_lgl()
 */
int igraph_write_graph_ncol(const igraph_t *graph, FILE *outstream,
                            const char *names, const char *weights) {
    igraph_eit_t it;
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
        if (weighttype != IGRAPH_ATTRIBUTE_NUMERIC) {
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
            ret = fprintf(outstream, "%li %li\n",
                          (long int) from,
                          (long int) to);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
            IGRAPH_EIT_NEXT(it);
        }
    } else if (weights == 0) {
        /* No weights, but use names */
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
            igraph_strvector_get(&nvec, from, &str1);
            igraph_strvector_get(&nvec, to, &str2);
            ret = fprintf(outstream, "%s %s\n", str1, str2);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
            IGRAPH_EIT_NEXT(it);
        }
        igraph_strvector_destroy(&nvec);
        IGRAPH_FINALLY_CLEAN(1);
    } else if (names == 0) {
        /* No names but weights */
        igraph_vector_t wvec;
        IGRAPH_VECTOR_INIT_FINALLY(&wvec, igraph_ecount(graph));
        IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, weights,
                     igraph_ess_all(IGRAPH_EDGEORDER_ID),
                     &wvec));
        while (!IGRAPH_EIT_END(it)) {
            igraph_integer_t edge = IGRAPH_EIT_GET(it);
            igraph_integer_t from, to;
            int ret1, ret2, ret3;
            igraph_edge(graph, edge, &from, &to);
            ret1 = fprintf(outstream, "%li %li ",
                           (long int)from, (long int)to);
            ret2 = igraph_real_fprintf_precise(outstream, VECTOR(wvec)[(long int)edge]);
            ret3 = fputc('\n', outstream);
            if (ret1 < 0 || ret2 < 0 || ret3 == EOF) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
            IGRAPH_EIT_NEXT(it);
        }
        igraph_vector_destroy(&wvec);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        /* Both names and weights */
        igraph_strvector_t nvec;
        igraph_vector_t wvec;
        IGRAPH_VECTOR_INIT_FINALLY(&wvec, igraph_ecount(graph));
        IGRAPH_CHECK(igraph_strvector_init(&nvec, igraph_vcount(graph)));
        IGRAPH_FINALLY(igraph_strvector_destroy, &nvec);
        IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, weights,
                     igraph_ess_all(IGRAPH_EDGEORDER_ID),
                     &wvec));
        IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, names,
                     igraph_vss_all(),
                     &nvec));
        while (!IGRAPH_EIT_END(it)) {
            igraph_integer_t edge = IGRAPH_EIT_GET(it);
            igraph_integer_t from, to;
            int ret = 0, ret2 = 0;
            char *str1, *str2;
            igraph_edge(graph, edge, &from, &to);
            igraph_strvector_get(&nvec, from, &str1);
            igraph_strvector_get(&nvec, to, &str2);
            ret = fprintf(outstream, "%s %s ", str1, str2);
            if (ret < 0) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
            ret = igraph_real_fprintf_precise(outstream, VECTOR(wvec)[(long int)edge]);
            ret2 = fputc('\n', outstream);
            if (ret < 0 || ret2 == EOF) {
                IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
            }
            IGRAPH_EIT_NEXT(it);
        }
        igraph_strvector_destroy(&nvec);
        igraph_vector_destroy(&wvec);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
