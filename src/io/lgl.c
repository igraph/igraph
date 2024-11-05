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

#include "io/lgl-header.h"
#include "io/parsers/lgl-parser.h"

int igraph_lgl_yylex_init_extra (igraph_i_lgl_parsedata_t *user_defined,
                                 void *scanner);
int igraph_lgl_yylex_destroy(void *scanner);
int igraph_lgl_yyparse(igraph_i_lgl_parsedata_t *context);
void igraph_lgl_yyset_in(FILE *in_str, void *yyscanner);

/* for IGRAPH_FINALLY, which assumes that destructor functions return void */
void igraph_lgl_yylex_destroy_wrapper (void *scanner ) {
    (void) igraph_lgl_yylex_destroy(scanner);
}

/**
 * \ingroup loadsave
 * \function igraph_read_graph_lgl
 * \brief Reads a graph from an <code>.lgl</code> file.
 *
 * The <code>.lgl</code> format is used by the Large Graph
 * Layout visualization software
 * (https://lgl.sourceforge.net), it can
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
 * \param names Logical value, if \c true the symbolic names of the
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
igraph_error_t igraph_read_graph_lgl(igraph_t *graph, FILE *instream,
                          igraph_bool_t names,
                          igraph_add_weights_t weights,
                          igraph_bool_t directed) {

    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t ws = IGRAPH_VECTOR_NULL;
    igraph_trie_t trie = IGRAPH_TRIE_NULL;
    igraph_vector_ptr_t name, weight;
    igraph_vector_ptr_t *pname = 0, *pweight = 0;
    igraph_attribute_record_t namerec, weightrec;
    const char *namestr = "name", *weightstr = "weight";
    igraph_i_lgl_parsedata_t context;

    IGRAPH_VECTOR_INIT_FINALLY(&ws, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_TRIE_INIT_FINALLY(&trie, names);

    context.has_weights = false;
    context.vector = &edges;
    context.weights = &ws;
    context.trie = &trie;
    context.errmsg[0] = '\0';
    context.igraph_errno = IGRAPH_SUCCESS;

    igraph_lgl_yylex_init_extra(&context, &context.scanner);
    IGRAPH_FINALLY(igraph_lgl_yylex_destroy_wrapper, context.scanner);

    igraph_lgl_yyset_in(instream, context.scanner);

    /* Use ENTER/EXIT to avoid destroying context.scanner before this function returns */
    IGRAPH_FINALLY_ENTER();
    int err = igraph_lgl_yyparse(&context);
    IGRAPH_FINALLY_EXIT();
    switch (err) {
    case 0: /* success */
        break;
    case 1: /* parse error */
        if (context.errmsg[0] != '\0') {
            IGRAPH_ERROR(context.errmsg, IGRAPH_PARSEERROR);
        } else if (context.igraph_errno != IGRAPH_SUCCESS) {
            IGRAPH_ERROR("", context.igraph_errno);
        } else {
            IGRAPH_ERROR("Cannot read LGL file.", IGRAPH_PARSEERROR);
        }
        break;
    case 2: /* out of memory */
        IGRAPH_ERROR("Cannot read LGL file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        break;
    default: /* must never reach here */
        /* Hint: This will usually be triggered if an IGRAPH_CHECK() is used in a Bison
         * action instead of an IGRAPH_YY_CHECK(), resulting in an igraph errno being
         * returned in place of a Bison error code.
         * TODO: What if future Bison versions introduce error codes other than 0, 1 and 2?
         */
        IGRAPH_FATALF("Parser returned unexpected error code (%d) when reading LGL file.", err);
    }

    /* Prepare attributes, if needed */

    if (names) {
        IGRAPH_CHECK(igraph_vector_ptr_init(&name, 1));
        IGRAPH_FINALLY(igraph_vector_ptr_destroy, &name);
        pname = &name;
        namerec.name = namestr;
        namerec.type = IGRAPH_ATTRIBUTE_STRING;
        namerec.value = igraph_i_trie_borrow_keys(&trie);
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

    /* Create graph */
    IGRAPH_CHECK(igraph_empty(graph, 0, directed));
    IGRAPH_FINALLY(igraph_destroy, graph);
    IGRAPH_CHECK(igraph_add_vertices(graph, igraph_trie_size(&trie), pname));
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
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&ws);
    igraph_lgl_yylex_destroy(context.scanner);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}


static igraph_error_t check_name(const char *name) {
    size_t len = 0;

    for (; *name != '\0'; name++, len++) {
        if (   *name <= 0x020 /* space or non-printable */
            || *name == 0x7f /* del */
            || *name == '#') {
            IGRAPH_ERRORF("The LGL format does not allow non-printable characters, spaces or '#' in vertex names. "
                          "Character code 0x%02X found.", IGRAPH_EINVAL,
                           *name);
        }
    }
    if (len == 0) {
        IGRAPH_ERROR("The LGL format does not support empty vertex names.", IGRAPH_EINVAL);
    }
    return IGRAPH_SUCCESS;
}


/**
 * \ingroup loadsave
 * \function igraph_write_graph_lgl
 * \brief Writes the graph to a file in <code>.lgl</code> format.
 *
 * <code>.lgl</code> is a format used by LGL, see \ref
 * igraph_read_graph_lgl() for details.
 *
 * </para><para>
 * Note that having multiple or loop edges in an
 * <code>.lgl</code> file breaks the  LGL software but \a igraph
 * does not check for this condition.
 *
 * \param graph The graph to write.
 * \param outstream The stream object to write to, it should be
 *        writable.
 * \param names The name of a string vertex attribute, if symbolic names
 *        are to be written to the file. Supply \c NULL to write vertex
 *        ids instead.
 * \param weights The name of a numerical edge attribute, which will be
 *        written as weights to the file. Supply \c NULL to skip writing
 *        edge weights.
 * \param isolates Logical, if \c true isolated vertices are also written
 *        to the file. If \c false they will be omitted.
 * \return Error code:
 *         \c IGRAPH_EFILE if there is an error
 *         writing the file.
 *
 * Time complexity: O(|E|), the number of edges if \p isolates is \c false,
 * O(|V|+|E|) otherwise. All file operations are expected to have
 * time complexity O(1).
 *
 * \sa \ref igraph_read_graph_lgl(), \ref igraph_write_graph_ncol()
 *
 * \example examples/simple/igraph_write_graph_lgl.c
 */
igraph_error_t igraph_write_graph_lgl(const igraph_t *graph, FILE *outstream,
                           const char *names, const char *weights,
                           igraph_bool_t isolates) {
    igraph_eit_t it;
    igraph_integer_t actvertex = -1;
    igraph_attribute_type_t nametype, weighttype;
    const igraph_integer_t vcount = igraph_vcount(graph), ecount = igraph_ecount(graph);

    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_FROM),
                                   &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);

    /* Check if we have the names attribute */
    if (names && !igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX,
            names)) {
        IGRAPH_WARNINGF("Names attribute '%s' does not exist.", names);
        names = NULL;
    }
    if (names) {
        IGRAPH_CHECK(igraph_i_attribute_gettype(graph, &nametype,
                                                IGRAPH_ATTRIBUTE_VERTEX, names));
        if (nametype != IGRAPH_ATTRIBUTE_STRING) {
            IGRAPH_WARNINGF("Ignoring names attribute '%s', unknown attribute type.", names);
            names = NULL;
        }
    }

    /* Check the weights as well */
    if (weights && !igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE, weights)) {
        IGRAPH_WARNINGF("Weights attribute '%s' does not exist.", weights);
        weights = NULL;
    }
    if (weights) {
        IGRAPH_CHECK(igraph_i_attribute_gettype(graph, &weighttype,
                                                IGRAPH_ATTRIBUTE_EDGE, weights));
        if (weighttype != IGRAPH_ATTRIBUTE_NUMERIC) {
            IGRAPH_WARNINGF("Ignoring weights attribute '%s', unknown attribute type.", weights);
            weights = NULL;
        }
    }

    if (names == NULL && weights == NULL) {
        /* No names, no weights */
        while (!IGRAPH_EIT_END(it)) {
            igraph_integer_t from, to;
            int ret;
            igraph_edge(graph, IGRAPH_EIT_GET(it), &from, &to);
            if (from == actvertex) {
                ret = fprintf(outstream, "%" IGRAPH_PRId "\n", to);
            } else {
                actvertex = from;
                ret = fprintf(outstream, "# %" IGRAPH_PRId "\n%" IGRAPH_PRId "\n", from, to);
            }
            if (ret < 0) {
                IGRAPH_ERROR("Writing LGL file failed.", IGRAPH_EFILE);
            }
            IGRAPH_EIT_NEXT(it);
        }
    } else if (weights == NULL) {
        /* No weights but use names */
        igraph_strvector_t nvec;
        IGRAPH_STRVECTOR_INIT_FINALLY(&nvec, vcount);
        IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, names,
                     igraph_vss_all(),
                     &nvec));
        while (!IGRAPH_EIT_END(it)) {
            igraph_integer_t edge = IGRAPH_EIT_GET(it);
            igraph_integer_t from, to;
            int ret = 0;
            const char *str1, *str2;
            igraph_edge(graph, edge, &from, &to);
            str2 = igraph_strvector_get(&nvec, to);
            IGRAPH_CHECK(check_name(str2));

            if (from == actvertex) {
                ret = fprintf(outstream, "%s\n", str2);
            } else {
                actvertex = from;
                str1 = igraph_strvector_get(&nvec, from);
                IGRAPH_CHECK(check_name(str1));
                ret = fprintf(outstream, "# %s\n%s\n", str1, str2);
            }
            if (ret < 0) {
                IGRAPH_ERROR("Writing LGL file failed.", IGRAPH_EFILE);
            }
            IGRAPH_EIT_NEXT(it);
        }
        igraph_strvector_destroy(&nvec);
        IGRAPH_FINALLY_CLEAN(1);
    } else if (names == NULL) {
        /* No names but weights */
        igraph_vector_t wvec;
        IGRAPH_VECTOR_INIT_FINALLY(&wvec, ecount);
        IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, weights,
                     igraph_ess_all(IGRAPH_EDGEORDER_ID),
                     &wvec));
        while (!IGRAPH_EIT_END(it)) {
            igraph_integer_t edge = IGRAPH_EIT_GET(it);
            igraph_integer_t from, to;
            int ret1, ret2, ret3;
            igraph_edge(graph, edge, &from, &to);
            if (from == actvertex) {
                ret1 = fprintf(outstream, "%" IGRAPH_PRId " ", to);
            } else {
                actvertex = from;
                ret1 = fprintf(outstream, "# %" IGRAPH_PRId "\n%" IGRAPH_PRId " ", from, to);
            }
            ret2 = igraph_real_fprintf_precise(outstream, VECTOR(wvec)[edge]);
            ret3 = fputc('\n', outstream);
            if (ret1 < 0 || ret2 < 0 || ret3 == EOF) {
                IGRAPH_ERROR("Writing LGL file failed.", IGRAPH_EFILE);
            }
            IGRAPH_EIT_NEXT(it);
        }
        igraph_vector_destroy(&wvec);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        /* Both names and weights */
        igraph_strvector_t nvec;
        igraph_vector_t wvec;
        IGRAPH_VECTOR_INIT_FINALLY(&wvec, ecount);
        IGRAPH_STRVECTOR_INIT_FINALLY(&nvec, vcount);
        IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, weights,
                     igraph_ess_all(IGRAPH_EDGEORDER_ID),
                     &wvec));
        IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, names,
                     igraph_vss_all(),
                     &nvec));
        while (!IGRAPH_EIT_END(it)) {
            igraph_integer_t edge = IGRAPH_EIT_GET(it);
            igraph_integer_t from, to;
            int ret = 0, ret2;
            const char *str1, *str2;
            igraph_edge(graph, edge, &from, &to);
            str2 = igraph_strvector_get(&nvec, to);
            IGRAPH_CHECK(check_name(str2));

            if (from == actvertex) {
                ret = fprintf(outstream, "%s ", str2);
            } else {
                actvertex = from;
                str1 = igraph_strvector_get(&nvec, from);
                IGRAPH_CHECK(check_name(str1));
                ret = fprintf(outstream, "# %s\n%s ", str1, str2);
            }
            if (ret < 0) {
                IGRAPH_ERROR("Writing LGL file failed.", IGRAPH_EFILE);
            }
            ret = igraph_real_fprintf_precise(outstream, VECTOR(wvec)[edge]);
            ret2 = fputc('\n', outstream);
            if (ret < 0 || ret2 == EOF) {
                IGRAPH_ERROR("Writing LGL file failed.", IGRAPH_EFILE);
            }
            IGRAPH_EIT_NEXT(it);
        }
        igraph_strvector_destroy(&nvec);
        igraph_vector_destroy(&wvec);
        IGRAPH_FINALLY_CLEAN(2);
    }

    if (isolates) {
        igraph_integer_t nov = vcount;
        igraph_integer_t i;
        int ret = 0;
        igraph_integer_t deg;
        igraph_strvector_t nvec;
        const char *str;

        IGRAPH_STRVECTOR_INIT_FINALLY(&nvec, 1);
        for (i = 0; i < nov; i++) {
            IGRAPH_CHECK(igraph_degree_1(graph, &deg, i, IGRAPH_ALL, IGRAPH_LOOPS));
            if (deg == 0) {
                if (names == NULL) {
                    ret = fprintf(outstream, "# %" IGRAPH_PRId "\n", i);
                } else {
                    IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, names,
                                 igraph_vss_1(i), &nvec));
                    str = igraph_strvector_get(&nvec, 0);
                    IGRAPH_CHECK(check_name(str));
                    ret = fprintf(outstream, "# %s\n", str);
                }
            }
            if (ret < 0) {
                IGRAPH_ERROR("Writing LGL file failed.", IGRAPH_EFILE);
            }
        }
        igraph_strvector_destroy(&nvec);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
