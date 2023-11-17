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

#include "io/dl-header.h"
#include "io/parsers/dl-parser.h"

int igraph_dl_yylex_init_extra (igraph_i_dl_parsedata_t* user_defined,
                                void* scanner);
int igraph_dl_yylex_destroy (void *scanner );
int igraph_dl_yyparse (igraph_i_dl_parsedata_t* context);
void igraph_dl_yyset_in  (FILE * in_str, void* yyscanner );

/* for IGRAPH_FINALLY, which assumes that destructor functions return void */
void igraph_dl_yylex_destroy_wrapper (void *scanner ) {
    (void) igraph_dl_yylex_destroy(scanner);
}

/**
 * \function igraph_read_graph_dl
 * \brief Reads a file in the DL format of UCINET.
 *
 * This is a simple textual file format used by UCINET. See
 * http://www.analytictech.com/networks/dataentry.htm for
 * examples. All the forms described here are supported by
 * igraph. Vertex names and edge weights are also supported and they
 * are added as attributes. (If an attribute handler is attached.)
 *
 * </para><para> Note the specification does not mention whether the
 * format is case sensitive or not. For igraph DL files are case
 * sensitive, i.e. \c Larry and \c larry are not the same.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param instream The stream to read the DL file from.
 * \param directed Logical scalar, whether to create a directed file.
 * \return Error code.
 *
 * Time complexity: linear in terms of the number of edges and
 * vertices, except for the matrix format, which is quadratic in the
 * number of vertices.
 *
 * \example examples/simple/igraph_read_graph_dl.c
 */

igraph_error_t igraph_read_graph_dl(igraph_t *graph, FILE *instream,
                         igraph_bool_t directed) {

    igraph_integer_t n, n2;
    const igraph_strvector_t *namevec = 0;
    igraph_vector_ptr_t name, weight;
    igraph_vector_ptr_t *pname = 0, *pweight = 0;
    igraph_attribute_record_t namerec, weightrec;
    const char *namestr = "name", *weightstr = "weight";
    igraph_i_dl_parsedata_t context;

    context.eof = 0;
    context.mode = 0;
    context.n = -1;
    context.from = 0;
    context.to = 0;
    context.errmsg[0] = '\0';
    context.igraph_errno = IGRAPH_SUCCESS;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&context.edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&context.weights, 0);
    IGRAPH_CHECK(igraph_strvector_init(&context.labels, 0));
    IGRAPH_FINALLY(igraph_strvector_destroy, &context.labels);
    IGRAPH_TRIE_INIT_FINALLY(&context.trie, /*names=*/ 1);

    igraph_dl_yylex_init_extra(&context, &context.scanner);
    IGRAPH_FINALLY(igraph_dl_yylex_destroy_wrapper, context.scanner);

    igraph_dl_yyset_in(instream, context.scanner);

    /* Use ENTER/EXIT to avoid destroying context.scanner before this function returns */
    IGRAPH_FINALLY_ENTER();
    int err = igraph_dl_yyparse(&context);
    IGRAPH_FINALLY_EXIT();
    switch (err) {
    case 0: /* success */
        break;
    case 1: /* parse error */
        if (context.errmsg[0] != 0) {
            IGRAPH_ERROR(context.errmsg, IGRAPH_PARSEERROR);
        } else if (context.igraph_errno != IGRAPH_SUCCESS) {
            IGRAPH_ERROR("", context.igraph_errno);
        } else {
            IGRAPH_ERROR("Cannot read DL file.", IGRAPH_PARSEERROR);
        }
        break;
    case 2: /* out of memory */
        IGRAPH_ERROR("Cannot read DL file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        break;
    default: /* must never reach here */
        /* Hint: This will usually be triggered if an IGRAPH_CHECK() is used in a Bison
         * action instead of an IGRAPH_YY_CHECK(), resulting in an igraph errno being
         * returned in place of a Bison error code.
         * TODO: What if future Bison versions introduce error codes other than 0, 1 and 2?
         */
        IGRAPH_FATALF("Parser returned unexpected error code (%d) when reading DL file.", err);
    }

    /* Extend the weight vector, if needed */
    n = igraph_vector_size(&context.weights);
    n2 = igraph_vector_int_size(&context.edges) / 2;
    if (n != 0) {
        IGRAPH_CHECK(igraph_vector_resize(&context.weights, n2));
        for (; n < n2; n++) {
            VECTOR(context.weights)[n] = IGRAPH_NAN;
        }
    }

    /* Check number of vertices */
    if (n2 > 0) {
        n = igraph_vector_int_max(&context.edges);
    } else {
        n = 0;
    }
    if (n >= context.n) {
        IGRAPH_WARNING("More vertices than specified in `DL' file");
        context.n = n;
    }

    /* Prepare attributes */

    /* Labels */
    if (igraph_strvector_size(&context.labels) != 0) {
        namevec = (const igraph_strvector_t*) &context.labels;
    } else if (igraph_trie_size(&context.trie) != 0) {
        namevec = igraph_i_trie_borrow_keys(&context.trie);
    }
    if (namevec) {
        IGRAPH_CHECK(igraph_vector_ptr_init(&name, 1));
        IGRAPH_FINALLY(igraph_vector_ptr_destroy, &name);
        pname = &name;
        namerec.name = namestr;
        namerec.type = IGRAPH_ATTRIBUTE_STRING;
        namerec.value = namevec;
        VECTOR(name)[0] = &namerec;
    }

    /* Weights */
    if (igraph_vector_size(&context.weights) != 0) {
        IGRAPH_CHECK(igraph_vector_ptr_init(&weight, 1));
        IGRAPH_FINALLY(igraph_vector_ptr_destroy, &weight);
        pweight = &weight;
        weightrec.name = weightstr;
        weightrec.type = IGRAPH_ATTRIBUTE_NUMERIC;
        weightrec.value = &context.weights;
        VECTOR(weight)[0] = &weightrec;
    }

    /* Create graph */
    IGRAPH_CHECK(igraph_empty(graph, 0, directed));
    IGRAPH_FINALLY(igraph_destroy, graph);
    IGRAPH_CHECK(igraph_add_vertices(graph, context.n, pname));
    IGRAPH_CHECK(igraph_add_edges(graph, &context.edges, pweight));

    if (pweight) {
        igraph_vector_ptr_destroy(pweight);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (pname) {
        igraph_vector_ptr_destroy(pname);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* don't destroy the graph itself but pop it from the finally stack */
    IGRAPH_FINALLY_CLEAN(1);

    igraph_trie_destroy(&context.trie);
    igraph_strvector_destroy(&context.labels);
    igraph_vector_int_destroy(&context.edges);
    igraph_vector_destroy(&context.weights);
    igraph_dl_yylex_destroy(context.scanner);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}
