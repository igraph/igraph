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
#include "igraph_interface.h"
#include "igraph_iterators.h"

#include "core/interruption.h"
#include "io/parse_utils.h"

/**
 * \section about_loadsave
 *
 * <para>These functions can write a graph to a file, or read a graph
 * from a file.</para>
 *
 * <para>They assume that the current locale uses a decimal point and not
 * a decimal comma. See \ref igraph_enter_safelocale() and
 * \ref igraph_exit_safelocale() for more information.</para>
 *
 * <para>Note that as \a igraph uses the traditional C streams, it is
 * possible to read/write files from/to memory, at least on GNU
 * operating systems supporting \quote non-standard\endquote streams.</para>
 */

/**
 * \ingroup loadsave
 * \function igraph_read_graph_edgelist
 * \brief Reads an edge list from a file and creates a graph.
 *
 * This format is simply a series of an even number of non-negative integers separated by
 * whitespace. The integers represent vertex IDs. Placing each edge (i.e. pair of integers)
 * on a separate line is not required, but it is recommended for readability.
 * Edges of directed graphs are assumed to be in "from, to" order.
 *
 * </para><para>
 * The largest vertex ID plus one, or the parameter \p n determines the vertex count,
 * whichever is larger. See \ref igraph_read_graph_ncol() for reading files where
 * vertices are specified by name instead of by a numerical vertex ID.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param instream Pointer to a stream, it should be readable.
 * \param n The number of vertices in the graph. If smaller than the
 *        largest integer in the file it will be ignored. It is thus
 *        safe to supply zero here.
 * \param directed Logical, if true the graph is directed, if false it
 *        will be undirected.
 * \return Error code:
 *         \c IGRAPH_PARSEERROR: if there is a
 *         problem reading the file, or the file is syntactically
 *         incorrect.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges. It is assumed that
 * reading an integer requires O(1) time.
 */
igraph_error_t igraph_read_graph_edgelist(igraph_t *graph, FILE *instream,
                               igraph_integer_t n, igraph_bool_t directed) {

    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_integer_t from, to;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, 100));

    for (;;) {
        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_CHECK(igraph_i_fskip_whitespace(instream));

        if (feof(instream)) break;

        IGRAPH_CHECK(igraph_i_fget_integer(instream, &from));
        IGRAPH_CHECK(igraph_i_fget_integer(instream, &to));

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
        /* Protect from very large memory allocations when fuzzing. */
#define IGRAPH_EDGELIST_MAX_VERTEX_COUNT (1L << 20)
        if (from > IGRAPH_EDGELIST_MAX_VERTEX_COUNT || to > IGRAPH_EDGELIST_MAX_VERTEX_COUNT) {
            IGRAPH_ERROR("Vertex count too large in edgelist file.", IGRAPH_EINVAL);
        }
#endif

        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup loadsave
 * \function igraph_write_graph_edgelist
 * \brief Writes the edge list of a graph to a file.
 *
 * </para><para>
 * Edges are represented as pairs of 0-based vertex indices.
 * One edge is written per line, separated by a single space.
 * For directed graphs edges are written in from, to order.
 *
 * \param graph The graph object to write.
 * \param outstream Pointer to a stream, it should be writable.
 * \return Error code:
 *         \c IGRAPH_EFILE if there is an error writing the
 *         file.
 *
 * Time complexity: O(|E|), the
 * number of edges in the  graph. It is assumed that writing an
 * integer to the file requires O(1)
 * time.
 */
igraph_error_t igraph_write_graph_edgelist(const igraph_t *graph, FILE *outstream) {

    igraph_eit_t it;

    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_FROM),
                                   &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);

    while (!IGRAPH_EIT_END(it)) {
        igraph_integer_t from, to;
        int ret;
        igraph_edge(graph, IGRAPH_EIT_GET(it), &from, &to);
        ret = fprintf(outstream, "%" IGRAPH_PRId " %" IGRAPH_PRId "\n",
                      from,
                      to);
        if (ret < 0) {
            IGRAPH_ERROR("Failed writing edgelist.", IGRAPH_EFILE);
        }
        IGRAPH_EIT_NEXT(it);
    }

    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
