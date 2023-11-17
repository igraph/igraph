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

#include <string.h>

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
/* Limit maximum vertex count when using a fuzzer, to avoid out-of-memory failure. */
#define IGRAPH_DIMACS_MAX_VERTEX_COUNT (1 << 20)
#define IGRAPH_DIMACS_MAX_EDGE_COUNT   (1 << 20)
#else
#define IGRAPH_DIMACS_MAX_VERTEX_COUNT INT32_MAX
#define IGRAPH_DIMACS_MAX_EDGE_COUNT   INT32_MAX
#endif

/**
 * \function igraph_read_graph_dimacs
 * \brief Read a graph in DIMACS format (deprecated alias).
 *
 * \deprecated-by igraph_read_graph_dimacs_flow 0.10.0
 */
igraph_error_t igraph_read_graph_dimacs(igraph_t *graph, FILE *instream,
                             igraph_strvector_t *problem,
                             igraph_vector_int_t *label,
                             igraph_integer_t *source,
                             igraph_integer_t *target,
                             igraph_vector_t *capacity,
                             igraph_bool_t directed) {
    return igraph_read_graph_dimacs_flow(
        graph, instream, problem, label, source, target, capacity, directed
    );
}

#define EXPECT(actual, expected) \
    do { \
        if ((actual) != (expected)) { \
            IGRAPH_ERROR("Reading DIMACS flow problem file failed.", IGRAPH_PARSEERROR); \
        } \
    } while (0)

#define CHECK_VID(vid) \
    do { \
        if (vid > IGRAPH_DIMACS_MAX_VERTEX_COUNT) { \
            IGRAPH_ERRORF("Vertex ID %" IGRAPH_PRId " too large in DIMACS file.", IGRAPH_PARSEERROR, vid); \
        } \
    } while(0)

/**
 * \function igraph_read_graph_dimacs_flow
 * \brief Read a graph in DIMACS format.
 *
 * This function reads the DIMACS file format, more specifically the
 * version for network flow problems, see the files at
 * http://archive.dimacs.rutgers.edu/pub/netflow/general-info/
 *
 * </para><para>
 * This is a line-oriented text file (ASCII) format. The first
 * character of each line defines the type of the line. If the first
 * character is \c c the line is a comment line and it is
 * ignored. There is one problem line (\c p in the file), it
 * must appear before any node and arc descriptor lines. The problem
 * line has three fields separated by spaces: the problem type
 * (\c max or \c edge), the number of vertices,
 * and number of edges in the graph. In MAX problems,
 * exactly two node identification lines are expected
 * (\c n), one for the source, and one for the target vertex.
 * These have two fields: the ID of the vertex and the type of the
 * vertex, either \c s ( = source) or \c t ( = target).
 * Arc lines start with \c a and have three fields: the source vertex,
 * the target vertex and the edge capacity. In EDGE problems,
 * there may be a node line (\c n) for each node. It specifies the
 * node index and an integer node label. Nodes for which no explicit
 * label was specified will use their index as label. In EDGE problems,
 * each edge is specified as an edge line (\c e).
 *
 * </para><para>
 * Within DIMACS files, vertex IDs are numbered from 1.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param instream The file to read from.
 * \param problem If not \c NULL, it will contain the problem type.
 * \param label If not \c NULL, node labels will be stored here for \c edge
 *    problems. Ignored for \c max problems.
 * \param source Pointer to an integer, the ID of the source node will
 *    be stored here. (The igraph vertex ID, which is one less than
 *    the actual number in the file.) It is ignored if \c NULL.
 * \param target Pointer to an integer, the (igraph) ID of the target
 *    node will be stored here. It is ignored if \c NULL.
 * \param capacity Pointer to an initialized vector, the capacity of
 *    the edges will be stored here if not \ NULL.
 * \param directed Boolean, whether to create a directed graph.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|+c), the number of vertices plus the
 * number of edges, plus the size of the file in characters.
 *
 * \sa \ref igraph_write_graph_dimacs()
 */
igraph_error_t igraph_read_graph_dimacs_flow(
        igraph_t *graph, FILE *instream,
        igraph_strvector_t *problem,
        igraph_vector_int_t *label,
        igraph_integer_t *source,
        igraph_integer_t *target,
        igraph_vector_t *capacity,
        igraph_bool_t directed) {

    igraph_vector_int_t edges;
    igraph_integer_t no_of_nodes = -1;
    igraph_integer_t no_of_edges = -1;
    igraph_integer_t tsource = -1;
    igraph_integer_t ttarget = -1;
    char prob[21];
    enum {
        PROBLEM_NONE,
        PROBLEM_EDGE,
        PROBLEM_MAX
    } problem_type = PROBLEM_NONE;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    if (capacity) {
        igraph_vector_clear(capacity);
    }

    while (!feof(instream)) {
        int read;
        char str[2];

        IGRAPH_ALLOW_INTERRUPTION();

        read = fscanf(instream, "%2c", str);
        if (feof(instream)) {
            break;
        }
        EXPECT(read, 1);
        switch (str[0]) {
            igraph_integer_t tmp, tmp2;
            igraph_integer_t from, to;
            igraph_real_t cap;

        case 'c':
            /* comment */
            break;

        case 'p':
            if (no_of_nodes != -1) {
                IGRAPH_ERROR("Reading DIMACS file failed, double 'p' line.",
                             IGRAPH_PARSEERROR);
            }
            read = fscanf(instream, "%20s %" IGRAPH_PRId " %" IGRAPH_PRId "", prob,
                          &no_of_nodes, &no_of_edges);
            EXPECT(read, 3);
            if (no_of_nodes > IGRAPH_DIMACS_MAX_VERTEX_COUNT) {
                IGRAPH_ERROR("Vertex count too large in DIMACS file.", IGRAPH_PARSEERROR);
            }
            if (no_of_nodes < 0) {
                IGRAPH_ERROR("Invalid (negative) vertex count in DIMACS file.", IGRAPH_PARSEERROR);
            }
            if (no_of_edges > IGRAPH_DIMACS_MAX_EDGE_COUNT) {
                IGRAPH_ERROR("Edge count too large in DIMACS file.", IGRAPH_PARSEERROR);
            }
            if (no_of_edges < 0) {
                IGRAPH_ERROR("Invalid (negative) edge count in DIMACS file.", IGRAPH_PARSEERROR);
            }
            if (!strcmp(prob, "edge")) {
                /* edge list */
                problem_type = PROBLEM_EDGE;
                if (label) {
                    IGRAPH_CHECK(igraph_vector_int_range(label, 1, no_of_nodes+1));
                }
            } else if (!strcmp(prob, "max")) {
                /* maximum flow problem */
                problem_type = PROBLEM_MAX;
                if (capacity) {
                    IGRAPH_CHECK(igraph_vector_reserve(capacity, no_of_edges));
                }
            } else {
                IGRAPH_ERROR("Unknown problem type, should be 'edge' or 'max'.",
                             IGRAPH_PARSEERROR);
            }
            if (problem) {
                igraph_strvector_clear(problem);
                IGRAPH_CHECK(igraph_strvector_push_back(problem, prob));
            }
            IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges * 2));
            break;

        case 'n':
            /* for MAX this is either the source or target vertex,
            for EDGE this is a vertex label */
            if (problem_type == PROBLEM_MAX) {
                str[0] = 'x';
                read = fscanf(instream, "%" IGRAPH_PRId " %1s", &tmp, str);
                EXPECT(read, 2);
                if (str[0] == 's') {
                    if (tsource != -1) {
                        IGRAPH_ERROR("Reading DIMACS file: multiple source vertex line.",
                                     IGRAPH_PARSEERROR);
                    } else {
                        tsource = tmp;
                    }
                } else if (str[0] == 't') {
                    if (ttarget != -1) {
                        IGRAPH_ERROR("Reading DIMACS file: multiple target vertex line.",
                                     IGRAPH_PARSEERROR);
                    } else {
                        ttarget = tmp;
                    }
                } else {
                    IGRAPH_ERROR("Invalid node descriptor line in DIMACS file.",
                                 IGRAPH_PARSEERROR);
                }
            } else { /* PROBLEM_EDGE */
                read = fscanf(instream, "%" IGRAPH_PRId " %" IGRAPH_PRId "", &tmp, &tmp2);
                EXPECT(read, 1);
                if (label) {
                    if (tmp < 0 || tmp >= no_of_nodes) {
                        IGRAPH_ERRORF("Invalid node index %" IGRAPH_PRId " in DIMACS file. "
                                      "Number of nodes was given as %" IGRAPH_PRId".",
                                      IGRAPH_PARSEERROR, tmp, no_of_nodes);
                    }
                    VECTOR(*label)[tmp] = tmp2;
                }
            }

            break;

        case 'a':
            /* This is valid only for MAX, a weighted edge */
            if (problem_type != PROBLEM_MAX) {
                IGRAPH_ERROR("'a' lines are allowed only in MAX problem files.",
                             IGRAPH_PARSEERROR);
            }
            read = fscanf(instream, "%" IGRAPH_PRId " %" IGRAPH_PRId " %lf", &from, &to, &cap);
            EXPECT(read, 3);
            CHECK_VID(from);
            CHECK_VID(to);
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from - 1));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to - 1));
            if (capacity) {
                IGRAPH_CHECK(igraph_vector_push_back(capacity, cap));
            }
            break;

        case 'e':
            /* Edge line, only in EDGE */
            if (problem_type != PROBLEM_EDGE) {
                IGRAPH_ERROR("'e' lines are allowed only in EDGE problem files.",
                             IGRAPH_PARSEERROR);
            }
            read = fscanf(instream, "%" IGRAPH_PRId " %" IGRAPH_PRId "", &from, &to);
            EXPECT(read, 2);
            CHECK_VID(from);
            CHECK_VID(to);
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from - 1));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to - 1));
            break;

        default:
            IGRAPH_ERROR("Unknown line type in DIMACS file.", IGRAPH_PARSEERROR);
        }

        /* Go to next line */
        while (!feof(instream) && getc(instream) != '\n') ;
    }

    if (source) {
        *source = tsource - 1;
    }
    if (target) {
        *target = ttarget - 1;
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_write_graph_dimacs
 * \brief Write a graph in DIMACS format (deprecated alias).
 *
 * \deprecated-by igraph_write_graph_dimacs_flow 0.10.0
 */
igraph_error_t igraph_write_graph_dimacs(const igraph_t *graph, FILE *outstream,
                              igraph_integer_t source, igraph_integer_t target,
                              const igraph_vector_t *capacity) {
    return igraph_write_graph_dimacs_flow(graph, outstream, source, target, capacity);
}

/**
 * \function igraph_write_graph_dimacs_flow
 * \brief Write a graph in DIMACS format.
 *
 * This function writes a graph to an output stream in DIMACS format,
 * describing a maximum flow problem.
 * See ftp://dimacs.rutgers.edu/pub/netflow/general-info/
 *
 * </para><para>
 * This file format is discussed in the documentation of \ref
 * igraph_read_graph_dimacs_flow(), see that for more information.
 *
 * \param graph The graph to write to the stream.
 * \param outstream The stream.
 * \param source Integer, the id of the source vertex for the maximum
 *     flow.
 * \param target Integer, the id of the target vertex.
 * \param capacity Pointer to an initialized vector containing the
 *     edge capacity values.
 * \return Error code.
 *
 * Time complexity: O(|E|), the number of edges in the graph.
 *
 * \sa \ref igraph_read_graph_dimacs_flow()
 */
igraph_error_t igraph_write_graph_dimacs_flow(const igraph_t *graph, FILE *outstream,
                              igraph_integer_t source, igraph_integer_t target,
                              const igraph_vector_t *capacity) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_eit_t it;
    igraph_integer_t i = 0;
    int ret, ret1, ret2, ret3;

    if (igraph_vector_size(capacity) != no_of_edges) {
        IGRAPH_ERROR("invalid capacity vector length", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID),
                                   &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);

    ret = fprintf(outstream,
                  "c created by igraph\np max %" IGRAPH_PRId " %" IGRAPH_PRId "\nn %" IGRAPH_PRId " s\nn %" IGRAPH_PRId " t\n",
                  no_of_nodes, no_of_edges, source + 1, target + 1);
    if (ret < 0) {
        IGRAPH_ERROR("Write error", IGRAPH_EFILE);
    }


    while (!IGRAPH_EIT_END(it)) {
        igraph_integer_t from, to;
        igraph_real_t cap;
        igraph_edge(graph, IGRAPH_EIT_GET(it), &from, &to);
        cap = VECTOR(*capacity)[i++];
        ret1 = fprintf(outstream, "a %" IGRAPH_PRId " %" IGRAPH_PRId " ",
                       from + 1, to + 1);
        ret2 = igraph_real_fprintf_precise(outstream, cap);
        ret3 = fputc('\n', outstream);
        if (ret1 < 0 || ret2 < 0 || ret3 == EOF) {
            IGRAPH_ERROR("Write error", IGRAPH_EFILE);
        }
        IGRAPH_EIT_NEXT(it);
    }

    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
