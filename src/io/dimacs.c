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

/**
 * \function igraph_read_graph_dimacs
 * \brief Read a graph in DIMACS format.
 *
 * This function reads the DIMACS file format, more specifically the
 * version for network flow problems, see the files at
 * ftp://dimacs.rutgers.edu/pub/netflow/general-info/
 *
 * </para><para>
 * This is a line-oriented text file (ASCII) format. The first
 * character of each line defines the type of the line. If the first
 * character is <code>c</code> the line is a comment line and it is
 * ignored. There is one problem line (<code>p</code> in the file, it
 * must appear before any node and arc descriptor lines. The problem
 * line has three fields separated by spaces: the problem type
 * (<code>min</code>, <code>max</code> or <code>asn</code>), the
 * number of vertices and number of edges in the graph.
 * Exactly two node identification lines are expected
 * (<code>n</code>), one for the source, one for the target vertex.
 * These have two fields: the id of the vertex and the type of the
 * vertex, either <code>s</code> (=source) or <code>t</code>
 * (=target). Arc lines start with <code>a</code> and have three
 * fields: the source vertex, the target vertex and the edge capacity.
 *
 * </para><para>
 * Vertex ids are numbered from 1.
 * \param graph Pointer to an uninitialized graph object.
 * \param instream The file to read from.
 * \param source Pointer to an integer, the id of the source node will
 *    be stored here. (The igraph vertex id, which is one less than
 *    the actual number in the file.) It is ignored if
 *    <code>NULL</code>.
 * \param target Pointer to an integer, the (igraph) id of the target
 *    node will be stored here. It is ignored if <code>NULL</code>.
 * \param capacity Pointer to an initialized vector, the capacity of
 *    the edges will be stored here if not <code>NULL</code>.
 * \param directed Boolean, whether to create a directed graph.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|+c), the number of vertices plus the
 * number of edges, plus the size of the file in characters.
 *
 * \sa \ref igraph_write_graph_dimacs()
 */
int igraph_read_graph_dimacs(igraph_t *graph, FILE *instream,
                             igraph_strvector_t *problem,
                             igraph_vector_t *label,
                             igraph_integer_t *source,
                             igraph_integer_t *target,
                             igraph_vector_t *capacity,
                             igraph_bool_t directed) {

    igraph_vector_t edges;
    long int no_of_nodes = -1;
    long int no_of_edges = -1;
    long int tsource = -1;
    long int ttarget = -1;
    char prob[21];
    char c;
    int problem_type = 0;

#define PROBLEM_EDGE  1
#define PROBLEM_MAX   2

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    if (capacity) {
        igraph_vector_clear(capacity);
    }

    while (!feof(instream)) {
        int read;
        char str[3];

        IGRAPH_ALLOW_INTERRUPTION();

        read = fscanf(instream, "%2c", str);
        if (feof(instream)) {
            break;
        }
        if (read != 1) {
            IGRAPH_ERROR("parsing dimacs file failed", IGRAPH_PARSEERROR);
        }
        switch (str[0]) {
            long int tmp, tmp2;
            long int from, to;
            igraph_real_t cap;

        case 'c':
            /* comment */
            break;

        case 'p':
            if (no_of_nodes != -1) {
                IGRAPH_ERROR("reading dimacs file failed, double 'p' line",
                             IGRAPH_PARSEERROR);
            }
            read = fscanf(instream, "%20s %li %li", prob,
                          &no_of_nodes, &no_of_edges);
            if (read != 3) {
                IGRAPH_ERROR("reading dimacs file failed", IGRAPH_PARSEERROR);
            }
            if (!strcmp(prob, "edge")) {
                /* edge list */
                problem_type = PROBLEM_EDGE;
                if (label) {
                    long int i;
                    IGRAPH_CHECK(igraph_vector_resize(label, no_of_nodes));
                    for (i = 0; i < no_of_nodes; i++) {
                        VECTOR(*label)[i] = i + 1;
                    }
                }
            } else if (!strcmp(prob, "max")) {
                /* maximum flow problem */
                problem_type = PROBLEM_MAX;
                if (capacity) {
                    IGRAPH_CHECK(igraph_vector_reserve(capacity, no_of_edges));
                }
            } else {
                IGRAPH_ERROR("Unknown problem type, should be 'edge' or 'max'",
                             IGRAPH_PARSEERROR);
            }
            if (problem) {
                igraph_strvector_clear(problem);
                IGRAPH_CHECK(igraph_strvector_add(problem, prob));
            }
            IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));
            break;

        case 'n':
            /* for MAX this is either the source or target vertex,
            for EDGE this is a vertex label */
            if (problem_type == PROBLEM_MAX) {
                str[0] = 'x';
                read = fscanf(instream, "%li %1s", &tmp, str);
                if (str[0] == 's') {
                    if (tsource != -1) {
                        IGRAPH_ERROR("reading dimacsfile: multiple source vertex line",
                                     IGRAPH_PARSEERROR);
                    } else {
                        tsource = tmp;
                    }
                } else if (str[0] == 't') {
                    if (ttarget != -1) {
                        IGRAPH_ERROR("reading dimacsfile: multiple target vertex line",
                                     IGRAPH_PARSEERROR);
                    } else {
                        ttarget = tmp;
                    }
                } else {
                    IGRAPH_ERROR("invalid node descriptor line in dimacs file",
                                 IGRAPH_PARSEERROR);
                }
            } else {
                read = fscanf(instream, "%li %li", &tmp, &tmp2);
                if (label) {
                    VECTOR(*label)[tmp] = tmp2;
                }
            }

            break;

        case 'a':
            /* This is valid only for MAX, a weighted edge */
            if (problem_type != PROBLEM_MAX) {
                IGRAPH_ERROR("'a' lines are allowed only in MAX problem files",
                             IGRAPH_PARSEERROR);
            }
            read = fscanf(instream, "%li %li %lf", &from, &to, &cap);
            if (read != 3) {
                IGRAPH_ERROR("reading dimacs file", IGRAPH_PARSEERROR);
            }
            IGRAPH_CHECK(igraph_vector_push_back(&edges, from - 1));
            IGRAPH_CHECK(igraph_vector_push_back(&edges, to - 1));
            if (capacity) {
                IGRAPH_CHECK(igraph_vector_push_back(capacity, cap));
            }
            break;

        case 'e':
            /* Edge line, only in EDGE */
            if (problem_type != PROBLEM_EDGE) {
                IGRAPH_ERROR("'e' lines are allowed only in EDGE problem files",
                             IGRAPH_PARSEERROR);
            }
            read = fscanf(instream, "%li %li", &from, &to);
            if (read != 2) {
                IGRAPH_ERROR("reading dimacs file", IGRAPH_PARSEERROR);
            }
            IGRAPH_CHECK(igraph_vector_push_back(&edges, from - 1));
            IGRAPH_CHECK(igraph_vector_push_back(&edges, to - 1));
            break;

        default:
            IGRAPH_ERROR("unknown line type in dimacs file", IGRAPH_PARSEERROR);
        }

        /* Go to next line */
        while (!feof(instream) && (c = (char) getc(instream)) != '\n') ;
    }

    if (source) {
        *source = (igraph_integer_t) tsource - 1;
    }
    if (target) {
        *target = (igraph_integer_t) ttarget - 1;
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               directed));
    igraph_vector_destroy(&edges);

    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_write_graph_dimacs
 * \brief Write a graph in DIMACS format.
 *
 * This function writes a graph to an output stream in DIMACS format,
 * describing a maximum flow problem.
 * See ftp://dimacs.rutgers.edu/pub/netflow/general-info/
 *
 * </para><para>
 * This file format is discussed in the documentation of \ref
 * igraph_read_graph_dimacs(), see that for more information.
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
 * \sa igraph_read_graph_dimacs()
 */
int igraph_write_graph_dimacs(const igraph_t *graph, FILE *outstream,
                              long int source, long int target,
                              const igraph_vector_t *capacity) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_eit_t it;
    long int i = 0;
    int ret, ret1, ret2, ret3;

    if (igraph_vector_size(capacity) != no_of_edges) {
        IGRAPH_ERROR("invalid capacity vector length", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID),
                                   &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);

    ret = fprintf(outstream,
                  "c created by igraph\np max %li %li\nn %li s\nn %li t\n",
                  no_of_nodes, no_of_edges, source + 1, target + 1);
    if (ret < 0) {
        IGRAPH_ERROR("Write error", IGRAPH_EFILE);
    }


    while (!IGRAPH_EIT_END(it)) {
        igraph_integer_t from, to;
        igraph_real_t cap;
        igraph_edge(graph, IGRAPH_EIT_GET(it), &from, &to);
        cap = VECTOR(*capacity)[i++];
        ret1 = fprintf(outstream, "a %li %li ",
                       (long int) from + 1, (long int) to + 1);
        ret2 = igraph_real_fprintf_precise(outstream, cap);
        ret3 = fputc('\n', outstream);
        if (ret1 < 0 || ret2 < 0 || ret3 == EOF) {
            IGRAPH_ERROR("Write error", IGRAPH_EFILE);
        }
        IGRAPH_EIT_NEXT(it);
    }

    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
