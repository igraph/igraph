/*
   IGraph library.
   Copyright (C) 2023  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "igraph_interface.h"
#include "igraph_foreign.h"

//TODO remove unnecessary
#include "igraph_components.h"
#include "igraph_constants.h"
#include "igraph_constructors.h"
#include "igraph_datatype.h"
#include "igraph_memory.h"
#include "igraph_structural.h"
#include "igraph_types.h"
//end TODO

#include <stdio.h>
#include <string.h>

/**
 * \function igraph_read_graph6
 * \brief Read a graph in graph6.
 *
 * Graph6 is a format that uses printable characters, see
 * http://users.cecs.anu.edu.au/~bdm/data/formats.txt
 * for details.
 *
 * </para><para>
 * \param graph Pointer to an uninitialized graph object.
 * \param instream The stream to read the graph6 file from.
 * \return Error code.
 *
 * Time complexity: should be proportional to the length of the file.
 *
 * \ref igraph_write_graph6() for writing GML files.
 *
 * \example examples/simple/graph6.c
 */

static igraph_integer_t igraph_read_graph6_number_and_skip(const char *str) {
    //TODO bigger numbers
    igraph_integer_t res = str[0] - 63;
    str++;
    return res;
}

static igraph_error_t igraph_from_nauty_digraph6(igraph_t *graph, const char *str) {

    igraph_integer_t no_of_nodes = igraph_read_graph6_number_and_skip(str);

    igraph_integer_t length = strlen(str);
    igraph_matrix_t adj_matrix;

    IGRAPH_CHECK(igraph_matrix_init(&adj_matrix, no_of_nodes, no_of_nodes));
    IGRAPH_FINALLY(igraph_matrix_destroy, &adj_matrix);

    igraph_integer_t bitIndex = 0;
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
            //TODO bigger numbers
            igraph_integer_t byte = 1 + bitIndex / 6;
            if (byte >= length) break;
            igraph_integer_t bit = 5 - bitIndex % 6;
            igraph_integer_t value = str[byte] - 63;
            igraph_integer_t is_set = value & (1 << bit);
            if (is_set) {
                MATRIX(adj_matrix, i, j) = 1;
            }
            bitIndex++;
        }
    }
    igraph_adjacency(graph, &adj_matrix, IGRAPH_ADJ_DIRECTED, IGRAPH_LOOPS_ONCE);
    igraph_matrix_destroy(&adj_matrix);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_from_nauty_graph6(igraph_t *graph, const char *str) {
    igraph_integer_t no_of_nodes = igraph_read_graph6_number_and_skip(str);

    igraph_integer_t length = strlen(str);
    igraph_matrix_t adj_matrix;

    IGRAPH_CHECK(igraph_matrix_init(&adj_matrix, no_of_nodes, no_of_nodes));
    IGRAPH_FINALLY(igraph_matrix_destroy, &adj_matrix);

    igraph_integer_t bitIndex = 0;
    for (igraph_integer_t j = 1; j < no_of_nodes; j++) {
        for (igraph_integer_t i = 0; i < j; i++) {
            //TODO bigger numbers
            igraph_integer_t byte = 1 + bitIndex / 6;
            if (byte >= length) break;
            igraph_integer_t bit = 5 - bitIndex % 6;
            igraph_integer_t value = str[byte] - 63;
            igraph_integer_t is_set = value & (1 << bit);
            if (is_set) {
                MATRIX(adj_matrix, i, j) = 1;
            }
            bitIndex++;
        }
    }
    igraph_adjacency(graph, &adj_matrix, IGRAPH_ADJ_UPPER, IGRAPH_NO_LOOPS);
    igraph_matrix_destroy(&adj_matrix);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_from_nauty(igraph_t *graph, const char *str) {
    if (str[0] == '&') {
        return igraph_from_nauty_digraph6(graph, str + 1);
    } else {
        return igraph_from_nauty_graph6(graph, str);
    }
}

static igraph_error_t igraph_to_nauty_graph6(const igraph_t *graph, char **res) {
    char *str_res = NULL;
    //TODO memory management
    str_res = IGRAPH_CALLOC(1000, char);

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    str_res[0] = no_of_nodes;

    igraph_integer_t bitIndex = 0;
    for (igraph_integer_t j = 1; j < no_of_nodes; j++) {
        for (igraph_integer_t i = 0; i < j; i++) {
            igraph_integer_t byte = 1 + bitIndex / 6;
            igraph_integer_t bit = 5 - bitIndex % 6;
            //TODO is this the best way to check if there is an edge between i and j?
            igraph_integer_t eid;
            igraph_get_eid(graph, &eid, i, j, false, false);
            if (eid != -1) {
                str_res[byte] |= (1 << bit);
            }
            bitIndex++;
        }
    }


    igraph_integer_t length = 1 + bitIndex / 6 + 1;
    for (igraph_integer_t i = 0; i < length; i++) {
        str_res[i] += 63;
    }

    *res = str_res;
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_to_nauty_digraph6(const igraph_t *graph, char **res) {
    char *str_res = NULL;
    //TODO memory management
    str_res = IGRAPH_CALLOC(1000, char);

    str_res[0] = '&';

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    str_res[1] = no_of_nodes;

    igraph_integer_t bitIndex = 0;
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
            igraph_integer_t byte = 2 + bitIndex / 6;
            igraph_integer_t bit = 5 - bitIndex % 6;
            //TODO is this the best way to check if there is an edge between i and j?
            igraph_integer_t eid;
            igraph_get_eid(graph, &eid, i, j, true, false);
            if (eid != -1) {
                str_res[byte] |= (1 << bit);
            }
            bitIndex++;
        }
    }


    igraph_integer_t length = 1 + bitIndex / 6 + 1;
    for (igraph_integer_t i = 1; i < length + 1; i++) {
        str_res[i] += 63;
    }

    *res = str_res;
    return IGRAPH_SUCCESS;
}


igraph_error_t igraph_to_nauty(const igraph_t *graph, char **res) {
    if (igraph_is_directed(graph)) {
        return igraph_to_nauty_digraph6(graph, res);
    } else {
        return igraph_to_nauty_graph6(graph, res);
    }
}
