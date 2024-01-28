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

static igraph_integer_t igraph_read_graph6_number_and_skip(char *str) {
    //TODO bigger numbers
    igraph_integer_t res = str[0] - 63;
    str++;
    return res;
}

igraph_error_t igraph_read_graph6(igraph_t *graph, FILE *infile) {
    // pipes and sockets don't support seeking.
    // https://stackoverflow.com/a/44894946 to read the whole file, or use some standard igraph solution?
    igraph_integer_t fsize = ftell(f);
    fseek(infile, 0, SEEK_END);
    fseek(infile, 0, SEEK_SET);

    char *instring = IGRAPH_CALLOC(fsize + 1);
    fread(string, fsize, 1, f);
    fclose(f);
    igraph_integer_t length = strlen(instring);

    igraph_integer_t nr_nodes = igraph_read_graph6_number_and_skip(str)
    igraph_matrix_t adj_matrix;

    IGRAPH_CHECK(igraph_matrix_init(adj_matrix, nr_nodes, nr_nodes));
    IGRAPH_FINALLY(igraph_matrix_destroy, adj_matrix);

    int bitIndex = 0;
    for (int j = 1; j < nr_nodes; j++) {
        for (int i = 0; i < j; i++) {
            //TODO bigger numbers
            int byte = 1 + bitIndex / 6;
            if (byte >= length) break;
            int bit = 5 - bitIndex % 6;
            int value = str[byte] - 63;
            int is_set = value & (1 << bit);
            if (is_set) {
                MATRIX(adj_matrix, i, j) = 1;
            }
            bitIndex++;
        }
    }
    igraph_adjacency(graph, adj_matrix, IGRAPH_ADJ_UPPER, IGRAPH_NO_LOOPS);
    igraph_matrix_destroy(&adj_matrix);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_write_graph6(const igraph_t *graph, FILE *outfile) {

    char str_res[100] = {0};

    bitIndex = 0;
    for (int j = 1; j < nr_nodes; j++) {
        for (int i = 0; i < j; i++) {
            int byte = 1 + bitIndex / 6;
            int bit = 5 - bitIndex % 6;
            if (adj_matrix[i][j]) {
                str_res[byte] |= (1 << bit);
            }
            bitIndex++;
        }
    }


    str_res[0] = nr_nodes;
    for (int i = 0; i < length; i++) {
        str_res[i] += 63;
    }
    printf("%s", str_res);

    return 0;
}
