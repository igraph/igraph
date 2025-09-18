/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

#include <igraph.h>

/* Create the condensation of a directed graph.
 * See https://en.wikipedia.org/wiki/Strongly_connected_component#Definitions
 * This example demonstrates how to write a basic igraph function, complete
 * with error handling. */
igraph_error_t condensation(const igraph_t *graph, igraph_t *cond) {
    igraph_vector_int_t membership;

    /* Data structures such as vector must be initialized in igraph before use. */
    IGRAPH_CHECK(igraph_vector_int_init(&membership, 0));

    /* Adding the initialized vector to the "finally" stack ensures that it will
     * be automatically destroyed if an error occurs. */
    IGRAPH_FINALLY(igraph_vector_int_destroy, &membership);

    /* Functions that return an error code can be wrapped in IGRAPH_CHECK to pass that error
     * up to the caller. */
    IGRAPH_CHECK(igraph_connected_components(graph, &membership, /* csize */ NULL, /* no */ NULL, IGRAPH_STRONG));

    /* To compute the condensation, we simply contract strongly connected components.
     * Since igraph_contract_vertices() modifies graphs in-place, we make a copy first. */
    IGRAPH_CHECK(igraph_copy(cond, graph));

    /* Since we are not done creating the condensation yet, we add 'cond' to the
     * "finally" stack, so that it will be destroyed if an error occurs. */
    IGRAPH_FINALLY(igraph_destroy, cond);

    /* Contract strongly connected components. */
    IGRAPH_CHECK(igraph_contract_vertices(cond, &membership, NULL));

    /* igraph_contract_vertices() preserves all edges, some of which become
     * parallel edges or self-loops after the contraction. We simplify these. */
    IGRAPH_CHECK(igraph_simplify(cond, /* remove_multiple */ true, /* remove_loops */ true, NULL));

    /* Data structures that are no longer needed must be explicitly destroyed.
     * If they were added to the "finally" stack, they must be removed explicitly,
     * in the opposite order to how they were added. IGRAPH_FINALLY_CLEAN removes
     * the indicated number of entries from the "finally" stack. We remove
     * 'membership' because it was destroyed, and 'cond' because the responsibility
     * to destroy it is now with the caller. */
    igraph_vector_int_destroy(&membership);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS; /* return with no error */
}

int main(void) {
    igraph_t graph, cond;

    /* Initialize the library. */
    igraph_setup();

    /* Create a random directed graph with mean degree 2 and compute its condensation. */
    igraph_erdos_renyi_game_gnm(&graph, 100, 200, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    condensation(&graph, &cond);

    printf("Number of vertices in the condensation: %" IGRAPH_PRId "\n", igraph_vcount(&cond));
    igraph_write_graph_edgelist(&cond, stdout);

    /* Destroy data structures that are no longer needed. */
    igraph_destroy(&graph);
    igraph_destroy(&cond);

    return 0;
}
