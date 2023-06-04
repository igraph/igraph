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

#include "internal/utils.h"

/**
 * \function igraph_i_matrix_subset_vertices
 * \brief Subsets a matrix whose rows/columns correspond to graph vertices.
 *
 * This is a convenience function to subset a matrix computed from a graph.
 * It takes a matrix whose rows and columns correspond to the vertices
 * of a graph, and subsets it in-place to retain only some of the vertices.
 *
 * \param m A square matrix with the same number of rows/columns as the vertex
 *    count of \p graph. It will be modified in-place, deleting rows \em not present
 *    in \p from and columns \em not present in \p to.
 * \param graph The corresponding graph. <code>m[u,v]</code> is assumed to contain
 *    a value associated with vertices \c u and \c v of \p graph, e.g. the graph
 *    distance between them, their similarity, etc.
 * \param from Vertex set, these rows of the matrix will be retained.
 * \param to Vertex set, these columns of the matrix will be retained.
 * \return Error code.
 *
 * Time complexity:
 * O(1) when taking all vertices,
 * O(|from|*|to|) otherwise where |from| and |to| denote the size
 * of the source and target vertex sets.
 */
igraph_error_t igraph_i_matrix_subset_vertices(
        igraph_matrix_t *m,
        const igraph_t *graph,
        igraph_vs_t from,
        igraph_vs_t to) {

    /* Assertion: the size of 'm' agrees with 'graph': */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t ncol = igraph_matrix_ncol(m);
    igraph_integer_t nrow = igraph_matrix_nrow(m);

    IGRAPH_ASSERT(nrow == no_of_nodes && nrow == ncol);

    /* When taking all vertices, nothing needs to be done: */

    if (igraph_vs_is_all(&from) && igraph_vs_is_all(&to)) {
        return IGRAPH_SUCCESS;
    }

    /* Otherwise, allocate a temporary matrix to copy the data into: */

    igraph_vit_t fromvit, tovit;
    igraph_matrix_t tmp;

    IGRAPH_CHECK(igraph_vit_create(graph, from, &fromvit));
    IGRAPH_FINALLY(igraph_vit_destroy, &fromvit);

    IGRAPH_CHECK(igraph_vit_create(graph, to, &tovit));
    IGRAPH_FINALLY(igraph_vit_destroy, &tovit);

    IGRAPH_MATRIX_INIT_FINALLY(&tmp, IGRAPH_VIT_SIZE(fromvit), IGRAPH_VIT_SIZE(tovit));

    for (igraph_integer_t j=0; ! IGRAPH_VIT_END(tovit); IGRAPH_VIT_NEXT(tovit), j++) {
        igraph_integer_t i;
        for (IGRAPH_VIT_RESET(fromvit), i=0; ! IGRAPH_VIT_END(fromvit); IGRAPH_VIT_NEXT(fromvit), i++) {
            MATRIX(tmp, i, j) = MATRIX(*m, IGRAPH_VIT_GET(fromvit), IGRAPH_VIT_GET(tovit));
        }
    }

    /* This is O(1) time */
    IGRAPH_CHECK(igraph_matrix_swap(m, &tmp));

    igraph_matrix_destroy(&tmp);
    igraph_vit_destroy(&tovit);
    igraph_vit_destroy(&fromvit);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}
