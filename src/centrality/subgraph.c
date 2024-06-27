/*
   IGraph library.
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

#include "igraph_centrality.h"

#include "igraph_blas.h"
#include "igraph_conversion.h"
#include "igraph_interface.h"
#include "igraph_lapack.h"

#include <math.h>

/**
 * \function igraph_subgraph_centrality
 * \brief The subgraph centrality of vertices.
 *
 * \experimental
 *
 * The subgraph centrality of vertex \c i is formally defined as
 * <code>(exp(A))_ii</code>,
 * where \c A is the adjacency matrix and \c exp denotes the matrix
 * exponential. This value reflects the number of closed walks that each
 * vertex participates in, downweighted by the factorial of the walk length.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * E Estrada, J A Rodríguez-Velázquez:
 * Subgraph Centrality in Complex Networks.
 * Physical Review E 71, no. 5 (2005): 056103.
 * https://doi.org/10.1103/PhysRevE.71.056103.
 *
 * \param graph The input graph. Edge directions are ignored.
 * \param res An initialized vector, the result will be stored here. It will be
 *   resized as needed.
 * \param loops If false, self-loops are ignored. Otherwise, the diagonal of
 *   the adjacency matrix is assumed to contain twice the number of self-loops.
 * \return Error code.
 *
 * Time complexity: Dominated by the DSYEVR LAPACK routine.
 *
 * </para><para>
 * Memory complexity: O(|V|^2) where |V| is the number of vertices.
 */
igraph_error_t igraph_subgraph_centrality(
        const igraph_t *graph, igraph_vector_t *res, igraph_bool_t loops) {

    const igraph_integer_t vcount = igraph_vcount(graph);
    const igraph_real_t eps = 1e-10;
    igraph_vector_t values;
    igraph_matrix_t adjmat, vectors;

    /* This function computes the matrix exponential by diagonalizing the matrix
     * first. This technique is only appropriate for symmetric matrices, therefore
     * we ignore edge directions. */

    IGRAPH_VECTOR_INIT_FINALLY(&values, vcount);
    IGRAPH_MATRIX_INIT_FINALLY(&vectors, vcount, vcount);
    IGRAPH_MATRIX_INIT_FINALLY(&adjmat, vcount, vcount);

    IGRAPH_CHECK(igraph_get_adjacency(
        graph, &adjmat, IGRAPH_GET_ADJACENCY_BOTH,
        /* weights */ NULL, loops ? IGRAPH_LOOPS_TWICE : IGRAPH_NO_LOOPS));

    /* Symmetrize if directed */
    if (igraph_is_directed(graph)) {
        for (igraph_integer_t i=0; i < vcount; i++) {
            for (igraph_integer_t j=0; j < vcount; j++) {
                MATRIX(adjmat, i, j) += MATRIX(adjmat, j, i);
            }
        }
    }

    IGRAPH_CHECK(igraph_lapack_dsyevr(
        &adjmat, IGRAPH_LAPACK_DSYEV_ALL,
        0, 0, 0, 0, 0, /* these parameters are unused with IGRAPH_LAPACK_DSYEV_ALL */
        eps, &values, &vectors, NULL));

    igraph_matrix_destroy(&adjmat);
    IGRAPH_FINALLY_CLEAN(1);

    for (igraph_integer_t j=0; j < vcount; j++) {
        for (igraph_integer_t i=0; i < vcount; i++) {
            MATRIX(vectors, i, j) *= MATRIX(vectors, i, j);
        }
    }

    for (igraph_integer_t i=0; i < vcount; i++) {
        VECTOR(values)[i] = exp(VECTOR(values)[i]);
    }

    IGRAPH_CHECK(igraph_vector_resize(res, vcount));
    IGRAPH_CHECK(igraph_blas_dgemv(false, /* alpha */ 1, &vectors, &values, /* beta */ 0, res));

    igraph_matrix_destroy(&vectors);
    igraph_vector_destroy(&values);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
