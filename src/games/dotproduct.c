/*
   igraph library.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_games.h"

#include "igraph_blas.h"
#include "igraph_constructors.h"
#include "igraph_random.h"

/**
 * \function igraph_dot_product_game
 * \brief Generates a random dot product graph.
 *
 * In this model, each vertex is represented by a latent
 * position vector. Probability of an edge between two vertices are given
 * by the dot product of their latent position vectors.
 *
 * </para><para>
 * See also Christine Leigh Myers Nickel: Random dot product graphs, a
 * model for social networks. Dissertation, Johns Hopkins University,
 * Maryland, USA, 2006.
 *
 * \param graph The output graph is stored here.
 * \param vecs A matrix in which each latent position vector is a
 *    column. The dot product of the latent position vectors should be
 *    in the [0,1] interval, otherwise a warning is given. For
 *    negative dot products, no edges are added; dot products that are
 *    larger than one always add an edge.
 * \param directed Should the generated graph be directed?
 * \return Error code.
 *
 * Time complexity: O(n*n*m), where n is the number of vertices,
 * and m is the length of the latent vectors.
 *
 * \sa \ref igraph_rng_sample_dirichlet(), \ref
 * igraph_rng_sample_sphere_volume(), \ref igraph_rng_sample_sphere_surface()
 * for functions to generate the latent vectors.
 */

igraph_error_t igraph_dot_product_game(igraph_t *graph, const igraph_matrix_t *vecs,
                            igraph_bool_t directed) {

    igraph_int_t nrow = igraph_matrix_nrow(vecs);
    igraph_int_t ncol = igraph_matrix_ncol(vecs);
    igraph_int_t i, j;
    igraph_vector_int_t edges;
    igraph_bool_t warned_neg = false, warned_big = false;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    for (i = 0; i < ncol; i++) {
        igraph_int_t from = directed ? 0 : i + 1;
        const igraph_vector_t v1 = igraph_vector_view(&MATRIX(*vecs, 0, i), nrow);
        for (j = from; j < ncol; j++) {
            igraph_real_t prob;
            const igraph_vector_t v2 = igraph_vector_view(&MATRIX(*vecs, 0, j), nrow);

            if (i == j) {
                continue;
            }
            IGRAPH_CHECK(igraph_blas_ddot(&v1, &v2, &prob));
            if (prob < 0 && ! warned_neg) {
                warned_neg = true;
                IGRAPH_WARNING("Negative connection probability in dot-product graph.");
            } else if (prob > 1 && ! warned_big) {
                warned_big = true;
                IGRAPH_WARNING("Greater than 1 connection probability in dot-product graph.");
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j));
            } else if (RNG_UNIF01() < prob) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j));
            }
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, ncol, directed));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
