/*
   IGraph library.
   Copyright (C) 2003-2023  The igraph development team <igraph@igraph.org>

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

#include "igraph_games.h"

#include "igraph_constructors.h"
#include "igraph_matrix.h"
#include "igraph_random.h"
#include "igraph_vector.h"

#include "core/interruption.h"
#include "math/safe_intop.h"

#include <float.h>      /* for DBL_EPSILON */
#include <math.h>       /* for sqrt and floor */

/**
 * \function igraph_sbm_game
 * \brief Sample from a stochastic block model.
 *
 * This function samples graphs from a stochastic block
 * model by (doing the equivalent of) Bernoulli
 * trials for each potential edge with the probabilities
 * given by the Bernoulli rate matrix, \p pref_matrix.
 * See Faust, K., &amp; Wasserman, S. (1992a). Blockmodels:
 * Interpretation and evaluation. Social Networks, 14, 5-–61.
 *
 * </para><para>
 * The order of the vertex IDs in the generated graph corresponds to
 * the \p block_sizes argument.
 *
 * \param graph The output graph. This should be a pointer to an
 *     uninitialized graph.
 * \param n Number of vertices.
 * \param pref_matrix The matrix giving the Bernoulli rates.
 *     This is a KxK matrix, where K is the number of groups.
 *     The probability of creating an edge between vertices from
 *     groups i and j is given by element (i,j).
 * \param block_sizes An integer vector giving the number of
 *     vertices in each group.
 * \param directed Boolean, whether to create a directed graph. If
 *     this argument is false, then \p pref_matrix must be symmetric.
 * \param loops Boolean, whether to create self-loops.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|+K^2), where |V| is the number of
 * vertices, |E| is the number of edges, and K is the number of
 * groups.
 *
 * \sa \ref igraph_erdos_renyi_game_gnp() for a simple Bernoulli graph.
 *
 */

igraph_error_t igraph_sbm_game(igraph_t *graph, igraph_integer_t n,
                    const igraph_matrix_t *pref_matrix,
                    const igraph_vector_int_t *block_sizes,
                    igraph_bool_t directed, igraph_bool_t loops) {

#define IGRAPH_CHECK_MAXEDGES() \
    do {if (maxedges > IGRAPH_MAX_EXACT_REAL) { \
        IGRAPH_ERROR("Too many vertices, overflow in maximum number of edges.", IGRAPH_EOVERFLOW); \
    }} while (0)

    igraph_integer_t no_blocks = igraph_matrix_nrow(pref_matrix);
    igraph_integer_t from, to, fromoff = 0;
    igraph_real_t minp, maxp;
    igraph_vector_int_t edges;

    /* ------------------------------------------------------------ */
    /* Check arguments                                              */
    /* ------------------------------------------------------------ */

    if (igraph_matrix_ncol(pref_matrix) != no_blocks) {
        IGRAPH_ERROR("Preference matrix is not square.",
                     IGRAPH_NONSQUARE);
    }

    if (no_blocks > 0) {
        igraph_matrix_minmax(pref_matrix, &minp, &maxp);
        if (minp < 0 || maxp > 1) {
            IGRAPH_ERROR("Connection probabilities must be in [0,1].", IGRAPH_EINVAL);
        }
    }

    if (!directed && !igraph_matrix_is_symmetric(pref_matrix)) {
        IGRAPH_ERROR("Preference matrix must be symmetric for undirected graphs.",
                     IGRAPH_EINVAL);
    }

    if (igraph_vector_int_size(block_sizes) != no_blocks) {
        IGRAPH_ERRORF("Block size vector length (%" IGRAPH_PRId ") does not agree with "
                      "preference matrix size (%" IGRAPH_PRId  ").", IGRAPH_EINVAL,
                      igraph_vector_int_size(block_sizes), no_blocks);
    }

    if (no_blocks > 0) {
        if (igraph_vector_int_min(block_sizes) < 0) {
            IGRAPH_ERRORF("Block sizes must be non-negative, but got %" IGRAPH_PRId ".",
                          IGRAPH_EINVAL, igraph_vector_int_min(block_sizes));
        }
    }

    if (igraph_vector_int_sum(block_sizes) != n) {
        IGRAPH_ERRORF("Sum of the block sizes (%" IGRAPH_PRId ") must equal the number of vertices (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, igraph_vector_int_sum(block_sizes), n);
    }

    /* Since the sum of the block sizes should equal the number of vertices,
     * and the block sizes are non-negative, the number of vertices is
     * guaranteed to be non-negative. This shouldn't be checked separately.
     */

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    RNG_BEGIN();

    for (from = 0; from < no_blocks; from++) {
        igraph_integer_t fromsize = VECTOR(*block_sizes)[from];
        igraph_integer_t start = directed ? 0 : from;
        igraph_integer_t i, tooff = 0;

        IGRAPH_ALLOW_INTERRUPTION();

        for (i = 0; i < start; i++) {
            tooff += VECTOR(*block_sizes)[i];
        }
        for (to = start; to < no_blocks; to++) {
            igraph_integer_t tosize = VECTOR(*block_sizes)[to];
            igraph_real_t prob = MATRIX(*pref_matrix, from, to);
            igraph_real_t maxedges;
            igraph_real_t last = RNG_GEOM(prob);  /* RNG_GEOM may return NaN so igraph_integer_t is not suitable */
            igraph_integer_t vfrom, vto;

            if (directed && loops) {
                maxedges = ((igraph_real_t) fromsize) * tosize;
                IGRAPH_CHECK_MAXEDGES();
                while (last < maxedges) {
                    vto = floor(last / fromsize);
                    vfrom = last - ((igraph_real_t) vto) * fromsize;
                    igraph_vector_int_push_back(&edges, fromoff + vfrom);
                    igraph_vector_int_push_back(&edges, tooff + vto);
                    last += RNG_GEOM(prob);
                    last += 1;
                }
            } else if (directed && !loops && from != to) {
                maxedges = ((igraph_real_t) fromsize) * tosize;
                IGRAPH_CHECK_MAXEDGES();
                while (last < maxedges) {
                    vto = floor(last / fromsize);
                    vfrom = last - ((igraph_real_t) vto) * fromsize;
                    igraph_vector_int_push_back(&edges, fromoff + vfrom);
                    igraph_vector_int_push_back(&edges, tooff + vto);
                    last += RNG_GEOM(prob);
                    last += 1;
                }
            } else if (directed && !loops && from == to) {
                maxedges = ((igraph_real_t) fromsize) * (fromsize - 1.0);
                IGRAPH_CHECK_MAXEDGES();
                while (last < maxedges) {
                    vto = floor(last / fromsize);
                    vfrom = last - ((igraph_real_t) vto) * fromsize;
                    if (vfrom == vto) {
                        vto = fromsize - 1;
                    }
                    igraph_vector_int_push_back(&edges, fromoff + vfrom);
                    igraph_vector_int_push_back(&edges, tooff + vto);
                    last += RNG_GEOM(prob);
                    last += 1;
                }
            } else if (!directed && loops && from != to) {
                maxedges = ((igraph_real_t) fromsize) * tosize;
                IGRAPH_CHECK_MAXEDGES();
                while (last < maxedges) {
                    vto = floor(last / fromsize);
                    vfrom = last - ((igraph_real_t) vto) * fromsize;
                    igraph_vector_int_push_back(&edges, fromoff + vfrom);
                    igraph_vector_int_push_back(&edges, tooff + vto);
                    last += RNG_GEOM(prob);
                    last += 1;
                }
            } else if (!directed && loops && from == to) {
                maxedges = ((igraph_real_t) fromsize) * (fromsize + 1.0) / 2.0;
                IGRAPH_CHECK_MAXEDGES();
                while (last < maxedges) {
                    vto = floor((sqrt(8 * last + 1) - 1) / 2);
                    vfrom = last - (((igraph_real_t) vto) * (vto + 1.0)) / 2.0;
                    igraph_vector_int_push_back(&edges, fromoff + vfrom);
                    igraph_vector_int_push_back(&edges, tooff + vto);
                    last += RNG_GEOM(prob);
                    last += 1;
                }
            } else if (!directed && !loops && from != to) {
                maxedges = ((igraph_real_t) fromsize) * tosize;
                IGRAPH_CHECK_MAXEDGES();
                while (last < maxedges) {
                    vto = floor(last / fromsize);
                    vfrom = last - ((igraph_real_t) vto) * fromsize;
                    igraph_vector_int_push_back(&edges, fromoff + vfrom);
                    igraph_vector_int_push_back(&edges, tooff + vto);
                    last += RNG_GEOM(prob);
                    last += 1;
                }
            } else { /*!directed && !loops && from==to */
                maxedges = ((igraph_real_t) fromsize) * (fromsize - 1.0) / 2.0;
                IGRAPH_CHECK_MAXEDGES();
                while (last < maxedges) {
                    vto = floor((sqrt(8 * last + 1) + 1) / 2);
                    vfrom = last - (((igraph_real_t) vto) * (vto - 1.0)) / 2.0;
                    igraph_vector_int_push_back(&edges, fromoff + vfrom);
                    igraph_vector_int_push_back(&edges, tooff + vto);
                    last += RNG_GEOM(prob);
                    last += 1;
                }
            }

            tooff += tosize;
        }
        fromoff += fromsize;
    }

    RNG_END();

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
#undef IGRAPH_CHECK_MAXEDGES
}

/**
 * \function igraph_hsbm_game
 * \brief Hierarchical stochastic block model.
 *
 * The function generates a random graph according to the hierarchical
 * stochastic block model.
 *
 * \param graph The generated graph is stored here.
 * \param n The number of vertices in the graph.
 * \param m The number of vertices per block. n/m must be integer.
 * \param rho The fraction of vertices per cluster,
 *        within a block. Must sum up to 1, and rho * m must be integer
 *        for all elements of rho.
 * \param C A square, symmetric numeric matrix, the Bernoulli rates for
 *        the clusters within a block. Its size must mach the size of the
 *        \p rho vector.
 * \param p The Bernoulli rate of connections between
 *        vertices in different blocks.
 * \return Error code.
 *
 * \sa \ref igraph_sbm_game() for the classic stochastic block model,
 * \ref igraph_hsbm_list_game() for a more general version.
 */

igraph_error_t igraph_hsbm_game(igraph_t *graph, igraph_integer_t n,
                     igraph_integer_t m, const igraph_vector_t *rho,
                     const igraph_matrix_t *C, igraph_real_t p) {

#define IGRAPH_CHECK_MAXEDGES() \
    do {if (maxedges > IGRAPH_MAX_EXACT_REAL) { \
        IGRAPH_ERROR("Too many vertices, overflow in maximum number of edges.", IGRAPH_EOVERFLOW); \
    }} while (0)
    igraph_integer_t b, i, k = igraph_vector_size(rho);
    igraph_vector_t csizes;
    igraph_real_t sq_dbl_epsilon = sqrt(DBL_EPSILON);
    igraph_integer_t no_blocks = n / m;
    igraph_vector_int_t edges;
    igraph_integer_t offset = 0;

    if (n < 1) {
        IGRAPH_ERROR("`n' must be positive for HSBM", IGRAPH_EINVAL);
    }
    if (m < 1) {
        IGRAPH_ERROR("`m' must be positive for HSBM", IGRAPH_EINVAL);
    }
    if (n % m) {
        IGRAPH_ERROR("`n' must be a multiple of `m' for HSBM", IGRAPH_EINVAL);
    }
    if (!igraph_vector_isininterval(rho, 0, 1)) {
        IGRAPH_ERROR("`rho' must be between zero and one for HSBM",
                     IGRAPH_EINVAL);
    }
    if (igraph_matrix_min(C) < 0 || igraph_matrix_max(C) > 1) {
        IGRAPH_ERROR("`C' must be between zero and one for HSBM", IGRAPH_EINVAL);
    }
    if (fabs(igraph_vector_sum(rho) - 1.0) > sq_dbl_epsilon) {
        IGRAPH_ERROR("`rho' must sum up to 1 for HSBM", IGRAPH_EINVAL);
    }
    if (igraph_matrix_nrow(C) != k || igraph_matrix_ncol(C) != k) {
        IGRAPH_ERROR("`C' dimensions must match `rho' dimensions in HSBM",
                     IGRAPH_EINVAL);
    }
    if (!igraph_matrix_is_symmetric(C)) {
        IGRAPH_ERROR("`C' must be a symmetric matrix", IGRAPH_EINVAL);
    }
    if (p < 0 || p > 1) {
        IGRAPH_ERROR("`p' must be a probability for HSBM", IGRAPH_EINVAL);
    }
    for (i = 0; i < k; i++) {
        igraph_real_t s = VECTOR(*rho)[i] * m;
        if (fabs(round(s) - s) > sq_dbl_epsilon) {
            IGRAPH_ERROR("`rho' * `m' is not integer in HSBM", IGRAPH_EINVAL);
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&csizes, k);
    for (i = 0; i < k; i++) {
        VECTOR(csizes)[i] = round(VECTOR(*rho)[i] * m);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    RNG_BEGIN();

    /* Block models first */

    for (b = 0; b < no_blocks; b++) {
        igraph_integer_t from, to, fromoff = 0;

        for (from = 0; from < k; from++) {
            igraph_integer_t fromsize = VECTOR(csizes)[from];
            igraph_integer_t i, tooff = 0;
            for (i = 0; i < from; i++) {
                tooff += VECTOR(csizes)[i];
            }
            for (to = from; to < k; to++) {
                igraph_integer_t tosize = VECTOR(csizes)[to];
                igraph_real_t prob = MATRIX(*C, from, to);
                igraph_real_t maxedges;
                igraph_real_t last = RNG_GEOM(prob);  /* RNG_GEOM may return NaN so igraph_integer_t is not suitable */
                if (from != to) {
                    maxedges = ((igraph_real_t) fromsize) * tosize;
                    IGRAPH_CHECK_MAXEDGES();
                    while (last < maxedges) {
                        igraph_integer_t vto = floor(last / fromsize);
                        igraph_integer_t vfrom = last - ((igraph_real_t) vto) * fromsize;
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, offset + fromoff + vfrom));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, offset + tooff + vto));
                        last += RNG_GEOM(prob);
                        last += 1;
                    }
                } else { /* from==to */
                    maxedges = ((igraph_real_t) fromsize) * (fromsize - 1.0) / 2.0;
                    IGRAPH_CHECK_MAXEDGES();
                    while (last < maxedges) {
                        igraph_integer_t vto = floor((sqrt(8 * last + 1) + 1) / 2);
                        igraph_integer_t vfrom = last - (((igraph_real_t) vto) * (vto - 1.0)) / 2.0;
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, offset + fromoff + vfrom));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, offset + tooff + vto));
                        last += RNG_GEOM(prob);
                        last += 1;
                    }
                }

                tooff += tosize;
            }
            fromoff += fromsize;
        }

        offset += m;
    }

    /* And now the rest, if not a special case */

    if (p == 1) {
        igraph_integer_t fromoff = 0, tooff = m;
        for (b = 0; b < no_blocks; b++) {
            igraph_integer_t fromsize = m;
            igraph_integer_t tosize = n - tooff;
            igraph_integer_t from, to;
            for (from = 0; from < fromsize; from++) {
                for (to = 0; to < tosize; to++) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, fromoff + from));
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, tooff + to));
                }
            }
            fromoff += m;
            tooff += m;
        }
    } else if (p > 0) {
        igraph_integer_t fromoff = 0, tooff = m;
        for (b = 0; b < no_blocks; b++) {
            igraph_integer_t fromsize = m;
            igraph_integer_t tosize = n - tooff;
            igraph_real_t maxedges = ((igraph_real_t) fromsize) * tosize;
            IGRAPH_CHECK_MAXEDGES();
            igraph_real_t last = RNG_GEOM(p);  /* RNG_GEOM may return NaN so igraph_integer_t is not suitable */
            while (last < maxedges) {
                igraph_integer_t vto = floor(last / fromsize);
                igraph_integer_t vfrom = last - ((igraph_real_t) vto) * fromsize;
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, fromoff + vfrom));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, tooff + vto));
                last += RNG_GEOM(p);
                last += 1;
            }

            fromoff += m;
            tooff += m;
        }
    }

    RNG_END();

    IGRAPH_CHECK(igraph_create(graph, &edges, n, /*directed=*/ 0));

    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&csizes);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
#undef IGRAPH_CHECK_MAXEDGES
}

/**
 * \function igraph_hsbm_list_game
 * \brief Hierarchical stochastic block model, more general version.
 *
 * The function generates a random graph according to the hierarchical
 * stochastic block model.
 *
 * \param graph The generated graph is stored here.
 * \param n The number of vertices in the graph.
 * \param mlist An integer vector of block sizes.
 * \param rholist A list of rho vectors (\c igraph_vector_t objects), one
 *        for each block.
 * \param Clist A list of square matrices (\c igraph_matrix_t objects),
 *        one for each block, specifying the Bernoulli rates of connections
 *        within the block.
 * \param p The Bernoulli rate of connections between
 *        vertices in different blocks.
 * \return Error code.
 *
 * \sa \ref igraph_sbm_game() for the classic stochastic block model,
 * \ref igraph_hsbm_game() for a simpler general version.
 */

igraph_error_t igraph_hsbm_list_game(igraph_t *graph, igraph_integer_t n,
                          const igraph_vector_int_t *mlist,
                          const igraph_vector_list_t *rholist,
                          const igraph_matrix_list_t *Clist,
                          igraph_real_t p) {

    igraph_integer_t i, no_blocks = igraph_vector_list_size(rholist);
    igraph_real_t sq_dbl_epsilon = sqrt(DBL_EPSILON);
    igraph_vector_int_t edges;
    igraph_vector_t csizes;
    igraph_integer_t b, offset = 0;

    if (n < 1) {
        IGRAPH_ERROR("`n' must be positive for HSBM.", IGRAPH_EINVAL);
    }
    if (no_blocks == 0) {
        IGRAPH_ERROR("`rholist' empty for HSBM.", IGRAPH_EINVAL);
    }
    if (igraph_matrix_list_size(Clist) != no_blocks &&
        igraph_vector_int_size(mlist) != no_blocks) {
        IGRAPH_ERROR("`rholist' must have same length as `Clist' and `m' "
                     "for HSBM.", IGRAPH_EINVAL);
    }
    if (p < 0 || p > 1) {
        IGRAPH_ERROR("`p' must be a probability for HSBM.", IGRAPH_EINVAL);
    }
    /* Checks for m's */
    if (igraph_vector_int_sum(mlist) != n) {
        IGRAPH_ERROR("`m' must sum up to `n' for HSBM.", IGRAPH_EINVAL);
    }
    if (igraph_vector_int_min(mlist) < 1) {
        IGRAPH_ERROR("`m' must be positive for HSBM.", IGRAPH_EINVAL);
    }
    /* Checks for the rhos */
    for (i = 0; i < no_blocks; i++) {
        const igraph_vector_t *rho = igraph_vector_list_get_ptr(rholist, i);
        if (!igraph_vector_isininterval(rho, 0, 1)) {
            IGRAPH_ERROR("`rho' must be between zero and one for HSBM.",
                         IGRAPH_EINVAL);
        }
        if (fabs(igraph_vector_sum(rho) - 1.0) > sq_dbl_epsilon) {
            IGRAPH_ERROR("`rho' must sum up to 1 for HSBM.", IGRAPH_EINVAL);
        }
    }
    /* Checks for the Cs */
    for (i = 0; i < no_blocks; i++) {
        const igraph_matrix_t *C = igraph_matrix_list_get_ptr(Clist, i);
        if (igraph_matrix_min(C) < 0 || igraph_matrix_max(C) > 1) {
            IGRAPH_ERROR("Bernoulli rates must be between zero and one for HSBM.",
                         IGRAPH_EINVAL);
        }
        if (!igraph_matrix_is_symmetric(C)) {
            IGRAPH_ERROR("Bernoulli rate matrices must be symmetric.", IGRAPH_EINVAL);
        }
    }
    /* Check that C and rho sizes match */
    for (i = 0; i < no_blocks; i++) {
        const igraph_vector_t *rho = igraph_vector_list_get_ptr(rholist, i);
        const igraph_matrix_t *C = igraph_matrix_list_get_ptr(Clist, i);
        igraph_integer_t k = igraph_vector_size(rho);
        if (igraph_matrix_nrow(C) != k || igraph_matrix_ncol(C) != k) {
            IGRAPH_ERROR("All Bernoulli rate matrix dimensions must match `rho' "
                    "dimensions in HSBM.",
                         IGRAPH_EINVAL);
        }
    }
    /* Check that rho * m is integer */
    for (i = 0; i < no_blocks; i++) {
        const igraph_vector_t *rho = igraph_vector_list_get_ptr(rholist, i);
        igraph_real_t m = VECTOR(*mlist)[i];
        igraph_integer_t j, k = igraph_vector_size(rho);
        for (j = 0; j < k; j++) {
            igraph_real_t s = VECTOR(*rho)[j] * m;
            if (fabs(round(s) - s) > sq_dbl_epsilon) {
                IGRAPH_ERROR("`rho' * `m' is not integer in HSBM.", IGRAPH_EINVAL);
            }
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&csizes, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    RNG_BEGIN();

    /* Block models first */

    for (b = 0; b < no_blocks; b++) {
        igraph_integer_t from, to, fromoff = 0;
        const igraph_vector_t *rho = igraph_vector_list_get_ptr(rholist, b);
        const igraph_matrix_t *C = igraph_matrix_list_get_ptr(Clist, b);
        igraph_integer_t m = VECTOR(*mlist)[b];
        igraph_integer_t k = igraph_vector_size(rho);

        IGRAPH_CHECK(igraph_vector_resize(&csizes, k));
        for (i = 0; i < k; i++) {
            VECTOR(csizes)[i] = round(VECTOR(*rho)[i] * m);
        }

        for (from = 0; from < k; from++) {
            igraph_integer_t fromsize = VECTOR(csizes)[from];
            igraph_integer_t i, tooff = 0;
            for (i = 0; i < from; i++) {
                tooff += VECTOR(csizes)[i];
            }
            for (to = from; to < k; to++) {
                igraph_integer_t tosize = VECTOR(csizes)[to];
                igraph_real_t prob = MATRIX(*C, from, to);
                igraph_real_t maxedges;
                igraph_real_t last = RNG_GEOM(prob);  /* RNG_GEOM may return NaN so igraph_integer_t is not suitable */
                if (from != to) {
                    maxedges = ((igraph_real_t) fromsize) * tosize;
                    while (last < maxedges) {
                        igraph_integer_t vto = floor(last / fromsize);
                        igraph_integer_t vfrom = last - ((igraph_real_t) vto) * fromsize;
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, offset + fromoff + vfrom));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, offset + tooff + vto));
                        last += RNG_GEOM(prob);
                        last += 1;
                    }
                } else { /* from==to */
                    maxedges = ((igraph_real_t) fromsize) * (fromsize - 1.0) / 2.0;
                    while (last < maxedges) {
                        igraph_integer_t vto = floor((sqrt(8 * last + 1) + 1) / 2);
                        igraph_integer_t vfrom = last - (((igraph_real_t) vto) * (vto - 1.0)) / 2.0;
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, offset + fromoff + vfrom));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, offset + tooff + vto));
                        last += RNG_GEOM(prob);
                        last += 1;
                    }
                }

                tooff += tosize;
            }
            fromoff += fromsize;
        }

        offset += m;
    }

    /* And now the rest, if not a special case */

    if (p == 1) {
        igraph_integer_t fromoff = 0, tooff = VECTOR(*mlist)[0];
        for (b = 0; b < no_blocks; b++) {
            igraph_integer_t fromsize = VECTOR(*mlist)[b];
            igraph_integer_t tosize = n - tooff;
            igraph_integer_t from, to;
            for (from = 0; from < fromsize; from++) {
                for (to = 0; to < tosize; to++) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, fromoff + from));
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, tooff + to));
                }
            }
            fromoff += fromsize;
            if (b + 1 < no_blocks) {
                tooff += VECTOR(*mlist)[b + 1];
            }
        }
    } else if (p > 0) {
        igraph_integer_t fromoff = 0, tooff = VECTOR(*mlist)[0];
        for (b = 0; b < no_blocks; b++) {
            igraph_integer_t fromsize = VECTOR(*mlist)[b];
            igraph_integer_t tosize = n - tooff;
            igraph_real_t maxedges = ((igraph_real_t) fromsize) * tosize;
            igraph_real_t last = RNG_GEOM(p);  /* RNG_GEOM may return NaN so igraph_integer_t is not suitable */
            while (last < maxedges) {
                igraph_integer_t vto = floor(last / fromsize);
                igraph_integer_t vfrom = last - ((igraph_real_t) vto) * fromsize;
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, fromoff + vfrom));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, tooff + vto));
                last += RNG_GEOM(p);
                last += 1;
            }

            fromoff += fromsize;
            if (b + 1 < no_blocks) {
                tooff += VECTOR(*mlist)[b + 1];
            }
        }
    }

    RNG_END();

    IGRAPH_CHECK(igraph_create(graph, &edges, n, /*directed=*/ 0));

    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&csizes);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
