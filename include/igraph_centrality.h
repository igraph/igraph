/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_CENTRALITY_H
#define IGRAPH_CENTRALITY_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_types.h"
#include "igraph_datatype.h"
#include "igraph_iterators.h"
#include "igraph_arpack.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Centrality                                         */
/* -------------------------------------------------- */

DECLDIR int igraph_closeness(const igraph_t *graph, igraph_vector_t *res,
                             const igraph_vs_t vids, igraph_neimode_t mode,
                             const igraph_vector_t *weights, igraph_bool_t normalized);
DECLDIR int igraph_closeness_estimate(const igraph_t *graph, igraph_vector_t *res,
                                      const igraph_vs_t vids, igraph_neimode_t mode,
                                      igraph_real_t cutoff,
                                      const igraph_vector_t *weights,
                                      igraph_bool_t normalized);

DECLDIR int igraph_betweenness(const igraph_t *graph, igraph_vector_t *res,
                               const igraph_vs_t vids, igraph_bool_t directed,
                               const igraph_vector_t *weights, igraph_bool_t nobigint);
DECLDIR int igraph_betweenness_estimate(const igraph_t *graph, igraph_vector_t *res,
                                        const igraph_vs_t vids, igraph_bool_t directed,
                                        igraph_real_t cutoff,
                                        const igraph_vector_t *weights,
                                        igraph_bool_t nobigint);
DECLDIR int igraph_edge_betweenness(const igraph_t *graph, igraph_vector_t *result,
                                    igraph_bool_t directed,
                                    const igraph_vector_t *weigths);
DECLDIR int igraph_edge_betweenness_estimate(const igraph_t *graph, igraph_vector_t *result,
        igraph_bool_t directed, igraph_real_t cutoff,
        const igraph_vector_t *weights);
DECLDIR int igraph_pagerank_old(const igraph_t *graph, igraph_vector_t *res,
                                const igraph_vs_t vids, igraph_bool_t directed,
                                igraph_integer_t niter, igraph_real_t eps,
                                igraph_real_t damping, igraph_bool_t old);

/**
 * \typedef igraph_pagerank_algo_t
 * \brief PageRank algorithm implementation
 *
 * Algorithms to calculate PageRank.
 * \enumval IGRAPH_PAGERANK_ALGO_POWER Use a simple power iteration,
 *   as it was implemented before igraph version 0.5.
 * \enumval IGRAPH_PAGERANK_ALGO_ARPACK Use the ARPACK library, this
 *   was the PageRank implementation in igraph from version 0.5, until
 *   version 0.7.
 * \enumval IGRAPH_PAGERANK_ALGO_PRPACK Use the PRPACK
 *   library. Currently this implementation is recommended.
 */

typedef enum {
    IGRAPH_PAGERANK_ALGO_POWER = 0,
    IGRAPH_PAGERANK_ALGO_ARPACK = 1,
    IGRAPH_PAGERANK_ALGO_PRPACK = 2
} igraph_pagerank_algo_t;

/**
 * \struct igraph_pagerank_power_options_t
 * \brief Options for the power method
 *
 * \member niter The number of iterations to perform, integer.
 * \member eps  The algorithm will consider the calculation as complete
 *        if the difference of values between iterations change
 *        less than this value for every vertex.
 */

typedef struct igraph_pagerank_power_options_t {
    igraph_integer_t niter;
    igraph_real_t eps;
} igraph_pagerank_power_options_t;

DECLDIR int igraph_pagerank(const igraph_t *graph, igraph_pagerank_algo_t algo,
                            igraph_vector_t *vector,
                            igraph_real_t *value, const igraph_vs_t vids,
                            igraph_bool_t directed, igraph_real_t damping,
                            const igraph_vector_t *weights, void *options);
DECLDIR int igraph_personalized_pagerank(const igraph_t *graph,
        igraph_pagerank_algo_t algo, igraph_vector_t *vector,
        igraph_real_t *value, const igraph_vs_t vids,
        igraph_bool_t directed, igraph_real_t damping,
        igraph_vector_t *reset,
        const igraph_vector_t *weights, void *options);
DECLDIR int igraph_personalized_pagerank_vs(const igraph_t *graph,
        igraph_pagerank_algo_t algo,
        igraph_vector_t *vector,
        igraph_real_t *value, const igraph_vs_t vids,
        igraph_bool_t directed, igraph_real_t damping,
        igraph_vs_t reset_vids,
        const igraph_vector_t *weights, void *options);

DECLDIR int igraph_eigenvector_centrality(const igraph_t *graph, igraph_vector_t *vector,
        igraph_real_t *value,
        igraph_bool_t directed, igraph_bool_t scale,
        const igraph_vector_t *weights,
        igraph_arpack_options_t *options);

DECLDIR int igraph_hub_score(const igraph_t *graph, igraph_vector_t *vector,
                             igraph_real_t *value, igraph_bool_t scale,
                             const igraph_vector_t *weights,
                             igraph_arpack_options_t *options);
DECLDIR int igraph_authority_score(const igraph_t *graph, igraph_vector_t *vector,
                                   igraph_real_t *value, igraph_bool_t scale,
                                   const igraph_vector_t *weights,
                                   igraph_arpack_options_t *options);

DECLDIR int igraph_constraint(const igraph_t *graph, igraph_vector_t *res,
                              igraph_vs_t vids, const igraph_vector_t *weights);

DECLDIR int igraph_strength(const igraph_t *graph, igraph_vector_t *res,
                            const igraph_vs_t vids, igraph_neimode_t mode,
                            igraph_bool_t loops, const igraph_vector_t *weights);

DECLDIR int igraph_convergence_degree(const igraph_t *graph, igraph_vector_t *result,
                                      igraph_vector_t *ins, igraph_vector_t *outs);

DECLDIR int igraph_sort_vertex_ids_by_degree(const igraph_t *graph,
        igraph_vector_t *outvids,
        igraph_vs_t vids,
        igraph_neimode_t mode,
        igraph_bool_t loops,
        igraph_order_t order,
        igraph_bool_t only_indices);

DECLDIR igraph_real_t igraph_centralization(const igraph_vector_t *scores,
        igraph_real_t theoretical_max,
        igraph_bool_t normalized);

DECLDIR int igraph_centralization_degree(const igraph_t *graph, igraph_vector_t *res,
        igraph_neimode_t mode, igraph_bool_t loops,
        igraph_real_t *centralization,
        igraph_real_t *theoretical_max,
        igraph_bool_t normalized);
DECLDIR int igraph_centralization_degree_tmax(const igraph_t *graph,
        igraph_integer_t nodes,
        igraph_neimode_t mode,
        igraph_bool_t loops,
        igraph_real_t *res);

DECLDIR int igraph_centralization_betweenness(const igraph_t *graph,
        igraph_vector_t *res,
        igraph_bool_t directed,
        igraph_bool_t nobigint,
        igraph_real_t *centralization,
        igraph_real_t *theoretical_max,
        igraph_bool_t normalized);
DECLDIR int igraph_centralization_betweenness_tmax(const igraph_t *graph,
        igraph_integer_t nodes,
        igraph_bool_t directed,
        igraph_real_t *res);

DECLDIR int igraph_centralization_closeness(const igraph_t *graph,
        igraph_vector_t *res,
        igraph_neimode_t mode,
        igraph_real_t *centralization,
        igraph_real_t *theoretical_max,
        igraph_bool_t normalized);
DECLDIR int igraph_centralization_closeness_tmax(const igraph_t *graph,
        igraph_integer_t nodes,
        igraph_neimode_t mode,
        igraph_real_t *res);

DECLDIR int igraph_centralization_eigenvector_centrality(
    const igraph_t *graph,
    igraph_vector_t *vector,
    igraph_real_t *value,
    igraph_bool_t directed,
    igraph_bool_t scale,
    igraph_arpack_options_t *options,
    igraph_real_t *centralization,
    igraph_real_t *theoretical_max,
    igraph_bool_t normalized);
DECLDIR int igraph_centralization_eigenvector_centrality_tmax(
    const igraph_t *graph,
    igraph_integer_t nodes,
    igraph_bool_t directed,
    igraph_bool_t scale,
    igraph_real_t *res);

__END_DECLS

#endif
