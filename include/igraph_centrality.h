/*
   igraph library.
   Copyright (C) 2009-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_CENTRALITY_H
#define IGRAPH_CENTRALITY_H

#include "igraph_decls.h"
#include "igraph_arpack.h"
#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_iterators.h"
#include "igraph_error.h"
#include "igraph_types.h"

IGRAPH_BEGIN_C_DECLS

/* -------------------------------------------------- */
/* Centrality                                         */
/* -------------------------------------------------- */

IGRAPH_EXPORT igraph_error_t igraph_closeness(const igraph_t *graph, igraph_vector_t *res,
                                   igraph_vector_int_t *reachable_count, igraph_bool_t *all_reachable,
                                   igraph_vs_t vids, igraph_neimode_t mode,
                                   const igraph_vector_t *weights, igraph_bool_t normalized);
IGRAPH_EXPORT igraph_error_t igraph_closeness_cutoff(const igraph_t *graph, igraph_vector_t *res,
                                          igraph_vector_int_t *reachable_count, igraph_bool_t *all_reachable,
                                          igraph_vs_t vids, igraph_neimode_t mode,
                                          const igraph_vector_t *weights,
                                          igraph_bool_t normalized,
                                          igraph_real_t cutoff);

IGRAPH_EXPORT igraph_error_t igraph_harmonic_centrality(const igraph_t *graph, igraph_vector_t *res,
                                             igraph_vs_t vids, igraph_neimode_t mode,
                                             const igraph_vector_t *weights,
                                             igraph_bool_t normalized);
IGRAPH_EXPORT igraph_error_t igraph_harmonic_centrality_cutoff(const igraph_t *graph, igraph_vector_t *res,
                                                    igraph_vs_t vids, igraph_neimode_t mode,
                                                    const igraph_vector_t *weights,
                                                    igraph_bool_t normalized,
                                                    igraph_real_t cutoff);

IGRAPH_EXPORT igraph_error_t igraph_betweenness(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_vector_t *res,
        igraph_vs_t vids,
        igraph_bool_t directed, igraph_bool_t normalized);

IGRAPH_EXPORT igraph_error_t igraph_betweenness_cutoff(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_vector_t *res,
        igraph_vs_t vids,
        igraph_bool_t directed, igraph_bool_t normalized,
        igraph_real_t cutoff);

IGRAPH_EXPORT igraph_error_t igraph_edge_betweenness(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_vector_t *res, igraph_es_t eids,
        igraph_bool_t directed, igraph_bool_t normalized);

IGRAPH_EXPORT igraph_error_t igraph_edge_betweenness_cutoff(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_vector_t *res, igraph_es_t eids,
        igraph_bool_t directed, igraph_bool_t normalized,
        igraph_real_t cutoff);

IGRAPH_EXPORT igraph_error_t igraph_betweenness_subset(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_vector_t *res,
        igraph_vs_t sources, igraph_vs_t targets,
        igraph_vs_t vids,
        igraph_bool_t directed, igraph_bool_t normalized);

IGRAPH_EXPORT igraph_error_t igraph_edge_betweenness_subset(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_vector_t *res,
        igraph_vs_t sources, igraph_vs_t targets,
        igraph_es_t eids,
        igraph_bool_t directed, igraph_bool_t normalized);

/**
 * \typedef igraph_pagerank_algo_t
 * \brief PageRank algorithm implementation.
 *
 * Algorithms to calculate PageRank.
 * \enumval IGRAPH_PAGERANK_ALGO_ARPACK Use the ARPACK library, this
 *   was the PageRank implementation in igraph from version 0.5, until
 *   version 0.7.
 * \enumval IGRAPH_PAGERANK_ALGO_PRPACK Use the PRPACK
 *   library. Currently this implementation is recommended.
 */

typedef enum {
    IGRAPH_PAGERANK_ALGO_ARPACK = 1,
    IGRAPH_PAGERANK_ALGO_PRPACK = 2
} igraph_pagerank_algo_t;

IGRAPH_EXPORT igraph_error_t igraph_pagerank(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_vector_t *vector, igraph_real_t *value,
        igraph_real_t damping, igraph_bool_t directed,
        igraph_vs_t vids,
        igraph_pagerank_algo_t algo,
        igraph_arpack_options_t *options);

IGRAPH_EXPORT igraph_error_t igraph_personalized_pagerank(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_vector_t *vector, igraph_real_t *value,
        const igraph_vector_t *reset,
        igraph_real_t damping, igraph_bool_t directed,
        igraph_vs_t vids,
        igraph_pagerank_algo_t algo,
        igraph_arpack_options_t *options);

IGRAPH_EXPORT igraph_error_t igraph_personalized_pagerank_vs(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_vector_t *vector, igraph_real_t *value,
        igraph_vs_t reset_vids, igraph_real_t damping,
        igraph_bool_t directed, igraph_vs_t vids,
        igraph_pagerank_algo_t algo,
        igraph_arpack_options_t *options);

IGRAPH_EXPORT igraph_error_t igraph_eigenvector_centrality(const igraph_t *graph, igraph_vector_t *vector,
                                                           igraph_real_t *value,
                                                           igraph_neimode_t mode,
                                                           const igraph_vector_t *weights,
                                                           igraph_arpack_options_t *options);

IGRAPH_EXPORT igraph_error_t igraph_hub_and_authority_scores(const igraph_t *graph, igraph_vector_t *hub_vector,
                                                             igraph_vector_t *authority_vector,
                                                             igraph_real_t *value,
                                                             const igraph_vector_t *weights,
                                                             igraph_arpack_options_t *options);

IGRAPH_EXPORT igraph_error_t igraph_constraint(const igraph_t *graph, igraph_vector_t *res,
                                    igraph_vs_t vids, const igraph_vector_t *weights);

IGRAPH_EXPORT igraph_error_t igraph_convergence_degree(const igraph_t *graph, igraph_vector_t *result,
                                            igraph_vector_t *ins, igraph_vector_t *outs);

IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_real_t igraph_centralization(const igraph_vector_t *scores,
                                                  igraph_real_t theoretical_max,
                                                  igraph_bool_t normalized);

IGRAPH_EXPORT igraph_error_t igraph_centralization_degree(const igraph_t *graph, igraph_vector_t *res,
                                               igraph_neimode_t mode, igraph_loops_t loops,
                                               igraph_real_t *centralization,
                                               igraph_real_t *theoretical_max,
                                               igraph_bool_t normalized);
IGRAPH_EXPORT igraph_error_t igraph_centralization_degree_tmax(const igraph_t *graph,
                                                    igraph_int_t nodes,
                                                    igraph_neimode_t mode,
                                                    igraph_loops_t loops,
                                                    igraph_real_t *res);

IGRAPH_EXPORT igraph_error_t igraph_centralization_betweenness(const igraph_t *graph,
                                                    igraph_vector_t *res,
                                                    igraph_bool_t directed,
                                                    igraph_real_t *centralization,
                                                    igraph_real_t *theoretical_max,
                                                    igraph_bool_t normalized);
IGRAPH_EXPORT igraph_error_t igraph_centralization_betweenness_tmax(const igraph_t *graph,
                                                         igraph_int_t nodes,
                                                         igraph_bool_t directed,
                                                         igraph_real_t *res);

IGRAPH_EXPORT igraph_error_t igraph_centralization_closeness(const igraph_t *graph,
                                                  igraph_vector_t *res,
                                                  igraph_neimode_t mode,
                                                  igraph_real_t *centralization,
                                                  igraph_real_t *theoretical_max,
                                                  igraph_bool_t normalized);
IGRAPH_EXPORT igraph_error_t igraph_centralization_closeness_tmax(const igraph_t *graph,
                                                       igraph_int_t nodes,
                                                       igraph_neimode_t mode,
                                                       igraph_real_t *res);

IGRAPH_EXPORT igraph_error_t igraph_centralization_eigenvector_centrality(const igraph_t *graph,
                                                                          igraph_vector_t *vector,
                                                                          igraph_real_t *value,
                                                                          igraph_neimode_t mode,
                                                                          igraph_arpack_options_t *options,
                                                                          igraph_real_t *centralization,
                                                                          igraph_real_t *theoretical_max,
                                                                          igraph_bool_t normalized);
IGRAPH_EXPORT igraph_error_t igraph_centralization_eigenvector_centrality_tmax(const igraph_t *graph,
                                                                               igraph_int_t nodes,
                                                                               igraph_neimode_t mode,
                                                                               igraph_real_t *res);

IGRAPH_END_C_DECLS

#endif
