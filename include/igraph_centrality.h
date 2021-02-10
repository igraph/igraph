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

IGRAPH_EXPORT int igraph_closeness(const igraph_t *graph, igraph_vector_t *res,
                                   igraph_vector_t *reachable_count, igraph_bool_t *all_reachable,
                                   const igraph_vs_t vids, igraph_neimode_t mode,
                                   const igraph_vector_t *weights, igraph_bool_t normalized);
IGRAPH_EXPORT int igraph_closeness_cutoff(const igraph_t *graph, igraph_vector_t *res,
                                          igraph_vector_t *reachable_count, igraph_bool_t *all_reachable,
                                          const igraph_vs_t vids, igraph_neimode_t mode,
                                          const igraph_vector_t *weights,
                                          igraph_bool_t normalized,
                                          igraph_real_t cutoff);

IGRAPH_EXPORT int igraph_harmonic_centrality(const igraph_t *graph, igraph_vector_t *res,
                                             const igraph_vs_t vids, igraph_neimode_t mode,
                                             const igraph_vector_t *weights,
                                             igraph_bool_t normalized);
IGRAPH_EXPORT int igraph_harmonic_centrality_cutoff(const igraph_t *graph, igraph_vector_t *res,
                                                    const igraph_vs_t vids, igraph_neimode_t mode,
                                                    const igraph_vector_t *weights,
                                                    igraph_bool_t normalized,
                                                    igraph_real_t cutoff);

IGRAPH_EXPORT int igraph_betweenness(const igraph_t *graph, igraph_vector_t *res,
                                     const igraph_vs_t vids, igraph_bool_t directed,
                                     const igraph_vector_t *weights);
IGRAPH_EXPORT int igraph_betweenness_cutoff(const igraph_t *graph, igraph_vector_t *res,
                                            const igraph_vs_t vids, igraph_bool_t directed,
                                            const igraph_vector_t *weights, igraph_real_t cutoff);
IGRAPH_EXPORT int igraph_edge_betweenness(const igraph_t *graph, igraph_vector_t *result,
                                          igraph_bool_t directed,
                                          const igraph_vector_t *weigths);
IGRAPH_EXPORT int igraph_edge_betweenness_cutoff(const igraph_t *graph, igraph_vector_t *result,
                                                 igraph_bool_t directed,
                                                 const igraph_vector_t *weights, igraph_real_t cutoff);

/**
 * \typedef igraph_pagerank_algo_t
 * \brief PageRank algorithm implementation
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

IGRAPH_EXPORT int igraph_pagerank(const igraph_t *graph, igraph_pagerank_algo_t algo,
                                  igraph_vector_t *vector,
                                  igraph_real_t *value, const igraph_vs_t vids,
                                  igraph_bool_t directed, igraph_real_t damping,
                                  const igraph_vector_t *weights, igraph_arpack_options_t *options);
IGRAPH_EXPORT int igraph_personalized_pagerank(const igraph_t *graph,
                                               igraph_pagerank_algo_t algo, igraph_vector_t *vector,
                                               igraph_real_t *value, const igraph_vs_t vids,
                                               igraph_bool_t directed, igraph_real_t damping,
                                               const igraph_vector_t *reset,
                                               const igraph_vector_t *weights, igraph_arpack_options_t *options);
IGRAPH_EXPORT int igraph_personalized_pagerank_vs(const igraph_t *graph,
                                                  igraph_pagerank_algo_t algo,
                                                  igraph_vector_t *vector,
                                                  igraph_real_t *value, const igraph_vs_t vids,
                                                  igraph_bool_t directed, igraph_real_t damping,
                                                  igraph_vs_t reset_vids,
                                                  const igraph_vector_t *weights, igraph_arpack_options_t *options);

IGRAPH_EXPORT int igraph_eigenvector_centrality(const igraph_t *graph, igraph_vector_t *vector,
                                                igraph_real_t *value,
                                                igraph_bool_t directed, igraph_bool_t scale,
                                                const igraph_vector_t *weights,
                                                igraph_arpack_options_t *options);

IGRAPH_EXPORT int igraph_hub_score(const igraph_t *graph, igraph_vector_t *vector,
                                   igraph_real_t *value, igraph_bool_t scale,
                                   const igraph_vector_t *weights,
                                   igraph_arpack_options_t *options);
IGRAPH_EXPORT int igraph_authority_score(const igraph_t *graph, igraph_vector_t *vector,
                                         igraph_real_t *value, igraph_bool_t scale,
                                         const igraph_vector_t *weights,
                                         igraph_arpack_options_t *options);

IGRAPH_EXPORT int igraph_constraint(const igraph_t *graph, igraph_vector_t *res,
                                    igraph_vs_t vids, const igraph_vector_t *weights);

IGRAPH_EXPORT int igraph_convergence_degree(const igraph_t *graph, igraph_vector_t *result,
                                            igraph_vector_t *ins, igraph_vector_t *outs);

IGRAPH_EXPORT igraph_real_t igraph_centralization(const igraph_vector_t *scores,
                                                  igraph_real_t theoretical_max,
                                                  igraph_bool_t normalized);

IGRAPH_EXPORT int igraph_centralization_degree(const igraph_t *graph, igraph_vector_t *res,
                                               igraph_neimode_t mode, igraph_bool_t loops,
                                               igraph_real_t *centralization,
                                               igraph_real_t *theoretical_max,
                                               igraph_bool_t normalized);
IGRAPH_EXPORT int igraph_centralization_degree_tmax(const igraph_t *graph,
                                                    igraph_integer_t nodes,
                                                    igraph_neimode_t mode,
                                                    igraph_bool_t loops,
                                                    igraph_real_t *res);

IGRAPH_EXPORT int igraph_centralization_betweenness(const igraph_t *graph,
                                                    igraph_vector_t *res,
                                                    igraph_bool_t directed,
                                                    igraph_real_t *centralization,
                                                    igraph_real_t *theoretical_max,
                                                    igraph_bool_t normalized);
IGRAPH_EXPORT int igraph_centralization_betweenness_tmax(const igraph_t *graph,
                                                         igraph_integer_t nodes,
                                                         igraph_bool_t directed,
                                                         igraph_real_t *res);

IGRAPH_EXPORT int igraph_centralization_closeness(const igraph_t *graph,
                                                  igraph_vector_t *res,
                                                  igraph_neimode_t mode,
                                                  igraph_real_t *centralization,
                                                  igraph_real_t *theoretical_max,
                                                  igraph_bool_t normalized);
IGRAPH_EXPORT int igraph_centralization_closeness_tmax(const igraph_t *graph,
                                                       igraph_integer_t nodes,
                                                       igraph_neimode_t mode,
                                                       igraph_real_t *res);

IGRAPH_EXPORT int igraph_centralization_eigenvector_centrality(
    const igraph_t *graph,
    igraph_vector_t *vector,
    igraph_real_t *value,
    igraph_bool_t directed,
    igraph_bool_t scale,
    igraph_arpack_options_t *options,
    igraph_real_t *centralization,
    igraph_real_t *theoretical_max,
    igraph_bool_t normalized);
IGRAPH_EXPORT int igraph_centralization_eigenvector_centrality_tmax(
    const igraph_t *graph,
    igraph_integer_t nodes,
    igraph_bool_t directed,
    igraph_bool_t scale,
    igraph_real_t *res);


/* Deprecated functions: */

IGRAPH_EXPORT IGRAPH_DEPRECATED int igraph_closeness_estimate(const igraph_t *graph, igraph_vector_t *res,
                                                              const igraph_vs_t vids, igraph_neimode_t mode,
                                                              igraph_real_t cutoff,
                                                              const igraph_vector_t *weights,
                                                              igraph_bool_t normalized);

IGRAPH_EXPORT IGRAPH_DEPRECATED int igraph_betweenness_estimate(const igraph_t *graph, igraph_vector_t *res,
                                                                const igraph_vs_t vids, igraph_bool_t directed,
                                                                igraph_real_t cutoff, const igraph_vector_t *weights);

IGRAPH_EXPORT IGRAPH_DEPRECATED int igraph_edge_betweenness_estimate(const igraph_t *graph, igraph_vector_t *result,
                                                                     igraph_bool_t directed, igraph_real_t cutoff,
                                                                     const igraph_vector_t *weights);

__END_DECLS

#endif
