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

#ifndef IGRAPH_COMMUNITY_H
#define IGRAPH_COMMUNITY_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_types.h"
#include "igraph_arpack.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* K-Cores                                            */
/* -------------------------------------------------- */

IGRAPH_EXPORT int igraph_coreness(const igraph_t *graph, igraph_vector_t *cores,
                                  igraph_neimode_t mode);

/* -------------------------------------------------- */
/* Community Structure                                */
/* -------------------------------------------------- */

/* TODO: cut.community */
/* TODO: edge.type.matrix */
/* TODO:  */

IGRAPH_EXPORT int igraph_community_optimal_modularity(const igraph_t *graph,
                                                      igraph_real_t *modularity,
                                                      igraph_vector_t *membership,
                                                      const igraph_vector_t *weights);

IGRAPH_EXPORT int igraph_community_spinglass(const igraph_t *graph,
                                             const igraph_vector_t *weights,
                                             igraph_real_t *modularity,
                                             igraph_real_t *temperature,
                                             igraph_vector_t *membership,
                                             igraph_vector_t *csize,
                                             igraph_integer_t spins,
                                             igraph_bool_t parupdate,
                                             igraph_real_t starttemp,
                                             igraph_real_t stoptemp,
                                             igraph_real_t coolfact,
                                             igraph_spincomm_update_t update_rule,
                                             igraph_real_t gamma,
                                       /* the rest is for the NegSpin implementation */
                                       igraph_spinglass_implementation_t implementation,
                                       /*                    igraph_matrix_t *adhesion, */
                                       /*                    igraph_matrix_t *normalised_adhesion, */
                                       /*                    igraph_real_t *polarization, */
                                       igraph_real_t lambda);

IGRAPH_EXPORT int igraph_community_spinglass_single(const igraph_t *graph,
                                                    const igraph_vector_t *weights,
                                                    igraph_integer_t vertex,
                                                    igraph_vector_t *community,
                                                    igraph_real_t *cohesion,
                                                    igraph_real_t *adhesion,
                                                    igraph_integer_t *inner_links,
                                                    igraph_integer_t *outer_links,
                                                    igraph_integer_t spins,
                                                    igraph_spincomm_update_t update_rule,
                                                    igraph_real_t gamma);

IGRAPH_EXPORT int igraph_community_walktrap(const igraph_t *graph,
                                            const igraph_vector_t *weights,
                                            int steps,
                                            igraph_matrix_t *merges,
                                            igraph_vector_t *modularity,
                                            igraph_vector_t *membership);

IGRAPH_EXPORT int igraph_community_infomap(const igraph_t * graph,
                                           const igraph_vector_t *e_weights,
                                           const igraph_vector_t *v_weights,
                                           int nb_trials,
                                           igraph_vector_t *membership,
                                           igraph_real_t *codelength);

IGRAPH_EXPORT int igraph_community_edge_betweenness(const igraph_t *graph,
                                                    igraph_vector_t *result,
                                                    igraph_vector_t *edge_betweenness,
                                                    igraph_matrix_t *merges,
                                                    igraph_vector_t *bridges,
                                                    igraph_vector_t *modularity,
                                                    igraph_vector_t *membership,
                                                    igraph_bool_t directed,
                                                    const igraph_vector_t *weights);
IGRAPH_EXPORT int igraph_community_eb_get_merges(const igraph_t *graph,
                                                 const igraph_bool_t directed,
                                                 const igraph_vector_t *edges,
                                                 const igraph_vector_t *weights,
                                                 igraph_matrix_t *merges,
                                                 igraph_vector_t *bridges,
                                                 igraph_vector_t *modularity,
                                                 igraph_vector_t *membership);

IGRAPH_EXPORT int igraph_community_fastgreedy(const igraph_t *graph,
                                              const igraph_vector_t *weights,
                                              igraph_matrix_t *merges,
                                              igraph_vector_t *modularity,
                                              igraph_vector_t *membership);

IGRAPH_EXPORT int igraph_community_to_membership(const igraph_matrix_t *merges,
                                                 igraph_integer_t nodes,
                                                 igraph_integer_t steps,
                                                 igraph_vector_t *membership,
                                                 igraph_vector_t *csize);
IGRAPH_EXPORT int igraph_le_community_to_membership(const igraph_matrix_t *merges,
                                                    igraph_integer_t steps,
                                                    igraph_vector_t *membership,
                                                    igraph_vector_t *csize);

IGRAPH_EXPORT int igraph_modularity(const igraph_t *graph,
                                    const igraph_vector_t *membership,
                                    const igraph_vector_t *weights,
                                    const igraph_real_t resolution,
                                    const igraph_bool_t directed,
                                    igraph_real_t *modularity);

IGRAPH_EXPORT int igraph_modularity_matrix(const igraph_t *graph,
                                           const igraph_vector_t *weights,
                                           const igraph_real_t resolution,
                                           igraph_matrix_t *modmat,
                                           igraph_bool_t directed);

IGRAPH_EXPORT int igraph_reindex_membership(igraph_vector_t *membership,
                                            igraph_vector_t *new_to_old,
                                            igraph_integer_t *nb_clusters);

typedef enum { IGRAPH_LEVC_HIST_SPLIT = 1,
               IGRAPH_LEVC_HIST_FAILED,
               IGRAPH_LEVC_HIST_START_FULL,
               IGRAPH_LEVC_HIST_START_GIVEN
             } igraph_leading_eigenvector_community_history_t;

/**
 * \typedef igraph_community_leading_eigenvector_callback_t
 * Callback for the leading eigenvector community finding method.
 *
 * The leading eigenvector community finding implementation in igraph
 * is able to call a callback function, after each eigenvalue
 * calculation. This callback function must be of \c
 * igraph_community_leading_eigenvector_callback_t type.
 * The following arguments are passed to the callback:
 * \param membership The actual membership vector, before recording
 *    the potential change implied by the newly found eigenvalue.
 * \param comm The id of the community that the algorithm tried to
 *    split in the last iteration. The community ids are indexed from
 *    zero here!
 * \param eigenvalue The eigenvalue the algorithm has just found.
 * \param eigenvector The eigenvector corresponding to the eigenvalue
 *    the algorithm just found.
 * \param arpack_multiplier A function that was passed to \ref
 *    igraph_arpack_rssolve() to solve the last eigenproblem.
 * \param arpack_extra The extra argument that was passed to the
 *    ARPACK solver.
 * \param extra Extra argument that as passed to \ref
 *    igraph_community_leading_eigenvector().
 *
 * \sa \ref igraph_community_leading_eigenvector(), \ref
 * igraph_arpack_function_t, \ref igraph_arpack_rssolve().
 */

typedef int igraph_community_leading_eigenvector_callback_t(
    const igraph_vector_t *membership,
    long int comm,
    igraph_real_t eigenvalue,
    const igraph_vector_t *eigenvector,
    igraph_arpack_function_t *arpack_multiplier,
    void *arpack_extra,
    void *extra);

IGRAPH_EXPORT int igraph_community_leading_eigenvector(const igraph_t *graph,
                                                       const igraph_vector_t *weights,
                                                       igraph_matrix_t *merges,
                                                       igraph_vector_t *membership,
                                                       igraph_integer_t steps,
                                                       igraph_arpack_options_t *options,
                                                       igraph_real_t *modularity,
                                                       igraph_bool_t start,
                                                       igraph_vector_t *eigenvalues,
                                                       igraph_vector_ptr_t *eigenvectors,
                                                       igraph_vector_t *history,
                                                       igraph_community_leading_eigenvector_callback_t *callback,
                                                       void *callback_extra);

IGRAPH_EXPORT int igraph_community_fluid_communities(const igraph_t *graph,
                                                     igraph_integer_t no_of_communities,
                                                     igraph_vector_t *membership,
                                                     igraph_real_t *modularity);

IGRAPH_EXPORT int igraph_community_label_propagation(const igraph_t *graph,
                                                     igraph_vector_t *membership,
                                                     const igraph_vector_t *weights,
                                                     const igraph_vector_t *initial,
                                                     igraph_vector_bool_t *fixed,
                                                     igraph_real_t *modularity);

IGRAPH_EXPORT int igraph_community_multilevel(const igraph_t *graph,
                                              const igraph_vector_t *weights,
                                              const igraph_real_t resolution,
                                              igraph_vector_t *membership,
                                              igraph_matrix_t *memberships,
                                              igraph_vector_t *modularity);

IGRAPH_EXPORT int igraph_community_leiden(const igraph_t *graph,
                                          const igraph_vector_t *edge_weights,
                                          const igraph_vector_t *node_weights,
                                          const igraph_real_t resolution_parameter,
                                          const igraph_real_t beta,
                                          const igraph_bool_t start,
                                          igraph_vector_t *membership,
                                          igraph_integer_t *nb_clusters,
                                          igraph_real_t *quality);
/* -------------------------------------------------- */
/* Community Structure Comparison                     */
/* -------------------------------------------------- */

IGRAPH_EXPORT int igraph_compare_communities(const igraph_vector_t *comm1,
                                             const igraph_vector_t *comm2,
                                             igraph_real_t* result,
                                             igraph_community_comparison_t method);
IGRAPH_EXPORT int igraph_split_join_distance(const igraph_vector_t *comm1,
                                             const igraph_vector_t *comm2,
                                             igraph_integer_t* distance12,
                                             igraph_integer_t* distance21);

__END_DECLS

#endif
