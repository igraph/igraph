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

#ifndef IGRAPH_GAMES_H
#define IGRAPH_GAMES_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_types.h"
#include "igraph_matrix.h"
#include "igraph_vector.h"
#include "igraph_datatype.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Constructors, games (=stochastic)                  */
/* -------------------------------------------------- */

IGRAPH_EXPORT int igraph_barabasi_game(igraph_t *graph, igraph_integer_t n,
                                       igraph_real_t power,
                                       igraph_integer_t m,
                                       const igraph_vector_t *outseq,
                                       igraph_bool_t outpref,
                                       igraph_real_t A,
                                       igraph_bool_t directed,
                                       igraph_barabasi_algorithm_t algo,
                                       const igraph_t *start_from);
IGRAPH_EXPORT int igraph_erdos_renyi_game(igraph_t *graph, igraph_erdos_renyi_t type,
                                          igraph_integer_t n, igraph_real_t p_or_m,
                                          igraph_bool_t directed, igraph_bool_t loops);
IGRAPH_EXPORT int igraph_erdos_renyi_game_gnp(igraph_t *graph, igraph_integer_t n, igraph_real_t p,
                                              igraph_bool_t directed, igraph_bool_t loops);
IGRAPH_EXPORT int igraph_erdos_renyi_game_gnm(igraph_t *graph, igraph_integer_t n, igraph_real_t m,
                                              igraph_bool_t directed, igraph_bool_t loops);
IGRAPH_EXPORT int igraph_degree_sequence_game(igraph_t *graph, const igraph_vector_t *out_deg,
                                              const igraph_vector_t *in_deg,
                                              igraph_degseq_t method);
IGRAPH_EXPORT int igraph_growing_random_game(igraph_t *graph, igraph_integer_t n,
                                             igraph_integer_t m, igraph_bool_t directed, igraph_bool_t citation);
IGRAPH_EXPORT int igraph_barabasi_aging_game(igraph_t *graph,
                                             igraph_integer_t nodes,
                                             igraph_integer_t m,
                                             const igraph_vector_t *outseq,
                                             igraph_bool_t outpref,
                                             igraph_real_t pa_exp,
                                             igraph_real_t aging_exp,
                                             igraph_integer_t aging_bin,
                                             igraph_real_t zero_deg_appeal,
                                             igraph_real_t zero_age_appeal,
                                             igraph_real_t deg_coef,
                                             igraph_real_t age_coef,
                                             igraph_bool_t directed);
IGRAPH_EXPORT int igraph_recent_degree_game(igraph_t *graph, igraph_integer_t n,
                                            igraph_real_t power,
                                            igraph_integer_t window,
                                            igraph_integer_t m,
                                            const igraph_vector_t *outseq,
                                            igraph_bool_t outpref,
                                            igraph_real_t zero_appeal,
                                            igraph_bool_t directed);
IGRAPH_EXPORT int igraph_recent_degree_aging_game(igraph_t *graph,
                                                  igraph_integer_t nodes,
                                                  igraph_integer_t m,
                                                  const igraph_vector_t *outseq,
                                                  igraph_bool_t outpref,
                                                  igraph_real_t pa_exp,
                                                  igraph_real_t aging_exp,
                                                  igraph_integer_t aging_bin,
                                                  igraph_integer_t window,
                                                  igraph_real_t zero_appeal,
                                                  igraph_bool_t directed);
IGRAPH_EXPORT int igraph_callaway_traits_game(igraph_t *graph, igraph_integer_t nodes,
                                              igraph_integer_t types, igraph_integer_t edges_per_step,
                                              const igraph_vector_t *type_dist,
                                              const igraph_matrix_t *pref_matrix,
                                              igraph_bool_t directed,
                                              igraph_vector_t *node_type_vec);
IGRAPH_EXPORT int igraph_establishment_game(igraph_t *graph, igraph_integer_t nodes,
                                            igraph_integer_t types, igraph_integer_t k,
                                            const igraph_vector_t *type_dist,
                                            const igraph_matrix_t *pref_matrix,
                                            igraph_bool_t directed,
                                            igraph_vector_t *node_type_vec);
IGRAPH_EXPORT int igraph_grg_game(igraph_t *graph, igraph_integer_t nodes,
                                  igraph_real_t radius, igraph_bool_t torus,
                                  igraph_vector_t *x, igraph_vector_t *y);
IGRAPH_EXPORT int igraph_preference_game(igraph_t *graph, igraph_integer_t nodes,
                                         igraph_integer_t types,
                                         const igraph_vector_t *type_dist,
                                         igraph_bool_t fixed_sizes,
                                         const igraph_matrix_t *pref_matrix,
                                         igraph_vector_t *node_type_vec,
                                         igraph_bool_t directed, igraph_bool_t loops);
IGRAPH_EXPORT int igraph_asymmetric_preference_game(igraph_t *graph, igraph_integer_t nodes,
                                                    igraph_integer_t out_types,
                                                    igraph_integer_t in_types,
                                                    const igraph_matrix_t *type_dist_matrix,
                                                    const igraph_matrix_t *pref_matrix,
                                                    igraph_vector_t *node_type_out_vec,
                                                    igraph_vector_t *node_type_in_vec,
                                                    igraph_bool_t loops);

IGRAPH_EXPORT int igraph_rewire_edges(igraph_t *graph, igraph_real_t prob,
                                      igraph_bool_t loops, igraph_bool_t multiple);
IGRAPH_EXPORT int igraph_rewire_directed_edges(igraph_t *graph, igraph_real_t prob,
                                               igraph_bool_t loops, igraph_neimode_t mode);

IGRAPH_EXPORT int igraph_watts_strogatz_game(igraph_t *graph, igraph_integer_t dim,
                                             igraph_integer_t size, igraph_integer_t nei,
                                             igraph_real_t p, igraph_bool_t loops,
                                             igraph_bool_t multiple);

IGRAPH_EXPORT int igraph_lastcit_game(igraph_t *graph,
                                      igraph_integer_t nodes, igraph_integer_t edges_per_node,
                                      igraph_integer_t agebins,
                                      const igraph_vector_t *preference, igraph_bool_t directed);

IGRAPH_EXPORT int igraph_cited_type_game(igraph_t *graph, igraph_integer_t nodes,
                                         const igraph_vector_t *types,
                                         const igraph_vector_t *pref,
                                         igraph_integer_t edges_per_step,
                                         igraph_bool_t directed);

IGRAPH_EXPORT int igraph_citing_cited_type_game(igraph_t *graph, igraph_integer_t nodes,
                                                const igraph_vector_t *types,
                                                const igraph_matrix_t *pref,
                                                igraph_integer_t edges_per_step,
                                                igraph_bool_t directed);

IGRAPH_EXPORT int igraph_forest_fire_game(igraph_t *graph, igraph_integer_t nodes,
                                          igraph_real_t fw_prob, igraph_real_t bw_factor,
                                          igraph_integer_t ambs, igraph_bool_t directed);


IGRAPH_EXPORT int igraph_simple_interconnected_islands_game(
    igraph_t *graph,
    igraph_integer_t islands_n,
    igraph_integer_t islands_size,
    igraph_real_t islands_pin,
    igraph_integer_t n_inter);

IGRAPH_EXPORT int igraph_static_fitness_game(igraph_t *graph, igraph_integer_t no_of_edges,
                                             const igraph_vector_t *fitness_out, const igraph_vector_t *fitness_in,
                                             igraph_bool_t loops, igraph_bool_t multiple);

IGRAPH_EXPORT int igraph_static_power_law_game(igraph_t *graph,
                                               igraph_integer_t no_of_nodes, igraph_integer_t no_of_edges,
                                               igraph_real_t exponent_out, igraph_real_t exponent_in,
                                               igraph_bool_t loops, igraph_bool_t multiple,
                                               igraph_bool_t finite_size_correction);

IGRAPH_EXPORT int igraph_k_regular_game(igraph_t *graph,
                                        igraph_integer_t no_of_nodes, igraph_integer_t k,
                                        igraph_bool_t directed, igraph_bool_t multiple);

IGRAPH_EXPORT int igraph_sbm_game(igraph_t *graph, igraph_integer_t n,
                                  const igraph_matrix_t *pref_matrix,
                                  const igraph_vector_int_t *block_sizes,
                                  igraph_bool_t directed, igraph_bool_t loops);

IGRAPH_EXPORT int igraph_hsbm_game(igraph_t *graph, igraph_integer_t n,
                                   igraph_integer_t m, const igraph_vector_t *rho,
                                   const igraph_matrix_t *C, igraph_real_t p);

IGRAPH_EXPORT int igraph_hsbm_list_game(igraph_t *graph, igraph_integer_t n,
                                        const igraph_vector_int_t *mlist,
                                        const igraph_vector_ptr_t *rholist,
                                        const igraph_vector_ptr_t *Clist,
                                        igraph_real_t p);

IGRAPH_EXPORT int igraph_correlated_game(const igraph_t *old_graph, igraph_t *new_graph,
                                         igraph_real_t corr, igraph_real_t p,
                                         const igraph_vector_t *permutation);

IGRAPH_EXPORT int igraph_correlated_pair_game(igraph_t *graph1, igraph_t *graph2,
                                              igraph_integer_t n, igraph_real_t corr, igraph_real_t p,
                                              igraph_bool_t directed,
                                              const igraph_vector_t *permutation);

IGRAPH_EXPORT int igraph_tree_game(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed,
                                   igraph_random_tree_t method);

IGRAPH_EXPORT int igraph_dot_product_game(igraph_t *graph, const igraph_matrix_t *vecs,
                                          igraph_bool_t directed);

IGRAPH_EXPORT int igraph_sample_sphere_surface(igraph_integer_t dim, igraph_integer_t n,
                                               igraph_real_t radius,
                                               igraph_bool_t positive,
                                               igraph_matrix_t *res);

IGRAPH_EXPORT int igraph_sample_sphere_volume(igraph_integer_t dim, igraph_integer_t n,
                                              igraph_real_t radius,
                                              igraph_bool_t positive,
                                              igraph_matrix_t *res);

IGRAPH_EXPORT int igraph_sample_dirichlet(igraph_integer_t n, const igraph_vector_t *alpha,
                                          igraph_matrix_t *res);

__END_DECLS

#endif
