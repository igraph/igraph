/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2021  The igraph development team

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

#ifndef IGRAPH_PATHS_H
#define IGRAPH_PATHS_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_vector_ptr.h"
#include "igraph_matrix.h"
#include "igraph_iterators.h"

__BEGIN_DECLS

IGRAPH_EXPORT int igraph_diameter(const igraph_t *graph, igraph_real_t *res,
                                  igraph_integer_t *from, igraph_integer_t *to,
                                  igraph_vector_t *path,
                                  igraph_bool_t directed, igraph_bool_t unconn);
IGRAPH_EXPORT int igraph_diameter_dijkstra(const igraph_t *graph,
                                           const igraph_vector_t *weights,
                                           igraph_real_t *pres,
                                           igraph_integer_t *pfrom,
                                           igraph_integer_t *pto,
                                           igraph_vector_t *path,
                                           igraph_bool_t directed,
                                           igraph_bool_t unconn);

IGRAPH_EXPORT int igraph_shortest_paths(const igraph_t *graph, igraph_matrix_t *res,
                                        const igraph_vs_t from, const igraph_vs_t to,
                                        igraph_neimode_t mode);
IGRAPH_EXPORT int igraph_get_shortest_paths(const igraph_t *graph,
                                            igraph_vector_ptr_t *vertices,
                                            igraph_vector_ptr_t *edges,
                                            igraph_integer_t from, const igraph_vs_t to,
                                            igraph_neimode_t mode,
                                            igraph_vector_long_t *predecessors,
                                            igraph_vector_long_t *inbound_edges);
IGRAPH_EXPORT int igraph_get_shortest_path(const igraph_t *graph,
                                           igraph_vector_t *vertices,
                                           igraph_vector_t *edges,
                                           igraph_integer_t from,
                                           igraph_integer_t to,
                                           igraph_neimode_t mode);

IGRAPH_EXPORT int igraph_get_all_shortest_paths(const igraph_t *graph,
                                                igraph_vector_ptr_t *res,
                                                igraph_vector_t *nrgeo,
                                                igraph_integer_t from, const igraph_vs_t to,
                                                igraph_neimode_t mode);
IGRAPH_EXPORT int igraph_shortest_paths_dijkstra(const igraph_t *graph,
                                                 igraph_matrix_t *res,
                                                 const igraph_vs_t from,
                                                 const igraph_vs_t to,
                                                 const igraph_vector_t *weights,
                                                 igraph_neimode_t mode);
IGRAPH_EXPORT int igraph_shortest_paths_bellman_ford(const igraph_t *graph,
                                                     igraph_matrix_t *res,
                                                     const igraph_vs_t from,
                                                     const igraph_vs_t to,
                                                     const igraph_vector_t *weights,
                                                     igraph_neimode_t mode);
IGRAPH_EXPORT int igraph_get_shortest_path_bellman_ford(const igraph_t *graph,
                                                        igraph_vector_t *vertices,
                                                        igraph_vector_t *edges,
                                                        igraph_integer_t from,
                                                        igraph_integer_t to,
                                                        const igraph_vector_t *weights,
                                                        igraph_neimode_t mode);
IGRAPH_EXPORT int igraph_get_shortest_paths_dijkstra(const igraph_t *graph,
                                                     igraph_vector_ptr_t *vertices,
                                                     igraph_vector_ptr_t *edges,
                                                     igraph_integer_t from,
                                                     igraph_vs_t to,
                                                     const igraph_vector_t *weights,
                                                     igraph_neimode_t mode,
                                                     igraph_vector_long_t *predecessors,
                                                     igraph_vector_long_t *inbound_edges);
IGRAPH_EXPORT int igraph_get_shortest_paths_bellman_ford(const igraph_t *graph,
                                                      igraph_vector_ptr_t *vertices,
                                                      igraph_vector_ptr_t *edges,
                                                      igraph_integer_t from,
                                                      igraph_vs_t to,
                                                      const igraph_vector_t *weights,
                                                      igraph_neimode_t mode,
                                                      igraph_vector_long_t *predecessors,
                                                      igraph_vector_long_t *inbound_edges);
IGRAPH_EXPORT int igraph_get_shortest_path_dijkstra(const igraph_t *graph,
                                                    igraph_vector_t *vertices,
                                                    igraph_vector_t *edges,
                                                    igraph_integer_t from,
                                                    igraph_integer_t to,
                                                    const igraph_vector_t *weights,
                                                    igraph_neimode_t mode);
IGRAPH_EXPORT int igraph_get_all_shortest_paths_dijkstra(const igraph_t *graph,
                                                         igraph_vector_ptr_t *res,
                                                         igraph_vector_t *nrgeo,
                                                         igraph_integer_t from, igraph_vs_t to,
                                                         const igraph_vector_t *weights,
                                                         igraph_neimode_t mode);
IGRAPH_EXPORT int igraph_shortest_paths_johnson(const igraph_t *graph,
                                                igraph_matrix_t *res,
                                                const igraph_vs_t from,
                                                const igraph_vs_t to,
                                                const igraph_vector_t *weights);

IGRAPH_EXPORT int igraph_average_path_length(const igraph_t *graph,
                                             igraph_real_t *res, igraph_real_t *unconn_pairs,
                                             igraph_bool_t directed, igraph_bool_t unconn);
IGRAPH_EXPORT int igraph_average_path_length_dijkstra(const igraph_t *graph,
                                                      igraph_real_t *res, igraph_real_t *unconn_pairs,
                                                      const igraph_vector_t *weights,
                                                      igraph_bool_t directed, igraph_bool_t unconn);
IGRAPH_EXPORT int igraph_path_length_hist(const igraph_t *graph, igraph_vector_t *res,
                                          igraph_real_t *unconnected, igraph_bool_t directed);

IGRAPH_EXPORT int igraph_global_efficiency(const igraph_t *graph, igraph_real_t *res,
                                           const igraph_vector_t *weights,
                                           igraph_bool_t directed);
IGRAPH_EXPORT int igraph_local_efficiency(const igraph_t *graph, igraph_vector_t *res,
                                          const igraph_vs_t vids,
                                          const igraph_vector_t *weights,
                                          igraph_bool_t directed, igraph_neimode_t mode);
IGRAPH_EXPORT int igraph_average_local_efficiency(const igraph_t *graph, igraph_real_t *res,
                                                  const igraph_vector_t *weights,
                                                  igraph_bool_t directed, igraph_neimode_t mode);

IGRAPH_EXPORT int igraph_eccentricity(const igraph_t *graph,
                                      igraph_vector_t *res,
                                      igraph_vs_t vids,
                                      igraph_neimode_t mode);

IGRAPH_EXPORT int igraph_radius(const igraph_t *graph, igraph_real_t *radius,
                                igraph_neimode_t mode);

IGRAPH_EXPORT int igraph_get_all_simple_paths(const igraph_t *graph,
                                              igraph_vector_int_t *res,
                                              igraph_integer_t from,
                                              const igraph_vs_t to,
                                              igraph_integer_t cutoff,
                                              igraph_neimode_t mode);

IGRAPH_EXPORT int igraph_random_walk(const igraph_t *graph, igraph_vector_t *walk,
                                     igraph_integer_t start, igraph_neimode_t mode,
                                     igraph_integer_t steps,
                                     igraph_random_walk_stuck_t stuck);

IGRAPH_EXPORT int igraph_random_edge_walk(const igraph_t *graph,
                                          const igraph_vector_t *weights,
                                          igraph_vector_t *edgewalk,
                                          igraph_integer_t start, igraph_neimode_t mode,
                                          igraph_integer_t steps,
                                          igraph_random_walk_stuck_t stuck);

__END_DECLS

#endif
