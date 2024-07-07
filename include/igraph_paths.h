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

#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_iterators.h"
#include "igraph_matrix.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_vector_list.h"

__BEGIN_DECLS

/**
 * \typedef igraph_astar_heuristic_func_t
 * \brief Distance estimator for A* algorithm.
 *
 * \ref igraph_get_shortest_path_astar() uses a heuristic based on a distance
 * estimate to the target vertex to guide its search, and determine
 * which vertex to try next. The heurstic function is expected to compute
 * an estimate of the distance between \p from and \p to. In order for
 * \ref igraph_get_shortest_path_astar() to find an exact shortest path,
 * the distance must not be overestimated, i.e. the heuristic function
 * must be \em admissible.
 *
 * \param result The result of the heuristic, i.e. the estimated distance.
 *    A lower value will mean this vertex will be a better candidate for
 *    exploration.
 * \param from The vertex ID of the candidate vertex will be passed here.
 * \param to The vertex ID of the endpoint of the path, i.e. the \c to parameter
 *    given to \ref igraph_get_shortest_path_astar(), will be passed here.
 * \param extra The \c extra argument that was passed to
 *    \ref igraph_get_shortest_path_astar().
 * \return Error code. Must return \c IGRAPH_SUCCESS if there were no errors.
 *    This can be used to break off the algorithm if something unexpected happens,
 *    like a failed memory allocation (\c IGRAPH_ENOMEM).
 *
 * \sa \ref igraph_get_shortest_path_astar()
 */
typedef igraph_error_t igraph_astar_heuristic_func_t(
            igraph_real_t *result,
            igraph_integer_t from, igraph_integer_t to,
            void *extra);

typedef enum {
    IGRAPH_FLOYD_WARSHALL_AUTOMATIC = 0,
    IGRAPH_FLOYD_WARSHALL_ORIGINAL = 1,
    IGRAPH_FLOYD_WARSHALL_TREE = 2
} igraph_floyd_warshall_algorithm_t;

IGRAPH_EXPORT igraph_error_t igraph_diameter(const igraph_t *graph, igraph_real_t *res,
                                  igraph_integer_t *from, igraph_integer_t *to,
                                  igraph_vector_int_t *vertex_path, igraph_vector_int_t *edge_path,
                                  igraph_bool_t directed, igraph_bool_t unconn);
IGRAPH_EXPORT igraph_error_t igraph_diameter_dijkstra(const igraph_t *graph,
                                           const igraph_vector_t *weights,
                                           igraph_real_t *res,
                                           igraph_integer_t *from,
                                           igraph_integer_t *to,
                                           igraph_vector_int_t *vertex_path,
                                           igraph_vector_int_t *edge_path,
                                           igraph_bool_t directed,
                                           igraph_bool_t unconn);

IGRAPH_EXPORT igraph_error_t igraph_distances_cutoff(const igraph_t *graph, igraph_matrix_t *res,
                                                     const igraph_vs_t from, const igraph_vs_t to,
                                                     igraph_neimode_t mode, igraph_real_t cutoff);
IGRAPH_EXPORT igraph_error_t igraph_distances(const igraph_t *graph, igraph_matrix_t *res,
                                              const igraph_vs_t from, const igraph_vs_t to,
                                              igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_distances_bellman_ford(const igraph_t *graph,
                                                     igraph_matrix_t *res,
                                                     const igraph_vs_t from,
                                                     const igraph_vs_t to,
                                                     const igraph_vector_t *weights,
                                                     igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_distances_dijkstra_cutoff(const igraph_t *graph,
                                                              igraph_matrix_t *res,
                                                              const igraph_vs_t from,
                                                              const igraph_vs_t to,
                                                              const igraph_vector_t *weights,
                                                              igraph_neimode_t mode,
                                                              igraph_real_t cutoff);
IGRAPH_EXPORT igraph_error_t igraph_distances_dijkstra(const igraph_t *graph,
                                                       igraph_matrix_t *res,
                                                       const igraph_vs_t from,
                                                       const igraph_vs_t to,
                                                       const igraph_vector_t *weights,
                                                       igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_distances_johnson(const igraph_t *graph,
                                                igraph_matrix_t *res,
                                                const igraph_vs_t from,
                                                const igraph_vs_t to,
                                                const igraph_vector_t *weights);
IGRAPH_EXPORT igraph_error_t igraph_distances_floyd_warshall(const igraph_t *graph,
                                                             igraph_matrix_t *res,
                                                             igraph_vs_t from,
                                                             igraph_vs_t to,
                                                             const igraph_vector_t *weights,
                                                             igraph_neimode_t mode,
                                                             igraph_floyd_warshall_algorithm_t method);
IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_shortest_paths(const igraph_t *graph, igraph_matrix_t *res,
                                        const igraph_vs_t from, const igraph_vs_t to,
                                        igraph_neimode_t mode);
IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_shortest_paths_bellman_ford(const igraph_t *graph,
                                                     igraph_matrix_t *res,
                                                     const igraph_vs_t from,
                                                     const igraph_vs_t to,
                                                     const igraph_vector_t *weights,
                                                     igraph_neimode_t mode);
IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_shortest_paths_dijkstra(const igraph_t *graph,
                                                 igraph_matrix_t *res,
                                                 const igraph_vs_t from,
                                                 const igraph_vs_t to,
                                                 const igraph_vector_t *weights,
                                                 igraph_neimode_t mode);
IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_shortest_paths_johnson(const igraph_t *graph,
                                                igraph_matrix_t *res,
                                                const igraph_vs_t from,
                                                const igraph_vs_t to,
                                                const igraph_vector_t *weights);

IGRAPH_EXPORT igraph_error_t igraph_get_shortest_paths(const igraph_t *graph,
                                            igraph_vector_int_list_t *vertices,
                                            igraph_vector_int_list_t *edges,
                                            igraph_integer_t from, const igraph_vs_t to,
                                            igraph_neimode_t mode,
                                            igraph_vector_int_t *parents,
                                            igraph_vector_int_t *inbound_edges);
IGRAPH_EXPORT igraph_error_t igraph_get_shortest_paths_bellman_ford(const igraph_t *graph,
                                                      igraph_vector_int_list_t *vertices,
                                                      igraph_vector_int_list_t *edges,
                                                      igraph_integer_t from,
                                                      igraph_vs_t to,
                                                      const igraph_vector_t *weights,
                                                      igraph_neimode_t mode,
                                                      igraph_vector_int_t *parents,
                                                      igraph_vector_int_t *inbound_edges);
IGRAPH_EXPORT igraph_error_t igraph_get_shortest_paths_dijkstra(const igraph_t *graph,
                                                     igraph_vector_int_list_t *vertices,
                                                     igraph_vector_int_list_t *edges,
                                                     igraph_integer_t from,
                                                     igraph_vs_t to,
                                                     const igraph_vector_t *weights,
                                                     igraph_neimode_t mode,
                                                     igraph_vector_int_t *parents,
                                                     igraph_vector_int_t *inbound_edges);

IGRAPH_EXPORT igraph_error_t igraph_get_shortest_path(const igraph_t *graph,
                                           igraph_vector_int_t *vertices,
                                           igraph_vector_int_t *edges,
                                           igraph_integer_t from,
                                           igraph_integer_t to,
                                           igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_get_shortest_path_bellman_ford(const igraph_t *graph,
                                                        igraph_vector_int_t *vertices,
                                                        igraph_vector_int_t *edges,
                                                        igraph_integer_t from,
                                                        igraph_integer_t to,
                                                        const igraph_vector_t *weights,
                                                        igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_get_shortest_path_dijkstra(const igraph_t *graph,
                                                    igraph_vector_int_t *vertices,
                                                    igraph_vector_int_t *edges,
                                                    igraph_integer_t from,
                                                    igraph_integer_t to,
                                                    const igraph_vector_t *weights,
                                                    igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_get_shortest_path_astar(const igraph_t *graph,
                                      igraph_vector_int_t *vertices,
                                      igraph_vector_int_t *edges,
                                      igraph_integer_t from,
                                      igraph_integer_t to,
                                      const igraph_vector_t *weights,
                                      igraph_neimode_t mode,
                                      igraph_astar_heuristic_func_t *heuristic,
                                      void *extra);

IGRAPH_EXPORT igraph_error_t igraph_get_all_shortest_paths(const igraph_t *graph,
                                                igraph_vector_int_list_t *vertices,
                                                igraph_vector_int_list_t *edges,
                                                igraph_vector_int_t *nrgeo,
                                                igraph_integer_t from, const igraph_vs_t to,
                                                igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_get_all_shortest_paths_dijkstra(const igraph_t *graph,
                                                         igraph_vector_int_list_t *vertices,
                                                         igraph_vector_int_list_t *edges,
                                                         igraph_vector_int_t *nrgeo,
                                                         igraph_integer_t from, igraph_vs_t to,
                                                         const igraph_vector_t *weights,
                                                         igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_average_path_length(const igraph_t *graph,
                                             igraph_real_t *res, igraph_real_t *unconn_pairs,
                                             igraph_bool_t directed, igraph_bool_t unconn);
IGRAPH_EXPORT igraph_error_t igraph_average_path_length_dijkstra(const igraph_t *graph,
                                                      igraph_real_t *res, igraph_real_t *unconn_pairs,
                                                      const igraph_vector_t *weights,
                                                      igraph_bool_t directed, igraph_bool_t unconn);
IGRAPH_EXPORT igraph_error_t igraph_path_length_hist(const igraph_t *graph, igraph_vector_t *res,
                                          igraph_real_t *unconnected, igraph_bool_t directed);

IGRAPH_EXPORT igraph_error_t igraph_global_efficiency(const igraph_t *graph, igraph_real_t *res,
                                           const igraph_vector_t *weights,
                                           igraph_bool_t directed);
IGRAPH_EXPORT igraph_error_t igraph_local_efficiency(const igraph_t *graph, igraph_vector_t *res,
                                          const igraph_vs_t vids,
                                          const igraph_vector_t *weights,
                                          igraph_bool_t directed, igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_average_local_efficiency(const igraph_t *graph, igraph_real_t *res,
                                                  const igraph_vector_t *weights,
                                                  igraph_bool_t directed, igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_eccentricity(const igraph_t *graph,
                                      igraph_vector_t *res,
                                      igraph_vs_t vids,
                                      igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_eccentricity_dijkstra(const igraph_t *graph,
                        const igraph_vector_t *weights,
                        igraph_vector_t *res,
                        igraph_vs_t vids,
                        igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_radius(const igraph_t *graph, igraph_real_t *radius,
                                igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_radius_dijkstra(const igraph_t *graph, const igraph_vector_t *weights,
                                                    igraph_real_t *radius, igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_graph_center(const igraph_t *graph,
                                 igraph_vector_int_t *res,
                                 igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_graph_center_dijkstra(
    const igraph_t *graph, const igraph_vector_t *weights,
    igraph_vector_int_t *res, igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_pseudo_diameter(const igraph_t *graph,
                                         igraph_real_t *diameter,
                                         igraph_integer_t vid_start,
                                         igraph_integer_t *from,
                                         igraph_integer_t *to,
                                         igraph_bool_t directed,
                                         igraph_bool_t unconn);
IGRAPH_EXPORT igraph_error_t igraph_pseudo_diameter_dijkstra(const igraph_t *graph,
                                                  const igraph_vector_t *weights,
                                                  igraph_real_t *diameter,
                                                  igraph_integer_t vid_start,
                                                  igraph_integer_t *from,
                                                  igraph_integer_t *to,
                                                  igraph_bool_t directed,
                                                  igraph_bool_t unconn);

IGRAPH_EXPORT igraph_error_t igraph_get_all_simple_paths(const igraph_t *graph,
                                              igraph_vector_int_t *res,
                                              igraph_integer_t from,
                                              const igraph_vs_t to,
                                              igraph_integer_t cutoff,
                                              igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_random_walk(const igraph_t *graph,
                                     const igraph_vector_t *weights,
                                     igraph_vector_int_t *vertices,
                                     igraph_vector_int_t *edges,
                                     igraph_integer_t start,
                                     igraph_neimode_t mode,
                                     igraph_integer_t steps,
                                     igraph_random_walk_stuck_t stuck);

IGRAPH_EXPORT igraph_error_t igraph_get_k_shortest_paths(const igraph_t *graph,
                                          const igraph_vector_t *weights,
                                          igraph_vector_int_list_t *vertex_paths,
                                          igraph_vector_int_list_t *edge_paths,
                                          igraph_integer_t k,
                                          igraph_integer_t from,
                                          igraph_integer_t to,
                                          igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_spanner(const igraph_t *graph,
                                igraph_vector_int_t *spanner,
                                igraph_real_t stretch,
                                const igraph_vector_t *weights);

IGRAPH_EXPORT igraph_error_t igraph_get_widest_paths(const igraph_t *graph,
                                             igraph_vector_int_list_t *vertices,
                                             igraph_vector_int_list_t *edges,
                                             igraph_integer_t from,
                                             igraph_vs_t to,
                                             const igraph_vector_t *weights,
                                             igraph_neimode_t mode,
                                             igraph_vector_int_t *parents,
                                             igraph_vector_int_t *inbound_edges);
IGRAPH_EXPORT igraph_error_t igraph_get_widest_path(const igraph_t *graph,
                                             igraph_vector_int_t *vertices,
                                             igraph_vector_int_t *edges,
                                             igraph_integer_t from,
                                             igraph_integer_t to,
                                             const igraph_vector_t *weights,
                                             igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_widest_path_widths_floyd_warshall(const igraph_t *graph,
                                                   igraph_matrix_t *res,
                                                   const igraph_vs_t from,
                                                   const igraph_vs_t to,
                                                   const igraph_vector_t *weights,
                                                   igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_widest_path_widths_dijkstra(const igraph_t *graph,
                                             igraph_matrix_t *res,
                                             const igraph_vs_t from,
                                             const igraph_vs_t to,
                                             const igraph_vector_t *weights,
                                             igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_voronoi(const igraph_t *graph,
                                            igraph_vector_int_t *membership,
                                            igraph_vector_t *distances,
                                            const igraph_vector_int_t *generators,
                                            const igraph_vector_t *weights,
                                            igraph_neimode_t mode,
                                            igraph_voronoi_tiebreaker_t tiebreaker);

IGRAPH_EXPORT igraph_error_t igraph_expand_path_to_pairs(igraph_vector_int_t *path);

IGRAPH_EXPORT igraph_error_t igraph_vertex_path_from_edge_path(
        const igraph_t *graph, igraph_integer_t start,
        const igraph_vector_int_t *edge_path, igraph_vector_int_t *vertex_path,
        igraph_neimode_t mode);

/* Deprecated functions: */

IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_random_edge_walk(const igraph_t *graph,
                                                            const igraph_vector_t *weights,
                                                            igraph_vector_int_t *edgewalk,
                                                            igraph_integer_t start,
                                                            igraph_neimode_t mode,
                                                            igraph_integer_t steps,
                                                            igraph_random_walk_stuck_t stuck);

__END_DECLS

#endif
