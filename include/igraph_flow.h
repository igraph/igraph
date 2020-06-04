/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_FLOW_H
#define IGRAPH_FLOW_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_types.h"
#include "igraph_datatype.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* MAximum flows, minimum cuts & such                 */
/* -------------------------------------------------- */

/**
 * \typedef igraph_maxflow_stats_t
 * A simple data type to return some statistics from the
 * push-relabel maximum flow solver.
 *
 * \param nopush The number of push operations performed.
 * \param norelabel The number of relabel operarions performed.
 * \param nogap The number of times the gap heuristics was used.
 * \param nogapnodes The total number of vertices that were
 *        omitted form further calculations because of the gap
 *        heuristics.
 * \param nobfs The number of times the reverse BFS was run to
 *        assign good values to the height function. This includes
 *        an initial run before the whole algorithm, so it is always
 *        at least one.
 */

typedef struct {
    int nopush, norelabel, nogap, nogapnodes, nobfs;
} igraph_maxflow_stats_t;

DECLDIR int igraph_maxflow(const igraph_t *graph, igraph_real_t *value,
                           igraph_vector_t *flow, igraph_vector_t *cut,
                           igraph_vector_t *partition, igraph_vector_t *partition2,
                           igraph_integer_t source, igraph_integer_t target,
                           const igraph_vector_t *capacity,
                           igraph_maxflow_stats_t *stats);
DECLDIR int igraph_maxflow_value(const igraph_t *graph, igraph_real_t *value,
                                 igraph_integer_t source, igraph_integer_t target,
                                 const igraph_vector_t *capacity,
                                 igraph_maxflow_stats_t *stats);

DECLDIR int igraph_st_mincut(const igraph_t *graph, igraph_real_t *value,
                             igraph_vector_t *cut, igraph_vector_t *partition,
                             igraph_vector_t *partition2,
                             igraph_integer_t source, igraph_integer_t target,
                             const igraph_vector_t *capacity);
DECLDIR int igraph_st_mincut_value(const igraph_t *graph, igraph_real_t *res,
                                   igraph_integer_t source, igraph_integer_t target,
                                   const igraph_vector_t *capacity);

DECLDIR int igraph_mincut_value(const igraph_t *graph, igraph_real_t *res,
                                const igraph_vector_t *capacity);
DECLDIR int igraph_mincut(const igraph_t *graph,
                          igraph_real_t *value,
                          igraph_vector_t *partition,
                          igraph_vector_t *partition2,
                          igraph_vector_t *cut,
                          const igraph_vector_t *capacity);

DECLDIR int igraph_st_vertex_connectivity(const igraph_t *graph,
        igraph_integer_t *res,
        igraph_integer_t source,
        igraph_integer_t target,
        igraph_vconn_nei_t neighbors);
DECLDIR int igraph_vertex_connectivity(const igraph_t *graph, igraph_integer_t *res,
                                       igraph_bool_t checks);

DECLDIR int igraph_st_edge_connectivity(const igraph_t *graph, igraph_integer_t *res,
                                        igraph_integer_t source,
                                        igraph_integer_t target);
DECLDIR int igraph_edge_connectivity(const igraph_t *graph, igraph_integer_t *res,
                                     igraph_bool_t checks);

DECLDIR int igraph_edge_disjoint_paths(const igraph_t *graph, igraph_integer_t *res,
                                       igraph_integer_t source,
                                       igraph_integer_t target);
DECLDIR int igraph_vertex_disjoint_paths(const igraph_t *graph, igraph_integer_t *res,
        igraph_integer_t source,
        igraph_integer_t target);

DECLDIR int igraph_adhesion(const igraph_t *graph, igraph_integer_t *res,
                            igraph_bool_t checks);
DECLDIR int igraph_cohesion(const igraph_t *graph, igraph_integer_t *res,
                            igraph_bool_t checks);

/* s-t cut listing related stuff */

DECLDIR int igraph_even_tarjan_reduction(const igraph_t *graph, igraph_t *graphbar,
        igraph_vector_t *capacity);

DECLDIR int igraph_residual_graph(const igraph_t *graph,
                                  const igraph_vector_t *capacity,
                                  igraph_t *residual,
                                  igraph_vector_t *residual_capacity,
                                  const igraph_vector_t *flow);
int igraph_i_residual_graph(const igraph_t *graph,
                            const igraph_vector_t *capacity,
                            igraph_t *residual,
                            igraph_vector_t *residual_capacity,
                            const igraph_vector_t *flow,
                            igraph_vector_t *tmp);

int igraph_i_reverse_residual_graph(const igraph_t *graph,
                                    const igraph_vector_t *capacity,
                                    igraph_t *residual,
                                    const igraph_vector_t *flow,
                                    igraph_vector_t *tmp);
DECLDIR int igraph_reverse_residual_graph(const igraph_t *graph,
        const igraph_vector_t *capacity,
        igraph_t *residual,
        const igraph_vector_t *flow);

DECLDIR int igraph_dominator_tree(const igraph_t *graph,
                                  igraph_integer_t root,
                                  igraph_vector_t *dom,
                                  igraph_t *domtree,
                                  igraph_vector_t *leftout,
                                  igraph_neimode_t mode);

DECLDIR int igraph_all_st_cuts(const igraph_t *graph,
                               igraph_vector_ptr_t *cuts,
                               igraph_vector_ptr_t *partition1s,
                               igraph_integer_t source,
                               igraph_integer_t target);

DECLDIR int igraph_all_st_mincuts(const igraph_t *graph, igraph_real_t *value,
                                  igraph_vector_ptr_t *cuts,
                                  igraph_vector_ptr_t *partition1s,
                                  igraph_integer_t source,
                                  igraph_integer_t target,
                                  const igraph_vector_t *capacity);

DECLDIR int igraph_gomory_hu_tree(const igraph_t *graph,
                                  igraph_t *tree,
                                  igraph_vector_t *flows,
                                  const igraph_vector_t *capacity);

__END_DECLS

#endif
