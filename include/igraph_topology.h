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

#ifndef IGRAPH_TOPOLOGY_H
#define IGRAPH_TOPOLOGY_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector_list.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Directed acyclic graphs                            */
/* -------------------------------------------------- */

IGRAPH_EXPORT igraph_error_t igraph_topological_sorting(
    const igraph_t *graph, igraph_vector_int_t *res, igraph_neimode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_is_dag(const igraph_t *graph, igraph_bool_t *res);
IGRAPH_DEPRECATED IGRAPH_EXPORT igraph_error_t igraph_transitive_closure_dag(const igraph_t *graph,
                                                igraph_t *closure);

/* -------------------------------------------------- */
/* Graph isomorphisms                                 */
/* -------------------------------------------------- */

/* Common functions */
IGRAPH_EXPORT igraph_error_t igraph_simplify_and_colorize(
    const igraph_t *graph, igraph_t *res,
    igraph_vector_int_t *vertex_color, igraph_vector_int_t *edge_color);

/* Generic interface */
IGRAPH_EXPORT igraph_error_t igraph_isomorphic(const igraph_t *graph1, const igraph_t *graph2,
                                    igraph_bool_t *iso);
IGRAPH_EXPORT igraph_error_t igraph_subisomorphic(const igraph_t *graph1, const igraph_t *graph2,
                                       igraph_bool_t *iso);

/* LAD */
IGRAPH_EXPORT igraph_error_t igraph_subisomorphic_lad(
    const igraph_t *pattern, const igraph_t *target, const igraph_vector_int_list_t *domains,
    igraph_bool_t *iso, igraph_vector_int_t *map, igraph_vector_int_list_t *maps,
    igraph_bool_t induced, igraph_integer_t time_limit
);

/* VF2 family*/
/**
 * \typedef igraph_isohandler_t
 * Callback type, called when an isomorphism was found
 *
 * See the details at the documentation of \ref
 * igraph_get_isomorphisms_vf2_callback().
 * \param map12 The mapping from the first graph to the second.
 * \param map21 The mapping from the second graph to the first, the
 *   inverse of \p map12 basically.
 * \param arg This extra argument was passed to \ref
 *   igraph_get_isomorphisms_vf2_callback() when it was called.
 * \return \c IGRAPH_SUCCESS to continue the search, \c IGRAPH_STOP to
 *   terminate the search. Any other return value is interpreted as an
 *   igraph error code, which will then abort the search and return the
 *   same error code from the caller function.
 */


typedef igraph_error_t igraph_isohandler_t(const igraph_vector_int_t *map12,
        const igraph_vector_int_t *map21, void *arg);

/**
 * \typedef igraph_isocompat_t
 * Callback type, called to check whether two vertices or edges are compatible
 *
 * VF2 (subgraph) isomorphism functions can be restricted by defining
 * relations on the vertices and/or edges of the graphs, and then checking
 * whether the vertices (edges) match according to these relations.
 *
 * </para><para>This feature is implemented by two callbacks, one for
 * vertices, one for edges. Every time igraph tries to match a vertex (edge)
 * of the first (sub)graph to a vertex of the second graph, the vertex
 * (edge) compatibility callback is called. The callback returns a
 * logical value, giving whether the two vertices match.
 *
 * </para><para>Both callback functions are of type \c igraph_isocompat_t.
 * \param graph1 The first graph.
 * \param graph2 The second graph.
 * \param g1_num The id of a vertex or edge in the first graph.
 * \param g2_num The id of a vertex or edge in the second graph.
 * \param arg Extra argument to pass to the callback functions.
 * \return Logical scalar, whether vertex (or edge) \p g1_num in \p graph1
 *    is compatible with vertex (or edge) \p g2_num in \p graph2.
 */

typedef igraph_bool_t igraph_isocompat_t(const igraph_t *graph1,
        const igraph_t *graph2,
        const igraph_integer_t g1_num,
        const igraph_integer_t g2_num,
        void *arg);

IGRAPH_EXPORT igraph_error_t igraph_isomorphic_vf2(const igraph_t *graph1, const igraph_t *graph2,
                                        const igraph_vector_int_t *vertex_color1,
                                        const igraph_vector_int_t *vertex_color2,
                                        const igraph_vector_int_t *edge_color1,
                                        const igraph_vector_int_t *edge_color2,
                                        igraph_bool_t *iso,
                                        igraph_vector_int_t *map12,
                                        igraph_vector_int_t *map21,
                                        igraph_isocompat_t *node_compat_fn,
                                        igraph_isocompat_t *edge_compat_fn,
                                        void *arg);
IGRAPH_EXPORT igraph_error_t igraph_count_isomorphisms_vf2(const igraph_t *graph1, const igraph_t *graph2,
                                                const igraph_vector_int_t *vertex_color1,
                                                const igraph_vector_int_t *vertex_color2,
                                                const igraph_vector_int_t *edge_color1,
                                                const igraph_vector_int_t *edge_color2,
                                                igraph_integer_t *count,
                                                igraph_isocompat_t *node_compat_fn,
                                                igraph_isocompat_t *edge_compat_fn,
                                                void *arg);
IGRAPH_EXPORT igraph_error_t igraph_get_isomorphisms_vf2(const igraph_t *graph1,
                                              const igraph_t *graph2,
                                              const igraph_vector_int_t *vertex_color1,
                                              const igraph_vector_int_t *vertex_color2,
                                              const igraph_vector_int_t *edge_color1,
                                              const igraph_vector_int_t *edge_color2,
                                              igraph_vector_int_list_t *maps,
                                              igraph_isocompat_t *node_compat_fn,
                                              igraph_isocompat_t *edge_compat_fn,
                                              void *arg);
IGRAPH_EXPORT igraph_error_t igraph_get_isomorphisms_vf2_callback(
    const igraph_t *graph1, const igraph_t *graph2,
    const igraph_vector_int_t *vertex_color1, const igraph_vector_int_t *vertex_color2,
    const igraph_vector_int_t *edge_color1, const igraph_vector_int_t *edge_color2,
    igraph_vector_int_t *map12, igraph_vector_int_t *map21,
    igraph_isohandler_t *isohandler_fn, igraph_isocompat_t *node_compat_fn,
    igraph_isocompat_t *edge_compat_fn, void *arg
);

/* Deprecated alias to igraph_get_isomorphisms_vf2_callback(), will be removed in 0.11 */
IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_isomorphic_function_vf2(
    const igraph_t *graph1, const igraph_t *graph2,
    const igraph_vector_int_t *vertex_color1, const igraph_vector_int_t *vertex_color2,
    const igraph_vector_int_t *edge_color1, const igraph_vector_int_t *edge_color2,
    igraph_vector_int_t *map12, igraph_vector_int_t *map21,
    igraph_isohandler_t *isohandler_fn, igraph_isocompat_t *node_compat_fn,
    igraph_isocompat_t *edge_compat_fn, void *arg
);

IGRAPH_EXPORT igraph_error_t igraph_subisomorphic_vf2(const igraph_t *graph1, const igraph_t *graph2,
                                           const igraph_vector_int_t *vertex_color1,
                                           const igraph_vector_int_t *vertex_color2,
                                           const igraph_vector_int_t *edge_color1,
                                           const igraph_vector_int_t *edge_color2,
                                           igraph_bool_t *iso,
                                           igraph_vector_int_t *map12,
                                           igraph_vector_int_t *map21,
                                           igraph_isocompat_t *node_compat_fn,
                                           igraph_isocompat_t *edge_compat_fn,
                                           void *arg);
IGRAPH_EXPORT igraph_error_t igraph_count_subisomorphisms_vf2(const igraph_t *graph1, const igraph_t *graph2,
                                                   const igraph_vector_int_t *vertex_color1,
                                                   const igraph_vector_int_t *vertex_color2,
                                                   const igraph_vector_int_t *edge_color1,
                                                   const igraph_vector_int_t *edge_color2,
                                                   igraph_integer_t *count,
                                                   igraph_isocompat_t *node_compat_fn,
                                                   igraph_isocompat_t *edge_compat_fn,
                                                   void *arg);
IGRAPH_EXPORT igraph_error_t igraph_get_subisomorphisms_vf2(const igraph_t *graph1,
                                                 const igraph_t *graph2,
                                                 const igraph_vector_int_t *vertex_color1,
                                                 const igraph_vector_int_t *vertex_color2,
                                                 const igraph_vector_int_t *edge_color1,
                                                 const igraph_vector_int_t *edge_color2,
                                                 igraph_vector_int_list_t *maps,
                                                 igraph_isocompat_t *node_compat_fn,
                                                 igraph_isocompat_t *edge_compat_fn,
                                                 void *arg);
IGRAPH_EXPORT igraph_error_t igraph_get_subisomorphisms_vf2_callback(
    const igraph_t *graph1, const igraph_t *graph2,
    const igraph_vector_int_t *vertex_color1, const igraph_vector_int_t *vertex_color2,
    const igraph_vector_int_t *edge_color1, const igraph_vector_int_t *edge_color2,
    igraph_vector_int_t *map12, igraph_vector_int_t *map21,
    igraph_isohandler_t *isohandler_fn, igraph_isocompat_t *node_compat_fn,
    igraph_isocompat_t *edge_compat_fn, void *arg
);

/* Deprecated alias to igraph_get_subisomorphisms_vf2_callback(), will be removed in 0.11 */
IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_subisomorphic_function_vf2(
    const igraph_t *graph1, const igraph_t *graph2,
    const igraph_vector_int_t *vertex_color1, const igraph_vector_int_t *vertex_color2,
    const igraph_vector_int_t *edge_color1, const igraph_vector_int_t *edge_color2,
    igraph_vector_int_t *map12, igraph_vector_int_t *map21,
    igraph_isohandler_t *isohandler_fn, igraph_isocompat_t *node_compat_fn,
    igraph_isocompat_t *edge_compat_fn, void *arg
);

/* BLISS family */
/**
 * \struct igraph_bliss_info_t
 * \brief Information about a Bliss run.
 *
 * Some secondary information found by the Bliss algorithm is stored
 * here. It is useful if you wany to study the internal working of the
 * algorithm.
 *
 * \member nof_nodes The number of nodes in the search tree.
 * \member nof_leaf_nodes The number of leaf nodes in the search tree.
 * \member nof_bad_nodes Number of bad nodes.
 * \member nof_canupdates Number of canrep updates.
 * \member nof_generators Number of generators of the automorphism group.
 * \member max_level Maximum level.
 * \member group_size The size of the automorphism group of the graph,
 *    given as a string. It should be deallocated via
 *    \ref igraph_free() if not needed any more.
 *
 * See https://users.aalto.fi/~tjunttil/bliss/
 * for details about the algorithm and these parameters.
 */
typedef struct igraph_bliss_info_t {
    unsigned long nof_nodes;
    unsigned long nof_leaf_nodes;
    unsigned long nof_bad_nodes;
    unsigned long nof_canupdates;
    unsigned long nof_generators;
    unsigned long max_level;
    char *group_size;
} igraph_bliss_info_t;

/**
 * \typedef igraph_bliss_sh_t
 * \brief Splitting heuristics for Bliss.
 *
 * \c IGRAPH_BLISS_FL provides good performance for many graphs, and is a reasonable
 * default choice. \c IGRAPH_BLISS_FSM is recommended for graphs that have some
 * combinatorial structure, and is the default of the Bliss library's command
 * line tool.
 *
 * \enumval IGRAPH_BLISS_F First non-singleton cell.
 * \enumval IGRAPH_BLISS_FL First largest non-singleton cell.
 * \enumval IGRAPH_BLISS_FS First smallest non-singleton cell.
 * \enumval IGRAPH_BLISS_FM First maximally non-trivially connected
 *      non-singleton cell.
 * \enumval IGRAPH_BLISS_FLM Largest maximally non-trivially connected
 *      non-singleton cell.
 * \enumval IGRAPH_BLISS_FSM Smallest maximally non-trivially
 *      connected non-singletion cell.
 */

typedef enum { IGRAPH_BLISS_F = 0, IGRAPH_BLISS_FL,
               IGRAPH_BLISS_FS, IGRAPH_BLISS_FM,
               IGRAPH_BLISS_FLM, IGRAPH_BLISS_FSM
             } igraph_bliss_sh_t;

IGRAPH_EXPORT igraph_error_t igraph_canonical_permutation(const igraph_t *graph, const igraph_vector_int_t *colors, igraph_vector_int_t *labeling,
                                               igraph_bliss_sh_t sh, igraph_bliss_info_t *info);
IGRAPH_EXPORT igraph_error_t igraph_isomorphic_bliss(const igraph_t *graph1, const igraph_t *graph2,
                                          const igraph_vector_int_t *colors1, const igraph_vector_int_t *colors2,
                                          igraph_bool_t *iso, igraph_vector_int_t *map12,
                                          igraph_vector_int_t *map21,
                                          igraph_bliss_sh_t sh,
                                          igraph_bliss_info_t *info1, igraph_bliss_info_t *info2);

IGRAPH_EXPORT igraph_error_t igraph_count_automorphisms(
        const igraph_t *graph, const igraph_vector_int_t *colors,
        igraph_bliss_sh_t sh, igraph_bliss_info_t *info);

IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_automorphisms(
        const igraph_t *graph, const igraph_vector_int_t *colors,
        igraph_bliss_sh_t sh, igraph_bliss_info_t *info);

IGRAPH_EXPORT igraph_error_t igraph_automorphism_group(
    const igraph_t *graph, const igraph_vector_int_t *colors,
    igraph_vector_int_list_t *generators, igraph_bliss_sh_t sh,
    igraph_bliss_info_t *info
);

/* Functions for small graphs (<= 4 vertices for directed graphs, <= 6 for undirected graphs) */
IGRAPH_EXPORT igraph_error_t igraph_isomorphic_small(const igraph_t *graph1, const igraph_t *graph2,
                                       igraph_bool_t *iso);
IGRAPH_EXPORT igraph_error_t igraph_isoclass(const igraph_t *graph, igraph_integer_t *isoclass);
IGRAPH_EXPORT igraph_error_t igraph_isoclass_subgraph(const igraph_t *graph, const igraph_vector_int_t *vids,
                                           igraph_integer_t *isoclass);
IGRAPH_EXPORT igraph_error_t igraph_isoclass_create(igraph_t *graph, igraph_integer_t size,
                                         igraph_integer_t number, igraph_bool_t directed);

IGRAPH_EXPORT igraph_error_t igraph_graph_count(igraph_integer_t n, igraph_bool_t directed, igraph_integer_t *count);

IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_isomorphic_34(
    const igraph_t *graph1, const igraph_t *graph2, igraph_bool_t *iso
);



__END_DECLS

#endif
