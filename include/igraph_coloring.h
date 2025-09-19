/*
   igraph library.
   Copyright (C) 2017-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_COLORING_H
#define IGRAPH_COLORING_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_error.h"

IGRAPH_BEGIN_C_DECLS

/**
 * \typedef igraph_coloring_greedy_t
 * \brief Ordering heuristics for greedy graph coloring.
 *
 * Ordering heuristics for \ref igraph_vertex_coloring_greedy().
 *
 * \enumval IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS
 *    Choose the vertex with largest number of already colored neighbors.
 * \enumval IGRAPH_COLORING_GREEDY_DSATUR
 *    Choose the vertex with largest number of unique colors in its neighborhood, i.e. its
 *    "saturation degree". When multiple vertices have the same saturation degree, choose
 *    the one with the most not yet colored neighbors. Added in igraph 0.10.4. This heuristic
 *    is known as "DSatur", and was proposed in
 *    Daniel Brélaz: New methods to color the vertices of a graph,
 *    Commun. ACM 22, 4 (1979), 251–256. https://doi.org/10.1145/359094.359101
 */
typedef enum {
    IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS = 0,
    IGRAPH_COLORING_GREEDY_DSATUR = 1
} igraph_coloring_greedy_t;

IGRAPH_EXPORT igraph_error_t igraph_vertex_coloring_greedy(const igraph_t *graph, igraph_vector_int_t *colors, igraph_coloring_greedy_t heuristic);

IGRAPH_EXPORT igraph_error_t igraph_is_vertex_coloring(const igraph_t *graph, const igraph_vector_int_t *types, igraph_bool_t *res);
IGRAPH_EXPORT igraph_error_t igraph_is_bipartite_coloring(const igraph_t *graph, const igraph_vector_bool_t *types, igraph_bool_t *res, igraph_neimode_t *mode);
IGRAPH_EXPORT igraph_error_t igraph_is_edge_coloring(const igraph_t *graph, const igraph_vector_int_t *types, igraph_bool_t *res);

IGRAPH_END_C_DECLS

#endif /* IGRAPH_COLORING_H */
