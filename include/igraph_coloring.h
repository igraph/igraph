/*
  Heuristic graph coloring algorithms.
  Copyright (C) 2017 Szabolcs Horvat <szhorvat@gmail.com>

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
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301 USA
*/

#ifndef IGRAPH_COLORING_H
#define IGRAPH_COLORING_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_error.h"

__BEGIN_DECLS

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

__END_DECLS

#endif /* IGRAPH_COLORING_H */
