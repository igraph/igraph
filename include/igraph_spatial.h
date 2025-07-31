/*
   IGraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_SPATIAL_H
#define IGRAPH_SPATIAL_H

#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_error.h"
#include "igraph_matrix.h"

IGRAPH_BEGIN_C_DECLS

/**
 * \typedef igraph_metric_t
 * \brief Metric functions for use with spatial computation.
 *
 * \enumval IGRAPH_METRIC_EUCLIDEAN The Euclidean distance, i.e. L2 metric.
 * \enumval IGRAPH_METRIC_MANHATTAN The Manhattan distance, i.e. L1 metric.
 */
typedef enum {
    IGRAPH_METRIC_EUCLIDEAN = 0,
    IGRAPH_METRIC_L2 = IGRAPH_METRIC_EUCLIDEAN,
    IGRAPH_METRIC_MANHATTAN = 1,
    IGRAPH_METRIC_L1 = IGRAPH_METRIC_MANHATTAN
} igraph_metric_t;

IGRAPH_EXPORT igraph_error_t igraph_nearest_neighbor_graph(
    igraph_t *graph,
    const igraph_matrix_t *points,
    igraph_metric_t metric,
    igraph_integer_t neighbors,
    igraph_real_t cutoff,
    igraph_bool_t directed
);

IGRAPH_END_C_DECLS

#endif
