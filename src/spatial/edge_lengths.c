/*
   igraph library.
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

#include "igraph_spatial.h"

#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_matrix.h"
#include "igraph_vector.h"

#include <math.h> /* sqrt, fabs */

/**
 * \function igraph_spatial_edge_lengths
 * \brief Edge lengths based on spatial vertex coordinates.
 *
 * \experimental
 *
 * The length of each edge is computed based on spatial coordinates. The length
 * can be employed by several igraph functions, such as \ref igraph_voronoi(),
 * \ref igraph_betweenness(), \ref igraph_closeness() and others.
 *
 * \param graph The graph whose edge lengths are to be computed.
 * \param lengths An initialized vector. Length will be stored here, in the
 *    order of edge IDs. It will be resized as needed.
 * \param points A matrix of vertex coordinates. Each row contains the
 *    coordinates of the corresponding vertex, in the order of vertex IDs.
 *    Arbitrary dimensional point sets are supported.
 * \param metric The distance metric to use. See \ref igraph_metric_t for
 *    valid values.
 * \return Error code.
 *
 * \sa \ref igraph_nearest_neighbor_graph() computes a k nearest neighbor graph
 * and \ref igraph_delaunay_graph() computes a Delaunay graph based on a set of
 * spatial points.
 *
 * Time complexity: O(|E| d) where |E| is the number of edges and d is the
 * dimensionality of the point set.
 */
igraph_error_t igraph_spatial_edge_lengths(
        const igraph_t *graph,
        igraph_vector_t *lengths,
        const igraph_matrix_t *points,
        igraph_metric_t metric) {

    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_int_t ecount = igraph_ecount(graph);
    const igraph_int_t dim = igraph_matrix_ncol(points);

    /* Validate input.
     *
     * We opt not to check for Inf/NaN, as these can safely propagate to the result. */

    if (igraph_matrix_nrow(points) != vcount) {
        IGRAPH_ERROR("Number of vertex coordinates must match the vertex count.", IGRAPH_EINVAL);
    }

    /* Special case: we allow zero columns when there are zero rows, i.e. no points.
     * Some languages cannot represent size-zero matrices and length-zero point
     * lists may translate as a 0-by-0 matrix to igraph. */
    if (dim == 0 && vcount > 0) {
        IGRAPH_ERROR("Vertex coordinates must not be zero-dimensional.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_resize(lengths, ecount));

    /* Validate distance metric. The switch statement is helpful because
     * compilers tend to show warnings when some enum values are missing. */
    switch (metric) {
        case IGRAPH_METRIC_EUCLIDEAN:
        case IGRAPH_METRIC_MANHATTAN:
            break;
        default:
            IGRAPH_ERROR("Invalid distance metric.", IGRAPH_EINVAL);
    }

    /* Compute edge lengths. */

    for (igraph_int_t eid=0; eid < ecount; eid++) {
        const igraph_int_t from = IGRAPH_FROM(graph, eid);
        const igraph_int_t to = IGRAPH_TO(graph, eid);
        igraph_real_t length = 0.0;

        for (igraph_int_t i=0; i < dim; i++) {
            igraph_real_t diff = MATRIX(*points, from, i) - MATRIX(*points, to, i);

            switch (metric) {
                case IGRAPH_METRIC_EUCLIDEAN:
                    length += diff*diff; break;
                case IGRAPH_METRIC_MANHATTAN:
                    length += fabs(diff); break;
            }
        }

        if (metric == IGRAPH_METRIC_EUCLIDEAN) {
            length = sqrt(length);
        }

        VECTOR(*lengths)[eid] = length;
    }

    return IGRAPH_SUCCESS;
}
