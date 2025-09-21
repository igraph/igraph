/*
   igraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_interface.h"
#include "igraph_nongraph.h"
#include "igraph_paths.h"

#include "core/interruption.h"

/**
 * \ingroup nongraph
 * \function igraph_running_mean
 * \brief Calculates the running mean of a vector.
 *
 * </para><para>
 * The running mean is defined by the mean of the
 * previous \p binwidth values.
 * \param data The vector containing the data.
 * \param res The vector containing the result. This should be
 *        initialized before calling this function and will be
 *        resized.
 * \param binwidth Integer giving the width of the bin for the running
 *        mean calculation.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the length of
 * the data vector.
 */

igraph_error_t igraph_running_mean(const igraph_vector_t *data, igraph_vector_t *res,
                        igraph_int_t binwidth) {

    double sum = 0;
    igraph_int_t i;

    /* Check */
    if (igraph_vector_size(data) < binwidth) {
        IGRAPH_ERRORF("Data vector length (%" IGRAPH_PRId ") smaller than bin width (%" IGRAPH_PRId ").", IGRAPH_EINVAL, igraph_vector_size(data), binwidth);
    }
    if (binwidth < 1) {
        IGRAPH_ERRORF("Bin width for running mean should be at least 1, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, binwidth);
    }

    /* Memory for result */

    IGRAPH_CHECK(igraph_vector_resize(res, (igraph_vector_size(data) - binwidth + 1)));

    /* Initial bin */
    for (i = 0; i < binwidth; i++) {
        sum += VECTOR(*data)[i];
    }

    VECTOR(*res)[0] = sum / binwidth;

    for (i = 1; i < igraph_vector_size(data) - binwidth + 1; i++) {
        IGRAPH_ALLOW_INTERRUPTION();
        sum -= VECTOR(*data)[i - 1];
        sum += VECTOR(*data)[ (i + binwidth - 1)];
        VECTOR(*res)[i] = sum / binwidth;
    }

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_expand_path_to_pairs
 * \brief Helper function to convert a sequence of vertex IDs describing a path into a "pairs" vector.
 *
 * </para><para>
 * This function is useful when you have a sequence of vertex IDs in a graph and
 * you would like to retrieve the IDs of the edges between them. The function
 * duplicates all but the first and the last elements in the vector, effectively
 * converting the path into a vector of vertex IDs that can be passed to
 * \ref igraph_get_eids().
 *
 * \param  path  the input vector. It will be modified in-place and it will be
 *         resized as needed. When the vector contains less than two vertex IDs,
 *         it will be cleared.
 * \return Error code: \c IGRAPH_ENOMEM if there is not enough memory to expand
 *         the vector.
 */
igraph_error_t igraph_expand_path_to_pairs(igraph_vector_int_t* path) {
    igraph_int_t no_of_vertices = igraph_vector_int_size(path);
    igraph_int_t i, j, no_of_items = (no_of_vertices - 1) * 2;

    if (no_of_vertices <= 1) {
        igraph_vector_int_clear(path);
    } else {
        IGRAPH_CHECK(igraph_vector_int_resize(path, no_of_items));

        i = no_of_vertices - 1;
        j = no_of_items - 1;
        VECTOR(*path)[j] = VECTOR(*path)[i];
        while (i > 1) {
            i--; j--;
            VECTOR(*path)[j] = VECTOR(*path)[i];
            j--;
            VECTOR(*path)[j] = VECTOR(*path)[i];
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_vertex_path_from_edge_path
 * \brief Converts a walk of edge IDs to the traversed vertex IDs.
 *
 * This function is useful when you have a sequence of edge IDs representing a
 * continuous walk in a graph and you would like to obtain the vertex IDs that
 * the walk traverses. The function is used implicitly by several shortest path
 * related functions to convert a path of edge IDs to the corresponding
 * representation that describes the path in terms of vertex IDs instead.
 *
 * </para><para>
 * The result will always contain one more vertex than the number of provided
 * edges. If no edges are given, the output will contain only the start vertex.
 *
 * </para><para>
 * The walk is allowed to traverse the same vertex more than once. It is
 * suitable for use on paths, cycles, or arbitrary walks.
 *
 * \param  graph The graph that the edge IDs refer to.
 * \param  start The start vertex of the path. If negative, it is determined
 *         automatically. This is only possible if the walk contains at least
 *         one edge. If only one edge is present in an undirected walk,
 *         one of its endpoints will be selected arbitrarily.
 * \param  edge_path The sequence of edge IDs that describe the path.
 * \param  vertex_path The sequence of vertex IDs traversed will be returned here.
 * \param  mode A constant specifying how edge directions are
 *         considered in directed graphs. \c IGRAPH_OUT follows edge
 *         directions, \c IGRAPH_IN follows the opposite directions,
 *         and \c IGRAPH_ALL ignores edge directions. This argument is
 *         ignored for undirected graphs.
 * \return Error code: \c IGRAPH_ENOMEM if there is not enough memory,
 *         \c IGRAPH_EINVVID if the start vertex is invalid,
 *         \c IGRAPH_EINVAL if the edge walk does not start at the given vertex
 *         or if there is at least one edge whose start vertex does not match
 *         the end vertex of the previous edge.
 *
 * Time complexity: O(n) where n is the length of the walk.
 */
igraph_error_t igraph_vertex_path_from_edge_path(
   const igraph_t *graph, igraph_int_t start,
   const igraph_vector_int_t *edge_path, igraph_vector_int_t *vertex_path,
   igraph_neimode_t mode
) {
    const igraph_int_t no_of_edges = igraph_vector_int_size(edge_path);

    igraph_vector_int_clear(vertex_path);
    IGRAPH_CHECK(igraph_vector_int_reserve(vertex_path, no_of_edges + 1));

    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    /* Detect start vertex automatically if necessary. */
    if (start < 0) {
        if (no_of_edges == 0) {
            IGRAPH_ERROR("The path must contain at least one edge in order to "
                         "determine its starting vertex automatically.", IGRAPH_EINVAL);
        }
        const igraph_int_t edge = VECTOR(*edge_path)[0];
        switch (mode) {
        case IGRAPH_OUT:
            start = IGRAPH_FROM(graph, edge);
            break;
        case IGRAPH_IN:
            start = IGRAPH_TO(graph, edge);
            break;
        case IGRAPH_ALL:
            if (no_of_edges > 1) {
                const igraph_int_t from = IGRAPH_FROM(graph, edge);
                const igraph_int_t to = IGRAPH_TO(graph, edge);
                const igraph_int_t next_edge = VECTOR(*edge_path)[1];
                if (to == IGRAPH_FROM(graph, next_edge) || to == IGRAPH_TO(graph, next_edge)) {
                    start = from;
                } else {
                    start = to;
                    /* If the walk is not continuous, this will be detected in the next stage. */
                }
            } else {
                /* There is a single undirected edge.
                 * Choose an arbitrary endpoint the start vertex. */
                start = IGRAPH_FROM(graph, edge);
            }
        }
    }

    if (start >= igraph_vcount(graph)) {
        IGRAPH_ERROR("Invalid start vertex.", IGRAPH_EINVVID);
    }

    for (igraph_int_t i = 0; i < no_of_edges; i++) {
        const igraph_int_t edge = VECTOR(*edge_path)[i];
        const igraph_int_t from = IGRAPH_FROM(graph, edge);
        const igraph_int_t to = IGRAPH_TO(graph, edge);
        igraph_bool_t next_edge_ok;
        igraph_int_t next_start;

        igraph_vector_int_push_back(vertex_path, start);  /* reserved */

        switch (mode) {
        case IGRAPH_OUT:
            next_edge_ok = from == start;
            next_start = to;
            break;

        case IGRAPH_IN:
            next_edge_ok = to == start;
            next_start = from;
            break;

        case IGRAPH_ALL:
            if (from == start) {
                next_edge_ok = true;
                next_start = to;
            } else if (to == start) {
                next_edge_ok = true;
                next_start = from;
            } else {
                next_edge_ok = false;
            }
            break;

        default:
            IGRAPH_ERROR("Invalid neighborhood mode.", IGRAPH_EINVMODE);
        }

        if (!next_edge_ok) {
            IGRAPH_ERROR("Edge IDs do not form a continuous path.", IGRAPH_EINVAL);
        }

        start = next_start;
    }

    igraph_vector_int_push_back(vertex_path, start);  /* reserved */

    return IGRAPH_SUCCESS;
}
