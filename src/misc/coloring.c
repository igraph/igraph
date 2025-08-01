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

#include "igraph_coloring.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"

#include "core/genheap.h"
#include "core/indheap.h"
#include "core/interruption.h"
#include "graph/internal.h" /* igraph_i_incident() */

/* COLORED_NEIGHBORS: Choose vertices based on the number of already coloured neighbours. */

static igraph_error_t igraph_i_vertex_coloring_greedy_cn(const igraph_t *graph, igraph_vector_int_t *colors) {
    igraph_integer_t i, vertex, maxdeg;
    igraph_integer_t vc = igraph_vcount(graph);
    igraph_2wheap_t cn; /* indexed heap storing number of already coloured neighbours */
    igraph_vector_int_t neighbors, nei_colors;

    IGRAPH_CHECK(igraph_vector_int_resize(colors, vc));
    igraph_vector_int_null(colors);

    /* Nothing to do for 0 or 1 vertices.
     * Remember that colours are integers starting from 0,
     * and the 'colors' vector is already 0-initialized above.
     */
    if (vc <= 1) {
        return IGRAPH_SUCCESS;
    }

    /* find maximum degree and a corresponding vertex */
    {
        igraph_vector_int_t degree;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&degree, 0);
        IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL, false));

        vertex = igraph_vector_int_which_max(&degree);
        maxdeg = VECTOR(degree)[vertex];

        igraph_vector_int_destroy(&degree);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&nei_colors, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&nei_colors, maxdeg));

    IGRAPH_VECTOR_INT_INIT_FINALLY(&neighbors, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&neighbors, maxdeg));

    /* two-way indexed heap holding number of already colored neighbors of yet-uncolored vertices */
    IGRAPH_CHECK(igraph_2wheap_init(&cn, vc));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &cn);
    for (i = 0; i < vc; ++i) {
        if (i != vertex) {
            igraph_2wheap_push_with_index(&cn, i, 0); /* should not fail since memory was already reserved */
        }
    }

    /* Within this loop, a color of 0 means "uncolored", and valid color indices start at 1.
     * At the beginning, all vertices are set as "uncolored", see the vector_int_fill() call above.
     * Colors will be decremented to start at 0 later. */
    while (true) {
        IGRAPH_CHECK(igraph_neighbors(graph, &neighbors, vertex, IGRAPH_ALL));
        igraph_integer_t nei_count = igraph_vector_int_size(&neighbors);

        /* Colour current vertex by finding the smallest available non-0 color.
         * Note that self-loops are effectively skipped as they merely prevent
         * the current vertex from being colored with the color value it presently
         * has, which is 0 (meaning uncolored). */
        {
            igraph_integer_t col;

            IGRAPH_CHECK(igraph_vector_int_resize(&nei_colors, nei_count));
            for (i = 0; i < nei_count; ++i) {
                VECTOR(nei_colors)[i] = VECTOR(*colors)[ VECTOR(neighbors)[i] ];
            }
            igraph_vector_int_sort(&nei_colors);

            i = 0;
            col = 0;
            do {
                while (i < nei_count && VECTOR(nei_colors)[i] == col) {
                    i++;
                }
                col++;
            } while (i < nei_count && VECTOR(nei_colors)[i] == col);

            VECTOR(*colors)[vertex] = col;
        }

        /* increment number of coloured neighbours for each neighbour of vertex */
        for (i = 0; i < nei_count; ++i) {
            igraph_integer_t idx = VECTOR(neighbors)[i];
            if (igraph_2wheap_has_elem(&cn, idx)) {
                igraph_2wheap_modify(&cn, idx, igraph_2wheap_get(&cn, idx) + 1);
            }
        }

        /* stop if no more vertices left to colour */
        if (igraph_2wheap_empty(&cn)) {
            break;
        }

        igraph_2wheap_delete_max_index(&cn, &vertex);

        IGRAPH_ALLOW_INTERRUPTION();
    }

    /* subtract 1 from each colour value, so that colours start at 0 */
    igraph_vector_int_add_constant(colors, -1);

    /* free data structures */
    igraph_vector_int_destroy(&neighbors);
    igraph_vector_int_destroy(&nei_colors);
    igraph_2wheap_destroy(&cn);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/* DSATUR: Choose vertices based on the number of adjacent colours, i.e. "saturation degree" */

typedef struct {
    igraph_integer_t saturation_degree; /* number of colors used by neighbors */
    igraph_integer_t edge_degree; /* degree in the subgraph induced by uncolored vertices */
} dsatur_t;

static int dsatur_t_compare(const void *left, const void *right) {
    const dsatur_t *left_d  = left;
    const dsatur_t *right_d = right;
    if (left_d->saturation_degree == right_d->saturation_degree) {
        if (left_d->edge_degree == right_d->edge_degree) {
            return 0;
        } else if (left_d->edge_degree > right_d->edge_degree) {
            return 1;
        } else {
            return -1;
        }
    }
    return left_d->saturation_degree > right_d->saturation_degree ? 1 : -1;
}

static igraph_bool_t dsatur_is_color_used_by_neighbour(
    const igraph_vector_int_t *colors, igraph_integer_t color,
    const igraph_vector_int_t *neighbors
) {
    igraph_integer_t nei_count = igraph_vector_int_size(neighbors);

    for (igraph_integer_t i=0; i < nei_count; i++) {
        igraph_integer_t nei = VECTOR(*neighbors)[i];
        if (VECTOR(*colors)[nei] == color) {
            return true;
        }
    }

    return false;
}

static void dsatur_update_heap(
    const igraph_adjlist_t *adjlist, igraph_gen2wheap_t *node_degrees_heap,
    const igraph_vector_int_t *neighbors, const igraph_vector_int_t *colors,
    igraph_integer_t color
) {
    igraph_gen2wheap_delete_max(node_degrees_heap);
    igraph_integer_t nei_count = igraph_vector_int_size(neighbors);
    for (igraph_integer_t i=0; i < nei_count; i++) {
        igraph_integer_t nei = VECTOR(*neighbors)[i];
        if (!igraph_gen2wheap_has_elem(node_degrees_heap, nei)) {
            continue;
        }
        dsatur_t deg_data = *((dsatur_t*) igraph_gen2wheap_get(node_degrees_heap, nei));
        if (!dsatur_is_color_used_by_neighbour(colors, color, igraph_adjlist_get(adjlist, nei))) {
            deg_data.saturation_degree++;
        }
        deg_data.edge_degree--;
        igraph_gen2wheap_modify(node_degrees_heap, nei, &deg_data);
    }
}

static igraph_integer_t dsatur_get_first_viable_color(const igraph_vector_int_t *used_colors_sorted) {
    igraph_integer_t color_count = igraph_vector_int_size(used_colors_sorted);
    igraph_integer_t i = 0;
    igraph_integer_t col = 0;
    while (i < color_count && VECTOR(*used_colors_sorted)[i] == col) {
        while (i < color_count && VECTOR(*used_colors_sorted)[i] == col) {
            i++;
        }
        col++;
    }
    return  col;
}

static igraph_error_t igraph_i_vertex_coloring_dsatur(
    const igraph_t *graph, igraph_vector_int_t *colors
) {
    igraph_integer_t vcount = igraph_vcount(graph);
    IGRAPH_CHECK(igraph_vector_int_resize(colors, vcount));

    if (vcount == 0) {
        return IGRAPH_SUCCESS;
    }

    if (vcount == 1) {
        VECTOR(*colors)[0] = 0;
        return IGRAPH_SUCCESS;
    }

    igraph_vector_int_fill(colors, -1); /* -1 as a color means uncolored */

    /* Multi-edges and self-loops are removed from the adjacency list in order to ensure the correct
     * updating of a vertex's neighbors' saturation degrees when that vertex is colored. */
    igraph_adjlist_t adjlist;
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    igraph_gen2wheap_t node_degrees_heap;
    IGRAPH_CHECK(igraph_gen2wheap_init(&node_degrees_heap, dsatur_t_compare, sizeof(dsatur_t), vcount));
    IGRAPH_FINALLY(igraph_gen2wheap_destroy, &node_degrees_heap);

    for (igraph_integer_t vertex = 0; vertex < vcount; vertex++) {
        dsatur_t dsatur;
        dsatur.saturation_degree = 0;
        dsatur.edge_degree = igraph_vector_int_size(igraph_adjlist_get(&adjlist, vertex));
        IGRAPH_CHECK(igraph_gen2wheap_push_with_index(&node_degrees_heap, vertex, &dsatur));
    }

    igraph_vector_int_t used_colors_sorted;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&used_colors_sorted, 0);

    while (! igraph_gen2wheap_empty(&node_degrees_heap)) {
        igraph_integer_t node_to_color = igraph_gen2wheap_max_index(&node_degrees_heap);
        igraph_vector_int_t *neighbors = igraph_adjlist_get(&adjlist, node_to_color);
        igraph_integer_t nei_count = igraph_vector_int_size(neighbors);
        igraph_vector_int_clear(&used_colors_sorted);
        for (igraph_integer_t i=0; i < nei_count; i++) {
            igraph_integer_t nei = VECTOR(*neighbors)[i];
            if (VECTOR(*colors)[nei] != -1) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&used_colors_sorted, VECTOR(*colors)[nei]));
            }
        }
        igraph_vector_int_sort(&used_colors_sorted);
        igraph_integer_t color = dsatur_get_first_viable_color(&used_colors_sorted);
        dsatur_update_heap(&adjlist, &node_degrees_heap, neighbors, colors, color);
        VECTOR(*colors)[node_to_color] = color;

        IGRAPH_ALLOW_INTERRUPTION();
    }

    igraph_vector_int_destroy(&used_colors_sorted);
    igraph_gen2wheap_destroy(&node_degrees_heap);
    igraph_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_vertex_coloring_greedy
 * \brief Computes a vertex coloring using a greedy algorithm.
 *
 * This function assigns a "color"—represented as a non-negative integer—to
 * each vertex of the graph in such a way that neighboring vertices never have
 * the same color. The obtained coloring is not necessarily minimal.
 *
 * </para><para>
 * Vertices are colored greedily, one by one, always choosing the smallest color
 * index that differs from that of already colored neighbors. Vertices are picked
 * in an order determined by the speified heuristic.
 * Colors are represented by non-negative integers 0, 1, 2, ...
 *
 * \param graph The input graph.
 * \param colors Pointer to an initialized integer vector. The vertex colors will be stored here.
 * \param heuristic The vertex ordering heuristic to use during greedy coloring.
 *    See \ref igraph_coloring_greedy_t for more information.
 *
 * \return Error code.
 *
 * \sa igraph_is_vertex_coloring() to check if a coloring is valid, i.e. if all
 * edges connect vertices of different colors.
 *
 * \example examples/simple/coloring.c
 */
igraph_error_t igraph_vertex_coloring_greedy(const igraph_t *graph, igraph_vector_int_t *colors, igraph_coloring_greedy_t heuristic) {
    switch (heuristic) {
    case IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS:
        return igraph_i_vertex_coloring_greedy_cn(graph, colors);
    case IGRAPH_COLORING_GREEDY_DSATUR:
        return igraph_i_vertex_coloring_dsatur(graph, colors);
    default:
        IGRAPH_ERROR("Invalid heuristic for greedy vertex coloring.", IGRAPH_EINVAL);
    }
}

/**
 * \function igraph_is_vertex_coloring
 * \brief Checks whether a vertex coloring is valid.
 *
 * \experimental
 *
 * This function checks whether the given vertex type/color assignment is a valid
 * vertex coloring, i.e., no two adjacent vertices have the same color.
 * Self-loops are ignored.
 *
 * \param graph The input graph.
 * \param types The vertex types/colors as an integer vector.
 * \param res Pointer to a boolean, the result is stored here.
 * \return Error code.
 *
 * Time complexity: O(|E|), linear in the number of edges.
 */
igraph_error_t igraph_is_vertex_coloring(
        const igraph_t *graph,
        const igraph_vector_int_t *types,
        igraph_bool_t *res) {

    const igraph_integer_t vcount = igraph_vcount(graph);
    const igraph_integer_t ecount = igraph_ecount(graph);
    int iter = 0;

    if (igraph_vector_int_size(types) != vcount) {
        IGRAPH_ERROR("Invalid vertex type vector length.", IGRAPH_EINVAL);
    }

    *res = true;

    for (igraph_integer_t e = 0; e < ecount; e++) {
        igraph_integer_t from = IGRAPH_FROM(graph, e);
        igraph_integer_t to = IGRAPH_TO(graph, e);

        /* Skip self-loops */
        if (from == to) {
            continue;
        }

        if (VECTOR(*types)[from] == VECTOR(*types)[to]) {
            *res = false;
            break;
        }

        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 10);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_bipartite_coloring
 * \brief Checks whether a bipartite vertex coloring is valid.
 *
 * \experimental
 *
 * This function checks whether the given vertex type assignment is a valid
 * bipartite coloring, i.e., no two adjacent vertices have the same type.
 * Additionally, for directed graphs, it determines the mode of edge directions.
 * Self-loops are ignored.
 *
 * \param graph The input graph.
 * \param types The vertex types as a boolean vector.
 * \param res Pointer to a boolean, the result is stored here.
 * \param mode Pointer to store the edge direction mode. Can be \c NULL if not needed.
 *   If all edges go from false to true vertices, \c IGRAPH_OUT is returned.
 *   If all edges go from true to false vertices, \c IGRAPH_IN is returned.
 *   If edges go in both directions or graph is undirected, \c IGRAPH_ALL is returned.
 * \return Error code.
 *
 * Time complexity: O(|E|), linear in the number of edges.
 *
 * \sa igraph_is_bipartite() to determine whether a graph is bipartite,
 * i.e. 2-colorable, and find such a coloring.
 */
igraph_error_t igraph_is_bipartite_coloring(
        const igraph_t *graph,
        const igraph_vector_bool_t *types,
        igraph_bool_t *res,
        igraph_neimode_t *mode) {

    const igraph_integer_t vcount = igraph_vcount(graph);
    const igraph_integer_t ecount = igraph_ecount(graph);
    int iter = 0;

    if (igraph_vector_bool_size(types) != vcount) {
        IGRAPH_ERROR("Invalid vertex type vector length.", IGRAPH_EINVAL);
    }

    *res = true;
    if (mode) {
        *mode = IGRAPH_ALL;
    }

    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_bool_t has_false_to_true = false;
    igraph_bool_t has_true_to_false = false;

    for (igraph_integer_t e = 0; e < ecount; e++) {
        igraph_integer_t from = IGRAPH_FROM(graph, e);
        igraph_integer_t to = IGRAPH_TO(graph, e);

        /* Skip self-loops */
        if (from == to) {
            continue;
        }

        igraph_bool_t from_type = VECTOR(*types)[from];
        igraph_bool_t to_type = VECTOR(*types)[to];

        /* Check if adjacent vertices have the same type */
        if (from_type == to_type) {
            *res = false;
            break;
        }

        /* Track edge directions for directed graphs */
        if (directed && mode) {
            if (!from_type && to_type) {
                has_false_to_true = true;
            } else if (from_type && !to_type) {
                has_true_to_false = true;
            }
        }

        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 10);
    }

    /* Determine the mode for directed graphs */
    if (*res && directed && mode) {
        if (has_false_to_true && !has_true_to_false) {
            *mode = IGRAPH_OUT;
        } else if (!has_false_to_true && has_true_to_false) {
            *mode = IGRAPH_IN;
        } else {
            *mode = IGRAPH_ALL;
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_edge_coloring
 * \brief Checks whether an edge coloring is valid.
 *
 * \experimental
 *
 * This function checks whether the given edge color assignment is a valid
 * edge coloring, i.e., no two adjacent edges have the same color.
 *
 * Note that this function does not consider self-edges (loops) as being
 * adjacent to themselves, so graphs with self-loops may still be considered
 * to have a valid edge coloring.
 *
 * \param graph The input graph.
 * \param types The edge colors as an integer vector.
 * \param res Pointer to a boolean, the result is stored here.
 * \return Error code.
 *
 * Time complexity: O(|V|*d*log(d)), where d is the maximum degree.
 */
igraph_error_t igraph_is_edge_coloring(
        const igraph_t *graph,
        const igraph_vector_int_t *types,
        igraph_bool_t *res) {

    const igraph_integer_t vcount = igraph_vcount(graph);
    const igraph_integer_t ecount = igraph_ecount(graph);
    igraph_vector_int_t edges, edge_colors;
    int iter = 0;

    if (igraph_vector_int_size(types) != ecount) {
        IGRAPH_ERROR("Invalid edge type vector length.", IGRAPH_EINVAL);
    }

    /* Be optimistic */
    *res = true;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edge_colors, 0);

    /* For each vertex, check that all incident edges have different colors */
    for (igraph_integer_t v = 0; v < vcount; v++) {
        IGRAPH_CHECK(igraph_i_incident(graph, &edges, v, IGRAPH_ALL, IGRAPH_LOOPS_ONCE));

        /* Get sorted edge color list */
        IGRAPH_CHECK(igraph_vector_int_index(types, &edge_colors, &edges));
        igraph_vector_int_sort(&edge_colors);

        /* Look for consecutive duplicates in edge color list */
        igraph_integer_t edge_color_count = igraph_vector_int_size(&edge_colors);
        for (igraph_integer_t i = 0; i < edge_color_count - 1; i++) {
            if (VECTOR(edge_colors)[i] == VECTOR(edge_colors)[i + 1]) {
                *res = false;
                goto done;
            }
        }

        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 7);
    }

done:
    igraph_vector_int_destroy(&edge_colors);
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
