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

/* COLORED_NEIGHBORS: Choose vertices based on the number of already coloured neighbours. */

static igraph_error_t igraph_i_vertex_coloring_greedy_cn(const igraph_t *graph, igraph_vector_int_t *colors) {
    igraph_integer_t i, vertex, maxdeg;
    igraph_integer_t vc = igraph_vcount(graph);
    igraph_2wheap_t cn; /* indexed heap storing number of already coloured neighbours */
    igraph_vector_int_t neigh_colors;
    igraph_adjlist_t adjlist;

    IGRAPH_CHECK(igraph_vector_int_resize(colors, vc));
    igraph_vector_int_fill(colors, 0);

    /* Nothing to do for 0 or 1 vertices.
     * Remember that colours are integers starting from 0,
     * and the 'colors' vector is already 0-initialized above.
     */
    if (vc <= 1) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    /* find maximum degree and a corresponding vertex */
    {
        igraph_vector_int_t degree;

        IGRAPH_CHECK(igraph_vector_int_init(&degree, 0));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &degree);
        IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL, 0));

        vertex = igraph_vector_int_which_max(&degree);
        maxdeg = VECTOR(degree)[vertex];

        igraph_vector_int_destroy(&degree);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_CHECK(igraph_vector_int_init(&neigh_colors, 0));
    IGRAPH_CHECK(igraph_vector_int_reserve(&neigh_colors, maxdeg));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &neigh_colors);

    IGRAPH_CHECK(igraph_2wheap_init(&cn, vc));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &cn);
    for (i = 0; i < vc; ++i)
        if (i != vertex) {
            igraph_2wheap_push_with_index(&cn, i, 0); /* should not fail since memory was already reserved */
        }

    while (1) {
        igraph_vector_int_t *neighbors = igraph_adjlist_get(&adjlist, vertex);
        igraph_integer_t neigh_count = igraph_vector_int_size(neighbors);

        /* colour current vertex */
        {
            igraph_integer_t col;

            IGRAPH_CHECK(igraph_vector_int_resize(&neigh_colors, neigh_count));
            for (i = 0; i < neigh_count; ++i) {
                VECTOR(neigh_colors)[i] = VECTOR(*colors)[ VECTOR(*neighbors)[i] ];
            }
            igraph_vector_int_sort(&neigh_colors);

            i = 0;
            col = 0;
            do {
                while (i < neigh_count && VECTOR(neigh_colors)[i] == col) {
                    i++;
                }
                col++;
            } while (i < neigh_count && VECTOR(neigh_colors)[i] == col);

            VECTOR(*colors)[vertex] = col;
        }

        /* increment number of coloured neighbours for each neighbour of vertex */
        for (i = 0; i < neigh_count; ++i) {
            igraph_integer_t idx = VECTOR(*neighbors)[i];
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
    igraph_vector_int_destroy(&neigh_colors);
    igraph_adjlist_destroy(&adjlist);
    igraph_2wheap_destroy(&cn);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/* DSATUR: Choose vertices based on the number of adjacent colours, i.e. "saturation degree" */

typedef struct {
    igraph_integer_t saturation_degree;
    igraph_integer_t edge_degree;
} dsatur_t;

static int dsatur_t_compare(const void* left, const void* right) {
    const dsatur_t *left_d  = (const dsatur_t*)left;
    const dsatur_t *right_d = (const dsatur_t*)right;
    if (left_d->saturation_degree == right_d->saturation_degree) {
        return left_d->edge_degree - right_d->edge_degree;
    }
    return left_d->saturation_degree - right_d->saturation_degree;
}

static igraph_bool_t is_color_used_by_neighbour(
    const igraph_vector_int_t *colors, igraph_integer_t color,
    const igraph_vector_int_t *neighbours
) {
    igraph_integer_t nbr_cnt = igraph_vector_int_size(neighbours);

    for (igraph_integer_t nbr_indx = 0; nbr_indx < nbr_cnt; nbr_indx++) {
        igraph_integer_t nbr = VECTOR(*neighbours)[nbr_indx];
        if (VECTOR(*colors)[nbr] == color) {
            return true;
        }
    }

    return false;
}

static void igraph_i_dsatur_update_heap(
    const igraph_adjlist_t *adjlist, igraph_gen2wheap_t *node_degrees_heap,
    igraph_vector_int_t *neighbours, igraph_vector_int_t *colors,
    igraph_integer_t color
) {
    igraph_gen2wheap_delete_max(node_degrees_heap);
    dsatur_t color_metadata;
    igraph_integer_t nbr_cnt = igraph_vector_int_size(neighbours);
    for (igraph_integer_t nbr_indx = 0; nbr_indx < nbr_cnt; nbr_indx++) {
        igraph_integer_t nbr = VECTOR(*neighbours)[nbr_indx];
        if (!igraph_gen2wheap_has_elem(node_degrees_heap, nbr)) {
            continue;
        }
        color_metadata = *((dsatur_t*)igraph_gen2wheap_get(node_degrees_heap, nbr));
        if (!is_color_used_by_neighbour(colors, color, igraph_adjlist_get(adjlist, nbr))) {
            color_metadata.saturation_degree++;
        }
        color_metadata.edge_degree--;
        igraph_gen2wheap_modify(node_degrees_heap, nbr, &color_metadata);
    }
}

static igraph_error_t igraph_i_get_first_viable_color(
    igraph_vector_int_t *used_colors
) {
    igraph_integer_t color_count = igraph_vector_int_size(used_colors);
    igraph_vector_int_sort(used_colors);
    igraph_integer_t i = 0;
    igraph_integer_t col = 0;
    while (i < color_count && VECTOR(*used_colors)[i] == col) {
        while (i < color_count && VECTOR(*used_colors)[i] == col) {
            i++;
            if (i == color_count) {
                break; /*loop second condition could read outside bounds without this*/
            }
        }
        col++;
        if (i == color_count) {
            break;
        }
    }
    return  col;
}

static igraph_error_t igraph_i_dsatur_main_loop(
    const igraph_t *graph, igraph_vector_int_t *colors,
    igraph_gen2wheap_t *node_degrees_heap, igraph_adjlist_t *adjlist,
    igraph_integer_t vc
) {
    igraph_integer_t vertices_colored = 0;
    igraph_vector_int_t used_colors;
    igraph_integer_t color;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&used_colors, 0);

    while (vertices_colored < vc) {
        igraph_integer_t node_to_color = igraph_gen2wheap_max_index(node_degrees_heap);
        igraph_vector_int_t *neighbours = igraph_adjlist_get(adjlist, node_to_color);
        igraph_integer_t neighbours_count = igraph_vector_int_size(neighbours);
        igraph_vector_int_clear(&used_colors);
        for (igraph_integer_t neighbor_index = 0; neighbor_index < neighbours_count; neighbor_index++) {
            igraph_integer_t neighbour = VECTOR(*neighbours)[neighbor_index];
            if (VECTOR(*colors)[neighbour] != -1) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&used_colors, VECTOR(*colors)[neighbour]));
            }
        }
        color = igraph_i_get_first_viable_color(&used_colors);
        igraph_i_dsatur_update_heap(adjlist, node_degrees_heap, neighbours, colors, color);
        VECTOR(*colors)[node_to_color] = color;

        vertices_colored++;
        IGRAPH_ALLOW_INTERRUPTION();
    }

    igraph_vector_int_destroy(&used_colors);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_vertex_coloring_dsatur(
    const igraph_t *graph, igraph_vector_int_t *colors
) {
    igraph_integer_t vc = igraph_vcount(graph);
    IGRAPH_CHECK(igraph_vector_int_resize(colors, vc));

    if (vc == 0) {
        return IGRAPH_SUCCESS;
    }


    if (vc == 1) {
        VECTOR(*colors)[0] = 0;
        return IGRAPH_SUCCESS;
    }

    igraph_vector_int_fill(colors, -1);   // -1 as a color means uncolored

    igraph_adjlist_t adjlist;
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    igraph_gen2wheap_t node_degrees_heap;
    IGRAPH_CHECK(igraph_gen2wheap_init(&node_degrees_heap, dsatur_t_compare, sizeof(dsatur_t), vc));
    IGRAPH_FINALLY(igraph_gen2wheap_destroy, &node_degrees_heap);

    dsatur_t igraph_color_metadata;
    for (igraph_integer_t vertex = 0; vertex < vc; vertex++) {
        igraph_color_metadata.saturation_degree = 0;
        igraph_color_metadata.edge_degree = igraph_vector_int_size(igraph_adjlist_get(&adjlist, vertex));
        IGRAPH_CHECK(igraph_gen2wheap_push_with_index(&node_degrees_heap, vertex, &igraph_color_metadata));
    }

    IGRAPH_CHECK(igraph_i_dsatur_main_loop(graph, colors, &node_degrees_heap, &adjlist, vc));

    igraph_gen2wheap_destroy(&node_degrees_heap);
    igraph_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_vertex_coloring_greedy
 * \brief Computes a vertex coloring using a greedy algorithm.
 *
 * </para><para>
 * This function assigns a "color"—represented as a non-negative integer—to
 * each vertex of the graph in such a way that neighboring vertices never have
 * the same color. The obtained coloring is not necessarily minimal.
 *
 * </para><para>
 * Vertices are colored one by one, choosing the smallest color index that
 * differs from that of already colored neighbors.
 * Colors are represented with non-negative integers 0, 1, 2, ...
 *
 * \param graph The input graph.
 * \param colors Pointer to an initialized integer vector. The vertex colors will be stored here.
 * \param heuristic The vertex ordering heuristic to use during greedy coloring. See \ref igraph_coloring_greedy_t
 *
 * \return Error code.
 *
 * \example examples/simple/igraph_coloring.c
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
