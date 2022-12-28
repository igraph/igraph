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

#include "igraph_interface.h"
#include "igraph_adjlist.h"

#include "core/indheap.h"
#include "core/interruption.h"

/* Greedy coloring heuristics */

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

/* DSatur coloring heuristics */

static igraph_integer_t igraph_i_dsatur_select_node(
    const igraph_t *graph, const igraph_vector_int_t *colors,
    const igraph_vector_int_t *saturation_degree, const igraph_vector_int_t *edge_degree,
    igraph_integer_t vc
) {
    igraph_integer_t most_saturated_node = -1, max_saturation = -1, max_degree = -1;

    for (igraph_integer_t node = 0; node < vc; node++) {
        //finding an uncolored node with max (saturation degree ,  degree)
        if (VECTOR(*colors)[node] != -1){
            continue;
        }
        igraph_integer_t saturation = VECTOR(*saturation_degree)[node];
        igraph_integer_t degree = VECTOR(*edge_degree)[node] ;
        if (( saturation < max_saturation) || (saturation == max_saturation && degree <= max_degree)) {
            continue;
        }
        max_saturation = saturation;
        max_degree = degree;
        most_saturated_node = node;
    }

    return most_saturated_node;
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

static void igraph_i_dsatur_update_saturation_degree_and_edge_degree(
    const igraph_adjlist_t *adjlist, igraph_vector_int_t *saturation_degree,
    igraph_vector_int_t *edge_degree, igraph_vector_int_t *neighbours, 
    igraph_vector_int_t *colors, igraph_integer_t color
) {
    igraph_integer_t nbr_cnt = igraph_vector_int_size(neighbours);
    for (igraph_integer_t nbr_indx = 0; nbr_indx < nbr_cnt; nbr_indx++){
        igraph_integer_t nbr = VECTOR(*neighbours)[nbr_indx];
        if (!is_color_used_by_neighbour(colors, color, igraph_adjlist_get(adjlist, nbr))) {
            VECTOR(*saturation_degree)[nbr]++;
        }
        VECTOR(*edge_degree)[nbr]--;
    }
}

static igraph_integer_t igraph_i_vector_bool_first_true(igraph_vector_bool_t *bool_vector){
    for(igraph_integer_t index = 0 ; index<igraph_vector_bool_size(bool_vector) ; index++){
        if(VECTOR(*bool_vector)[index]){
            return index;
        }
    } 
    return -1;
}
static igraph_error_t igraph_i_dsatur_main_loop(
    const igraph_t *graph, igraph_vector_int_t *colors,
    igraph_vector_int_t *saturation_degree, igraph_vector_int_t *edge_degree,
    igraph_adjlist_t *adjlist, igraph_integer_t vc
) {
    igraph_integer_t vertices_colored = 0;
    igraph_vector_bool_t viable_colors;
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&viable_colors, vc);
    while (vertices_colored < vc) {
        igraph_integer_t node_to_color = igraph_i_dsatur_select_node(graph, colors, saturation_degree, edge_degree, vc);
        igraph_vector_int_t *neighbours = igraph_adjlist_get(adjlist, node_to_color);
        igraph_integer_t neighbours_count = igraph_vector_int_size(neighbours);
        igraph_vector_bool_fill(&viable_colors, true);
        for(igraph_integer_t neighbor_index = 0 ; neighbor_index<neighbours_count ; neighbor_index++){
            igraph_integer_t neighbour = VECTOR(*neighbours)[neighbor_index]; 
            if(VECTOR(*colors)[neighbour] != -1){
                VECTOR(viable_colors)[VECTOR(*colors)[neighbour]] = false;
            }
        }
        igraph_integer_t color = igraph_i_vector_bool_first_true(&viable_colors);
        igraph_i_dsatur_update_saturation_degree_and_edge_degree(adjlist, saturation_degree, edge_degree, neighbours, colors, color);
        VECTOR(*colors)[node_to_color] = color;

        vertices_colored++;
    }
    igraph_vector_bool_destroy(&viable_colors);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
static igraph_error_t igraph_i_vertex_coloring_dsatur(
    const igraph_t *graph, igraph_vector_int_t *colors
) {
    igraph_vector_int_fill(colors, -1);   // -1 as a color means uncolored
    igraph_integer_t vc = igraph_vcount(graph);

    igraph_adjlist_t adjlist;
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    igraph_vector_int_t saturation_degree;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&saturation_degree, vc );


    igraph_vector_int_t edge_degree;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edge_degree, vc );
    IGRAPH_CHECK(igraph_degree(graph, &edge_degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS ) );

    IGRAPH_CHECK(igraph_i_dsatur_main_loop(graph, colors, &saturation_degree, &edge_degree, &adjlist, vc));

    igraph_vector_int_destroy(&saturation_degree);
    igraph_vector_int_destroy(&edge_degree);
    igraph_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(3);

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
