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

/* Heuristic: The next vertex to color will be the one with the most already-colored neighbors. */
static igraph_error_t igraph_i_vertex_coloring_greedy_cn(const igraph_t *graph, igraph_vector_int_t *colors) {
    igraph_integer_t i, vertex, maxdeg;
    igraph_integer_t vc = igraph_vcount(graph);
    igraph_2wheap_t cn; /* indexed heap storing number of already coloured neighbours */
    igraph_vector_int_t neighbors, nei_colors;

    IGRAPH_CHECK(igraph_vector_int_resize(colors, vc));
    igraph_vector_int_fill(colors, 0);

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
        IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL, 0));

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

        /* Colour current vertex by finding smallest available non-0 color.
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
    default:
        IGRAPH_ERROR("Invalid heuristic for greedy vertex coloring.", IGRAPH_EINVAL);
    }
}
