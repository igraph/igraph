/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2021 The igraph development team

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

#include "igraph_games.h"

#include "igraph_constructors.h"
#include "igraph_random.h"

#include "core/interruption.h"

/**
 * \function igraph_grg_game
 * \brief Generates a geometric random graph.
 *
 * A geometric random graph is created by dropping points (i.e. vertices)
 * randomly on the unit square and then connecting all those pairs
 * which are less than \c radius apart in Euclidean distance.
 *
 * </para><para>
 * Original code contributed by Keith Briggs, thanks Keith.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param nodes The number of vertices in the graph.
 * \param radius The radius within which the vertices will be connected.
 * \param torus Logical constant. If true, periodic boundary conditions
 *        will be used, i.e. the vertices are assumed to be on a torus
 *        instead of a square.
 * \param x An initialized vector or \c NULL. If not \c NULL, the points'
 *          x coordinates will be returned here.
 * \param y An initialized vector or \c NULL. If not \c NULL, the points'
 *          y coordinates will be returned here.
 * \return Error code.
 *
 * Time complexity: TODO, less than O(|V|^2+|E|).
 *
 * \example examples/simple/igraph_grg_game.c
 */
int igraph_grg_game(igraph_t *graph, igraph_integer_t nodes,
                    igraph_real_t radius, igraph_bool_t torus,
                    igraph_vector_t *x, igraph_vector_t *y) {

    long int i;
    igraph_vector_t myx, myy, *xx = &myx, *yy = &myy, edges;
    igraph_real_t r2 = radius * radius;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, nodes));

    if (x) {
        xx = x;
        IGRAPH_CHECK(igraph_vector_resize(xx, nodes));
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(xx, nodes);
    }
    if (y) {
        yy = y;
        IGRAPH_CHECK(igraph_vector_resize(yy, nodes));
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(yy, nodes);
    }

    RNG_BEGIN();

    for (i = 0; i < nodes; i++) {
        VECTOR(*xx)[i] = RNG_UNIF01();
        VECTOR(*yy)[i] = RNG_UNIF01();
    }

    RNG_END();

    igraph_vector_sort(xx);

    if (!torus) {
        for (i = 0; i < nodes; i++) {
            igraph_real_t xx1 = VECTOR(*xx)[i];
            igraph_real_t yy1 = VECTOR(*yy)[i];
            long int j = i + 1;
            igraph_real_t dx, dy;

            IGRAPH_ALLOW_INTERRUPTION();

            while ( j < nodes && (dx = VECTOR(*xx)[j] - xx1) < radius) {
                dy = VECTOR(*yy)[j] - yy1;
                if (dx * dx + dy * dy < r2) {
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
                }
                j++;
            }
        }
    } else {
        for (i = 0; i < nodes; i++) {
            igraph_real_t xx1 = VECTOR(*xx)[i];
            igraph_real_t yy1 = VECTOR(*yy)[i];
            long int j = i + 1;
            igraph_real_t dx, dy;

            IGRAPH_ALLOW_INTERRUPTION();

            while ( j < nodes && (dx = VECTOR(*xx)[j] - xx1) < radius) {
                dy = fabs(VECTOR(*yy)[j] - yy1);
                if (dx > 0.5) {
                    dx = 1 - dx;
                }
                if (dy > 0.5) {
                    dy = 1 - dy;
                }
                if (dx * dx + dy * dy < r2) {
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
                }
                j++;
            }
            if (j == nodes) {
                j = 0;
                while (j < i && (dx = 1 - xx1 + VECTOR(*xx)[j]) < radius &&
                       xx1 - VECTOR(*xx)[j] >= radius) {
                    dy = fabs(VECTOR(*yy)[j] - yy1);
                    if (dy > 0.5) {
                        dy = 1 - dy;
                    }
                    if (dx * dx + dy * dy < r2) {
                        IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                        IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
                    }
                    j++;
                }
            }
        }
    }

    if (!y) {
        igraph_vector_destroy(yy);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (!x) {
        igraph_vector_destroy(xx);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, IGRAPH_UNDIRECTED));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
