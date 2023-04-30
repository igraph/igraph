/* -*- mode: C -*-  */
/*
   IGraph R package.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_types.h"

#include "core/grid.h"

#include <math.h>

/* internal function */

static void igraph_i_2dgrid_which(
    igraph_2dgrid_t *grid, igraph_real_t xc, igraph_real_t yc,
    igraph_integer_t *x, igraph_integer_t *y
) {
    if (xc <= grid->minx) {
        *x = 0;
    } else if (xc >= grid->maxx) {
        *x = grid->stepsx - 1;
    } else {
        *x = floor((xc - (grid->minx)) / (grid->deltax));
    }

    if (yc <= grid->miny) {
        *y = 0;
    } else if (yc >= grid->maxy) {
        *y = grid->stepsy - 1;
    } else {
        *y = floor((yc - (grid->miny)) / (grid->deltay));
    }
}

igraph_error_t igraph_2dgrid_init(igraph_2dgrid_t *grid, igraph_matrix_t *coords,
                       igraph_real_t minx, igraph_real_t maxx, igraph_real_t deltax,
                       igraph_real_t miny, igraph_real_t maxy, igraph_real_t deltay) {
    igraph_integer_t no_of_points;

    IGRAPH_ASSERT(minx <= maxx);
    IGRAPH_ASSERT(miny <= maxy);
    IGRAPH_ASSERT(deltax > 0 && deltay > 0);
    IGRAPH_ASSERT(isfinite(minx) && isfinite(maxx) && isfinite(miny) && isfinite(maxy));
    IGRAPH_ASSERT(isfinite(deltax) && isfinite(deltay));

    grid->coords = coords;
    grid->minx = minx;
    grid->maxx = maxx;
    grid->deltax = deltax;
    grid->miny = miny;
    grid->maxy = maxy;
    grid->deltay = deltay;

    grid->stepsx = ceil((maxx - minx) / deltax);
    grid->stepsy = ceil((maxy - miny) / deltay);

    no_of_points = igraph_matrix_nrow(coords);

    IGRAPH_CHECK(igraph_matrix_int_init(&grid->startidx, grid->stepsx, grid->stepsy));
    IGRAPH_FINALLY(igraph_matrix_int_destroy, &grid->startidx);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&grid->next, no_of_points);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&grid->prev, no_of_points);

    igraph_vector_int_fill(&grid->prev, 0);
    igraph_vector_int_fill(&grid->next, 0);

    grid->massx = 0;
    grid->massy = 0;
    grid->vertices = 0;

    IGRAPH_FINALLY_CLEAN(3);
    return IGRAPH_SUCCESS;
}

void igraph_2dgrid_destroy(igraph_2dgrid_t *grid) {
    igraph_matrix_int_destroy(&grid->startidx);
    igraph_vector_int_destroy(&grid->next);
    igraph_vector_int_destroy(&grid->prev);
}

void igraph_2dgrid_add(igraph_2dgrid_t *grid, igraph_integer_t elem,
                       igraph_real_t xc, igraph_real_t yc) {
    igraph_integer_t x, y;
    igraph_integer_t first;

    MATRIX(*grid->coords, elem, 0) = xc;
    MATRIX(*grid->coords, elem, 1) = yc;

    /* add to cell */
    igraph_i_2dgrid_which(grid, xc, yc, &x, &y);
    first = MATRIX(grid->startidx, x, y);
    VECTOR(grid->prev)[elem] = 0;
    VECTOR(grid->next)[elem] = first;
    if (first != 0) {
        VECTOR(grid->prev)[first - 1] = elem + 1;
    }
    MATRIX(grid->startidx, x, y) = elem + 1;

    grid->massx += xc;
    grid->massy += yc;
    grid->vertices += 1;
}

void igraph_2dgrid_add2(igraph_2dgrid_t *grid, igraph_integer_t elem) {
    igraph_integer_t x, y;
    igraph_integer_t first;
    igraph_real_t xc, yc;

    xc = MATRIX(*grid->coords, elem, 0);
    yc = MATRIX(*grid->coords, elem, 1);

    /* add to cell */
    igraph_i_2dgrid_which(grid, xc, yc, &x, &y);
    first = MATRIX(grid->startidx, x, y);
    VECTOR(grid->prev)[elem] = 0;
    VECTOR(grid->next)[elem] = first;
    if (first != 0) {
        VECTOR(grid->prev)[first - 1] = elem + 1;
    }
    MATRIX(grid->startidx, x, y) = elem + 1;

    grid->massx += xc;
    grid->massy += yc;
    grid->vertices += 1;
}

void igraph_2dgrid_move(igraph_2dgrid_t *grid, igraph_integer_t elem,
                        igraph_real_t xc, igraph_real_t yc) {
    igraph_integer_t oldx, oldy;
    igraph_integer_t newx, newy;
    igraph_real_t oldxc = MATRIX(*grid->coords, elem, 0);
    igraph_real_t oldyc = MATRIX(*grid->coords, elem, 1);
    igraph_integer_t first;

    xc = oldxc + xc; yc = oldyc + yc;

    igraph_i_2dgrid_which(grid, oldxc, oldyc, &oldx, &oldy);
    igraph_i_2dgrid_which(grid, xc, yc, &newx, &newy);
    if (oldx != newx || oldy != newy) {
        /* remove from this cell */
        if (VECTOR(grid->prev)[elem] != 0) {
            VECTOR(grid->next) [ VECTOR(grid->prev)[elem] - 1 ] =
                VECTOR(grid->next)[elem];
        } else {
            MATRIX(grid->startidx, oldx, oldy) = VECTOR(grid->next)[elem];
        }
        if (VECTOR(grid->next)[elem] != 0) {
            VECTOR(grid->prev)[ VECTOR(grid->next)[elem] - 1 ] =
                VECTOR(grid->prev)[elem];
        }

        /* add to this cell */
        first = MATRIX(grid->startidx, newx, newy);
        VECTOR(grid->prev)[elem] = 0;
        VECTOR(grid->next)[elem] = first;
        if (first != 0) {
            VECTOR(grid->prev)[first - 1] = elem + 1;
        }
        MATRIX(grid->startidx, newx, newy) = elem + 1;
    }

    grid->massx += -oldxc + xc;
    grid->massy += -oldyc + yc;

    MATRIX(*grid->coords, elem, 0) = xc;
    MATRIX(*grid->coords, elem, 1) = yc;

}

void igraph_2dgrid_getcenter(const igraph_2dgrid_t *grid,
                             igraph_real_t *massx, igraph_real_t *massy) {
    *massx = (grid->massx) / (grid->vertices);
    *massy = (grid->massy) / (grid->vertices);
}

igraph_bool_t igraph_2dgrid_in(const igraph_2dgrid_t *grid, igraph_integer_t elem) {
    return VECTOR(grid->next)[elem] != -1;
}

igraph_real_t igraph_2dgrid_sq_dist(const igraph_2dgrid_t *grid,
                                    igraph_integer_t e1, igraph_integer_t e2) {
    igraph_real_t x = MATRIX(*grid->coords, e1, 0) - MATRIX(*grid->coords, e2, 0);
    igraph_real_t y = MATRIX(*grid->coords, e1, 1) - MATRIX(*grid->coords, e2, 1);

    return x * x + y * y;
}

void igraph_2dgrid_reset(igraph_2dgrid_t *grid, igraph_2dgrid_iterator_t *it) {
    /* Search for the first cell containing a vertex */
    it->x = 0; it->y = 0; it->vid = MATRIX(grid->startidx, 0, 0);
    while ( it->vid == 0 && (it->x < grid->stepsx - 1 || it->y < grid->stepsy - 1)) {
        it->x += 1;
        if (it->x == grid->stepsx) {
            it->x = 0; it->y += 1;
        }
        it->vid = MATRIX(grid->startidx, it->x, it->y);
    }
}

igraph_integer_t igraph_2dgrid_next(igraph_2dgrid_t *grid,
                                    igraph_2dgrid_iterator_t *it) {
    igraph_integer_t ret = it->vid;

    if (ret == 0) {
        return 0;
    }

    /* First neighbor */
    it->ncells = -1;
    if (it->x != grid->stepsx - 1) {
        it->ncells += 1;
        it->nx[it->ncells] = it->x + 1;
        it->ny[it->ncells] = it->y;
    }
    if (it->y != grid->stepsy - 1) {
        it->ncells += 1;
        it->nx[it->ncells] = it->x;
        it->ny[it->ncells] = it->y + 1;
    }
    if (it->ncells == 1) {
        it->ncells += 1;
        it->nx[it->ncells] = it->x + 1;
        it->ny[it->ncells] = it->y + 1;
    }
    it->ncells += 1;
    it->nx[it->ncells] = it->x;
    it->ny[it->ncells] = it->y;

    it->nei = VECTOR(grid->next) [ ret - 1 ];
    while (it->ncells > 0 && it->nei == 0 ) {
        it->ncells -= 1;
        it->nei = MATRIX(grid->startidx, it->nx[it->ncells], it->ny[it->ncells]);
    }

    /* Next vertex */
    it->vid = VECTOR(grid->next)[ it->vid - 1 ];
    while ( (it->x < grid->stepsx - 1 || it->y < grid->stepsy - 1) &&
            it->vid == 0) {
        it->x += 1;
        if (it->x == grid->stepsx) {
            it->x = 0; it->y += 1;
        }
        it->vid = MATRIX(grid->startidx, it->x, it->y);
    }

    return ret;
}

igraph_integer_t igraph_2dgrid_next_nei(igraph_2dgrid_t *grid,
                                        igraph_2dgrid_iterator_t *it) {
    igraph_integer_t ret = it->nei;

    if (it->nei != 0) {
        it->nei = VECTOR(grid->next) [ ret - 1 ];
    }
    while (it->ncells > 0 && it->nei == 0 ) {
        it->ncells -= 1;
        it->nei = MATRIX(grid->startidx, it->nx[it->ncells], it->ny[it->ncells]);
    }

    return ret;
}
