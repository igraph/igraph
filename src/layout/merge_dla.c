/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2020  The igraph development team

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

#include "igraph_layout.h"

#include "igraph_progress.h"
#include "igraph_random.h"

#include "core/interruption.h"
#include "core/math.h"
#include "layout/merge_grid.h"
#include "layout/layout_internal.h"

/**
 * \function igraph_layout_merge_dla
 * \brief Merges multiple layouts by using a DLA algorithm.
 *
 * \experimental
 *
 * First each layout is covered by a circle. Then the layout of the
 * largest graph is placed at the origin. Then the other layouts are
 * placed by the DLA algorithm, larger ones first and smaller ones
 * last.
 *
 * \param thegraphs Pointer vector containing the graph objects of
 *        which the layouts will be merged.
 * \param coords List of matrices with the 2D layouts of the graphs in \p thegraphs.
 * \param res Pointer to an initialized matrix object, the result will
 *        be stored here. It will be resized if needed.
 * \return Error code.
 *
 * Added in version 0.2.
 *
 * </para><para>
 * Time complexity: TODO.
 */

igraph_error_t igraph_layout_merge_dla(
    const igraph_vector_ptr_t *thegraphs, const igraph_matrix_list_t *coords,
    igraph_matrix_t *res
) {
    igraph_integer_t coords_len = igraph_matrix_list_size(coords);
    igraph_vector_t sizes;
    igraph_vector_t x, y, r;
    igraph_vector_t nx, ny, nr;
    igraph_integer_t allnodes = 0;
    igraph_integer_t i, j;
    igraph_integer_t actg;
    igraph_i_layout_mergegrid_t grid;
    igraph_integer_t jpos = 0;
    igraph_real_t minx, maxx, miny, maxy;
    igraph_real_t area = 0;
    igraph_real_t maxr = 0;
    igraph_integer_t respos;

    /* Graphs are currently not used, only the coordinates */
    IGRAPH_UNUSED(thegraphs);

    IGRAPH_VECTOR_INIT_FINALLY(&sizes, coords_len);
    IGRAPH_VECTOR_INIT_FINALLY(&x, coords_len);
    IGRAPH_VECTOR_INIT_FINALLY(&y, coords_len);
    IGRAPH_VECTOR_INIT_FINALLY(&r, coords_len);
    IGRAPH_VECTOR_INIT_FINALLY(&nx, coords_len);
    IGRAPH_VECTOR_INIT_FINALLY(&ny, coords_len);
    IGRAPH_VECTOR_INIT_FINALLY(&nr, coords_len);

    RNG_BEGIN();

    for (i = 0; i < coords_len; i++) {
        igraph_matrix_t *mat = igraph_matrix_list_get_ptr(coords, i);
        igraph_integer_t size = igraph_matrix_nrow(mat);

        if (igraph_matrix_ncol(mat) != 2) {
            IGRAPH_ERROR("igraph_layout_merge_dla works for 2D layouts only",
                         IGRAPH_EINVAL);
        }

        IGRAPH_ALLOW_INTERRUPTION();
        allnodes += size;
        VECTOR(sizes)[i] = size;
        VECTOR(r)[i] = pow(size, .75);
        area += VECTOR(r)[i] * VECTOR(r)[i];
        if (VECTOR(r)[i] > maxr) {
            maxr = VECTOR(r)[i];
        }

        igraph_i_layout_sphere_2d(mat,
                                  igraph_vector_get_ptr(&nx, i),
                                  igraph_vector_get_ptr(&ny, i),
                                  igraph_vector_get_ptr(&nr, i));
    }
    igraph_vector_order2(&sizes); /* largest first */

    /* 0. create grid */
    minx = miny = -sqrt(5 * area);
    maxx = maxy = sqrt(5 * area);
    igraph_i_layout_mergegrid_init(&grid, minx, maxx, 200,
                                   miny, maxy, 200);
    IGRAPH_FINALLY(igraph_i_layout_mergegrid_destroy, &grid);

    /*   fprintf(stderr, "Ok, starting DLA\n"); */

    /* 1. place the largest  */
    actg = VECTOR(sizes)[jpos++];
    igraph_i_layout_merge_place_sphere(&grid, 0, 0, VECTOR(r)[actg], actg);

    IGRAPH_PROGRESS("Merging layouts via DLA", 0.0, NULL);
    while (jpos < coords_len) {
        IGRAPH_ALLOW_INTERRUPTION();
        /*     fprintf(stderr, "comp: %li", jpos); */
        IGRAPH_PROGRESS("Merging layouts via DLA", (100.0 * jpos) / coords_len, NULL);

        actg = VECTOR(sizes)[jpos++];
        /* 2. random walk, TODO: tune parameters */
        igraph_i_layout_merge_dla(&grid, actg,
                                  igraph_vector_get_ptr(&x, actg),
                                  igraph_vector_get_ptr(&y, actg),
                                  VECTOR(r)[actg], 0, 0,
                                  maxx, maxx + 5);

        /* 3. place sphere */
        igraph_i_layout_merge_place_sphere(&grid, VECTOR(x)[actg], VECTOR(y)[actg],
                                           VECTOR(r)[actg], actg);
    }
    IGRAPH_PROGRESS("Merging layouts via DLA", 100.0, NULL);

    /* Create the result */
    IGRAPH_CHECK(igraph_matrix_resize(res, allnodes, 2));
    respos = 0;
    for (i = 0; i < coords_len; i++) {
        igraph_matrix_t *mat = igraph_matrix_list_get_ptr(coords, i);
        igraph_integer_t size = igraph_matrix_nrow(mat);
        igraph_real_t xx = VECTOR(x)[i];
        igraph_real_t yy = VECTOR(y)[i];
        igraph_real_t rr = VECTOR(r)[i] / VECTOR(nr)[i];
        IGRAPH_ALLOW_INTERRUPTION();
        if (VECTOR(nr)[i] == 0) {
            rr = 1;
        }
        for (j = 0; j < size; j++) {
            MATRIX(*res, respos, 0) = rr * (MATRIX(*mat, j, 0) - VECTOR(nx)[i]);
            MATRIX(*res, respos, 1) = rr * (MATRIX(*mat, j, 1) - VECTOR(ny)[i]);
            MATRIX(*res, respos, 0) += xx;
            MATRIX(*res, respos, 1) += yy;
            ++respos;
        }
    }

    RNG_END();

    igraph_i_layout_mergegrid_destroy(&grid);
    igraph_vector_destroy(&sizes);
    igraph_vector_destroy(&x);
    igraph_vector_destroy(&y);
    igraph_vector_destroy(&r);
    igraph_vector_destroy(&nx);
    igraph_vector_destroy(&ny);
    igraph_vector_destroy(&nr);
    IGRAPH_FINALLY_CLEAN(8);
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_layout_sphere_2d(igraph_matrix_t *coords,
                              igraph_real_t *x, igraph_real_t *y,
                              igraph_real_t *r) {
    igraph_integer_t nodes = igraph_matrix_nrow(coords);
    igraph_integer_t i;
    igraph_real_t xmin, xmax, ymin, ymax;

    xmin = xmax = MATRIX(*coords, 0, 0);
    ymin = ymax = MATRIX(*coords, 0, 1);
    for (i = 1; i < nodes; i++) {

        if (MATRIX(*coords, i, 0) < xmin) {
            xmin = MATRIX(*coords, i, 0);
        } else if (MATRIX(*coords, i, 0) > xmax) {
            xmax = MATRIX(*coords, i, 0);
        }

        if (MATRIX(*coords, i, 1) < ymin) {
            ymin = MATRIX(*coords, i, 1);
        } else if (MATRIX(*coords, i, 1) > ymax) {
            ymax = MATRIX(*coords, i, 1);
        }

    }

    *x = (xmin + xmax) / 2;
    *y = (ymin + ymax) / 2;
    *r = sqrt((xmax - xmin)*(xmax - xmin) +  (ymax - ymin)*(ymax - ymin)) / 2;

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_layout_sphere_3d(igraph_matrix_t *coords,
                              igraph_real_t *x, igraph_real_t *y,
                              igraph_real_t *z, igraph_real_t *r) {
    igraph_integer_t nodes = igraph_matrix_nrow(coords);
    igraph_integer_t i;
    igraph_real_t xmin, xmax, ymin, ymax, zmin, zmax;

    xmin = xmax = MATRIX(*coords, 0, 0);
    ymin = ymax = MATRIX(*coords, 0, 1);
    zmin = zmax = MATRIX(*coords, 0, 2);
    for (i = 1; i < nodes; i++) {

        if (MATRIX(*coords, i, 0) < xmin) {
            xmin = MATRIX(*coords, i, 0);
        } else if (MATRIX(*coords, i, 0) > xmax) {
            xmax = MATRIX(*coords, i, 0);
        }

        if (MATRIX(*coords, i, 1) < ymin) {
            ymin = MATRIX(*coords, i, 1);
        } else if (MATRIX(*coords, i, 1) > ymax) {
            ymax = MATRIX(*coords, i, 1);
        }

        if (MATRIX(*coords, i, 2) < zmin) {
            zmin = MATRIX(*coords, i, 2);
        } else if (MATRIX(*coords, i, 2) > zmax) {
            zmax = MATRIX(*coords, i, 2);
        }

    }

    *x = (xmin + xmax) / 2;
    *y = (ymin + ymax) / 2;
    *z = (zmin + zmax) / 2;
    *r = sqrt( (xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin) +
               (zmax - zmin) * (zmax - zmin) ) / 2;

    return IGRAPH_SUCCESS;
}

#define DIST(x,y) (sqrt(pow((x)-cx,2)+pow((y)-cy,2)))

igraph_error_t igraph_i_layout_merge_dla(igraph_i_layout_mergegrid_t *grid,
                              igraph_integer_t actg, igraph_real_t *x, igraph_real_t *y, igraph_real_t r,
                              igraph_real_t cx, igraph_real_t cy, igraph_real_t startr,
                              igraph_real_t killr) {
    igraph_integer_t sp = -1;
    igraph_real_t angle, len;

    /* The graph is not used, only its coordinates */
    IGRAPH_UNUSED(actg);

    while (sp < 0) {
        /* start particle */
        do {
            angle = RNG_UNIF(0, 2 * M_PI);
            len = RNG_UNIF(.5 * startr, startr);
            *x = cx + len * cos(angle);
            *y = cy + len * sin(angle);
            sp = igraph_i_layout_mergegrid_get_sphere(grid, *x, *y, r);
        } while (sp >= 0);

        while (sp < 0 && DIST(*x, *y) < killr) {
            igraph_real_t nx, ny;
            angle = RNG_UNIF(0, 2 * M_PI);
            len = RNG_UNIF(0, startr / 100);
            nx = *x + len * cos(angle);
            ny = *y + len * sin(angle);
            sp = igraph_i_layout_mergegrid_get_sphere(grid, nx, ny, r);
            if (sp < 0) {
                *x = nx; *y = ny;
            }
        }
    }

    return IGRAPH_SUCCESS;
}
