/* -*- mode: C -*-  */
/*
   IGraph package.
   Copyright (C) 2006-2021 The igraph development team

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

#include "layout/merge_grid.h"

#include "igraph_memory.h"


static igraph_error_t igraph_i_layout_mergegrid_which(igraph_i_layout_mergegrid_t *grid,
                                    igraph_real_t xc, igraph_real_t yc,
                                    igraph_integer_t *x, igraph_integer_t *y) {
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

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_layout_mergegrid_init(igraph_i_layout_mergegrid_t *grid,
                                              igraph_real_t minx, igraph_real_t maxx, igraph_integer_t stepsx,
                                              igraph_real_t miny, igraph_real_t maxy, igraph_integer_t stepsy) {
    grid->minx = minx;
    grid->maxx = maxx;
    grid->stepsx = stepsx;
    grid->deltax = (maxx - minx) / stepsx;
    grid->miny = miny;
    grid->maxy = maxy;
    grid->stepsy = stepsy;
    grid->deltay = (maxy - miny) / stepsy;

    grid->data = IGRAPH_CALLOC(stepsx * stepsy, igraph_integer_t);
    if (grid->data == 0) {
        IGRAPH_ERROR("Cannot create grid", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    return IGRAPH_SUCCESS;
}

void igraph_i_layout_mergegrid_destroy(igraph_i_layout_mergegrid_t *grid) {
    IGRAPH_FREE(grid->data);
}

#define MAT(i,j) (grid->data[(grid->stepsy)*(j)+(i)])
#define DIST2(x2,y2) (sqrt(pow(x-(x2),2)+pow(y-(y2), 2)))

igraph_error_t igraph_i_layout_merge_place_sphere(igraph_i_layout_mergegrid_t *grid,
                                       igraph_real_t x, igraph_real_t y, igraph_real_t r,
                                       igraph_integer_t id) {
    igraph_integer_t cx, cy;
    igraph_integer_t i, j;

    igraph_i_layout_mergegrid_which(grid, x, y, &cx, &cy);

    MAT(cx, cy) = id + 1;

#define DIST(i,j) (DIST2(grid->minx+(cx+(i))*grid->deltax, \
                         grid->miny+(cy+(j))*grid->deltay))

    for (i = 0; cx + i < grid->stepsx && DIST(i, 0) < r; i++) {
        for (j = 0; cy + j < grid->stepsy && DIST(i, j) < r; j++) {
            MAT(cx + i, cy + j) = id + 1;
        }
    }

#undef DIST
#define DIST(i,j) (DIST2(grid->minx+(cx+(i))*grid->deltax, \
                         grid->miny+(cy-(j)+1)*grid->deltay))

    for (i = 0; cx + i < grid->stepsx && DIST(i, 0) < r; i++) {
        for (j = 1; cy - j > 0 && DIST(i, j) < r; j++) {
            MAT(cx + i, cy - j) = id + 1;
        }
    }

#undef DIST
#define DIST(i,j) (DIST2(grid->minx+(cx-(i)+1)*grid->deltax, \
                         grid->miny+(cy+(j))*grid->deltay))

    for (i = 1; cx - i > 0 && DIST(i, 0) < r; i++) {
        for (j = 0; cy + j < grid->stepsy && DIST(i, j) < r; j++) {
            MAT(cx - i, cy + j) = id + 1;
        }
    }

#undef DIST
#define DIST(i,j) (DIST2(grid->minx+(cx-(i)+1)*grid->deltax, \
                         grid->miny+(cy-(j)+1)*grid->deltay))

    for (i = 1; cx - i > 0 && DIST(i, 0) < r; i++) {
        for (j = 1; cy - j > 0 && DIST(i, j) < r; j++) {
            MAT(cx - i, cy - j) = id + 1;
        }
    }

#undef DIST
#undef DIST2

    return IGRAPH_SUCCESS;
}

igraph_integer_t igraph_i_layout_mergegrid_get(igraph_i_layout_mergegrid_t *grid,
                                       igraph_real_t x, igraph_real_t y) {
    igraph_integer_t cx, cy;
    igraph_integer_t res;

    if (x <= grid->minx || x >= grid->maxx ||
        y <= grid->miny || y >= grid->maxy) {
        res = -1;
    } else {
        igraph_i_layout_mergegrid_which(grid, x, y, &cx, &cy);
        res = MAT(cx, cy) - 1;
    }

    return res;
}

#define DIST2(x2,y2) (sqrt(pow(x-(x2),2)+pow(y-(y2), 2)))

igraph_integer_t igraph_i_layout_mergegrid_get_sphere(igraph_i_layout_mergegrid_t *grid,
        igraph_real_t x, igraph_real_t y, igraph_real_t r) {
    igraph_integer_t cx, cy;
    igraph_integer_t i, j;
    igraph_integer_t ret;

    if (x - r <= grid->minx || x + r >= grid->maxx ||
        y - r <= grid->miny || y + r >= grid->maxy) {
        ret = -1;
    } else {
        igraph_i_layout_mergegrid_which(grid, x, y, &cx, &cy);

        ret = MAT(cx, cy) - 1;

#define DIST(i,j) (DIST2(grid->minx+(cx+(i))*grid->deltax, \
                         grid->miny+(cy+(j))*grid->deltay))

        for (i = 0; ret < 0 && cx + i < grid->stepsx && DIST(i, 0) < r; i++) {
            for (j = 0; ret < 0 && cy + j < grid->stepsy && DIST(i, j) < r; j++) {
                ret = MAT(cx + i, cy + j) - 1;
            }
        }

#undef DIST
#define DIST(i,j) (DIST2(grid->minx+(cx+(i))*grid->deltax, \
                         grid->miny+(cy-(j)+1)*grid->deltay))

        for (i = 0; ret < 0 && cx + i < grid->stepsx && DIST(i, 0) < r; i++) {
            for (j = 1; ret < 0 && cy - j > 0 && DIST(i, j) < r; j++) {
                ret = MAT(cx + i, cy - j) - 1;
            }
        }

#undef DIST
#define DIST(i,j) (DIST2(grid->minx+(cx-(i)+1)*grid->deltax, \
                         grid->miny+(cy+(j))*grid->deltay))

        for (i = 1; ret < 0 && cx - i > 0 && DIST(i, 0) < r; i++) {
            for (j = 0; ret < 0 && cy + j < grid->stepsy && DIST(i, j) < r; j++) {
                ret = MAT(cx - i, cy + j) - 1;
            }
        }

#undef DIST
#define DIST(i,j) (DIST2(grid->minx+(cx-(i)+1)*grid->deltax, \
                         grid->miny+(cy-(j)+1)*grid->deltay))

        for (i = 1; ret < 0 && cx + i > 0 && DIST(i, 0) < r; i++) {
            for (j = 1; ret < 0 && cy + i > 0 && DIST(i, j) < r; j++) {
                ret = MAT(cx - i, cy - j) - 1;
            }
        }

#undef DIST

    }

    return ret;
}

/* int print_grid(igraph_i_layout_mergegrid_t *grid) { */
/*   igraph_integer_t i,j; */

/*   for (i=0; i<grid->stepsx; i++) { */
/*     for (j=0; j<grid->stepsy; j++) { */
/*       printf("%li ", MAT(i,j)-1); */
/*     } */
/*     printf("\n"); */
/*   } */
/* } */
