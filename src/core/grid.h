/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2020  The igraph development team

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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_CORE_GRID_H
#define IGRAPH_CORE_GRID_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_matrix.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/**
 * 2d grid containing points
 */

typedef struct igraph_2dgrid_t {
    igraph_matrix_t *coords;          /* The current coordinates in the grid */
    igraph_real_t minx, maxx, deltax; /* Minimum and maximum X coordinates and X spacing */
    igraph_real_t miny, maxy, deltay; /* Minimum and maximum Y coordinates and Y spacing */
    igraph_integer_t stepsx, stepsy;  /* Number of cells in the X and Y directions */
    igraph_matrix_int_t startidx;     /* startidx[i, j] is the index of an arbitrary point in that grid cell, plus one; zero means "empty cell" */
    igraph_vector_int_t next;         /* next[i] is the index of the point following point i in the same cell, plus one; zero means "last point" */
    igraph_vector_int_t prev;         /* prev[i] is the index of the point preceding point i in the same cell, plus one; zero means "first point" */
    igraph_real_t massx, massy;       /* The sum of the coordinates */
    igraph_integer_t vertices;        /* Number of active vertices  */
} igraph_2dgrid_t;

igraph_error_t igraph_2dgrid_init(igraph_2dgrid_t *grid, igraph_matrix_t *coords,
                       igraph_real_t minx, igraph_real_t maxx, igraph_real_t deltax,
                       igraph_real_t miny, igraph_real_t maxy, igraph_real_t deltay);
void igraph_2dgrid_destroy(igraph_2dgrid_t *grid);
void igraph_2dgrid_add(igraph_2dgrid_t *grid, igraph_integer_t elem,
                       igraph_real_t xc, igraph_real_t yc);
void igraph_2dgrid_add2(igraph_2dgrid_t *grid, igraph_integer_t elem);
void igraph_2dgrid_move(igraph_2dgrid_t *grid, igraph_integer_t elem,
                        igraph_real_t xc, igraph_real_t yc);
void igraph_2dgrid_getcenter(const igraph_2dgrid_t *grid,
                             igraph_real_t *massx, igraph_real_t *massy);
igraph_bool_t igraph_2dgrid_in(const igraph_2dgrid_t *grid, igraph_integer_t elem);

typedef struct igraph_2dgrid_iterator_t {
    igraph_integer_t vid, x, y;
    igraph_integer_t nei;
    igraph_integer_t nx[4], ny[4], ncells;
} igraph_2dgrid_iterator_t;

void igraph_2dgrid_reset(igraph_2dgrid_t *grid, igraph_2dgrid_iterator_t *it);
igraph_integer_t igraph_2dgrid_next(igraph_2dgrid_t *grid,
                                    igraph_2dgrid_iterator_t *it);
igraph_integer_t igraph_2dgrid_next_nei(igraph_2dgrid_t *grid,
                                        igraph_2dgrid_iterator_t *it);

__END_DECLS

#endif
