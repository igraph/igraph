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

#ifndef IGRAPH_LAYOUT_MERGE_GRID_H
#define IGRAPH_LAYOUT_MERGE_GRID_H

#include "igraph_decls.h"
#include "igraph_types.h"

__BEGIN_DECLS

/* A type of grid used for merging layouts; each cell is owned by exactly one graph */

typedef struct igraph_i_layout_mergegrid_t {
    long int *data;
    long int stepsx, stepsy;
    igraph_real_t minx, maxx, deltax;
    igraph_real_t miny, maxy, deltay;
} igraph_i_layout_mergegrid_t;

IGRAPH_PRIVATE_EXPORT int igraph_i_layout_mergegrid_init(igraph_i_layout_mergegrid_t *grid,
                                                         igraph_real_t minx, igraph_real_t maxx, long int stepsx,
                                                         igraph_real_t miny, igraph_real_t maxy, long int stepsy);

IGRAPH_PRIVATE_EXPORT void igraph_i_layout_mergegrid_destroy(igraph_i_layout_mergegrid_t *grid);

IGRAPH_PRIVATE_EXPORT int igraph_i_layout_merge_place_sphere(igraph_i_layout_mergegrid_t *grid,
                                                             igraph_real_t x, igraph_real_t y, igraph_real_t r,
                                                             long int id);

long int igraph_i_layout_mergegrid_get(igraph_i_layout_mergegrid_t *grid,
                                       igraph_real_t x, igraph_real_t y);

long int igraph_i_layout_mergegrid_get_sphere(igraph_i_layout_mergegrid_t *g,
                                              igraph_real_t x, igraph_real_t y, igraph_real_t r);

__END_DECLS

#endif
