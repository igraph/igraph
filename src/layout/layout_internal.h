/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2021 The igraph development team

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

#ifndef IGRAPH_LAYOUT_INTERNAL_H
#define IGRAPH_LAYOUT_INTERNAL_H

#include "igraph_types.h"

#include "layout/merge_grid.h"

__BEGIN_DECLS

IGRAPH_PRIVATE_EXPORT int igraph_i_layout_merge_dla(igraph_i_layout_mergegrid_t *grid,
                                                    long int actg, igraph_real_t *x, igraph_real_t *y, igraph_real_t r,
                                                    igraph_real_t cx, igraph_real_t cy, igraph_real_t startr,
                                                    igraph_real_t killr);

IGRAPH_PRIVATE_EXPORT int igraph_i_layout_sphere_2d(igraph_matrix_t *coords, igraph_real_t *x,
                                                    igraph_real_t *y, igraph_real_t *r);

IGRAPH_PRIVATE_EXPORT int igraph_i_layout_sphere_3d(igraph_matrix_t *coords, igraph_real_t *x,
                                                    igraph_real_t *y, igraph_real_t *z,
                                                    igraph_real_t *r);

IGRAPH_PRIVATE_EXPORT float igraph_i_layout_point_segment_dist2(float v_x, float v_y,
                                                                float u1_x, float u1_y,
                                                                float u2_x, float u2_y);

IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_i_layout_segments_intersect(float p0_x, float p0_y,
                                                                       float p1_x, float p1_y,
                                                                       float p2_x, float p2_y,
                                                                       float p3_x, float p3_y);

__END_DECLS

#endif
