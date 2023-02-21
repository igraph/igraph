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

#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_matrix.h"

#include "layout/merge_grid.h"

__BEGIN_DECLS

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_i_layout_merge_dla(igraph_i_layout_mergegrid_t *grid,
                                                    igraph_integer_t actg, igraph_real_t *x, igraph_real_t *y, igraph_real_t r,
                                                    igraph_real_t cx, igraph_real_t cy, igraph_real_t startr,
                                                    igraph_real_t killr);

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_i_layout_sphere_2d(igraph_matrix_t *coords, igraph_real_t *x,
                                                    igraph_real_t *y, igraph_real_t *r);

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_i_layout_sphere_3d(igraph_matrix_t *coords, igraph_real_t *x,
                                                    igraph_real_t *y, igraph_real_t *z,
                                                    igraph_real_t *r);

IGRAPH_PRIVATE_EXPORT igraph_real_t igraph_i_layout_point_segment_dist2(igraph_real_t v_x, igraph_real_t v_y,
                                                                igraph_real_t u1_x, igraph_real_t u1_y,
                                                                igraph_real_t u2_x, igraph_real_t u2_y);

IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_i_layout_segments_intersect(igraph_real_t p0_x, igraph_real_t p0_y,
                                                                       igraph_real_t p1_x, igraph_real_t p1_y,
                                                                       igraph_real_t p2_x, igraph_real_t p2_y,
                                                                       igraph_real_t p3_x, igraph_real_t p3_y);

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_i_umap_fit_ab(igraph_real_t min_dist,
                                                          igraph_real_t *a_p,
                                                          igraph_real_t *b_p);

igraph_error_t igraph_i_layout_random_bounded(
        const igraph_t *graph, igraph_matrix_t *res,
        const igraph_vector_t *minx, const igraph_vector_t *maxx,
        const igraph_vector_t *miny, const igraph_vector_t *maxy);

igraph_error_t igraph_i_layout_random_bounded_3d(
        const igraph_t *graph, igraph_matrix_t *res,
        const igraph_vector_t *minx, const igraph_vector_t *maxx,
        const igraph_vector_t *miny, const igraph_vector_t *maxy,
        const igraph_vector_t *minz, const igraph_vector_t *maxz);

__END_DECLS

#endif
