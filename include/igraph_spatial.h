/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/
#ifndef IGRAPH_SPATIAL_H
#define IGRAPH_SPATIAL_H

#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_error.h"
#include "igraph_matrix.h"

typedef enum {
    IGRAPH_L2_METRIC
} igraph_metric_t;

IGRAPH_EXPORT igraph_error_t igraph_nearest_neighbor_graph(
    igraph_t *graph,
    igraph_matrix_t *points,
    igraph_metric_t metric,
    igraph_integer_t neighbors,
    igraph_real_t cutoff
);

#endif
