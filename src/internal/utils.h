/*
   IGraph library.
   Copyright (C) 2008-2021  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef IGRAPH_INTERNAL_UTILS_H
#define IGRAPH_INTERNAL_UTILS_H

#include "igraph_datatype.h"
#include "igraph_iterators.h"
#include "igraph_matrix.h"

igraph_error_t igraph_i_matrix_subset_vertices(
        igraph_matrix_t *m,
        const igraph_t *graph,
        igraph_vs_t from,
        igraph_vs_t to);

#endif /* IGRAPH_INTERNAL_UTILS_H */
