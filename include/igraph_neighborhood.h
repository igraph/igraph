/*
   igraph library.
   Copyright (C) 2009-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_NEIGHBORHOOD_H
#define IGRAPH_NEIGHBORHOOD_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_graph_list.h"
#include "igraph_iterators.h"
#include "igraph_vector_list.h"

IGRAPH_BEGIN_C_DECLS

IGRAPH_EXPORT igraph_error_t igraph_neighborhood_size(const igraph_t *graph, igraph_vector_int_t *res,
                                           igraph_vs_t vids, igraph_int_t order,
                                           igraph_neimode_t mode, igraph_int_t mindist);
IGRAPH_EXPORT igraph_error_t igraph_neighborhood(const igraph_t *graph, igraph_vector_int_list_t *res,
                                      igraph_vs_t vids, igraph_int_t order,
                                      igraph_neimode_t mode, igraph_int_t mindist);
IGRAPH_EXPORT igraph_error_t igraph_neighborhood_graphs(const igraph_t *graph, igraph_graph_list_t *res,
                                             igraph_vs_t vids, igraph_int_t order,
                                             igraph_neimode_t mode,
                                             igraph_int_t mindist);

IGRAPH_END_C_DECLS

#endif
