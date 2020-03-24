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

#ifndef IGRAPH_EULERIAN_H
#define IGRAPH_EULERIAN_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"
#include "igraph_datatype.h"
#include "igraph_iterators.h"
#include "igraph_attributes.h"
#include "igraph_sparsemat.h"

__BEGIN_DECLS

DECLDIR int igraph_is_eulerian(igraph_t *graph);
DECLDIR int is_eulerian_undirected(igraph_t *graph);
DECLDIR int is_eulerian_directed(igraph_t *graph);
DECLDIR igraph_bool_t check_if_bridge(igraph_t *g, igraph_integer_t start, igraph_integer_t v, long edge);
DECLDIR int print_euler_undirected_implementation(igraph_integer_t start, igraph_t *g, igraph_vector_t *path);
DECLDIR int igraph_euler_path_undirected(igraph_t *graph, igraph_vector_t *path);
DECLDIR int eulerian_path_directed_implementation(igraph_t *graph, igraph_integer_t *start_node, igraph_vector_t *outgoing_list, igraph_vector_t *res);
DECLDIR int igraph_eulerian_path_directed(igraph_t *graph, igraph_vector_t *res);
DECLDIR int igraph_eulerian_paths(igraph_t *graph, igraph_vector_t *res);

__END_DECLS

#endif