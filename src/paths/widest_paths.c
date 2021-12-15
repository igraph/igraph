/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2021 The igraph development team

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


#include "igraph_paths.h"

#include "igraph_adjlist.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_stack.h"
#include "igraph_qsort.h"

#include "core/indheap.h"
#include "core/interruption.h"

#include <string.h>   /* memset */


igraph_error_t igraph_get_widest_paths(const igraph_t *graph,
                                       igraph_vector_ptr_t *vertices,
                                       igraph_vector_ptr_t *edges,
                                       igraph_integer_t from,
                                       igraph_vs_t to,
                                       const igraph_vector_t *weights,
                                       igraph_neimode_t mode,
                                       igraph_vector_int_t *predecessors,
                                       igraph_vector_int_t *inbound_edges) {

    // TODO: implement a modified dijkstra to calculate widest paths from a source to every other node

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_get_widest_path(const igraph_t *graph,
                                      igraph_vector_int_t *vertices,
                                      igraph_vector_int_t *edges,
                                      igraph_integer_t from,
                                      igraph_integer_t to,
                                      const igraph_vector_t *weights,
                                      igraph_neimode_t mode) {

    // TODO: calls igraph_get_widest_paths and returns the one specific path to "to"

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_get_all_shortest_paths(const igraph_t *graph,
                                            igraph_vector_ptr_t *vertices,
                                            igraph_vector_ptr_t *edges,
                                            igraph_vector_int_t *nrgeo,
                                            igraph_integer_t from, igraph_vs_t to,
                                            const igraph_vector_t *weights,
                                            igraph_neimode_t mode) {

    // TODO: implement a modified dijkstra that stores all paths

    return IGRAPH_SUCCESS;
}



igraph_error_t igraph_widest_paths_floyd_warshalls(const igraph_t *graph,
                                   igraph_matrix_t *res,
                                   const igraph_vs_t from,
                                   const igraph_vs_t to,
                                   const igraph_vector_t *weights,
                                   igraph_neimode_t mode) {

    // TODO: implement floyd warshalls to calculate widest paths between all nodes

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_widest_paths_dijkstra(const igraph_t *graph,
                                   igraph_matrix_t *res,
                                   const igraph_vs_t from,
                                   const igraph_vs_t to,
                                   const igraph_vector_t *weights,
                                   igraph_neimode_t mode) {

    // TODO: calls igraph_get_widest_path for each node

    return IGRAPH_SUCCESS;
}
