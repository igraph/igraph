/*
   IGraph library.
   Copyright (C) 2016-2022  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_CLIQUER_H
#define IGRAPH_CLIQUER_H

#include "igraph_decls.h"
#include "igraph_cliques.h"

__BEGIN_DECLS

igraph_error_t igraph_i_cliquer_cliques(const igraph_t *graph, igraph_vector_int_list_t *res,
                             igraph_integer_t min_size, igraph_integer_t max_size);

igraph_error_t igraph_i_cliquer_histogram(const igraph_t *graph, igraph_vector_t *hist,
                               igraph_integer_t min_size, igraph_integer_t max_size);

igraph_error_t igraph_i_cliquer_callback(const igraph_t *graph,
                              igraph_integer_t min_size, igraph_integer_t max_size,
                              igraph_clique_handler_t *cliquehandler_fn, void *arg);

igraph_error_t igraph_i_weighted_cliques(const igraph_t *graph,
                              const igraph_vector_t *vertex_weights, igraph_vector_int_list_t *res,
                              igraph_real_t min_weight, igraph_real_t max_weight, igraph_bool_t maximal);

igraph_error_t igraph_i_largest_weighted_cliques(const igraph_t *graph,
                                      const igraph_vector_t *vertex_weights, igraph_vector_int_list_t *res);

igraph_error_t igraph_i_weighted_clique_number(const igraph_t *graph,
                                    const igraph_vector_t *vertex_weights, igraph_real_t *res);

__END_DECLS

#endif // IGRAPH_CLIQUER_H
