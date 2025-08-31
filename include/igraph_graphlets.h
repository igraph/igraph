/*
   igraph library.
   Copyright (C) 2013-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_GRAPHLETS_H
#define IGRAPH_GRAPHLETS_H

#include "igraph_decls.h"

#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_vector_list.h"

IGRAPH_BEGIN_C_DECLS

IGRAPH_EXPORT igraph_error_t igraph_graphlets_candidate_basis(const igraph_t *graph,
                                                   const igraph_vector_t *weights,
                                                   igraph_vector_int_list_t *cliques,
                                                   igraph_vector_t *thresholds);

IGRAPH_EXPORT igraph_error_t igraph_graphlets_project(const igraph_t *graph,
                                           const igraph_vector_t *weights,
                                           const igraph_vector_int_list_t *cliques,
                                           igraph_vector_t *Mu, igraph_bool_t startMu,
                                           igraph_int_t niter);

IGRAPH_EXPORT igraph_error_t igraph_graphlets(const igraph_t *graph,
                                   const igraph_vector_t *weights,
                                   igraph_vector_int_list_t *cliques,
                                   igraph_vector_t *Mu, igraph_int_t niter);

IGRAPH_END_C_DECLS

#endif
