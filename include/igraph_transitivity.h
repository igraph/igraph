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

#ifndef IGRAPH_TRANSITIVITY_H
#define IGRAPH_TRANSITIVITY_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_constants.h"
#include "igraph_error.h"
#include "igraph_iterators.h"

IGRAPH_BEGIN_C_DECLS

IGRAPH_EXPORT igraph_error_t igraph_transitivity_undirected(const igraph_t *graph,
                                                 igraph_real_t *res,
                                                 igraph_transitivity_mode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_transitivity_local_undirected(const igraph_t *graph,
                                                       igraph_vector_t *res,
                                                       igraph_vs_t vids,
                                                       igraph_transitivity_mode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_transitivity_avglocal_undirected(const igraph_t *graph,
                                                          igraph_real_t *res,
                                                          igraph_transitivity_mode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_transitivity_barrat(const igraph_t *graph,
                                             igraph_vector_t *res,
                                             igraph_vs_t vids,
                                             const igraph_vector_t *weights,
                                             igraph_transitivity_mode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_ecc(const igraph_t *graph,
                                        igraph_vector_t *res,
                                        igraph_es_t eids,
                                        igraph_int_t k,
                                        igraph_bool_t offset,
                                        igraph_bool_t normalize);

IGRAPH_END_C_DECLS

#endif
