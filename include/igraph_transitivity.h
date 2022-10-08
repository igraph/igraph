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

#ifndef IGRAPH_TRANSITIVITY_H
#define IGRAPH_TRANSITIVITY_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_constants.h"
#include "igraph_error.h"
#include "igraph_iterators.h"

__BEGIN_DECLS

IGRAPH_EXPORT igraph_error_t igraph_transitivity_undirected(const igraph_t *graph,
                                                 igraph_real_t *res,
                                                 igraph_transitivity_mode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_transitivity_local_undirected(const igraph_t *graph,
                                                       igraph_vector_t *res,
                                                       const igraph_vs_t vids,
                                                       igraph_transitivity_mode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_transitivity_avglocal_undirected(const igraph_t *graph,
                                                          igraph_real_t *res,
                                                          igraph_transitivity_mode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_transitivity_barrat(const igraph_t *graph,
                                             igraph_vector_t *res,
                                             const igraph_vs_t vids,
                                             const igraph_vector_t *weights,
                                             const igraph_transitivity_mode_t mode);
IGRAPH_EXPORT igraph_error_t igraph_ecc(const igraph_t *graph,
                                        igraph_vector_t *res,
                                        igraph_es_t eids,
                                        igraph_integer_t k,
                                        igraph_bool_t offset,
                                        igraph_bool_t normalize);

__END_DECLS

#endif
