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

#ifndef IGRAPH_SCAN_H
#define IGRAPH_SCAN_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_constants.h"
#include "igraph_error.h"
#include "igraph_vector_list.h"

IGRAPH_BEGIN_C_DECLS

IGRAPH_EXPORT igraph_error_t igraph_local_scan_0(const igraph_t *graph, igraph_vector_t *res,
                                      const igraph_vector_t *weights, igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_local_scan_0_them(const igraph_t *us, const igraph_t *them,
                                           igraph_vector_t *res,
                                           const igraph_vector_t *weights_them,
                                           igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_local_scan_1_ecount(const igraph_t *graph, igraph_vector_t *res,
                                             const igraph_vector_t *weights,
                                             igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_local_scan_1_ecount_them(const igraph_t *us, const igraph_t *them,
                                                  igraph_vector_t *res,
                                                  const igraph_vector_t *weights,
                                                  igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_local_scan_k_ecount(const igraph_t *graph, igraph_int_t k,
                                             igraph_vector_t *res,
                                             const igraph_vector_t *weights,
                                             igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_local_scan_k_ecount_them(const igraph_t *us, const igraph_t *them,
                                                  igraph_int_t k, igraph_vector_t *res,
                                                  const igraph_vector_t *weights_them,
                                                  igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_local_scan_neighborhood_ecount(const igraph_t *graph,
                                                        igraph_vector_t *res,
                                                        const igraph_vector_t *weights,
                                                        const igraph_vector_int_list_t *neighborhoods);
IGRAPH_EXPORT igraph_error_t igraph_local_scan_subset_ecount(const igraph_t *graph,
                                                        igraph_vector_t *res,
                                                        const igraph_vector_t *weights,
                                                        const igraph_vector_int_list_t *neighborhoods);
IGRAPH_END_C_DECLS

#endif
