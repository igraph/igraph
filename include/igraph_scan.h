/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_SCAN_H
#define IGRAPH_SCAN_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_arpack.h"
#include "igraph_constants.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

DECLDIR int igraph_local_scan_0(const igraph_t *graph, igraph_vector_t *res,
                                const igraph_vector_t *weights, igraph_neimode_t mode);

DECLDIR int igraph_local_scan_0_them(const igraph_t *us, const igraph_t *them,
                                     igraph_vector_t *res,
                                     const igraph_vector_t *weigths_them,
                                     igraph_neimode_t mode);

DECLDIR int igraph_local_scan_1_ecount(const igraph_t *graph, igraph_vector_t *res,
                                       const igraph_vector_t *weights,
                                       igraph_neimode_t mode);

DECLDIR int igraph_local_scan_1_ecount_them(const igraph_t *us, const igraph_t *them,
        igraph_vector_t *res,
        const igraph_vector_t *weights,
        igraph_neimode_t mode);

DECLDIR int igraph_local_scan_k_ecount(const igraph_t *graph, int k,
                                       igraph_vector_t *res,
                                       const igraph_vector_t *weights,
                                       igraph_neimode_t mode);

DECLDIR int igraph_local_scan_k_ecount_them(const igraph_t *us, const igraph_t *them,
        int k, igraph_vector_t *res,
        const igraph_vector_t *weights_them,
        igraph_neimode_t mode);

DECLDIR int igraph_local_scan_neighborhood_ecount(const igraph_t *graph,
        igraph_vector_t *res,
        const igraph_vector_t *weights,
        const igraph_vector_ptr_t *neighborhoods);

__END_DECLS

#endif
