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

#ifndef IGRAPH_MIXING_H
#define IGRAPH_MIXING_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_datatype.h"
#include "igraph_vector.h"

__BEGIN_DECLS

DECLDIR int igraph_assortativity_nominal(const igraph_t *graph,
        const igraph_vector_t *types,
        igraph_real_t *res,
        igraph_bool_t directed);

DECLDIR int igraph_assortativity(const igraph_t *graph,
                                 const igraph_vector_t *types1,
                                 const igraph_vector_t *types2,
                                 igraph_real_t *res,
                                 igraph_bool_t directed);

DECLDIR int igraph_assortativity_degree(const igraph_t *graph,
                                        igraph_real_t *res,
                                        igraph_bool_t directed);

__END_DECLS

#endif
