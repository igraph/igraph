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

#ifndef IGRAPH_NEIGHBORHOOD_H
#define IGRAPH_NEIGHBORHOOD_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_iterators.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

DECLDIR int igraph_neighborhood_size(const igraph_t *graph, igraph_vector_t *res,
                                     igraph_vs_t vids, igraph_integer_t order,
                                     igraph_neimode_t mode, igraph_integer_t mindist);
DECLDIR int igraph_neighborhood(const igraph_t *graph, igraph_vector_ptr_t *res,
                                igraph_vs_t vids, igraph_integer_t order,
                                igraph_neimode_t mode, igraph_integer_t mindist);
DECLDIR int igraph_neighborhood_graphs(const igraph_t *graph, igraph_vector_ptr_t *res,
                                       igraph_vs_t vids, igraph_integer_t order,
                                       igraph_neimode_t mode,
                                       igraph_integer_t mindist);

__END_DECLS

#endif
