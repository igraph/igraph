/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_PRPACK
#define IGRAPH_PRPACK

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_datatype.h"
#include "igraph_iterators.h"

#include "igraph_interface.h"

__BEGIN_DECLS

int igraph_i_personalized_pagerank_prpack(const igraph_t *graph, igraph_vector_t *vector,
                                          igraph_real_t *value, const igraph_vs_t vids,
                                          igraph_bool_t directed, igraph_real_t damping,
                                          const igraph_vector_t *reset,
                                          const igraph_vector_t *weights);

__END_DECLS

#endif

