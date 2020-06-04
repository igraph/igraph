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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_GRAPHLETS_H
#define IGRAPH_GRAPHLETS_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_vector_ptr.h"
#include "igraph_interface.h"

__BEGIN_DECLS

DECLDIR int igraph_graphlets_candidate_basis(const igraph_t *graph,
        const igraph_vector_t *weights,
        igraph_vector_ptr_t *cliques,
        igraph_vector_t *thresholds);

DECLDIR int igraph_graphlets_project(const igraph_t *graph,
                                     const igraph_vector_t *weights,
                                     const igraph_vector_ptr_t *cliques,
                                     igraph_vector_t *Mu, igraph_bool_t startMu,
                                     int niter);

DECLDIR int igraph_graphlets(const igraph_t *graph,
                             const igraph_vector_t *weights,
                             igraph_vector_ptr_t *cliques,
                             igraph_vector_t *Mu, int niter);

__END_DECLS

#endif
