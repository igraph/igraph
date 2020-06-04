/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_FLOW_INTERNAL_H
#define IGRAPH_FLOW_INTERNAL_H

#include "igraph_types.h"
#include "igraph_marked_queue.h"
#include "igraph_estack.h"
#include "igraph_datatype.h"

typedef int igraph_provan_shier_pivot_t(const igraph_t *graph,
                                        const igraph_marked_queue_t *S,
                                        const igraph_estack_t *T,
                                        long int source,
                                        long int target,
                                        long int *v,
                                        igraph_vector_t *Isv,
                                        void *arg);

#endif

