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

#ifndef IGRAPH_SEPARATORS_H
#define IGRAPH_SEPARATORS_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#include "igraph_constants.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_vector_ptr.h"
#include "igraph_datatype.h"
#include "igraph_iterators.h"

__BEGIN_DECLS

int igraph_is_separator(const igraph_t *graph, 
			const igraph_vs_t candidate,
			igraph_bool_t *res);

int igraph_all_minimal_st_separators(const igraph_t *graph, 
				     igraph_vector_ptr_t *separators);

int igraph_is_minimal_separator(const igraph_t *graph,
				const igraph_vs_t candidate, 
				igraph_bool_t *res);

int igraph_minimum_size_separators(const igraph_t *graph,
				   igraph_vector_ptr_t *separators);

__END_DECLS

#endif
