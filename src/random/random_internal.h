/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2021 The igraph development team

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

#ifndef IGRAPH_RANDOM_INTERNAL_H
#define IGRAPH_RANDOM_INTERNAL_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_vector.h"

__BEGIN_DECLS

igraph_error_t igraph_random_sample_real(
        igraph_vector_t *res, igraph_real_t l, igraph_real_t h,
        igraph_integer_t length);

__END_DECLS

#endif
