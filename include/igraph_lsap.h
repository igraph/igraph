/*
   IGraph library.
   Copyright (C) 2014-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_LSAP_H
#define IGRAPH_LSAP_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_matrix.h"
#include "igraph_vector.h"
#include "igraph_types.h"

IGRAPH_BEGIN_C_DECLS

IGRAPH_EXPORT igraph_error_t igraph_solve_lsap(const igraph_matrix_t *c, igraph_integer_t n,
                      igraph_vector_int_t *p);

IGRAPH_END_C_DECLS

#endif
