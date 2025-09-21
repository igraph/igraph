/*
  igraph library.
  Copyright (C) 2021-2025 The igraph development team <igraph@igraph.org>

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  this program. If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef IGRAPH_RANDOM_INTERNAL_H
#define IGRAPH_RANDOM_INTERNAL_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_vector.h"

IGRAPH_BEGIN_C_DECLS

igraph_error_t igraph_i_random_sample_real(
        igraph_vector_t *res, igraph_real_t l, igraph_real_t h,
        igraph_int_t length);

igraph_uint_t igraph_i_get_random_seed(void);

IGRAPH_END_C_DECLS

#endif
