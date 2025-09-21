/*
   igraph library.
   Copyright (C) 2003-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_SAMPLING_H
#define IGRAPH_SAMPLING_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_matrix.h"
#include "igraph_random.h"

IGRAPH_BEGIN_C_DECLS

IGRAPH_EXPORT igraph_error_t igraph_rng_sample_sphere_surface(
    igraph_rng_t* rng, igraph_int_t dim, igraph_int_t n, igraph_real_t radius,
    igraph_bool_t positive, igraph_matrix_t *res
);

IGRAPH_EXPORT igraph_error_t igraph_rng_sample_sphere_volume(
    igraph_rng_t* rng, igraph_int_t dim, igraph_int_t n, igraph_real_t radius,
    igraph_bool_t positive, igraph_matrix_t *res
);

IGRAPH_EXPORT igraph_error_t igraph_rng_sample_dirichlet(
    igraph_rng_t* rng, igraph_int_t n, const igraph_vector_t *alpha, igraph_matrix_t *res
);

IGRAPH_END_C_DECLS

#endif
