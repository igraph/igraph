/*
   IGraph library.
   Copyright (C) 2022 The igraph development team

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

#include "math/safe_intop.h"

igraph_error_t igraph_i_safe_add(igraph_integer_t a, igraph_integer_t b, igraph_integer_t *res) {
    IGRAPH_SAFE_ADD(a, b, res);
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_safe_mult(igraph_integer_t a, igraph_integer_t b, igraph_integer_t *res) {
    IGRAPH_SAFE_MULT(a, b, res);
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_safe_vector_int_sum(const igraph_vector_int_t *vec, igraph_integer_t *res) {
    igraph_integer_t i, n = igraph_vector_int_size(vec);
    igraph_integer_t sum = 0;
    for (i=0; i < n; ++i) {
        IGRAPH_SAFE_ADD(sum, VECTOR(*vec)[i], &sum);
    }
    *res = sum;
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_safe_vector_int_prod(const igraph_vector_int_t *vec, igraph_integer_t *res) {
    igraph_integer_t i, n = igraph_vector_int_size(vec);
    igraph_integer_t prod = 1;
    for (i=0; i < n; ++i) {
        IGRAPH_SAFE_MULT(prod, VECTOR(*vec)[i], &prod);
    }
    *res = prod;
    return IGRAPH_SUCCESS;
}
