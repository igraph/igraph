/*
   igraph library.
   Copyright (C) 2007-2021  The igraph development team <igraph@igraph.org>

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

#include "centrality/centrality_internal.h"

igraph_bool_t igraph_i_vector_mostly_negative(const igraph_vector_t *vector) {
    /* Many of the centrality measures correspond to the eigenvector of some
     * matrix. When v is an eigenvector, c*v is also an eigenvector, therefore
     * it may happen that all the scores in the eigenvector are negative, in which
     * case we want to negate them since the centrality scores should be positive.
     * However, since ARPACK is not always stable, sometimes it happens that
     * *some* of the centrality scores are small negative numbers. This function
     * helps distinguish between the two cases; it should return true if most of
     * the values are relatively large negative numbers, in which case we should
     * negate the eigenvector.
     */
    igraph_int_t n = igraph_vector_size(vector);
    igraph_real_t mi, ma;

    if (n == 0) {
        return false;
    }

    igraph_vector_minmax(vector, &mi, &ma);

    if (mi >= 0) {
        return false;
    }
    if (ma <= 0) {
        return true;
    }

    /* is the most negative value larger in magnitude than the most positive? */
    return (-mi/ma > 1);
}

/* Normalizes a vector of real numbers such that the largest value, as well as
 * the largest value by magnitude, are 1.0. This is used by functions that
 * produce eigenvector-like centrality values, scaling the largest centrality
 * to 1.0. */
void igraph_i_vector_scale_by_max_abs(igraph_vector_t *vec) {
    const igraph_int_t n = igraph_vector_size(vec);
    igraph_real_t amax = 0;
    igraph_int_t which = 0;

    for (igraph_int_t i = 0; i < n; i++) {
        igraph_real_t tmp;
        tmp = fabs(VECTOR(*vec)[i]);
        if (tmp > amax) {
            amax = tmp;
            which = i;
        }
    }
    if (amax != 0) {
        igraph_vector_scale(vec, 1 / VECTOR(*vec)[which]);
    }
}
