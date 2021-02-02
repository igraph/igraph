/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_PSUMTREE_H
#define IGRAPH_PSUMTREE_H

#include "igraph_decls.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/*
 * Defines a partial prefix sum tree which is handy for drawing random numbers
 * from a dynamic discrete distribution. The first part (0,...,offset - 1) of
 * the vector v contains the prefixes of the values contained in the latter part
 * (offset, offset + size - 1) of vector v.
 *
 * More precisely: the part between (offset, offset + size - 1) of vector v
 * contains the values (not necessarily probabilities) corresponding to the
 * individual items. For the part in front of it, it holds that the value at
 * index i (zero-based) is the sum of values at index (2*i + 1) and index
 * (2*i + 2). The item at index zero contains the sum of all values in the
 * slice between (offset, offset + size - 1).
 *
 * In order for the partial prefix sum tree to be suitable for drawing values
 * from a dynamic discrete distribution, it must hold that all the values in
 * the tree must be non-negative. Since this is the only use-case that we use
 * this data structure for, this condition is ensured when the tree is updated.
 */

typedef struct {
    igraph_vector_t v;
    long int size;
    long int offset;
} igraph_psumtree_t;

IGRAPH_EXPORT int igraph_psumtree_init(igraph_psumtree_t *t, long int size);
IGRAPH_EXPORT void igraph_psumtree_reset(igraph_psumtree_t *t);
IGRAPH_EXPORT void igraph_psumtree_destroy(igraph_psumtree_t *t);
IGRAPH_EXPORT igraph_real_t igraph_psumtree_get(const igraph_psumtree_t *t, long int idx);
IGRAPH_EXPORT long int igraph_psumtree_size(const igraph_psumtree_t *t);
IGRAPH_EXPORT int igraph_psumtree_search(const igraph_psumtree_t *t, long int *idx,
                                         igraph_real_t elem);
IGRAPH_EXPORT int igraph_psumtree_update(igraph_psumtree_t *t, long int idx,
                                         igraph_real_t new_value);
IGRAPH_EXPORT igraph_real_t igraph_psumtree_sum(const igraph_psumtree_t *t);

__END_DECLS

#endif
