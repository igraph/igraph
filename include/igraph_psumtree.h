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
 */

typedef struct {
    igraph_vector_t v;
    long int size;
    long int offset;
} igraph_psumtree_t;

DECLDIR int igraph_psumtree_init(igraph_psumtree_t *t, long int size);
DECLDIR void igraph_psumtree_reset(igraph_psumtree_t *t);
DECLDIR void igraph_psumtree_destroy(igraph_psumtree_t *t);
DECLDIR igraph_real_t igraph_psumtree_get(const igraph_psumtree_t *t, long int idx);
DECLDIR long int igraph_psumtree_size(const igraph_psumtree_t *t);
DECLDIR int igraph_psumtree_search(const igraph_psumtree_t *t, long int *idx,
                                   igraph_real_t elem);
DECLDIR int igraph_psumtree_update(igraph_psumtree_t *t, long int idx,
                                   igraph_real_t new_value);
DECLDIR igraph_real_t igraph_psumtree_sum(const igraph_psumtree_t *t);

__END_DECLS

#endif
