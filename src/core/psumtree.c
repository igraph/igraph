/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   Copyright (C) 2006 Elliot Paquette <Elliot.Paquette05@kzoo.edu>
   Kalamazoo College, 1200 Academy st, Kalamazoo, MI

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

#include "igraph_types.h"
#include "igraph_psumtree.h"
#include "igraph_error.h"
#include "config.h"

#include <math.h>
#include <stdio.h>

static double igraph_i_log2(double f) {
    return log(f) / log(2.0);
}

int igraph_psumtree_init(igraph_psumtree_t *t, long int size) {
    t->size = size;
    t->offset = (long int) (pow(2, ceil(igraph_i_log2(size))) - 1);
    IGRAPH_CHECK(igraph_vector_init((igraph_vector_t *)t, t->offset + t->size));
    return 0;
}

void igraph_psumtree_reset(igraph_psumtree_t *t) {
    igraph_vector_fill(&(t->v), 0);
}

void igraph_psumtree_destroy(igraph_psumtree_t *t) {
    igraph_vector_destroy((igraph_vector_t *)t);
}

igraph_real_t igraph_psumtree_get(const igraph_psumtree_t *t, long int idx) {
    const igraph_vector_t *tree = &t->v;
    return VECTOR(*tree)[t->offset + idx];
}

int igraph_psumtree_search(const igraph_psumtree_t *t, long int *idx,
                           igraph_real_t search) {
    const igraph_vector_t *tree = &t->v;
    long int i = 1;
    long int size = igraph_vector_size(tree);

    while ( 2 * i + 1 <= size) {
        if ( search <= VECTOR(*tree)[i * 2 - 1] ) {
            i <<= 1;
        } else {
            search -= VECTOR(*tree)[i * 2 - 1];
            i <<= 1;
            i += 1;
        }
    }
    if (2 * i <= size) {
        i = 2 * i;
    }

    *idx = i - t->offset - 1;
    return IGRAPH_SUCCESS;
}

int igraph_psumtree_update(igraph_psumtree_t *t, long int idx,
                           igraph_real_t new_value) {
    const igraph_vector_t *tree = &t->v;
    igraph_real_t difference;

    idx = idx + t->offset + 1;
    difference = new_value - VECTOR(*tree)[idx - 1];

    while ( idx >= 1 ) {
        VECTOR(*tree)[idx - 1] += difference;
        idx >>= 1;
    }
    return IGRAPH_SUCCESS;
}

long int igraph_psumtree_size(const igraph_psumtree_t *t) {
    return t->size;
}

igraph_real_t igraph_psumtree_sum(const igraph_psumtree_t *t) {
    return VECTOR(t->v)[0];
}
