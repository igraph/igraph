/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2009-2020  Gabor Csardi <csardi.gabor@gmail.com>

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

#ifndef IGRAPH_REACHABILITY_H
#define IGRAPH_REACHABILITY_H

#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_vector.h"

#include "limits.h"

__BEGIN_DECLS

#define BITMASK(b) (1 << ((b) % IGRAPH_INTEGER_SIZE))
#define BITSLOT(b) ((b) / IGRAPH_INTEGER_SIZE)
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b) ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb) ((nb + IGRAPH_INTEGER_SIZE - 1) / IGRAPH_INTEGER_SIZE)

IGRAPH_EXPORT igraph_integer_t igraph_bitset_popcount(igraph_vector_int_t* bitset);

IGRAPH_EXPORT igraph_error_t igraph_reachability_directed(
        const igraph_t *graph,
        igraph_vector_int_t *membership,
        igraph_vector_int_t *csize,
        igraph_integer_t *no_of_components,
        igraph_vector_int_list_t *reach);

__END_DECLS

#endif // IGRAPH_GRAPHICALITY_H
