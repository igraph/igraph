/*
   IGraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

#include "igraph_bitset_list.h"
#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_vector.h"

__BEGIN_DECLS

IGRAPH_EXPORT igraph_error_t igraph_reachability(
            const igraph_t *graph,
            igraph_vector_int_t *membership,
            igraph_vector_int_t *csize,
            igraph_integer_t *no_of_components,
            igraph_bitset_list_t *reach,
            igraph_neimode_t mode);

IGRAPH_EXPORT igraph_error_t igraph_count_reachable(
            const igraph_t *graph,
            igraph_vector_int_t *counts,
            igraph_neimode_t mode);


IGRAPH_EXPORT igraph_error_t igraph_transitive_closure(const igraph_t *graph,
                                                       igraph_t* closure);

__END_DECLS

#endif // IGRAPH_REACHABILITY_H
