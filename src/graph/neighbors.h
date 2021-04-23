/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_NEIGHBORS_H
#define IGRAPH_NEIGHBORS_H

#include "igraph_constants.h"

__BEGIN_DECLS

IGRAPH_PRIVATE_EXPORT int igraph_i_neighbors(const igraph_t *graph,
                       igraph_vector_t *neis,
                       igraph_integer_t pnode,
                       igraph_neimode_t mode,
                       igraph_loops_t loops,
                       igraph_multiple_t multiple);

IGRAPH_PRIVATE_EXPORT int igraph_i_incident(const igraph_t *graph,
                       igraph_vector_t *eids,
                       igraph_integer_t pnode,
                       igraph_neimode_t mode,
                       igraph_loops_t loops,
                       igraph_multiple_t multiple);

__END_DECLS

#endif /* IGRAPH_NEIGHBORS_H */
