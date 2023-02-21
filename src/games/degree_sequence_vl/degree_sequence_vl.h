/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_GAMES_DEGREE_SEQUENCE_VL_H
#define IGRAPH_GAMES_DEGREE_SEQUENCE_VL_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_vector.h"

__BEGIN_DECLS

igraph_error_t igraph_degree_sequence_game_vl(igraph_t *graph,
                                   const igraph_vector_int_t *out_seq,
                                   const igraph_vector_int_t *in_seq);

__END_DECLS

#endif /* IGRAPH_GAMES_DEGREE_SEQUENCE_VL_H */
