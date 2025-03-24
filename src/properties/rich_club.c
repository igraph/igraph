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

#include "igraph_structural.h"
#include "igraph_interface.h"

igraph_error_t igraph_rich_club_density_sequence(const igraph_t *graph,
                                                 const igraph_vector_int_t *vertex_order,
                                                 igraph_vector_t *res) {
    /*
       Takes a vertex ordering (vertex_order) as input and sequentially computes the density of
       the graph, removing vertices in order.
    */

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_rich_club_coefficient(const igraph_t *graph, igraph_vector_t *res) {
    /* Computes the rich-club coefficient as a function of degree. */

    return IGRAPH_SUCCESS;
}
