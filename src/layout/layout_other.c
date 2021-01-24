/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2020  The igraph development team

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_layout.h"

/**
 * \section about_layouts
 *
 * <para>Layout generator functions (or at least most of them) try to place the
 * vertices and edges of a graph on a 2D plane or in 3D space in a way
 * which visually pleases the human eye.</para>
 *
 * <para>They take a graph object and a number of parameters as arguments
 * and return an \type igraph_matrix_t, in which each row gives the
 * coordinates of a vertex.</para>
 */

igraph_integer_t igraph_layout_springs(const igraph_t *graph, igraph_matrix_t *res,
                          igraph_real_t mass, igraph_real_t equil, igraph_real_t k,
                          igraph_real_t repeqdis, igraph_real_t kfr, igraph_bool_t repulse) {

    IGRAPH_UNUSED(graph); IGRAPH_UNUSED(res); IGRAPH_UNUSED(mass);
    IGRAPH_UNUSED(equil); IGRAPH_UNUSED(k); IGRAPH_UNUSED(repeqdis);
    IGRAPH_UNUSED(kfr); IGRAPH_UNUSED(repulse);
    IGRAPH_ERROR("Springs layout not implemented", IGRAPH_UNIMPLEMENTED);
    /* TODO */
}
