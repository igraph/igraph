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

#ifndef IGRAPH_MISC_CONVERSION_INTERNAL_H
#define IGRAPH_MISC_CONVERSION_INTERNAL_H

#include "igraph_sparsemat.h"
#include "igraph_types.h"

int igraph_i_normalize_sparsemat(igraph_sparsemat_t *sparsemat,
                                 igraph_bool_t column_wise);

#endif
