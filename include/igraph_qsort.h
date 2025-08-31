/*
   igraph library.
   Copyright (C) 2011-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_QSORT_H
#define IGRAPH_QSORT_H

#include "igraph_decls.h"

#include <stddef.h>

IGRAPH_BEGIN_C_DECLS

IGRAPH_EXPORT void igraph_qsort(void *base, size_t nel, size_t width,
                                int (*compar)(const void *, const void *));
IGRAPH_EXPORT void igraph_qsort_r(void *base, size_t nel, size_t width, void *thunk,
                                  int (*compar)(void *, const void *, const void *));

IGRAPH_END_C_DECLS

#endif
