/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA 02139, USA

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

#ifndef IGRAPH_QSORT_H
#define IGRAPH_QSORT_H

#include "igraph_decls.h"

#include <stddef.h>

__BEGIN_DECLS

DECLDIR void igraph_qsort(void *base, size_t nel, size_t width,
                          int (*compar)(const void *, const void *));
DECLDIR void igraph_qsort_r(void *base, size_t nel, size_t width, void *thunk,
                            int (*compar)(void *, const void *, const void *));

__END_DECLS

#endif
