/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#ifndef IGRAPH_HEAP_H
#define IGRAPH_HEAP_H

#include "igraph_decls.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Heap                                               */
/* -------------------------------------------------- */

/**
 * Heap data type.
 * \ingroup internal
 */

#define BASE_IGRAPH_REAL
#define HEAP_TYPE_MAX
#include "igraph_pmt.h"
#include "igraph_heap_pmt.h"
#include "igraph_pmt_off.h"
#undef HEAP_TYPE_MAX
#define HEAP_TYPE_MIN
#include "igraph_pmt.h"
#include "igraph_heap_pmt.h"
#include "igraph_pmt_off.h"
#undef HEAP_TYPE_MIN
#undef BASE_IGRAPH_REAL

#define BASE_LONG
#define HEAP_TYPE_MAX
#include "igraph_pmt.h"
#include "igraph_heap_pmt.h"
#include "igraph_pmt_off.h"
#undef HEAP_TYPE_MAX
#define HEAP_TYPE_MIN
#include "igraph_pmt.h"
#include "igraph_heap_pmt.h"
#include "igraph_pmt_off.h"
#undef HEAP_TYPE_MIN
#undef BASE_LONG

#define BASE_CHAR
#define HEAP_TYPE_MAX
#include "igraph_pmt.h"
#include "igraph_heap_pmt.h"
#include "igraph_pmt_off.h"
#undef HEAP_TYPE_MAX
#define HEAP_TYPE_MIN
#include "igraph_pmt.h"
#include "igraph_heap_pmt.h"
#include "igraph_pmt_off.h"
#undef HEAP_TYPE_MIN
#undef BASE_CHAR

#define IGRAPH_HEAP_NULL { 0,0,0 }

__END_DECLS

#endif
