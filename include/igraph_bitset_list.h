/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2024  The igraph development team

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

#ifndef IGRAPH_BITSET_LIST_H
#define IGRAPH_BITSET_LIST_H

#include "igraph_bitset.h"
#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* List of graphs                                     */
/* -------------------------------------------------- */

#define BITSET_LIST
#define BASE_BITSET
#include "igraph_pmt.h"
#include "igraph_typed_list_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_BITSET
#undef BITSET_LIST

#define IGRAPH_BITSET_LIST_INIT_FINALLY(v, size) \
    do { IGRAPH_CHECK(igraph_bitset_list_init(v, size)); \
        IGRAPH_FINALLY(igraph_bitset_list_destroy, v); } while (0)

__END_DECLS

#endif
