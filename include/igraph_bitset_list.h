/*
   igraph library.
   Copyright (C) 2024-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_BITSET_LIST_H
#define IGRAPH_BITSET_LIST_H

#include "igraph_bitset.h"
#include "igraph_decls.h"
#include "igraph_error.h"

IGRAPH_BEGIN_C_DECLS

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

IGRAPH_END_C_DECLS

#endif
