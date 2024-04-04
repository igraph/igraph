/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2022  The igraph development team
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

#include "igraph_bitset_list.h"

#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_types.h"

#define BITSET_LIST
#define BASE_BITSET
#define CUSTOM_INIT_DESTROY
#include "igraph_pmt.h"
#include "typed_list.pmt"
#include "igraph_pmt_off.h"
#undef CUSTOM_INIT_DESTROY
#undef BASE_BITSET
#undef BITSET_LIST

static igraph_error_t igraph_i_bitset_list_init_item(
    const igraph_bitset_list_t* list, igraph_bitset_t* item
) {
    return igraph_bitset_init(item, 0);
}

static igraph_error_t igraph_i_bitset_list_copy_item(
    igraph_bitset_t* dest, const igraph_bitset_t* source
) {
    return igraph_bitset_init_copy(dest, source);
}

static void igraph_i_bitset_list_destroy_item(igraph_bitset_t* item) {
    igraph_bitset_destroy(item);
}
